"""Object-based coherence rule: label connected patches (8-connected) of non-permanent RAPID water, take each patch's
MEDIAN event coherence (pixels with coherence only), and set the whole patch to 0 when the median > THR.
Env PCT (default 50) uses that percentile instead of the median; env MIN_KM2 (default 0) only removes patches at least that large.
Patches with no coherence are kept. Permanent water and other classes are untouched.
usage: coh_rule_obj.py <rapid.tif> <vv_linear.tif> <coh.tif (b2 = event)> <out_prefix> <nonperm classes> <thr> [thr ...]
Per threshold writes <out_prefix>_obj<thr>.tif (final map), <out_prefix>_obj<thr>_overview.png, _zooms.png, _tiles.csv."""
import sys, os, numpy as np, pandas as pd
from osgeo import gdal
from scipy import ndimage
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
gdal.UseExceptions()
RF, VF, CF, OUT, NP = sys.argv[1:6]; THRS = [float(x) for x in sys.argv[6:]]; NONP = [int(x) for x in NP.split(',')]
PCT = float(os.environ.get('PCT', 50)); MINKM2 = float(os.environ.get('MIN_KM2', 0))
r = gdal.Open(RF); gt = r.GetGeoTransform(); X, Y = r.RasterXSize, r.RasterYSize; m = r.ReadAsArray()
c = gdal.Open(CF); ev = c.GetRasterBand(2).ReadAsArray(); vds = gdal.Open(VF)
b = (gt[0], gt[3] + Y * gt[5], gt[0] + X * gt[1], gt[3])
lat = gt[3] + gt[5] * (np.arange(Y) + .5); lon = gt[0] + gt[1] * (np.arange(X) + .5)
px = (gt[1] * 111.32) ** 2 * np.cos(np.radians(lat))[:, None]; A = lambda s: float((s * px).sum())
nonp = np.isin(m, NONP)
lab, n = ndimage.label(nonp, structure=np.ones((3, 3), bool)); print(f'{n} patches, {A(nonp):.1f} km2 non-permanent', flush=True)
# per-patch median coherence over pixels that have coherence
sel = nonp & np.isfinite(ev); li = lab[sel]; vals = ev[sel]
o = np.lexsort((vals, li)); li, vals = li[o], vals[o]
starts = np.r_[0, np.flatnonzero(np.diff(li)) + 1]; ends = np.r_[starts[1:], len(li)]
med = np.full(n + 1, np.nan, np.float32)
if PCT == 50: med[li[starts]] = 0.5 * (vals[starts + (ends - starts - 1) // 2] + vals[starts + (ends - starts) // 2])
else: med[li[starts]] = vals[starts + np.floor((ends - starts - 1) * PCT / 100).astype(int)]
size = np.bincount(lab.ravel(), minlength=n + 1); size[0] = 0
medpx = med[lab]   # per-pixel value of its patch median (NaN outside patches / no coherence)
pa = (gt[1] * 111.32) ** 2 * np.cos(np.radians(np.mean(lat)))
bigenough = size * pa >= MINKM2
print(f'patches with coherence: {np.isfinite(med[1:]).sum()}; patch median coh p10/50/90 (area-weighted): ' +
      '/'.join(f'{q:.2f}' for q in np.percentile(medpx[nonp & np.isfinite(medpx)], [10, 50, 90])), flush=True)
for T in THRS:
    rmlab = np.isfinite(med) & (med > T) & bigenough; rmlab[0] = False; rm = rmlab[lab] & nonp
    pix = nonp & np.isfinite(ev) & (ev > T)
    print(f'\nTHR {T:g} (patch p{PCT:g}, min {MINKM2:g} km2): removed {A(rm):.1f} km2 ({100*A(rm)/A(nonp):.1f}%), {rmlab.sum()} of {n} patches; '
          f'pixel rule would remove {A(pix):.1f} km2 ({100*A(pix)/A(nonp):.1f}%)')
    for k in NONP:
        s = m == k
        if s.any(): print(f'  class {k:3d} {A(s):7.1f} km2  removed {100*A(s&rm)/A(s):5.1f}%  (pixel rule {100*A(s&pix)/A(s):5.1f}%)')
    big = np.argsort(size)[::-1][:8]
    print('  largest patches: ' + '; '.join(f'{size[i]*pa:.1f} km2 med {med[i]:.2f}{" REMOVED" if rmlab[i] else ""}' for i in big if size[i] > 0))
    fo = m.copy(); fo[rm] = 0; tag = f'{OUT}_obj{T:g}' + (f'_p{PCT:g}' if PCT != 50 else '') + (f'_min{MINKM2:g}' if MINKM2 else '')
    d = gdal.GetDriverByName('GTiff').Create(f'{tag}.tif', X, Y, 1, gdal.GDT_Byte, ['COMPRESS=DEFLATE', 'TILED=YES'])
    d.SetGeoTransform(gt); d.SetProjection(r.GetProjection()); d.GetRasterBand(1).WriteArray(fo)
    nd = r.GetRasterBand(1).GetNoDataValue(); d.GetRasterBand(1).SetNoDataValue(nd) if nd is not None else None; d = None
    TT = 0.25; rows = []
    for la in np.arange(np.floor(b[1] / TT) * TT, b[3], TT):
        ri = (lat >= la) & (lat < la + TT)
        if not ri.any(): continue
        for lo in np.arange(np.floor(b[0] / TT) * TT, b[2], TT):
            ci = (lon >= lo) & (lon < lo + TT)
            if not ci.any(): continue
            sub = np.ix_(ri, ci); nn = float((nonp[sub] * px[ri]).sum())
            if nn < 0.2: continue
            rows.append(dict(lon=lo, lat=la, nonperm_km2=round(nn, 2), removed_km2=round(float((rm[sub] * px[ri]).sum()), 2)))
    t = pd.DataFrame(rows); t['removed_pct'] = (100 * t.removed_km2 / t.nonperm_km2).round(1); t = t.sort_values('nonperm_km2', ascending=False)
    t.to_csv(f'{tag}_tiles.csv', index=False)
    st = max(1, max(X, Y) // 1600); h, w = Y // st, X // st
    blk = lambda a: a[:h*st, :w*st].reshape(h, st, w, st).sum((1, 3))
    vv = vds.ReadAsArray(buf_xsize=w, buf_ysize=h)
    g = np.clip((10 * np.log10(np.clip(vv, 1e-6, None)) + 22) / 20, 0, 1); g[~np.isfinite(g)] = 0; rgb = np.dstack([g, g, g]) * 0.6 + 0.2
    mk250 = m == 250; mk100 = m == 100; mkoth = (m > 100) & (m < 200)   # mask labels: 250 black, 100 white, 101-199 (106-115) pink
    for cnt, col in ((blk(m == 1), (0.15, 0.3, 1)), (blk(mk100), (1, 1, 1)), (blk(mkoth), (1, 0.2, 0.8)), (blk(mk250), (0, 0, 0)), (blk(nonp & ~rm), (0, 0.9, 1)), (blk(rm), (1, 0.15, 0.15))):
        a = np.clip(cnt / (st * st) * 6, 0, 1)[..., None] * 0.95; rgb = rgb * (1 - a) + a * np.array(col)
    fig, ax = plt.subplots(figsize=(18, 18 * h / w + 1), constrained_layout=True)
    ax.imshow(rgb, extent=(b[0], b[2], b[1], b[3]))
    for _, q in t.head(6).iterrows(): ax.add_patch(plt.Rectangle((q.lon, q.lat), TT, TT, fill=False, ec='yellow', lw=1.2))
    ax.set_title(f'patch rule: red = patch removed (p{PCT:g} event coh > {T:g}, patch >= {MINKM2:g} km2), cyan kept; blue permanent; black 250; white 100; pink 106-115; boxes = zooms')
    fig.savefig(f'{tag}_overview.png', dpi=90); plt.close(fig)
    fig, axs = plt.subplots(2, 6, figsize=(30, 10.5), constrained_layout=True)
    for j, (_, q) in enumerate(t.head(6).iterrows()):
        c0, c1 = max(int((q.lon - gt[0]) / gt[1]), 0), min(int((q.lon + TT - gt[0]) / gt[1]), X)
        r0, r1 = max(int((q.lat + TT - gt[3]) / gt[5]), 0), min(int((q.lat - gt[3]) / gt[5]), Y)
        v = vds.ReadAsArray(c0, r0, c1 - c0, r1 - r0)
        gg = np.clip((10 * np.log10(np.clip(v, 1e-6, None)) + 22) / 20, 0, 1); gg[~np.isfinite(gg)] = 0; im = np.dstack([gg] * 3)
        mm = m[r0:r1, c0:c1]; nn = np.isin(mm, NONP); rr = rm[r0:r1, c0:c1]
        for msk, col in ((mm == 1, (0.15, 0.3, 1)), (mm == 100, (1, 1, 1)), ((mm > 100) & (mm < 200), (1, 0.2, 0.8)), (mm == 250, (0, 0, 0)), (nn & ~rr, (0, 0.9, 1)), (rr, (1, 0.15, 0.15))):
            im[msk] = 0.35 * im[msk] + 0.65 * np.array(col)
        e = (gt[0] + c0 * gt[1], gt[0] + c1 * gt[1], gt[3] + r1 * gt[5], gt[3] + r0 * gt[5])
        axs[0, j].imshow(gg, cmap='gray', extent=e); axs[0, j].set_title(f'VV dB  {q.lon:.2f}, {q.lat:.2f}')
        axs[1, j].imshow(im, extent=e); axs[1, j].set_title(f'{q.nonperm_km2} km2 non-perm; red removed {q.removed_pct}%')
    fig.savefig(f'{tag}_zooms.png', dpi=65); plt.close(fig); print('  ->', f'{tag}.tif, _overview.png, _zooms.png', flush=True)
