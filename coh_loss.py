"""Mosaic HyP3 burst coherence onto the RAPID grid and compute coherence loss.

  python coh_loss.py <project_dir> <RAPID_reference.tif> <out_prefix>

Writes <out_prefix>_coh.tif (float32, NaN nodata):
    band 1  baseline coherence (pre2, pre1)
    band 2  event coherence    (pre1, co)
    band 3  coherence loss = baseline - event
and <out_prefix>_loss.png, a quick look of band 3.
"""
import glob
import os
import sys

import numpy as np
import rasterio
from osgeo import gdal

gdal.UseExceptions()
project, ref_path, prefix = sys.argv[1:4]

with rasterio.open(ref_path) as ref:
    rows, cols = ref.height, ref.width
    crs_wkt, bounds = ref.crs.to_wkt(), tuple(ref.bounds)
    profile = ref.profile.copy()
    ref_invalid = ref.read_masks(1) == 0


def mosaic(kind):
    files = sorted(glob.glob(os.path.join(project, kind, "*", "*_corr.tif")))
    print(f"{kind}: {len(files)} burst coherence files")
    # Bilinear onto the RAPID grid; burst overlaps take the later file.
    ds = gdal.Warp("", files, format="MEM", dstSRS=crs_wkt, outputBounds=bounds,
                   width=cols, height=rows, resampleAlg="bilinear",
                   srcNodata=0, dstNodata=np.nan, outputType=gdal.GDT_Float32,
                   multithread=True)
    arr = ds.ReadAsArray()
    arr[ref_invalid] = np.nan
    return arr


base = mosaic("baseline")
event = mosaic("event")
loss = base - event

profile.update(driver="GTiff", count=3, dtype="float32", nodata=np.nan,
               compress="deflate", predictor=3, tiled=True, blockxsize=256,
               blockysize=256, BIGTIFF="IF_SAFER")
with rasterio.open(f"{prefix}_coh.tif", "w", **profile) as dst:
    for i, (arr, desc) in enumerate(((base, "baseline_coherence_pre2_pre1"),
                                     (event, "event_coherence_pre1_co"),
                                     (loss, "coherence_loss_baseline_minus_event")), 1):
        dst.write(arr.astype(np.float32), i)
        dst.set_band_description(i, desc)

ok = np.isfinite(loss)
for name, arr in (("baseline", base), ("event", event), ("loss", loss)):
    q = np.nanpercentile(arr[ok], [5, 25, 50, 75, 95])
    print(f"{name:8s} p5/25/50/75/95: " + " ".join(f"{v:.2f}" for v in q))
print(f"coverage of RAPID data: {100 * ok.sum() / (~ref_invalid).sum():.1f}%")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

step = max(1, max(rows, cols) // 3000)
fig, axes = plt.subplots(1, 3, figsize=(21, 8), constrained_layout=True)
extent = (bounds[0], bounds[2], bounds[1], bounds[3])
for ax, arr, title, cmap, vmin, vmax in (
        (axes[0], base, "Baseline coherence (pre2, pre1)", "gray", 0, 1),
        (axes[1], event, "Event coherence (pre1, co)", "gray", 0, 1),
        (axes[2], loss, "Coherence loss (baseline - event)", "RdBu_r", -0.6, 0.6)):
    im = ax.imshow(arr[::step, ::step], cmap=cmap, vmin=vmin, vmax=vmax, extent=extent,
                   interpolation="nearest")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.7)
fig.savefig(f"{prefix}_loss.png", dpi=110)
print(f"-> {prefix}_coh.tif, {prefix}_loss.png")
