"""Where snow is physically possible, from 2 m temperature downscaled to a DEM.

  python snow_possible.py <S1_NAME> <REFERENCE_RASTER> <OUT_TIF> [days=30] [tmin_c=5] [--temp-only]

--temp-only fetches the temperature into <OUT_TIF stem>_temp.npz and stops; the full
run reuses it.

Window: the days before the S1 acquisition up to the acquisition time. Temperature is
NOAA GFS 0.25 deg 6-hourly minimum 2 m temperature (TMIN, a few hours behind real
time; AWS archive from 2021-01-01). Windows starting before the archive use ERA5
hourly 2 m temperature (Google ARCO copy, ~5 days behind real time) instead. The
temperature is moved from the model orography to sea level with a 6.5 K/km lapse
rate, the minimum is taken, and the result is moved to the Copernicus DEM (GLO-90)
height.

Output values on the reference grid:
    1    snow possible: downscaled minimum temperature <= tmin_c
    0    snow impossible
    255  reference is nodata
The last line printed is SNOW_POSSIBLE or SNOW_IMPOSSIBLE for the whole scene; when it
is SNOW_IMPOSSIBLE the S2 snow mask is not needed.
"""
import os
import re
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
import urllib.error
import urllib.request
from datetime import datetime, timedelta

import numpy as np
import planetary_computer
import pystac_client
import rasterio
import xarray as xr
import zarr
from osgeo import gdal
from rasterio.warp import transform_bounds

gdal.UseExceptions()
ERA5 = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{c:%Y%m%d}/{c:%H}/{sub}gfs.t{c:%H}z.pgrb2.0p25.f{f:03d}"
GFS_START = datetime(2021, 1, 1)
LAPSE = 0.0065  # K per m
G = 9.80665
MARGIN = 0.5  # degrees around the scene, so bilinear has neighbours

temp_only = "--temp-only" in sys.argv
argv = [a for a in sys.argv if a != "--temp-only"]
name, ref_path, out_path = argv[1:4]
days = int(argv[4]) if len(argv) > 4 else 30
tmin_c = float(argv[5]) if len(argv) > 5 else 5.0
temp_cache = os.path.splitext(out_path)[0] + "_temp.npz"
started = time.perf_counter()

with rasterio.open(ref_path) as ref:
    rows, cols = ref.height, ref.width
    crs_wkt, bounds = ref.crs.to_wkt(), tuple(ref.bounds)
    w, s, e, n = transform_bounds(ref.crs, "EPSG:4326", *ref.bounds, densify_pts=21)
    profile = ref.profile.copy()
    ref_invalid = ref.read_masks(1) == 0

when = datetime.strptime(re.search(r"_(\d{8}T\d{6})_", name).group(1), "%Y%m%dT%H%M%S")
start = when - timedelta(days=days)


def fetch_era5():
    """Sea-level min temperature from ERA5, its grid, and the last hour ERA5 has."""
    # One global chunk per hour: fetch the hours concurrently, without dask.
    zarr.config.set({"async.concurrency": 64})
    ds = xr.open_zarr(ERA5, chunks=None, storage_options=dict(token="anon"))
    box = dict(latitude=slice(n + MARGIN, s - MARGIN))
    # ERA5 longitudes run 0..360 in 0.25 steps; contiguous slices (two across 0/360)
    lo, hi = (w - MARGIN) % 360, (e + MARGIN) % 360
    i0, i1 = int(np.floor(lo / 0.25)), int(np.ceil(hi / 0.25)) + 1
    parts = [slice(i0, i1)] if i0 < i1 else [slice(i0, 1440), slice(0, i1)]
    t2m = [ds["2m_temperature"].sel(time=slice(start, when), **box).isel(longitude=p)
           for p in parts]
    t2m = t2m[0] if len(t2m) == 1 else xr.concat(t2m, "longitude")
    arr = t2m.values
    have = ~np.isnan(arr).all(axis=(1, 2))  # hours past the ERA5 release are NaN
    times = t2m.time.values[have]
    print(f"ERA5: {have.sum()} of {have.size} hours {start:%Y-%m-%d %H:%M} .. "
          f"{when:%Y-%m-%d %H:%M}, {t2m.sizes['latitude']}x{t2m.sizes['longitude']} cells")
    # static field stored on the time axis: one hour is enough
    oro = xr.concat([ds["geopotential_at_surface"].sel(time=when, method="nearest").sel(**box)
                     .isel(longitude=p) for p in parts], "longitude").values / G
    lons = (t2m.longitude.values + 180) % 360 - 180
    order = np.argsort(lons)  # put a dateline-crossing box in increasing longitude
    t0 = (np.nanmin(arr[have], axis=0) - 273.15 + LAPSE * oro) if have.any() else \
        np.full(oro.shape, np.inf)
    last = (times[-1].astype("datetime64[s]").astype(datetime) if have.any()
            else start)
    return t0[:, order], t2m.latitude.values, lons[order], last


def gfs_field(cycle, fcst, element):
    """One GRIB message of a GFS file, by byte range from its .idx, as a GDAL dataset."""
    for sub in ("atmos/", ""):  # the atmos/ level appeared in March 2021
        url = GFS.format(c=cycle, f=fcst, sub=sub)
        try:
            idx = urllib.request.urlopen(url + ".idx", timeout=60).read().decode().splitlines()
            break
        except urllib.error.HTTPError:
            continue
    else:
        raise KeyError(f"no GFS file for {cycle:%Y-%m-%d %H}z f{fcst:03d}")
    for i, line in enumerate(idx):
        if element in line:
            first = int(line.split(":")[1])
            last = int(idx[i + 1].split(":")[1]) - 1 if i + 1 < len(idx) else ""
            break
    else:
        raise KeyError(f"{element} not in {url}.idx")
    req = urllib.request.Request(url, headers={"Range": f"bytes={first}-{last}"})
    tmp = tempfile.NamedTemporaryFile(suffix=".grb2", delete=False)
    tmp.write(urllib.request.urlopen(req, timeout=120).read())
    tmp.close()
    ds = gdal.Translate("", tmp.name, format="MEM")  # -180..180, 90..-90, 0.25 deg
    os.unlink(tmp.name)
    return ds.ReadAsArray().astype(np.float32)


def box_grid():
    """0.25 deg cell centres around the scene, the GFS / ERA5 grid."""
    lats = np.arange(np.ceil((n + MARGIN) * 4), np.floor((s - MARGIN) * 4) - 1, -1) / 4
    lons = np.arange(np.floor((w - MARGIN) * 4), np.ceil((e + MARGIN) * 4) + 1) / 4
    return lats, (lons + 180) % 360 - 180


def fetch_gfs(lats, lons, after):
    """Sea-level min temperature from GFS 6-hourly TMIN for (after, when]."""
    crop = np.ix_(np.rint((90 - lats) * 4).astype(int),
                  np.rint((lons + 180) * 4).astype(int) % 1440)

    def window(cycle):
        # TMIN of f006 covers (cycle, cycle+6h]; fall back to the previous cycle's f012
        for cyc, f in ((cycle, 6), (cycle - timedelta(hours=6), 12)):
            try:
                tmin = gfs_field(cyc, f, ":TMIN:2 m above ground:")[crop]
                oro = gfs_field(cyc, f, ":HGT:surface:")[crop]
                return tmin + LAPSE * oro
            except (urllib.error.URLError, KeyError):
                continue
        print(f"GFS: nothing for {cycle:%Y-%m-%d %H}z")
        return None

    first = datetime(after.year, after.month, after.day, after.hour // 6 * 6)
    cycles = [first + timedelta(hours=6 * i)
              for i in range(int((when - first) / timedelta(hours=6)) + 1)]
    with ThreadPoolExecutor(16) as pool:
        fields = [f for f in pool.map(window, cycles) if f is not None]
    print(f"GFS TMIN: {len(fields)} of {len(cycles)} 6-hour windows "
          f"{cycles[0]:%Y-%m-%d %H}z .. {cycles[-1]:%Y-%m-%d %H}z")
    return np.min(fields, axis=0) if fields else np.full(crop[0].shape[:1] + (lons.size,), np.inf)


if not os.path.exists(temp_cache):
    if start >= GFS_START:
        lats, lons = box_grid()
        t0, sources = fetch_gfs(lats, lons, start), "GFS"
    else:
        t0, lats, lons, last = fetch_era5()
        sources = "ERA5"
        if when - last >= timedelta(hours=1):
            print(f"ERA5 ends {last:%Y-%m-%d %H:%M}; GFS for the rest")
            t0 = np.fmin(t0, fetch_gfs(lats, lons, last))
            sources = "ERA5+GFS"
    np.savez(temp_cache, t0=t0, lats=lats, lons=lons, sources=sources)
if temp_only:
    print(f"-> {temp_cache} ({time.perf_counter() - started:.0f} s)")
    sys.exit()
cache = np.load(temp_cache)
t0, lats, lons = cache["t0"], cache["lats"], cache["lons"]
print(f"min 2m temperature at sea level over the box ({cache['sources']}): "
      f"{np.nanmin(t0):.1f} .. {np.nanmax(t0):.1f} C")


def to_grid(arr):
    """Bilinear from the 0.25 deg cell centres onto the reference grid."""
    src = gdal.GetDriverByName("MEM").Create("", arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    src.SetGeoTransform((lons[0] - 0.125, 0.25, 0, lats[0] + 0.125, 0, -0.25))
    src.SetProjection("EPSG:4326")
    src.GetRasterBand(1).WriteArray(arr.astype(np.float32))
    out = gdal.Warp("", src, format="MEM", dstSRS=crs_wkt, outputBounds=bounds,
                    width=cols, height=rows, resampleAlg="bilinear")
    return out.ReadAsArray()


catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
)
dem_hrefs = ["/vsicurl/" + it.assets["data"].href for it in
             catalog.search(collections=["cop-dem-glo-90"], bbox=(w, s, e, n)).items()]
print(f"Copernicus DEM GLO-90 tiles: {len(dem_hrefs)}")
if dem_hrefs:
    dem = gdal.Warp("", dem_hrefs, format="MEM", dstSRS=crs_wkt, outputBounds=bounds,
                    width=cols, height=rows, resampleAlg="bilinear", dstNodata=np.nan,
                    outputType=gdal.GDT_Float32, multithread=True).ReadAsArray()
else:  # open ocean: no DEM tile, sea level
    dem = np.zeros((rows, cols), np.float32)
dem = np.where(np.isfinite(dem), dem, 0)

tmin_dem = to_grid(t0) - LAPSE * dem
possible = tmin_dem <= tmin_c
out = np.where(ref_invalid, 255, possible).astype(np.uint8)
profile.update(driver="GTiff", count=1, dtype="uint8", nodata=255, compress="deflate",
               tiled=True, blockxsize=256, blockysize=256)
with rasterio.open(out_path, "w", **profile) as dst:
    dst.write(out, 1)
    dst.update_tags(tmin_c=tmin_c, days=days, lapse_k_per_m=LAPSE)

valid = ~ref_invalid
print(f"downscaled min temperature in scene: {np.nanmin(tmin_dem[valid]):.1f} .. "
      f"{np.nanmax(tmin_dem[valid]):.1f} C; snow possible: "
      f"{100 * possible[valid].mean():.2f}% of {valid.sum()} px")
print(f"-> {out_path} ({time.perf_counter() - started:.0f} s)")
print("SNOW_POSSIBLE" if possible[valid].any() else "SNOW_IMPOSSIBLE")
