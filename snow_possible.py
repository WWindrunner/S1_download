"""Where snow or ice is physically possible at the S1 time, from model weather and
surface state.

  python snow_possible.py <S1_NAME> <REFERENCE_RASTER> <OUT_TIF> [t_c=0] [--fetch-only] [--gfs]

--fetch-only fetches the model fields into <OUT_TIF stem>_model.npz and stops; the full
run reuses it. --gfs uses GFS even when ERA5 is available (for checks against ERA5).

Model fields at the hour closest to the S1 acquisition. ERA5 hourly (Google ARCO copy)
when it has that hour: 2 m temperature, snow depth, lake ice depth where the lake cover
is > 5%, and sea ice cover. Otherwise (recent scenes, ERA5 lags by days) and from
2021-01-01 on, NOAA GFS 0.25 deg hourly output of the cycle before the acquisition,
forecast hour = acquisition time rounded to the hour (AWS archive, a few hours behind
real time): 2 m temperature (TMP), snow depth (SNOD) and ice thickness (ICETK). Snow counts from 1 cm depth (ERA5
water equivalent / snow density): ERA5 has 1-1.5 mm of snow on the flooded Red River
valley of the scene below, which would otherwise make it snow-possible.
ERA5 is preferred because its lake model has ice on more lakes: on the Red River scene
below GFS has no ice on many small lakes east of the valley that are frozen.

Real-time and research product: a run within about a week of the S1 time uses GFS
(real-time product). Running the same command again once ERA5 has the hour replaces the
GFS cache with ERA5 (research product); then rerun snow_latest_clear.py (if not done yet
and now SNOW_POSSIBLE) and Post_masks/apply_masks.py. RAPID itself is not rerun. The
output tif has the tags source (ERA5 / GFS) and product (research / realtime); the
second-to-last line printed is "PRODUCT research" or "PRODUCT realtime <date>", the date
from which a rerun gets ERA5.
The 2 m temperature is moved from the model orography to sea level with a 6.5 K/km
lapse rate and then to the Copernicus DEM (GLO-90) height. Snow and ice are taken per
0.25 deg cell (nearest); ice is grown by one cell, because the model grid misses small
lakes (it sees Upper Red Lake and Leech Lake, Minnesota, as land).

Lakes can still be frozen when the air is well above 0 C (Red River scene
S1A_IW_GRDH_1SDV_20220501T001528_20220501T001603_043010_05229E_121A: 2.5-5 C over
Lake of the Woods and Upper / Lower Red Lake, all frozen in Planet imagery, 0.4-0.6 m
of model ice), so the temperature alone is not enough.

Output values on the reference grid:
    1    snow or ice possible: downscaled 2 m temperature <= t_c, or the model has
         snow on the ground or ice in the cell (or a neighbouring cell for ice)
    0    snow and ice impossible
    255  reference is nodata
The last line printed is SNOW_POSSIBLE or SNOW_IMPOSSIBLE for the whole scene; when it
is SNOW_IMPOSSIBLE the S2 snow mask is not needed.
"""
import os
import re
import sys
import tempfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timedelta

import numpy as np
import planetary_computer
from scipy.ndimage import maximum_filter
import pystac_client
import rasterio
import xarray as xr
from osgeo import gdal
from rasterio.warp import transform_bounds

gdal.UseExceptions()
ERA5 = "gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3"
GFS = "https://noaa-gfs-bdp-pds.s3.amazonaws.com/gfs.{c:%Y%m%d}/{c:%H}/{sub}gfs.t{c:%H}z.pgrb2.0p25.f{f:03d}"
GFS_START = datetime(2021, 1, 1)
LAPSE = 0.0065  # K per m
G = 9.80665
MARGIN = 0.5  # degrees around the scene, so bilinear has neighbours

fetch_only = "--fetch-only" in sys.argv
force_gfs = "--gfs" in sys.argv
argv = [a for a in sys.argv if not a.startswith("--")]
name, ref_path, out_path = argv[1:4]
t_c = float(argv[4]) if len(argv) > 4 else 0.0
model_cache = os.path.splitext(out_path)[0] + "_model.npz"
LAKE_COVER = 0.05  # ERA5 lake ice depth is also set on cells with almost no lake
SNOW_MIN = 0.01  # m; ERA5 leaves traces (~1 mm) on the flooded Red River valley
started = time.perf_counter()

with rasterio.open(ref_path) as ref:
    rows, cols = ref.height, ref.width
    crs_wkt, bounds = ref.crs.to_wkt(), tuple(ref.bounds)
    w, s, e, n = transform_bounds(ref.crs, "EPSG:4326", *ref.bounds, densify_pts=21)
    profile = ref.profile.copy()
    ref_invalid = ref.read_masks(1) == 0

when = datetime.strptime(re.search(r"_(\d{8}T\d{6})_", name).group(1), "%Y%m%dT%H%M%S")
hour = (when + timedelta(minutes=30)).replace(minute=0, second=0)  # nearest full hour


def fetch_era5():
    """Sea-level 2 m temperature, snow depth (m) and ice (m or cover) from ERA5 at
    the hour nearest the S1 time, and their grid."""
    ds = xr.open_zarr(ERA5, chunks=None, storage_options=dict(token="anon"))
    # the time axis runs to 2050; ERA5T (preliminary, ~6 days behind) ends at this date
    stop = ds.attrs.get("valid_time_stop_era5t", ds.attrs["valid_time_stop"])
    if hour >= datetime.fromisoformat(stop) + timedelta(days=1):
        return None
    box = dict(latitude=slice(n + MARGIN, s - MARGIN))
    # ERA5 longitudes run 0..360 in 0.25 steps; contiguous slices (two across 0/360)
    lo, hi = (w - MARGIN) % 360, (e + MARGIN) % 360
    i0, i1 = int(np.floor(lo / 0.25)), int(np.ceil(hi / 0.25)) + 1
    parts = [slice(i0, i1)] if i0 < i1 else [slice(i0, 1440), slice(0, i1)]
    def get(var):
        return xr.concat([ds[var].sel(time=hour).sel(**box).isel(longitude=p)
                          for p in parts], "longitude")
    t2m = get("2m_temperature")
    arr = t2m.values
    if np.isnan(arr).all():
        return None  # not released yet
    print(f"ERA5 2 m temperature at {hour:%Y-%m-%d %H}:00, "
          f"{t2m.sizes['latitude']}x{t2m.sizes['longitude']} cells")
    # static field stored on the time axis: one hour is enough
    oro = get("geopotential_at_surface").values / G
    swe = np.nan_to_num(get("snow_depth").values)  # m of water equivalent
    snow = np.where(swe > 0, swe * 1000 / get("snow_density").values, 0)
    lake_ice = np.where(get("lake_cover").values > LAKE_COVER,
                        np.nan_to_num(get("lake_ice_depth").values), 0)
    ice = np.maximum(lake_ice, np.nan_to_num(get("sea_ice_cover").values))
    lons = (t2m.longitude.values + 180) % 360 - 180
    order = np.argsort(lons)  # put a dateline-crossing box in increasing longitude
    t0 = arr - 273.15 + LAPSE * oro
    return t0[:, order], snow[:, order], ice[:, order], t2m.latitude.values, lons[order]


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


def fetch_gfs(lats, lons):
    """Sea-level 2 m temperature, snow depth (m) and ice thickness (m) from GFS at the
    hour nearest the S1 time."""
    crop = np.ix_(np.rint((90 - lats) * 4).astype(int),
                  np.rint((lons + 180) * 4).astype(int) % 1440)
    cycle = datetime(hour.year, hour.month, hour.day, hour.hour // 6 * 6)
    # hourly output: the latest cycle before the hour, then the one before it
    for cyc in (cycle, cycle - timedelta(hours=6)):
        f = int((hour - cyc) / timedelta(hours=1))
        try:
            tmp = gfs_field(cyc, f, ":TMP:2 m above ground:")[crop]
            oro = gfs_field(cyc, f, ":HGT:surface:")[crop]
            snow = gfs_field(cyc, f, ":SNOD:surface:")[crop]
            ice = gfs_field(cyc, f, ":ICETK:surface:")[crop]
        except (urllib.error.URLError, KeyError):
            continue
        print(f"GFS 2 m temperature, snow depth, ice thickness: "
              f"{cyc:%Y-%m-%d %H}z f{f:03d} = {hour:%Y-%m-%d %H}:00")
        valid = lambda a: np.where(np.isfinite(a) & (a < 1e3), a, 0)  # missing on land
        # GDAL's GRIB driver already gives deg C
        return tmp + LAPSE * oro, valid(snow), valid(ice)
    raise SystemExit(f"no GFS fields for {hour:%Y-%m-%d %H}:00")


ERA5T_LAG = timedelta(days=7)  # ERA5T is about 6 days behind real time
cached = str(np.load(model_cache)["sources"]) if os.path.exists(model_cache) else None
# fetch when there is no cache, or the cache is real-time GFS (upgrade once ERA5 is out)
if cached is None or (cached == "GFS" and not force_gfs):
    era5 = None if force_gfs else fetch_era5()
    if era5 is not None:
        if cached == "GFS":
            print("ERA5 is out now: replacing the real-time GFS fields")
        (t0, snow, ice, lats, lons), sources = era5, "ERA5"
    elif cached == "GFS":
        sources = None  # ERA5 not out yet: keep the GFS cache
    elif hour >= GFS_START:
        lats, lons = box_grid()
        (t0, snow, ice), sources = fetch_gfs(lats, lons), "GFS"
    else:
        raise SystemExit(f"ERA5 has no fields for {hour:%Y-%m-%d %H}:00")
    if sources:
        np.savez(model_cache, t0=t0, snow=snow, ice=ice, lats=lats, lons=lons,
                 sources=sources, hour=str(hour))
if fetch_only:
    print(f"-> {model_cache} ({time.perf_counter() - started:.0f} s)")
    sys.exit()
cache = np.load(model_cache)
t0, snow, ice, lats, lons = (cache[k] for k in ("t0", "snow", "ice", "lats", "lons"))
print(f"model fields over the box ({cache['sources']}, {cache['hour']}): 2 m temperature at "
      f"sea level {np.nanmin(t0):.1f} .. {np.nanmax(t0):.1f} C; cells with snow >= "
      f"{SNOW_MIN} m {(snow >= SNOW_MIN).sum()}, with ice {(ice > 0).sum()} of {snow.size}")


def to_grid(arr, resample="bilinear"):
    """From the 0.25 deg cell centres onto the reference grid."""
    src = gdal.GetDriverByName("MEM").Create("", arr.shape[1], arr.shape[0], 1, gdal.GDT_Float32)
    src.SetGeoTransform((lons[0] - 0.125, 0.25, 0, lats[0] + 0.125, 0, -0.25))
    src.SetProjection("EPSG:4326")
    src.GetRasterBand(1).WriteArray(arr.astype(np.float32))
    out = gdal.Warp("", src, format="MEM", dstSRS=crs_wkt, outputBounds=bounds,
                    width=cols, height=rows, resampleAlg=resample)
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

t_dem = to_grid(t0) - LAPSE * dem
cold = t_dem <= t_c
frozen_cells = (snow >= SNOW_MIN) | maximum_filter(ice > 0, size=3)
frozen = to_grid(frozen_cells, "near") > 0
possible = cold | frozen
out = np.where(ref_invalid, 255, possible).astype(np.uint8)
profile.update(driver="GTiff", count=1, dtype="uint8", nodata=255, compress="deflate",
               tiled=True, blockxsize=256, blockysize=256)
research = str(cache["sources"]) == "ERA5"
with rasterio.open(out_path, "w", **profile) as dst:
    dst.write(out, 1)
    dst.update_tags(t_c=t_c, model_time=str(cache["hour"]), source=str(cache["sources"]),
                    product="research" if research else "realtime", lapse_k_per_m=LAPSE)

valid = ~ref_invalid
print(f"downscaled 2 m temperature in scene: {np.nanmin(t_dem[valid]):.1f} .. "
      f"{np.nanmax(t_dem[valid]):.1f} C; of {valid.sum()} px: <= {t_c} C "
      f"{100 * cold[valid].mean():.2f}%, model snow / ice {100 * frozen[valid].mean():.2f}%, "
      f"snow or ice possible {100 * possible[valid].mean():.2f}%")
print(f"-> {out_path} ({time.perf_counter() - started:.0f} s)")
print("PRODUCT research" if research else f"PRODUCT realtime {hour + ERA5T_LAG:%Y-%m-%d}")
print("SNOW_POSSIBLE" if possible[valid].any() else "SNOW_IMPOSSIBLE")
