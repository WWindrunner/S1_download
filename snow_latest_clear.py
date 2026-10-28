"""Snow mask from the latest clear Sentinel-2 SCL observation, with its age.

  python snow_latest_clear.py <S1_NAME> <REFERENCE_RASTER> <OUT_TIF> [days=15]

Window: the S1 acquisition day plus the preceding days-1 days (same-day S2
counts, even when acquired after S1). Per pixel, walk back from the S1 day to
the first day with a clear SCL observation and take snow (SCL 11) from it.

Output values on the reference grid:
    0      latest clear observation is not snow
    1..15  snow; day of that observation counted back from the S1 day
           (1 = same day, 2 = one day earlier, ...). Smaller = more trustable
    200    observed only as cloud / cloud shadow / cirrus: snow undetermined
    255    no data: no S2 observation in the window, or reference is nodata
"""
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timedelta

import numpy as np
import planetary_computer
import pystac_client
import rasterio
from osgeo import gdal
from rasterio.warp import transform_bounds

gdal.UseExceptions()
NO_SNOW, CLOUD, NODATA = 0, 200, 255
CLOUD_CLASSES = (3, 8, 9, 10)  # cloud shadow, cloud medium/high, cirrus
NOT_CLEAR = (0, 1) + CLOUD_CLASSES  # plus nodata and saturated/defective
SNOW = 11

name, ref_path, out_path = sys.argv[1:4]
days = int(sys.argv[4]) if len(sys.argv) > 4 else 15
started = time.perf_counter()

with rasterio.open(ref_path) as ref:
    rows, cols = ref.height, ref.width
    crs_wkt, bounds = ref.crs.to_wkt(), tuple(ref.bounds)
    bbox = transform_bounds(ref.crs, "EPSG:4326", *ref.bounds, densify_pts=21)
    profile = ref.profile.copy()
    ref_values = ref.read(1)
    ref_invalid = ref.read_masks(1) == 0
    if np.issubdtype(ref_values.dtype, np.floating):
        ref_invalid |= ~np.isfinite(ref_values)
    del ref_values

day = datetime.strptime(re.search(r"_(\d{8})T\d{6}_", name).group(1), "%Y%m%d")
start = day - timedelta(days=days - 1)
catalog = pystac_client.Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1",
    modifier=planetary_computer.sign_inplace,
)
items = catalog.search(
    collections=["sentinel-2-l2a"], bbox=bbox,
    datetime=f"{start:%Y-%m-%d}/{day:%Y-%m-%d}T23:59:59Z",
).items()
by_date = defaultdict(list)
for item in items:
    by_date[item.datetime.strftime("%Y-%m-%d")].append(item.assets["SCL"].href)
print(f"{sum(map(len, by_date.values()))} S2 SCL assets over {len(by_date)} dates, "
      f"{start:%Y-%m-%d} .. {day:%Y-%m-%d}")

mask = np.full((rows, cols), NODATA, np.uint8)
seen_cloud = np.zeros((rows, cols), bool)
print("date        value  newly-decided%  snow-of-new%  undecided%")
for date in sorted(by_date, reverse=True):
    value = (day - datetime.strptime(date, "%Y-%m-%d")).days + 1
    # Warp handles sources in different UTM zones; nearest keeps SCL classes.
    ds = gdal.Warp("", ["/vsicurl/" + h for h in by_date[date]], format="MEM",
                   dstSRS=crs_wkt, outputBounds=bounds, width=cols, height=rows,
                   resampleAlg="near", srcNodata=0, dstNodata=0, multithread=True)
    scl = ds.ReadAsArray()
    ds = None
    undecided = mask == NODATA
    new = undecided & ~np.isin(scl, NOT_CLEAR)
    seen_cloud |= undecided & np.isin(scl, CLOUD_CLASSES)
    mask[new] = np.where(scl[new] == SNOW, value, NO_SNOW)
    n_new = new.sum()
    print(f"{date}  {value:5d}  {100 * n_new / mask.size:14.2f}  "
          f"{100 * (mask[new] == value).sum() / max(1, n_new):12.2f}  "
          f"{100 * (mask == NODATA).mean():10.2f}", flush=True)
    del scl, new, undecided

mask[(mask == NODATA) & seen_cloud] = CLOUD
mask[ref_invalid] = NODATA

profile.update(driver="GTiff", count=1, dtype="uint8", nodata=NODATA,
               compress="deflate", tiled=True, blockxsize=256, blockysize=256,
               BIGTIFF="IF_SAFER")
with rasterio.open(out_path, "w", **profile) as dst:
    dst.write(mask, 1)
    dst.set_band_description(
        1, "0 no snow; 1-15 snow, day of latest clear S2 obs counted back from "
           "S1 day (1 = same day); 200 cloud only; 255 nodata")

valid = ~ref_invalid
counts = {v: int(((mask == v) & valid).sum()) for v in np.unique(mask[valid])}
print(f"value counts within reference data ({valid.sum()} px): {counts}")
print(f"-> {out_path} ({time.perf_counter() - started:.0f} s)")
