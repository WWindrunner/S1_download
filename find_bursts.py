"""Find SLC bursts of an S1 GRD scene plus the two previous acquisitions of each burst.

  python find_bursts.py <S1_GRD_NAME> <out.json> [lookback_days=40] [candidates.tif]

With candidates.tif (EPSG:4326, 1 = candidate pixel), only bursts whose footprint
contains at least one candidate pixel are kept.
"""
import json
import sys
from datetime import timedelta

import asf_search as asf
from dateutil.parser import isoparse
from shapely.geometry import shape

name, out = sys.argv[1:3]
lookback = int(sys.argv[3]) if len(sys.argv) > 3 else 40
cand_path = sys.argv[4] if len(sys.argv) > 4 else None
grd = [g for g in asf.granule_search([name])
       if g.properties["processingLevel"].startswith("GRD")][0]
gp = grd.properties
t_start, t_stop = isoparse(gp["startTime"]), isoparse(gp["stopTime"])
print("GRD", name, "path", gp["pathNumber"], gp["flightDirection"], t_start, t_stop)

pad = timedelta(seconds=10)
co = asf.search(dataset=asf.DATASET.SLC_BURST, relativeOrbit=gp["pathNumber"],
                polarization="VV", start=t_start - pad, end=t_stop + pad,
                intersectsWith=shape(grd.geometry).wkt)
print(f"{len(co)} co-event VV bursts")
if cand_path:
    import numpy as np
    import rasterio
    from rasterio.features import geometry_mask
    with rasterio.open(cand_path) as src:
        cand = src.read(1) == 1
        transform = src.transform
    co = [b for b in co if (cand & geometry_mask([b.geometry], cand.shape, transform,
                                                 invert=True)).any()]
    print(f"{len(co)} bursts contain candidate pixels")

pairs = []
for b in sorted(co, key=lambda b: b.properties["burst"]["fullBurstID"]):
    bid = b.properties["burst"]["fullBurstID"]
    t = isoparse(b.properties["startTime"])
    stack = asf.search(dataset=asf.DATASET.SLC_BURST, fullBurstID=bid, polarization="VV",
                       start=t - timedelta(days=lookback), end=t - timedelta(hours=1))
    stack = sorted(stack, key=lambda s: s.properties["startTime"], reverse=True)[:2]
    rec = dict(burst=bid, co=b.properties["sceneName"], co_date=f"{t:%Y-%m-%d}")
    for k, s in enumerate(stack, 1):
        rec[f"pre{k}"] = s.properties["sceneName"]
        rec[f"pre{k}_date"] = s.properties["startTime"][:10]
    pairs.append(rec)
    print(bid, rec["co_date"], "<-", rec.get("pre1_date"), rec.get("pre1", "")[-10:-6],
          "<-", rec.get("pre2_date"), rec.get("pre2", "")[-10:-6])
json.dump(pairs, open(out, "w"), indent=1)
