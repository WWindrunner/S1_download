"""Historical flood-warning Sentinel-1 search with a year/month cache.

Each UTC day is searched once and cached under ROOT/cache/YYYY/MM/:
  YYYY-MM-DD.json            candidate table (scene names with overlap metrics)
  YYYY-MM-DD_<source>.zip    the downloaded GloFAS layer, retained for rebuilds

Warning source per day:
  >= 2019-11-02  GloFAS RapidFloodMapping (warning areas), as in the daily workflow
  <  2019-11-02  GloFAS ReportingPoints (alert points; first available 2015-12-17)
A day with no RapidFloodMapping file falls back to ReportingPoints.

The cache keeps every candidate above loose floors (regions > 1000 km2, any
warning-pixel overlap; points with ThresGroup >= 2 within 25 km), so the final
selection thresholds can be changed later with `merge` without searching again.
`merge` writes ROOT/results/YYYY.json: a sorted JSON list of unique dual-polarization
S1 names (1SDV/1SDH; single-pol 1SSV/1SSH scenes are excluded).

Usage:
  python S1_archive_search.py search --start 2016-01-01 --end 2016-01-31 --root DIR
  python S1_archive_search.py merge --root DIR --years 2016 2025
"""
import argparse
import datetime
import io
import json
import math
import os
import sys
import tempfile
import zipfile
from time import perf_counter

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.errors import WindowError
from rasterio.features import geometry_window, rasterize
from rasterio.windows import Window
from requests.exceptions import HTTPError
from shapely.geometry import box, mapping
from shapely.ops import nearest_points

import Sentinel_1_ESA_search_download_process_chain_v4 as chain

CACHE_VERSION = 1
GLOFAS_API = "https://european-flood.emergency.copernicus.eu/api/fms/download/glofas"
CATALOGUE = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"
RFM_START = datetime.date(2019, 11, 2)

# Cache floors: the loosest selections that can still be applied later.
WINDOW_SIZE = 10
CACHE_MIN_REGION_KM2 = 1000
CACHE_MIN_POINT_TG = 2
CACHE_MAX_POINT_DIST_KM = 25.0

# Default final selection (chosen from the 2016-2025 sampling test).
DEFAULT_MIN_REGION_KM2 = 5000
DEFAULT_MIN_OVERLAP_KM2 = 50.0
DEFAULT_POINT_TG = 3
DEFAULT_POINT_PROB = 50.0
DEFAULT_POINT_DIST_KM = 10.0

# GloFAS and CDSE throttle parallel clients (HTTP 504/429); retry patiently.
REQUEST_TIMEOUT = 300
REQUEST_ATTEMPTS = 8


class MissingSource(Exception):
    """GloFAS has no file for this layer and date."""


def _get(url):
    return chain._get_search_data(url, timeout=REQUEST_TIMEOUT, attempts=REQUEST_ATTEMPTS)


def day_paths(root, day):
    directory = os.path.join(root, "cache", f"{day:%Y}", f"{day:%m}")
    return directory, os.path.join(directory, f"{day.isoformat()}.json")


def _atomic_write(path, data, mode="w"):
    directory = os.path.dirname(path)
    os.makedirs(directory, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=directory, prefix=".tmp_", suffix=os.path.splitext(path)[1])
    try:
        with os.fdopen(fd, mode, **({} if "b" in mode else {"encoding": "utf-8"})) as stream:
            stream.write(data)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def fetch_layer(root, day, layer, tag):
    """Return the cached GloFAS zip path, downloading it if absent or corrupt."""
    directory, _ = day_paths(root, day)
    path = os.path.join(directory, f"{day.isoformat()}_{tag}.zip")
    if os.path.isfile(path) and zipfile.is_zipfile(path):
        return path
    url = f"{GLOFAS_API}/{layer}/{day:%Y%m%d}T00:00Z"
    try:
        response = _get(url)
    except HTTPError as exc:
        response = exc.response
        if response is not None and response.status_code == 400 and "No files" in response.text:
            raise MissingSource(f"{layer} {day}") from exc
        raise
    if not zipfile.is_zipfile(io.BytesIO(response.content)):
        raise ValueError(f"GloFAS {layer} {day} did not return a zip file")
    _atomic_write(path, response.content, "wb")
    return path


def catalogue_products(bounds, day):
    """All Sentinel-1 IW GRDH products starting on `day` that intersect `bounds`."""
    left, bottom, right, top = bounds
    polygon = f"{left} {top}, {left} {bottom}, {right} {bottom}, {right} {top}, {left} {top}"
    start = f"{day.isoformat()}T00:00:00.000Z"
    end = f"{(day + datetime.timedelta(days=1)).isoformat()}T00:00:00.000Z"
    url = (CATALOGUE
           + f"$filter=OData.CSC.Intersects(area=geography'SRID=4326;POLYGON(({polygon}))') "
           "and Collection/Name eq 'SENTINEL-1' "
           "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' "
           "and att/OData.CSC.StringAttribute/Value eq 'IW_GRDH_1S') "
           f"and ContentDate/Start ge {start} and ContentDate/Start lt {end}")
    products = {}
    while url:
        payload = _get(url).json()
        for product in payload["value"]:
            products.setdefault(product["Id"], product)
        url = payload.get("@odata.nextLink") or payload.get("@OData.nextLink")
    return products


def overlap_km2(product, feature):
    """Full count (no early stop) of this region's original warning pixels in the footprint."""
    footprint = chain._product_footprint(product)
    extent = box(*footprint.bounds).intersection(
        box(feature["minx"], feature["miny"], feature["maxx"], feature["maxy"]))
    if extent.is_empty or extent.area == 0:
        return 0.0
    with rasterio.open(feature["warning_raster"]) as src:
        try:
            window = geometry_window(src, [mapping(extent)])
        except WindowError:
            return 0.0
        transform = src.transform
        pixel_km2 = abs(transform.a * transform.e - transform.b * transform.d) * 110 * 110
        count = 0
        row_stop = int(window.row_off + window.height)
        col_stop = int(window.col_off + window.width)
        for row in range(int(window.row_off), row_stop, 512):
            for col in range(int(window.col_off), col_stop, 512):
                tile = Window(col, row, min(512, col_stop - col),
                                               min(512, row_stop - row))
                selected = src.read(1, window=tile) == int(feature["warning_region"])
                if not selected.any():
                    continue
                cover = rasterize([(mapping(footprint), 1)], out_shape=selected.shape,
                                  transform=src.window_transform(tile), dtype="uint8",
                                  all_touched=False)
                count += int(np.count_nonzero(selected & (cover == 1)))
    return count * pixel_km2


def search_areas(zip_path, day, workdir):
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(workdir)
    shapefile = os.path.join(workdir, f"FloodMaskMerged{day:%Y%m%d}00.shp")
    regions = chain.simplify_flood_warning_shp_from_ESA(shapefile, WINDOW_SIZE, CACHE_MIN_REGION_KM2)
    candidates = []
    for group in chain._group_sentinel_warning_regions(regions):
        products = catalogue_products(group["bounds"], day)
        for feature in group["features"]:
            for product in products.values():
                area = overlap_km2(product, feature)
                if area > 0:
                    candidates.append({
                        "name": product["Name"].removesuffix(".SAFE"),
                        "region": int(feature["warning_region"]),
                        "region_km2": round(float(feature["area_km2"]), 3),
                        "overlap_km2": round(area, 3)})
    return candidates, {"regions": len(regions)}


def _haversine_km(p, q):
    lat1, lat2 = math.radians(p.y), math.radians(q.y)
    a = (math.sin((lat2 - lat1) / 2) ** 2
         + math.cos(lat1) * math.cos(lat2) * math.sin(math.radians(q.x - p.x) / 2) ** 2)
    return 6371.0 * 2 * math.asin(math.sqrt(min(1.0, a)))


def search_points(zip_path, day, workdir):
    with zipfile.ZipFile(zip_path) as archive:
        archive.extractall(workdir)
    shapefiles = [f for f in os.listdir(workdir) if f.lower().endswith(".shp")]
    if len(shapefiles) != 1:
        raise ValueError(f"Expected one reporting-point shapefile, found {shapefiles}")
    points = gpd.read_file(os.path.join(workdir, shapefiles[0]))
    points = points.set_crs(4326) if points.crs is None else points.to_crs(4326)
    points["ThresGroup"] = pd.to_numeric(points["ThresGroup"], errors="coerce")
    points["Prob"] = pd.to_numeric(points["Prob"], errors="coerce")
    alerts = points[points["ThresGroup"] >= CACHE_MIN_POINT_TG].reset_index(drop=True)
    pad = CACHE_MAX_POINT_DIST_KM / 111.0
    features = pd.DataFrame({
        "minx": [max(-180.0, g.x - pad / max(math.cos(math.radians(g.y)), 0.2)) for g in alerts.geometry],
        "maxx": [min(180.0, g.x + pad / max(math.cos(math.radians(g.y)), 0.2)) for g in alerts.geometry],
        "miny": [max(-90.0, g.y - pad) for g in alerts.geometry],
        "maxy": [min(90.0, g.y + pad) for g in alerts.geometry],
        "index": range(len(alerts))})
    candidates = []
    for group in chain._group_sentinel_warning_regions(features):
        products = catalogue_products(group["bounds"], day)
        footprints = {pid: chain._product_footprint(p) for pid, p in products.items()}
        for feature in group["features"]:
            alert = alerts.iloc[int(feature["index"])]
            point = alert.geometry
            for pid, footprint in footprints.items():
                distance = 0.0 if footprint.contains(point) else _haversine_km(
                    point, nearest_points(point, footprint)[1])
                if distance <= CACHE_MAX_POINT_DIST_KM:
                    candidates.append({
                        "name": products[pid]["Name"].removesuffix(".SAFE"),
                        "point": str(alert.get("PointID")),
                        "tg": int(alert["ThresGroup"]),
                        "prob": float(alert["Prob"]),
                        "dist_km": round(distance, 3)})
    return candidates, {"points": len(points), "alert_points": len(alerts)}


def _cache_valid(saved, day):
    return (isinstance(saved, dict) and saved.get("version") == CACHE_VERSION
            and saved.get("date") == day.isoformat()
            and saved.get("source") in ("areas", "points", "none")
            and isinstance(saved.get("candidates"), list))


def read_day(root, day):
    _, path = day_paths(root, day)
    try:
        with open(path, encoding="utf-8") as stream:
            saved = json.load(stream)
    except (OSError, ValueError):
        return None
    return saved if _cache_valid(saved, day) else None


def search_day(root, day, refresh=False):
    """Search one day into the cache; return the cached record."""
    if not refresh:
        saved = read_day(root, day)
        if saved is not None:
            return saved
    started = perf_counter()
    plans = ([("areas", "RapidFloodMapping", "rfm", search_areas)] if day >= RFM_START else [])
    plans.append(("points", "ReportingPoints", "points", search_points))
    record = None
    for source, layer, tag, search in plans:
        try:
            zip_path = fetch_layer(root, day, layer, tag)
        except MissingSource:
            print(f"{day}: no GloFAS {layer} file.", flush=True)
            continue
        with tempfile.TemporaryDirectory(prefix="s1_archive_") as workdir:
            candidates, info = search(zip_path, day, workdir)
        record = {"source": source, "candidates": candidates, **info}
        break
    if record is None:
        record = {"source": "none", "candidates": []}
    record = {"version": CACHE_VERSION, "date": day.isoformat(), **record,
              "params": {"window_size": WINDOW_SIZE, "min_region_km2": CACHE_MIN_REGION_KM2,
                         "min_point_tg": CACHE_MIN_POINT_TG,
                         "max_point_dist_km": CACHE_MAX_POINT_DIST_KM},
              "seconds": round(perf_counter() - started, 1)}
    _, path = day_paths(root, day)
    _atomic_write(path, json.dumps(record, indent=1) + "\n")
    names = {c["name"] for c in record["candidates"]}
    print(f"{day}: source={record['source']}, {len(names)} candidate scenes, "
          f"{record['seconds']:.0f}s.", flush=True)
    return record


def search_range(root, start, end, refresh=False):
    failed = []
    day = start
    while day <= end:
        try:
            search_day(root, day, refresh)
        except Exception as exc:
            failed.append(day.isoformat())
            print(f"{day}: search failed: {exc!r}", file=sys.stderr, flush=True)
        day += datetime.timedelta(days=1)
    print(f"Searched {start} to {end}; failed days: {', '.join(failed) or 'none'}.", flush=True)
    return failed


def select_names(record, min_region_km2=DEFAULT_MIN_REGION_KM2,
                 min_overlap_km2=DEFAULT_MIN_OVERLAP_KM2, point_tg=DEFAULT_POINT_TG,
                 point_prob=DEFAULT_POINT_PROB, point_dist_km=DEFAULT_POINT_DIST_KM):
    if record["source"] == "areas":
        keep = (c for c in record["candidates"]
                if c["region_km2"] > min_region_km2 and c["overlap_km2"] >= min_overlap_km2)
    elif record["source"] == "points":
        keep = (c for c in record["candidates"]
                if c["tg"] >= point_tg and c["prob"] >= point_prob and c["dist_km"] <= point_dist_km)
    else:
        keep = ()
    # Dual-polarization only (VV+VH or HH+HV); single-pol 1SSV/1SSH is out of scope.
    return {c["name"] for c in keep if c["name"].startswith("S1") and "_1SS" not in c["name"]}


def merge_year(root, year, allow_incomplete=False, **thresholds):
    """Write ROOT/results/YEAR.json from the cache; return (count, missing_days)."""
    names, missing = set(), []
    day = datetime.date(year, 1, 1)
    while day.year == year:
        record = read_day(root, day)
        if record is None:
            missing.append(day.isoformat())
        else:
            names |= select_names(record, **thresholds)
        day += datetime.timedelta(days=1)
    if missing and not allow_incomplete:
        return None, missing
    _atomic_write(os.path.join(root, "results", f"{year}.json"),
                  json.dumps(sorted(names), indent=2) + "\n")
    return len(names), missing


def _date(value):
    return datetime.date.fromisoformat(value)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    search = sub.add_parser("search", help="Search days into the cache (resumable).")
    search.add_argument("--start", type=_date, required=True, help="First UTC day, inclusive.")
    search.add_argument("--end", type=_date, required=True, help="Last UTC day, inclusive.")
    search.add_argument("--root", required=True)
    search.add_argument("--refresh", action="store_true", help="Search cached days again.")
    merge = sub.add_parser("merge", help="Write results/YYYY.json from the cache.")
    merge.add_argument("--root", required=True)
    merge.add_argument("--years", type=int, nargs=2, metavar=("FIRST", "LAST"), required=True)
    merge.add_argument("--allow-incomplete", action="store_true",
                       help="Write years that have uncached days (they are still reported).")
    merge.add_argument("--min-region-km2", type=float, default=DEFAULT_MIN_REGION_KM2)
    merge.add_argument("--min-overlap-km2", type=float, default=DEFAULT_MIN_OVERLAP_KM2)
    merge.add_argument("--point-tg", type=int, default=DEFAULT_POINT_TG)
    merge.add_argument("--point-prob", type=float, default=DEFAULT_POINT_PROB)
    merge.add_argument("--point-dist-km", type=float, default=DEFAULT_POINT_DIST_KM)
    args = parser.parse_args(argv)

    if args.command == "search":
        if args.start > args.end:
            parser.error("--start must be on or before --end")
        return 1 if search_range(args.root, args.start, args.end, args.refresh) else 0

    if (args.min_region_km2 < CACHE_MIN_REGION_KM2 or args.point_tg < CACHE_MIN_POINT_TG
            or args.point_dist_km > CACHE_MAX_POINT_DIST_KM):
        parser.error(f"Thresholds must be within the cache floors: region > {CACHE_MIN_REGION_KM2} km2, "
                     f"ThresGroup >= {CACHE_MIN_POINT_TG}, distance <= {CACHE_MAX_POINT_DIST_KM} km")
    thresholds = dict(min_region_km2=args.min_region_km2, min_overlap_km2=args.min_overlap_km2,
                      point_tg=args.point_tg, point_prob=args.point_prob,
                      point_dist_km=args.point_dist_km)
    incomplete = False
    for year in range(args.years[0], args.years[1] + 1):
        count, missing = merge_year(args.root, year, args.allow_incomplete, **thresholds)
        if missing:
            incomplete = True
            shown = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
            print(f"{year}: {len(missing)} uncached days ({shown}).", flush=True)
        if count is None:
            print(f"{year}: not written (incomplete).", flush=True)
        else:
            print(f"{year}: {count} scenes -> results/{year}.json", flush=True)
    return 1 if incomplete else 0


if __name__ == "__main__":
    sys.exit(main())
