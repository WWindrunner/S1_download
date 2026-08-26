"""Generate snow and persistent-cloud masks from Sentinel-2 SCL data."""

import argparse
from collections import defaultdict
from datetime import datetime, timedelta
import os
import re
import shutil
import sys
import time
import urllib.request

import numpy as np
import planetary_computer
import rasterio
from eodag import EODataAccessGateway, setup_logging
from osgeo import gdal
from rasterio.warp import transform_bounds

gdal.UseExceptions()

MASK_NODATA = 255
CLOUD_CLASSES = (3, 8, 9, 10)
S1_TIME_PATTERN = re.compile(r"_([0-9]{8}T[0-9]{6})_")
S2_DATE_PATTERN = re.compile(r"MSIL2A?_([0-9]{8})T")


def parse_s1_time(product_name):
    match = S1_TIME_PATTERN.search(product_name)
    if match is None:
        raise ValueError(f"Cannot parse acquisition time from S1 name: {product_name}")
    return datetime.strptime(match.group(1), "%Y%m%dT%H%M%S")


def download_scl_asset(href, destination):
    signed_href = planetary_computer.sign_url(href)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    partial = f"{destination}.part"
    try:
        with urllib.request.urlopen(signed_href, timeout=300) as response:
            with open(partial, "wb") as output:
                shutil.copyfileobj(response, output)
        os.replace(partial, destination)
    finally:
        if os.path.exists(partial):
            os.remove(partial)


def write_mask(path, data, profile, description):
    output_profile = profile.copy()
    output_profile.update(
        driver="GTiff", count=1, dtype="uint8", nodata=MASK_NODATA,
        compress="deflate", tiled=True, blockxsize=256, blockysize=256,
        BIGTIFF="IF_SAFER",
    )
    with rasterio.open(path, "w", **output_profile) as dst:
        dst.write(data.astype(np.uint8), 1)
        dst.set_band_description(1, description)


def generate_s2_masks(product_name, reference_path, output_dir):
    started = time.perf_counter()
    reference_path = os.path.abspath(os.path.expanduser(reference_path))
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    if not os.path.isfile(reference_path):
        raise FileNotFoundError(f"Gamma0 reference not found: {reference_path}")
    os.makedirs(output_dir, exist_ok=True)
    workspace = os.path.join(output_dir, "snow_temp")
    os.makedirs(workspace, exist_ok=True)
    snow_out = os.path.join(output_dir, f"{product_name}_ice.tif")
    cloud_out = os.path.join(output_dir, f"{product_name}_cloud.tif")

    with rasterio.open(reference_path) as reference:
        if reference.crs is None:
            raise ValueError(f"Reference raster has no CRS: {reference_path}")
        profile = reference.profile.copy()
        rows, cols = reference.height, reference.width
        bounds_wgs84 = transform_bounds(
            reference.crs, "EPSG:4326", *reference.bounds, densify_pts=21
        )
        reference_projection = reference.crs.to_wkt()
        reference_bounds = tuple(reference.bounds)

    end_dt = parse_s1_time(product_name)
    start_dt = end_dt - timedelta(days=15)
    geom = {
        "lonmin": bounds_wgs84[0], "latmin": bounds_wgs84[1],
        "lonmax": bounds_wgs84[2], "latmax": bounds_wgs84[3],
    }
    os.environ["EODAG__PLANETARY_COMPUTER__DOWNLOAD__OUTPUT_DIR"] = workspace
    setup_logging(0)
    results = EODataAccessGateway().search_all(
        provider="planetary_computer", collection="S2_MSI_L2A",
        start=start_dt.strftime("%Y-%m-%d"), end=end_dt.strftime("%Y-%m-%d"),
        geom=geom,
    )
    print(f"Sentinel-2 search completed: {len(results)} products found.")

    files_by_date = defaultdict(list)
    for product in results:
        title = product.properties.get("title", product.properties.get("id", ""))
        date_match = S2_DATE_PATTERN.search(title)
        asset = product.assets.get("SCL_20m")
        href = asset.get("href") if asset is not None else None
        if date_match is None or not isinstance(href, str) or not href:
            print(f"Skipping Sentinel-2 product without usable SCL/date: {title}")
            continue
        scl_path = os.path.join(workspace, title, "SCL_20m.tif")
        if not os.path.isfile(scl_path):
            try:
                download_scl_asset(href, scl_path)
            except Exception as exc:
                print(f"SCL download failed for {title}: {exc}", file=sys.stderr)
                continue
        files_by_date[date_match.group(1)].append(scl_path)

    snow_any = np.zeros((rows, cols), dtype=bool)
    valid_day_count = np.zeros((rows, cols), dtype=np.uint16)
    cloud_day_count = np.zeros((rows, cols), dtype=np.uint16)

    for date, scl_paths in sorted(files_by_date.items()):
        daily_snow = np.zeros((rows, cols), dtype=bool)
        daily_valid = np.zeros((rows, cols), dtype=bool)
        daily_noncloud = np.zeros((rows, cols), dtype=bool)
        for index, scl_path in enumerate(scl_paths):
            aligned_path = os.path.join(workspace, f"aligned_{date}_{index}.tif")
            aligned = gdal.Warp(
                aligned_path, scl_path, dstSRS=reference_projection,
                outputBounds=reference_bounds, width=cols, height=rows,
                resampleAlg=gdal.GRA_NearestNeighbour, srcNodata=0, dstNodata=0,
                multithread=True,
                creationOptions=["TILED=YES", "COMPRESS=DEFLATE"],
            )
            if aligned is None:
                print(f"Failed to align SCL raster: {scl_path}", file=sys.stderr)
                continue
            aligned = None
            with rasterio.open(aligned_path) as src:
                arr = src.read(1)
            os.remove(aligned_path)
            valid = arr != 0
            cloud = valid & np.isin(arr, CLOUD_CLASSES)
            daily_snow |= valid & (arr == 11)
            daily_valid |= valid
            daily_noncloud |= valid & ~cloud

        daily_cloud = daily_valid & ~daily_noncloud
        snow_any |= daily_snow
        valid_day_count += daily_valid.astype(np.uint16)
        cloud_day_count += daily_cloud.astype(np.uint16)

    observed = valid_day_count > 0
    snow_mask = np.full((rows, cols), MASK_NODATA, dtype=np.uint8)
    cloud_mask = np.full((rows, cols), MASK_NODATA, dtype=np.uint8)
    snow_mask[observed] = snow_any[observed].astype(np.uint8)
    cloud_mask[observed] = (
        cloud_day_count[observed] == valid_day_count[observed]
    ).astype(np.uint8)
    write_mask(snow_out, snow_mask, profile, "snow_or_ice_seen_in_15_day_window")
    write_mask(cloud_out, cloud_mask, profile, "cloud_on_all_valid_observation_days")
    shutil.rmtree(workspace)
    print(f"Snow mask: {snow_out}")
    print(f"Cloud mask: {cloud_out}")
    print(f"Sentinel-2 masks completed in {time.perf_counter() - started:.2f} seconds.")
    return snow_out, cloud_out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("s1_product_name")
    parser.add_argument(
        "reference_raster",
        help=(
            "Gamma0 reference raster, or the legacy output root when "
            "OUTPUT_DIR is omitted."
        ),
    )
    parser.add_argument("output_dir", nargs="?")
    args = parser.parse_args()
    legacy_cli = args.output_dir is None
    if legacy_cli:
        output_dir = os.path.join(
            os.path.abspath(os.path.expanduser(args.reference_raster)),
            args.s1_product_name,
        )
        reference_path = os.path.join(output_dir, "Gamma0_VV.tif")
    else:
        output_dir = args.output_dir
        reference_path = args.reference_raster

    generate_s2_masks(
        args.s1_product_name, reference_path, output_dir
    )


if __name__ == "__main__":
    main()
