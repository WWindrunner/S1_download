import os
from eodag import EODataAccessGateway,setup_logging,SearchResult
from eodag.utils import ProgressCallback
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from concurrent.futures import ThreadPoolExecutor
from pyroSAR import identify
from datetime import timedelta
import glob
from datetime import datetime, timedelta
import math
import time
import re
import sys
from collections import defaultdict
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
import shutil
import traceback
import urllib.request

import planetary_computer


if len(sys.argv) != 3:
    print(
        "Usage: python Snow_detect.py "
        "<S1_PRODUCT_NAME> <OUTPUT_FOLDER>"
    )
    sys.exit(1)

S1name = sys.argv[1]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))
safe_dir = os.path.join(folder, S1name)

print(f"Path: {safe_dir}")

if not os.path.isdir(safe_dir):
    raise FileNotFoundError(f"S1 product folder not found: {safe_dir}")

tifs = glob.glob(os.path.join(safe_dir, "*VV.tif"))
s1_file = tifs[0]
workspace = os.path.join(safe_dir, "snow_temp")
os.makedirs(workspace, exist_ok=True)
snow_out = os.path.join(safe_dir, f"{S1name}_ice.tif")


def cleanup_intermediate_files():
    if os.path.isdir(workspace):
        shutil.rmtree(workspace)

    cleanup_patterns = [
        "*.zip",
        "*incidenceAngleFromEllipsoid.tif",
        "*localIncidenceAngle.tif",
        "*_manifest.safe",
        "*_proc.xml",
        "*_gamma0-elp.tif",
    ]

    for pattern in cleanup_patterns:
        for file_path in glob.glob(os.path.join(safe_dir, pattern)):
            os.remove(file_path)
            print(f"Removed intermediate file: {file_path}")

    dem_tiles = os.path.join(safe_dir, "dem_tiles")
    if os.path.isdir(dem_tiles):
        shutil.rmtree(dem_tiles)
        print(f"Removed intermediate directory: {dem_tiles}")


def download_scl_asset(scl_href, scl_file):
    """Download one Planetary Computer asset without EODAG's auth headers."""
    signed_href = planetary_computer.sign_url(scl_href)
    os.makedirs(os.path.dirname(scl_file), exist_ok=True)
    partial_file = f"{scl_file}.part"
    try:
        with urllib.request.urlopen(signed_href, timeout=300) as response:
            with open(partial_file, "wb") as output:
                shutil.copyfileobj(response, output)
        os.replace(partial_file, scl_file)
    finally:
        if os.path.exists(partial_file):
            os.remove(partial_file)

    return scl_file


def create_nodata_ice_raster():
    s1_ds = gdal.Open(s1_file)
    driver = gdal.GetDriverByName("GTiff")
    out_ds = driver.Create(
        snow_out,
        s1_ds.RasterXSize,
        s1_ds.RasterYSize,
        1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"]
    )
    out_ds.SetGeoTransform(s1_ds.GetGeoTransform())
    out_ds.SetProjection(s1_ds.GetProjection())
    band = out_ds.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    band.Fill(-9999)
    band.FlushCache()
    out_ds = None
    s1_ds = None


zip_path = glob.glob(os.path.join(safe_dir, "*.zip"))[0]
scene = identify(zip_path)
os.environ["EODAG__PLANETARY_COMPUTER__DOWNLOAD__OUTPUT_DIR"] = os.path.abspath(workspace)
dag = EODataAccessGateway()
setup_logging(0)
end_dt = datetime.fromisoformat(scene.start)
start_dt = end_dt - timedelta(days=15)
end = end_dt.strftime("%Y-%m-%d")
start = start_dt.strftime("%Y-%m-%d")
bbox = scene.bbox()
geom = {
    "lonmin": bbox.extent['xmin'],
    "latmin": bbox.extent['ymin'],
    "lonmax": bbox.extent['xmax'],
    "latmax": bbox.extent['ymax'],
}
search_results = dag.search_all(
    provider="planetary_computer",
    collection="S2_MSI_L2A",
    start=start,
    end=end,
    geom=geom
)
n = len(search_results)
print(n)
if not search_results:
    print(
        "WARNING: No Sentinel-2 products found for the search period and "
        "area. Creating a nodata ice raster."
    )
    create_nodata_ice_raster()
    cleanup_intermediate_files()
    print("\nFinished with no Sentinel-2 observations:")
    print(snow_out)
    sys.exit(0)

for product in search_results:
    product_name = product.properties.get(
        "title",
        product.properties.get("id", "unknown_product")
    )
    scl_file = os.path.join(
        workspace,
        product_name,
        "SCL_20m.tif"
    )
    print(scl_file)
    if os.path.exists(scl_file):
        print(f"Already exists, skip: {product_name}")
        continue

    scl_asset = product.assets.get("SCL_20m")
    if scl_asset is None:
        print(f"Download skipped for {product_name}: SCL_20m asset is missing")
        print(f"Available assets: {list(product.assets.keys())}")
        continue

    scl_href = scl_asset.get("href")
    if not isinstance(scl_href, str) or not scl_href:
        print(
            f"Download skipped for {product_name}: "
            f"SCL_20m has an invalid href: {scl_href!r}"
        )
        continue

    product_start_time = time.perf_counter()
    try:
        downloaded_path = download_scl_asset(scl_href, scl_file)
        elapsed = time.perf_counter() - product_start_time
        print(f"Downloaded: {downloaded_path}")
        print(f"Download time: {elapsed:.2f} seconds")
    except Exception as e:
        elapsed = time.perf_counter() - product_start_time
        print(
            f"Download failed for {product_name} after {elapsed:.2f} seconds: "
            f"{type(e).__name__}: {e}"
        )
        traceback.print_exc()
        time.sleep(10)
        continue

mosaic_dir = os.path.join(workspace, "Daily_Mosaic_S1grid")
os.makedirs(mosaic_dir, exist_ok=True)
s1_ds = gdal.Open(s1_file)
s1_proj = s1_ds.GetProjection()
s1_gt = s1_ds.GetGeoTransform()
s1_cols = s1_ds.RasterXSize
s1_rows = s1_ds.RasterYSize
xmin = s1_gt[0]
ymax = s1_gt[3]
pixel_x = s1_gt[1]
pixel_y = s1_gt[5]
xmax = xmin + s1_cols * pixel_x
ymin = ymax + s1_rows * pixel_y
s1_bounds = [xmin, ymin, xmax, ymax]
s1_ds=None
date_dict = defaultdict(list)
products = sorted(glob.glob(os.path.join(workspace, "S2*_MSIL2A_*")))

for product in products:
    name=os.path.basename(product)
    m=re.search(r"MSIL2A_(\d{8})T", name)
    if m is None:
        continue
    date=m.group(1)
    tifs=glob.glob(os.path.join(product, "**", "*.tif"), recursive=True)
    date_dict[date].extend(tifs)

for date in sorted(date_dict.keys()):
    tif_list=date_dict[date]
    outfile=os.path.join(mosaic_dir, f"S2_{date}_S1grid.tif")
    if os.path.exists(outfile):
        print("exist skip")
        continue

    gdal.Warp(
        outfile,
        tif_list,
        # Sentinel-1 CRS
        dstSRS=s1_proj,
        # Sentinel-1 extent
        outputBounds=s1_bounds,
        # Sentinel-1 resolution
        xRes=abs(pixel_x),
        yRes=abs(pixel_y),
        # 分类数据必须nearest
        resampleAlg=
        gdal.GRA_NearestNeighbour,
        dstNodata=0,
        multithread=True,
        creationOptions=[
            "COMPRESS=LZW",
            "TILED=YES",
            "BIGTIFF=IF_SAFER"
        ]
    )

daily_files=sorted(glob.glob(os.path.join(mosaic_dir,"*.tif")))
snow_count=np.zeros((s1_rows,s1_cols),dtype=np.uint16)
valid_count=np.zeros((s1_rows,s1_cols),dtype=np.uint16)

for i,tif in enumerate(daily_files):
    ds=gdal.Open(tif)
    arr=(ds.GetRasterBand(1).ReadAsArray())
    nodata=(ds.GetRasterBand(1).GetNoDataValue())

    if nodata is not None:
        valid=arr!=nodata
    else:
        valid=arr!=0

    snow=((arr==11)&valid)
    snow_count += snow.astype(np.uint16)
    valid_count += valid.astype(np.uint16)
    ds=None

snow_frequency=np.zeros((s1_rows,s1_cols),dtype=np.float32)
idx=valid_count>0

snow_frequency[idx]=(snow_count[idx]/valid_count[idx])

driver=gdal.GetDriverByName("GTiff")

out_ds=driver.Create(
    snow_out,
    s1_cols,
    s1_rows,
    1,
    gdal.GDT_Float32,
    options=[
        "COMPRESS=LZW",
        "TILED=YES",
        "BIGTIFF=IF_SAFER"
    ]
)
out_ds.SetGeoTransform(s1_gt)
out_ds.SetProjection(s1_proj)
band=out_ds.GetRasterBand(1)
band.WriteArray(snow_frequency)
band.SetNoDataValue(-9999)
band.FlushCache()
out_ds=None

cleanup_intermediate_files()


print("\nFinished:")
print(snow_out)





