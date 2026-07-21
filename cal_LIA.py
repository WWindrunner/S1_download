from pyroSAR import identify
from pystac_client import Client
import planetary_computer as pc
import requests
import os
from osgeo import gdal
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
import shutil
# -*- coding: utf-8 -*-
import rasterio
import zipfile
import xml.etree.ElementTree as ET
import sys


# -------------------------
# Parse input arguments
# -------------------------
if len(sys.argv) != 3:
    print(
        "Usage: python cal_LIA.py "
        "<S1_PRODUCT_NAME> <OUTPUT_FOLDER>"
    )
    sys.exit(1)

S1name = sys.argv[1]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))
path = os.path.join(folder, S1name)

print(f"Path: {path}")

if not os.path.isdir(path):
    raise FileNotFoundError(f"S1 product folder not found: {path}")

#path = "/tank/data/SFS/xinyis/zhao89/data/Sentinel_1/specific_images/S1A_IW_GRDH_1SDV_20231104T091539_20231104T091604_051065_06285B_6CA3/"

zip_path = glob.glob(os.path.join(path, "*.zip"))[0]
print(zip_path)
scene = identify(zip_path)

S1_angle_path = glob.glob(os.path.join(path, f"*incidenceAngleFromEllipsoid.tif"))[0]
dem_dir = os.path.join(path, "dem_tiles")
os.makedirs(dem_dir, exist_ok=True)


SAFE_file = glob.glob(os.path.join(path, f"*_manifest.safe"))[0]

dem_merge = os.path.join(path, f"DEM_merged.tif")
dem_merge_resample = os.path.join(path, f"DEM_merged_res.tif")



bbox = scene.bbox()
bbox = [bbox.extent['xmin'],bbox.extent['ymin'],bbox.extent['xmax'],bbox.extent['ymax']]
# print(box_list)




# =========================
# 1. bbox（你已经有）
# =========================
# bbox = [83.262001, 27.47629, 86.07119, 29.404852]

# =========================
# 2. STAC search Copernicus DEM
# =========================
catalog = Client.open(
    "https://planetarycomputer.microsoft.com/api/stac/v1"
)

search = catalog.search(
    collections=["cop-dem-glo-30"],
    bbox=bbox
)

items = list(search.items())

if len(items) == 0:
    raise ValueError("No DEM tiles found")

print(f"Found {len(items)} DEM tiles")

# =========================


local_files = []

for i, item in enumerate(items):

    # ⭐ 必须 sign（否则 PublicAccessNotPermitted）
    item = pc.sign(item)

    url = item.assets["data"].href

    out_path = os.path.join(dem_dir, f"dem_{i}.tif")

    print("Downloading:", url)

    r = requests.get(url, stream=True)
    r.raise_for_status()

    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    local_files.append(out_path)

print("Download complete:", len(local_files), "tiles")

# =========================
# 4. 合并 DEM（可选）
# =========================
# rasters = [rioxarray.open_rasterio(f) for f in local_files]

# dem = rioxarray.concat(rasters, dim="band").squeeze()

# dem.rio.to_raster("DEM_merged.tif")

# print("DONE -> DEM_merged.tif")





tif_list = glob.glob(os.path.join(dem_dir, "*.tif"))


gdal.Warp(
    dem_merge,
    tif_list,
    format="GTiff",
    options=gdal.WarpOptions(
        multithread=True,
        resampleAlg="bilinear",
        creationOptions=["TILED=YES", "COMPRESS=LZW"]
    )
)

print("Done -> DEM_merged.tif")






if os.path.exists(dem_dir):
    shutil.rmtree(dem_dir)

with rasterio.open(S1_angle_path) as s1_src:
    s1 = s1_src.read(1).astype(np.float32)
    s1_transform = s1_src.transform
    s1_crs = s1_src.crs
    out_shape = s1.shape
    profile = s1_src.profile.copy()

with rasterio.open(dem_merge) as dem_src:
    dem = dem_src.read(1)

    dem_reproj = np.empty(out_shape, dtype=np.float32)

    reproject(
        source=dem,
        destination=dem_reproj,
        src_transform=dem_src.transform,
        src_crs=dem_src.crs,
        dst_transform=s1_transform,
        dst_crs=s1_crs,
        resampling=Resampling.bilinear
    )


# 更新关键参数
profile.update(
    dtype="float32",
    count=1,
    compress="lzw",
    tiled=True,
    BIGTIFF="YES"
)

with rasterio.open(dem_merge_resample, "w", **profile) as dst:
    dst.write(dem_reproj.astype(np.float32), 1)


os.remove(dem_merge)











dem_resample_path = os.path.join(path, f"DEM_merged_res.tif")
LIA_path = os.path.join(path, f"{S1name}_LIA.tif")


with rasterio.open(dem_resample_path) as src:
    dem = src.read(1).astype(np.float32)
    transform = src.transform
    pixel_size_x = transform.a
    pixel_size_y = -transform.e

# =========================
# 2. slope & aspect（核心替代 terrain.products）
# =========================

dzdx = np.gradient(dem, axis=1) / pixel_size_x
dzdy = np.gradient(dem, axis=0) / pixel_size_y

slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))

aspect = np.arctan2(-dzdx, dzdy)
aspect = np.degrees(aspect)
aspect = (aspect + 90) % 360

# =========================
# 3. Sentinel-1 incidence angle
# =========================

inc = s1
inc_rad = np.deg2rad(inc)

# =========================
# 4. orbit direction → azimuth
# =========================




tree = ET.parse(SAFE_file)
root = tree.getroot()

# 不管 namespace，直接匹配 tag ending
for elem in root.iter():
    if elem.tag.endswith("pass"):
        print(elem.text)
        pass_dir = elem.text  # or DESCENDING from metadata
# pass_dir = "DESCENDING"  # or DESCENDING from metadata

AZIMUTH = 102 if pass_dir == "ASCENDING" else 282
azi_rad = np.deg2rad(AZIMUTH)

# =========================
# 5. LIA formula
# =========================

inc_sin = np.sin(inc_rad)
inc_cos = np.cos(inc_rad)
slope_sin = np.sin(slope)
slope_cos = np.cos(slope)

cos_lia = (
    inc_cos * slope_cos +
    inc_sin * slope_sin * np.cos(np.deg2rad(aspect) - azi_rad)
)

# clamp numerical issues
cos_lia = np.clip(cos_lia, -1, 1)

lia = np.degrees(np.arccos(cos_lia))

# =========================
# 6. save result
# =========================
with rasterio.open(LIA_path, "w",
                   driver="GTiff",
                   height=lia.shape[0],
                   width=lia.shape[1],
                   count=1,
                   dtype="float32",
                   crs=src.crs,
                   transform=src.transform) as dst:
    dst.write(lia.astype(np.float32), 1)

print(f"LIA saved -> {LIA_path}")

os.remove(dem_resample_path)
