import glob
import os
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import planetary_computer as pc
import rasterio
import requests
from osgeo import gdal
from pyroSAR import identify
from pystac_client import Client
from rasterio.warp import reproject, Resampling

start_time = time.perf_counter()
OUTPUT_NODATA = -9999.0


def pixel_spacing_meters(transform, crs, rows):
    """Return per-row x/y pixel spacing in meters for north-up rasters."""
    if transform.b != 0 or transform.d != 0:
        raise ValueError("Rotated raster transforms are not supported")

    if crs is None:
        raise ValueError("The incidence-angle raster has no CRS")

    if crs.is_projected:
        _, unit_factor = crs.linear_units_factor
        x_spacing = np.full(rows, abs(transform.a) * unit_factor)
        y_spacing = np.full(rows, abs(transform.e) * unit_factor)
        return x_spacing, y_spacing

    if crs.to_epsg() != 4326:
        raise ValueError(
            f"Unsupported geographic CRS for slope calculation: {crs}"
        )

    # EPSG:4326 coordinates are angular while Copernicus DEM heights are
    # meters. Convert the longitude/latitude pixel sizes to WGS84 ground
    # distances at the centre latitude of every raster row.
    row_centres = np.arange(rows, dtype=np.float64) + 0.5
    lat = transform.f + row_centres * transform.e
    lat_rad = np.deg2rad(lat)
    meters_per_degree_lat = (
        111132.92
        - 559.82 * np.cos(2 * lat_rad)
        + 1.175 * np.cos(4 * lat_rad)
        - 0.0023 * np.cos(6 * lat_rad)
    )
    meters_per_degree_lon = (
        111412.84 * np.cos(lat_rad)
        - 93.5 * np.cos(3 * lat_rad)
        + 0.118 * np.cos(5 * lat_rad)
    )
    x_spacing = abs(transform.a) * meters_per_degree_lon
    y_spacing = abs(transform.e) * meters_per_degree_lat

    if np.any(x_spacing <= 0) or np.any(y_spacing <= 0):
        raise ValueError("Invalid ground pixel spacing in the output grid")

    return x_spacing, y_spacing

if len(sys.argv) != 3:
    print(
        "Usage: python cal_LIA.py "
        "<S1_PRODUCT_NAME> <OUTPUT_FOLDER>"
    )
    sys.exit(1)

S1name = sys.argv[1]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))
path = os.path.join(folder, S1name)

# print(f"Path: {path}")

if not os.path.isdir(path):
    raise FileNotFoundError(f"S1 product folder not found: {path}")

zip_path = glob.glob(os.path.join(path, "*.zip"))[0]
# print(zip_path)
scene = identify(zip_path)

S1_angle_path = glob.glob(os.path.join(path, f"*incidenceAngleFromEllipsoid.tif"))[0]
dem_dir = os.path.join(path, "dem_tiles")
os.makedirs(dem_dir, exist_ok=True)

SAFE_file = glob.glob(os.path.join(path, f"*_manifest.safe"))[0]

dem_merge = os.path.join(path, f"DEM_merged.tif")
dem_merge_resample = os.path.join(path, f"DEM_merged_res.tif")

bbox = scene.bbox()
bbox = [
    bbox.extent["xmin"],
    bbox.extent["ymin"],
    bbox.extent["xmax"],
    bbox.extent["ymax"],
]

# Search for Copernicus DEM tiles covering the Sentinel-1 scene.
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

print(f"DEM search completed: {len(items)} tiles found.")

local_files = []

for i, item in enumerate(items):
    item = pc.sign(item)

    url = item.assets["data"].href

    out_path = os.path.join(dem_dir, f"dem_{i}.tif")

    # print("Downloading:", url)

    r = requests.get(url, stream=True)
    r.raise_for_status()

    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)

    local_files.append(out_path)

print(f"DEM download completed: {len(local_files)} tiles.")

tif_list = glob.glob(os.path.join(dem_dir, "*.tif"))
gdal.Warp(
    dem_merge,
    tif_list,
    format="GTiff",
    options=gdal.WarpOptions(
        multithread=True,
        resampleAlg="bilinear",
        creationOptions=["TILED=YES", "COMPRESS=LZW"],
    ),
)

# print("DEM mosaic completed.")

with rasterio.open(S1_angle_path) as s1_src:
    s1_masked = s1_src.read(1, masked=True).astype(np.float32)
    s1 = s1_masked.filled(np.nan)
    s1_transform = s1_src.transform
    s1_crs = s1_src.crs
    out_shape = s1.shape
    profile = s1_src.profile.copy()

with rasterio.open(dem_merge) as dem_src:
    dem_masked = dem_src.read(1, masked=True).astype(np.float32)
    dem = dem_masked.filled(np.nan)

    dem_reproj = np.full(out_shape, np.nan, dtype=np.float32)

    reproject(
        source=dem,
        destination=dem_reproj,
        src_transform=dem_src.transform,
        src_crs=dem_src.crs,
        dst_transform=s1_transform,
        dst_crs=s1_crs,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        init_dest_nodata=True,
        resampling=Resampling.bilinear,
    )


profile.update(
    dtype="float32",
    count=1,
    compress="lzw",
    tiled=True,
    BIGTIFF="YES",
    nodata=np.nan,
)

with rasterio.open(dem_merge_resample, "w", **profile) as dst:
    dst.write(dem_reproj.astype(np.float32), 1)

dem_resample_path = os.path.join(path, f"DEM_merged_res.tif")
LIA_path = os.path.join(path, f"{S1name}_LIA.tif")


with rasterio.open(dem_resample_path) as src:
    dem = src.read(1, masked=True).astype(np.float32).filled(np.nan)
    transform = src.transform
    output_crs = src.crs

# Derive terrain slope and aspect from the resampled DEM.
pixel_size_x, pixel_size_y = pixel_spacing_meters(
    transform,
    output_crs,
    dem.shape[0],
)
dzdx = np.gradient(dem, axis=1) / pixel_size_x[:, np.newaxis]
dzdy = np.gradient(dem, axis=0) / pixel_size_y[:, np.newaxis]

slope = np.arctan(np.sqrt(dzdx**2 + dzdy**2))

aspect = np.arctan2(-dzdx, dzdy)
aspect = np.degrees(aspect)
aspect = aspect % 360

inc = s1
inc_rad = np.deg2rad(inc)

tree = ET.parse(SAFE_file)
root = tree.getroot()

# Read the orbit direction without depending on the XML namespace.
for elem in root.iter():
    if elem.tag.endswith("pass"):
        # print(elem.text)
        pass_dir = elem.text

AZIMUTH = 282 if pass_dir == "ASCENDING" else 102
azi_rad = np.deg2rad(AZIMUTH)

inc_sin = np.sin(inc_rad)
inc_cos = np.cos(inc_rad)
slope_sin = np.sin(slope)
slope_cos = np.cos(slope)

cos_lia = inc_cos * slope_cos + inc_sin * slope_sin * np.cos(
    np.deg2rad(aspect) - azi_rad
)

cos_lia = np.clip(cos_lia, -1, 1)

lia = np.degrees(np.arccos(cos_lia))
valid = np.isfinite(dem) & np.isfinite(inc) & np.isfinite(lia)
lia = np.where(valid, lia, OUTPUT_NODATA)

with rasterio.open(
    LIA_path,
    "w",
    driver="GTiff",
    height=lia.shape[0],
    width=lia.shape[1],
    count=1,
    dtype="float32",
    crs=output_crs,
    transform=transform,
    nodata=OUTPUT_NODATA,
    compress="lzw",
    tiled=True,
    BIGTIFF="IF_SAFER",
) as dst:
    dst.write(lia.astype(np.float32), 1)

elapsed = time.perf_counter() - start_time
print(f"Local incidence angle completed in {elapsed:.2f} seconds.")
