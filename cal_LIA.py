import argparse
import glob
import os
import shutil
import time
import xml.etree.ElementTree as ET
import zipfile

import numpy as np
import planetary_computer as pc
import rasterio
import requests
from osgeo import gdal
from pystac_client import Client
from rasterio.warp import reproject, Resampling, transform_bounds

start_time = time.perf_counter()
MASK_NODATA = 255
LIA_THRESHOLD_DEGREES = 50.0


def normalize_orbit_direction(value):
    if value is None:
        return None
    value = str(value).strip().upper()
    aliases = {"A": "ASCENDING", "D": "DESCENDING"}
    value = aliases.get(value, value)
    if value not in {"ASCENDING", "DESCENDING"}:
        return None
    return value


def direction_from_xml(xml_content):
    root = ET.fromstring(xml_content)
    for elem in root.iter():
        if elem.tag.endswith("pass"):
            direction = normalize_orbit_direction(elem.text)
            if direction:
                return direction
    return None


def direction_from_local_metadata(scene_dir):
    manifests = sorted(glob.glob(os.path.join(scene_dir, "*manifest.safe")))
    for manifest in manifests:
        with open(manifest, "rb") as source:
            direction = direction_from_xml(source.read())
        if direction:
            print(f"Orbit direction read from manifest: {direction}")
            return direction

    for archive_path in sorted(glob.glob(os.path.join(scene_dir, "*.zip"))):
        try:
            with zipfile.ZipFile(archive_path) as archive:
                manifest_names = [
                    name for name in archive.namelist()
                    if name.lower().endswith("manifest.safe")
                ]
                for name in manifest_names:
                    direction = direction_from_xml(archive.read(name))
                    if direction:
                        print(f"Orbit direction read from ZIP: {direction}")
                        return direction
        except (OSError, zipfile.BadZipFile, ET.ParseError) as exc:
            print(f"Could not read orbit direction from {archive_path}: {exc}")
    return None


def direction_from_catalog(product_name):
    safe_name = product_name if product_name.endswith(".SAFE") else f"{product_name}.SAFE"
    escaped_name = safe_name.replace("'", "''")
    url = (
        "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
        f"?$filter=Name eq '{escaped_name}'&$expand=Attributes"
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    products = response.json().get("value", [])
    for product in products:
        for attribute in product.get("Attributes", []):
            if str(attribute.get("Name", "")).lower() == "orbitdirection":
                direction = normalize_orbit_direction(attribute.get("Value"))
                if direction:
                    print(f"Orbit direction read from Copernicus catalogue: {direction}")
                    return direction
    return None


def resolve_orbit_direction(product_name, scene_dir, user_direction=None):
    direction = direction_from_local_metadata(scene_dir)
    if direction:
        return direction
    try:
        direction = direction_from_catalog(product_name)
    except (requests.RequestException, ValueError) as exc:
        print(f"Copernicus orbit-direction lookup failed: {exc}")
        direction = None
    if direction:
        return direction
    direction = normalize_orbit_direction(user_direction)
    if direction:
        print(f"Using user-provided orbit direction: {direction}")
        return direction
    raise ValueError(
        "Orbit direction is required for LIA. Provide manifest.safe, the "
        "original ZIP, catalogue access, or --orbit-direction ASCENDING/DESCENDING."
    )


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

parser = argparse.ArgumentParser(
    description="Generate a binary LIA > 50 degree exclusion mask."
)
parser.add_argument("s1_product_name")
parser.add_argument("work_dir", help="Writable per-scene preprocessing directory.")
parser.add_argument(
    "--incidence-angle",
    help=(
        "Read-only ellipsoid incidence-angle raster. If omitted, use the "
        "legacy layout <WORK_DIR>/<S1_PRODUCT_NAME>/ and discover it there."
    ),
)
parser.add_argument(
    "--metadata-dir",
    help="Read-only directory that may contain manifest.safe or the S1 ZIP.",
)
parser.add_argument(
    "--orbit-direction", choices=["ASCENDING", "DESCENDING", "A", "D"],
    help="Fallback when orbit direction cannot be read or queried.",
)
cli_args = parser.parse_args()

S1name = cli_args.s1_product_name
path = os.path.abspath(os.path.expanduser(cli_args.work_dir))
legacy_cli = cli_args.incidence_angle is None
if legacy_cli:
    path = os.path.join(path, S1name)
    incidence_candidates = sorted(
        glob.glob(os.path.join(path, "*incidenceAngleFromEllipsoid.tif"))
    )
    if len(incidence_candidates) != 1:
        raise FileNotFoundError(
            "Legacy invocation expected exactly one ellipsoid incidence-angle "
            f"raster in {path}, found {len(incidence_candidates)}"
        )
    cli_args.incidence_angle = incidence_candidates[0]
metadata_dir = (
    os.path.abspath(os.path.expanduser(cli_args.metadata_dir))
    if cli_args.metadata_dir else path
)
os.makedirs(path, exist_ok=True)

S1_angle_path = os.path.abspath(os.path.expanduser(cli_args.incidence_angle))
if not os.path.isfile(S1_angle_path):
    raise FileNotFoundError(
        f"Ellipsoid incidence-angle raster not found: {S1_angle_path}"
    )
dem_dir = os.path.join(path, "dem_tiles")
os.makedirs(dem_dir, exist_ok=True)

dem_merge = os.path.join(path, f"DEM_merged.tif")
dem_merge_resample = os.path.join(path, f"DEM_merged_res.tif")

with rasterio.open(S1_angle_path) as angle_source:
    if angle_source.crs is None:
        raise ValueError(f"Incidence-angle raster has no CRS: {S1_angle_path}")
    bbox = transform_bounds(
        angle_source.crs, "EPSG:4326", *angle_source.bounds, densify_pts=21
    )

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
    blockxsize=256,
    blockysize=256,
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

pass_dir = resolve_orbit_direction(
    S1name, metadata_dir, cli_args.orbit_direction
)
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
lia_mask = np.full(lia.shape, MASK_NODATA, dtype=np.uint8)
lia_mask[valid] = (lia[valid] > LIA_THRESHOLD_DEGREES).astype(np.uint8)

with rasterio.open(
    LIA_path,
    "w",
    driver="GTiff",
    height=lia.shape[0],
    width=lia.shape[1],
    count=1,
    dtype="uint8",
    crs=output_crs,
    transform=transform,
    nodata=MASK_NODATA,
    compress="deflate",
    tiled=True,
    blockxsize=256,
    blockysize=256,
    BIGTIFF="IF_SAFER",
) as dst:
    dst.write(lia_mask, 1)
    dst.set_band_description(1, "local_incidence_angle_greater_than_50_degrees")
    dst.update_tags(
        mask_semantics="0=keep,1=exclude,255=nodata",
        lia_threshold_degrees=LIA_THRESHOLD_DEGREES,
        orbit_direction=pass_dir,
        radar_azimuth_degrees=AZIMUTH,
    )

shutil.rmtree(dem_dir)
for temporary_path in (dem_merge, dem_merge_resample):
    if os.path.isfile(temporary_path):
        os.remove(temporary_path)

elapsed = time.perf_counter() - start_time
print(f"Local incidence angle completed in {elapsed:.2f} seconds.")
