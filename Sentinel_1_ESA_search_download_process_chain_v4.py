#import packages
import os
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")
import json
from shapely.geometry import shape, MultiPolygon, LineString, mapping, box
from time import perf_counter
from shapely.wkt import loads, dumps
import xml.etree.ElementTree as ET
import re
import requests
import datetime
import tarfile
from io import BytesIO
import tempfile
import shutil
import matplotlib.gridspec as gridspec
import argparse
import subprocess
import sys
import zipfile
import io
import rasterio
from rasterio.features import rasterize, geometry_window
from rasterio.windows import Window
from rasterio.errors import WindowError
from rasterio.transform import from_origin
import numpy as np
from scipy.ndimage import maximum_filter
from rasterio.features import shapes
from tqdm import tqdm
from sys import stdout
import glob
os.environ['PATH'] += ':/shared/stormcenter/zby3135/Software/snap/bin/'
import shutil
import pdb

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MIN_FLOOD_OVERLAP_KM2 = 50.0

def log_in(username,password):
    auth_url = "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token"
    auth_data = {
        "client_id": "cdse-public",
        "grant_type": "password",
        "username": username,
        "password": password
    }
    response = requests.post(auth_url, data=auth_data)
    access_token = response.json().get("access_token")

    if not access_token:
        print("Wrong username and password")
        exit()
    return access_token

def download_flood_warning_shp_from_ESA(year,month,day, directory):
    glofas_date = datetime.datetime(year, month, day).strftime("%Y%m%dT00:00Z")
    url = f"https://european-flood.emergency.copernicus.eu/api/fms/download/glofas/RapidFloodMapping/{glofas_date}"
    response = requests.get(url, timeout=120)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as z:
        z.extractall(path=directory)
    print("GloFAS data ready. Process Flood warning data:")

def simplify_flood_warning_shp_from_ESA(shapefile,window_size,area_thresholds):
    # Rasterize
    base_filename = os.path.splitext(shapefile)[0]
    gdf = gpd.read_file(shapefile)
    if gdf.empty:
        return gdf
    pixel_size = 1/111
    minx, miny, maxx, maxy = gdf.total_bounds
    width = int((maxx - minx) / pixel_size)
    height = int((maxy - miny) / pixel_size)
    transform = from_origin(minx, maxy, pixel_size, pixel_size)
    simplified_shapes = ((geom, 1) for geom in gdf.geometry)
    print("1. Rasterizing......", end="", flush=True)
    raster = rasterize(
        simplified_shapes,
        out_shape=(height, width),
        transform=transform,
        fill=0,
        dtype='uint8'
    )
    out_tif = f"{base_filename}_simplified.tif"
    with rasterio.open(
        out_tif,
        'w',
        driver='GTiff',
        height=height,
        width=width,
        count=1,
        dtype=raster.dtype,
        crs=gdf.crs,
        tiled=True, blockxsize=512, blockysize=512, compress='deflate',
        transform=transform
    ) as dst:
        dst.write(raster, 1)
    print("finished!")

    # Window filter
    print("2. Window filtering......", end="", flush=True)
    with rasterio.open(out_tif) as src:
        data = src.read(1)
        profile = src.profile
    filtered = maximum_filter(data, size=window_size, mode='constant', cval=0)
    result = np.where(filtered == 1, 1, data)
    with rasterio.open(f"{base_filename}_simplified_filtered.tif", 'w', **profile) as dst:
        dst.write(result, 1)
    print("finished!")

    # Window filter vectorize
    print("3. Window filter vectorize......", end="", flush=True)
    raster_path = f"{base_filename}_simplified_filtered.tif"
    with rasterio.open(raster_path) as src:
        image = src.read(1)
        mask = image != src.nodata
        transform = src.transform
        crs = src.crs
    results = []
    for geom, value in shapes(image, mask=mask, transform=transform):
        results.append({
            "geometry": shape(geom),
            "value": value
        })
    gdf = gpd.GeoDataFrame(results, crs=crs)
    gdf.to_file(f"{base_filename}_simplified_filtered.shp")
    print("finished!")

    # Feature area filter
    print("4. Selecting flood regions by area......", end="", flush=True)
    gdf = gpd.read_file(f"{base_filename}_simplified_filtered.shp")
    gdf_filtered = gdf[gdf["value"] == 1]
    gdf_filtered['area_km2'] = gdf_filtered.geometry.area*110*110
    gdf_filtered.to_file(f"{base_filename}_simplified_filtered_flood.shp")
    gdf_filtered_over_10000 = gdf_filtered[gdf_filtered["area_km2"] > area_thresholds]
    # gdf_filtered_over_10000 = gdf_filtered.loc[[gdf_filtered["area_km2"].idxmax()]]
    gdf_filtered_over_10000.to_file(f"{base_filename}_simplified_filtered_flood_over_{area_thresholds}.shp")
    gdf_filtered_over_10000["minx"] = gdf_filtered_over_10000.bounds.minx
    gdf_filtered_over_10000["miny"] = gdf_filtered_over_10000.bounds.miny
    gdf_filtered_over_10000["maxx"] = gdf_filtered_over_10000.bounds.maxx
    gdf_filtered_over_10000["maxy"] = gdf_filtered_over_10000.bounds.maxy
    print(f"finished! Retained {len(gdf_filtered_over_10000)} regions.", flush=True)
    # Keep original warning pixels only, tagged by their retained expanded region.
    warning_raster = os.path.abspath(f"{base_filename}_warning_regions.tif")
    gdf_filtered_over_10000 = gdf_filtered_over_10000.copy()
    gdf_filtered_over_10000["warning_region"] = np.arange(1, len(gdf_filtered_over_10000) + 1)
    gdf_filtered_over_10000["warning_raster"] = warning_raster
    if not gdf_filtered_over_10000.empty:
        print("5. Preparing original warning pixels by region...", flush=True)
        started = perf_counter()
        _write_warning_region_raster(out_tif, gdf_filtered_over_10000, warning_raster)
        print(f"Warning pixels ready in {perf_counter() - started:.1f}s.", flush=True)
    return gdf_filtered_over_10000


def _write_warning_region_raster(original_path, regions, output_path):
    """Write region IDs on original warning pixels using bounded raster tiles.

    Region polygons come from the expanded binary raster; rasterizing them back
    to that same grid preserves connected-region membership without processing
    the detailed original warning polygons. Zero denotes non-warning pixels or
    regions rejected by the existing area filter.
    """
    geometries = list(regions.geometry)
    bounds = np.array([geometry.bounds for geometry in geometries])
    ids = regions["warning_region"].to_numpy()
    with rasterio.open(original_path) as src:
        profile = src.profile.copy()
        profile.update(dtype="uint32", nodata=0, tiled=True, blockxsize=512,
                       blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
        with rasterio.open(output_path, "w", **profile) as dst:
            total = ((src.width + 511) // 512) * ((src.height + 511) // 512)
            for number, (_, window) in enumerate(dst.block_windows(1), 1):
                left, bottom, right, top = src.window_bounds(window)
                candidates = np.flatnonzero(
                    (bounds[:, 0] < right) & (bounds[:, 2] > left)
                    & (bounds[:, 1] < top) & (bounds[:, 3] > bottom)
                )
                labels = np.zeros((int(window.height), int(window.width)), dtype="uint32")
                if len(candidates):
                    original = src.read(1, window=window) == 1
                    if original.any():
                        rasterize([(geometries[i], int(ids[i])) for i in candidates],
                                  out=labels, transform=src.window_transform(window),
                                  all_touched=False)
                        labels[~original] = 0
                dst.write(labels, 1, window=window)
                if number % 100 == 0 or number == total:
                    print(f"   Warning raster tiles: {number}/{total}", flush=True)


def _has_raster_flood_overlap(product, feature):
    """Count original warning pixels for this region inside a SAR footprint."""
    footprint = _product_footprint(product)
    if footprint.is_empty:
        return False
    region_box = box(feature["minx"], feature["miny"], feature["maxx"], feature["maxy"])
    extent = box(*footprint.bounds).intersection(region_box)
    if extent.is_empty or extent.area == 0:
        return False
    with rasterio.open(feature["warning_raster"]) as src:
        try:
            window = geometry_window(src, [mapping(extent)])
        except WindowError:
            return False
        transform = src.transform
        pixel_area_km2 = abs(transform.a * transform.e - transform.b * transform.d) * 110 * 110
        count = 0
        # Bound memory even for very large or multipart scene footprints.
        row_stop = int(window.row_off + window.height)
        col_stop = int(window.col_off + window.width)
        for row in range(int(window.row_off), row_stop, 512):
            for col in range(int(window.col_off), col_stop, 512):
                tile = Window(col, row, min(512, col_stop - col), min(512, row_stop - row))
                selected = src.read(1, window=tile) == int(feature["warning_region"])
                if not selected.any():
                    continue
                coverage = rasterize([(mapping(footprint), 1)], out_shape=selected.shape,
                                     transform=src.window_transform(tile), dtype="uint8",
                                     all_touched=False)
                count += int(np.count_nonzero(selected & (coverage == 1)))
                if count * pixel_area_km2 >= MIN_FLOOD_OVERLAP_KM2:
                    return True
    return False


def _warning_geometry(feature):
    """Prefer original warnings; accept legacy geometry/bounds-only features."""
    warning_wkt = feature.get("warning_wkt")
    if isinstance(warning_wkt, str) and warning_wkt:
        return loads(warning_wkt)
    geometry = feature.get("geometry")
    if geometry is not None:
        return geometry
    return box(feature["minx"], feature["miny"], feature["maxx"], feature["maxy"])


def _product_footprint(product):
    """Read a catalogue footprint in either supported representation."""
    footprint = product.get("GeoFootprint")
    if footprint:
        footprint = shape(footprint)
    else:
        footprint = product.get("Footprint")
        if not footprint:
            raise ValueError(f"Missing footprint for Sentinel-1 product {product.get('Id')}")
        footprint = footprint.split(";", 1)[-1].strip().strip("'")
        footprint = loads(footprint)
    return footprint


def _has_flood_overlap(product, warning_geometry):
    """Preserve geometric filtering for legacy callers without raster metadata."""
    overlap = _product_footprint(product).intersection(warning_geometry)
    return not overlap.is_empty and overlap.area * 110 * 110 >= MIN_FLOOD_OVERLAP_KM2


def search_sentinel_with_shape_extent_and_data(df,year,month,day,feature):
    """Keep the legacy signature; filter candidates by >=50 km2 warning overlap.

    The catalogue bounding box is only a coarse prefilter. Original warning
    pixels belonging to this retained region determine final acceptance.
    Legacy bounds-only features still work, using their rectangle as the AOI.
    """
    use_raster = bool(feature.get("warning_raster"))
    warning_geometry = None if use_raster else _warning_geometry(feature)
    if not use_raster and warning_geometry.is_empty:
        return df

    start_date = f"{year}-{str(month).zfill(2)}-{str(day).zfill(2)}T00:00:00.000Z"
    next_day = datetime.date(year, month, day) + datetime.timedelta(days=1)
    end_date = f"{next_day.isoformat()}T00:00:00.000Z"
    coords = [
        (feature["minx"], feature["maxy"]),  # ⌈ leftup
        (feature["minx"], feature["miny"]),  # ⌊ leftdown
        (feature["maxx"], feature["miny"]),  # ⌋ rightdown
        (feature["maxx"], feature["maxy"]),  # ⌉ rightup
        (feature["minx"], feature["maxy"])   # ⌈lefup
    ]
    polygon_str = ", ".join([f"{x} {y}" for x, y in coords])
    query_url = (
        "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"
        f"$filter=OData.CSC.Intersects(area=geography'SRID=4326;POLYGON(({polygon_str}))') "
        "and Collection/Name eq 'SENTINEL-1' "
        "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq 'productType' "
        "and att/OData.CSC.StringAttribute/Value eq 'IW_GRDH_1S') "
        f"and ContentDate/Start ge {start_date} "
        f"and ContentDate/Start lt {end_date}"
    )
    while query_url:
        response = requests.get(query_url, timeout=120)
        response.raise_for_status()
        response_json = response.json()
        df2 = pd.DataFrame.from_dict([
            product for product in response_json['value']
            if (_has_raster_flood_overlap(product, feature) if use_raster
                else _has_flood_overlap(product, warning_geometry))
        ])
        if not df2.empty:
            df = pd.concat([df, df2], ignore_index=True)
        query_url = response_json.get('@odata.nextLink') or response_json.get('@OData.nextLink')
    return df


def search_flood_images_by_date_range(start_date, end_date, work_directory,
                                     window_size=10, area_thresholds=1000):
    """Return product names without .SAFE and count for daily warning-area SAR scenes.

    Dates are YYYY-MM-DD strings, inclusive in UTC. Each day's GloFAS mask
    selects that same day's Sentinel-1 IW_GRDH_1S acquisitions. Warning files
    and intermediate masks are saved under work_directory; SAR is not downloaded.
    Any missing warning archive or failed catalogue request raises an exception
    rather than returning an incomplete count as a successful result.
    """
    start = datetime.date.fromisoformat(start_date)
    end = datetime.date.fromisoformat(end_date)
    if start > end:
        raise ValueError("start_date must be on or before end_date")
    if not isinstance(window_size, int) or isinstance(window_size, bool) or window_size < 1:
        raise ValueError("window_size must be a positive integer")
    if not np.isfinite(area_thresholds) or area_thresholds < 0:
        raise ValueError("area_thresholds must be finite and nonnegative")

    products = {}
    current = start
    while current <= end:
        directory = os.path.join(work_directory, current.isoformat(), "ESA_flood_waring")
        os.makedirs(directory, exist_ok=True)
        print(f"Searching flood-warning images for {current.isoformat()}", flush=True)
        try:
            download_flood_warning_shp_from_ESA(
                current.year, current.month, current.day, directory
            )
            shapefile = os.path.join(directory, f"FloodMaskMerged{current:%Y%m%d}00.shp")
            regions = simplify_flood_warning_shp_from_ESA(shapefile, window_size, area_thresholds)
            daily = pd.DataFrame(columns=["Id", "Name"])
            for _, feature in regions.iterrows():
                daily = search_sentinel_with_shape_extent_and_data(
                    daily, current.year, current.month, current.day, feature
                )
            for product in daily[["Id", "Name"]].to_dict("records"):
                products.setdefault(product["Id"], product)
        except Exception as exc:
            raise RuntimeError(f"Flood image search failed for {current.isoformat()}: {exc}") from exc
        current += datetime.timedelta(days=1)

    names = [product["Name"].removesuffix(".SAFE") for product in products.values()]
    return {"names": names, "count": len(names)}


def download_Sentinel_with_ids_names(ids,name,output_dir,access_token):
    download_url = f"https://download.dataspace.copernicus.eu/odata/v1/Products({ids})/$value"
    download_headers = {
        "Authorization": f"Bearer {access_token}"
    }
    output_filename = f"{name}.zip"
    print(f"\nDownloading {output_filename}......", end="", flush=True)
    stdout.flush()
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(output_dir, exist_ok=True)
    with requests.get(download_url, headers=download_headers, stream=True) as r:
        if r.status_code == 200:
            total_size = int(r.headers.get('content-length', 0))
            chunk_size = 8192

            with open(output_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    f.write(chunk)


            print(f"finished!")
        else:
            print(f"Error ({r.status_code}): {r.text}")

def process_snentinel_images(file,processed_path):
    from pyroSAR.snap.util import geocode, identify, sub_parametrize

    target_resolution = 20
    terrain_flat_bool = True
    remove_therm_noise_bool = True
    fileid = identify(file)
    corners = fileid.getCorners()
    subsetnode = sub_parametrize(fileid, geometry=corners)

    geocode(
        infile=file,
        outdir=processed_path,
        spacing=int(target_resolution),
        polarizations=['VV','VH'],
        refarea='gamma0',
        t_srs=4326,
        scaling='linear',
        clean_edges=True,
        terrainFlattening=terrain_flat_bool,
        removeS1ThermalNoise=remove_therm_noise_bool,
        export_extra=['incidenceAngleFromEllipsoid', 'localIncidenceAngle', 'DEM'],
        demName='Copernicus 30m Global DEM',
        nodataValueAtSea=False,
        allow_RES_OSV=True
    )
    print

def incidence_process(VV_VH_incidence_path):
    from osgeo import gdal
    ds1 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VV_gamma0-rtc.tif")[0])
    band1 = ds1.GetRasterBand(1).ReadAsArray().astype(float)

    ds2 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VH_gamma0-rtc.tif")[0])
    band2 = ds2.GetRasterBand(1).ReadAsArray().astype(float)

    # 打开第二个影像
    ds3 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*incidenceAngleFromEllipsoid.tif")[0])
    band3 = ds3.GetRasterBand(1).ReadAsArray().astype(float)

    # 避免除零
    cos_band = np.cos(band3*3.1415926/180)
    # cos_band[cos_band == 0] = np.nan

    # 进行逐像素除法
    result1 = band1 / cos_band

    # 创建输出文件
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(f"{VV_VH_incidence_path}/Gamma0_VV.tif",
                           ds1.RasterXSize,
                           ds1.RasterYSize,
                           1,
                           gdal.GDT_Float32)

    # 设置地理信息
    out_ds.SetGeoTransform(ds1.GetGeoTransform())
    out_ds.SetProjection(ds1.GetProjection())

    # 写入数据
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result1)
    out_band.SetNoDataValue(np.nan)

    result2 = band2 / cos_band
    # 创建输出文件
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(f"{VV_VH_incidence_path}/Gamma0_VH.tif",
                           ds1.RasterXSize,
                           ds1.RasterYSize,
                           1,
                           gdal.GDT_Float32)

    # 设置地理信息
    out_ds.SetGeoTransform(ds1.GetGeoTransform())
    out_ds.SetProjection(ds1.GetProjection())

    # 写入数据
    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result2)
    out_band.SetNoDataValue(np.nan)

    # 关闭文件
    out_ds.FlushCache()
    out_ds = None
    ds1 = None
    ds2 = None

    pattern = os.path.join(VV_VH_incidence_path, 'S1A*')

    for file in glob.glob(pattern):
        try:
            os.remove(file)

        except Exception as e:
            a=1


def run_new_processing_chain(product_name, output_dir, desert_mask_vrt):
    """Run the maintained per-product workflow and its ancillary masks."""
    s1name = product_name.removesuffix(".SAFE")
    product_dir = os.path.join(output_dir, s1name)
    stages = [
        (
            "Sentinel-1 download and preprocessing",
            "Sentinel_1_specific_name_download_process.py",
            [s1name, output_dir, username, password],
        ),
        (
            "desert mask",
            "Desert_mask.py",
            [s1name, output_dir, desert_mask_vrt],
        ),
    ]

    for label, script_name, arguments in stages:
        print(f"Running {label} for {s1name}...", flush=True)
        subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, script_name), *arguments],
            check=True,
        )

    incidence_paths = glob.glob(
        os.path.join(product_dir, "*incidenceAngleFromEllipsoid.tif")
    )
    if len(incidence_paths) != 1:
        raise RuntimeError(
            "Expected exactly one ellipsoid incidence-angle raster, found "
            f"{len(incidence_paths)}"
        )
    ancillary_stages = [
        (
            "local incidence angle",
            "cal_LIA.py",
            [
                s1name,
                product_dir,
                "--incidence-angle",
                incidence_paths[0],
                "--metadata-dir",
                product_dir,
            ],
        ),
        (
            "snow and cloud masks",
            "Snow_detect.py",
            [s1name, os.path.join(product_dir, "Gamma0_VV.tif"), product_dir],
        ),
    ]
    for label, script_name, arguments in ancillary_stages:
        print(f"Running {label} for {s1name}...", flush=True)
        subprocess.run(
            [sys.executable, os.path.join(SCRIPT_DIR, script_name), *arguments],
            check=True,
        )

    required_outputs = [
        "Gamma0_VV.tif",
        "Gamma0_VH.tif",
        f"{s1name}_desert.tif",
        f"{s1name}_LIA.tif",
        f"{s1name}_ice.tif",
        f"{s1name}_cloud.tif",
    ]
    missing = [
        filename
        for filename in required_outputs
        if not os.path.isfile(os.path.join(product_dir, filename))
    ]
    if missing:
        raise RuntimeError("Missing expected outputs: " + ", ".join(missing))

    # Match execute.sh: only clean intermediates after every stage succeeds.
    for directory_name in ("snow_temp", "dem_tiles"):
        shutil.rmtree(os.path.join(product_dir, directory_name), ignore_errors=True)
    # Keep *_DEM.tif exported by SNAP for LIA reuse and inspection.
    cleanup_patterns = [
        "*.zip",
        "*incidenceAngleFromEllipsoid.tif",
        "*localIncidenceAngle.tif",
        "*_manifest.safe",
        "*_proc.xml",
        "*_gamma0-rtc.tif",
        "DEM_merged.tif",
        "DEM_merged_res.tif",
    ]
    for pattern in cleanup_patterns:
        for filename in glob.glob(os.path.join(product_dir, pattern)):
            if os.path.isfile(filename):
                os.remove(filename)

    print(f"Completed {s1name}; intermediate files cleaned.", flush=True)






def main():
    global username, password
    parser = argparse.ArgumentParser(
        description="Process today's GloFAS flood warnings with Sentinel-1 data."
    )
    parser.add_argument(
        "--desert-mask-vrt",
        help="Path to the global desert-mask VRT used by Desert_mask.py.",
    )
    parser.add_argument("--search-only", action="store_true", help="List flood-area SAR products without downloading SAR.")
    parser.add_argument("--start-date", help="First UTC date, YYYY-MM-DD (inclusive).")
    parser.add_argument("--end-date", help="Last UTC date, YYYY-MM-DD (inclusive).")
    parser.add_argument("--work-directory", default=os.path.join(SCRIPT_DIR, "data"))
    parser.add_argument("--output-json", help="Save product names without .SAFE and count to this JSON file.")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--area-thresholds", type=float, default=1000)
    args = parser.parse_args()
    if args.search_only:
        if not args.start_date or not args.end_date:
            parser.error("--search-only requires --start-date and --end-date")
        result = search_flood_images_by_date_range(
            args.start_date, args.end_date, args.work_directory,
            window_size=args.window_size, area_thresholds=args.area_thresholds,
        )
        output = json.dumps(result, indent=2)
        if args.output_json:
            with open(args.output_json, "w", encoding="utf-8") as stream:
                stream.write(output + "\n")
        print(output)
        return
    if args.start_date or args.end_date or args.output_json:
        parser.error("Date-range and JSON output options require --search-only")
    if not args.desert_mask_vrt:
        parser.error("--desert-mask-vrt is required for processing")
    username = os.environ.get("CDSE_USERNAME", "")
    password = os.environ.get("CDSE_PASSWORD", "")
    if not username or not password:
        raise RuntimeError(
            "Set CDSE_USERNAME and CDSE_PASSWORD before running the flood-warning "
            "workflow."
        )
    desert_mask_vrt = os.path.abspath(os.path.expanduser(args.desert_mask_vrt))
    if not os.path.isfile(desert_mask_vrt):
        parser.error(f"Desert mask VRT not found: {desert_mask_vrt}")


    # 获取当前时间
    now =  datetime.datetime.now()

    # 年、月、日
    year = now.year
    month = now.month
    day = now.day

    # year = 2025
    # month = 10
    # day = 29

    project_root = os.path.abspath(
        os.path.expanduser(os.environ.get("RAPID_PROJECT_DIR", SCRIPT_DIR))
    )
    workfolder = os.path.join(
        project_root,
        "data",
        f"{year}-{str(month).zfill(2)}-{str(day).zfill(2)}",
    )
    os.makedirs(workfolder, exist_ok=True)


    window_size = args.window_size
    area_thresholds = args.area_thresholds

    # Main program

    # path set up
    flood_waring_directory = os.path.join(workfolder, "ESA_flood_waring")
    os.makedirs(flood_waring_directory, exist_ok=True)
    Sentinel_process_dir = os.path.join(workfolder, "Sentinel_1")
    os.makedirs(Sentinel_process_dir, exist_ok=True)
    processed_images_dir = os.path.join(Sentinel_process_dir, "processed_images")
    os.makedirs(processed_images_dir, exist_ok=True)
    glofas_date_name = datetime.datetime(year, month, day).strftime("%Y%m%d")
    shapefile = os.path.join(
        flood_waring_directory,
        f"FloodMaskMerged{glofas_date_name}00.shp",
    )
    Processed_Sentinel_1_data_path_filename = os.path.join(
        workfolder,
        f"Processed_Sentinel_1_data_path_{glofas_date_name}.txt",
    )
    if os.path.exists(Processed_Sentinel_1_data_path_filename):
        os.remove(Processed_Sentinel_1_data_path_filename)

    # download ESA flood warning data
    download_flood_warning_shp_from_ESA(year,month,day, flood_waring_directory)
    gdf = simplify_flood_warning_shp_from_ESA(shapefile,window_size,area_thresholds)
    df = pd.DataFrame(columns=["Id", "Name"])

    # Search Sentinel-1 images with extent and date
    print("\nSearching images......", end="", flush=True)
    for idx, feature in gdf.iterrows():
        df = search_sentinel_with_shape_extent_and_data(df,year,month,day,feature)
    print(f"finished!")
    if not df.empty:
        df = df.drop_duplicates(subset=["Name"]).reset_index(drop=True)
    print(f"Found {len(df)} images")
    for idx, row in df[['Id','Name']].iterrows():
        print(f"No.{idx+1}", row['Id'], row['Name'])

    # download images and process images to gamm0
    if len(df)==0:
        print(f"No images on {glofas_date_name}")
        # shutil.rmtree(workfolder)
    else:

        for idx, row in df.iterrows():
            name = row["Name"]
            try:
                start_time_each_image = datetime.datetime.now()
                run_new_processing_chain(name, processed_images_dir, desert_mask_vrt)
                interval_time = datetime.datetime.now()
                print(f"Finished in {interval_time - start_time_each_image}")

            except Exception as e:
                print(f"Error processing {name}: {e}")
                continue
        with open(Processed_Sentinel_1_data_path_filename, "a", encoding="utf-8") as f:
            f.write(processed_images_dir + os.sep + "\n")

    total_time = datetime.datetime.now()
    print(f"\n{year}-{str(month).zfill(2)}-{str(day).zfill(2)} triger finished in " +str(total_time - now))


if __name__ == "__main__":
    main()
