import datetime
import glob
import os
import shutil
import sys
from sys import stdout

import numpy as np
import pandas as pd
import requests
from osgeo import gdal
from pyroSAR.snap.util import geocode

os.environ['PATH'] += ':/tank/data/SFS/xinyis/shared/apps/esa-snap/bin'


def log_in(username, password):
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

def search_sentinel_with_S1name(S1name):
    S1name_with_safe = S1name + ".SAFE"
    query_url = (
        "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"
        f"$filter=Name eq '{S1name_with_safe}'"
    )
    response = requests.get(query_url)
    response_json = response.json()
    df = pd.DataFrame.from_dict(response_json['value'])
    
    return df


def download_Sentinel_with_ids_names(ids, name, output_dir, access_token):
    download_url = f"https://download.dataspace.copernicus.eu/odata/v1/Products({ids})/$value"
    download_headers = {
        "Authorization": f"Bearer {access_token}"
    }
    output_filename = f"{name}.zip"
    print("Downloading Sentinel-1 product...", end="", flush=True)
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


            print(" completed.")
        else:
            print(f"Sentinel-1 download failed with status {r.status_code}.")

def process_snentinel_images(file, processed_path):
    target_resolution = 20
    terrain_flat_bool = False
    remove_therm_noise_bool = True

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
        export_extra=['incidenceAngleFromEllipsoid', 'localIncidenceAngle'],
        demName='ACE30',
        nodataValueAtSea=False,
    )


def incidence_process(VV_VH_incidence_path):
    ds1 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VV_gamma0-elp.tif")[0])
    band1 = ds1.GetRasterBand(1).ReadAsArray().astype(float)

    ds2 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VH_gamma0-elp.tif")[0])
    band2 = ds2.GetRasterBand(1).ReadAsArray().astype(float)

    ds3 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*incidenceAngleFromEllipsoid.tif")[0])
    band3 = ds3.GetRasterBand(1).ReadAsArray().astype(float)

    # Normalize backscatter using the ellipsoid incidence angle.
    cos_band = np.cos(np.deg2rad(band3))

    result1 = band1 / cos_band

    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(f"{VV_VH_incidence_path}/Gamma0_VV.tif",
                           ds1.RasterXSize,
                           ds1.RasterYSize,
                           1,
                           gdal.GDT_Float32)

    out_ds.SetGeoTransform(ds1.GetGeoTransform())
    out_ds.SetProjection(ds1.GetProjection())

    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result1)
    out_band.SetNoDataValue(np.nan)

    result2 = band2 / cos_band
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(f"{VV_VH_incidence_path}/Gamma0_VH.tif",
                           ds1.RasterXSize,
                           ds1.RasterYSize,
                           1,
                           gdal.GDT_Float32)

    out_ds.SetGeoTransform(ds1.GetGeoTransform())
    out_ds.SetProjection(ds1.GetProjection())

    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(result2)
    out_band.SetNoDataValue(np.nan)

    out_ds.FlushCache()
    out_ds = None
    ds1 = None
    ds2 = None

if len(sys.argv) != 5:
    print(
        "Usage: python Sentinel_1_specific_name_download_process.py "
        "<S1_PRODUCT_NAME1,S1_PRODUCT_NAME2,...> <OUTPUT_FOLDER> "
        "<USERNAME> <PASSWORD>"
    )
    sys.exit(1)

S1names = [name.strip() for name in sys.argv[1].split(",") if name.strip()]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))
username = sys.argv[3]
password = sys.argv[4]

if not username or not password:
    print("Copernicus Data Space username and password are required.")
    sys.exit(1)

os.makedirs(folder, exist_ok=True)

print(f"Processing: {S1names}")
# print(f"Output folder: {folder}")

for S1name in S1names:
    workfolder = os.path.join(folder, S1name)
    os.makedirs(workfolder, exist_ok=True)

    df = search_sentinel_with_S1name(S1name)
    access_token = log_in(username, password)

    if len(df) == 0:
        print("No images found")
        shutil.rmtree(workfolder)
    else:
        for idx, row in df.iterrows():
            try:
                start_time_each_image = datetime.datetime.now()
                ids = row["Id"]
                name = row["Name"]
                Sentinel_ori_dir = workfolder
                Sentinel_1_GRD_file = os.path.join(workfolder, f'{name}.zip')
                if not os.path.exists(Sentinel_1_GRD_file):
                    access_token = log_in(username, password)
                    download_Sentinel_with_ids_names(ids, name, Sentinel_ori_dir, access_token)
                else:
                    # print("Product already downloaded.")
                    pass

                print("Preprocessing SAR imagery...", end="", flush=True)
                process_snentinel_images(Sentinel_1_GRD_file, workfolder)
                incidence_process(workfolder)

                interval_time = datetime.datetime.now()
                print(f" completed in {interval_time - start_time_each_image}.")

            except Exception as e:
                print(f"Error processing {name}: {e}")
                continue
