#import packages
import os 
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings("ignore")
import json
from shapely.geometry import shape, MultiPolygon, LineString, mapping
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
import zipfile
import io
import wget
import rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
import numpy as np
from scipy.ndimage import maximum_filter
from rasterio.features import shapes
# from tqdm import tqdm
from sys import stdout
import glob
os.environ['PATH'] += ':/tank/data/SFS/xinyis/shared/apps/esa-snap/bin'
#os.system('cls' if os.name == 'nt' else 'clear')
import shutil
from pyroSAR.snap.util import geocode,ID,identify,sub_parametrize
# from pyroSAR.snap.util import geocode
import pdb
from osgeo import gdal

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

def search_sentinel_with_S1name(S1name):
    S1name_with_safe = S1name + '.SAFE'
    query_url = f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Name eq '{S1name_with_safe}'"
    # print(query_url)
    response = requests.get(query_url)
    response_json = response.json()
    df = pd.DataFrame.from_dict(response_json['value'])
    
    return df


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

    target_resolution = 20
    terrain_flat_bool = False
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
        export_extra=['incidenceAngleFromEllipsoid', 'localIncidenceAngle'],
        demName='ACE30',
        nodataValueAtSea=False
    )

def incidence_process(VV_VH_incidence_path):
    ds1 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VV_gamma0-elp.tif")[0])
    band1 = ds1.GetRasterBand(1).ReadAsArray().astype(float)

    ds2 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VH_gamma0-elp.tif")[0])
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
            #os.remove(file)
            pass
        except Exception as e:
            a=1



    
# parameters
username = "maopuxu@uwm.edu"
password = "Xmp*20021109"


# acquire the current time
now =  datetime.datetime.now()

import sys

# -------------------------
# Parse input arguments
# -------------------------
if len(sys.argv) != 3:
    print(
        "Usage: python Sentinel_1_specific_name_download_process.py "
        "<S1_PRODUCT_NAME1,S1_PRODUCT_NAME2,...> <OUTPUT_FOLDER>"
    )
    sys.exit(1)

S1names = [name.strip() for name in sys.argv[1].split(",") if name.strip()]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))

os.makedirs(folder, exist_ok=True)

print(f"Processing: {S1names}")
print(f"Output folder: {folder}")

for S1name in S1names:
    workfolder = os.path.join(folder, S1name)
    os.makedirs(workfolder, exist_ok=True) 

    df = search_sentinel_with_S1name(S1name)
    # print(df)
    access_token = log_in(username,password)
  
    # download images and process images to gamm0
    if len(df)==0:
        print(f"No images on ")
        shutil.rmtree(workfolder)
    else:       
        # for idx, row in df.iloc[:3].iterrows():
        for idx, row in df.iterrows():
            try:
                start_time_each_image = datetime.datetime.now()
                ids = row["Id"]
                name = row["Name"]
                Sentinel_ori_dir = workfolder
                Sentinel_1_GRD_file = os.path.join(workfolder, f'{name}.zip')
                if not os.path.exists(Sentinel_1_GRD_file):
                    access_token = log_in(username,password)
                    download_Sentinel_with_ids_names(ids, name, Sentinel_ori_dir, access_token)
                else:
                    print("Downloaded!")    
            
                print("Processing to Gamma0......", end="", flush=True )
                process_snentinel_images(Sentinel_1_GRD_file,workfolder)
            
                incidence_process(workfolder)            
            
                interval_time = datetime.datetime.now()
                print(f"finished in " +str(interval_time - start_time_each_image))

            except Exception as e:
                print(f"Error processing {name}: {e}")
                continue

        
