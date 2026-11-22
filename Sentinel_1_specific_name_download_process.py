import argparse
from contextlib import contextmanager
import datetime
import glob
import os
import shutil
import sys
from sys import stdout

import numpy as np
import pandas as pd
import requests

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
        raise RuntimeError("Copernicus Data Space authentication failed")
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


# CDSE allows 4 concurrent download connections per user (documentation, Quotas).
CDSE_MAX_DOWNLOADS = int(os.environ.get("CDSE_MAX_DOWNLOADS", "4"))


@contextmanager
def cdse_download_slot():
    """Hold one of CDSE_MAX_DOWNLOADS lock files shared by all of this user's processes,
    so parallel scenes and array tasks stay within the CDSE connection limit."""
    import fcntl
    import time

    lock_dir = os.path.expanduser(os.environ.get("CDSE_LOCK_DIR", "~/.cache/s1_download/cdse_slots"))
    os.makedirs(lock_dir, exist_ok=True)
    while True:
        for i in range(CDSE_MAX_DOWNLOADS):
            handle = open(os.path.join(lock_dir, f"slot{i}"), "w")
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                handle.close()
                continue
            try:
                yield
            finally:
                fcntl.flock(handle, fcntl.LOCK_UN)
                handle.close()
            return
        time.sleep(10)


@contextmanager
def user_lock(name):
    """Exclusive lock shared by all of this user's processes on the machine."""
    import fcntl

    lock_dir = os.path.expanduser(os.environ.get("CDSE_LOCK_DIR", "~/.cache/s1_download/cdse_slots"))
    os.makedirs(lock_dir, exist_ok=True)
    with open(os.path.join(lock_dir, f"{name}.lock"), "w") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def prime_snap_version_cache(ExamineSnap):
    """Fill pyroSAR's SNAP module-version cache from <SNAP>/*/update_tracking/*.xml.

    These are the installed-module records the launcher itself reports, so the
    result matches `snap --modules --list` without starting SNAP. Returns False
    (cache untouched) if the records or the microwave toolbox are not found.
    """
    import xml.etree.ElementTree as ET

    if ExamineSnap._version_dict is not None:
        return True
    try:
        root = os.path.dirname(os.path.dirname(os.path.realpath(ExamineSnap().path)))
    except Exception:
        return False
    versions = {}
    for record in glob.glob(os.path.join(root, "*", "update_tracking", "*.xml")):
        try:
            module = ET.parse(record).getroot()
            release = module.find("module_version")
            versions[module.get("codename")] = {"version": release.get("specification_version"),
                                                "state": "Installed"}
        except Exception:
            continue
    # Clusters without update_tracking (e.g. snap/, which holds snap-core): the
    # version is the module JAR's OpenIDE-Module-Specification-Version.
    import zipfile
    for config in glob.glob(os.path.join(root, "*", "config", "Modules", "*.xml")):
        try:
            module = ET.parse(config).getroot()
            code = module.get("name")
            jar = module.find("param[@name='jar']")
            if code in versions or jar is None:
                continue
            cluster = os.path.dirname(os.path.dirname(os.path.dirname(config)))
            with zipfile.ZipFile(os.path.join(cluster, jar.text)) as archive:
                manifest = archive.read("META-INF/MANIFEST.MF").decode("utf-8", "replace")
            for line in manifest.splitlines():
                if line.startswith("OpenIDE-Module-Specification-Version:"):
                    versions[code] = {"version": line.split(":", 1)[1].strip(), "state": "Installed"}
                    break
        except Exception:
            continue
    if "eu.esa.microwavetbx.microwavetbx.kit" not in versions:
        return False
    ExamineSnap._version_dict = versions
    return True


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

            # Write to .part first: an interrupted download must not look complete.
            partial_path = output_path + ".part"
            with open(partial_path, 'wb') as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    f.write(chunk)
            os.replace(partial_path, output_path)


            print(" completed.")
        else:
            raise RuntimeError(f"Sentinel-1 download failed with status {r.status_code}.")

def process_snentinel_images(file, processed_path):
    from pyroSAR.snap.auxil import ExamineSnap
    from pyroSAR.snap.util import geocode

    # pyroSAR reads the SNAP version through the SNAP GUI launcher
    # (snap --nogui --modules --list --refresh): one instance per user, a network
    # refresh, and on compute nodes it can stall for hours. Fill pyroSAR's
    # version cache from the installed-module records instead; fall back to the
    # launcher only if those cannot be read.
    if not prime_snap_version_cache(ExamineSnap):
        with user_lock("snap_version"):
            ExamineSnap().get_version("microwavetbx")

    target_resolution = 20
    terrain_flat_bool = True
    remove_therm_noise_bool = True

    # Optional per-process SNAP limits, so several scenes can share a node:
    # SNAP_THREADS -> gpt -q; SNAP_HEAP (e.g. 24G) -> Java -Xmx (the SNAP install sets 528G).
    gpt_args = ["-q", os.environ["SNAP_THREADS"]] if os.environ.get("SNAP_THREADS") else None
    if os.environ.get("SNAP_HEAP"):
        os.environ["_JAVA_OPTIONS"] = f"-Xmx{os.environ['SNAP_HEAP']}"

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
        gpt_args=gpt_args,
    )


def incidence_process(VV_VH_incidence_path):
    from osgeo import gdal

    ds1 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VV_gamma0-rtc.tif")[0])
    band1 = ds1.GetRasterBand(1).ReadAsArray().astype(float)

    ds2 = gdal.Open(glob.glob(f"{VV_VH_incidence_path}/*VH_gamma0-rtc.tif")[0])
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

def main(argv=None):
    parser = argparse.ArgumentParser(description="Download and process Sentinel-1 scenes.")
    parser.add_argument("names", help="Comma-separated Sentinel-1 product names")
    parser.add_argument("output_folder")
    parser.add_argument("username", nargs="?", default="")
    parser.add_argument("password", nargs="?", default="")
    parser.add_argument("--download-source", choices=("cdse", "asf"), default="cdse",
                        help="cdse (CDSE + SNAP, Copernicus 30 m DEM; default) or asf (HyP3 RTC)")
    args = parser.parse_args(argv)
    S1names = [name.strip().removesuffix(".SAFE") for name in args.names.split(",") if name.strip()]
    if not S1names:
        parser.error("At least one product name is required")
    folder = os.path.abspath(os.path.expanduser(args.output_folder))
    if args.download_source == "asf":
        from Sentinel_1_specific_name_ASF import download_and_process
        return 1 if download_and_process(S1names, folder) else 0
    username = args.username or os.environ.get("CDSE_USERNAME", "")
    password = args.password or os.environ.get("CDSE_PASSWORD", "")
    if not username or not password:
        parser.error("Copernicus Data Space username and password are required")
    os.makedirs(folder, exist_ok=True)
    print(f"Processing: {S1names}")
    failures = 0
    for S1name in S1names:
        workfolder = os.path.join(folder, S1name)
        os.makedirs(workfolder, exist_ok=True)

        df = search_sentinel_with_S1name(S1name)
        access_token = log_in(username, password)

        if len(df) == 0:
            print("No images found")
            failures += 1
        else:
            for idx, row in df.iterrows():
                try:
                    start_time_each_image = datetime.datetime.now()
                    ids = row["Id"]
                    name = row["Name"]
                    Sentinel_ori_dir = workfolder
                    Sentinel_1_GRD_file = os.path.join(workfolder, f'{name}.zip')
                    if not os.path.exists(Sentinel_1_GRD_file):
                        # Log in inside the slot: a token expires after 10 min of waiting.
                        with cdse_download_slot():
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
                    failures += 1
                    continue

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
