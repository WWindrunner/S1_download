"""ASF HyP3 RTC adapter for the maintained Sentinel-1 processing chain.

Authentication uses Earthdata credentials in ~/.netrc. Importing this module
never submits jobs. HyP3 local incidence angles are radians, not ellipsoid angles;
the ellipsoid incidence angle is interpolated from the geolocation grid of the GRD
annotation (the source of SNAP's incidenceAngleFromEllipsoid), read remotely from
the ASF GRD zip. Gamma0 is left as HyP3 RTC; RAPID normalizes it with this angle
(--inc_angle_path).
"""
import argparse
from pathlib import Path
import shutil
import tempfile
import zipfile


def extract_products(archives, directory):
    """Extract only recognized rasters; reject missing or ambiguous layers."""
    suffixes = {"vv": "_vv.tif", "vh": "_vh.tif",
                "angle": "_inc_map.tif", "dem": "_dem.tif", "ls": "_ls_map.tif"}
    products = {}
    for archive in archives:
        with zipfile.ZipFile(archive) as package:
            for member in package.infolist():
                for key, suffix in suffixes.items():
                    if member.filename.lower().endswith(suffix):
                        if key in products:
                            raise ValueError(f"Multiple ASF rasters for {key}")
                        # Do not use archive paths as filesystem paths.
                        target = Path(directory) / (key + ".tif")
                        with package.open(member) as source, target.open("wb") as dest:
                            shutil.copyfileobj(source, dest)
                        products[key] = target
    missing = set(suffixes) - products.keys()
    if missing:
        raise ValueError(f"Missing ASF RTC layers: {sorted(missing)}")
    return products


def annotation_points(granule):
    """(lon, lat, ellipsoid incidence deg) of the GRD annotation geolocation grid."""
    import netrc
    import xml.etree.ElementTree as ET

    import asf_search as asf
    import numpy as np

    def read():
        user, _, password = netrc.netrc().authenticators("urs.earthdata.nasa.gov")
        session = asf.ASFSession().auth_with_creds(user, password)
        products = [p for p in asf.granule_search([granule])
                    if p.properties["processingLevel"].startswith("GRD")]
        if not products:
            raise ValueError(f"GRD not found in the ASF catalogue: {granule}")
        # Only the annotation (~1.6 MB) is read from the GRD zip; one polarization suffices.
        with products[0].remotezip(session) as package:
            names = [n for n in package.namelist() if "/annotation/" in n and n.endswith(".xml")
                     and "/calibration/" not in n and "/rfi/" not in n]
            if not names:
                raise ValueError(f"No annotation XML in the GRD zip of {granule}")
            return package.read(sorted(names)[0])

    root = ET.fromstring(retry_transient(read, f"ASF annotation for {granule}"))
    points = np.array([[float(p.find(k).text) for k in ("longitude", "latitude", "incidenceAngle")]
                       for p in root.iter("geolocationGridPoint")])
    if len(points) < 3:
        raise ValueError(f"Annotation geolocation grid too small for {granule}")
    return points


def ellipsoid_incidence(points, transform, width, height, step=50):
    """Interpolate annotation points onto a north-up EPSG:4326 grid (degrees)."""
    import numpy as np
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    # The angle is smooth: interpolate every `step` pixels, then bilinear.
    rows = np.unique(np.append(np.arange(0, height, step), height - 1))
    cols = np.unique(np.append(np.arange(0, width, step), width - 1))
    lon = transform.c + (cols + 0.5) * transform.a
    lat = transform.f + (rows + 0.5) * transform.e
    glon, glat = np.meshgrid(lon, lat)
    nodes = LinearNDInterpolator(points[:, :2], points[:, 2])(glon, glat)
    gap = np.isnan(nodes)  # outside the point hull, e.g. scene edges
    nodes[gap] = NearestNDInterpolator(points[:, :2], points[:, 2])(glon[gap], glat[gap])
    # separable linear interpolation: along columns, then along rows
    across = np.array([np.interp(np.arange(width), cols, row) for row in nodes], dtype="float32")
    out = np.empty((height, width), dtype="float32")
    for c0 in range(0, width, 4096):
        part = across[:, c0:c0 + 4096]
        out[:, c0:c0 + 4096] = np.array([np.interp(np.arange(height), rows, col)
                                         for col in part.T], dtype="float32").T
    return out


def normalize_products(products, product_dir, granule, incidence_points=None):
    """Write float32 gamma0/DEM/angles and uint8 LIA on one WGS84 grid.

    Gamma0 is not normalized here: the ellipsoid incidence angle is written for
    RAPID's --inc_angle_path. The LIA mask excludes LIA > 50
    degrees and HyP3 layover-shadow shadow flags (8, 16), so shadow is excluded even
    where the local incidence angle is folded below 90 degrees.
    """
    import numpy as np
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

    if incidence_points is None:
        incidence_points = annotation_points(granule)

    product_dir = Path(product_dir)
    with rasterio.open(products["vv"]) as source:
        transform, width, height = calculate_default_transform(
            source.crs, "EPSG:4326", source.width, source.height,
            *source.bounds, resolution=20 / 111320.0)
    profile = dict(driver="GTiff", width=width, height=height, count=1,
                   crs="EPSG:4326", transform=transform, dtype="float32",
                   nodata=float("nan"), compress="deflate", tiled=True)
    names = {"vv": "Gamma0_VV.tif", "vh": "Gamma0_VH.tif",
             "dem": f"{granule}_DEM.tif", "angle": f"{granule}_localIncidenceAngle.tif"}
    with rasterio.open(products["ls"]) as source:
        ls = np.zeros((height, width), dtype="uint8")  # 0 = not tested
        reproject(source.read(1), ls, src_transform=source.transform, src_crs=source.crs,
                  dst_transform=transform, dst_crs=profile["crs"], src_nodata=0, dst_nodata=0,
                  resampling=Resampling.nearest)
        shadow = (ls & (8 | 16)) > 0
    inc = ellipsoid_incidence(incidence_points, transform, width, height)
    for key, filename in names.items():
        with rasterio.open(products[key]) as source:
            if source.crs is None or source.count != 1:
                raise ValueError(f"Invalid ASF raster: {products[key]}")
            values = source.read(1, masked=True).astype("float32").filled(np.nan)
            if key == "angle":
                values[(values < 0) | (values > np.pi)] = np.nan
                values = np.rad2deg(values)
            destination = np.full((height, width), np.nan, dtype="float32")
            reproject(values, destination, src_transform=source.transform,
                      src_crs=source.crs, src_nodata=np.nan,
                      dst_transform=transform, dst_crs=profile["crs"],
                      dst_nodata=np.nan, resampling=Resampling.bilinear)
        if not np.isfinite(destination).any():
            raise ValueError(f"ASF raster has no valid pixels: {products[key]}")
        if key == "vv":
            valid_vv = np.isfinite(destination)
        with rasterio.open(product_dir / filename, "w", **profile) as target:
            target.write(destination, 1)
            target.update_tags(download_source="asf", processing="HyP3 RTC")
        if key == "angle":
            mask = np.full(destination.shape, 255, dtype="uint8")
            valid = np.isfinite(destination)
            mask[valid] = ((destination[valid] > 50) | shadow[valid]).astype("uint8")
            with rasterio.open(product_dir / f"{granule}_LIA.tif", "w",
                               **dict(profile, dtype="uint8", nodata=255)) as target:
                target.write(mask, 1)
                target.set_band_description(1, "local_incidence_angle_greater_than_50_degrees_or_shadow")
                target.update_tags(mask_semantics="0=keep,1=exclude,255=nodata",
                                   lia_threshold_degrees=50, download_source="asf",
                                   shadow="HyP3 ls_map values 8 or 16 also excluded")
    inc[~valid_vv] = np.nan
    with rasterio.open(product_dir / f"{granule}_incidenceAngleFromEllipsoid.tif", "w", **profile) as target:
        target.write(inc, 1)
        target.set_band_description(1, "ellipsoid_incidence_angle_degrees")
        target.update_tags(download_source="asf", source="S1 annotation geolocationGrid")


RTC_OPTIONS = dict(radiometry="gamma0", scale="power", resolution=20,
                   include_inc_map=True, include_dem=True)


def output_names(granule):
    return ["Gamma0_VV.tif", "Gamma0_VH.tif", f"{granule}_DEM.tif",
            f"{granule}_localIncidenceAngle.tif", f"{granule}_LIA.tif",
            f"{granule}_incidenceAngleFromEllipsoid.tif"]


def outputs_complete(product_dir, granule):
    return all((Path(product_dir) / name).is_file() and (Path(product_dir) / name).stat().st_size > 0
               for name in output_names(granule))


def retry_transient(action, label, attempts=4, wait=30):
    """Retry HyP3 5xx/unavailable and network errors with backoff; other errors fail at once."""
    import time

    import requests
    from hyp3_sdk.exceptions import ServerError, ServiceUnavailableError

    for attempt in range(1, attempts + 1):
        try:
            return action()
        except (ServerError, ServiceUnavailableError, requests.RequestException) as exc:
            if attempt == attempts:
                raise
            delay = wait * 2 ** (attempt - 1)
            print(f"{label}: transient error ({type(exc).__name__}: {exc}); "
                  f"retry {attempt}/{attempts - 1} in {delay} s", flush=True)
            time.sleep(delay)


def find_or_submit(hyp3, granule):
    """Reuse a pending or unexpired successful job with identical options; otherwise submit.

    Each retry searches again first, so a submission that reached HyP3 is never duplicated.
    """
    return retry_transient(lambda: _find_or_submit(hyp3, granule), f"ASF submission for {granule}")


def _find_or_submit(hyp3, granule):
    reusable = []
    for job in hyp3.find_jobs(name=granule, job_type="RTC_GAMMA"):
        params = job.job_parameters or {}
        if (params.get("granules") == [granule] and all(params.get(k) == v for k, v in RTC_OPTIONS.items())
                and not job.failed() and not job.expired()):
            reusable.append(job)
    if reusable:
        job = max(reusable, key=lambda j: j.request_time)
        print(f"Reusing ASF RTC job {job.job_id} ({job.status_code}) for {granule}", flush=True)
        return job
    job = list(hyp3.submit_rtc_job(granule=granule, name=granule, **RTC_OPTIONS))[0]
    print(f"Submitted ASF RTC job {job.job_id} for {granule}", flush=True)
    return job


def download_and_convert(job, granule, product_dir):
    # Keep downloads on failure, and publish outputs only after conversion succeeds.
    staging = Path(tempfile.mkdtemp(prefix="asf_", dir=product_dir))
    retry_transient(lambda: job.download_files(location=str(staging)), f"ASF download for {granule}")
    products = extract_products(sorted(staging.glob("*.zip")), staging)
    normalize_products(products, staging, granule)
    for filename in output_names(granule):
        (staging / filename).replace(Path(product_dir) / filename)
    shutil.rmtree(staging)
    print(f"ASF processing completed: {product_dir}", flush=True)


def download_and_process(granules, output_dir, workers=4, poll_seconds=60, timeout=6 * 3600):
    """Run all HyP3 RTC jobs concurrently; download each as soon as it succeeds.

    A failure is recorded for its granule only. Returns {granule: error message}.
    """
    from concurrent.futures import ThreadPoolExecutor
    import time

    import hyp3_sdk as sdk

    output_dir = Path(output_dir).expanduser().resolve()
    failures = {}
    pending = {}  # job_id -> (granule, job)
    hyp3 = None
    for granule in dict.fromkeys(g.removesuffix(".SAFE") for g in granules):
        try:
            if not granule or Path(granule).name != granule or granule in (".", ".."):
                raise ValueError(f"Invalid granule name: {granule}")
            product_dir = output_dir / granule
            product_dir.mkdir(parents=True, exist_ok=True)
            if outputs_complete(product_dir, granule):
                print(f"ASF outputs already present, skipping: {granule}", flush=True)
                continue
            hyp3 = hyp3 or sdk.HyP3()
            job = find_or_submit(hyp3, granule)
            pending[job.job_id] = (granule, job)
        except Exception as exc:
            failures[granule] = f"{type(exc).__name__}: {exc}"
            print(f"ASF submission failed for {granule}: {failures[granule]}", flush=True)

    futures = {}
    deadline = time.monotonic() + timeout
    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        while pending:
            try:
                refreshed = {job.job_id: job for job in hyp3.refresh(sdk.Batch([j for _, j in pending.values()]))}
            except Exception as exc:  # Transient API errors: keep polling.
                print(f"ASF status check failed, retrying: {exc}", flush=True)
                refreshed = {}
            for job_id, job in refreshed.items():
                granule = pending[job_id][0]
                pending[job_id] = (granule, job)
                if job.succeeded():
                    del pending[job_id]
                    futures[granule] = pool.submit(download_and_convert, job, granule, output_dir / granule)
                elif job.failed():
                    del pending[job_id]
                    failures[granule] = f"ASF RTC job {job_id} failed"
                    print(f"{failures[granule]} for {granule}", flush=True)
            if not pending:
                break
            if time.monotonic() > deadline:
                for job_id, (granule, _) in pending.items():
                    failures[granule] = f"ASF RTC job {job_id} still running after {timeout} s; rerun to resume"
                break
            print(f"ASF jobs: {len(pending)} processing, {len(futures)} downloaded/downloading, "
                  f"{len(failures)} failed", flush=True)
            time.sleep(poll_seconds)
        for granule, future in futures.items():
            try:
                future.result()
            except Exception as exc:
                failures[granule] = f"{type(exc).__name__}: {exc}"
                print(f"ASF download/conversion failed for {granule}: {failures[granule]}", flush=True)
    print(f"ASF batch finished: {len(failures)} failed", flush=True)
    return failures


def main(argv=None):
    parser = argparse.ArgumentParser(description="Download ASF HyP3 RTC products using Earthdata .netrc authentication")
    parser.add_argument("names", help="Comma-separated Sentinel-1 product names")
    parser.add_argument("output_folder")
    parser.add_argument("--workers", type=int, default=4, help="Concurrent downloads/conversions")
    args = parser.parse_args(argv)
    names = [name.strip() for name in args.names.split(",") if name.strip()]
    if not names:
        parser.error("At least one product name is required")
    return 1 if download_and_process(names, args.output_folder, args.workers) else 0


if __name__ == "__main__":
    raise SystemExit(main())
