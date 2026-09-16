"""ASF HyP3 RTC adapter for the maintained Sentinel-1 processing chain.

Authentication uses Earthdata credentials in ~/.netrc. Importing this module
never submits jobs. HyP3 local incidence angles are radians, not ellipsoid angles.
"""
import argparse
from pathlib import Path
import shutil
import tempfile
import zipfile


def extract_products(archives, directory):
    """Extract only recognized rasters; reject missing or ambiguous layers."""
    suffixes = {"vv": "_vv.tif", "vh": "_vh.tif",
                "angle": "_inc_map.tif", "dem": "_dem.tif"}
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


def normalize_products(products, product_dir, granule):
    """Write float32 gamma0/DEM and uint8 LIA on one WGS84 grid."""
    import numpy as np
    import rasterio
    from rasterio.warp import calculate_default_transform, reproject, Resampling

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
        with rasterio.open(product_dir / filename, "w", **profile) as target:
            target.write(destination, 1)
            target.update_tags(download_source="asf", processing="HyP3 RTC")
        if key == "angle":
            mask = np.full(destination.shape, 255, dtype="uint8")
            valid = np.isfinite(destination)
            mask[valid] = (destination[valid] > 50).astype("uint8")
            with rasterio.open(product_dir / f"{granule}_LIA.tif", "w",
                               **dict(profile, dtype="uint8", nodata=255)) as target:
                target.write(mask, 1)
                target.set_band_description(1, "local_incidence_angle_greater_than_50_degrees")
                target.update_tags(mask_semantics="0=keep,1=exclude,255=nodata",
                                   lia_threshold_degrees=50, download_source="asf")


def download_and_process(granules, output_dir):
    import hyp3_sdk as sdk

    hyp3 = sdk.HyP3()
    for granule in granules:
        granule = granule.removesuffix(".SAFE")
        if not granule or Path(granule).name != granule or granule in (".", ".."):
            raise ValueError(f"Invalid granule name: {granule}")
        product_dir = Path(output_dir).expanduser().resolve() / granule
        product_dir.mkdir(parents=True, exist_ok=True)
        jobs = hyp3.submit_rtc_job(
            granule=granule, name=granule, radiometry="gamma0", scale="power",
            resolution=20, include_inc_map=True, include_dem=True)
        jobs = hyp3.watch(jobs)
        completed = list(jobs)
        if len(completed) != 1 or not completed[0].succeeded:
            raise RuntimeError(f"ASF RTC job did not succeed for {granule}")
        # Keep downloads on failure, and publish outputs only after conversion succeeds.
        staging = Path(tempfile.mkdtemp(prefix="asf_", dir=product_dir))
        completed[0].download_files(location=str(staging))
        products = extract_products(sorted(staging.glob("*.zip")), staging)
        normalize_products(products, staging, granule)
        for filename in ("Gamma0_VV.tif", "Gamma0_VH.tif", f"{granule}_DEM.tif",
                         f"{granule}_localIncidenceAngle.tif", f"{granule}_LIA.tif"):
            (staging / filename).replace(product_dir / filename)
        shutil.rmtree(staging)
        print(f"ASF processing completed: {product_dir}")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Download ASF HyP3 RTC products using Earthdata .netrc authentication")
    parser.add_argument("names", help="Comma-separated Sentinel-1 product names")
    parser.add_argument("output_folder")
    args = parser.parse_args(argv)
    names = [name.strip() for name in args.names.split(",") if name.strip()]
    if not names:
        parser.error("At least one product name is required")
    download_and_process(names, args.output_folder)


if __name__ == "__main__":
    main()
