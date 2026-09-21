"""Download/cache Copernicus DEM and align ellipsoidal heights to a SAR grid."""
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.warp import Resampling, transform_bounds

from nisar_common import download_http, raster_profile, validate_raster

GEOID_URL = "https://cdn.proj.org/us_nga_egm08_25.tif"


def dem_tiles(reference_path, cache_dir):
    import planetary_computer as pc
    from pystac_client import Client

    with rasterio.open(reference_path) as ref:
        bbox = transform_bounds(ref.crs, "EPSG:4326", *ref.bounds)
    catalog = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")
    items = sorted(catalog.search(collections=["cop-dem-glo-30"], bbox=bbox).items(), key=lambda item: item.id)
    if not items:
        raise ValueError("No Copernicus DEM tiles found")
    paths = []
    for item in items:
        if Path(item.id).name != item.id or any(c in item.id for c in '/\\:'):
            raise ValueError("Invalid DEM tile ID")
        path = Path(cache_dir) / "cop-dem-glo-30" / (item.id + ".tif")
        if path.exists():
            try:
                validate_raster(path)
            except (OSError, ValueError):
                download_http(pc.sign(item.assets["data"].href), path)
        else:
            download_http(pc.sign(item.assets["data"].href), path)
        validate_raster(path)
        paths.append(path)
    return paths


def prepare_dem(reference_path, output_path, cache_dir, dem_path=None,
                height_reference="egm2008", geoid_path=None):
    """Write float32 WGS84 ellipsoidal heights: h = H(EGM2008) + N."""
    if height_reference not in ("egm2008", "ellipsoid"):
        raise ValueError("DEM height reference must be egm2008 or ellipsoid")
    if dem_path is None and height_reference != "egm2008":
        raise ValueError("Downloaded Copernicus DEM uses EGM2008; --ellipsoid requires a local DEM")
    sources = [Path(dem_path)] if dem_path else dem_tiles(reference_path, cache_dir)
    if height_reference == "egm2008":
        if geoid_path is None:
            geoid_path = Path(cache_dir) / "us_nga_egm08_25.tif"
            if not geoid_path.exists():
                download_http(GEOID_URL, geoid_path)
        validate_raster(geoid_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    partial = output_path.with_name(output_path.name + ".part")
    with ExitStack() as stack:
        ref = stack.enter_context(rasterio.open(reference_path))
        profile = raster_profile(ref)
        vrt_args = dict(crs=ref.crs, transform=ref.transform, width=ref.width, height=ref.height,
                        dtype="float32", nodata=np.nan, resampling=Resampling.bilinear,
                        warp_mem_limit=128)
        vrts = []
        for path in sources:
            src = stack.enter_context(rasterio.open(path))
            if src.crs is None or src.count != 1:
                raise ValueError(f"Invalid DEM: {path}")
            vrts.append(stack.enter_context(WarpedVRT(src, **vrt_args)))
        geoid = None
        if height_reference == "egm2008":
            src = stack.enter_context(rasterio.open(geoid_path))
            geoid = stack.enter_context(WarpedVRT(src, **vrt_args))
        dst = stack.enter_context(rasterio.open(partial, "w", **profile))
        for _, window in dst.block_windows(1):
            values = np.full((int(window.height), int(window.width)), np.nan, dtype=np.float32)
            for vrt in vrts:
                tile = vrt.read(1, window=window, masked=True).filled(np.nan)
                take = ~np.isfinite(values) & np.isfinite(tile)
                values[take] = tile[take]
            if geoid is not None:
                values += geoid.read(1, window=window, masked=True).filled(np.nan)
            dst.write(values, 1, window=window)
        dst.update_tags(height_reference="WGS84_ellipsoid", units="metres",
                        source_height_reference=height_reference,
                        geoid_model="EGM2008" if geoid is not None else "none")
    validate_raster(partial, reference_path)
    partial.replace(output_path)
    return str(output_path)
