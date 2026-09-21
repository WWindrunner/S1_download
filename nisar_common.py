"""Shared I/O and product metadata for the NISAR GCOV workflow."""
import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np
import rasterio


def local_path(value):
    """Absolute path; allow long NISAR product names on Windows as well as Linux."""
    path = str(Path(value).expanduser().resolve())
    if os.name == "nt" and not path.startswith("\\\\?\\"):
        path = "\\\\?\\UNC\\" + path[2:] if path.startswith("\\\\") else "\\\\?\\" + path
    return Path(path)


def scene_name(value):
    name = str(value).removesuffix(".h5")
    if (not name or Path(name).name != name or any(c in name for c in '/\\:*?')
            or not name.startswith("NISAR_") or "_GCOV_" not in name):
        raise ValueError("Expected an exact NISAR GCOV product name (optionally ending in .h5)")
    return name


def text_value(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def gcov_root(handle, band="LSAR"):
    path = f"science/{band}/GCOV"
    if path not in handle:
        raise ValueError(f"Missing {path}; this workflow requires a {band} GCOV product, not GUNW")
    return handle[path]


def read_crs(group):
    projection = group["projection"]
    try:
        return rasterio.crs.CRS.from_epsg(int(projection[()]))
    except (TypeError, ValueError):
        for key in ("spatial_ref", "crs_wkt"):
            if key in projection.attrs:
                return rasterio.crs.CRS.from_wkt(text_value(projection.attrs[key]))
        raise ValueError(f"Cannot read projection in {group.name}")


def read_values(dataset, selection=()):
    values = np.asarray(dataset[selection], dtype=np.float32)
    # Ignore HDF5's implicit zero fill: only explicit product attributes apply.
    for key in ("_FillValue", "missing_value"):
        if key in dataset.attrs:
            for value in np.asarray(dataset.attrs[key]).ravel():
                values[values == value] = np.nan
    values[~np.isfinite(values)] = np.nan
    return values


def inspect_h5(path, band="LSAR", frequency="A"):
    import h5py

    with h5py.File(path, "r") as handle:
        root = gcov_root(handle, band)
        grid = root[f"grids/frequency{frequency}"]
        pols = [p for p in ("HH", "HV", "VV", "VH", "RH", "RV") if p + p in grid]
        if not pols:
            raise ValueError(f"No diagonal backscatter terms in {grid.name}")
        for pol in pols:
            if grid[pol + pol].shape != (len(grid["yCoordinates"]), len(grid["xCoordinates"])):
                raise ValueError(f"Invalid dimensions for {pol}")
        identification = handle.get(f"science/{band}/identification")
        metadata = {}
        for key in ("zeroDopplerStartTime", "orbitPassDirection", "productType"):
            if identification is not None and key in identification:
                metadata[key] = text_value(identification[key][()])
        if metadata.get("productType", "GCOV") != "GCOV":
            raise ValueError("Only GCOV backscatter products are supported")
        if "zeroDopplerStartTime" not in metadata:
            match = re.search(r"_(\d{8})T(\d{6})_", Path(path).name)
            if match:
                from datetime import datetime
                metadata["zeroDopplerStartTime"] = datetime.strptime(
                    "".join(match.groups()), "%Y%m%d%H%M%S").isoformat() + "Z"
        metadata.update(sensor="nisar", product_type="GCOV", band=band,
                        frequency=frequency, polarizations=pols,
                        radiometry="gamma0", scale="linear_power",
                        source_crs=read_crs(grid).to_string(),
                        source_pixel_spacing=[float(grid["xCoordinateSpacing"][()]),
                                              float(grid["yCoordinateSpacing"][()])])
        return metadata


def raster_profile(reference, **overrides):
    profile = dict(driver="GTiff", width=reference.width, height=reference.height,
                   transform=reference.transform, crs=reference.crs, count=1,
                   dtype="float32", nodata=np.nan, tiled=True, blockxsize=512,
                   blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
    profile.update(overrides)
    return profile


def validate_raster(path, reference=None, require_valid=True):
    with rasterio.open(path) as src:
        if src.count != 1 or src.crs is None:
            raise ValueError(f"Invalid single-band georeferenced raster: {path}")
        if reference is not None:
            with rasterio.open(reference) as ref:
                if (src.crs != ref.crs or src.shape != ref.shape
                        or not src.transform.almost_equals(ref.transform)):
                    raise ValueError(f"Grid mismatch: {path}")
        valid = False
        for _, window in src.block_windows(1):
            values = src.read(1, window=window, masked=True)
            valid |= bool(np.any(np.isfinite(values.compressed())))
        if require_valid and not valid:
            raise ValueError(f"No valid pixels: {path}")


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".part")
    partial.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    partial.replace(path)


def file_stamp(path):
    path = Path(path).resolve()
    return [str(path), path.stat().st_size, path.stat().st_mtime_ns]


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def download_http(url, destination, session=None, expected_size=None):
    """Publish only a complete response; retain .part on failure for inspection."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    own_session = session is None
    session = session or requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=Retry(
        total=3, backoff_factor=1, status_forcelist=(429, 500, 502, 503, 504))))
    try:
        with session.get(url, stream=True, timeout=(30, 300)) as response:
            response.raise_for_status()
            with partial.open("wb") as output:
                for chunk in response.iter_content(1024 * 1024):
                    if chunk:
                        output.write(chunk)
            size = partial.stat().st_size
            length = response.headers.get("Content-Length")
            if not size or (length and not response.headers.get("Content-Encoding") and size != int(length)):
                raise IOError("Incomplete download (HTTP size mismatch)")
            if expected_size is not None and size != int(expected_size):
                raise IOError("Incomplete download (catalogue size mismatch)")
        partial.replace(destination)
    finally:
        if own_session:
            session.close()
    return destination
