"""Interpolate NISAR radar metadata at DEM heights and calculate terrain LIA."""
import argparse
from contextlib import ExitStack
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window
from scipy.interpolate import RegularGridInterpolator

from nisar_common import gcov_root, raster_profile, read_crs, read_values, text_value, validate_raster
from sar_dem import prepare_dem


def cube_interpolators(root):
    radar = root["metadata/radarGrid"]
    axes = tuple(np.asarray(radar[key][:], dtype=float) for key in
                 ("heightAboveEllipsoid", "yCoordinates", "xCoordinates"))
    for axis in axes:
        if (axis.ndim != 1 or len(axis) < 2 or not np.isfinite(axis).all()
                or not (np.all(np.diff(axis) > 0) or np.all(np.diff(axis) < 0))):
            raise ValueError("Radar cube axes must be finite and strictly monotonic")
    fields = {}
    for key in ("incidenceAngle", "losUnitVectorX", "losUnitVectorY"):
        if key not in radar:
            raise ValueError(f"NISAR LIA needs radarGrid/{key}; an incidence angle alone is insufficient")
        dataset = radar[key]
        values = read_values(dataset)
        if values.shape != tuple(len(a) for a in axes):
            raise ValueError(f"Unexpected radar cube shape for {key}: {values.shape}")
        if key == "incidenceAngle":
            units = text_value(dataset.attrs.get("units", "degrees")).lower()
            if "rad" in units:
                values = np.rad2deg(values)
            elif "deg" not in units:
                raise ValueError(f"Unknown incidence-angle units: {units}")
            values[(values < 0) | (values > 90)] = np.nan
        else:
            values[np.abs(values) > 1] = np.nan
        fields[key] = RegularGridInterpolator(axes, values, bounds_error=False, fill_value=np.nan)
    # The metadata cube has its own projection, independent of the image grid.
    crs = read_crs(radar)
    return fields, crs


def local_incidence(dem, transform, incidence, east, north):
    """Angle between terrain normal and ground-to-sensor LOS, on a WGS84 grid."""
    if min(dem.shape) < 2 or transform.b != 0 or transform.d != 0 or transform.e >= 0:
        raise ValueError("LIA needs a north-up grid with at least 2 rows and columns")
    lat = transform.f + (np.arange(dem.shape[0]) + 0.5) * transform.e
    lat = np.deg2rad(lat)
    metres_lat = 111132.92 - 559.82 * np.cos(2 * lat) + 1.175 * np.cos(4 * lat) - 0.0023 * np.cos(6 * lat)
    metres_lon = 111412.84 * np.cos(lat) - 93.5 * np.cos(3 * lat) + 0.118 * np.cos(5 * lat)
    dx = transform.a * metres_lon[:, None]
    dy = transform.e * metres_lat[:, None]  # negative: raster rows go south
    dzde = np.gradient(dem, axis=1) / dx
    dzdn = np.gradient(dem, axis=0) / dy
    horizontal = np.hypot(east, north)
    angle = np.deg2rad(incidence)
    # Use the interpolated incidence for inclination and LOS components for azimuth.
    scale = np.divide(np.sin(angle), horizontal, out=np.full_like(angle, np.nan), where=horizontal > 1e-8)
    lx, ly, lz = east * scale, north * scale, np.cos(angle)
    overhead = (horizontal <= 1e-8) & (np.abs(incidence) < 1e-5)
    lx[overhead], ly[overhead] = 0, 0
    cosine = (-dzde * lx - dzdn * ly + lz) / np.sqrt(1 + dzde**2 + dzdn**2)
    return np.rad2deg(np.arccos(np.clip(cosine, -1, 1))).astype(np.float32)


def interpolate_angles(h5_file, reference_path, dem_path, output_dir, scene,
                       band="LSAR", threshold=50.0):
    import h5py
    from pyproj import Transformer

    if not np.isfinite(threshold) or not 0 <= threshold <= 180:
        raise ValueError("LIA threshold must be between 0 and 180 degrees")
    validate_raster(dem_path, reference_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {"incidence": output_dir / f"{scene}_incidenceAngle.tif",
               "local_angle": output_dir / f"{scene}_localIncidenceAngle.tif",
               "lia": output_dir / f"{scene}_LIA.tif"}
    partials = {key: path.with_name(path.name + ".part") for key, path in outputs.items()}
    with ExitStack() as stack:
        handle = stack.enter_context(h5py.File(h5_file, "r"))
        fields, cube_crs = cube_interpolators(gcov_root(handle, band))
        reference = stack.enter_context(rasterio.open(reference_path))
        dem = stack.enter_context(rasterio.open(dem_path))
        if reference.crs.to_epsg() != 4326:
            raise ValueError("NISAR angle processing currently expects the converter's EPSG:4326 grid")
        if dem.tags().get("height_reference") != "WGS84_ellipsoid":
            raise ValueError("DEM must carry height_reference=WGS84_ellipsoid; prepare it before interpolation")
        transformer = Transformer.from_crs(reference.crs, cube_crs, always_xy=True)
        writers = {key: stack.enter_context(rasterio.open(path, "w", **raster_profile(
            reference, **(dict(dtype="uint8", nodata=255) if key == "lia" else {}))))
            for key, path in partials.items()}
        for _, window in reference.block_windows(1):
            # One-pixel halo keeps slope estimates continuous across tile boundaries.
            row0 = max(0, int(window.row_off) - 1)
            col0 = max(0, int(window.col_off) - 1)
            row1 = min(reference.height, int(window.row_off + window.height) + 1)
            col1 = min(reference.width, int(window.col_off + window.width) + 1)
            halo = Window(col0, row0, col1 - col0, row1 - row0)
            heights = dem.read(1, window=halo, masked=True).filled(np.nan)
            affine = reference.window_transform(halo)
            xx, yy = np.meshgrid(affine.c + (np.arange(heights.shape[1]) + 0.5) * affine.a,
                                 affine.f + (np.arange(heights.shape[0]) + 0.5) * affine.e)
            x, y = transformer.transform(xx, yy)
            points = np.column_stack((heights.ravel(), np.asarray(y).ravel(), np.asarray(x).ravel()))
            values = {key: field(points).reshape(heights.shape).astype(np.float32)
                      for key, field in fields.items()}
            angle = values["incidenceAngle"]
            lia = local_incidence(heights, affine, angle, values["losUnitVectorX"], values["losUnitVectorY"])
            crop = (slice(int(window.row_off) - row0, int(window.row_off + window.height) - row0),
                    slice(int(window.col_off) - col0, int(window.col_off + window.width) - col0))
            angle, lia = angle[crop].copy(), lia[crop].copy()
            valid_sar = np.isfinite(reference.read(1, window=window, masked=True).filled(np.nan))
            angle[~valid_sar], lia[~valid_sar] = np.nan, np.nan
            mask = np.full(angle.shape, 255, dtype=np.uint8)
            valid = np.isfinite(lia)
            mask[valid] = (lia[valid] > threshold).astype(np.uint8)
            for key, array in (("incidence", angle), ("local_angle", lia), ("lia", mask)):
                writers[key].write(array, 1, window=window)
        writers["incidence"].update_tags(units="degrees", angle_type="ellipsoid_incidence",
                                         interpolation="linear_xyz_no_extrapolation")
        writers["local_angle"].update_tags(units="degrees", angle_type="local_incidence",
                                           method="DEM_surface_normal_and_NISAR_LOS")
        writers["lia"].update_tags(mask_semantics="0=keep,1=exclude,255=nodata",
                                   lia_threshold_degrees=threshold)
    for key, partial in partials.items():
        validate_raster(partial, reference_path)
    for key, path in outputs.items():
        partials[key].replace(path)
    return {key: str(path) for key, path in outputs.items()}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5_file")
    parser.add_argument("reference_path")
    parser.add_argument("output_dir")
    parser.add_argument("--band", choices=("LSAR", "SSAR"), default="LSAR")
    parser.add_argument("--dem", help="Optional local DEM, otherwise download Copernicus GLO-30")
    parser.add_argument("--dem-height-reference", choices=("egm2008", "ellipsoid"), default="egm2008")
    parser.add_argument("--geoid", help="EGM2008 undulation raster; otherwise cached PROJ grid is downloaded")
    parser.add_argument("--cache-dir")
    parser.add_argument("--lia-threshold", type=float, default=50.0)
    args = parser.parse_args(argv)
    scene = Path(args.h5_file).stem
    dem = prepare_dem(args.reference_path, Path(args.output_dir) / f"{scene}_DEM.tif",
                      args.cache_dir or Path(args.output_dir) / ".cache", args.dem,
                      args.dem_height_reference, args.geoid)
    interpolate_angles(args.h5_file, args.reference_path, dem, args.output_dir, scene,
                       args.band, args.lia_threshold)


if __name__ == "__main__":
    main()
