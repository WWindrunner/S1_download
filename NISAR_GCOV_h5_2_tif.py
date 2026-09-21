"""Export NISAR GCOV diagonal covariance terms as linear gamma0 GeoTIFFs."""
import argparse
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine, array_bounds
from rasterio.warp import calculate_default_transform, reproject, Resampling

from nisar_common import gcov_root, inspect_h5, read_crs, read_values, validate_raster


def parse_resolution(value):
    if value is None or str(value).lower() == "auto":
        return None
    resolution = float(value)
    if not np.isfinite(resolution) or resolution <= 0:
        raise ValueError("Resolution must be auto or a positive number of metres")
    return resolution


def source_grid(group):
    x = group["xCoordinates"][:]
    y = group["yCoordinates"][:]
    dx = float(group["xCoordinateSpacing"][()])
    dy = float(group["yCoordinateSpacing"][()])
    if (dx <= 0 or dy >= 0 or not np.allclose(np.diff(x), dx)
            or not np.allclose(np.diff(y), dy)):
        raise ValueError("Expected a regular north-up GCOV grid")
    # NISAR coordinates locate pixel centres, affine origins locate outer corners.
    return Affine(dx, 0, x[0] - dx / 2, 0, dy, y[0] - dy / 2), len(x), len(y)


def convert_gcov(h5_file, output_dir, frequency="A", band="LSAR", resolution="auto", polarizations=None):
    import h5py

    metadata = inspect_h5(h5_file, band, frequency)
    pols = polarizations or metadata["polarizations"]
    missing = set(pols) - set(metadata["polarizations"])
    if missing:
        raise ValueError(f"Requested polarizations absent from product: {sorted(missing)}")
    spacing = parse_resolution(resolution)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    with h5py.File(h5_file, "r") as handle:
        grid = gcov_root(handle, band)[f"grids/frequency{frequency}"]
        transform, width, height = source_grid(grid)
        crs = read_crs(grid)
        options = {} if spacing is None else {"resolution": spacing / 111320.0}
        target_transform, target_width, target_height = calculate_default_transform(
            crs, "EPSG:4326", width, height, *array_bounds(height, width, transform), **options)
        profile = dict(driver="GTiff", width=width, height=height, count=1, crs=crs,
                       transform=transform, dtype="float32", nodata=np.nan,
                       tiled=True, blockxsize=512, blockysize=512, compress="deflate", BIGTIFF="IF_SAFER")
        for pol in pols:
            native = output_dir / f".{pol}_native.tif"
            with rasterio.open(native, "w", **profile) as dst:
                for _, window in dst.block_windows(1):
                    slices = window.toslices()
                    values = read_values(grid[pol + pol], slices)
                    values[values <= 0] = np.nan
                    if "mask" in grid:
                        mask = grid["mask"][slices]
                        values[(mask == 0) | (mask == 255)] = np.nan
                    dst.write(values, 1, window=window)
            output = output_dir / f"Gamma0_{pol}.tif"
            partial = output.with_name(output.name + ".part")
            target_profile = dict(profile, crs="EPSG:4326", transform=target_transform,
                                  width=target_width, height=target_height)
            with rasterio.open(native) as src, rasterio.open(partial, "w", **target_profile) as dst:
                reproject(rasterio.band(src, 1), rasterio.band(dst, 1),
                          resampling=Resampling.bilinear, warp_mem_limit=256)
                dst.update_tags(sensor="NISAR", product_type="GCOV", frequency=frequency,
                                band=band, polarization=pol, radiometry="gamma0", scale="linear_power",
                                resolution_mode="auto" if spacing is None else str(spacing))
            validate_raster(partial)
            partial.replace(output)
            native.unlink()
            outputs[pol] = str(output)
            print(f"GCOV {pol}: {output}")
    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("h5_file")
    parser.add_argument("output_dir")
    parser.add_argument("--frequency", choices=("A", "B"), default="A")
    parser.add_argument("--band", choices=("LSAR", "SSAR"), default="LSAR")
    parser.add_argument("--resolution", default="auto", help="auto (source-derived WGS84 grid), or nominal metres")
    parser.add_argument("--polarizations", nargs="+")
    convert_gcov(**vars(parser.parse_args(argv)))


if __name__ == "__main__":
    main()
