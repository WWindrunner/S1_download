"""Align an external categorical desert mask to a SAR reference grid."""
import argparse
from pathlib import Path
import sys
import time

import numpy as np
import rasterio
from rasterio.vrt import WarpedVRT
from rasterio.warp import Resampling

from nisar_common import raster_profile, validate_raster


def generate_desert_mask(product_name, reference_path, output_dir, desert_mask_vrt, strict=False):
    started = time.perf_counter()
    output = Path(output_dir) / f"{product_name}_desert.tif"
    output.parent.mkdir(parents=True, exist_ok=True)
    partial = output.with_name(output.name + ".part")
    with rasterio.open(reference_path) as reference:
        if reference.crs is None or reference.transform.b != 0 or reference.transform.d != 0:
            raise ValueError("Reference must have a CRS and an unrotated grid")
        try:
            with rasterio.open(desert_mask_vrt) as source:
                if source.count != 1 or source.crs is None:
                    raise ValueError("Desert source must be single-band with a CRS")
                nodata = source.nodata if source.nodata is not None else -9999
                profile = raster_profile(reference, nodata=nodata)
                with WarpedVRT(source, crs=reference.crs, transform=reference.transform,
                               width=reference.width, height=reference.height, dtype="float32",
                               nodata=nodata, resampling=Resampling.nearest) as aligned, rasterio.open(
                                   partial, "w", **profile) as dst:
                    for _, window in dst.block_windows(1):
                        values = aligned.read(1, window=window, masked=True).filled(nodata)
                        valid = np.isfinite(reference.read(1, window=window, masked=True).filled(np.nan))
                        values[~valid] = nodata
                        dst.write(values, 1, window=window)
        except Exception:
            if strict:
                raise
            print("Desert mask failed; creating the legacy all-nodata fallback.", file=sys.stderr)
            with rasterio.open(partial, "w", **raster_profile(reference, nodata=-9999)) as dst:
                for _, window in dst.block_windows(1):
                    dst.write(np.full((int(window.height), int(window.width)), -9999, dtype="float32"), 1, window=window)
    validate_raster(partial, reference_path, require_valid=False)
    partial.replace(output)
    print(f"Desert mask completed in {time.perf_counter() - started:.2f} seconds.")
    return str(output)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("product_name")
    parser.add_argument("output_folder")
    parser.add_argument("desert_mask_vrt")
    parser.add_argument("--reference-raster", help="Default: <OUTPUT_FOLDER>/<PRODUCT>/Gamma0_VV.tif")
    parser.add_argument("--strict", action="store_true", help="Stop on source errors instead of writing a nodata fallback")
    args = parser.parse_args(argv)
    directory = Path(args.output_folder).expanduser().resolve() / args.product_name
    reference = args.reference_raster or directory / "Gamma0_VV.tif"
    generate_desert_mask(args.product_name, reference, directory, args.desert_mask_vrt, args.strict)


if __name__ == "__main__":
    main()
