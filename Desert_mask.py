import os
import sys
import traceback

from osgeo import gdal


gdal.UseExceptions()


if len(sys.argv) != 4:
    print(
        "Usage: python Desert_mask.py "
        "<S1_PRODUCT_NAME> <OUTPUT_FOLDER> <DESERT_MASK_VRT>"
    )
    sys.exit(1)

S1name = sys.argv[1]
folder = os.path.abspath(os.path.expanduser(sys.argv[2]))
desert_mask_vrt = os.path.abspath(os.path.expanduser(sys.argv[3]))
product_dir = os.path.join(folder, S1name)

if not os.path.isdir(product_dir):
    raise FileNotFoundError(f"S1 product folder not found: {product_dir}")

reference_path = os.path.join(product_dir, "Gamma0_VV.tif")
if not os.path.isfile(reference_path):
    raise FileNotFoundError(f"Reference raster not found: {reference_path}")

output_path = os.path.join(product_dir, f"{S1name}_desert.tif")


def create_nodata_desert_raster(reference_ds):
    driver = gdal.GetDriverByName("GTiff")
    output_ds = driver.Create(
        output_path,
        reference_ds.RasterXSize,
        reference_ds.RasterYSize,
        1,
        gdal.GDT_Float32,
        options=["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"],
    )
    if output_ds is None:
        raise RuntimeError(f"Failed to create empty desert mask: {output_path}")
    output_ds.SetGeoTransform(reference_ds.GetGeoTransform())
    output_ds.SetProjection(reference_ds.GetProjection())
    band = output_ds.GetRasterBand(1)
    band.SetNoDataValue(-9999)
    band.Fill(-9999)
    band.FlushCache()
    output_ds = None


# The reference raster defines the output grid and fallback raster.
reference_ds = gdal.Open(reference_path)
if reference_ds is None:
    raise RuntimeError(f"Failed to open reference raster: {reference_path}")

try:
    if not os.path.isfile(desert_mask_vrt):
        raise FileNotFoundError(f"Desert mask VRT not found: {desert_mask_vrt}")

    reference_gt = reference_ds.GetGeoTransform()
    reference_projection = reference_ds.GetProjection()
    reference_cols = reference_ds.RasterXSize
    reference_rows = reference_ds.RasterYSize

    if not reference_projection:
        raise ValueError(f"Reference raster has no CRS: {reference_path}")

    if reference_gt[2] != 0 or reference_gt[4] != 0:
        raise ValueError("Rotated reference rasters are not supported")

    x_end = reference_gt[0] + reference_cols * reference_gt[1]
    y_end = reference_gt[3] + reference_rows * reference_gt[5]
    reference_bounds = [
        min(reference_gt[0], x_end),
        min(reference_gt[3], y_end),
        max(reference_gt[0], x_end),
        max(reference_gt[3], y_end),
    ]

    desert_ds = gdal.Open(desert_mask_vrt)
    if desert_ds is None:
        raise RuntimeError(f"Failed to open desert mask: {desert_mask_vrt}")
    if desert_ds.RasterCount != 1:
        raise ValueError(
            "Desert mask must contain exactly one band, "
            f"found {desert_ds.RasterCount}"
        )

    source_nodata = desert_ds.GetRasterBand(1).GetNoDataValue()
    warp_kwargs = {
        "format": "GTiff",
        "dstSRS": reference_projection,
        "outputBounds": reference_bounds,
        "width": reference_cols,
        "height": reference_rows,
        "resampleAlg": gdal.GRA_NearestNeighbour,
        "multithread": True,
        "creationOptions": ["COMPRESS=LZW", "TILED=YES", "BIGTIFF=IF_SAFER"],
    }
    if source_nodata is not None:
        warp_kwargs["srcNodata"] = source_nodata
        warp_kwargs["dstNodata"] = source_nodata

    print(f"Reference raster: {reference_path}")
    print(f"Desert mask: {desert_mask_vrt}")
    print(f"Output raster: {output_path}")

    output_ds = gdal.Warp(
        output_path,
        desert_ds,
        options=gdal.WarpOptions(**warp_kwargs),
    )
    if output_ds is None:
        raise RuntimeError(f"Failed to create desert mask: {output_path}")

    output_ds.FlushCache()
    output_ds = None
    desert_ds = None
    print(f"Desert mask saved -> {output_path}")
except Exception:
    print(
        "ERROR: Desert mask calculation failed. Creating an empty nodata mask ",
        "so processing can continue.",
        file=sys.stderr,
    )
    traceback.print_exc()
    create_nodata_desert_raster(reference_ds)
    print(f"Empty desert mask saved -> {output_path}")
finally:
    reference_ds = None
