# Sentinel-1 Download and Processing

This project downloads specified Sentinel-1 GRD products and generates the
SAR and ancillary raster products used by the downstream workflow.

## Processing workflow

For each Sentinel-1 product name, the workflow runs the following steps:

1. `Sentinel_1_specific_name_download_process.py`
   - Searches for and downloads the Sentinel-1 GRD product from Copernicus
     Data Space.
   - Uses SNAP through pyroSAR to process VV and VH backscatter.
   - Produces the incidence-angle rasters required by the next step.
2. `Desert_mask.py`
   - Clips and reprojects the external global desert-mask VRT to the
     Sentinel-1 grid.
   - Preserves the mask classes using nearest-neighbour resampling.
3. `cal_LIA.py`
   - Reuses the matching `*_DEM.tif` exported by SNAP on the Sentinel-1 grid.
   - Checks its CRS, dimensions, transform, and invalid elevations before
     calculating local incidence angle (LIA). A mismatched grid raises an error.
   - If the SNAP DEM is absent, downloads Copernicus DEM tiles from Microsoft
     Planetary Computer and resamples them to the Sentinel-1 grid.
   - Writes a uint8 exclusion mask: 0 for LIA up to 50 degrees, 1 for LIA
     greater than 50 degrees, and 255 for nodata.
4. `Snow_detect.py`
   - Searches for Sentinel-2 L2A products from the 15 days before the
     Sentinel-1 acquisition.
   - Uses the Sentinel-2 Scene Classification Layer (SCL) to calculate snow
     occurrence and persistent cloud on the Sentinel-1 grid.
   - Combines overlapping Sentinel-2 tiles by acquisition day before counting
     valid/cloud observations, avoiding duplicate observations on tile seams.

After all four stages succeed, `execute.sh` removes downloaded ZIP files, DEM
tiles, and SNAP/intermediate rasters. If a stage fails, its intermediate files
are retained for inspection or retry.

## Environment

The current Conda environment was exported to `s1pro.yml`:

```bash
conda env create -f s1pro.yml
conda activate s1pro
```

The workflow also requires ESA SNAP. Make sure the SNAP `bin` directory is
available on `PATH`. The current server-specific SNAP path is set in
`Sentinel_1_specific_name_download_process.py` and may need to be changed on a
different system.

Copernicus Data Space credentials are passed to the download script as command
line arguments. Set the empty `username` and `password` variables in
`execute.sh` before running the workflow.

## Running the workflow

Edit `execute.sh` before submission:

- Set `path` to the output root directory.
- Set `desert_mask_vrt` to the external global desert-mask VRT.
- Set `username` and `password` to your Copernicus Data Space credentials.
- Add one or more Sentinel-1 product names to the `s1names` array.
- Update the Conda, project, SNAP, Slurm partition, and log paths when running
  on a different system.

Submit the workflow with:

```bash
sbatch execute.sh
```

The four stages can also be run manually:

```bash
python Sentinel_1_specific_name_download_process.py \
  <S1_PRODUCT_NAME> <OUTPUT_FOLDER> <USERNAME> <PASSWORD>
python Desert_mask.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER> <DESERT_MASK_VRT>
python cal_LIA.py <S1_PRODUCT_NAME> <PRODUCT_FOLDER> \
  --incidence-angle <ELLIPSOID_INCIDENCE_ANGLE> \
  --metadata-dir <PRODUCT_FOLDER>
python Snow_detect.py \
  <S1_PRODUCT_NAME> <REFERENCE_GAMMA0> <PRODUCT_FOLDER>
```

The legacy two-argument forms remain supported:

```bash
python cal_LIA.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER>
python Snow_detect.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER>
```

They discover inputs under `<OUTPUT_FOLDER>/<S1_PRODUCT_NAME>/`. Both the
legacy and explicit calling forms use the same scene-prefixed output names.

## Search flood-warning SAR images over a date range

Run in the project's Conda environment:

```bash
python Sentinel_1_ESA_search_download_process_chain_v4.py --search-only \
  --start-date 2025-10-01 --end-date 2025-10-07 \
  --work-directory ./data --output-json flood_images.json
```

Both dates are inclusive in UTC. For each day, the search downloads that day's
GloFAS warning shapefile and queries that day's Sentinel-1 `IW_GRDH_1S` products
overlapping original warning pixels inside retained regions by at least
50 km? per scene per region. Bounding boxes serve as catalogue prefilters;
final filtering rasterizes each candidate SAR footprint onto the warning grid
and counts original warning pixels tagged with the current region ID. Expanded
pixels and other regions inside the same bounding box do not contribute.
It follows catalogue
pagination and deduplicates by CDSE product UUID across regions and days.
No SAR downloads, CDSE credentials, or desert-mask VRT are required in this mode.
Warning files and intermediate rasters are saved in the work directory.
Missing warning data or request failures stop the search with the affected date;
they are not silently counted as zero images.

The JSON result contains `names` (product names with the trailing `.SAFE` removed)
and `count` (product count); UUIDs are used internally for deduplication only.
The same result is
available from Python without triggering the processing workflow on import:

```python
from Sentinel_1_ESA_search_download_process_chain_v4 import search_flood_images_by_date_range

result = search_flood_images_by_date_range("2025-10-01", "2025-10-07", "./data")
print(result["names"], result["count"])
```

Defaults match the daily warning workflow: raster size `1/111` degrees,
`--window-size 10` (10 by 10 maximum filter), and `--area-thresholds 1000`
(strictly greater than 1000 square kilometres after expansion). Area is currently
approximated by geometry area in square degrees multiplied by `110 * 110`,
so it is not an equal-area measurement. All input warning geometries are used;
the code applies no warning probability, discharge, return-period or severity
attribute threshold. Matches are candidate scenes over warning areas, not
confirmation that flooding is visible in the SAR imagery.

The overlap threshold is `MIN_FLOOD_OVERLAP_KM2 = 50.0` in the workflow script;
its area conversion remains `geometry.area * 110 * 110`. Both the daily workflow
and date-range search use this filter, with existing function signatures and
commands unchanged. The main workflow returns region features carrying
`warning_raster` and `warning_region` for pixel filtering. Legacy callers without
this metadata retain their existing geometry / `warning_wkt` / bounds fallback.
A missing product footprint raises an error rather than silently accepting an
unverified match.

Step 5 creates a tiled, compressed `*_warning_regions.tif`: 0 denotes background
and positive values identify retained regions, only on the original pre-expansion
warning pixels. It rasterizes the existing expanded region polygons in bounded
512-by-512 tiles; detailed original warning polygons are no longer indexed,
clipped or unioned. Candidate scenes read only local windows from this raster,
in bounded tiles, and can stop once the area threshold is reached. Only stage-level elapsed times are logged, including download/extraction,
shapefile reading, raster processing, and catalogue search with overlap filtering.
Each searched date ends with its unique image count, total elapsed time and
cumulative image count; per-tile and per-region progress is omitted. Intermediate files remain on disk as before.

Overlap is now a raster approximation using `all_touched=False`, consistent with
the original warning rasterization. The existing area convention gives each
`1/111`-degree pixel about 0.982 km?, so 50 warning pixels fall below the threshold
and 51 pass. Narrow features and boundary pixels can differ from exact polygon
intersection; no latitude correction or threshold change is introduced.

## Output structure

Each Sentinel-1 product is written to its own subdirectory:

```text
<OUTPUT_FOLDER>/
+-- <S1_PRODUCT_NAME>/
    +-- Gamma0_VV.tif
    +-- Gamma0_VH.tif
    +-- <S1_PRODUCT_NAME>_desert.tif
    +-- <S1_PRODUCT_NAME>_LIA.tif
    +-- <S1_PRODUCT_NAME>_ice.tif
    +-- <S1_PRODUCT_NAME>_cloud.tif
```

The `*_DEM.tif` exported by SNAP is retained after successful completion for
reuse and inspection. Its elevations are used as exported, without applying
another geoid correction. Other intermediate SNAP output files are removed.

## Download channel

All existing commands continue to use Copernicus Data Space + SNAP by default.
Use `--download-source cdse` to select that path explicitly, or
`--download-source asf` for ASF HyP3 RTC. Catalogue searching remains on CDSE
and does not require authentication. ASF authentication uses Earthdata
`~/.netrc`, independently of the positional CDSE credentials.

```bash
# Existing invocation, unchanged:
python Sentinel_1_specific_name_download_process.py "$SCENE" ./data "$CDSE_USERNAME" "$CDSE_PASSWORD"
# ASF does not require the two CDSE credential arguments:
python Sentinel_1_specific_name_download_process.py "$SCENE" ./data --download-source asf
# Direct ASF adapter (also accepts comma-separated scene names):
python Sentinel_1_specific_name_ASF.py "$SCENE" ./data
# Daily warning workflow:
python Sentinel_1_ESA_search_download_process_chain_v4.py --desert-mask-vrt /path/desert.vrt --download-source asf
# Shell/Slurm entry points (default DOWNLOAD_SOURCE=cdse):
DOWNLOAD_SOURCE=asf bash execute.sh
sbatch --export=ALL,DOWNLOAD_SOURCE=asf execute_flood_warning.sh
```

Install the optional ASF dependency into the processing environment with
`python -m pip install hyp3-sdk`. The existing SNAP environment remains usable
without it. Configure `~/.netrc` with your Earthdata account:

```text
machine urs.earthdata.nasa.gov
  login YOUR_EARTHDATA_USERNAME
  password YOUR_EARTHDATA_PASSWORD
```

On Linux use `chmod 600 ~/.netrc`. ASF mode submits HyP3 RTC jobs using the
account's available processing credits and waits for completion; each invocation
submits new jobs. Failed jobs or incomplete VV/VH/angle/DEM packages stop that
scene. Download/conversion intermediates are retained on failure.

Both channels produce `Gamma0_VV.tif`, `Gamma0_VH.tif` and the same scene-prefixed
`_desert.tif`, `_LIA.tif`, `_ice.tif`, `_cloud.tif` interface in the scene directory.
ASF requests 20 m linear gamma0 RTC and reprojects layers onto a shared WGS84
(EPSG:4326) grid at 20/111320 degrees. Gamma0 and DEM are float32 with NaN nodata.
The DEM is retained as `<SCENE>_DEM.tif`. The ancillary mask workflow uses this
grid. Exact grid extents and pixel values are not guaranteed to match SNAP:
HyP3 and SNAP use different processing algorithms, and the existing SNAP path
also applies its original ellipsoid-angle normalization, which is preserved.
ASF gamma0 is used directly without an additional cosine correction.

The [HyP3 product guide](https://hyp3-docs.asf.alaska.edu/guides/rtc_product_guide/)
defines `_inc_map.tif` as **local incidence angle in radians**. The adapter converts
it to degrees and writes the LIA mask directly: 0 for <=50 degrees, 1 for >50,
255 for nodata. It must not be passed to `cal_LIA.py` as an ellipsoid angle.
The ASF branch therefore skips that calculation; desert, snow and cloud stages
still run. Successful workflow cleanup removes the intermediate local-angle
raster, while retaining the DEM. The standalone ASF adapter produces the SAR,
DEM, angle and LIA products; the full chain adds the other masks.
