# Sentinel-1 and NISAR Download and Processing

This project downloads specified Sentinel-1 GRD or NISAR GCOV products and
generates backscatter, DEM, incidence-angle and ancillary-mask rasters.
`SAR_download_process.py` is the unified entry point; `execute.sh` calls it.
The existing Sentinel-1 Python commands and flood-warning workflow remain available.

## Unified entry point

Sensor/product and download provider are independent parameters:

| Sensor | Product type | Download source | Processing |
| --- | --- | --- | --- |
| `s1` (default) | `GRD` (default for s1) | `cdse` (default for s1) | CDSE download + SNAP |
| `s1` | `GRD` | `asf` | Submit Sentinel-1 HyP3 RTC jobs |
| `nisar` | `GCOV` (default for nisar) | `asf` (default for nisar) | Download existing HDF5, export backscatter and calculate masks |

Unsupported combinations are rejected before searching/downloading.
**GCOV** means Geocoded Polarimetric Covariance: its diagonal terms contain
terrain-corrected gamma0 backscatter power. **GUNW** means Geocoded Unwrapped
Interferogram: it contains interferometric products from paired acquisitions,
including unwrapped phase. GUNW cannot be passed to the GCOV converter, and is
not supported by this backscatter workflow. The original NISAR search example's
`processingLevel="GUNW"` has therefore been replaced with `"GCOV"`.
See the official [GCOV](https://nisar-docs.asf.alaska.edu/gcov/) and
[GUNW](https://nisar-docs.asf.alaska.edu/gunw/) descriptions.

```bash
# Search GCOV over a rectangle; dates include the entire final UTC day.
python SAR_download_process.py --sensor nisar --product-type GCOV \
  --download-source asf --bbox -81.5 25.0 -80.5 26.5 \
  --start-date 2026-07-01 --end-date 2026-07-07 \
  --search-only --output-json nisar_scenes.json

# Download and fully process exact products (space- or comma-separated names).
python SAR_download_process.py --sensor nisar --names "$NISAR_SCENE" \
  --output-dir ./data --desert-mask-vrt /path/to/desert.vrt

# Process an already downloaded H5 without accessing the ASF catalogue.
python SAR_download_process.py --sensor nisar --h5 /path/to/NISAR_PRODUCT.h5 \
  --output-dir ./data --desert-mask-vrt /path/to/desert.vrt --resolution 30

# Sentinel-1 uses the same entry point and existing processing adapters.
export CDSE_USERNAME="your_username"
export CDSE_PASSWORD="your_password"
python SAR_download_process.py --sensor s1 --product-type GRD \
  --download-source cdse --names "$S1_SCENE" \
  --output-dir ./data --desert-mask-vrt /path/to/desert.vrt
```

The H5 basename must be its original NISAR GCOV product name. Select either
`--names`, `--h5`, or a spatial/time search. Spatial inputs are mutually exclusive:
`--bbox WEST SOUTH EAST NORTH`, `--polygon-wkt "POLYGON(...)"`, or `--extent-file`.
NISAR uses the vector geometry (SHP/GeoJSON etc.) or the raster's WGS84 bounding
box; Sentinel-1 uses the bounding box for catalogue searching, as before.
Search selects intersecting full products; it does not crop output rasters to the AOI.
NISAR accepts UTC ISO timestamps as well as dates; an explicit end timestamp is
exclusive. Sentinel-1 continues to use inclusive `YYYY-MM-DD` dates.
`--names --search-only` writes the supplied names without querying their existence.

Use `--download-only` for NISAR H5 downloads without conversion or masks; this
mode and `--search-only` do not require a desert VRT. Full processing always
generates desert, LIA, snow and cloud masks and requires `--desert-mask-vrt`.
ASF downloads use Earthdata `~/.netrc` (see below); NISAR downloads do not submit
HyP3 jobs or consume HyP3 processing credits.

## NISAR processing and defaults

1. Search ASF GCOV or resolve exact product names, and download one H5 per scene.
2. Read actual diagonal polarizations from `science/LSAR/GCOV/grids/frequencyA`.
   Export `Gamma0_HH.tif`, `Gamma0_HV.tif`, or the other polarizations present.
   `--band LSAR` and `--frequency A` are defaults; SSAR and frequency B may be
   selected when present in the input product. `--polarizations HH HV` selects
   specific layers and fails clearly if any are absent.
3. Align a DEM to the backscatter grid and convert EGM2008 heights to WGS84
   ellipsoid heights. By default, download Copernicus GLO-30 tiles and the
   [PROJ EGM2008 undulation grid](https://cdn.proj.org/) (about 77 MB, cached).
4. Interpolate incidence angle and LOS components in the metadata cube's own
   `(heightAboveEllipsoid, yCoordinates, xCoordinates)` coordinates. Use the
   DEM surface normal and the ground-to-sensor LOS direction to calculate LIA.
5. Generate categorical desert, LIA exclusion, snow and persistent-cloud masks
   on the same reference grid.

The default **`--resolution auto`** retains the original converter's approach:
Rasterio derives an EPSG:4326 grid from the source CRS, dimensions and bounds.
It is **not a fixed 20 m or 30 m**. GCOV source pixel spacing depends on bandwidth;
frequency A commonly has 10 m or 20 m spacing, while frequency B can have 80 m
spacing. Source spacing and the actual output transform are written to
`metadata.json`. Pixel-centre coordinates are correctly converted to outer
pixel bounds before warping.

`--resolution 20` or `--resolution 30` selects nominal metre spacing using
`resolution / 111320` degrees in both axes, consistent with the existing ASF
Sentinel-1 output convention. East-west ground spacing varies with latitude;
this is not a constant-metre projected grid. All NISAR layers use the one
backscatter reference grid. Continuous rasters use bilinear resampling and
categorical masks use nearest-neighbour resampling.

Backscatter remains **linear gamma0 power**. No dB conversion, gamma-to-sigma
factor, or additional division by cosine is applied. Source mask values 0 and
255, nonpositive power and explicit fill values are excluded. H5 backscatter
is processed in blocks; angle interpolation and DEM processing also use blocks.

The interpolated `incidenceAngle` is relative to the ellipsoid normal, not the
terrain normal. LIA is computed from DEM gradients and NISAR LOS X/Y, using
interpolated incidence for inclination. Missing LOS data causes an error;
Sentinel-1's fixed azimuth assumptions are not used. The cube is not extrapolated:
out-of-range coordinates/heights produce nodata, while exact end height levels
are included. One-pixel halos avoid block-edge slope seams.
See the official [metadata](https://nisar-docs.asf.alaska.edu/metadata/) and
[angle/LOS definitions](https://nisar-docs.asf.alaska.edu/glossary/).

The LIA mask is uint8: 0 = keep (LIA <= 50 degrees), 1 = exclude (> 50),
255 = nodata. Set `--lia-threshold` to change the default 50 degrees for NISAR.
Snow/cloud retain the existing 15-day Sentinel-2 SCL algorithm and 0/1/255
encoding; no Sentinel-2 observations yield unknown (255), not a clear-sky mask.
Desert classes are retained from the supplied source, with source nodata (or
-9999 if absent). All three mask routines respect invalid SAR pixels.

To reuse local auxiliary data:

```bash
python SAR_download_process.py --sensor nisar --h5 /path/to/NISAR_PRODUCT.h5 \
  --output-dir ./data --desert-mask-vrt /path/to/desert.vrt \
  --dem /path/to/dem.tif --dem-height-reference egm2008 \
  --geoid /path/to/us_nga_egm08_25.tif --cache-dir /path/to/shared-cache
```

Use `--dem-height-reference ellipsoid` only for a local DEM already expressed
in WGS84 ellipsoidal heights. Otherwise the conversion is `h = H + N`, where
H is EGM2008 orthometric height and N is geoid undulation. The saved `_DEM.tif`
is float32 ellipsoidal height in metres, tagged with its height reference.

```text
data/
  .cache/                           # shared DEM tiles and EGM2008 grid
  <NISAR_PRODUCT>/
    <NISAR_PRODUCT>.h5              # retained; --h5 uses its existing location
    Gamma0_HH.tif                  # actual polarizations, never renamed to VV/VH
    Gamma0_HV.tif
    <NISAR_PRODUCT>_DEM.tif
    <NISAR_PRODUCT>_incidenceAngle.tif
    <NISAR_PRODUCT>_localIncidenceAngle.tif
    <NISAR_PRODUCT>_LIA.tif
    <NISAR_PRODUCT>_desert.tif
    <NISAR_PRODUCT>_ice.tif
    <NISAR_PRODUCT>_cloud.tif
    metadata.json
    processing_status.json
```

`metadata.json` records the selected band/polarizations, radiometry, source and
output spacing, acquisition time, height reference and active output paths.
Changing selected polarizations does not delete older outputs; use the active
paths in this manifest rather than globbing every TIFF in the scene folder.
RAPID's current VV/VH interface and classification parameters still require
separate adaptation/validation before using NISAR HH/HV.

## Retry and checkpoints

The unified runner processes scenes independently, continues after failures,
and exits nonzero if any selected scene fails. HTTP requests have timeouts and
bounded retries; incomplete files retain `.part` names and are restarted on
retry (byte-range continuation is not implemented). NISAR H5 downloads are
checked against available catalogue byte counts and opened/validated as GCOV
before publication. Existing complete H5 files are reused.

For both sensors, `processing_status.json` records stages. A stage is reused
only when its configuration, input file sizes/mtimes, implementation hashes,
and output stamps match and its raster outputs can be read and validated.
`--force` reruns processing stages; it does not redownload a valid H5. Use it
when updating the underlying data referenced by an unchanged desert VRT, or
when refreshing Sentinel-2 observations/remote auxiliary datasets.
Do not run two workers against the same scene directory simultaneously.

Successful unified runs retain source downloads, DEMs and angle products for
inspection/reuse. Per-stage temporary rasters are removed after success;
failed-stage intermediates remain. The flood-warning runner retains its existing
Sentinel-1 cleanup rules; its optional NISAR branch uses the NISAR checkpoints
and retains H5, DEM and angle products.

Individual NISAR tools remain usable: `NISAR_extent_time_download.py` for search
and H5 downloads, `NISAR_GCOV_h5_2_tif.py H5 OUTPUT_DIR` for conversion, and
`NISAR_incidence_angle_interpolate.py H5 REFERENCE_TIF OUTPUT_DIR` for DEM/angles.
`NISAR_specific_name_download_process.py` accepts the unified CLI arguments and
defaults to `--sensor nisar`. Run each with `--help` for its options.

## Validation

Run the offline tests in the configured environment:

```bash
python -m unittest discover -s tests -v
bash -n execute.sh
```

The NISAR tests use actual small HDF5/GeoTIFF fixtures to check pixel-centre
placement, gamma0 preservation, metadata-cube interpolation, EGM2008 conversion,
slope/LOS direction, nodata, download publication, batch failures and checkpoint
reuse. ASF/Earthdata and Sentinel-2 services are mocked. These checks do not
replace a full-scene run with authenticated downloads, real auxiliary data and
the server's GDAL/SNAP environment.

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

The unified `execute.sh` retains downloads and processing rasters for stage
reuse. The legacy flood-warning chain removes intermediate products only after
all four stages succeed. Failed-stage files remain available for inspection.

## Environment

The current Conda environment was exported to `s1pro.yml`:

```bash
conda env create -f s1pro.yml
conda activate s1pro-rtc
```

The workflow also requires ESA SNAP. Make sure the SNAP `bin` directory is
available on `PATH`. The current server-specific SNAP path is set in
`Sentinel_1_specific_name_download_process.py` and may need to be changed on a
different system.

The environment now explicitly includes `h5py`, `pyproj`, `asf-search` and
`hyp3-sdk`. Existing user-specified dependency pins are retained. SNAP is only
needed for the Sentinel-1 CDSE branch. The environment file has not been solved
or installed as a complete environment on this workstation; validate it on the
processing server. To update an existing environment, use
`conda env update -n YOUR_ENV -f s1pro.yml`.

Set `CDSE_USERNAME` and `CDSE_PASSWORD` environment variables for the unified
Sentinel-1 CDSE runner. The legacy direct download command still accepts its
positional credential arguments. NISAR/ASF uses Earthdata credentials instead.

## Running the workflow

Edit `execute.sh` before submission:

- Set `path` to the output root directory.
- Set `desert_mask_vrt` to the external global desert-mask VRT.
- Export `CDSE_USERNAME` and `CDSE_PASSWORD` for Sentinel-1 CDSE.
- Set `SENSOR`, `PRODUCT_TYPE` and `DOWNLOAD_SOURCE`, or use their defaults.
- Add one or more exact SAR product names to `product_names`. The ordinary
  `execute.sh` submission does not perform spatial/time searches; those searches
  remain in `execute_flood_warning.sh` for the daily flood-warning workflow.
- Update the Conda, project, SNAP, Slurm partition, and log paths when running
  on a different system.

Submit the workflow with:

```bash
sbatch execute.sh
# NISAR GCOV with source-derived output spacing:
sbatch --export=ALL,SENSOR=nisar,PRODUCT_TYPE=GCOV,DOWNLOAD_SOURCE=asf execute.sh
# Explicit nominal 30 m spacing:
sbatch --export=ALL,SENSOR=nisar,RESOLUTION=30 execute.sh
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

## Daily flood-warning sensor selection

`execute_flood_warning.sh` defaults to `SENSOR=s1`. Without any new settings,
existing daily submissions still search Sentinel-1 and use CDSE + SNAP;
`DOWNLOAD_SOURCE=asf` still selects Sentinel-1 HyP3 RTC as before.
Set the single new sensor variable to search and process NISAR GCOV instead:

```bash
# Existing daily job, unchanged:
sbatch execute_flood_warning.sh
# Existing Sentinel-1 ASF selection, unchanged:
sbatch --export=ALL,DOWNLOAD_SOURCE=asf execute_flood_warning.sh
# NISAR GCOV: ASF is selected automatically when DOWNLOAD_SOURCE is unset.
sbatch --export=ALL,SENSOR=nisar execute_flood_warning.sh
# Explicit source (also overrides any inherited DOWNLOAD_SOURCE=cdse):
sbatch --export=ALL,SENSOR=nisar,DOWNLOAD_SOURCE=asf execute_flood_warning.sh
```

Alternatively, change `sensor="${SENSOR:-s1}"` in the shell script to use your
preferred default. The existing server paths, Conda activation, logging,
Slurm settings, warning download and region processing are unchanged. Update
that activated environment with the new NISAR dependencies in `s1pro.yml`
before choosing NISAR. NISAR requires Earthdata `.netrc`, not CDSE credentials;
`SENSOR=nisar` combined with `DOWNLOAD_SOURCE=cdse` is rejected.

The Python entry point also accepts an optional sensor:

```bash
python Sentinel_1_ESA_search_download_process_chain_v4.py \
  --sensor nisar --desert-mask-vrt /path/to/desert.vrt
```

Both sensors use the original daily date, GloFAS warning regions, expansion and
area settings, and the >=50 km2 overlap test against original warning pixels.
The only catalogue switch is from Sentinel-1 IW_GRDH_1S/CDSE to NISAR GCOV/ASF.
Scenes are deduplicated before processing. No look-back window or retry on later
catalogue publication has been added; if that day's GCOV products are not yet
available, the job reports zero images.

Sentinel-1 outputs and its `Processed_Sentinel_1_data_path_YYYYMMDD.txt` file
keep their existing locations and behavior. The NISAR branch writes under
`data/YYYY-MM-DD/NISAR/processed_images/` and writes its directory to
`Processed_NISAR_data_path_YYYYMMDD.txt` when at least one scene succeeds.
It invokes the full GCOV pipeline, including desert, LIA, snow and cloud masks,
with the existing NISAR defaults (frequency A, LSAR, resolution auto, LIA 50
degrees). Scene failures retain intermediates, allow subsequent scenes to run,
and give the NISAR daily job a nonzero exit status. Existing downstream readers
of the Sentinel-1 path file are unaffected; RAPID HH/HV adaptation remains
separate from this download/preprocessing change.

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

All previous commands and positional Python arguments remain valid. Add
`--sensor nisar` to the date-range CLI, or `sensor="nisar"` to the Python
function, to search GCOV using the same warnings and overlap rule. The return
schema remains `{"names": [...], "count": N}`; NISAR names have no `.h5` suffix.
For example:

```bash
python Sentinel_1_ESA_search_download_process_chain_v4.py --sensor nisar \
  --search-only --start-date 2026-07-01 --end-date 2026-07-07 \
  --work-directory ./data --output-json nisar_flood_images.json
```

The existing `search_sentinel_with_shape_extent_and_data` function is unchanged.
`run_new_processing_chain(name, output_dir, desert_vrt)` and its fourth
positional `download_source` argument retain the old Sentinel-1 behavior;
passing `sensor="nisar"` selects the full NISAR chain with ASF as its default.

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
another geoid correction. The unified runner also retains SNAP intermediates
for checkpoints; the legacy flood-warning runner removes them after success.

## Download channel

Existing Sentinel-1 commands continue to use Copernicus Data Space + SNAP by default.
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
still run. The legacy flood-warning cleanup removes the intermediate local-angle
raster, while retaining the DEM; the unified runner retains both. The standalone ASF adapter produces the SAR,
DEM, angle and LIA products; the full chain adds the other masks.
