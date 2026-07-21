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
   - Downloads Copernicus DEM tiles from Microsoft Planetary Computer.
   - Resamples the DEM to the Sentinel-1 grid and calculates local incidence
     angle (LIA).
4. `Snow_detect.py`
   - Searches for Sentinel-2 L2A products from the 15 days before the
     Sentinel-1 acquisition.
   - Uses the Sentinel-2 Scene Classification Layer (SCL) to calculate snow
     occurrence on the Sentinel-1 grid.
   - Removes downloaded ZIP files and intermediate incidence-angle rasters
     after successful completion.

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

Copernicus Data Space credentials are currently configured in
`Sentinel_1_specific_name_download_process.py`.

## Running the workflow

Edit `execute.sh` before submission:

- Set `path` to the output root directory.
- Set `desert_mask_vrt` to the external global desert-mask VRT.
- Add one or more Sentinel-1 product names to the `s1names` array.
- Update the Conda, project, SNAP, Slurm partition, and log paths when running
  on a different system.

Submit the workflow with:

```bash
sbatch execute.sh
```

The four stages can also be run manually:

```bash
python Sentinel_1_specific_name_download_process.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER>
python Desert_mask.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER> <DESERT_MASK_VRT>
python cal_LIA.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER>
python Snow_detect.py <S1_PRODUCT_NAME> <OUTPUT_FOLDER>
```

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
```

Additional SNAP output files may also be present in the product directory.
