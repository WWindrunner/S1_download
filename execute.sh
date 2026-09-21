#!/bin/bash

#SBATCH --partition=HydroIntel
#SBATCH --mem=30G
#SBATCH --output=/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src/download_and_process_%j.out

source /tank/data/SFS/xinyis/zhao89/software/conda/bin/activate
conda activate "${CONDA_ENV:-s1pro-rtc}" || exit 1
cd /tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src || exit 1

path="/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/past_events/20260616_mask"
desert_mask_vrt="/path/to/global_desert_mask.vrt"
extent_file="/path/to/search_extent.tif"
start_date="2025-11-15"
end_date="2025-11-20"
xmin="-95.0"
xmax="-90.0"
ymin="40.0"
ymax="45.0"

# Sensor/product and download provider are separate settings.
# s1/GRD: cdse (SNAP) or asf (HyP3 RTC); nisar/GCOV: asf (existing H5).
sensor="${SENSOR:-s1}"
if [ "$sensor" = "nisar" ]; then
    product_type="${PRODUCT_TYPE:-GCOV}"
    download_source="${DOWNLOAD_SOURCE:-asf}"
else
    product_type="${PRODUCT_TYPE:-GRD}"
    download_source="${DOWNLOAD_SOURCE:-cdse}"
fi

# Required only for Sentinel-1 CDSE; ASF uses Earthdata ~/.netrc.
export CDSE_USERNAME="${CDSE_USERNAME:-}"
export CDSE_PASSWORD="${CDSE_PASSWORD:-}"

# Set exact product names to bypass spatial/time searching.
product_names=()
# product_names=("S1A_IW_GRDH_1SDV_20240301T020957_20240301T021022_052782_066322_5F93")
# product_names=("NISAR_L2_PR_GCOV_024_156_D_069_2005_QPDH_A_20260707T004257_20260707T004331_P05023_N_F_J_001")

args=(--sensor "$sensor" --product-type "$product_type"
      --download-source "$download_source" --output-dir "$path"
      --desert-mask-vrt "$desert_mask_vrt")
if [ "${#product_names[@]}" -gt 0 ]; then
    args+=(--names "${product_names[@]}")
else
    args+=(--start-date "$start_date" --end-date "$end_date")
    if [ -f "$extent_file" ]; then
        args+=(--extent-file "$extent_file")
    else
        args+=(--bbox "$xmin" "$ymin" "$xmax" "$ymax")
    fi
fi
if [ "$sensor" = "nisar" ]; then
    # auto preserves source-derived WGS84 spacing; 20/30 are nominal metres.
    args+=(--resolution "${RESOLUTION:-auto}" --frequency "${FREQUENCY:-A}")
fi

# Extra CLI switches may be passed to bash/sbatch (e.g. --search-only or --force).
python SAR_download_process.py "${args[@]}" "$@"
exit $?
