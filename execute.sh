#!/bin/bash

#SBATCH --partition=HydroIntel
#SBATCH --mem=30G
#SBATCH --output=/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src/download_and_process_%j.out

source /tank/data/SFS/xinyis/zhao89/software/conda/bin/activate
conda activate "${CONDA_ENV:-s1pro-rtc}" || exit 1
cd /tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src || exit 1

path="/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/past_events/20260616_mask"
desert_mask_vrt="/path/to/global_desert_mask.vrt"

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

# Set one or more exact SAR product names. Spatial/time searches belong to the
# daily flood-warning workflow in execute_flood_warning.sh.
product_names=()
# product_names=("S1A_IW_GRDH_1SDV_20240301T020957_20240301T021022_052782_066322_5F93")
# product_names=("NISAR_L2_PR_GCOV_024_156_D_069_2005_QPDH_A_20260707T004257_20260707T004331_P05023_N_F_J_001")

if [ "${#product_names[@]}" -eq 0 ]; then
    echo "Set at least one exact SAR product name in product_names before submitting."
    exit 1
fi

args=(--sensor "$sensor" --product-type "$product_type"
      --download-source "$download_source" --output-dir "$path"
      --desert-mask-vrt "$desert_mask_vrt"
      --names "${product_names[@]}")
if [ "$sensor" = "nisar" ]; then
    # auto preserves source-derived WGS84 spacing; 20/30 are nominal metres.
    args+=(--resolution "${RESOLUTION:-auto}" --frequency "${FREQUENCY:-A}")
fi

# Extra processing switches may be passed to bash/sbatch (e.g. --force).
python SAR_download_process.py "${args[@]}" "$@"
exit $?
