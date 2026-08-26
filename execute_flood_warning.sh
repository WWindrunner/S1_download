#!/bin/bash

#SBATCH --partition=HydroIntel
#SBATCH --mem=30G
#SBATCH --output=/dev/null

set -o pipefail

: "${CDSE_USERNAME:?Export CDSE_USERNAME before submitting this job}"
: "${CDSE_PASSWORD:?Export CDSE_PASSWORD before submitting this job}"

project_dir="${RAPID_PROJECT_DIR:-/shared/stormcenter/zby3135/RAPID}"
python_bin="${S1PRO_PYTHON:-/shared/stormcenter/zby3135/Software/conda/envs/s1pro/bin/python}"
data_dir="$project_dir/data"

# Set this to the actual global desert-mask VRT on the server.
desert_mask_vrt="${DESERT_MASK_VRT:-/path/to/global_desert_mask.vrt}"

export RAPID_PROJECT_DIR="$project_dir"

mkdir -p "$data_dir"
log_file="$data_dir/Daily_trigger_output_$(date +%Y-%m-%d).txt"
exec >> "$log_file" 2>&1

if [ ! -x "$python_bin" ]; then
    echo "Python executable not found or not executable: $python_bin"
    exit 1
fi

if [ ! -f "$desert_mask_vrt" ]; then
    echo "Desert mask VRT not found: $desert_mask_vrt"
    echo "Update desert_mask_vrt in $0 before scheduling this job."
    exit 1
fi

cd "$project_dir" || exit 1
"$python_bin" Sentinel_1_ESA_search_download_process_chain_v4.py \
    --desert-mask-vrt "$desert_mask_vrt"
status=$?
exit "$status"
