#!/bin/bash

#SBATCH --partition=priority
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --job-name=daily_S1_download
#SBATCH --time=10:00:00

source ~/.bashrc
source /shared/stormcenter/Linzq25/E001_MAE_Bathymetry/miniconda3/etc/profile.d/conda.sh

conda activate /shared/stormcenter/Linzq25/E001_MAE_Bathymetry/miniconda3/envs/s1pro-snap12

SNAP_HOME="/gpfs/sharedfs1/manoslab/CREST_app/apps/snap12"
export PATH="$SNAP_HOME/bin:$PATH"

set -o pipefail

: "${CDSE_USERNAME:?Export CDSE_USERNAME before submitting this job}"
: "${CDSE_PASSWORD:?Export CDSE_PASSWORD before submitting this job}"

project_dir="/shared/stormcenter/zby3135/RAPID/S1_download"
data_dir="/shared/stormcenter/zby3135/RAPID/data"
desert_mask_vrt="/shared/stormcenter/Shen/retrieval/RAPID/desert_mask/desert.vrt"

export RAPID_PROJECT_DIR="/shared/stormcenter/zby3135/RAPID"

mkdir -p "$data_dir"
log_file="$data_dir/Daily_trigger_output_$(date +%Y-%m-%d).txt"
exec >> "$log_file" 2>&1

if [ ! -f "$desert_mask_vrt" ]; then
    echo "Desert mask VRT not found: $desert_mask_vrt"
    echo "Update desert_mask_vrt in $0 before scheduling this job."
    exit 1
fi

cd "$project_dir" || exit 1
python Sentinel_1_ESA_search_download_process_chain_v4.py --desert-mask-vrt "$desert_mask_vrt"
status=$?
exit "$status"
