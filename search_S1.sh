#!/bin/bash

#SBATCH --output=/shared/stormcenter/Shen/retrieval/RAPID/S1_download/log-%j.out
#SBATCH --partition=priority
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --job-name=S1_search
#SBATCH --time=10:00:00

set -e
# Slurm runs a spool copy of this script; use the submission directory there.
project_directory="${PROJECT_DIRECTORY:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}"
cd "$project_directory"

source /shared/stormcenter/Linzq25/E001_MAE_Bathymetry/miniconda3/etc/profile.d/conda.sh
conda activate /shared/stormcenter/Linzq25/E001_MAE_Bathymetry/miniconda3/envs/s1pro-snap12

start_date="${START_DATE:-2025-06-01}"
end_date="${END_DATE:-2025-07-01}"
work_directory="${WORK_DIRECTORY:-/shared/stormcenter/Shen/retrieval/RAPID/S1_download/temp}"
output_json="${OUTPUT_JSON:-${work_directory}/result_${start_date//-/}_${end_date//-/}.json}"

mkdir -p "$work_directory" "$(dirname "$output_json")"

python Sentinel_1_ESA_search_download_process_chain_v4.py \
    --search-only \
    --clean-search \
    --names-only \
    --start-date "$start_date" \
    --end-date "$end_date" \
    --work-directory "$work_directory" \
    --output-json "$output_json" \
    "$@"
