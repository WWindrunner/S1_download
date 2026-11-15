#!/bin/bash

#SBATCH --output=/tank/data/SFS/xinyis/shared/src/RAPID_python/RAPID_2.0_Python/S1_download/logs/log-%j.out
#SBATCH --partition=HydroIntel
#SBATCH --mem=50G
#SBATCH --job-name=S1_search
#SBATCH --time=10:00:00

set -e
# Slurm runs a spool copy of this script; use the submission directory there.
project_directory="${PROJECT_DIRECTORY:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}"
cd "$project_directory"

source /tank/data/SFS/xinyis/shared/maopuxu/miniconda3/etc/profile.d/conda.sh
conda activate s1pro

start_date="${START_DATE:-2019-01-01}"
end_date="${END_DATE:-2020-01-01}"
work_directory="${WORK_DIRECTORY:-/tank/data/SFS/xinyis/shared/src/RAPID_python/RAPID_2.0_Python/S1_download/output_search}"
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
