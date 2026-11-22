#!/bin/bash
# One Slurm array task searches one calendar month into the shared cache.
# Submit through submit_S1_archive.sh, which sets the array range and the merge job.

#SBATCH --output=/tank/data/SFS/xinyis/shared/src/RAPID_python/RAPID_2.0_Python/S1_download/logs/archive-%A_%a.out
#SBATCH --partition=HydroIntel
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --job-name=S1_archive
#SBATCH --time=12:00:00

set -e
project_directory="${PROJECT_DIRECTORY:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}"
cd "$project_directory"

source /tank/data/SFS/xinyis/shared/maopuxu/miniconda3/etc/profile.d/conda.sh
conda activate s1pro
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1

start_month="${START_MONTH:-2016-01}"
root="${ARCHIVE_ROOT:-$project_directory/output_search}"
first=$(date -u -d "${start_month}-01 +${SLURM_ARRAY_TASK_ID} month" +%F)
last=$(date -u -d "$first +1 month -1 day" +%F)

python S1_archive_search.py search --start "$first" --end "$last" --root "$root" "$@"
