#!/bin/bash
# Build results/YYYY.json from the cache. Runs after the search array; also safe
# to resubmit alone with other thresholds (see S1_archive_search.py merge --help).

#SBATCH --output=/tank/data/SFS/xinyis/shared/src/RAPID_python/RAPID_2.0_Python/S1_download/logs/archive-merge-%j.out
#SBATCH --partition=HydroIntel
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=S1_archive_merge
#SBATCH --time=01:00:00

project_directory="${PROJECT_DIRECTORY:-${SLURM_SUBMIT_DIR:-$(dirname "${BASH_SOURCE[0]}")}}"
cd "$project_directory"

source /tank/data/SFS/xinyis/shared/maopuxu/miniconda3/etc/profile.d/conda.sh
conda activate s1pro
export PYTHONDONTWRITEBYTECODE=1

root="${ARCHIVE_ROOT:-$project_directory/output_search}"
python S1_archive_search.py merge --root "$root" \
    --years "${FIRST_YEAR:-2016}" "${LAST_YEAR:-2025}" "$@"
