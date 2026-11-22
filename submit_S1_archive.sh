#!/bin/bash
# Submit the month-by-month archive search and a dependent merge job.
#   bash submit_S1_archive.sh                    # 2016-01 .. 2025-12
#   START_MONTH=2020-01 END_MONTH=2020-12 bash submit_S1_archive.sh
# Rerunning resumes: cached days are skipped. PARALLEL limits concurrent months;
# GloFAS and the CDSE catalogue start rejecting requests well below 32 clients.
set -e
cd "$(dirname "${BASH_SOURCE[0]}")"
start_month="${START_MONTH:-2016-01}"
end_month="${END_MONTH:-2025-12}"
parallel="${PARALLEL:-6}"
months=$(( ( $(date -u -d "${end_month}-01" +%Y) - $(date -u -d "${start_month}-01" +%Y) ) * 12
           + 10#$(date -u -d "${end_month}-01" +%m) - 10#$(date -u -d "${start_month}-01" +%m) + 1 ))
mkdir -p logs
array=$(sbatch --parsable --array=0-$((months - 1))%${parallel} \
    --export=ALL,START_MONTH="$start_month" S1_archive_search.sh)
merge=$(sbatch --parsable --dependency=afterany:${array} \
    --export=ALL,FIRST_YEAR="${start_month%-*}",LAST_YEAR="${end_month%-*}" S1_archive_merge.sh)
echo "Search array: ${array} (${months} months, ${parallel} at a time); merge job: ${merge}"
