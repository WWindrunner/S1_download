#!/bin/bash

#SBATCH --partition=HydroIntel
#SBATCH --mem=30G
#SBATCH --output=/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src/download_and_process_%j.out


source /tank/data/SFS/xinyis/zhao89/software/conda/bin/activate
#source /home/uwm/maopuxu/miniconda3/bin/activate
conda activate s1pro

cd /tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src

#s1name="S1A_IW_GRDH_1SDV_20240301T020957_20240301T021022_052782_066322_5F93"
#path="/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/past_events/20260616_mask"

#python Sentinel_1_specific_name_download_process.py $s1name $path
#python cal_LIA.py $s1name $path

path="/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/past_events/20260616_mask"
desert_mask_vrt="/path/to/global_desert_mask.vrt"
extent_file="/path/to/search_extent.tif"
start_date="2025-11-15"
end_date="2025-11-20"

# Manual product-name input (kept as a reference):
# s1names=(
# "S1A_IW_GRDH_1SDV_20240301T020957_20240301T021022_052782_066322_5F93"
# )

search_output=$(python Sentinel_1_search_by_extent.py \
    "$extent_file" \
    "$start_date" \
    "$end_date" \
    --names-only)
search_status=$?

if [ "$search_status" -ne 0 ]; then
    echo "Sentinel-1 image search failed."
    exit "$search_status"
fi

if [ -z "$search_output" ]; then
    echo "No Sentinel-1 images found for the specified extent and date range."
    exit 0
fi

mapfile -t s1names <<< "$search_output"
echo "Found ${#s1names[@]} Sentinel-1 images."

for s1name in "${s1names[@]}"; do
    echo "Processing $s1name"

    python Sentinel_1_specific_name_download_process.py "$s1name" "$path"
    python Desert_mask.py "$s1name" "$path" "$desert_mask_vrt"
    python cal_LIA.py "$s1name" "$path"
    python Snow_detect.py "$s1name" "$path"
done
