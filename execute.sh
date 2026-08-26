#!/bin/bash

#SBATCH --partition=HydroIntel
#SBATCH --mem=30G
#SBATCH --output=/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src/download_and_process_%j.out

shopt -s nullglob

source /tank/data/SFS/xinyis/zhao89/software/conda/bin/activate
conda activate s1pro

cd /tank/data/SFS/xinyis/FS650/maopuxu/lab_2/src

path="/tank/data/SFS/xinyis/FS650/maopuxu/lab_2/past_events/20260616_mask"
desert_mask_vrt="/path/to/global_desert_mask.vrt"
extent_file="/path/to/search_extent.tif"
start_date="2025-11-15"
end_date="2025-11-20"
xmin="-95.0"
xmax="-90.0"
ymin="40.0"
ymax="45.0"

# Copernicus Data Space credentials.
username=""
password=""

if [ -z "$username" ] || [ -z "$password" ]; then
    echo "Set the Copernicus Data Space username and password in execute.sh."
    exit 1
fi

# Manual product-name input (kept as a reference):
# s1names=(
# "S1A_IW_GRDH_1SDV_20240301T020957_20240301T021022_052782_066322_5F93"
# )

search_output=$(python Sentinel_1_search_by_extent.py \
    "$extent_file" \
    "$start_date" \
    "$end_date" \
    "$xmin" \
    "$xmax" \
    "$ymin" \
    "$ymax" \
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
printf '  %s\n' "${s1names[@]}"

cleanup_product_intermediates() {
    local product_dir="$path/$1"

    if [ ! -d "$product_dir" ]; then
        echo "Cleanup skipped: product directory not found: $product_dir"
        return
    fi

    rm -rf -- \
        "$product_dir/snow_temp" \
        "$product_dir/dem_tiles"
    rm -f -- \
        "$product_dir"/*.zip \
        "$product_dir"/*incidenceAngleFromEllipsoid.tif \
        "$product_dir"/*localIncidenceAngle.tif \
        "$product_dir"/*_manifest.safe \
        "$product_dir"/*_proc.xml \
        "$product_dir"/*_gamma0-elp.tif \
        "$product_dir"/DEM_merged.tif \
        "$product_dir"/DEM_merged_res.tif
}

for s1name in "${s1names[@]}"; do
    echo "Processing $s1name"
    product_dir="$path/$s1name"

    if ! python Sentinel_1_specific_name_download_process.py \
        "$s1name" \
        "$path" \
        "$username" \
        "$password"; then
        echo "Sentinel-1 processing failed for $s1name; intermediates retained."
        continue
    fi

    if ! python Desert_mask.py "$s1name" "$path" "$desert_mask_vrt"; then
        echo "Desert-mask processing failed for $s1name; intermediates retained."
        continue
    fi

    incidence_angles=("$product_dir"/*incidenceAngleFromEllipsoid.tif)
    if [ "${#incidence_angles[@]}" -ne 1 ]; then
        echo "Expected exactly one incidence-angle raster for $s1name; found ${#incidence_angles[@]}."
        continue
    fi

    if ! python cal_LIA.py \
        "$s1name" \
        "$product_dir" \
        --incidence-angle "${incidence_angles[0]}" \
        --metadata-dir "$product_dir"; then
        echo "LIA processing failed for $s1name; intermediates retained."
        continue
    fi

    if ! python Snow_detect.py \
        "$s1name" \
        "$product_dir/Gamma0_VV.tif" \
        "$product_dir"; then
        echo "Snow/cloud processing failed for $s1name; intermediates retained."
        continue
    fi

    cleanup_product_intermediates "$s1name"
    echo "Completed and cleaned intermediate files for $s1name."
done
