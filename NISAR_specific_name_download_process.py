"""Run a complete NISAR GCOV scene through backscatter, DEM, LIA and masks."""
from pathlib import Path

import rasterio

from NISAR_GCOV_h5_2_tif import convert_gcov, parse_resolution
from NISAR_incidence_angle_interpolate import interpolate_angles
from nisar_common import file_stamp, inspect_h5, local_path, scene_name, write_json
from sar_dem import prepare_dem
from sar_stages import SceneStages, implementation_stamp


def process_scene(h5_file, output_root, desert_mask_vrt, frequency="A", band="LSAR",
                  resolution="auto", polarizations=None, dem=None,
                  dem_height_reference="egm2008", geoid=None, cache_dir=None,
                  lia_threshold=50.0, force=False):
    from Desert_mask import generate_desert_mask
    from Snow_detect import generate_s2_masks

    h5_file = local_path(h5_file)
    name = scene_name(h5_file.name)
    metadata = inspect_h5(h5_file, band, frequency)
    if "zeroDopplerStartTime" not in metadata:
        raise ValueError("NISAR acquisition time is required for snow/cloud masks")
    directory = local_path(output_root) / name
    cache_dir = Path(cache_dir or Path(output_root) / ".cache").resolve()
    desert_mask_vrt = Path(desert_mask_vrt).resolve()
    if not desert_mask_vrt.is_file():
        raise FileNotFoundError(f"Desert mask VRT not found: {desert_mask_vrt}")
    spacing = parse_resolution(resolution)
    code = implementation_stamp("nisar_common.py", "sar_stages.py")
    stages = SceneStages(directory, force)
    config = dict(source=file_stamp(h5_file), band=band, frequency=frequency,
                  resolution=spacing, polarizations=polarizations, common_code=code)
    gamma = stages.run("backscatter", dict(config, code=implementation_stamp("NISAR_GCOV_h5_2_tif.py")),
                       lambda: convert_gcov(h5_file, directory, frequency, band, resolution, polarizations).values())
    reference = gamma[0]
    dem_output = directory / f"{name}_DEM.tif"
    dem_inputs = dict(reference=file_stamp(reference), source=file_stamp(dem) if dem else "cop-dem-glo-30",
                      height_reference=dem_height_reference,
                      geoid=file_stamp(geoid) if geoid else "PROJ_us_nga_egm08_25",
                      code=implementation_stamp("sar_dem.py"), common_code=code)
    stages.run("dem", dem_inputs,
               lambda: [prepare_dem(reference, dem_output, cache_dir, dem, dem_height_reference, geoid)], reference)
    angle_inputs = dict(h5=file_stamp(h5_file), dem=file_stamp(dem_output), reference=file_stamp(reference),
                        band=band, threshold=lia_threshold,
                        code=implementation_stamp("NISAR_incidence_angle_interpolate.py"), common_code=code)
    angles = stages.run("angles_and_lia", angle_inputs,
                        lambda: interpolate_angles(h5_file, reference, dem_output, directory, name,
                                                   band, lia_threshold).values(), reference)
    desert = stages.run("desert", dict(reference=file_stamp(reference), source=file_stamp(desert_mask_vrt),
                                        code=implementation_stamp("Desert_mask.py"), common_code=code),
                         lambda: [generate_desert_mask(name, reference, directory, desert_mask_vrt, strict=True)],
                         reference, require_valid=False)
    snow = stages.run("snow_cloud", dict(reference=file_stamp(reference), time=metadata["zeroDopplerStartTime"],
                                          code=implementation_stamp("Snow_detect.py"), common_code=code),
                       lambda: generate_s2_masks(name, reference, directory, metadata["zeroDopplerStartTime"], strict=True),
                       reference, require_valid=False)
    with rasterio.open(reference) as ref:
        metadata["output_grid"] = dict(crs=ref.crs.to_string(), transform=list(ref.transform)[:6],
                                       width=ref.width, height=ref.height,
                                       resolution_degrees=list(ref.res))
    metadata.update(scene=name, source_h5=str(h5_file), resolution_mode="auto" if spacing is None else spacing,
                    dem_height_reference="WGS84_ellipsoid", lia_threshold_degrees=lia_threshold,
                    polarizations=[Path(p).stem.removeprefix("Gamma0_") for p in gamma],
                    outputs=dict(backscatter=gamma, dem=str(dem_output), angles=angles, desert=desert, snow_cloud=snow))
    write_json(directory / "metadata.json", metadata)
    stages.finish()
    return metadata


def main(argv=None):
    # Keep one parser and one set of validation rules for all entry points.
    from SAR_download_process import main as run
    import sys
    return run(["--sensor", "nisar", *(sys.argv[1:] if argv is None else argv)])


if __name__ == "__main__":
    raise SystemExit(main())
