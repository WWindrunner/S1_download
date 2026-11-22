"""Unified Sentinel-1 / NISAR search, download and ancillary-mask workflow."""
import argparse
import os
from pathlib import Path
import signal
import subprocess
import sys

from nisar_common import file_stamp, inspect_h5, local_path, scene_name, write_json
from sar_search import date_interval, search_polygon
from sar_stages import SceneStages, implementation_stamp

ROOT = Path(__file__).resolve().parent
DEFAULT_PARALLEL_SNAP_HEAP = "24G"


def process_sentinel(name, args):
    directory = local_path(args.output_dir) / name
    stages = SceneStages(directory, args.force)
    source = args.download_source

    def command(script, *values, timeout=None):
        # Own process group, so a timeout also stops grandchildren (SNAP gpt / Java).
        cmd = [sys.executable, str(ROOT / script), *map(str, values)]
        proc = subprocess.Popen(cmd, start_new_session=True)
        try:
            returncode = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            for sig, grace in ((signal.SIGTERM, 60), (signal.SIGKILL, None)):
                try:
                    os.killpg(proc.pid, sig)
                except ProcessLookupError:
                    break
                try:
                    proc.wait(timeout=grace)
                    break
                except subprocess.TimeoutExpired:
                    continue
            raise TimeoutError(f"{script} exceeded {timeout} s and was stopped") from None
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd)

    def backscatter():
        # CDSE credentials reach the child through the environment, never the command
        # line, which error messages and logs would print.
        command("Sentinel_1_specific_name_download_process.py", name, args.output_dir,
                "--download-source", source, timeout=args.scene_timeout or None)
        paths = [directory / "Gamma0_VV.tif", directory / "Gamma0_VH.tif"]
        if source == "asf":
            paths += [directory / f"{name}_DEM.tif", directory / f"{name}_localIncidenceAngle.tif",
                      directory / "incidenceAngleFromEllipsoid.tif", directory / f"{name}_ls_map.tif"]
        else:
            angles = list(directory.glob("*incidenceAngleFromEllipsoid.tif"))
            if len(angles) != 1:
                raise ValueError("Expected exactly one Sentinel-1 ellipsoid incidence raster")
            paths += angles + list(directory.glob("*_DEM.tif"))
        return paths

    stages.run("backscatter", dict(sensor="s1", scene=name, source=source,
                                   code=implementation_stamp("Sentinel_1_specific_name_download_process.py",
                                                             "Sentinel_1_specific_name_ASF.py")), backscatter)
    reference = directory / "Gamma0_VV.tif"
    # Both sources have the ellipsoid incidence angle and a DEM on the Gamma0 grid.
    angle = next(directory.glob("*incidenceAngleFromEllipsoid.tif"))

    def lia():
        command("cal_LIA.py", name, directory, "--incidence-angle", angle, "--metadata-dir", directory)
        return [directory / f"{name}_LIA.tif"]

    stages.run("lia", dict(angle=file_stamp(angle), dem=[file_stamp(p) for p in directory.glob("*_DEM.tif")],
                           code=implementation_stamp("cal_LIA.py")), lia, reference)
    from Desert_mask import generate_desert_mask
    from Snow_detect import generate_s2_masks

    stages.run("desert", dict(reference=file_stamp(reference), source=file_stamp(args.desert_mask_vrt),
                              code=implementation_stamp("Desert_mask.py", "nisar_common.py")),
               lambda: [generate_desert_mask(name, reference, directory, args.desert_mask_vrt, strict=True)],
               reference, require_valid=False)
    stages.run("snow_cloud", dict(reference=file_stamp(reference), code=implementation_stamp("Snow_detect.py")),
               lambda: generate_s2_masks(name, reference, directory, strict=True), reference, require_valid=False)
    stages.finish()


def _process_sentinel_worker(name, args):
    """Run one scene in a worker; return errors as plain text, because some exceptions
    (e.g. HTTPError holding an open response) cannot be sent back between processes."""
    try:
        process_sentinel(name, args)
    except Exception as exc:
        raise RuntimeError(f"{type(exc).__name__}: {exc}") from None


def process_sentinel_parallel(names, args, asf_failures):
    """Process Sentinel-1 scenes in args.workers processes; one failure affects only its scene."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # SNAP limits reach the gpt child processes through the environment.
    cpus = len(os.sched_getaffinity(0))
    os.environ["SNAP_THREADS"] = str(args.snap_threads or max(1, cpus // args.workers))
    os.environ["SNAP_HEAP"] = args.snap_heap or DEFAULT_PARALLEL_SNAP_HEAP
    print(f"Parallel: {args.workers} scenes at a time, SNAP {os.environ['SNAP_THREADS']} threads / "
          f"{os.environ['SNAP_HEAP']} heap each ({cpus} CPUs available)", flush=True)
    failures = []
    for name in names:
        if name in asf_failures:
            failures.append(name)
            print(f"Failed {name}: {asf_failures[name]}; intermediates retained.", file=sys.stderr)
    todo = [name for name in names if name not in asf_failures]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(_process_sentinel_worker, name, args): name for name in todo}
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
                print(f"Completed: {name}", flush=True)
            except Exception as exc:
                failures.append(name)
                print(f"Failed {name}: {exc}; intermediates retained.", file=sys.stderr, flush=True)
    print(f"Finished: {len(names) - len(failures)} succeeded, {len(failures)} failed")
    return 1 if failures else 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sensor", choices=("s1", "nisar"), default="s1")
    parser.add_argument("--product-type", type=str.upper, help="Default: GRD for s1, GCOV for nisar; GUNW is unsupported")
    parser.add_argument("--download-source", choices=("cdse", "asf"),
                        help="s1: cdse (CDSE + SNAP, Copernicus 30 m DEM; default) or asf (HyP3 RTC); "
                             "nisar: asf (default)")
    parser.add_argument("--output-dir", default="data")
    inputs = parser.add_mutually_exclusive_group()
    inputs.add_argument("--names", nargs="+", help="Exact product names (space- or comma-separated)")
    inputs.add_argument("--h5", nargs="+", help="Local NISAR H5 paths; bypass download")
    area = parser.add_mutually_exclusive_group()
    area.add_argument("--extent-file")
    area.add_argument("--bbox", nargs=4, type=float, metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    area.add_argument("--polygon-wkt", help="EPSG:4326 polygon; s1 searches its bounding box")
    parser.add_argument("--start-date")
    parser.add_argument("--end-date")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--search-only", action="store_true")
    mode.add_argument("--download-only", action="store_true", help="NISAR: download H5 without processing")
    parser.add_argument("--output-json", help="Write selected product names and count")
    parser.add_argument("--desert-mask-vrt")
    parser.add_argument("--resolution", default="auto", help="NISAR: auto (default), or nominal metres (20, 30, ...)")
    parser.add_argument("--frequency", choices=("A", "B"), default="A")
    parser.add_argument("--band", choices=("LSAR", "SSAR"), default="LSAR")
    parser.add_argument("--polarizations", nargs="+", choices=("HH", "HV", "VV", "VH", "RH", "RV"))
    parser.add_argument("--dem", help="NISAR: optional local DEM")
    parser.add_argument("--dem-height-reference", choices=("egm2008", "ellipsoid"), default="egm2008")
    parser.add_argument("--geoid", help="NISAR: local EGM2008 geoid undulation raster")
    parser.add_argument("--cache-dir", help="NISAR: shared DEM/geoid cache (default OUTPUT_DIR/.cache)")
    parser.add_argument("--lia-threshold", type=float, default=50.0)
    parser.add_argument("--force", action="store_true", help="Reprocess stages, retaining complete downloads")
    parser.add_argument("--workers", type=int, default=1,
                        help="s1: scenes processed in parallel (separate processes). CDSE downloads stay "
                             "limited to CDSE_MAX_DOWNLOADS (default 4) across all processes of the user.")
    parser.add_argument("--snap-threads", type=int,
                        help="s1 CDSE: SNAP gpt threads per scene. Default with --workers > 1: "
                             "available CPUs / workers")
    parser.add_argument("--scene-timeout", type=int, default=4 * 3600,
                        help="s1: seconds allowed for one scene's download + backscatter processing "
                             "before it is stopped and reported as failed (default 14400; 0 disables)")
    parser.add_argument("--snap-heap",
                        help="s1 CDSE: SNAP Java heap per scene, e.g. 24G. Default with --workers > 1: "
                             f"{DEFAULT_PARALLEL_SNAP_HEAP}; otherwise the SNAP install setting")
    return parser


def validate_args(parser, args):
    args.product_type = args.product_type or ("GCOV" if args.sensor == "nisar" else "GRD")
    args.download_source = args.download_source or ("asf" if args.sensor == "nisar" else "cdse")
    supported = {("s1", "GRD", "cdse"), ("s1", "GRD", "asf"), ("nisar", "GCOV", "asf")}
    if (args.sensor, args.product_type, args.download_source) not in supported:
        parser.error("Supported combinations: s1/GRD/cdse, s1/GRD/asf, nisar/GCOV/asf. GUNW needs a separate interferometric workflow.")
    if args.sensor != "nisar" and (args.h5 or args.download_only or args.resolution != "auto" or args.dem
                                   or args.geoid or args.polarizations or args.band != "LSAR" or args.frequency != "A"
                                   or args.dem_height_reference != "egm2008" or args.lia_threshold != 50):
        parser.error("H5, download-only, resolution, DEM, frequency, band, polarization and LIA options are NISAR-only")
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    if args.scene_timeout < 0:
        parser.error("--scene-timeout must be 0 (disabled) or a positive number of seconds")
    if args.sensor != "s1" and (args.workers != 1 or args.snap_threads or args.snap_heap):
        parser.error("--workers, --snap-threads and --snap-heap are Sentinel-1 options")
    if args.h5 and (args.search_only or args.download_only):
        parser.error("Local --h5 is for processing; it cannot be combined with search/download-only")
    spatial = args.extent_file or args.bbox or args.polygon_wkt
    if (args.names or args.h5) and (spatial or args.start_date or args.end_date):
        parser.error("Choose exact names/local H5 OR a spatial/time search")
    if not (args.names or args.h5) and not (spatial and args.start_date and args.end_date):
        parser.error("Supply --names, --h5, or a spatial extent with --start-date and --end-date")
    if args.sensor == "nisar":
        from NISAR_GCOV_h5_2_tif import parse_resolution
        try:
            parse_resolution(args.resolution)
        except ValueError as exc:
            parser.error(str(exc))
        if not 0 <= args.lia_threshold <= 180:
            parser.error("LIA threshold must be between 0 and 180")
        if args.dem_height_reference == "ellipsoid" and not args.dem:
            parser.error("--dem-height-reference ellipsoid requires a local --dem")
    if not args.search_only and not args.download_only:
        if not args.desert_mask_vrt or not Path(args.desert_mask_vrt).is_file():
            parser.error("Full processing requires an existing --desert-mask-vrt")
        if args.sensor == "s1" and args.download_source == "cdse":
            if not os.environ.get("CDSE_USERNAME") or not os.environ.get("CDSE_PASSWORD"):
                parser.error("Set CDSE_USERNAME and CDSE_PASSWORD for Sentinel-1 CDSE downloads")
    args.output_dir = str(Path(args.output_dir).expanduser().resolve())


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    products = {}
    local = {}
    try:
        if args.h5:
            local = {scene_name(Path(p).name): Path(p).expanduser().resolve() for p in args.h5}
            for path in local.values():
                inspect_h5(path, args.band, args.frequency)
            names = list(local)
        elif args.names:
            names = list(dict.fromkeys(n.strip() for value in args.names for n in value.split(",") if n.strip()))
            if args.sensor == "nisar":
                names = list(dict.fromkeys(scene_name(n) for n in names))
            else:
                names = list(dict.fromkeys(n.removesuffix(".SAFE") for n in names))
                if any(Path(n).name != n or any(c in n for c in '/\\:*?') or n in (".", "..") for n in names):
                    raise ValueError("Expected product names, not paths")
            if not names:
                raise ValueError("At least one product name is required")
        else:
            polygon = search_polygon(args.extent_file, args.bbox, args.polygon_wkt)
            date_interval(args.start_date, args.end_date)
            if args.sensor == "nisar":
                from NISAR_extent_time_download import search_nisar, product_name
                results = search_nisar(args.start_date, args.end_date, polygon.wkt)
                products = {product_name(p): p for p in results}
                names = list(products)
            else:
                from Sentinel_1_search_by_extent import search_sentinel_1_product_names
                west, south, east, north = polygon.bounds
                names = search_sentinel_1_product_names("none", args.start_date, args.end_date, west, east, south, north)
        print(f"Selected {len(names)} {args.sensor}/{args.product_type} products")
        for name in names:
            print(name)
        if args.output_json:
            write_json(args.output_json, dict(sensor=args.sensor, product_type=args.product_type,
                                              names=names, count=len(names)))
        if args.search_only or not names:
            return 0
    except Exception as exc:
        print(f"Product selection failed: {exc}", file=sys.stderr)
        return 1

    failures = []
    # Run every ASF RTC job concurrently up front; each scene's backscatter stage then
    # finds its outputs. A failed job only excludes its own scene.
    asf_failures = {}
    if args.sensor == "s1" and args.download_source == "asf":
        try:
            from Sentinel_1_specific_name_ASF import download_and_process
            asf_failures = download_and_process(names, args.output_dir)
        except Exception as exc:
            asf_failures = {name: f"ASF batch failed: {exc}" for name in names}
    if args.sensor == "s1" and args.workers > 1:
        return process_sentinel_parallel(names, args, asf_failures)
    if args.snap_threads:
        os.environ["SNAP_THREADS"] = str(args.snap_threads)
    if args.snap_heap:
        os.environ["SNAP_HEAP"] = args.snap_heap
    for name in names:
        try:
            if name in asf_failures:
                raise RuntimeError(asf_failures[name])
            if args.sensor == "s1":
                process_sentinel(name, args)
            else:
                h5 = local.get(name)
                if h5 is None:
                    from NISAR_extent_time_download import download_product, find_product
                    product = products.get(name) or find_product(name)
                    h5 = download_product(product, args.output_dir, args.band, args.frequency)
                if not args.download_only:
                    from NISAR_specific_name_download_process import process_scene
                    process_scene(h5, args.output_dir, args.desert_mask_vrt,
                                  frequency=args.frequency, band=args.band, resolution=args.resolution,
                                  polarizations=args.polarizations, dem=args.dem,
                                  dem_height_reference=args.dem_height_reference, geoid=args.geoid,
                                  cache_dir=args.cache_dir, lia_threshold=args.lia_threshold, force=args.force)
            print(f"Completed: {name}")
        except Exception as exc:
            failures.append(name)
            print(f"Failed {name}: {exc}; intermediates retained.", file=sys.stderr)
    print(f"Finished: {len(names) - len(failures)} succeeded, {len(failures)} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
