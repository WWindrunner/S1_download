"""Search and download exact NISAR GCOV products from ASF/Earthdata."""
import argparse
from datetime import datetime
from pathlib import Path
from urllib.parse import unquote, urlparse

from nisar_common import download_http, inspect_h5, local_path, scene_name, write_json
from sar_search import date_interval, search_polygon


def product_name(product):
    properties = product.properties
    value = properties.get("sceneName") or properties.get("fileName")
    return scene_name(value)


def search_nisar(start_date, end_date, polygon, product_type="GCOV"):
    import asf_search as asf

    if product_type != "GCOV":
        raise ValueError("The backscatter workflow supports GCOV; GUNW is an interferometric product")
    start, end = date_interval(start_date, end_date)
    results = asf.search(dataset="NISAR", processingLevel="GCOV", intersectsWith=polygon,
                         start=start.isoformat(), end=end.isoformat())
    # ASF may include acquisitions on the end boundary; use [start, end).
    unique = {}
    for product in results:
        acquired = datetime.fromisoformat(product.properties["startTime"].replace("Z", "+00:00"))
        if start <= acquired < end:
            unique[product_name(product)] = product
    return [unique[name] for name in sorted(unique)]


def find_product(name):
    import asf_search as asf

    name = scene_name(name)
    matches = [p for p in asf.granule_search([name]) if product_name(p) == name]
    if len(matches) != 1:
        raise ValueError(f"Expected one exact ASF GCOV match for {name}, found {len(matches)}")
    return matches[0]


def download_product(product, output_root, band="LSAR", frequency="A"):
    import asf_search as asf

    name = product_name(product)
    folder = local_path(output_root) / name
    target = folder / (name + ".h5")
    properties = product.properties
    urls = [properties.get("url"), *(properties.get("additionalUrls") or [])]
    urls = list(dict.fromkeys(url for url in urls if isinstance(url, str)
                             and Path(unquote(urlparse(url).path)).name == target.name))
    if len(urls) != 1:
        raise ValueError(f"Expected one primary HDF5 download URL for {name}, found {len(urls)}")
    sizes = properties.get("bytes")
    expected_size = sizes.get(target.name, {}).get("bytes") if isinstance(sizes, dict) else None
    if target.exists():
        try:
            if expected_size is not None and target.stat().st_size != int(expected_size):
                raise ValueError("Cached H5 size mismatch")
            inspect_h5(target, band, frequency)
            print(f"Reusing downloaded H5: {target}")
            return target
        except (OSError, ValueError, KeyError):
            print(f"Replacing incomplete/invalid H5: {target}")
    # ASFSession handles Earthdata redirect authentication using ~/.netrc.
    candidate = target.with_name(target.name + ".download")
    with asf.ASFSession() as session:
        download_http(urls[0], candidate, session=session, expected_size=expected_size)
    inspect_h5(candidate, band, frequency)
    candidate.replace(target)
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="data")
    parser.add_argument("--start-date")
    parser.add_argument("--end-date")
    parser.add_argument("--names", nargs="+", help="Exact names instead of a spatial search")
    area = parser.add_mutually_exclusive_group()
    area.add_argument("--polygon-wkt")
    area.add_argument("--extent-file")
    area.add_argument("--bbox", nargs=4, type=float, metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    parser.add_argument("--product-type", default="GCOV", choices=("GCOV",))
    parser.add_argument("--band", choices=("LSAR", "SSAR"), default="LSAR")
    parser.add_argument("--frequency", choices=("A", "B"), default="A")
    parser.add_argument("--search-only", action="store_true")
    parser.add_argument("--output-json")
    args = parser.parse_args(argv)
    if args.names:
        results = [find_product(name) for name in dict.fromkeys(args.names)]
    else:
        if not args.start_date or not args.end_date:
            parser.error("Spatial search requires --start-date and --end-date")
        polygon = search_polygon(args.extent_file, args.bbox, args.polygon_wkt)
        results = search_nisar(args.start_date, args.end_date, polygon.wkt, args.product_type)
    names = [product_name(p) for p in results]
    print(f"Found {len(names)} GCOV products")
    for name in names:
        print(name)
    if args.output_json:
        write_json(args.output_json, dict(names=names, count=len(names)))
    if not args.search_only:
        for product in results:
            download_product(product, args.output_dir, args.band, args.frequency)


if __name__ == "__main__":
    main()
