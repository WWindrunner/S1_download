import argparse
import os
from datetime import date, timedelta

import geopandas as gpd
import rasterio
import requests
from rasterio.warp import transform_bounds


def read_search_extent(input_path):
    """Return the bounding box of a shapefile or raster in EPSG:4326."""
    input_path = os.path.abspath(os.path.expanduser(input_path))
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input extent file not found: {input_path}")

    extension = os.path.splitext(input_path)[1].lower()

    if extension == ".shp":
        gdf = gpd.read_file(input_path)
        if gdf.empty:
            raise ValueError(f"Input shapefile contains no features: {input_path}")
        if gdf.crs is None:
            raise ValueError(f"Input shapefile has no CRS: {input_path}")

        minx, miny, maxx, maxy = gdf.to_crs(epsg=4326).total_bounds
    elif extension in {".tif", ".tiff"}:
        with rasterio.open(input_path) as src:
            if src.crs is None:
                raise ValueError(f"Input raster has no CRS: {input_path}")

            minx, miny, maxx, maxy = transform_bounds(
                src.crs,
                "EPSG:4326",
                *src.bounds,
                densify_pts=21,
            )
    else:
        raise ValueError(
            "Unsupported input format. Please provide a .shp, .tif, "
            f"or .tiff file: {input_path}"
        )

    return {
        "minx": float(minx),
        "miny": float(miny),
        "maxx": float(maxx),
        "maxy": float(maxy),
    }


def parse_search_dates(start_date, end_date):
    """Validate YYYY-MM-DD dates and return an inclusive date range."""
    try:
        start = date.fromisoformat(start_date)
        end = date.fromisoformat(end_date)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Dates must use the YYYY-MM-DD format, for example 2025-11-15."
        ) from exc

    if start > end:
        raise ValueError("start_date cannot be later than end_date.")

    return start, end


def search_sentinel_1_product_names(extent_file, start_date, end_date):
    """
    Search Sentinel-1 IW GRDH products for an extent and inclusive date range.

    Parameters
    ----------
    extent_file : str
        Path to a .shp, .tif, or .tiff file defining the search extent.
    start_date : str
        First acquisition date to include, formatted as YYYY-MM-DD.
    end_date : str
        Last acquisition date to include, formatted as YYYY-MM-DD.

    Returns
    -------
    list[str]
        Unique Sentinel-1 product names without the trailing ``.SAFE``.
    """
    bounds = read_search_extent(extent_file)
    start, end = parse_search_dates(start_date, end_date)
    start_datetime = f"{start.isoformat()}T00:00:00.000Z"
    end_datetime = f"{(end + timedelta(days=1)).isoformat()}T00:00:00.000Z"

    polygon = (
        f"{bounds['minx']} {bounds['maxy']}, "
        f"{bounds['minx']} {bounds['miny']}, "
        f"{bounds['maxx']} {bounds['miny']}, "
        f"{bounds['maxx']} {bounds['maxy']}, "
        f"{bounds['minx']} {bounds['maxy']}"
    )

    query_url = (
        "https://catalogue.dataspace.copernicus.eu/odata/v1/Products?"
        f"$filter=OData.CSC.Intersects(area=geography'SRID=4326;"
        f"POLYGON(({polygon}))') "
        "and Collection/Name eq 'SENTINEL-1' "
        "and Attributes/OData.CSC.StringAttribute/any(att:att/Name eq "
        "'productType' and att/OData.CSC.StringAttribute/Value eq "
        f"'IW_GRDH_1S') and ContentDate/Start ge {start_datetime} "
        f"and ContentDate/Start lt {end_datetime}"
    )

    product_names = []
    while query_url:
        response = requests.get(query_url, timeout=60)
        response.raise_for_status()
        response_json = response.json()

        for product in response_json.get("value", []):
            name = product.get("Name")
            if name:
                product_names.append(name.removesuffix(".SAFE"))

        query_url = response_json.get("@odata.nextLink")

    return list(dict.fromkeys(product_names))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Return Sentinel-1 IW GRDH product names intersecting an input "
            "raster or shapefile during an inclusive date range."
        )
    )
    parser.add_argument(
        "extent_file",
        help="Path to the input .shp, .tif, or .tiff extent file.",
    )
    parser.add_argument(
        "start_date",
        help="First acquisition date to include (YYYY-MM-DD).",
    )
    parser.add_argument(
        "end_date",
        help="Last acquisition date to include (YYYY-MM-DD).",
    )
    args = parser.parse_args()

    names = search_sentinel_1_product_names(
        args.extent_file,
        args.start_date,
        args.end_date,
    )
    print(f"Search dates: {args.start_date} to {args.end_date}")
    print(f"Found {len(names)} Sentinel-1 products:")
    for name in names:
        print(name)
