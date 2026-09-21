"""Spatial and UTC date inputs shared by the command-line workflows."""
from datetime import datetime, timedelta, timezone
from pathlib import Path


def date_interval(start, end):
    def parse(value):
        result = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if result.tzinfo is None:
            result = result.replace(tzinfo=timezone.utc)
        return result.astimezone(timezone.utc)

    lower, upper = parse(start), parse(end)
    if len(end) == 10:
        upper += timedelta(days=1)
    if lower >= upper:
        raise ValueError("Start must precede end; date-only end includes that entire UTC day")
    return lower, upper


def search_polygon(extent_file=None, bbox=None, polygon_wkt=None):
    from shapely.geometry import box
    from shapely.wkt import loads

    if polygon_wkt:
        geometry = loads(polygon_wkt)
    elif extent_file and str(extent_file).lower() != "none":
        path = Path(extent_file)
        if path.suffix.lower() in (".tif", ".tiff"):
            import rasterio
            from rasterio.warp import transform_bounds
            with rasterio.open(path) as src:
                if src.crs is None:
                    raise ValueError("Extent raster has no CRS")
                geometry = box(*transform_bounds(src.crs, "EPSG:4326", *src.bounds))
        else:
            import geopandas as gpd
            from shapely.ops import unary_union
            frame = gpd.read_file(path)
            if frame.empty or frame.crs is None:
                raise ValueError("Extent vector is empty or has no CRS")
            geometry = unary_union(frame.to_crs(4326).geometry)
    elif bbox is not None:
        west, south, east, north = map(float, bbox)
        if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
            raise ValueError("bbox must be WEST SOUTH EAST NORTH in EPSG:4326")
        geometry = box(west, south, east, north)
    else:
        raise ValueError("Supply --extent-file, --bbox, or --polygon-wkt")
    if geometry.is_empty or not geometry.is_valid or geometry.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Search geometry must be a valid Polygon or MultiPolygon")
    west, south, east, north = geometry.bounds
    if not (-180 <= west < east <= 180 and -90 <= south < north <= 90):
        raise ValueError("Search coordinates must be longitude/latitude")
    return geometry
