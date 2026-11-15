"""Archive search cache, source fallback and year merge, without external services."""
import datetime
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile
from unittest.mock import Mock, patch

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.transform import from_origin
from requests.exceptions import HTTPError
from shapely.geometry import Point, box, mapping

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import S1_archive_search as archive  # noqa: E402


def _zip_bytes(files):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as stream:
        for name, data in files.items():
            stream.writestr(name, data)
    return buffer.getvalue()


def _missing():
    response = Mock(status_code=400, text='{"message":"No files were found for that service"}')
    return HTTPError("400 Client Error", response=response)


def _record(day, source, candidates):
    return {"version": archive.CACHE_VERSION, "date": day, "source": source,
            "candidates": candidates}


class ArchiveSearchTests(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.root)

    def test_cache_layout_is_year_month(self):
        directory, path = archive.day_paths(self.root, datetime.date(2020, 1, 5))
        self.assertEqual(Path(directory), Path(self.root) / "cache" / "2020" / "01")
        self.assertEqual(Path(path).name, "2020-01-05.json")

    def test_source_selection_fallback_cache_reuse_and_refresh(self):
        calls = []

        def fake_get(url):
            calls.append(url)
            if "RapidFloodMapping/20191102" in url:
                raise _missing()
            return Mock(content=_zip_bytes({"x.txt": "x"}))

        areas = Mock(return_value=([{"name": "S1A_area", "region": 1, "region_km2": 6000.0,
                                     "overlap_km2": 80.0}], {"regions": 1}))
        points = Mock(return_value=([{"name": "S1B_point", "point": "7", "tg": 3,
                                      "prob": 90.0, "dist_km": 0.0}], {"points": 5}))
        plans = {"areas": areas, "points": points}
        with patch.object(archive, "_get", side_effect=fake_get), \
                patch.object(archive, "search_areas", areas), \
                patch.object(archive, "search_points", points):
            self.assertEqual(archive.search_day(self.root, datetime.date(2019, 11, 1))["source"], "points")
            self.assertEqual(archive.search_day(self.root, datetime.date(2019, 11, 2))["source"], "points")
            self.assertEqual(archive.search_day(self.root, datetime.date(2019, 11, 3))["source"], "areas")
            self.assertEqual(plans["areas"].call_count, 1)
            self.assertEqual(plans["points"].call_count, 2)
            # Cached days are reused without any request; the zip is kept beside the JSON.
            requests = len(calls)
            archive.search_day(self.root, datetime.date(2019, 11, 3))
            self.assertEqual(len(calls), requests)
            directory, _ = archive.day_paths(self.root, datetime.date(2019, 11, 3))
            self.assertTrue(os.path.isfile(os.path.join(directory, "2019-11-03_rfm.zip")))
            # Refresh searches again but reuses the cached zip.
            archive.search_day(self.root, datetime.date(2019, 11, 3), refresh=True)
            self.assertEqual(len(calls), requests)
            self.assertEqual(plans["areas"].call_count, 2)

    def test_day_without_any_source_is_cached_as_none(self):
        with patch.object(archive, "_get", side_effect=lambda url: (_ for _ in ()).throw(_missing())):
            record = archive.search_day(self.root, datetime.date(2016, 1, 1))
        self.assertEqual((record["source"], record["candidates"]), ("none", []))
        self.assertIsNotNone(archive.read_day(self.root, datetime.date(2016, 1, 1)))

    def test_failed_day_is_not_cached_and_range_continues(self):
        with patch.object(archive, "search_day",
                          side_effect=[None, RuntimeError("429"), None]) as search:
            failed = archive.search_range(self.root, datetime.date(2020, 1, 1),
                                          datetime.date(2020, 1, 3))
        self.assertEqual(failed, ["2020-01-02"])
        self.assertEqual(search.call_count, 3)

    def test_select_names_applies_thresholds(self):
        areas = _record("2020-01-01", "areas", [
            {"name": "S1A_big", "region": 1, "region_km2": 6000, "overlap_km2": 50},
            {"name": "S1A_small_region", "region": 2, "region_km2": 4000, "overlap_km2": 900},
            {"name": "S1A_small_overlap", "region": 1, "region_km2": 6000, "overlap_km2": 49.9}])
        self.assertEqual(archive.select_names(areas), {"S1A_big"})
        self.assertEqual(archive.select_names(areas, min_region_km2=1000),
                         {"S1A_big", "S1A_small_region"})
        points = _record("2018-01-01", "points", [
            {"name": "S1A_ok", "point": "1", "tg": 3, "prob": 50, "dist_km": 10},
            {"name": "S1A_tg2", "point": "2", "tg": 2, "prob": 100, "dist_km": 0},
            {"name": "S1A_far", "point": "3", "tg": 3, "prob": 100, "dist_km": 11},
            {"name": "S1A_lowprob", "point": "4", "tg": 3, "prob": 49, "dist_km": 0}])
        self.assertEqual(archive.select_names(points), {"S1A_ok"})
        polarizations = _record("2020-01-01", "areas", [
            {"name": f"S1A_IW_GRDH_{mode}_X", "region": 1, "region_km2": 6000, "overlap_km2": 60}
            for mode in ("1SDV", "1SDH", "1SSV", "1SSH")])
        self.assertEqual(archive.select_names(polarizations),
                         {"S1A_IW_GRDH_1SDV_X", "S1A_IW_GRDH_1SDH_X"})
        self.assertEqual(archive.select_names(_record("2016-01-01", "none", [])), set())

    def test_merge_writes_sorted_unique_names_and_refuses_incomplete_years(self):
        def write(day, record):
            _, path = archive.day_paths(self.root, day)
            archive._atomic_write(path, json.dumps(record))

        day = datetime.date(2017, 1, 1)
        while day.year == 2017:
            write(day, _record(day.isoformat(), "none", []))
            day += datetime.timedelta(days=1)
        hit = {"name": "S1B_z", "point": "1", "tg": 3, "prob": 80, "dist_km": 0}
        write(datetime.date(2017, 3, 1), _record("2017-03-01", "points",
                                                 [hit, dict(hit, name="S1A_a"), hit]))
        self.assertEqual(archive.merge_year(self.root, 2017), (2, []))
        with open(os.path.join(self.root, "results", "2017.json")) as stream:
            self.assertEqual(json.load(stream), ["S1A_a", "S1B_z"])
        count, missing = archive.merge_year(self.root, 2018)
        self.assertIsNone(count)
        self.assertEqual(len(missing), 365)
        self.assertFalse(os.path.exists(os.path.join(self.root, "results", "2018.json")))

    def test_merge_cli_rejects_thresholds_below_cache_floor(self):
        with self.assertRaises(SystemExit):
            archive.main(["merge", "--root", self.root, "--years", "2017", "2017",
                          "--min-region-km2", "500"])

    def test_point_distances_use_footprints(self):
        footprint = box(10.0, 0.0, 11.0, 1.0)
        product = {"Id": "p", "Name": "S1A_scene.SAFE", "GeoFootprint": mapping(footprint)}
        with tempfile.TemporaryDirectory() as workdir:
            frame = gpd.GeoDataFrame({"PointID": [1, 2, 3], "ThresGroup": [3, 2, 1],
                                      "Prob": [90, 60, 100]},
                                     geometry=[Point(10.5, 0.5), Point(11.1, 0.5), Point(10.5, 0.5)],
                                     crs=4326)
            frame.to_file(os.path.join(workdir, "pts.shp"))
            zip_path = os.path.join(self.root, "pts.zip")
            with zipfile.ZipFile(zip_path, "w") as stream:
                for name in os.listdir(workdir):
                    stream.write(os.path.join(workdir, name), name)
            with tempfile.TemporaryDirectory() as extract, \
                    patch.object(archive, "catalogue_products", return_value={"p": product}):
                candidates, info = archive.search_points(zip_path, datetime.date(2017, 1, 1), extract)
        self.assertEqual(info, {"points": 3, "alert_points": 2})
        by_point = {c["point"]: c for c in candidates}
        self.assertEqual(set(by_point), {"1", "2"})
        self.assertEqual(by_point["1"]["dist_km"], 0.0)
        self.assertAlmostEqual(by_point["2"]["dist_km"], 11.1, delta=0.3)
        self.assertEqual(by_point["1"]["name"], "S1A_scene")

    def test_area_overlap_counts_only_this_region(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "regions.tif")
            data = np.zeros((100, 100), dtype="uint32")
            data[:, :50] = 1
            data[:, 50:] = 2
            transform = from_origin(0, 1, 0.01, 0.01)
            with rasterio.open(path, "w", driver="GTiff", height=100, width=100, count=1,
                               dtype="uint32", crs="EPSG:4326", transform=transform) as dst:
                dst.write(data, 1)
            feature = {"minx": 0, "miny": 0, "maxx": 1, "maxy": 1,
                       "warning_raster": path, "warning_region": 1}
            product = {"Id": "p", "GeoFootprint": mapping(box(0.25, 0.0, 0.75, 1.0))}
            # 25 columns x 100 rows of region 1, at (0.01 deg * 110 km)^2 per pixel.
            self.assertAlmostEqual(archive.overlap_km2(product, feature),
                                   2500 * (0.01 * 110) ** 2, places=3)


if __name__ == "__main__":
    unittest.main()
