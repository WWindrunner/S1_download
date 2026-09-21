"""Exercise search logic without requiring SNAP/GDAL or external services."""
import ast
import argparse
import contextlib
import datetime
import io
import json
import os
from pathlib import Path
import tempfile
import sys
from time import perf_counter
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
from shapely.geometry import box, Polygon, mapping, shape
import rasterio
from rasterio.features import rasterize, geometry_window
from rasterio.windows import Window
from rasterio.errors import WindowError
from rasterio.transform import from_origin
from shapely.wkt import loads


SOURCE = Path(__file__).resolve().parents[1] / "Sentinel_1_ESA_search_download_process_chain_v4.py"


class FloodSearchTests(unittest.TestCase):
    def setUp(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        self.ns = dict(datetime=datetime, os=os, np=np, pd=pd, requests=Mock(),
                       box=box, shape=shape, loads=loads, MIN_FLOOD_OVERLAP_KM2=50.0,
                       rasterio=rasterio, rasterize=rasterize, geometry_window=geometry_window,
                       Window=Window, WindowError=WindowError, mapping=mapping, perf_counter=perf_counter)
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(SOURCE), "exec"), self.ns)

    def test_pagination_and_midnight_boundaries(self):
        responses = [Mock(), Mock()]
        responses[0].json.return_value = {
            "value": [{"Id": "a", "Name": "scene-a"}], "@odata.nextLink": "https://example.test/page2"
        }
        responses[1].json.return_value = {"value": [{"Id": "b", "Name": "scene-b"}]}
        for response in responses:
            for product in response.json.return_value["value"]:
                product["GeoFootprint"] = mapping(box(0, 0, 1, 1))
        responses[0].json.return_value["value"].append({
            "Id": "too-small", "Name": "rejected", "GeoFootprint": mapping(box(0, 0, 0.001, 0.001))
        })
        self.ns["requests"].get.side_effect = responses
        result = self.ns["search_sentinel_with_shape_extent_and_data"](
            pd.DataFrame(columns=["Id", "Name"]), 2024, 2, 29,
            dict(minx=0, miny=0, maxx=1, maxy=1),
        )
        self.assertEqual(result.Id.tolist(), ["a", "b"])
        calls = self.ns["requests"].get.call_args_list
        self.assertIn("ContentDate/Start ge 2024-02-29T00:00:00.000Z", calls[0].args[0])
        self.assertIn("ContentDate/Start lt 2024-03-01T00:00:00.000Z", calls[0].args[0])
        self.assertEqual(calls[1].args[0], "https://example.test/page2")
        for response in responses:
            response.raise_for_status.assert_called_once()

    def test_original_polygon_excludes_bbox_and_expansion_only_matches(self):
        warning = Polygon([(0, 0), (2, 0), (0, 2)])
        feature = dict(minx=-1, miny=-1, maxx=3, maxy=3,
                       geometry=box(-1, -1, 3, 3), warning_wkt=warning.wkt)
        geometry = self.ns["_warning_geometry"](feature)
        check = self.ns["_has_flood_overlap"]
        self.assertFalse(check({"GeoFootprint": mapping(box(1.5, 1.5, 2, 2))}, geometry))
        self.assertFalse(check({"GeoFootprint": mapping(box(-0.5, -0.5, -0.1, -0.1))}, geometry))
        self.assertTrue(check({"GeoFootprint": mapping(box(0, 0, 0.2, 0.2))}, geometry))

    def test_overlap_threshold_holes_and_footprint_formats(self):
        warning = box(0, 0, 2, 2)
        check = self.ns["_has_flood_overlap"]
        for area, expected in [(49, False), (50, True), (51, True)]:
            footprint = box(0, 0, 1, area / 12100)
            self.assertEqual(check({"GeoFootprint": mapping(footprint)}, warning), expected)
            self.assertEqual(check({"Footprint": f"geography'SRID=4326;{footprint.wkt}'"}, warning), expected)
        with_hole = Polygon(warning.exterior.coords, [box(0.2, 0.2, 1.8, 1.8).exterior.coords])
        self.assertFalse(check({"GeoFootprint": mapping(box(0.5, 0.5, 1, 1))}, with_hole))
        self.assertFalse(check({"GeoFootprint": mapping(box(2, 0, 3, 1))}, warning))
        with self.assertRaisesRegex(ValueError, "Missing footprint"):
            check({"Id": "missing"}, warning)

    def test_inclusive_range_and_duplicate_regions(self):
        download = self.ns["download_flood_warning_shp_from_ESA"] = Mock()
        self.ns["simplify_flood_warning_shp_from_ESA"] = Mock(return_value=pd.DataFrame([{"x": 1}, {"x": 2}]))

        def search(df, year, month, day, feature):
            name = f"scene-{day}.SAFE" if day != 1 else "scene-1"
            return pd.concat([df, pd.DataFrame([{"Id": str(day), "Name": name}])])

        self.ns["search_sentinel_with_shape_extent_and_data"] = search
        with tempfile.TemporaryDirectory() as directory:
            result = self.ns["search_flood_images_by_date_range"]("2024-02-28", "2024-03-01", directory)
        self.assertEqual(result, {"names": ["scene-28", "scene-29", "scene-1"], "count": 3})
        self.assertEqual(download.call_count, 3)

    def write_raster(self, path, data, transform):
        with rasterio.open(path, "w", driver="GTiff", height=data.shape[0],
                           width=data.shape[1], count=1, dtype=data.dtype,
                           crs="EPSG:4326", transform=transform) as dst:
            dst.write(data, 1)

    def test_warning_region_raster_excludes_expansion_and_other_regions(self):
        transform = from_origin(0, 20 / 111, 1 / 111, 1 / 111)
        original = np.ones((20, 20), dtype="uint8")
        original[2:5, 2:5] = 0
        outer = box(0, 0, 20 / 111, 20 / 111)
        inner = box(5 / 111, 5 / 111, 15 / 111, 15 / 111)
        ring = Polygon(outer.exterior.coords, [inner.exterior.coords])
        regions = pd.DataFrame({"geometry": [ring, inner], "warning_region": [1, 2]})
        with tempfile.TemporaryDirectory() as directory:
            original_path = os.path.join(directory, "original.tif")
            output_path = os.path.join(directory, "labels.tif")
            self.write_raster(original_path, original, transform)
            self.ns["_write_warning_region_raster"](original_path, regions, output_path)
            with rasterio.open(output_path) as src:
                labels = src.read(1)
                self.assertEqual(src.transform, transform)
            self.assertTrue(np.all(labels[2:5, 2:5] == 0))
            self.assertTrue(np.all(labels[5:15, 5:15] == 2))
            feature = dict(minx=0, miny=0, maxx=20/111, maxy=20/111,
                           warning_raster=output_path, warning_region=1)
            product = {"GeoFootprint": mapping(inner)}
            self.assertFalse(self.ns["_has_raster_flood_overlap"](product, feature))
            feature["warning_region"] = 2
            self.assertTrue(self.ns["_has_raster_flood_overlap"](product, feature))

    def test_pixel_threshold_and_catalogue_integration(self):
        transform = from_origin(0, 20 / 111, 1 / 111, 1 / 111)
        footprint = box(0, 0, 20 / 111, 20 / 111)
        product = {"Id": "a", "Name": "scene.SAFE", "GeoFootprint": mapping(footprint)}
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "labels.tif")
            feature = dict(minx=0, miny=0, maxx=20/111, maxy=20/111,
                           warning_raster=path, warning_region=1)
            for count, expected in [(0, False), (50, False), (51, True)]:
                labels = np.zeros((20, 20), dtype="uint32")
                labels.flat[:count] = 1
                self.write_raster(path, labels, transform)
                self.assertEqual(self.ns["_has_raster_flood_overlap"](product, feature), expected)
            # Confirm the real search dispatch uses pixels, not legacy geometry.
            self.ns["_warning_geometry"] = Mock(side_effect=AssertionError("unexpected geometry path"))
            response = Mock()
            response.json.return_value = {"value": [product]}
            self.ns["requests"].get.return_value = response
            result = self.ns["search_sentinel_with_shape_extent_and_data"](
                pd.DataFrame(columns=["Id", "Name"]), 2025, 4, 1, feature)
            self.assertEqual(result.Id.tolist(), ["a"])
            self.assertFalse(self.ns["_has_raster_flood_overlap"](
                {"GeoFootprint": mapping(box(10, 10, 11, 11))}, feature))
            # A footprint hole excludes warning pixels, even with a large bbox.
            hole = Polygon(footprint.exterior.coords,
                           [box(0.0001, 0.0001, 19.99/111, 19.99/111).exterior.coords])
            self.assertFalse(self.ns["_has_raster_flood_overlap"](
                {"GeoFootprint": mapping(hole)}, feature))

    def test_raster_tiles_preserve_counts_across_boundary(self):
        transform = from_origin(0, 1 / 111, 1 / 111, 1 / 111)
        original = np.zeros((1, 540), dtype="uint8")
        original[0, 490:540] = 1
        original[0, 0] = 1
        footprint = box(0, 0, 540 / 111, 1 / 111)
        regions = pd.DataFrame({"geometry": [footprint], "warning_region": [1]})
        with tempfile.TemporaryDirectory() as directory:
            original_path = os.path.join(directory, "original.tif")
            output_path = os.path.join(directory, "labels.tif")
            self.write_raster(original_path, original, transform)
            self.ns["_write_warning_region_raster"](original_path, regions, output_path)
            with rasterio.open(output_path) as src:
                np.testing.assert_array_equal(src.read(1), original)
            feature = dict(minx=0, miny=0, maxx=540/111, maxy=1/111,
                           warning_raster=output_path, warning_region=1)
            self.assertTrue(self.ns["_has_raster_flood_overlap"](
                {"Footprint": f"geography'SRID=4326;{footprint.wkt}'"}, feature))

    def test_empty_invalid_and_failed_searches(self):
        self.ns["download_flood_warning_shp_from_ESA"] = Mock()
        self.ns["simplify_flood_warning_shp_from_ESA"] = Mock(return_value=pd.DataFrame())
        search = self.ns["search_flood_images_by_date_range"]
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(search("2024-01-01", "2024-01-01", directory), {"names": [], "count": 0})
            for start, end in [("2024-02-30", "2024-03-01"), ("2024-03-02", "2024-03-01")]:
                with self.assertRaises(ValueError):
                    search(start, end, directory)
            self.ns["download_flood_warning_shp_from_ESA"].side_effect = OSError("unavailable")
            with self.assertRaisesRegex(RuntimeError, "2024-01-01"):
                search("2024-01-01", "2024-01-02", directory)

    def test_nisar_search_uses_original_pixels_and_same_overlap_threshold(self):
        import NISAR_extent_time_download as nisar

        footprint = box(0, 0, 20 / 111, 20 / 111)
        candidate = Mock(properties={"sceneName": "NISAR_L2_PR_GCOV_TEST"})
        candidate.geojson.return_value = {"geometry": mapping(footprint)}
        rejected = Mock(properties={"sceneName": "NISAR_L2_PR_GCOV_OUTSIDE"})
        rejected.geojson.return_value = {"geometry": mapping(box(10, 10, 11, 11))}
        self.ns["_warning_geometry"] = Mock(side_effect=AssertionError("must use original pixels"))
        with tempfile.TemporaryDirectory() as directory, patch.object(
                nisar, "search_nisar", return_value=[candidate, rejected]) as query:
            path = os.path.join(directory, "labels.tif")
            feature = dict(minx=0, miny=0, maxx=20/111, maxy=20/111,
                           warning_raster=path, warning_region=1)
            for count, expected in ((50, []), (51, ["NISAR_L2_PR_GCOV_TEST"])):
                labels = np.zeros((20, 20), dtype="uint32")
                labels.flat[:count] = 1
                self.write_raster(path, labels, from_origin(0, 20/111, 1/111, 1/111))
                results = self.ns["search_nisar_with_shape_extent_and_data"](
                    pd.DataFrame(columns=["Id", "Name"]), 2026, 7, 7, feature)
                self.assertEqual(results.Name.tolist(), expected)
            self.assertEqual(query.call_args.args[:2], ("2026-07-07", "2026-07-07"))
            self.assertTrue(loads(query.call_args.args[2]).equals_exact(footprint, 1e-12))
            candidate.geojson.return_value = {"geometry": None}
            with self.assertRaisesRegex(ValueError, "Missing footprint"):
                self.ns["search_nisar_with_shape_extent_and_data"](
                    pd.DataFrame(columns=["Id", "Name"]), 2026, 7, 7, feature)
        self.ns["requests"].get.assert_not_called()

    def test_nisar_date_range_keeps_schema_and_deduplicates_regions_and_days(self):
        import NISAR_extent_time_download as nisar

        self.ns["download_flood_warning_shp_from_ESA"] = Mock()
        region = dict(minx=0, miny=0, maxx=1, maxy=1, geometry=box(0, 0, 1, 1))
        self.ns["simplify_flood_warning_shp_from_ESA"] = Mock(return_value=pd.DataFrame([region, region]))
        candidate = Mock(properties={"sceneName": "NISAR_L2_PR_GCOV_TEST"})
        candidate.geojson.return_value = {"geometry": mapping(box(0, 0, 1, 1))}
        self.ns["search_sentinel_with_shape_extent_and_data"] = Mock(side_effect=AssertionError("unexpected CDSE search"))
        with tempfile.TemporaryDirectory() as directory, patch.object(
                nisar, "search_nisar", return_value=[candidate]) as query:
            result = self.ns["search_flood_images_by_date_range"](
                "2026-07-07", "2026-07-08", directory, 10, 1000, sensor="nisar")
        self.assertEqual(result, {"names": ["NISAR_L2_PR_GCOV_TEST"], "count": 1})
        self.assertEqual(query.call_count, 4)
        self.assertEqual(query.call_args.args[:2], ("2026-07-08", "2026-07-08"))

    def configure_daily_main(self):
        self.ns.update(argparse=argparse, json=json, sys=sys, SCRIPT_DIR=str(SOURCE.parent))
        self.ns["download_flood_warning_shp_from_ESA"] = Mock()
        self.ns["simplify_flood_warning_shp_from_ESA"] = Mock(return_value=pd.DataFrame([{"region": 1}]))
        self.ns["run_new_processing_chain"] = Mock()

    def test_daily_main_default_and_nisar_routing_preserve_legacy_layout(self):
        for sensor in ("s1", "nisar"):
            with self.subTest(sensor=sensor), tempfile.TemporaryDirectory() as directory:
                self.configure_daily_main()
                name = "scene.SAFE" if sensor == "s1" else "NISAR_L2_PR_GCOV_TEST"
                data = pd.DataFrame([dict(Id="a", Name=name), dict(Id="a", Name=name)])
                s1 = self.ns["search_sentinel_with_shape_extent_and_data"] = Mock(return_value=data)
                nisar = self.ns["search_nisar_with_shape_extent_and_data"] = Mock(return_value=data)
                vrt = Path(directory) / "desert.vrt"
                vrt.touch()
                env = dict(RAPID_PROJECT_DIR=directory)
                argv = [str(SOURCE), "--desert-mask-vrt", str(vrt)]
                if sensor == "s1":
                    env.update(CDSE_USERNAME="user", CDSE_PASSWORD="password")
                else:
                    argv += ["--sensor", "nisar"]
                with patch.dict(os.environ, env, clear=True), patch.object(sys, "argv", argv):
                    self.assertEqual(self.ns["main"](), None if sensor == "s1" else 0)
                selected, other = (s1, nisar) if sensor == "s1" else (nisar, s1)
                selected.assert_called_once()
                other.assert_not_called()
                label = "Sentinel_1" if sensor == "s1" else "NISAR"
                now = datetime.datetime.now()
                day_root = Path(directory) / "data" / now.strftime("%Y-%m-%d")
                output = day_root / label / "processed_images"
                manifest = day_root / f"Processed_{label}_data_path_{now:%Y%m%d}.txt"
                self.assertEqual(manifest.read_text().strip(), str(output) + os.sep)
                process = self.ns["run_new_processing_chain"]
                if sensor == "s1":
                    process.assert_called_once_with(name, str(output), str(vrt), "cdse")
                else:
                    process.assert_called_once_with(name, str(output), str(vrt), "asf", sensor="nisar")
                    self.assertFalse((day_root / "Sentinel_1").exists())
                    self.assertEqual(self.ns["username"], "")

    def test_nisar_daily_failure_does_not_publish_success_path(self):
        self.configure_daily_main()
        self.ns["search_nisar_with_shape_extent_and_data"] = Mock(
            return_value=pd.DataFrame([dict(Id="a", Name="NISAR_L2_PR_GCOV_TEST")]))
        self.ns["run_new_processing_chain"].side_effect = RuntimeError("download failed")
        with tempfile.TemporaryDirectory() as directory:
            vrt = Path(directory) / "desert.vrt"
            vrt.touch()
            with patch.dict(os.environ, dict(RAPID_PROJECT_DIR=directory), clear=True), patch.object(
                    sys, "argv", [str(SOURCE), "--sensor", "nisar", "--desert-mask-vrt", str(vrt)]):
                self.assertEqual(self.ns["main"](), 1)
            self.assertEqual(list(Path(directory).rglob("Processed_NISAR_data_path_*.txt")), [])

    def test_nisar_search_only_needs_no_credentials_or_desert_vrt(self):
        self.configure_daily_main()
        search = self.ns["search_flood_images_by_date_range"] = Mock(return_value={"names": [], "count": 0})
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {}, clear=True):
            argv = [str(SOURCE), "--sensor", "nisar", "--search-only", "--start-date", "2026-07-07",
                    "--end-date", "2026-07-08", "--work-directory", directory]
            with patch.object(sys, "argv", argv):
                self.ns["main"]()
            search.assert_called_once_with("2026-07-07", "2026-07-08", directory,
                                           window_size=10, area_thresholds=1000, sensor="nisar")
            self.ns["run_new_processing_chain"].assert_not_called()
            search.reset_mock()
            with (patch.object(sys, "argv", argv + ["--download-source", "cdse"]),
                  contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit)):
                self.ns["main"]()
            search.assert_not_called()


if __name__ == "__main__":
    unittest.main()
