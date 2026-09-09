"""Exercise search logic without requiring SNAP/GDAL or external services."""
import ast
import datetime
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import Mock

import numpy as np
import pandas as pd
from shapely.geometry import box, Polygon, mapping, shape
from shapely.wkt import loads


SOURCE = Path(__file__).resolve().parents[1] / "Sentinel_1_ESA_search_download_process_chain_v4.py"


class FloodSearchTests(unittest.TestCase):
    def setUp(self):
        tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
        self.ns = dict(datetime=datetime, os=os, np=np, pd=pd, requests=Mock(),
                       box=box, shape=shape, loads=loads, MIN_FLOOD_OVERLAP_KM2=50.0)
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
            return pd.concat([df, pd.DataFrame([{"Id": str(day), "Name": f"scene-{day}"}])])

        self.ns["search_sentinel_with_shape_extent_and_data"] = search
        with tempfile.TemporaryDirectory() as directory:
            result = self.ns["search_flood_images_by_date_range"]("2024-02-28", "2024-03-01", directory)
        self.assertEqual(result["ids"], ["28", "29", "1"])
        self.assertEqual(result["count"], 3)
        self.assertEqual(len(result["images"]), 3)
        self.assertEqual(download.call_count, 3)

    def test_empty_invalid_and_failed_searches(self):
        self.ns["download_flood_warning_shp_from_ESA"] = Mock()
        self.ns["simplify_flood_warning_shp_from_ESA"] = Mock(return_value=pd.DataFrame())
        search = self.ns["search_flood_images_by_date_range"]
        with tempfile.TemporaryDirectory() as directory:
            self.assertEqual(search("2024-01-01", "2024-01-01", directory), {"ids": [], "count": 0, "images": []})
            for start, end in [("2024-02-30", "2024-03-01"), ("2024-03-02", "2024-03-01")]:
                with self.assertRaises(ValueError):
                    search(start, end, directory)
            self.ns["download_flood_warning_shp_from_ESA"].side_effect = OSError("unavailable")
            with self.assertRaisesRegex(RuntimeError, "2024-01-01"):
                search("2024-01-01", "2024-01-02", directory)


if __name__ == "__main__":
    unittest.main()
