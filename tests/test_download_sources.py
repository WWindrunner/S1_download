"""Offline integration checks for download routing and the ASF raster contract."""
import ast
import glob
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch
import zipfile

import numpy as np
import rasterio
from rasterio.transform import from_origin

import Sentinel_1_specific_name_ASF as asf
import Sentinel_1_specific_name_download_process as downloader

ROOT = Path(__file__).resolve().parents[1]


class DownloadTests(unittest.TestCase):
    def test_legacy_cli_still_runs_snap(self):
        frame = Mock()
        frame.__len__ = Mock(return_value=1)
        frame.iterrows.return_value = [(0, {"Id": "id", "Name": "scene.SAFE"})]
        with tempfile.TemporaryDirectory() as folder, patch.object(
            downloader, "search_sentinel_with_S1name", return_value=frame
        ), patch.object(downloader, "log_in", return_value="token"), patch.object(
            downloader, "download_Sentinel_with_ids_names"
        ) as download, patch.object(downloader, "process_snentinel_images") as snap, patch.object(
            downloader, "incidence_process"
        ) as incidence:
            self.assertEqual(downloader.main(["scene", folder, "user", "pass"]), 0)
            download.assert_called_once()
            snap.assert_called_once_with(os.path.join(folder, "scene", "scene.SAFE.zip"), os.path.join(folder, "scene"))
            incidence.assert_called_once()
            snap.side_effect = RuntimeError("SNAP failed")
            self.assertEqual(downloader.main(["scene", folder, "user", "pass"]), 1)

    def test_asf_cli_does_not_require_cdse_credentials(self):
        with patch.object(asf, "download_and_process") as process:
            self.assertEqual(downloader.main(["one.SAFE,two", "data", "--download-source", "asf"]), 0)
            process.assert_called_once_with(["one", "two"], os.path.abspath("data"))
        with self.assertRaises(SystemExit):
            downloader.main(["one", "data"])

    def test_archive_mapping_and_missing_layers(self):
        with tempfile.TemporaryDirectory() as folder:
            archive = Path(folder) / "product.zip"
            with zipfile.ZipFile(archive, "w") as package:
                for suffix in ("VV", "VH", "dem", "inc_map", "ls_map", "rgb"):
                    package.writestr(f"../../nested/product_{suffix}.tif", b"raster")
            products = asf.extract_products([archive], folder)
            self.assertEqual(set(products), {"vv", "vh", "dem", "angle"})
            self.assertTrue(all(p.parent == Path(folder) for p in products.values()))
            with self.assertRaises(ValueError):
                asf.extract_products([archive, archive], folder)
            with self.assertRaises(ValueError):
                asf.extract_products([], folder)

    def test_raster_interface_and_radian_threshold(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            spacing = 20 / 111320.0
            profile = dict(driver="GTiff", width=3, height=2, count=1,
                           dtype="float32", crs="EPSG:4326",
                           transform=from_origin(-95, 40, spacing, spacing), nodata=-9999)
            products = {}
            for key in ("vv", "vh", "dem", "angle"):
                values = np.ones((2, 3), dtype="float32")
                if key == "angle":
                    values[:] = [[np.deg2rad(40), np.deg2rad(60), -9999]] * 2
                products[key] = root / (key + ".tif")
                with rasterio.open(products[key], "w", **profile) as target:
                    target.write(values, 1)
            asf.normalize_products(products, root, "scene")
            with rasterio.open(root / "Gamma0_VV.tif") as vv, rasterio.open(
                root / "Gamma0_VH.tif"
            ) as vh, rasterio.open(root / "scene_DEM.tif") as dem, rasterio.open(
                root / "scene_LIA.tif"
            ) as lia:
                for raster in (vh, dem, lia):
                    self.assertEqual(raster.transform, vv.transform)
                    self.assertEqual(raster.shape, vv.shape)
                    self.assertEqual(raster.crs.to_epsg(), 4326)
                self.assertEqual(vv.dtypes, ("float32",))
                self.assertTrue(np.isnan(vv.nodata))
                np.testing.assert_allclose(vv.read(1)[:2, :3], 1)
                self.assertEqual(lia.dtypes, ("uint8",))
                self.assertEqual(lia.nodata, 255)
                np.testing.assert_array_equal(lia.read(1)[:2, :3], [[0, 1, 255]] * 2)
                self.assertTrue(np.all(lia.read(1)[~np.isfinite(vv.read(1))] == 255))

    def test_failed_hyp3_job_stops_before_download(self):
        sdk = Mock()
        sdk.HyP3.return_value.watch.return_value = [Mock(succeeded=False)]
        with tempfile.TemporaryDirectory() as folder, patch.dict(sys.modules, hyp3_sdk=sdk):
            with self.assertRaises(RuntimeError):
                asf.download_and_process(["scene"], folder)
            sdk.HyP3.return_value.submit_rtc_job.assert_called_once_with(
                granule="scene", name="scene", radiometry="gamma0", scale="power",
                resolution=20, include_inc_map=True, include_dem=True)

    def test_full_chain_routes_ancillary_steps(self):
        tree = ast.parse((ROOT / "Sentinel_1_ESA_search_download_process_chain_v4.py").read_text(encoding="utf-8"))
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "run_new_processing_chain")
        for channel in ("cdse", "asf"):
            with self.subTest(channel=channel), tempfile.TemporaryDirectory() as folder:
                scene = Path(folder) / "scene"
                scene.mkdir()
                for name in ("Gamma0_VV.tif", "Gamma0_VH.tif", "scene_desert.tif", "scene_LIA.tif",
                             "scene_ice.tif", "scene_cloud.tif", "scene_DEM.tif"):
                    (scene / name).touch()
                if channel == "cdse":
                    (scene / "scene_incidenceAngleFromEllipsoid.tif").touch()
                runner = Mock()
                ns = dict(os=os, glob=glob, shutil=shutil, sys=sys, subprocess=runner,
                          SCRIPT_DIR=str(ROOT), username="", password="")
                exec(compile(ast.Module(body=[function], type_ignores=[]), "chain", "exec"), ns)
                if channel == "cdse":
                    ns["run_new_processing_chain"]("scene.SAFE", folder, "desert.vrt")
                else:
                    ns["run_new_processing_chain"]("scene.SAFE", folder, "desert.vrt", channel)
                commands = [call.args[0] for call in runner.run.call_args_list]
                self.assertEqual(commands[0][-2:], ["--download-source", channel])
                scripts = [Path(command[1]).name for command in commands]
                self.assertEqual("cal_LIA.py" in scripts, channel == "cdse")
                self.assertIn("Snow_detect.py", scripts)
                self.assertTrue((scene / "scene_DEM.tif").exists())


if __name__ == "__main__":
    unittest.main()
