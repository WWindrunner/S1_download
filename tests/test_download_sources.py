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

from hyp3_sdk import exceptions as hyp3_exceptions

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
            self.assertEqual(downloader.main(["scene", folder, "user", "pass", "--download-source", "cdse"]), 0)
            download.assert_called_once()
            snap.assert_called_once_with(os.path.join(folder, "scene", "scene.SAFE.zip"), os.path.join(folder, "scene"))
            incidence.assert_called_once()
            snap.side_effect = RuntimeError("SNAP failed")
            self.assertEqual(downloader.main(["scene", folder, "user", "pass", "--download-source", "cdse"]), 1)

    def test_asf_cli_does_not_require_cdse_credentials(self):
        with patch.object(asf, "download_and_process", return_value={}) as process:
            # ASF is the default source.
            self.assertEqual(downloader.main(["one.SAFE,two", "data"]), 0)
            process.assert_called_once_with(["one", "two"], os.path.abspath("data"))
            process.return_value = {"two": "failed"}
            self.assertEqual(downloader.main(["one,two", "data", "--download-source", "asf"]), 1)
        with self.assertRaises(SystemExit):
            downloader.main(["one", "data", "--download-source", "cdse"])

    def test_archive_mapping_and_missing_layers(self):
        with tempfile.TemporaryDirectory() as folder:
            archive = Path(folder) / "product.zip"
            with zipfile.ZipFile(archive, "w") as package:
                for suffix in ("VV", "VH", "dem", "inc_map", "ls_map", "rgb"):
                    package.writestr(f"../../nested/product_{suffix}.tif", b"raster")
            products = asf.extract_products([archive], folder)
            self.assertEqual(set(products), {"vv", "vh", "dem", "angle", "ls"})
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
            for key in ("vv", "vh", "dem", "angle", "ls"):
                values = np.ones((2, 3), dtype="float32")
                if key == "angle":
                    values[:] = [[np.deg2rad(40), np.deg2rad(60), -9999]] * 2
                if key == "ls":  # second row: shadow (1 + 8 + 16) at a folded 40 degrees
                    values[:] = [[1, 1, 1], [25, 1, 1]]
                products[key] = root / (key + ".tif")
                with rasterio.open(products[key], "w", **dict(profile, nodata=0) if key == "ls" else profile) as target:
                    target.write(values, 1)
            # constant 60 degree ellipsoid incidence around the scene
            points = np.array([[-96, 41, 60], [-94, 41, 60], [-96, 39, 60], [-94, 39, 60]], dtype=float)
            asf.normalize_products(products, root, "scene", incidence_points=points)
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
                np.testing.assert_allclose(vv.read(1)[:2, :3], 1)  # not normalized here
                np.testing.assert_allclose(vh.read(1)[:2, :3], 1)
                np.testing.assert_allclose(dem.read(1)[:2, :3], 1)
                self.assertEqual(lia.dtypes, ("uint8",))
                self.assertEqual(lia.nodata, 255)
                np.testing.assert_array_equal(lia.read(1)[:2, :3], [[0, 1, 255], [1, 1, 255]])
                self.assertTrue(np.all(lia.read(1)[~np.isfinite(vv.read(1))] == 255))
            with rasterio.open(root / "incidenceAngleFromEllipsoid.tif") as inc:
                self.assertEqual(inc.transform, vv.transform)
                np.testing.assert_allclose(inc.read(1)[:2, :3], 60, rtol=1e-6)

    def fake_hyp3(self, outcomes, existing=()):
        """HyP3 mock: each granule's job runs one poll, then takes its outcome."""
        sdk = Mock()
        sdk.Batch = list
        hyp3 = sdk.HyP3.return_value
        events = []

        def make_job(granule, status="RUNNING", **params):
            job = Mock(job_id="job-" + granule, status_code=status, request_time=0,
                       job_parameters=dict(asf.RTC_OPTIONS, granules=[granule], **params))
            job.succeeded.side_effect = lambda: job.status_code == "SUCCEEDED"
            job.failed.side_effect = lambda: job.status_code == "FAILED"
            job.expired.return_value = False
            return job

        def submit(granule, **kwargs):
            events.append(("submit", granule))
            if outcomes[granule] == "submit-error":
                raise RuntimeError("submission rejected")
            return [make_job(granule)]

        def refresh(jobs):
            for job in jobs:
                job.status_code = outcomes[job.job_parameters["granules"][0]]
            return jobs

        hyp3.submit_rtc_job.side_effect = submit
        hyp3.refresh.side_effect = refresh
        hyp3.find_jobs.side_effect = lambda name, job_type: [j for j in existing
                                                              if j.job_parameters["granules"] == [name]]
        return sdk, events, make_job

    def test_asf_jobs_run_concurrently_and_failures_are_isolated(self):
        outcomes = {"a": "SUCCEEDED", "b": "FAILED", "c": "SUCCEEDED", "d": "submit-error", "e": "SUCCEEDED"}
        sdk, events, _ = self.fake_hyp3(outcomes)

        def convert(job, granule, product_dir):
            events.append(("download", granule))
            if granule == "c":
                raise ValueError("bad archive")

        with tempfile.TemporaryDirectory() as folder, patch.dict(sys.modules, {"hyp3_sdk": sdk, "hyp3_sdk.exceptions": hyp3_exceptions}), \
                patch.object(asf, "download_and_convert", side_effect=convert), patch("time.sleep"):
            failures = asf.download_and_process(list(outcomes), folder)
        # Every job is submitted before any download starts.
        self.assertEqual([e for e in events if e[0] == "submit"], [("submit", g) for g in outcomes])
        self.assertLess(max(i for i, e in enumerate(events) if e[0] == "submit"),
                        min(i for i, e in enumerate(events) if e[0] == "download"))
        self.assertEqual(sorted(g for kind, g in events if kind == "download"), ["a", "c", "e"])
        self.assertEqual(sorted(failures), ["b", "c", "d"])
        sdk.HyP3.return_value.submit_rtc_job.assert_any_call(
            granule="a", name="a", radiometry="gamma0", scale="power",
            resolution=20, include_inc_map=True, include_dem=True)

    def test_asf_rerun_skips_outputs_and_reuses_matching_jobs(self):
        sdk, events, make_job = self.fake_hyp3({"done": "SUCCEEDED", "running": "SUCCEEDED", "other": "SUCCEEDED"})
        existing = [make_job("running"), make_job("other", resolution=30)]
        sdk.HyP3.return_value.find_jobs.side_effect = lambda name, job_type: [
            j for j in existing if j.job_parameters["granules"] == [name]]
        with tempfile.TemporaryDirectory() as folder, patch.dict(sys.modules, {"hyp3_sdk": sdk, "hyp3_sdk.exceptions": hyp3_exceptions}), \
                patch.object(asf, "download_and_convert") as convert, patch("time.sleep"):
            (Path(folder) / "done").mkdir()
            for name in asf.output_names("done"):
                (Path(folder) / "done" / name).write_bytes(b"x")
            self.assertEqual(asf.download_and_process(["done", "running", "other"], folder), {})
        # Matching running job is reused; a job with different options is not.
        self.assertEqual(events, [("submit", "other")])
        self.assertEqual(sorted(c.args[1] for c in convert.call_args_list), ["other", "running"])

    def test_transient_submit_error_is_retried_without_duplicates(self):
        sdk, events, make_job = self.fake_hyp3({"flaky": "SUCCEEDED", "bad": "SUCCEEDED"})
        hyp3 = sdk.HyP3.return_value
        submitted = []

        def submit(granule, **kwargs):
            events.append(("submit", granule))
            if granule == "bad":
                raise hyp3_exceptions.HyP3Error("400 invalid granule")
            if len(events) == 1:
                # The first request reaches HyP3 but the gateway times out.
                submitted.append(make_job(granule))
                raise hyp3_exceptions.ServerError()
            raise AssertionError("must reuse the job created by the timed-out request")

        hyp3.submit_rtc_job.side_effect = submit
        hyp3.find_jobs.side_effect = lambda name, job_type: [j for j in submitted
                                                              if j.job_parameters["granules"] == [name]]
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(sys.modules, {"hyp3_sdk": sdk, "hyp3_sdk.exceptions": hyp3_exceptions}), \
                patch.object(asf, "download_and_convert") as convert, patch("time.sleep"):
            failures = asf.download_and_process(["flaky", "bad"], folder)
        # Client errors are not retried; the flaky granule is processed once.
        self.assertEqual(events, [("submit", "flaky"), ("submit", "bad")])
        self.assertEqual(list(failures), ["bad"])
        self.assertEqual([c.args[1] for c in convert.call_args_list], ["flaky"])

    def test_drivers_default_sentinel_to_asf(self):
        import SAR_download_process as driver
        parser = driver.build_parser()
        args = parser.parse_args(["--names", "scene", "--search-only"])
        driver.validate_args(parser, args)
        self.assertEqual(args.download_source, "asf")
        text = (ROOT / "Sentinel_1_ESA_search_download_process_chain_v4.py").read_text(encoding="utf-8")
        self.assertIn('choices=("cdse", "asf"), default="asf"', text)
        for script in ("execute.sh", "execute_flood_warning.sh"):
            self.assertNotIn("DOWNLOAD_SOURCE:-cdse", (ROOT / script).read_text(encoding="utf-8"))

    def test_driver_asf_failure_does_not_stop_other_scenes(self):
        import SAR_download_process as driver
        with tempfile.TemporaryDirectory() as folder:
            vrt = Path(folder) / "desert.vrt"
            vrt.touch()
            with patch("Sentinel_1_specific_name_ASF.download_and_process",
                       return_value={"bad": "ASF RTC job failed"}) as batch, \
                    patch.object(driver, "process_sentinel") as scene:
                status = driver.main(["--names", "good", "bad", "good2", "--output-dir", folder,
                                      "--desert-mask-vrt", str(vrt)])
            batch.assert_called_once()
            self.assertEqual(status, 1)
            self.assertEqual([c.args[0] for c in scene.call_args_list], ["good", "good2"])

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
                    ns["run_new_processing_chain"]("scene.SAFE", folder, "desert.vrt", "cdse")
                else:
                    ns["run_new_processing_chain"]("scene.SAFE", folder, "desert.vrt", channel)
                commands = [call.args[0] for call in runner.run.call_args_list]
                self.assertEqual(commands[0][-2:], ["--download-source", channel])
                scripts = [Path(command[1]).name for command in commands]
                self.assertEqual("cal_LIA.py" in scripts, channel == "cdse")
                self.assertIn("Snow_detect.py", scripts)
                self.assertTrue((scene / "scene_DEM.tif").exists())

    def test_daily_nisar_adapter_uses_gcov_chain_and_retains_inputs(self):
        tree = ast.parse((ROOT / "Sentinel_1_ESA_search_download_process_chain_v4.py").read_text(encoding="utf-8"))
        function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "run_new_processing_chain")
        runner = Mock()
        ns = dict(os=os, sys=sys, subprocess=runner, SCRIPT_DIR=str(ROOT))
        exec(compile(ast.Module(body=[function], type_ignores=[]), "chain", "exec"), ns)
        with tempfile.TemporaryDirectory() as directory:
            name = "NISAR_L2_PR_GCOV_TEST"
            scene = Path(directory) / name
            scene.mkdir()
            original = scene / (name + ".h5")
            original.touch()
            # Omitted download_source selects ASF; no CDSE credential globals are needed.
            ns["run_new_processing_chain"](name, directory, "desert.vrt", sensor="nisar")
            command = runner.run.call_args.args[0]
            self.assertEqual(Path(command[1]).name, "NISAR_specific_name_download_process.py")
            self.assertEqual(command[2:], ["--names", name, "--output-dir", directory,
                                          "--desert-mask-vrt", "desert.vrt", "--download-source", "asf"])
            self.assertTrue(runner.run.call_args.kwargs["check"])
            self.assertTrue(original.exists())
            with self.assertRaises(ValueError):
                ns["run_new_processing_chain"](name, directory, "desert.vrt", "cdse", sensor="nisar")
            runner.run.assert_called_once()


if __name__ == "__main__":
    unittest.main()
