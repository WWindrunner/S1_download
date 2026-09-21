"""Offline numerical and routing checks using small, real HDF5/GeoTIFF fixtures."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

import h5py
import numpy as np
import rasterio
from rasterio.transform import from_origin

from Desert_mask import generate_desert_mask
from NISAR_GCOV_h5_2_tif import convert_gcov, source_grid
from NISAR_incidence_angle_interpolate import cube_interpolators, interpolate_angles, local_incidence
from NISAR_extent_time_download import search_nisar, download_product
from NISAR_specific_name_download_process import process_scene
from nisar_common import download_http, local_path
from sar_dem import prepare_dem
from SAR_download_process import main

SCENE = "NISAR_L2_PR_GCOV_024_156_D_069_2005_QPDH_A_20260707T004257_20260707T004331_P05023_N_F_J_001"


def write_tif(path, values, transform=None, tags=None, nodata=np.nan):
    values = np.asarray(values, dtype=np.float32)
    with rasterio.open(path, "w", driver="GTiff", crs="EPSG:4326", count=1,
                       width=values.shape[1], height=values.shape[0], dtype="float32",
                       transform=transform or from_origin(-80, 25, .0002, .0002), nodata=nodata) as dst:
        dst.write(values, 1)
        if tags:
            dst.update_tags(**tags)


def make_h5(path):
    with h5py.File(path, "w") as handle:
        grid = handle.create_group("science/LSAR/GCOV/grids/frequencyA")
        x = -80 + (np.arange(6) + .5) * .0002
        y = 25 - (np.arange(5) + .5) * .0002
        for key, data in dict(xCoordinates=x, yCoordinates=y, xCoordinateSpacing=.0002,
                              yCoordinateSpacing=-.0002, projection=4326).items():
            grid[key] = data
        for pol, power in (("HHHH", 2), ("HVHV", 0.25)):
            grid[pol] = np.full((5, 6), power, dtype=np.float32)
            grid[pol].attrs["_FillValue"] = np.nan
        mask = np.ones((5, 6), dtype=np.uint8)
        mask[0, 0], mask[0, 1] = 0, 255
        grid["mask"] = mask
        grid["rtcGammaToSigmaFactor"] = 7.0  # Must not be multiplied into gamma0.
        identification = handle.create_group("science/LSAR/identification")
        identification["productType"] = np.bytes_("GCOV")
        identification["zeroDopplerStartTime"] = np.bytes_("2026-07-07T00:42:57Z")
        radar = handle.create_group("science/LSAR/GCOV/metadata/radarGrid")
        # Deliberately different horizontal extent/spacing from the image grid.
        radar["xCoordinates"] = [-80.01, -79.99]
        radar["yCoordinates"] = [25.01, 24.99]
        radar["heightAboveEllipsoid"] = [0., 100., 200.]
        radar["projection"] = 4326
        z, yy, xx = np.meshgrid(radar["heightAboveEllipsoid"][:], radar["yCoordinates"][:],
                               radar["xCoordinates"][:], indexing="ij")
        angles = 40 + .01 * z + 100 * (xx + 80) + 20 * (yy - 25)
        radar["incidenceAngle"] = angles.astype(np.float32)
        radar["incidenceAngle"].attrs["units"] = np.bytes_("degrees")
        radar["losUnitVectorX"] = np.sin(np.deg2rad(angles)).astype(np.float32)
        radar["losUnitVectorY"] = np.zeros_like(angles, dtype=np.float32)


class NisarTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.temp.name = str(local_path(self.temp.name))
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.h5 = self.root / (SCENE + ".h5")
        make_h5(self.h5)

    def converted(self, **kwargs):
        return convert_gcov(self.h5, self.root / "converted", **kwargs)

    def test_pixel_centres_linear_power_mask_and_common_grid(self):
        with h5py.File(self.h5) as handle:
            transform, width, height = source_grid(handle["science/LSAR/GCOV/grids/frequencyA"])
        self.assertEqual((width, height), (6, 5))
        self.assertAlmostEqual(transform.c, -80)
        self.assertAlmostEqual(transform.f, 25)
        outputs = self.converted()
        with rasterio.open(outputs["HH"]) as hh, rasterio.open(outputs["HV"]) as hv:
            self.assertTrue(hh.transform.almost_equals(transform))
            self.assertEqual(hh.transform, hv.transform)
            self.assertEqual(hh.shape, (5, 6))
            self.assertTrue(np.isnan(hh.read(1)[0, :2]).all())
            np.testing.assert_allclose(hh.read(1)[1:], 2)
            np.testing.assert_allclose(hv.read(1)[1:], .25)
        with self.assertRaises(ValueError):
            self.converted(polarizations=["VV"])

    def test_explicit_resolution(self):
        outputs = self.converted(resolution=30)
        with rasterio.open(outputs["HH"]) as src:
            np.testing.assert_allclose(src.res, (30 / 111320, 30 / 111320))

    def test_cube_uses_own_coordinates_exact_levels_and_no_extrapolation(self):
        with h5py.File(self.h5) as handle:
            fields, crs = cube_interpolators(handle["science/LSAR/GCOV"])
        self.assertEqual(crs.to_epsg(), 4326)
        values = fields["incidenceAngle"]([[0, 25, -80], [100, 25, -80], [200, 25, -80],
                                            [-1, 25, -80], [201, 25, -80], [100, 26, -80]])
        np.testing.assert_allclose(values[:3], [40, 41, 42])
        self.assertTrue(np.isnan(values[3:]).all())

    def test_egm2008_conversion_and_nodata(self):
        reference = self.root / "ref.tif"
        dem = self.root / "orthometric.tif"
        geoid = self.root / "geoid.tif"
        write_tif(reference, np.ones((5, 6)))
        heights = np.full((5, 6), 70., dtype=np.float32)
        heights[2, 3] = np.nan
        write_tif(dem, heights)
        write_tif(geoid, np.full((5, 6), 30.))
        target = prepare_dem(reference, self.root / "dem.tif", self.root, dem, "egm2008", geoid)
        with rasterio.open(target) as src:
            self.assertEqual(src.tags()["height_reference"], "WGS84_ellipsoid")
            expected = heights + 30
            np.testing.assert_allclose(src.read(1), expected, equal_nan=True)

    def test_flat_lia_and_reference_nodata(self):
        reference = self.converted()["HH"]
        dem = self.root / "ellipsoid.tif"
        write_tif(dem, np.full((5, 6), 100.), tags={"height_reference": "WGS84_ellipsoid"})
        outputs = interpolate_angles(self.h5, reference, dem, self.root, SCENE, threshold=41)
        with (rasterio.open(outputs["incidence"]) as inc, rasterio.open(outputs["local_angle"]) as lia,
              rasterio.open(outputs["lia"]) as mask):
            np.testing.assert_allclose(lia.read(1), inc.read(1), equal_nan=True, atol=1e-4)
            self.assertTrue((mask.read(1)[0, :2] == 255).all())
            valid = np.isfinite(lia.read(1))
            np.testing.assert_array_equal(mask.read(1)[valid], (lia.read(1)[valid] > 41).astype(np.uint8))

    def test_lia_slope_direction(self):
        transform = from_origin(-80, 25, .0002, .0002)
        rows, cols = 5, 6
        lat = np.deg2rad(transform.f + (np.arange(rows) + .5) * transform.e)
        metres_lon = 111412.84 * np.cos(lat) - 93.5 * np.cos(3 * lat) + .118 * np.cos(5 * lat)
        # Slope rises eastwards; an eastward sensor sees a larger local angle.
        dem = np.tile(np.arange(cols), (rows, 1)) * (transform.a * metres_lon[2]) * np.tan(np.deg2rad(10))
        inc = np.full((rows, cols), 40., dtype=np.float32)
        east = np.full_like(inc, np.sin(np.deg2rad(40)))
        north = np.zeros_like(inc)
        angle = local_incidence(dem, transform, inc, east, north)
        np.testing.assert_allclose(angle[2], 50, atol=.001)
        angle = local_incidence(-dem, transform, inc, east, north)
        np.testing.assert_allclose(angle[2], 30, atol=.001)
        # Northward slope uses a negative row spacing (north-up raster).
        metres_lat = 111132.92 - 559.82 * np.cos(2*lat[2]) + 1.175*np.cos(4*lat[2]) - .0023*np.cos(6*lat[2])
        dem = np.tile(np.arange(rows)[:, None], (1, cols)) * transform.e * metres_lat * np.tan(np.deg2rad(10))
        angle = local_incidence(dem, transform, inc, north, east)
        np.testing.assert_allclose(angle[2], 50, atol=.001)

    def test_desert_preserves_classes_and_missing_source_fails_strict(self):
        reference = self.converted()["HH"]
        source = self.root / "desert.tif"
        write_tif(source, [[0, 1, 2], [2, 1, 0]], transform=from_origin(-80, 25, .0004, .0005), nodata=-9999)
        result = generate_desert_mask(SCENE, reference, self.root, source, strict=True)
        with rasterio.open(result) as src:
            self.assertTrue(set(np.unique(src.read(1))).issubset({-9999., 0., 1., 2.}))
        with self.assertRaises(Exception):
            generate_desert_mask(SCENE, reference, self.root, self.root / "missing.tif", strict=True)

    def test_workflow_resume_changed_threshold_and_failed_stage(self):
        import Snow_detect
        dem = self.root / "dem.tif"
        desert = self.root / "desert.tif"
        write_tif(dem, np.full((5, 6), 100.))
        write_tif(desert, np.zeros((5, 6)))

        def snow(name, reference, directory, *args, **kwargs):
            paths = []
            with rasterio.open(reference) as ref:
                for suffix in ("ice", "cloud"):
                    path = Path(directory) / f"{name}_{suffix}.tif"
                    write_tif(path, np.zeros(ref.shape), ref.transform)
                    paths.append(str(path))
            return paths

        kwargs = dict(h5_file=self.h5, output_root=self.root / "out", desert_mask_vrt=desert,
                      dem=dem, dem_height_reference="ellipsoid")
        with patch.object(Snow_detect, "generate_s2_masks", side_effect=RuntimeError("network")):
            with self.assertRaises(RuntimeError):
                process_scene(**kwargs)
        folder = self.root / "out" / SCENE
        state = json.loads((folder / "processing_status.json").read_text())
        self.assertEqual(state["stages"]["snow_cloud"]["status"], "failed")
        self.assertFalse(state["complete"])
        gamma_stamp = (folder / "Gamma0_HH.tif").stat().st_mtime_ns
        with patch.object(Snow_detect, "generate_s2_masks", side_effect=snow) as mocked:
            process_scene(**kwargs)
            self.assertEqual(mocked.call_count, 1)
            process_scene(**kwargs)
            self.assertEqual(mocked.call_count, 1)
            process_scene(**kwargs, lia_threshold=30)
            self.assertEqual(mocked.call_count, 1)
        self.assertEqual(gamma_stamp, (folder / "Gamma0_HH.tif").stat().st_mtime_ns)
        with rasterio.open(folder / f"{SCENE}_LIA.tif") as mask:
            self.assertEqual(mask.tags()["lia_threshold_degrees"], "30")
        self.assertTrue(json.loads((folder / "processing_status.json").read_text())["complete"])

    def test_search_deduplicates_and_excludes_end_boundary(self):
        def product(suffix, acquired):
            return Mock(properties={"sceneName": SCENE + suffix, "startTime": acquired})
        first = product("", "2026-07-07T00:00:00Z")
        asf = Mock()
        asf.search.return_value = [first, first, product("_next", "2026-07-08T00:00:00Z")]
        with patch.dict(sys.modules, asf_search=asf):
            results = search_nisar("2026-07-07", "2026-07-07", "POLYGON((-80 25,-79 25,-79 26,-80 25))")
        self.assertEqual(results, [first])
        self.assertEqual(asf.search.call_args.kwargs["processingLevel"], "GCOV")

    def test_download_size_mismatch_does_not_publish(self):
        response = Mock()
        response.__enter__ = Mock(return_value=response)
        response.__exit__ = Mock(return_value=False)
        response.headers = {"Content-Length": "10"}
        response.iter_content.return_value = [b"short"]
        session = Mock()
        session.get.return_value = response
        target = self.root / "download.h5"
        with self.assertRaises(IOError):
            download_http("https://example.test/data.h5", target, session)
        self.assertFalse(target.exists())
        self.assertTrue(target.with_name("download.h5.part").exists())

    def test_complete_h5_is_reused_and_primary_asset_selected(self):
        target_dir = self.root / "downloads" / SCENE
        target_dir.mkdir(parents=True)
        target = target_dir / self.h5.name
        target.write_bytes(self.h5.read_bytes())
        product = Mock(properties=dict(sceneName=SCENE, url=f"https://example.test/{self.h5.name}",
                                       additionalUrls=["https://example.test/metadata.xml"],
                                       bytes={self.h5.name: {"bytes": target.stat().st_size}}))
        asf = Mock()
        with patch.dict(sys.modules, asf_search=asf):
            self.assertEqual(download_product(product, self.root / "downloads"), target)
        asf.ASFSession.assert_not_called()

    def test_invalid_download_is_not_published_as_h5(self):
        product = Mock(properties=dict(sceneName=SCENE, url=f"https://example.test/{self.h5.name}"))
        asf = Mock()
        session = asf.ASFSession.return_value
        session.__enter__ = Mock(return_value=session)
        session.__exit__ = Mock(return_value=False)

        def bad_download(url, path, **kwargs):
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"not an HDF5 product")

        with patch.dict(sys.modules, asf_search=asf), patch(
                "NISAR_extent_time_download.download_http", side_effect=bad_download):
            with self.assertRaises(OSError):
                download_product(product, self.root / "downloads")
        self.assertFalse((self.root / "downloads" / SCENE / self.h5.name).exists())

    def test_snow_explicit_acquisition_time_and_no_observations(self):
        import Snow_detect
        reference = self.converted()["HH"]
        eodag = Mock()
        eodag.EODataAccessGateway.return_value.search_all.return_value = []
        osgeo = Mock()
        with patch.dict(sys.modules, eodag=eodag, osgeo=osgeo):
            paths = Snow_detect.generate_s2_masks("opaque_scene", reference, self.root,
                                                  "2026-07-07T01:00:00+01:00", strict=True)
        search = eodag.EODataAccessGateway.return_value.search_all.call_args.kwargs
        self.assertEqual(search["end"], "2026-07-07")
        self.assertEqual(search["start"], "2026-06-22")
        for path in paths:
            with rasterio.open(path) as src:
                self.assertTrue((src.read(1) == 255).all())

    def test_batch_continues_and_returns_failure(self):
        import NISAR_extent_time_download as downloads
        import NISAR_specific_name_download_process as processing
        desert = self.root / "desert.tif"
        write_tif(desert, np.zeros((5, 6)))
        with (patch.object(downloads, "find_product", side_effect=[RuntimeError("missing"), Mock()]) as find,
              patch.object(downloads, "download_product", return_value=self.h5),
              patch.object(processing, "process_scene") as process,
              contextlib.redirect_stderr(io.StringIO())):
            code = main(["--sensor", "nisar", "--names", SCENE + "_missing", SCENE,
                         "--output-dir", str(self.root), "--desert-mask-vrt", str(desert)])
        self.assertEqual(code, 1)
        self.assertEqual(find.call_count, 2)
        process.assert_called_once()

    def test_cli_rejects_unsupported_combinations_before_network(self):
        for args in (["--sensor", "nisar", "--download-source", "cdse"],
                     ["--sensor", "nisar", "--product-type", "GUNW"],
                     ["--sensor", "s1", "--resolution", "30"]):
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                main([*args, "--names", SCENE, "--search-only"])
        report = self.root / "selection.json"
        self.assertEqual(main(["--sensor", "nisar", "--names", SCENE, "--search-only", "--output-json", str(report)]), 0)
        self.assertEqual(json.loads(report.read_text())["count"], 1)


if __name__ == "__main__":
    unittest.main()
