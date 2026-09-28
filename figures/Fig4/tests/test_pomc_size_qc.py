"""Independent checks of the regional size rule and its auditability."""
import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/shared"))
from pomc_size_qc import filter_pomc_by_local_dapi


class RegionalSizeQC(unittest.TestCase):
    def test_each_region_uses_its_own_reference_and_preserves_ids(self):
        dapi = np.zeros((20, 20), np.int32)
        regions = np.ones_like(dapi); regions[:, 10:] = 2
        dapi[1:3, 1:3] = 1  # ARC area 4
        dapi[6:8, 1:5] = 2  # ARC area 8 -> regional mean 6
        dapi[1:5, 12:16] = 3  # ME area 16
        pomc = np.zeros_like(dapi)
        pomc[10:12, 1:4] = 10  # equality, keep
        pomc[14:15, 1:6] = 20  # area 5, reject
        pomc[10:13, 12:17] = 30  # area 15, reject in ME
        pomc[15:19, 12:16] = 40  # equality, keep
        original = pomc.copy()
        out, audit = filter_pomc_by_local_dapi(pomc, dapi, regions)
        self.assertEqual(set(np.unique(out)), {0, 10, 40})
        self.assertEqual(audit["objects_removed"], 2)
        np.testing.assert_array_equal(pomc, original)

    def test_full_object_area_is_not_clipped_to_region(self):
        dapi = np.zeros((10, 10), np.int32); dapi[1:3, 4:8] = 1
        regions = np.ones_like(dapi); regions[:, 6:] = 2
        pomc = np.zeros_like(dapi); pomc[5:7, 4:8] = 27
        out, audit = filter_pomc_by_local_dapi(pomc, dapi, regions)
        self.assertEqual(audit["objects"][0]["mean_dapi_area_px"], 8)
        self.assertTrue(np.any(out == 27))

    def test_missing_local_reference_does_not_invent_threshold(self):
        dapi = np.zeros((6, 6), np.int32)
        pomc = dapi.copy(); pomc[2, 2] = 5
        out, audit = filter_pomc_by_local_dapi(pomc, dapi, np.ones_like(dapi))
        self.assertEqual(out[2, 2], 5)
        self.assertEqual(audit["objects"][0]["decision"], "retain_no_local_dapi_reference")

    def test_mismatched_geometry_fails(self):
        with self.assertRaises(ValueError):
            filter_pomc_by_local_dapi(np.zeros((2, 2), int), np.zeros((3, 2), int), np.zeros((2, 2), int))


if __name__ == "__main__":
    unittest.main()
