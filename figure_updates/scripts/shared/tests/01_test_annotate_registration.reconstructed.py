#!/usr/bin/env python3
"""Synthetic regression tests for POMC stitch/HIL receipt binding."""
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import tifffile


SCRIPT = Path(__file__).resolve().parents[1] / "01_annotate_regions.py"
SPEC = importlib.util.spec_from_file_location("paper_annotate_regions", SCRIPT)
ANNOTATE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ANNOTATE
SPEC.loader.exec_module(ANNOTATE)


class RegistrationReceiptTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="fig4_hil_test_", dir="/tmp")
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def _payload(**extra):
        payload = {
            "animal": "Animal1",
            "section": "Sample1_S01",
            "reviewer": "synthetic reviewer",
            "session": "synthetic session",
            "polygons": {
                "ARC": [[[3, 3], [15, 3], [15, 11], [3, 11]]],
                "ME": [],
                "VMN": [],
            },
        }
        payload.update(extra)
        return payload

    def _pomc_tree(self, status="manual_review_required", method=None):
        reconstructed = self.root / "reconstructed"
        section = reconstructed / "Animal1" / "Sample1_S01"
        section.mkdir(parents=True)
        dapi = section / "reconstructed_DAPI.tif"
        image = np.arange(768, dtype=np.uint16).reshape(24, 32)
        tifffile.imwrite(dapi, image, photometric="minisblack", metadata=None)
        registration = section / "stitch_registration.json"
        requires_review = status == "manual_review_required"
        method = method or {
            "manual_review_required": "anatomical_y_fallback",
            "verified_overlap": "direct_ncc",
            "single_half_native": "single_half_native",
        }[status]
        record = {
            "schema_version": ANNOTATE.STITCH_REGISTRATION_SCHEMA,
            "algorithm": ANNOTATE.STITCH_REGISTRATION_ALGORITHM,
            "sample": "Sample1",
            "animal_id": "Animal1",
            "section": "S01",
            "section_index": 1,
            "tile": "",
            "status": status,
            "requires_manual_review": requires_review,
            "method": method,
            "fallback_reason": "synthetic low-confidence overlap" if requires_review else None,
            "transform": {
                "model": "integer_translation_crop_right_medial",
                "right_origin_x_px": 32,
                "right_shift_y_px": 0,
                "right_medial_crop_px": 0,
                "overlap_width_px": 0,
                "gap_px": 0,
                "orientation": "left_then_right",
                "resampling": False,
            },
            "output_shape_px": [24, 32],
            "input": {
                "left_dapi_path": "synthetic_left.tif",
                "right_dapi_path": "synthetic_right.tif",
                "left_dapi_sha256": "1" * 64,
                "right_dapi_sha256": "2" * 64,
            },
            "output": {
                "dapi_path": str(dapi),
                "dapi_sha256": ANNOTATE.sha256(dapi),
            },
            "qc": {
                "raw_overlap_ncc": None if requires_review else 0.94,
                "gradient_overlap_ncc": None if requires_review else 0.88,
                "sift_inlier_count": 0,
                "sift_residual_p95_px": None,
                "overlap_area_px": 0,
            },
        }
        registration.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8")
        return reconstructed, registration, record

    def _app(self, raw, dataset="pomc"):
        analysis = self.root / "analysis"
        output = self.root / "output"
        analysis.mkdir(exist_ok=True)
        return ANNOTATE.App(dataset, raw, analysis, output), output

    def test_manual_registration_requires_explicit_confirmation_and_binds_hash(self):
        raw, registration, registration_record = self._pomc_tree()
        app, output = self._app(raw)
        with self.assertRaisesRegex(ValueError, "explicit image-by-image confirmation"):
            app.save(self._payload())

        saved = app.save(self._payload(
            registered_composite_reviewed=True,
            stitch_registration_confirmation=ANNOTATE.STITCH_MANUAL_CONFIRMATION,
        ))
        self.assertTrue(saved["pair_valid"])
        self.assertTrue(app.rows()[0]["accepted"])
        receipt = json.loads(Path(saved["json_path"]).read_text(encoding="utf-8"))
        self.assertEqual(receipt["schema_version"], 4)
        self.assertTrue(receipt["registered_composite_reviewed"])
        self.assertEqual(
            receipt["stitch_registration_sha256"], ANNOTATE.sha256(registration)
        )

        # Any registration JSON byte change must invalidate the bound HIL receipt.
        registration_record["fallback_reason"] = "synthetic registration revised"
        registration.write_text(
            json.dumps(registration_record, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        self.assertFalse(ANNOTATE.accepted(output, app.sections[0]))
        self.assertFalse(app.rows()[0]["accepted"])

    def test_verified_transform_escalated_by_seam_qc_requires_confirmation(self):
        raw, _, _ = self._pomc_tree(
            status="manual_review_required", method="sift_ransac_raw_ncc")

        app, _ = self._app(raw)
        with self.assertRaisesRegex(ValueError, "explicit image-by-image confirmation"):
            app.save(self._payload())
        saved = app.save(self._payload(
            registered_composite_reviewed=True,
            stitch_registration_confirmation=ANNOTATE.STITCH_MANUAL_CONFIRMATION,
        ))
        receipt = json.loads(Path(saved["json_path"]).read_text(encoding="utf-8"))
        self.assertEqual(receipt["stitch_registration_method"], "sift_ransac_raw_ncc")
        self.assertTrue(receipt["registered_composite_reviewed"])

    def test_verified_overlap_persists_automated_registration_confirmation(self):
        raw, _, _ = self._pomc_tree(status="verified_overlap")
        app, _ = self._app(raw)
        saved = app.save(self._payload())
        receipt = json.loads(Path(saved["json_path"]).read_text(encoding="utf-8"))
        self.assertEqual(receipt["stitch_registration_status"], "verified_overlap")
        self.assertEqual(
            receipt["stitch_registration_confirmation"],
            "automated_overlap_qc_verified",
        )
        self.assertFalse(receipt["registered_composite_reviewed"])
        self.assertTrue(app.rows()[0]["accepted"])

    def test_legacy_algorithm_key_fails_closed(self):
        raw, registration, record = self._pomc_tree(status="verified_overlap")
        record["algorithm_version"] = record.pop("algorithm")
        registration.write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")

        app, _ = self._app(raw)
        self.assertFalse(app.sections[0].registration_valid)
        with self.assertRaisesRegex(ValueError, "no valid stitch_registration"):
            app.save(self._payload())

    def test_gfap_save_remains_schema_three_and_needs_no_registration(self):
        raw = self.root / "gfap_raw"
        section = raw / "Animal1" / "Sample1_S01"
        section.mkdir(parents=True)
        dapi = section / "Image_Sample1_S01_DAPI.tif"
        tifffile.imwrite(
            dapi, np.arange(768, dtype=np.uint16).reshape(24, 32),
            photometric="minisblack", metadata=None,
        )
        app, _ = self._app(raw, dataset="gfap")
        saved = app.save(self._payload())
        receipt = json.loads(Path(saved["json_path"]).read_text(encoding="utf-8"))
        self.assertEqual(receipt["schema_version"], 3)
        self.assertNotIn("stitch_registration_sha256", receipt)
        self.assertTrue(app.rows()[0]["accepted"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
