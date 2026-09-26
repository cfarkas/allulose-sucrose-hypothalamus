#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
import unittest
import zipfile
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "release"))
from build_archives import local_input_path


def input_path(relative):
    return local_input_path(ROOT, relative)

OFFICE_RECEIPT_RELATIVE = Path(
    "manuscript_sources/PUBLIC_RELEASE_OFFICE_SANITIZATION_20260902.json"
)
REVISION_RECEIPT_RELATIVE = Path(
    "manuscript_sources/PUBLIC_RELEASE_MANUSCRIPT_FIGS3_REVISION_20260904.json"
)
CONSUMPTION_RECEIPT_RELATIVE = Path(
    "manuscript_sources/PUBLIC_RELEASE_SINGLE_BOTTLE_CONSUMPTION_REVISION_20260904.json"
)
OFFICE_RECEIPT = input_path(OFFICE_RECEIPT_RELATIVE)
REVISION_RECEIPT = input_path(REVISION_RECEIPT_RELATIVE)
CONSUMPTION_RECEIPT = input_path(CONSUMPTION_RECEIPT_RELATIVE)

sys.path.insert(0, str(input_path("manuscript_sources")))
from caps_en import CAPS_EN  # noqa: E402


EXPECTED_SANITIZED_STILL_CURRENT = {
    "manuscript_sources/TESIS_FINAL_NS4_source_20260830.docx": (
        17222762,
        "676316082c324982ee53f7bc4026467a9028629cdd7b68368d095fcdb76d342b",
    ),
    "FigS4/raw_data/UCSC_MICE.xlsx": (
        29205,
        "9759e888643b82c3f4edac5d0a2eb3d76d1dcfcabe700a1851ed47002ab1d657",
    ),
}

EXPECTED_REVISION_CURRENT = {
    "TESIS_FINAL_NS(4)_revision_cientifica_metodos_figuras_final.docx": (
        300232507,
        "2200f59069a67030b87767dd177be56572ad4099726cb77ad0504ee565417f74",
    ),
    "TESIS_FINAL_NS(4)_revision_cientifica_metodos_figuras_final.pdf": (
        18674779,
        "ee28827f6387dac101b6534a6e74752b480291db0c6770ef2ebad43c7f328755",
    ),
}

EXPECTED_REMOVED_PARTS = {
    "word/comments.xml",
    "word/commentsExtended.xml",
    "word/commentsExtensible.xml",
    "word/commentsIds.xml",
    "word/people.xml",
}

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


class PublicOfficeSanitizationReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.office = json.loads(OFFICE_RECEIPT.read_text(encoding="utf-8"))
        cls.revision = json.loads(REVISION_RECEIPT.read_text(encoding="utf-8"))
        cls.consumption = json.loads(
            CONSUMPTION_RECEIPT.read_text(encoding="utf-8")
        )

    def test_receipts_bind_every_current_public_file(self) -> None:
        office_entries = {
            item["path"]: item
            for group in ("docx_files", "xlsx_files", "unchanged_companion_files")
            for item in self.office[group]
        }
        for relative, (size, expected_hash) in EXPECTED_SANITIZED_STILL_CURRENT.items():
            path = input_path(relative)
            self.assertEqual(path.stat().st_size, size)
            self.assertEqual(sha256(path), expected_hash)
            identity = office_entries[relative].get("after", office_entries[relative])
            self.assertEqual(identity["size_bytes"], size)
            self.assertEqual(identity["sha256"], expected_hash)

        revision_entries = {
            self.consumption["docx"]["path"]: self.consumption["docx"]["after"],
            self.consumption["pdf"]["path"]: self.consumption["pdf"]["after"],
        }
        self.assertEqual(set(revision_entries), set(EXPECTED_REVISION_CURRENT))
        for relative, (size, expected_hash) in EXPECTED_REVISION_CURRENT.items():
            path = input_path(relative)
            self.assertEqual(path.stat().st_size, size)
            self.assertEqual(sha256(path), expected_hash)
            self.assertEqual(revision_entries[relative]["size_bytes"], size)
            self.assertEqual(revision_entries[relative]["sha256"], expected_hash)

        predecessor_docx = office_entries[self.revision["docx"]["path"]]["after"]
        predecessor_pdf = office_entries[self.revision["pdf"]["path"]]
        self.assertEqual(self.revision["docx"]["before"], predecessor_docx)
        for key in ("size_bytes", "sha256"):
            self.assertEqual(self.revision["pdf"]["before"][key], predecessor_pdf[key])

        self.assertEqual(
            self.consumption["docx"]["before"], self.revision["docx"]["after"]
        )
        self.assertEqual(
            self.consumption["pdf"]["before"], self.revision["pdf"]["after"]
        )

    def test_transition_and_validation_contracts_are_complete(self) -> None:
        self.assertTrue(
            self.office["authoritative_for_original_sanitization_transition"]
        )
        self.assertFalse(self.office["authoritative_for_public_live_office_identity"])
        self.assertFalse(
            self.revision["authoritative_for_public_live_manuscript_identity"]
        )
        self.assertTrue(
            self.consumption["authoritative_for_public_live_manuscript_identity"]
        )
        self.assertEqual(
            self.office["public_live_manuscript_identity_successor"],
            REVISION_RECEIPT_RELATIVE.as_posix(),
        )
        self.assertEqual(
            self.revision["public_live_manuscript_identity_successor"],
            CONSUMPTION_RECEIPT_RELATIVE.as_posix(),
        )
        documents = {item["role"]: item for item in self.office["docx_files"]}
        self.assertEqual(
            documents["live_thesis_docx"]["transformation_counts"][
                "accepted_insertion_wrappers"
            ],
            0,
        )
        source_counts = documents["manuscript_build_source_docx"][
            "transformation_counts"
        ]
        self.assertEqual(source_counts["accepted_insertion_wrappers"], 25)
        self.assertEqual(source_counts["deleted_revision_elements"], 15)
        self.assertEqual(source_counts["property_change_elements"], 129)
        for document in documents.values():
            self.assertEqual(set(document["removed_parts"]), EXPECTED_REMOVED_PARTS)
            for key in (
                "accepted_revisions_visible_text_unchanged",
                "member_inventory_equal_except_removed_parts",
                "media_members_byte_identical",
                "zip_crc_pass",
            ):
                self.assertIs(document["validation"][key], True)
        revision_validation = self.revision["docx"]["validation"]
        for key in (
            "zip_crc_pass",
            "member_inventory_unchanged",
            "all_other_members_byte_identical",
            "embedded_figure_byte_identical_to_canonical_png",
            "figure_dimensions_unchanged",
            "comment_and_people_parts_absent",
            "comment_markers_absent",
            "tracked_revision_markup_absent",
            "revision_session_metadata_absent",
            "last_modified_by_absent",
        ):
            self.assertIs(revision_validation[key], True)

    def test_live_docx_contains_segmentation_only_figure_s3(self) -> None:
        document_path = input_path(self.revision["docx"]["path"])
        figure = self.revision["canonical_figure"]
        with zipfile.ZipFile(document_path) as archive:
            self.assertIsNone(archive.testzip())
            self.assertFalse(EXPECTED_REMOVED_PARTS.intersection(archive.namelist()))
            root = ElementTree.fromstring(archive.read("word/document.xml"))
            paragraphs = [
                "".join(node.text or "" for node in paragraph.iter(W + "t"))
                for paragraph in root.iter(W + "p")
            ]
            starts = [
                index
                for index, text in enumerate(paragraphs)
                if text.startswith("Figure S3.")
            ]
            self.assertEqual(len(starts), 1)
            caption = paragraphs[starts[0] : starts[0] + 4]
            self.assertEqual(caption, [CAPS_EN[key] for key in (671, 672, 673, 674)])
            joined = "\n".join(caption).casefold()
            self.assertIsNone(re.search(r"\bregistered\s+multi-organ\b", joined))
            self.assertNotIn("registered independently", joined)
            self.assertIn("m26-015", joined)
            self.assertIn("closest to the median among water animals", joined)
            embedded = archive.read(self.revision["docx"]["figure_member"])
            self.assertEqual(sha256_bytes(embedded), figure["sha256"])
        self.assertEqual(sha256(ROOT / figure["path"]), figure["sha256"])

    @unittest.skipUnless(shutil.which("pdftotext"), "pdftotext is unavailable")
    def test_live_pdf_contains_the_same_figure_s3_caption(self) -> None:
        pdf = input_path(self.revision["pdf"]["path"])
        text = subprocess.check_output(
            ["pdftotext", "-layout", str(pdf), "-"], text=True
        )
        start = text.index("Figure S3. Multi-organ H&E computer-vision analysis.")
        end = text.index("Artículo — Figure S4", start)
        caption = " ".join(text[start:end].split()).casefold()
        self.assertIn("m26-015", caption)
        self.assertIn("closest to the median among water animals", caption)
        self.assertIsNone(re.search(r"\bregistered\s+multi-organ\b", caption))
        self.assertNotIn("registered independently", caption)

    def test_historical_receipts_still_contain_the_bridged_hashes(self) -> None:
        hil = json.loads(
            (
                input_path("manuscript_sources/TESIS_FINAL_NS4_hil_refresh_receipt_20260901.json")
            ).read_text(encoding="utf-8")
        )
        completion = json.loads(
            (
                input_path("manuscript_sources/TESIS_FINAL_NS4_completion_receipt_20260830.json")
            ).read_text(encoding="utf-8")
        )
        bridges = {
            item["receipt_field"]: item
            for item in self.office["historical_receipt_bridges"]
        }
        self.assertEqual(
            hil["live_promotion"]["final_live_docx_sha256"],
            bridges["live_promotion.final_live_docx_sha256"]["historical_sha256"],
        )
        self.assertEqual(
            completion["source_sha256"], bridges["source_sha256"]["historical_sha256"]
        )

    def test_public_receipts_have_no_private_payload_or_machine_local_path(
        self,
    ) -> None:
        expected_privacy = {
            "comment_text_included": False,
            "reviewer_or_user_identities_included": False,
            "account_identifiers_included": False,
            "machine_local_source_paths_included": False,
            "digests_of_removed_sensitive_values_included": False,
        }
        for path, receipt in (
            (OFFICE_RECEIPT, self.office),
            (REVISION_RECEIPT, self.revision),
            (CONSUMPTION_RECEIPT, self.consumption),
        ):
            serialized = path.read_text(encoding="utf-8")
            self.assertNotIn("sensitive_value_sha256", serialized)
            self.assertIsNone(
                re.search(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+", serialized)
            )
            self.assertIsNone(re.search(r'"(?:/home/|/media/|/tmp/)', serialized))
            self.assertEqual(receipt["privacy"], expected_privacy)

    def test_readmes_point_to_both_receipts(self) -> None:
        for receipt_name in (
            OFFICE_RECEIPT_RELATIVE.name,
            REVISION_RECEIPT_RELATIVE.name,
            CONSUMPTION_RECEIPT_RELATIVE.name,
        ):
            self.assertIn(
                receipt_name, (ROOT / "README.txt").read_text(encoding="utf-8")
            )
            self.assertIn(
                receipt_name,
                input_path("manuscript_sources/README.txt").read_text(encoding="utf-8"),
            )


if __name__ == "__main__":
    unittest.main()
