#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path, PurePosixPath


RELEASE = Path(__file__).resolve().parents[1]
PAPER = RELEASE.parent
sys.path.insert(0, str(RELEASE))

import build_archives as archives  # noqa: E402
import finalize_archives as finalizer  # noqa: E402


class ArchiveFinalizerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paper-finalizer-test-")
        self.base = Path(self.temporary.name)
        self.paper = self.base / "Paper"
        self.paper.mkdir(mode=0o750)
        os.chmod(self.paper, 0o750)
        canonical = PAPER / "scripts" / "utilities" / "01_validate_bundle.py"
        self.write(
            "scripts/utilities/01_validate_bundle.py", canonical.read_bytes(), 0o755
        )
        self.write("README.txt", b"prior software-core payload\n")
        self.write("Fig2/raw_data/main.bin", b"main data\n")
        self.write("FigS3/exports/supplementary.bin", b"supplementary data\n")
        self.write("FigS3/exports/m26-001/1_L0_rgb.tif", b"WSI data\n")

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write(self, relative: str, payload: bytes, mode: int = 0o644) -> Path:
        target = self.paper.joinpath(*PurePosixPath(relative).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        os.chmod(target, mode)
        return target

    def assert_archive_hashes(self, root: Path, manifest: dict) -> None:
        for record in archives.RECORDS:
            for shard in manifest["records"][record]["shards"]:
                path = root / shard["relative_archive_path"]
                self.assertEqual(archives.sha256_file(path), shard["sha256"])

    def test_hard_linked_prior_manifest_is_not_modified(self) -> None:
        prior = archives.build_archives(
            self.paper, self.base / "archive-generation-1", jobs=2
        )
        prior_manifest_path = prior / "manifest.json"
        prior_bytes = prior_manifest_path.read_bytes()
        prior_stat = prior_manifest_path.stat()
        prior_document = json.loads(prior_bytes)

        stage = self.base / ".archive-generation-2.building"
        shutil.copytree(prior, stage, copy_function=os.link)
        self.assertTrue(os.path.samefile(prior_manifest_path, stage / "manifest.json"))

        self.write("README.txt", b"changed software-core payload\n")
        output = self.base / "archive-generation-2"
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            observed = finalizer.finalize_archives(self.paper, stage, output, jobs=2)

        self.assertEqual(observed, output)
        self.assertFalse(stage.exists())
        self.assertEqual(prior_manifest_path.read_bytes(), prior_bytes)
        self.assertEqual(prior_manifest_path.stat().st_ino, prior_stat.st_ino)
        self.assertEqual(prior_manifest_path.stat().st_mtime_ns, prior_stat.st_mtime_ns)
        self.assertFalse(
            os.path.samefile(prior_manifest_path, output / "manifest.json")
        )

        current_bytes = (output / "manifest.json").read_bytes()
        self.assertNotEqual(current_bytes, prior_bytes)
        current_document = json.loads(current_bytes)
        archives.validate_manifest(prior_document)
        archives.validate_manifest(current_document)
        self.assertNotEqual(
            prior_document["tree"]["sha256"], current_document["tree"]["sha256"]
        )
        self.assert_archive_hashes(prior, prior_document)
        self.assert_archive_hashes(output, current_document)
        self.assertIn("[SHARD] rebuilt: software_core/", captured.getvalue())
        self.assertIn("[PASS] Final archive set:", captured.getvalue())

    def test_unexpected_staged_file_fails_without_renaming(self) -> None:
        prior = archives.build_archives(self.paper, self.base / "prior", jobs=1)
        stage = self.base / ".next.building"
        shutil.copytree(prior, stage, copy_function=os.link)
        (stage / "unreviewed.txt").write_text("unexpected\n", encoding="utf-8")
        output = self.base / "next"

        with self.assertRaisesRegex(finalizer.FinalizationError, "unexpected"):
            finalizer.finalize_archives(self.paper, stage, output, jobs=1)

        self.assertTrue(stage.is_dir())
        self.assertFalse(output.exists())

    def test_explicit_hash_can_authorize_legacy_rows_only(self) -> None:
        prior = archives.build_archives(
            self.paper, self.base / "archive-generation-legacy", jobs=2
        )
        legacy_path = prior / "manifest.json"
        legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
        legacy["tree"]["directories"].append(
            {"path": "FigS3/registration", "mode": "0755", "empty": True}
        )
        legacy["tree"]["directories"].sort(key=lambda row: row["path"])
        legacy["tree"]["directory_count"] = len(legacy["tree"]["directories"])
        legacy["tree"]["sha256"] = archives.tree_digest(
            legacy["tree"]["root_mode"],
            legacy["tree"]["directories"],
            legacy["tree"]["files"],
        )
        legacy_path.write_bytes(archives.canonical_json(legacy))
        legacy_sha256 = archives.sha256_file(legacy_path)
        with self.assertRaisesRegex(archives.ArchiveError, "excluded directory"):
            archives.validate_manifest(legacy)

        stage = self.base / ".archive-generation-current.building"
        shutil.copytree(prior, stage, copy_function=os.link)
        output = self.base / "archive-generation-current"
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            finalizer.finalize_archives(
                self.paper,
                stage,
                output,
                jobs=2,
                trusted_legacy_manifest_sha256=legacy_sha256,
            )

        current = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        archives.validate_manifest(current)
        for record in archives.RECORDS:
            prior_zip = next((prior / record).glob("*.zip"))
            current_zip = next((output / record).glob("*.zip"))
            self.assertTrue(os.path.samefile(prior_zip, current_zip))
        self.assertIn("SHA-256-bound shard rows are eligible", captured.getvalue())

        wrong_stage = self.base / ".wrong-hash.building"
        shutil.copytree(prior, wrong_stage, copy_function=os.link)
        with self.assertRaisesRegex(
            finalizer.FinalizationError, "legacy manifest SHA-256 differs"
        ):
            finalizer.finalize_archives(
                self.paper,
                wrong_stage,
                self.base / "wrong-hash-output",
                jobs=1,
                trusted_legacy_manifest_sha256="0" * 64,
            )


if __name__ == "__main__":
    unittest.main()
