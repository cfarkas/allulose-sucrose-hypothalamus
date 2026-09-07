#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import stat
import sys
import tempfile
import unittest
import zipfile
from dataclasses import replace
from pathlib import Path, PurePosixPath


RELEASE = Path(__file__).resolve().parents[1]
PAPER = RELEASE.parent
sys.path.insert(0, str(RELEASE))

import build_archives as archive  # noqa: E402


class ArchiveBuilderTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paper-archive-test-")
        self.base = Path(self.temporary.name)
        self.paper = self.base / "Paper"
        self.paper.mkdir(mode=0o750)
        os.chmod(self.paper, 0o750)
        self._make_fixture()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write(self, relative: str, payload: bytes, mode: int = 0o644) -> Path:
        target = self.paper.joinpath(*PurePosixPath(relative).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        os.chmod(target, mode)
        return target

    def _make_fixture(self) -> None:
        self.write("README.txt", b"public paper\n")
        self.write("README2.txt", b"root-only excluded\n")
        self.write("run.sh", b"#!/bin/sh\nexit 0\n", 0o755)
        self.write("nested/README2.txt", b"nested name remains public\n")
        self.write("FigS3/render.py", b"print('histology code')\n")
        self.write("FigS3/exports/tile-a.bin", bytes(range(200)) * 3)
        self.write("FigS3/exports/m26-001/1_L0_rgb.tif", b"W" * 100)
        self.write("FigS3/registration/obsolete_alignment.bin", b"excluded experiment")
        self.write("FigS5/channel_tiffs/acquisition-s5.bin", b"S" * 50)
        self.write("FigS6/raw_data/predictions.csv", b"cell,prediction\nc1,Rod-like\n")
        self.write("FigS6/01_make_figure_s6.py", b"print(1)\n")
        self.write("Fig2/raw_data/acquisition-a.bin", b"A" * 600)
        self.write("Fig2/raw_data/acquisition-b.bin", b"B" * 600)
        empty = self.paper / "empty" / "kept"
        empty.mkdir(parents=True)
        os.chmod(empty, 0o711)
        self.write("analyses/secret.bin", b"excluded")
        self.write("legacy/secret.bin", b"excluded")
        self.write("literature_webscrap_30_08_2026/metadata.json", b"excluded")
        self.write(
            "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot/reference/atlas.bin",
            b"unused registration snapshot excluded",
        )
        self.write(
            "Fig4/figure_bundle_v3/raw_data/required-acquisition.bin",
            b"nearby data remains public",
        )
        self.write(
            "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot_backup/reference.bin",
            b"same-prefix sibling remains public",
        )
        self.write("nested/__pycache__/secret.pyc", b"excluded")
        self.write("nested/.mypy_cache/secret.json", b"excluded")
        self.write("nested/obsolete.py.orig", b"excluded")
        self.write("nested/conflict.py.rej", b"excluded")
        self.write(
            "literature/source/full_text_open_access.pdf",
            b"third-party article copy excluded",
        )
        self.write("apotome_rebuild_old/secret.bin", b"excluded")

    def build(self, name: str, *, jobs: int = 1) -> Path:
        return archive.build_archives(
            self.paper,
            self.base / name,
            payload_limit=700,
            max_archive_bytes=5000,
            compression_level=6,
            jobs=jobs,
        )

    def test_canonical_exclusions_match(self) -> None:
        canonical = PAPER / "scripts/utilities/01_validate_bundle.py"
        digest = archive.assert_canonical_exclusions(canonical)
        self.assertEqual(digest, archive.sha256_file(canonical))
        spec = importlib.util.spec_from_file_location(
            "canonical_bundle_test", canonical
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        expected = {
            path.relative_to(self.paper).as_posix()
            for path in module.bundle_entries(self.paper)
        }
        observed = {
            path.relative_to(self.paper).as_posix()
            for path in archive.public_walk(self.paper)
        }
        self.assertEqual(observed, expected)

    def test_inventory_exclusions_empty_dirs_modes_hashes_and_records(self) -> None:
        inventory, contract = archive.inventory_tree(self.paper)
        self.assertIsNone(contract)
        files = {item.path: item for item in inventory.files}
        self.assertNotIn("README2.txt", files)
        self.assertNotIn("analyses/secret.bin", files)
        self.assertNotIn("legacy/secret.bin", files)
        self.assertNotIn("literature_webscrap_30_08_2026/metadata.json", files)
        self.assertNotIn(
            "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot/reference/atlas.bin",
            files,
        )
        self.assertNotIn("FigS3/registration/obsolete_alignment.bin", files)
        self.assertIn("Fig4/figure_bundle_v3/raw_data/required-acquisition.bin", files)
        self.assertIn(
            "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot_backup/reference.bin",
            files,
        )
        self.assertNotIn("nested/__pycache__/secret.pyc", files)
        self.assertNotIn("nested/.mypy_cache/secret.json", files)
        self.assertNotIn("nested/obsolete.py.orig", files)
        self.assertNotIn("nested/conflict.py.rej", files)
        self.assertNotIn("literature/source/full_text_open_access.pdf", files)
        self.assertNotIn("apotome_rebuild_old/secret.bin", files)
        self.assertIn("nested/README2.txt", files)
        self.assertEqual(files["run.sh"].mode, 0o755)
        self.assertEqual(
            files["run.sh"].sha256,
            hashlib.sha256(b"#!/bin/sh\nexit 0\n").hexdigest(),
        )
        self.assertEqual(files["FigS3/render.py"].record, "software_core")
        self.assertEqual(
            files["FigS3/exports/tile-a.bin"].record,
            "supplementary_figure_data",
        )
        self.assertEqual(
            files["FigS3/exports/m26-001/1_L0_rgb.tif"].record,
            "figs3_wsi",
        )
        self.assertEqual(
            files["Fig2/raw_data/acquisition-a.bin"].record,
            "main_figure_data",
        )
        self.assertEqual(
            files["FigS5/channel_tiffs/acquisition-s5.bin"].record,
            "supplementary_figure_data",
        )
        self.assertEqual(
            archive.assign_record("TESIS_FINAL_NS.docx", 300_000_000),
            "software_core",
        )
        self.assertEqual(files["FigS6/raw_data/predictions.csv"].record, "supplementary_figure_data")
        self.assertEqual(files["FigS6/01_make_figure_s6.py"].record, "software_core")
        for sample in range(1, 25):
            self.assertEqual(
                archive.assign_record(
                    f"FigS3/exports/m26-{sample:03d}/1_L0_rgb.tif", 12
                ),
                "figs3_wsi",
            )
        for near_match in (
            "FigS3/exports/m26-000/1_L0_rgb.tif",
            "FigS3/exports/m26-025/1_L0_rgb.tif",
            "FigS3/exports/m26-001/1_L0_rgb.tiff",
            "FigS3/exports/m26-001/1_L2_rgb.tif",
            "FigS3/exports/m26-001/nested/1_L0_rgb.tif",
        ):
            self.assertEqual(
                archive.assign_record(near_match, 12), "supplementary_figure_data"
            )
        for component in (
            "exports",
            "histoplus",
            "organ_segmentation",
            "raw",
            "raw_data",
        ):
            self.assertEqual(
                archive.assign_record(f"FigS3/{component}/artifact.bin", 12),
                "supplementary_figure_data",
            )
        with self.assertRaisesRegex(archive.ArchiveError, "Unknown figure data root"):
            archive.assign_record("Fig6/raw_data/future.bin", 12)
        with self.assertRaises(archive.ArchiveError):
            archive.assign_record("invalid.bin", -1)
        directories = {item.path: item for item in inventory.directories}
        self.assertTrue(directories["empty/kept"].empty)
        self.assertEqual(directories["empty/kept"].mode, 0o711)

    def test_build_is_deterministic_and_manifest_is_reconstructable(self) -> None:
        first = self.build("release-one")
        second = self.build("release-two", jobs=3)
        first_files = sorted(
            path.relative_to(first) for path in first.rglob("*") if path.is_file()
        )
        second_files = sorted(
            path.relative_to(second) for path in second.rglob("*") if path.is_file()
        )
        self.assertEqual(first_files, second_files)
        for relative in first_files:
            self.assertEqual(
                (first / relative).read_bytes(), (second / relative).read_bytes()
            )

        manifest = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
        archive.validate_manifest(manifest)
        excluded = json.loads(json.dumps(manifest))
        row = excluded["tree"]["files"][0]
        row["path"] = "literature_webscrap_30_08_2026/private.json"
        row["member"] = "Paper/" + row["path"]
        excluded["tree"]["sha256"] = archive.tree_digest(
            excluded["tree"]["root_mode"],
            excluded["tree"]["directories"],
            excluded["tree"]["files"],
        )
        with self.assertRaisesRegex(archive.ArchiveError, "excluded file"):
            archive.validate_manifest(excluded)
        tampered = json.loads(json.dumps(manifest))
        tampered["records"]["main_figure_data"]["shards"][0]["member_count"] += 1
        with self.assertRaises(archive.ArchiveError):
            archive.validate_manifest(tampered)
        wrong_record = json.loads(json.dumps(manifest))
        wsi_row = next(
            row
            for row in wrong_record["tree"]["files"]
            if row["path"] == "FigS3/exports/m26-001/1_L0_rgb.tif"
        )
        wsi_row["record"] = "supplementary_figure_data"
        with self.assertRaisesRegex(archive.ArchiveError, "classification differs"):
            archive.validate_manifest(wrong_record)
        self.assertEqual(set(manifest["records"]), set(archive.RECORDS))
        self.assertGreater(manifest["records"]["main_figure_data"]["shard_count"], 1)
        for record in archive.RECORDS:
            for shard in manifest["records"][record]["shards"]:
                self.assertLessEqual(shard["uncompressed_bytes"], 700)
                self.assertLessEqual(shard["compressed_bytes"], 5000)
                path = first / shard["relative_archive_path"]
                self.assertEqual(shard["sha256"], archive.sha256_file(path))
                with zipfile.ZipFile(path) as zipped:
                    for info in zipped.infolist():
                        self.assertEqual(info.compress_type, zipfile.ZIP_DEFLATED)
                        self.assertEqual(info.date_time, archive.ZIP_DATETIME)

        rebuilt = self.base / "rebuilt" / "Paper"
        rebuilt.mkdir(parents=True)
        for row in manifest["tree"]["directories"]:
            target = rebuilt.joinpath(*PurePosixPath(row["path"]).parts)
            target.mkdir(parents=True, exist_ok=True)
        rows_by_shard = {}
        for row in manifest["tree"]["files"]:
            rows_by_shard.setdefault((row["record"], row["shard"]), []).append(row)
        for (record, shard_name), rows in rows_by_shard.items():
            with zipfile.ZipFile(first / record / shard_name) as zipped:
                for row in rows:
                    target = rebuilt.joinpath(*PurePosixPath(row["path"]).parts)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(zipped.read(row["member"]))
                    os.chmod(target, int(row["mode"], 8))
        for row in sorted(
            manifest["tree"]["directories"], key=lambda item: item["path"], reverse=True
        ):
            os.chmod(
                rebuilt.joinpath(*PurePosixPath(row["path"]).parts), int(row["mode"], 8)
            )
        os.chmod(rebuilt, int(manifest["tree"]["root_mode"], 8))
        for row in manifest["tree"]["files"]:
            target = rebuilt.joinpath(*PurePosixPath(row["path"]).parts)
            self.assertEqual(target.stat().st_size, row["size_bytes"])
            self.assertEqual(archive.sha256_file(target), row["sha256"])
            self.assertEqual(stat.S_IMODE(target.stat().st_mode), int(row["mode"], 8))
        self.assertTrue((rebuilt / "empty/kept").is_dir())
        self.assertEqual(stat.S_IMODE((rebuilt / "empty/kept").stat().st_mode), 0o711)

    def test_fail_closed_path_symlink_duplicate_and_oversize(self) -> None:
        with self.assertRaisesRegex(archive.ArchiveError, "positive integer"):
            archive.build_archives(self.paper, self.base / "invalid-jobs", jobs=0)
        for unsafe in ("../escape", "/absolute", "a/../../b", "a\\b"):
            with self.assertRaises(archive.ArchiveError):
                archive.safe_relative(unsafe)
        inventory, _ = archive.inventory_tree(self.paper)
        duplicate = list(inventory.files) + [
            replace(inventory.files[0], record="software_core")
        ]
        with self.assertRaises(archive.ArchiveError):
            archive.validate_inventory(duplicate)
        too_large = replace(inventory.files[0], size=701)
        with self.assertRaises(archive.ArchiveError):
            archive.plan_shards([too_large], payload_limit=700)
        link = self.paper / "public-link"
        root_link = self.base / "Paper-link"
        try:
            link.symlink_to(self.paper / "README.txt")
            root_link.symlink_to(self.paper, target_is_directory=True)
        except OSError as exc:
            self.skipTest(f"symlinks unavailable: {exc}")
        with self.assertRaises(archive.ArchiveError):
            archive.inventory_tree(root_link)
        with self.assertRaises(archive.ArchiveError):
            archive.inventory_tree(self.paper)


if __name__ == "__main__":
    unittest.main()
