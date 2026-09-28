#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import io
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
import urllib.parse
import zipfile
from pathlib import Path, PurePosixPath
from unittest import mock


RELEASE = Path(__file__).resolve().parents[1]
PAPER = RELEASE.parent
sys.path.insert(0, str(RELEASE))

import build_archives as builder  # noqa: E402
import reconstruct_paper as reconstructor  # noqa: E402


class FakeResponse(io.BytesIO):
    def __init__(self, payload: bytes, url: str):
        super().__init__(payload)
        self.status = 200
        self.headers = {
            "Content-Length": str(len(payload)),
            "Content-Encoding": "identity",
        }
        self._url = url

    def geturl(self) -> str:
        return self._url


class FakeHTTPSOpener:
    def __init__(self, payloads: dict[str, bytes]):
        self.payloads = payloads
        self.calls: list[tuple[str, float]] = []

    def __call__(self, request, *, timeout: float):
        url = request.full_url
        self.calls.append((url, timeout))
        if url not in self.payloads:
            raise AssertionError(f"Unexpected URL: {url}")
        return FakeResponse(self.payloads[url], url)


class FailOnceHTTPSOpener(FakeHTTPSOpener):
    def __init__(self, payloads: dict[str, bytes]):
        super().__init__(payloads)
        self.failed = False

    def __call__(self, request, *, timeout: float):
        url = request.full_url
        self.calls.append((url, timeout))
        if not self.failed:
            self.failed = True
            raise OSError("temporary transport failure")
        if url not in self.payloads:
            raise AssertionError(f"Unexpected URL: {url}")
        return FakeResponse(self.payloads[url], url)


class ReconstructorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(prefix="paper-reconstruct-test-")
        self.base = Path(self.temporary.name)
        self.paper = self.base / "Paper"
        self.paper.mkdir(mode=0o750)
        os.chmod(self.paper, 0o750)
        self._make_fixture()
        self.archives = builder.build_archives(
            self.paper,
            self.base / "archive-set",
            contract_script=self.paper
            / "scripts"
            / "utilities"
            / "01_validate_bundle.py",
            payload_limit=20_000,
            max_archive_bytes=50_000,
            compression_level=6,
        )
        self.manifest_path = self.archives / "manifest.json"
        self.manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def write(self, relative: str, payload: bytes, mode: int = 0o644) -> Path:
        target = self.paper.joinpath(*PurePosixPath(relative).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        os.chmod(target, mode)
        return target

    def _make_fixture(self) -> None:
        canonical = PAPER / "scripts" / "utilities" / "01_validate_bundle.py"
        self.write(
            "scripts/utilities/01_validate_bundle.py",
            canonical.read_bytes(),
            0o755,
        )
        self.write("README.txt", b"fixture paper\n")
        self.write("run.sh", b"#!/bin/sh\nexit 0\n", 0o755)
        self.write("FigS3/render.py", b"print('histology code')\n")
        self.write("FigS5/channel_tiffs/supplement.bin", b"supplement" * 20)
        self.write("FigS3/exports/tile-a.bin", bytes(range(200)) * 3)
        self.write("FigS3/exports/m26-001/1_L0_rgb.tif", b"WSI" * 200)
        self.write("Fig2/raw_data/acquisition-a.bin", b"A" * 12_000)
        self.write("Fig2/raw_data/acquisition-b.bin", b"B" * 12_000)
        empty = self.paper / "empty" / "kept"
        empty.mkdir(parents=True)
        os.chmod(empty, 0o711)
        self.write("README2.txt", b"excluded root file\n")
        self.write("analyses/private.bin", b"excluded directory\n")

    def shard_rows(self):
        for record in builder.RECORDS:
            for shard in self.manifest["records"][record]["shards"]:
                yield record, shard

    def segmented_url_map(self, segment_bytes: int = 257):
        archives = {}
        payloads = {}
        for archive_index, (_, shard) in enumerate(self.shard_rows(), start=1):
            relative = shard["relative_archive_path"]
            canonical = self.archives.joinpath(
                *PurePosixPath(relative).parts
            ).read_bytes()
            segments = []
            for segment_index, offset in enumerate(
                range(0, len(canonical), segment_bytes), start=1
            ):
                payload = canonical[offset : offset + segment_bytes]
                name = f"{shard['name']}.segment{segment_index:03d}"
                quoted = urllib.parse.quote(name)
                url = (
                    f"https://zenodo.org/records/{2000 + archive_index}/files/"
                    f"{quoted}?download=1"
                )
                segments.append(
                    {
                        "name": name,
                        "size_bytes": len(payload),
                        "sha256": hashlib.sha256(payload).hexdigest(),
                        "url": url,
                    }
                )
                payloads[url] = payload
            archives[relative] = {
                "size_bytes": shard["compressed_bytes"],
                "sha256": shard["sha256"],
                "segments": segments,
            }
        canonical_manifest = builder.canonical_json(self.manifest)
        return (
            {
                "schema": reconstructor.URL_MAP_SCHEMA_V2,
                "manifest": {
                    "size_bytes": len(canonical_manifest),
                    "sha256": hashlib.sha256(canonical_manifest).hexdigest(),
                },
                "archives": archives,
            },
            payloads,
        )

    def assert_exact_tree(self, rebuilt: Path) -> None:
        expected_directories = {
            row["path"]: int(row["mode"], 8)
            for row in self.manifest["tree"]["directories"]
        }
        expected_files = {row["path"]: row for row in self.manifest["tree"]["files"]}
        observed_directories = set()
        observed_files = set()
        for directory, dirnames, filenames in os.walk(rebuilt):
            base = Path(directory)
            for name in dirnames:
                observed_directories.add((base / name).relative_to(rebuilt).as_posix())
            for name in filenames:
                observed_files.add((base / name).relative_to(rebuilt).as_posix())
        self.assertEqual(observed_directories, set(expected_directories))
        self.assertEqual(observed_files, set(expected_files))
        self.assertEqual(
            stat.S_IMODE(rebuilt.stat().st_mode),
            int(self.manifest["tree"]["root_mode"], 8),
        )
        for relative, mode in expected_directories.items():
            target = rebuilt.joinpath(*PurePosixPath(relative).parts)
            self.assertEqual(stat.S_IMODE(target.stat().st_mode), mode)
        for relative, row in expected_files.items():
            target = rebuilt.joinpath(*PurePosixPath(relative).parts)
            self.assertEqual(target.stat().st_size, row["size_bytes"])
            self.assertEqual(builder.sha256_file(target), row["sha256"])
            self.assertEqual(
                stat.S_IMODE(target.stat().st_mode),
                int(row["mode"], 8),
            )
        self.assertFalse((rebuilt / "README2.txt").exists())
        self.assertFalse((rebuilt / "analyses").exists())

    def test_local_cli_reconstructs_exact_tree(self) -> None:
        output = self.base / "rebuilt-local"
        environment = os.environ.copy()
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        completed = subprocess.run(
            [
                sys.executable,
                str(RELEASE / "reconstruct_paper.py"),
                "--manifest",
                str(self.manifest_path),
                "--archive-root",
                str(self.archives),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("[PASS] exact reconstructed Paper tree:", completed.stdout)
        self.assert_exact_tree(output)

    def test_https_url_map_downloads_then_reconstructs(self) -> None:
        urls: dict[str, str] = {}
        payloads: dict[str, bytes] = {}
        for index, (_, shard) in enumerate(self.shard_rows(), start=1):
            relative = shard["relative_archive_path"]
            quoted = urllib.parse.quote(shard["name"])
            url = (
                f"https://zenodo.org/records/{1000 + index}/files/"
                f"{quoted}?download=1"
            )
            urls[relative] = url
            payloads[url] = self.archives.joinpath(
                *PurePosixPath(relative).parts
            ).read_bytes()
        url_map = self.base / "zenodo-urls.json"
        url_map.write_text(
            json.dumps(
                {
                    "schema": reconstructor.URL_MAP_SCHEMA,
                    "archives": urls,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        opener = FakeHTTPSOpener(payloads)
        downloads, output = reconstructor.reconstruct(
            self.manifest_path,
            self.base / "rebuilt-remote",
            url_map_path=url_map,
            download_dir=self.base / "downloads",
            timeout_seconds=9,
            opener=opener,
        )
        self.assertEqual(downloads, self.base / "downloads")
        self.assertEqual({url for url, _ in opener.calls}, set(payloads))
        self.assertTrue(all(timeout == 9 for _, timeout in opener.calls))
        for relative, url in urls.items():
            downloaded = downloads.joinpath(*PurePosixPath(relative).parts)
            self.assertEqual(downloaded.read_bytes(), payloads[url])
        self.assert_exact_tree(output)

    def test_https_v2_segments_reassemble_canonical_archives(self) -> None:
        document, payloads = self.segmented_url_map()
        url_map = self.base / "zenodo-segment-urls.json"
        url_map.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        opener = FakeHTTPSOpener(payloads)
        downloads, output = reconstructor.reconstruct(
            self.manifest_path,
            self.base / "rebuilt-segmented",
            url_map_path=url_map,
            download_dir=self.base / "segment-downloads",
            timeout_seconds=11,
            opener=opener,
        )
        self.assertEqual({url for url, _ in opener.calls}, set(payloads))
        self.assertTrue(all(timeout == 11 for _, timeout in opener.calls))
        expected_paths = {
            shard["relative_archive_path"] for _, shard in self.shard_rows()
        }
        observed_paths = {
            path.relative_to(downloads).as_posix()
            for path in downloads.rglob("*")
            if path.is_file()
        }
        self.assertEqual(observed_paths, expected_paths)
        for _, shard in self.shard_rows():
            relative = shard["relative_archive_path"]
            rebuilt_archive = downloads.joinpath(*PurePosixPath(relative).parts)
            source_archive = self.archives.joinpath(*PurePosixPath(relative).parts)
            self.assertEqual(rebuilt_archive.read_bytes(), source_archive.read_bytes())
        self.assert_exact_tree(output)

    def test_https_v2_retries_one_failed_segment_then_succeeds(self) -> None:
        document, payloads = self.segmented_url_map(segment_bytes=10_000_000)
        url_map = self.base / "zenodo-retry-urls.json"
        url_map.write_text(json.dumps(document), encoding="utf-8")
        opener = FailOnceHTTPSOpener(payloads)
        with mock.patch.object(reconstructor.time, "sleep") as sleeper:
            downloads, output = reconstructor.reconstruct(
                self.manifest_path,
                self.base / "rebuilt-after-retry",
                url_map_path=url_map,
                download_dir=self.base / "retry-downloads",
                timeout_seconds=13,
                opener=opener,
            )
        called_urls = [url for url, _ in opener.calls]
        first_url = called_urls[0]
        self.assertEqual(called_urls.count(first_url), 2)
        sleeper.assert_called_once_with(reconstructor.SEGMENT_RETRY_INITIAL_SECONDS)
        self.assertFalse(
            any(path.suffix == ".segment" for path in downloads.rglob("*"))
        )
        self.assert_exact_tree(output)

    def test_v2_interruption_resumes_archives_segments_and_mid_append(self) -> None:
        document, payloads = self.segmented_url_map()
        archive_items = list(document["archives"].items())
        self.assertGreaterEqual(len(archive_items), 2)
        active_index = next(
            index
            for index, (_, row) in enumerate(archive_items)
            if index > 0 and len(row["segments"]) >= 2
        )
        prior_archives = archive_items[:active_index]
        active_relative, active_archive = archive_items[active_index]
        failed_segment = active_archive["segments"][1]
        failed_url = failed_segment["url"]

        url_map = self.base / "resume-segment-urls.json"
        url_map.write_text(json.dumps(document), encoding="utf-8")
        unavailable = dict(payloads)
        unavailable.pop(failed_url)
        downloads = self.base / "resumable-downloads"
        output = self.base / "resumed-output"
        with mock.patch.object(reconstructor.time, "sleep"), self.assertRaises(
            reconstructor.ReconstructionError
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=url_map,
                download_dir=downloads,
                opener=FakeHTTPSOpener(unavailable),
            )

        stage = downloads.with_name(f".{downloads.name}.downloading")
        self.assertFalse(downloads.exists())
        self.assertTrue((stage / reconstructor.RESUME_RECEIPT_NAME).is_file())
        for relative, _ in prior_archives:
            self.assertTrue(
                stage.joinpath(*PurePosixPath(relative).parts).is_file()
            )
        partial = stage.joinpath(
            *PurePosixPath(active_relative + ".part").parts
        )
        segment_cache = stage.joinpath(
            *PurePosixPath(active_relative + ".segment").parts
        )
        self.assertTrue(partial.is_file())
        with partial.open("ab") as handle:
            handle.write(payloads[failed_url][:17])
        segment_cache.write_bytes(payloads[failed_url])

        opener = FakeHTTPSOpener(payloads)
        observed_downloads, observed_output = reconstructor.reconstruct(
            self.manifest_path,
            output,
            url_map_path=url_map,
            download_dir=downloads,
            opener=opener,
        )
        skipped_urls = {
            segment["url"]
            for _, prior_archive in prior_archives
            for segment in prior_archive["segments"]
        }
        skipped_urls.update(
            segment["url"] for segment in active_archive["segments"][:2]
        )
        self.assertFalse(skipped_urls & {url for url, _ in opener.calls})
        self.assertEqual(observed_downloads, downloads)
        self.assertEqual(observed_output, output)
        self.assertFalse(stage.exists())
        self.assertFalse((downloads / reconstructor.RESUME_RECEIPT_NAME).exists())
        self.assertFalse(
            any(
                path.name.endswith((".part", ".segment"))
                for path in downloads.rglob("*")
            )
        )
        self.assert_exact_tree(output)

    def test_v2_resume_binding_change_and_corrupt_archive_fail_closed(self) -> None:
        document, payloads = self.segmented_url_map()
        archive_items = list(document["archives"].items())
        first_relative, first_archive = archive_items[0]
        _, second_archive = archive_items[1]
        failed_url = second_archive["segments"][0]["url"]
        url_map = self.base / "bound-resume-urls.json"
        url_map.write_text(json.dumps(document), encoding="utf-8")
        unavailable = dict(payloads)
        unavailable.pop(failed_url)
        downloads = self.base / "bound-resume-downloads"
        output = self.base / "bound-resume-output"
        with mock.patch.object(reconstructor.time, "sleep"), self.assertRaises(
            reconstructor.ReconstructionError
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=url_map,
                download_dir=downloads,
                opener=FakeHTTPSOpener(unavailable),
            )

        stage = downloads.with_name(f".{downloads.name}.downloading")
        canonical = stage.joinpath(*PurePosixPath(first_relative).parts)
        with canonical.open("r+b") as handle:
            first = handle.read(1)
            handle.seek(0)
            handle.write(bytes([first[0] ^ 0xFF]))
        opener = FakeHTTPSOpener(payloads)
        with self.assertRaisesRegex(
            reconstructor.ReconstructionError,
            "Retained canonical archive size/SHA-256 differs",
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=url_map,
                download_dir=downloads,
                opener=opener,
            )
        self.assertEqual(opener.calls, [])
        self.assertTrue(stage.exists())
        self.assertFalse(downloads.exists())
        self.assertFalse(output.exists())

        changed = json.loads(json.dumps(document))
        changed_segment = changed["archives"][first_relative]["segments"][-1]
        changed_segment["url"] += "&mirror=changed"
        changed_map = self.base / "changed-bound-resume-urls.json"
        changed_map.write_text(json.dumps(changed), encoding="utf-8")
        with self.assertRaisesRegex(
            reconstructor.ReconstructionError,
            "different manifest or URL map",
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=changed_map,
                download_dir=downloads,
                opener=opener,
            )
        self.assertEqual(opener.calls, [])

    def test_v1_failure_still_removes_transient_download_stage(self) -> None:
        urls = {}
        payloads = {}
        for index, (_, shard) in enumerate(self.shard_rows(), start=1):
            relative = shard["relative_archive_path"]
            url = f"https://zenodo.org/records/{3000 + index}/files/{shard['name']}"
            urls[relative] = url
            payloads[url] = self.archives.joinpath(
                *PurePosixPath(relative).parts
            ).read_bytes()
        payloads.pop(next(iter(urls.values())))
        url_map = self.base / "failed-v1-urls.json"
        url_map.write_text(
            json.dumps(
                {"schema": reconstructor.URL_MAP_SCHEMA_V1, "archives": urls}
            ),
            encoding="utf-8",
        )
        downloads = self.base / "failed-v1-downloads"
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.reconstruct(
                self.manifest_path,
                self.base / "failed-v1-output",
                url_map_path=url_map,
                download_dir=downloads,
                opener=FakeHTTPSOpener(payloads),
            )
        self.assertFalse(downloads.exists())
        self.assertFalse(
            any(
                path.name.startswith(f".{downloads.name}.downloading.")
                for path in self.base.iterdir()
            )
        )

    def test_existing_output_refuses_before_v2_download_state_is_created(self) -> None:
        document, payloads = self.segmented_url_map()
        url_map = self.base / "preflight-v2-urls.json"
        url_map.write_text(json.dumps(document), encoding="utf-8")
        downloads = self.base / "preflight-downloads"
        output = self.base / "existing-output"
        output.mkdir()
        opener = FakeHTTPSOpener(payloads)
        with self.assertRaisesRegex(
            reconstructor.ReconstructionError,
            "reconstructed destination must not already exist",
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=url_map,
                download_dir=downloads,
                opener=opener,
            )
        self.assertEqual(opener.calls, [])
        self.assertFalse(downloads.exists())
        self.assertFalse(
            downloads.with_name(f".{downloads.name}.downloading").exists()
        )

    def test_v2_manifest_binding_and_segment_integrity_fail_closed(self) -> None:
        bad_binding, _ = self.segmented_url_map()
        bad_binding["manifest"]["sha256"] = "0" * 64
        bad_binding_path = self.base / "bad-binding.json"
        bad_binding_path.write_text(json.dumps(bad_binding), encoding="utf-8")
        with self.assertRaisesRegex(
            reconstructor.ReconstructionError, "exact canonical release manifest"
        ):
            reconstructor.load_url_map(bad_binding_path, self.manifest)

        bad_segment, payloads = self.segmented_url_map()
        first_archive = next(iter(bad_segment["archives"].values()))
        first_url = first_archive["segments"][0]["url"]
        payloads[first_url] = payloads[first_url] + b"tampered"
        bad_segment_path = self.base / "bad-segment.json"
        bad_segment_path.write_text(json.dumps(bad_segment), encoding="utf-8")
        downloads = self.base / "failed-segment-downloads"
        output = self.base / "failed-segment-output"
        with mock.patch.object(reconstructor.time, "sleep"), self.assertRaises(
            reconstructor.ReconstructionError
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=bad_segment_path,
                download_dir=downloads,
                opener=FakeHTTPSOpener(payloads),
            )
        self.assertFalse(downloads.exists())
        self.assertFalse(output.exists())

    def test_v2_wrong_segment_order_fails_final_archive_hash(self) -> None:
        document, payloads = self.segmented_url_map()
        first_archive = next(iter(document["archives"].values()))
        self.assertGreater(len(first_archive["segments"]), 1)
        first_archive["segments"][0], first_archive["segments"][1] = (
            first_archive["segments"][1],
            first_archive["segments"][0],
        )
        url_map = self.base / "wrong-segment-order.json"
        url_map.write_text(json.dumps(document), encoding="utf-8")
        downloads = self.base / "failed-order-downloads"
        output = self.base / "failed-order-output"
        with self.assertRaisesRegex(
            reconstructor.ReconstructionError,
            "Concatenated archive size/SHA-256 differs",
        ):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                url_map_path=url_map,
                download_dir=downloads,
                opener=FakeHTTPSOpener(payloads),
            )
        self.assertFalse(downloads.exists())
        self.assertFalse(output.exists())

    def test_tampered_archive_and_bad_url_map_fail_closed(self) -> None:
        forbidden_parent = self.archives / "must-not-be-created"
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.reconstruct(
                self.manifest_path,
                forbidden_parent / "Paper",
                archive_root=self.archives,
            )
        self.assertFalse(forbidden_parent.exists())

        _, shard = next(self.shard_rows())
        archive_path = self.archives.joinpath(
            *PurePosixPath(shard["relative_archive_path"]).parts
        )
        with archive_path.open("r+b") as handle:
            original = handle.read(1)
            handle.seek(0)
            handle.write(bytes([original[0] ^ 0xFF]))
        output = self.base / "must-not-exist"
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                archive_root=self.archives,
            )
        self.assertFalse(output.exists())
        self.assertFalse(
            any(
                path.name.startswith(".must-not-exist.reconstructing.")
                for path in self.base.iterdir()
            )
        )

        expected = {
            row["relative_archive_path"]: (
                "https://zenodo.org/records/123/files/" + row["name"]
            )
            for _, row in self.shard_rows()
        }
        first = next(iter(expected))
        expected[first] = "https://example.com/<REPLACE_ME>"
        bad_map = self.base / "bad-urls.json"
        bad_map.write_text(
            json.dumps(
                {
                    "schema": reconstructor.URL_MAP_SCHEMA,
                    "archives": expected,
                }
            ),
            encoding="utf-8",
        )
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.load_url_map(bad_map, self.manifest)

    def test_extra_traversal_member_and_archive_symlink_are_rejected(self) -> None:
        record, shard = next(self.shard_rows())
        archive_path = self.archives.joinpath(
            *PurePosixPath(shard["relative_archive_path"]).parts
        )
        with zipfile.ZipFile(archive_path, mode="a") as zipped:
            zipped.writestr("Paper/../escape.txt", b"escape")
        old_size = shard["compressed_bytes"]
        shard["compressed_bytes"] = archive_path.stat().st_size
        shard["sha256"] = hashlib.sha256(archive_path.read_bytes()).hexdigest()
        self.manifest["records"][record]["compressed_archive_bytes"] += (
            shard["compressed_bytes"] - old_size
        )
        self.manifest_path.write_bytes(builder.canonical_json(self.manifest))

        output = self.base / "extra-member-output"
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.reconstruct(
                self.manifest_path,
                output,
                archive_root=self.archives,
            )
        self.assertFalse(output.exists())
        self.assertFalse((self.base / "escape.txt").exists())

        backup = self.base / "archive-backup.zip"
        archive_path.rename(backup)
        try:
            archive_path.symlink_to(backup)
        except OSError as exc:
            self.skipTest(f"symlinks unavailable: {exc}")
        symlink_output = self.base / "symlink-output"
        with self.assertRaises(reconstructor.ReconstructionError):
            reconstructor.reconstruct(
                self.manifest_path,
                symlink_output,
                archive_root=self.archives,
            )
        self.assertFalse(symlink_output.exists())


if __name__ == "__main__":
    unittest.main()
