"""Portable preflight regressions; no large data or network required."""
import contextlib
import hashlib
import io
import os
import sys
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import check_reproduction as check
import reconstruct_paper as reconstruct


class PreflightTests(unittest.TestCase):
    def test_different_numerical_thread_profile_blocks_replay(self):
        with patch.dict(os.environ, {"OMP_NUM_THREADS": "16", "MKL_NUM_THREADS": "16", "OPENBLAS_NUM_THREADS": "16"}, clear=True), contextlib.redirect_stdout(io.StringIO()):
            report = check.Report()
            check.check_replay_threads(report, Path("/clone"))
            self.assertEqual(report.failures, 3)

    def test_reference_numerical_thread_profile_passes(self):
        with patch.dict(os.environ, {"OMP_NUM_THREADS": "64", "MKL_NUM_THREADS": "64", "OPENBLAS_NUM_THREADS": "64"}, clear=True), contextlib.redirect_stdout(io.StringIO()):
            report = check.Report()
            check.check_replay_threads(report, Path("/clone"))
            self.assertEqual(report.failures, 0)

    def test_fresh_clone_exact_bytes_and_separate_scratch(self):
        required, recommended = check.storage_budget(240253512366, 319524773963, "download", False, False)
        self.assertEqual(required, 559778286329)
        self.assertEqual(recommended, 700 * check.GB)

    def test_shared_scratch_is_added_not_counted_twice(self):
        required, recommended = check.storage_budget(240253512366, 319524773963, "download", False, True)
        self.assertEqual(required, 559778286329)
        self.assertEqual(recommended, 800 * check.GB)

    def test_existing_archives_credit_only_archive_payload(self):
        required, recommended = check.storage_budget(240253512366, 319524773963, "download", True, False)
        self.assertEqual(required, 319524773963)
        self.assertEqual(recommended, 700 * check.GB - 240253512366)

    def test_reproduce_still_reserves_new_work_space(self):
        required, recommended = check.storage_budget(240253512366, 319524773963, "reproduce", True, True)
        self.assertEqual(required, 0)
        self.assertEqual(recommended, 800 * check.GB - 559778286329)

    def test_filesystem_probe_leaves_no_files(self):
        with tempfile.TemporaryDirectory() as directory:
            check.filesystem_probe(Path(directory))
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_low_disk_fails_before_download(self):
        manifest = {"records": {"a": {"shards": [{"compressed_bytes": 240253512366}]}},
                    "tree": {"payload_bytes": 319524773963, "files": [], "directories": []}}
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            report = check.Report()
            with patch.object(check.shutil, "disk_usage", return_value=SimpleNamespace(free=10 * check.GB)):
                check.check_storage(report, Path(directory), Path('/tmp'), manifest, "download", None)
            self.assertGreater(report.failures, 0)

    def test_full_existing_paper_volume_blocks_even_if_clone_has_space(self):
        manifest = {"records": {"a": {"shards": [{"compressed_bytes": 240253512366}]}},
                    "tree": {"payload_bytes": 319524773963, "files": [], "directories": []}}
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory)
            paper = root / "Paper"
            paper.mkdir()
            report = check.Report()
            def space(path):
                return SimpleNamespace(free=0 if path == paper else 10_000 * check.GB)
            with patch.object(check.shutil, "disk_usage", side_effect=space):
                check.check_storage(report, root, Path('/tmp'), manifest, "reproduce", None)
            self.assertGreater(report.failures, 0)

    def test_explicit_invalid_python_is_not_silently_replaced(self):
        with patch.dict(os.environ, {"PAPER_PYTHON": "/missing/main/python", "PAPER_S3_PYTHON": "/missing/s3/python"}, clear=True):
            found = check.find_environments(None)
        self.assertEqual(found, {"PAPER_PYTHON": None, "PAPER_S3_PYTHON": None})

    def test_missing_paper_fails_without_launching(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(check.subprocess, "run") as launch, contextlib.redirect_stdout(io.StringIO()):
            report = check.Report()
            check.check_paper(report, Path(directory), {"tree": {"files": []}}, {})
            self.assertEqual(report.failures, 1)
            launch.assert_not_called()

    def test_changed_paper_script_blocks_launcher(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(check.subprocess, "run") as launch, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory)
            (root / "Paper").mkdir()
            (root / "Paper/run.sh").write_text("changed bytes")
            manifest = {"tree": {"files": [{"path": "run.sh", "sha256": hashlib.sha256(b"original bytes").hexdigest(), "mode": "0755"}]}}
            report = check.Report()
            check.check_paper(report, root, manifest, {})
            self.assertEqual(report.failures, 1)
            launch.assert_not_called()

    def test_failed_shipped_check_is_not_reported_ready(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(check, "check_hil_polygons"), patch.object(check.subprocess, "run", return_value=SimpleNamespace(returncode=1)) as launch, contextlib.redirect_stdout(io.StringIO()):
            root = Path(directory)
            (root / "Paper").mkdir()
            report = check.Report()
            check.check_paper(report, root, {"tree": {"files": []}}, {"PAPER_PYTHON": "/main/python", "PAPER_S3_PYTHON": "/s3/python"})
            self.assertEqual(report.failures, 1)
            self.assertEqual(launch.call_args.kwargs["cwd"], root / "Paper")
            self.assertEqual(launch.call_args.kwargs["env"]["PAPER_REVIEW_PYTHON"], "/s3/python")

    def test_zenodo_head_size_mismatch_fails(self):
        response = SimpleNamespace(status=200, headers={"Content-Length": "99"}, geturl=lambda: "https://zenodo.org/records/123/files/test")
        with patch.object(check.urllib.request, "urlopen") as opener:
            opener.return_value.__enter__.return_value = response
            status, _ = check.probe_zenodo({"url": response.geturl(), "size_bytes": 100})
        self.assertEqual(status, "FAIL")
        self.assertEqual(opener.call_args.args[0].method, "HEAD")

    def test_existing_interpreter_missing_scientific_module_fails(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            specification = Path(directory) / "environment.yml"
            specification.write_text("name: fixture\n")
            report = check.Report()
            check.check_scientific_environment(report, sys.executable, specification,
                                               ["nonexistent_paper_preflight_test_module"])
            self.assertEqual(report.failures, 1)

    def test_missing_pinned_package_is_visible_as_version_drift(self):
        with tempfile.TemporaryDirectory() as directory, contextlib.redirect_stdout(io.StringIO()):
            specification = Path(directory) / "environment.yml"
            specification.write_text("  - nonexistent-paper-preflight-distribution==1.0\n")
            report = check.Report()
            check.check_scientific_environment(report, sys.executable, specification, ["json"])
            self.assertEqual(report.failures, 0)
            self.assertEqual(report.warnings, 1)

    def test_hil_pixel_mismatch_blocks_ready_status(self):
        result = SimpleNamespace(returncode=0, stdout='{"opencv":"4.12.0","masks":[{"condition":"water","differing_pixels":12},{"condition":"allulose","differing_pixels":31}]}')
        with patch.object(check, "run", return_value=result), contextlib.redirect_stdout(io.StringIO()) as output:
            report = check.Report()
            check.check_hil_polygons(report, "/main/python", Path("/Paper"))
        self.assertEqual(report.failures, 2)
        self.assertIn("requirements-replay-compatibility.txt", output.getvalue())

    def test_exact_hil_masks_pass_runtime_comparison(self):
        result = SimpleNamespace(returncode=0, stdout='{"opencv":"4.10.0","masks":[{"condition":"water","differing_pixels":0},{"condition":"allulose","differing_pixels":0}]}')
        with patch.object(check, "run", return_value=result), contextlib.redirect_stdout(io.StringIO()):
            report = check.Report()
            check.check_hil_polygons(report, "/main/python", Path("/Paper"))
        self.assertEqual(report.failures, 0)

    def test_non_linux_host_cannot_pass(self):
        root = Path(__file__).resolve().parents[1]
        with patch.object(check.platform, "system", return_value="Darwin"), patch.object(check, "filesystem_probe"), patch.object(check, "check_storage"), patch.object(check, "find_conda", return_value=None), patch.object(check, "find_environments", return_value={}), contextlib.redirect_stdout(io.StringIO()) as output:
            status = check.main(["--root", str(root)])
        self.assertEqual(status, 2)
        self.assertIn("requires 64-bit Linux", output.getvalue())


class VerboseDownloadTests(unittest.TestCase):
    def test_verified_download_reports_actual_bytes_and_hash(self):
        payload = b"small verified HTTPS fixture"
        digest = hashlib.sha256(payload).hexdigest()
        url = "https://zenodo.org/records/123/files/chunk"
        class Response(io.BytesIO):
            status = 200
            headers = {"Content-Length": str(len(payload))}
            def geturl(self):
                return url
        with tempfile.TemporaryDirectory() as directory, self.assertLogs(reconstruct.LOG, level="INFO") as output:
            path = Path(directory) / "chunk"
            reconstruct.download_one(url, path, expected_size=len(payload), expected_sha256=digest,
                                     timeout_seconds=1, opener=lambda *a, **kw: Response(payload))
            self.assertEqual(path.read_bytes(), payload)
        self.assertTrue(any("VERIFIED DOWNLOAD" in line and digest in line for line in output.output))

    def test_corrupt_download_never_logs_verified(self):
        url = "https://zenodo.org/records/123/files/chunk"
        class Response(io.BytesIO):
            status = 200
            headers = {"Content-Length": "3"}
            def geturl(self):
                return url
        with tempfile.TemporaryDirectory() as directory, self.assertLogs(reconstruct.LOG, level="INFO") as output:
            with self.assertRaises(reconstruct.ReconstructionError):
                reconstruct.download_one(url, Path(directory) / "chunk", expected_size=3,
                    expected_sha256=hashlib.sha256(b"yes").hexdigest(), timeout_seconds=1,
                    opener=lambda *a, **kw: Response(b"bad"))
        self.assertFalse(any("VERIFIED DOWNLOAD" in line for line in output.output))


if __name__ == "__main__":
    unittest.main()
