"""Exercise orchestration and failure boundaries without downloading the paper."""
import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import reproduce


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.runner = reproduce.Runner(self.root, environ={})
        self.conda = patch.object(reproduce, "find_conda", return_value="/test/conda")
        self.conda.start()
        self.addCleanup(self.conda.stop)
        self.stdout = contextlib.redirect_stdout(io.StringIO())
        self.stdout.__enter__()
        self.addCleanup(self.stdout.__exit__, None, None, None)

    def paper_fixture(self):
        rows = []
        for relative in (*reproduce.SETUP_FILES, "reproduce_all_figures.sh"):
            path = self.root / "Paper" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# small immutable test fixture\n")
            path.chmod(0o755 if relative.endswith(".sh") else 0o644)
            rows.append({"path": relative, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                         "mode": "0755" if relative.endswith(".sh") else "0644"})
        return {"tree": {"files": rows}}

    def test_low_space_stops_before_download_or_setup(self):
        with patch.object(self.runner, "command", side_effect=reproduce.ReproductionError("Insufficient space")) as command, patch.object(self.runner, "reconstruct") as data, patch.object(self.runner, "setup_software") as software:
            with self.assertRaisesRegex(reproduce.ReproductionError, "Insufficient space"):
                self.runner.execute()
            self.assertIn("--system-only", command.call_args.args[1])
            self.assertIn("--strict", command.call_args.args[1])
            data.assert_not_called()
            software.assert_not_called()

    def test_missing_conda_stops_before_large_download(self):
        with patch.object(reproduce, "find_conda", return_value=None), patch.object(self.runner, "command") as command:
            with self.assertRaisesRegex(reproduce.ReproductionError, "SETUP.md"):
                self.runner.execute()
            command.assert_not_called()

    def test_computer_check_does_not_download_or_install(self):
        with patch.object(self.runner, "command", return_value=""), patch.object(self.runner, "reconstruct") as data, patch.object(self.runner, "setup_software") as software:
            self.runner.execute(check_only=True)
            data.assert_not_called()
            software.assert_not_called()

    def test_existing_data_computer_check_does_not_require_installed_science(self):
        (self.root / "Paper").mkdir()
        with patch.object(self.runner, "command", return_value="") as command:
            self.runner.check_computer()
            arguments = command.call_args.args[1]
            self.assertEqual(arguments[arguments.index("--stage") + 1], "reproduce")
            self.assertIn("--system-only", arguments)

    def test_new_clone_uses_public_zenodo_map(self):
        def download(*args, **kwargs):
            self.manifest = self.paper_fixture()
            return ""
        with patch.object(self.runner, "command", side_effect=download) as command, patch.object(reproduce, "load_manifest", side_effect=lambda _: self.manifest):
            self.runner.reconstruct()
            arguments = command.call_args.args[1]
            self.assertIn("--url-map", arguments)
            self.assertEqual(arguments[arguments.index("--download-dir") + 1], self.root / "zenodo-archives")
            self.assertIn("--verbose", arguments)

    def test_completed_archives_are_reused_without_network_download(self):
        archives = self.root / "zenodo-archives"
        archives.mkdir()
        with patch.object(self.runner, "command", return_value=""):
            self.runner.check_computer()
        def extract(*args, **kwargs):
            self.manifest = self.paper_fixture()
            return ""
        with patch.object(self.runner, "command", side_effect=extract) as command, patch.object(reproduce, "load_manifest", side_effect=lambda _: self.manifest):
            self.runner.reconstruct()
            arguments = command.call_args.args[1]
            self.assertNotIn("--url-map", arguments)
            self.assertEqual(arguments[arguments.index("--archive-root") + 1], archives)

    def test_existing_paper_is_authenticated_without_reconstruction(self):
        manifest = self.paper_fixture()
        with patch.object(reproduce, "load_manifest", return_value=manifest), patch.object(self.runner, "command") as command:
            self.runner.reconstruct()
            command.assert_not_called()

    def test_modified_installer_specification_blocks_setup(self):
        manifest = self.paper_fixture()
        (self.root / "Paper" / reproduce.SETUP_FILES[0]).write_text("untrusted changes")
        with patch.object(reproduce, "load_manifest", return_value=manifest), patch.object(self.runner, "check_computer"), patch.object(self.runner, "setup_software") as software:
            with self.assertRaisesRegex(reproduce.ReproductionError, "specification is missing or changed"):
                self.runner.execute()
            software.assert_not_called()

    def test_modified_paper_script_blocks_setup(self):
        manifest = self.paper_fixture()
        (self.root / "Paper/reproduce_all_figures.sh").write_text("untrusted changes")
        with patch.object(reproduce, "load_manifest", return_value=manifest), patch.object(self.runner, "check_computer"), patch.object(self.runner, "setup_software") as software:
            with self.assertRaisesRegex(reproduce.ReproductionError, "programs are missing or changed"):
                self.runner.execute()
            software.assert_not_called()

    def test_explicit_interpreters_are_used_without_installing_into_them(self):
        self.runner.env.update(PAPER_PYTHON=sys.executable, PAPER_S3_PYTHON=sys.executable)
        with patch.object(reproduce, "find_conda", return_value=None), patch.object(self.runner, "command", return_value=""), patch.object(self.runner, "install_environment") as install:
            self.runner.check_computer()
            self.runner.setup_software()
            install.assert_not_called()
        self.assertEqual(self.runner.env["PAPER_PYTHON_CZI"], sys.executable)
        self.assertEqual(self.runner.env["PAPER_REVIEW_PYTHON"], sys.executable)
        for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
            self.assertEqual(self.runner.env[name], "64")

    def test_invalid_explicit_interpreter_is_not_silently_ignored(self):
        self.runner.env["PAPER_PYTHON"] = "/missing/python"
        with self.assertRaisesRegex(reproduce.ReproductionError, "PAPER_PYTHON points to a missing program"):
            self.runner.check_computer()

    def test_software_installation_stays_inside_clone_and_applies_correction(self):
        self.paper_fixture()
        (self.root / "requirements-replay-compatibility.txt").write_text("opencv-python-headless==4.10.0.84\n")
        self.runner.conda = "/test/conda"
        with patch.object(self.runner, "command", return_value="") as command:
            self.runner.setup_software()
        calls = [call.args[1] for call in command.call_args_list]
        self.assertEqual(len(calls), 4)
        self.assertEqual(calls[0][calls[0].index("--prefix") + 1], self.root / ".paper-envs/main")
        self.assertEqual(calls[1][-1], self.root / "requirements-replay-compatibility.txt")
        self.assertEqual(calls[2][calls[2].index("--prefix") + 1], self.root / ".paper-envs/s3")
        self.assertEqual(calls[3][-1], self.root / "Paper" / reproduce.SETUP_FILES[2])
        self.assertEqual(self.runner.env["CONDA_PKGS_DIRS"], str(self.root / ".paper-envs/cache/conda"))

    def test_interrupted_main_installation_is_updated_and_retried(self):
        self.paper_fixture()
        correction = self.root / "requirements-replay-compatibility.txt"
        correction.write_text("opencv-python-headless==4.10.0.84\n")
        history = self.root / ".paper-envs/main/conda-meta/history"
        history.parent.mkdir(parents=True)
        history.touch()
        self.runner.conda = "/test/conda"
        with patch.object(self.runner, "command", side_effect=reproduce.ReproductionError("connection lost")) as command:
            with self.assertRaisesRegex(reproduce.ReproductionError, "connection lost"):
                self.runner.install_environment("main", [self.root / "Paper" / reproduce.SETUP_FILES[0], correction])
        self.assertEqual(command.call_args.args[1][1:3], ["env", "update"])
        self.assertFalse((self.root / ".paper-envs/main-prepared.json").exists())

    def test_completed_installation_is_reused_but_changed_specification_is_retried(self):
        self.paper_fixture()
        correction = self.root / "requirements-replay-compatibility.txt"
        correction.write_text("opencv-python-headless==4.10.0.84\n")
        specifications = [self.root / "Paper" / reproduce.SETUP_FILES[0], correction]
        python = self.root / ".paper-envs/main/bin/python"
        python.parent.mkdir(parents=True)
        python.touch()
        history = self.root / ".paper-envs/main/conda-meta/history"
        history.parent.mkdir()
        history.touch()
        self.runner.conda = "/test/conda"
        with patch.object(self.runner, "command", return_value="") as command:
            self.runner.install_environment("main", specifications)
            self.assertEqual(command.call_count, 2)
            command.reset_mock()
            self.runner.install_environment("main", specifications)
            command.assert_not_called()
            correction.write_text(correction.read_text() + "# updated instructions\n")
            self.runner.install_environment("main", specifications)
            self.assertEqual(command.call_count, 2)

    def test_concurrent_launch_in_same_folder_is_rejected(self):
        import fcntl
        with (self.root / ".reproduce.lock").open("a") as lock, patch.object(reproduce.Runner, "execute") as execute, contextlib.redirect_stderr(io.StringIO()) as error:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            status = reproduce.main(["--root", str(self.root), "--check-only"])
            self.assertEqual(status, 2)
            self.assertIn("already running", error.getvalue())
            execute.assert_not_called()

    def test_preparation_runs_final_readiness_but_no_figures(self):
        with patch.object(self.runner, "check_computer"), patch.object(self.runner, "reconstruct"), patch.object(self.runner, "setup_software"), patch.object(self.runner, "readiness") as ready, patch.object(self.runner, "figures") as figures:
            self.runner.execute(prepare_only=True)
            ready.assert_called_once()
            figures.assert_not_called()

    def test_readiness_failure_prevents_figure_run(self):
        with patch.object(self.runner, "check_computer"), patch.object(self.runner, "reconstruct"), patch.object(self.runner, "setup_software"), patch.object(self.runner, "readiness", side_effect=reproduce.ReproductionError("bad mask")), patch.object(self.runner, "figures") as figures:
            with self.assertRaisesRegex(reproduce.ReproductionError, "bad mask"):
                self.runner.execute()
            figures.assert_not_called()

    def test_zero_exit_without_all_ten_completion_is_not_success(self):
        with patch.object(self.runner, "command", return_value="Only a preflight passed\n"):
            with self.assertRaisesRegex(reproduce.ReproductionError, "did not confirm completion"):
                self.runner.figures()

    def test_all_ten_run_uses_original_launcher_and_confirmation(self):
        with patch.object(self.runner, "command", return_value=reproduce.SUCCESS_MARKER + "\n") as command:
            self.runner.figures()
        self.assertEqual(command.call_args.args[1], ["bash", "./reproduce_all_figures.sh", "--output", self.root / "Paper", "--force"])
        self.assertEqual(command.call_args.kwargs["cwd"], self.root / "Paper")

    def test_old_success_in_appended_log_cannot_certify_current_run(self):
        self.runner.logs.mkdir()
        (self.runner.logs / "run.log").write_text(reproduce.SUCCESS_MARKER + "\n")
        output = self.runner.command("Test child", [sys.executable, "-c", "print('incomplete current run')"], "run.log")
        self.assertNotIn(reproduce.SUCCESS_MARKER, output)
        self.assertIn("incomplete current run", output)

    def test_real_nonzero_child_exit_is_preserved(self):
        with self.assertRaisesRegex(reproduce.ReproductionError, "exit 7"):
            self.runner.command("Test child", [sys.executable, "-c", "raise SystemExit(7)"], "failure.log")


if __name__ == "__main__":
    unittest.main()
