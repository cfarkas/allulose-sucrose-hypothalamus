"""Check update installation boundaries and the eleven-figure completion contract."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import figure_updates as updates
import reproduce
from check_reproduction import Report, check_paper_sources


def sha(data):
    return hashlib.sha256(data).hexdigest()


class FigureUpdateTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.paper = self.root / "Paper"
        self.paper.mkdir()
        self.script = self.paper / "reproduce_all_figures.sh"
        self.script.write_bytes(b"# archived launcher\n")
        self.script.chmod(0o755)
        self.base = {"tree": {"files": [{"path": self.script.name, "sha256": sha(self.script.read_bytes()), "mode": "0755"}]}}
        (self.root / "manifest.json").write_text(json.dumps(self.base))
        rows = []
        for name, data in [(self.script.name, b"# eleven-figure launcher\n"), ("FigS6/raw_data/references.csv", b"cell,label\nc1,Rod-like\n")]:
            path = self.root / updates.PAYLOAD / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            mode = "0755" if name.endswith(".sh") else "0644"
            path.chmod(int(mode, 8))
            rows.append({"path": name, "size": len(data), "sha256": sha(data), "mode": mode})
        self.update = {"schema": updates.SCHEMA, "base_manifest_sha256": updates.digest(self.root / "manifest.json"), "files": rows}
        self.write_manifest()

    def write_manifest(self):
        (self.root / updates.MANIFEST).write_text(json.dumps(self.update))

    def test_install_preserves_archived_files_and_base_manifest_and_is_repeatable(self):
        before = (self.root / "manifest.json").read_bytes()
        result = updates.apply_updates(self.root)
        self.assertEqual((Path(result["backup"]) / self.script.name).read_bytes(), b"# archived launcher\n")
        self.assertEqual(self.script.read_bytes(), b"# eleven-figure launcher\n")
        self.assertEqual((self.root / "manifest.json").read_bytes(), before)
        self.assertEqual(updates.apply_updates(self.root)["status"], "already_installed")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertTrue(check_paper_sources(Report(), self.root, self.base))

    def test_unknown_local_edit_stops_before_any_installation(self):
        self.script.write_text("# user changes\n")
        with self.assertRaisesRegex(ValueError, "Preserving locally changed"):
            updates.apply_updates(self.root)
        self.assertFalse((self.paper / "FigS6").exists())
        self.assertEqual(self.script.read_text(), "# user changes\n")

    def test_bad_source_hash_stops_before_any_installation(self):
        (self.root / updates.PAYLOAD / "FigS6/raw_data/references.csv").write_text("changed")
        with self.assertRaisesRegex(ValueError, "absent or changed"):
            updates.apply_updates(self.root)
        self.assertEqual(self.script.read_bytes(), b"# archived launcher\n")
        self.assertFalse((self.paper / updates.RECEIPT).exists())

    def test_destination_link_stops_before_any_installation(self):
        (self.paper / "FigS6").symlink_to(self.root, target_is_directory=True)
        with self.assertRaisesRegex(ValueError, "Linked update"):
            updates.apply_updates(self.root)
        self.assertEqual(self.script.read_bytes(), b"# archived launcher\n")

    def test_traversal_and_duplicate_entries_are_rejected(self):
        self.update["files"].append(dict(self.update["files"][0]))
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            updates.load_updates(self.root)
        self.update["files"][-1]["path"] = "../outside.py"
        self.write_manifest()
        with self.assertRaisesRegex(ValueError, "Unsafe"):
            updates.load_updates(self.root)

    def test_installed_inputs_and_programs_are_checked_before_replay(self):
        updates.apply_updates(self.root)
        (self.paper / "FigS6/raw_data/references.csv").write_text("changed labels")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(check_paper_sources(Report(), self.root, self.base))
        (self.paper / "FigS6/raw_data/references.csv").write_bytes(b"cell,label\nc1,Rod-like\n")
        self.script.write_text("# changed program\n")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(check_paper_sources(Report(), self.root, self.base))

    def test_missing_or_changed_install_receipt_does_not_bypass_checks(self):
        updates.apply_updates(self.root)
        receipt = self.paper / updates.RECEIPT
        receipt.write_text('{"manifest_sha256":"changed"}')
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(check_paper_sources(Report(), self.root, self.base))
        receipt.unlink()
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(check_paper_sources(Report(), self.root, self.base))

    def test_archived_mode_does_not_claim_success_for_updated_tree(self):
        updates.apply_updates(self.root)
        runner = reproduce.Runner(self.root)
        with self.assertRaisesRegex(reproduce.ReproductionError, "--figure-updates"):
            runner.reconstruct()

    def test_microglial_flag_is_optional_and_forwarded_only_when_requested(self):
        default = reproduce.Runner(self.root, figure_updates=True)
        with patch.object(default, "command", return_value=reproduce.UPDATED_SUCCESS_MARKER + "\n") as command, contextlib.redirect_stdout(io.StringIO()):
            default.figures()
        self.assertNotIn("--do_microglial_choices", command.call_args.args[1])
        independent = reproduce.Runner(self.root, do_microglial_choices=True)
        with patch.object(independent, "command", return_value=reproduce.UPDATED_SUCCESS_MARKER + "\n") as command, contextlib.redirect_stdout(io.StringIO()):
            independent.figures()
        self.assertIn("--do_microglial_choices", command.call_args.args[1])
        self.assertTrue(independent.figure_updates)

    def test_beginner_command_defaults_to_current_choices(self):
        with patch.object(reproduce.Runner, "execute") as execute, patch.object(reproduce.Runner, "__init__", return_value=None) as init, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(reproduce.main(["--root", str(self.root), "--check-only"]), 0)
        self.assertTrue(init.call_args.kwargs["figure_updates"])
        self.assertFalse(init.call_args.kwargs["do_microglial_choices"])

    def test_review_url_is_visible_while_child_is_waiting(self):
        log = self.root / "progress.log"
        log.write_text("ordinary log line\nMICROGLIAL_REVIEW_URL: http://127.0.0.1:45678\n")
        runner = reproduce.Runner(self.root)
        with contextlib.redirect_stdout(io.StringIO()) as out:
            offset = runner.show_microglial_progress(log, 0)
            runner.show_microglial_progress(log, offset)
        self.assertEqual(out.getvalue().count("http://127.0.0.1:45678"), 1)
        self.assertNotIn("ordinary log", out.getvalue())

    def test_updated_mode_requires_eleven_figure_completion(self):
        runner = reproduce.Runner(self.root, figure_updates=True)
        with patch.object(runner, "command", return_value=reproduce.SUCCESS_MARKER + "\n"):
            with self.assertRaisesRegex(reproduce.ReproductionError, "did not confirm completion"):
                runner.figures()
        with patch.object(runner, "command", return_value=reproduce.UPDATED_SUCCESS_MARKER + "\n"), contextlib.redirect_stdout(io.StringIO()) as output:
            runner.figures()
        self.assertIn("All 11 figures reproduced and verified", output.getvalue())


if __name__ == "__main__":
    unittest.main()
