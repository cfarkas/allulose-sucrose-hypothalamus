"""The public scripts mirror must resolve the same Figure S3 data directory."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest


class FigureS3Paths(unittest.TestCase):
    def test_mirror_resolves_figure_data_in_relocated_tree(self):
        root = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory() as directory:
            paper = Path(directory) / "relocated_paper"
            expected = paper / "FigS3"
            for location in (expected, paper / "scripts" / "FigS3"):
                location.mkdir(parents=True)
                target = location / "figs3_common.py"
                shutil.copyfile(root / "FigS3/figs3_common.py", target)
                spec = importlib.util.spec_from_file_location("relocated_s3", target)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                self.assertEqual(module.PAPER_DIR, paper)
                self.assertEqual(module.FIGS3_DIR, expected)


if __name__ == "__main__":
    unittest.main()
