import hashlib
import importlib.util
import sys
import unittest
from pathlib import Path
from unittest.mock import patch


PAPER_ROOT = Path(__file__).resolve().parents[2]


def load_module(name: str, relative_path: str):
    path = PAPER_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class CleanRoomProvenancePathTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fig4_analyzer = load_module(
            "cleanroom_test_fig4_analyzer", "Fig4/02_analyze_pomc_cfos.py"
        )
        cls.fig4_renderer = load_module(
            "cleanroom_test_fig4_renderer", "Fig4/03_make_figure_4_pomc_cfos.py"
        )
        cls.fig5_analyzer = load_module(
            "cleanroom_test_fig5_analyzer", "Fig5/01_analyze_gfap_iba1_microglia.py"
        )
        cls.fig5_renderer = load_module(
            "cleanroom_test_fig5_renderer",
            "Fig5/02_make_figure_5_gfap_iba1_microglia.py",
        )

    def test_legacy_absolute_paths_are_never_resolved(self):
        legacy_fig4 = Path(
            "/unmounted-legacy/Paper/Fig4/hil_review/prepared/"
            "reconstructed_sections/FR6-1/FR6-1_S01/reconstructed_DAPI.tif"
        )
        current_fig4 = (
            PAPER_ROOT
            / "Fig4/analysis/reconstructed_sections/FR6-1/FR6-1_S01/"
            "reconstructed_DAPI.tif"
        )
        legacy_registration = legacy_fig4.with_name("stitch_registration.json")
        current_registration = current_fig4.with_name("stitch_registration.json")
        legacy_fig5 = Path(
            "/unmounted-legacy/Paper/Fig5/raw_data/FR6-1/XY01/FR6-1_XY01_CH1.tif"
        )
        current_fig5 = PAPER_ROOT / "Fig5/raw_data/FR6-1/XY01/FR6-1_XY01_CH1.tif"

        original_resolve = Path.resolve

        def guarded_resolve(path, *args, **kwargs):
            if str(path).startswith("/unmounted-legacy/"):
                raise AssertionError(f"legacy provenance path was probed: {path}")
            return original_resolve(path, *args, **kwargs)

        with patch.object(Path, "resolve", guarded_resolve):
            self.assertEqual(
                self.fig4_analyzer._receipt_path_identity(
                    legacy_fig4, current_fig4
                ),
                "relocated",
            )
            self.assertTrue(
                self.fig4_analyzer._registration_receipt_file_matches(
                    "stitch_registration.json", current_registration
                )
            )
            self.assertFalse(
                self.fig4_analyzer._registration_receipt_file_matches(
                    legacy_registration, current_registration
                )
            )
            normalized_fig4 = self.fig4_renderer._resolve_provenance_path(
                legacy_fig4,
                folder=current_fig4.parent,
                root=PAPER_ROOT,
            )
            self.assertEqual(normalized_fig4, legacy_fig4)
            self.assertTrue(
                self.fig4_renderer._same_reconstructed_section(
                    normalized_fig4, current_fig4
                )
            )
            self.assertEqual(
                self.fig5_analyzer._gi_resolve_relocated_figure_path(
                    legacy_fig5, repository_root=PAPER_ROOT
                ),
                current_fig5,
            )
            self.assertTrue(
                self.fig5_analyzer._gi_review_provenance_path_matches(
                    legacy_fig5,
                    current_fig5,
                    repository_root=PAPER_ROOT,
                )
            )
            self.assertTrue(
                self.fig5_renderer._review_provenance_path_matches(
                    legacy_fig5,
                    current_fig5,
                    repository_root=PAPER_ROOT,
                )
            )

    def test_active_front_and_mirror_scripts_match(self):
        pairs = (
            (
                "Fig4/02_analyze_pomc_cfos.py",
                "scripts/Fig4/02_analyze_pomc_cfos.py",
            ),
            (
                "Fig4/03_make_figure_4_pomc_cfos.py",
                "scripts/Fig4/03_make_figure_4_pomc_cfos.py",
            ),
            (
                "Fig5/01_analyze_gfap_iba1_microglia.py",
                "scripts/Fig5/01_analyze_gfap_iba1_microglia.py",
            ),
            (
                "Fig5/02_make_figure_5_gfap_iba1_microglia.py",
                "scripts/Fig5/02_make_figure_5_gfap_iba1_microglia.py",
            ),
        )
        for front, mirror in pairs:
            with self.subTest(front=front, mirror=mirror):
                self.assertEqual(sha256(PAPER_ROOT / front), sha256(PAPER_ROOT / mirror))


if __name__ == "__main__":
    unittest.main()
