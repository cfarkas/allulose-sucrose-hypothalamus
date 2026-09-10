#!/usr/bin/env python3
"""Create uncluttered ARC/ME panels A-J and the combined publication figure.

Panels A, B and C are the isolated DAPI, c-FOS and POMC channels of the corrected
DAPI-registered field, each in the tint it contributes to the merge; their clipped
sum is exactly Panel E. Panel D is a channel-separated Cellpose cartoon row.
Panel E is that three-channel merge with its deterministic POMC-positive-cell
magnification inset. Panel F shows DAPI/c-FOS/NPY-GFP microscopy plus real POMC
intensity strictly masked to the registered POMC ROI interiors; its smaller
inset shows native POMC/NPY-GFP microscopy intensities inside the same ROIs,
using the same channel display settings, and adds dashed white tissue
and ventricular guides. Panels G-H show compact
animal-level results with vertical primary p-value blocks. Detailed raw ANOVA/MWU
methods and cage sensitivities are written to the legend and deterministic TXT/CSV
outputs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import shutil
import sys
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# Stabilize PDF CreationDate so command-level reruns are byte-reproducible.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt, gaussian_filter, label as ndi_label
from scipy.stats import f_oneway, mannwhitneyu, rankdata
from skimage.measure import label as component_label, regionprops
from skimage.morphology import binary_dilation, disk, remove_small_objects
from skimage.segmentation import find_boundaries
import tifffile


HERE = Path(__file__).resolve()
_PAPER_ROOT_OVERRIDE = os.environ.get("FIG4_PAPER_ROOT", "").strip()
PAPER_ROOT = (
    Path(_PAPER_ROOT_OVERRIDE).expanduser().resolve()
    if _PAPER_ROOT_OVERRIDE
    else next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                                  and (path / "README.txt").is_file()))), None)
)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")
PROJECT_ROOT = PAPER_ROOT
FIG4_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig4"
FIG4_FIGURE_ROOT = PAPER_ROOT / "Fig4"
FIG4_HIL_REVIEW_ROOT = FIG4_FIGURE_ROOT / "hil_review"
SHARED_SCRIPTS = PAPER_ROOT / "scripts" / "shared"
if str(SHARED_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SHARED_SCRIPTS))
from spanish_panel_text import translate_figure_texts_to_spanish
from pomc_size_qc import filter_pomc_by_local_dapi


def _first_existing(*candidates: Path) -> Path:
    """First candidate that exists, else the first, so defaults name real paths."""
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    return next(path for path in candidates if path is not None)


def _latest_accepted_hil_root() -> Optional[Path]:
    """Newest shipped human-final accepted review, if the tree carries one."""
    if not FIG4_HIL_REVIEW_ROOT.is_dir():
        return None
    roots = sorted(
        path for path in FIG4_HIL_REVIEW_ROOT.glob("human_final_*")
        if path.is_dir() and any(path.glob("*/*/*.json"))
    )
    return roots[-1] if roots else None


# The analyzer writes analyses/Fig4/results/human_final_run; results/final_run is
# the older working copy and is used only when it is the one that exists.
DEFAULT_ANALYSIS_ROOT = _first_existing(
    FIG4_ANALYSIS_ROOT / "results" / "human_final_run",
    FIG4_ANALYSIS_ROOT / "results" / "final_run",
)
DEFAULT_MANUAL_REGION_ROOT = _first_existing(
    FIG4_ANALYSIS_ROOT / "hil" / "accepted_annotations", _latest_accepted_hil_root(),
)
# Generated, so it is absent in a freshly unpacked tree and the documented command
# runs without --force. Promotion into Fig4/ stays a separate, deliberate copy.
DEFAULT_FIGURE_OUTPUT = FIG4_ANALYSIS_ROOT / "figure" / "human_final_render"
EXPECTED_POMC_QUANTIFIABLE_SECTIONS = 76
DEFAULT_REPRESENTATIVE_SAMPLE = "FR722"
DEFAULT_REPRESENTATIVE_SECTION = 2

# Figure 4 must never render from the former vertical-shift-only stitch. These
# values are a cross-script contract with 02_analyze_pomc_cfos.py. Bumping the
# analyzer algorithm therefore deliberately makes this maker fail closed until
# its registration audit is updated too.
STITCH_REGISTRATION_SCHEMA = "fig4_stitch_registration_v1"
STITCH_REGISTRATION_ALGORITHM = "dapi_integer_overlap_v1"
STITCH_REGISTRATION_STATUSES = {
    "verified_overlap", "manual_review_required", "single_half_native",
}
STITCH_REGISTRATION_PAIR_METHODS = {
    "sift_ransac_raw_ncc", "direct_ncc", "anatomical_y_fallback",
}

CONDITION_ORDER = ["Water", "Sucrose", "Allulose"]
CONDITION_COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
# ME moves off orange to brick so it cannot be read as amber POMC signal.
REGION_COLORS = {"ARC": "#00bcd4", "ME": "#c2410c", "VMN": "#be3eae"}
PHENOTYPE_COLORS = {
    "DAPI only": "#9e9ea8", "c-FOS+": "#ff05d1",
    "POMC+": "#b86e00", "c-FOS∧POMC": "#ffffff",
}
# Channel colours are chosen for separability in an additive merge on black.
# Blue for DAPI, magenta for c-FOS and green for NPY are fixed paper-wide, which
# leaves POMC as the only free channel, and amber is the one well-separated hue
# left: it differs from green in red, from magenta in green and blue, and from
# blue in everything. POMC was green until 2026-08-27, which collided with NPY.
# Cyan was considered and rejected because it shares the green axis with NPY,
# and the representative field carries both markers.
#
# Two shades are needed because the same channel is drawn on opposite
# backgrounds: bright amber reads on the black microscopy field, darker amber
# reads as a panel title on white and as a cartoon cell on the light tissue
# fill. A single value would be invisible on one of them.
CHANNEL_COLORS = {"DAPI": "#4f86ff", "POMC": "#b86e00", "c-FOS": "#ff37d4"}
POMC_MERGE_LABEL_COLOR = "#ffb43d"
# Green means NPY-GFP across the paper; Figure 3 uses this exact value.
NPY_GFP_COLOR = "#00e85e"
# Panel F shows the real registered POMC intensity only within accepted POMC
# ROI interiors, with exactly the Panel C channel display settings. Screen
# compositing avoids the additional clipping caused by adding POMC to the base.
# Main field and inset share intensity layers; ROIs never supply image color.
# This bright green is used only for text; image color comes from microscopy.
NPY_GFP_INSET_COLOR = "#00ff4f"
# Exact per-channel contributions summed by raw_composite. Panels A-C show these
# three views in isolation, so the merged field in Panel E is their clipped sum.
COMPOSITE_CHANNEL_TINTS = {
    "DAPI": (0.12, 0.28, 0.95),
    "c-FOS": (1.00, 0.02, 0.82),
    "POMC": (1.00, 0.65, 0.12),
}
COMPOSITE_CHANNEL_ORDER = ("DAPI", "c-FOS", "POMC")
# Panel F adds the independently normalized registered NPY-GFP plane with a
# screen blend. Red is exactly zero in this tint, so the red component that
# makes amber POMC conspicuous can never be reduced or overwritten by NPY-GFP.
NPY_GFP_COMPOSITE_TINT = (0.00, 1.00, 0.18)
# Panel letters run in reading order across the assembled figure. The Cellpose
# cartoon row sits directly under the A/B/C single-channel row it depicts, and
# its three cartoons are ordered DAPI, c-FOS, POMC to match them, so the merged
# field and the accepted ROIs follow as E and F. The suffixes name the content,
# so a re-lettering does not silently rename a panel to another panel's file.
PANEL_FILE_SUFFIXES = {
    "A": "dapi", "B": "cfos", "C": "pomc",
    "D": "cellpose_cartoon", "E": "microscopy", "F": "microscopy_npy_gfp",
    "G": "cfos_dapi", "H": "pomc_activation",
    "I": "spatial_cfos", "J": "spatial_cfos_pomc",
}
# Cartoon order in panel D, matching panels A, B and C.
CARTOON_MARKER_ORDER = ("DAPI", "c-FOS", "POMC")
CELLPOSE_MARKERS = {"DAPI": "dapi", "POMC": "pomc", "c-FOS": "cfos", "NPY": "npy"}
# Panels I and J are rendered by Fig4/05_analyze_spatial_distributions.py, which
# writes them under these names.
SPATIAL_PANEL_NAMES = ("Figure4_Spatial_cFOS_occurrence.png",
                       "Figure4_Spatial_cFOS_POMC_occurrence.png")
SPANISH_PANEL_TEXT = {
    "POMC (within ROIs)": "POMC (en ROIs)",
    "Tissue boundary": "Límite tisular",
    "Cellpose segmentation cartoon of the representative field":
        "Segmentación Cellpose del campo representativo",
    "c-FOS/DAPI per animal": "c-FOS/DAPI por animal",
    "POMC activation per animal": "Activación de POMC por animal",
    "c-FOS/DAPI ratio": "Razón c-FOS/DAPI",
    "c-FOS∧POMC /\nPOMC ratio": "Razón c-FOS∧POMC /\nPOMC",
    "Mean ± SD": "Media ± DE",
    "mean ± SD": "media ± DE",
    "one-way ANOVA": "ANOVA de una vía",
    "Pairwise exact MWU": "MWU exacta por pares",
    "exact MWU": "MWU exacta",
    "segmentation": "segmentación",
    "animal": "animal",
}


def spanish_spatial_variants(paths: Sequence[Path]) -> tuple[Path, Path]:
    variants = tuple(
        path.with_name(f"{path.stem}_spanish{path.suffix}") for path in paths
    )
    missing = [str(path) for path in variants if not path.is_file()]
    if missing:
        raise SystemExit(
            "Spanish spatial panels not found: " + ", ".join(missing)
            + ". Re-run Fig4/05_analyze_spatial_distributions.py at 600 dpi."
        )
    return variants


def resolve_spatial_panels(args) -> tuple[Path, Path]:
    """Locate the two spatial occurrence panels drawn as I and J.

    Explicit --panel-i/--panel-j win; otherwise they are read from --spatial-dir,
    which defaults to a "spatial" directory beside the analysis this figure is
    built from. Missing panels are a hard error rather than a silently shorter
    figure, because the panel letters are cited by the manuscript.
    """
    explicit = (args.panel_i, args.panel_j)
    if any(path is not None for path in explicit):
        if not all(path is not None for path in explicit):
            raise SystemExit("Pass both --panel-i and --panel-j, or neither")
        resolved = tuple(Path(path).expanduser().resolve() for path in explicit)
    else:
        spatial_dir = (Path(args.spatial_dir).expanduser().resolve()
                       if args.spatial_dir is not None
                       else Path(args.analysis_dir).expanduser().resolve() / "spatial")
        resolved = tuple(spatial_dir / name for name in SPATIAL_PANEL_NAMES)
    missing = [str(path) for path in resolved if not path.is_file()]
    if missing:
        raise SystemExit(
            "Spatial occurrence panels not found: " + ", ".join(missing)
            + ". Run Fig4/05_analyze_spatial_distributions.py first, or pass "
              "--spatial-dir / --panel-i / --panel-j."
        )
    return resolved
TISSUE_FILL = "#e8edf2"
TISSUE_EDGE = "#354052"


def _prepare_generated_output_dir(
    output: Path,
    *,
    force: bool,
    protected_paths: Sequence[Path],
    preserved_names: Sequence[str] = (),
) -> Path:
    """Refuse mixed publication generations, or replace one safe output tree."""
    requested = Path(output).expanduser()
    if requested.exists() and requested.is_symlink():
        raise SystemExit(f"Refusing publication output through a symbolic link: {requested}")
    resolved = requested.resolve()
    filesystem_root = Path(resolved.anchor).resolve()
    protected = [Path(path).expanduser().resolve() for path in protected_paths]
    if resolved == filesystem_root:
        raise SystemExit(f"Refusing filesystem root as publication output: {resolved}")
    for path in protected:
        if resolved == path:
            raise SystemExit(f"Refusing publication output equal to protected path: {resolved}")
        if resolved in path.parents:
            raise SystemExit(
                f"Refusing publication output that contains protected path {path}: {resolved}"
            )
    if resolved.exists() and not resolved.is_dir():
        raise SystemExit(f"Publication output exists and is not a directory: {resolved}")
    nonempty = resolved.exists() and any(resolved.iterdir())
    if nonempty and not force:
        raise SystemExit(
            f"Publication output already exists and is not empty: {resolved}\n"
            "Use --force to replace only this validated output directory."
        )
    if force and resolved.exists():
        preserved = set(preserved_names)
        for child in resolved.iterdir():
            if child.name in preserved:
                continue
            if child.is_symlink() or not child.is_dir():
                child.unlink()
            else:
                shutil.rmtree(child)
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Generate individual A-E PDFs and the combined ARC/ME figure.")
    p.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    p.add_argument(
        "--manual-region-dir", type=Path, default=DEFAULT_MANUAL_REGION_ROOT,
        help=("Root containing accepted native POMC ARC/ME/VMN TIFF+JSON pairs. "
              "Complete 76/76 accepted coverage is mandatory; each receipt may "
              "contain any nonempty subset of ARC, ME, and VMN."),
    )
    p.add_argument(
        "--require-manual-regions", "--human-in-the-loop", "--human_in_the_loop",
        dest="require_manual_regions", action="store_true", default=True,
        help="Compatibility flag; complete native POMC HIL provenance is always required.",
    )
    p.add_argument("--outdir", type=Path, default=DEFAULT_FIGURE_OUTPUT)
    p.add_argument(
        "--force", action="store_true",
        help=("Replace the designated nonempty publication output directory after strict path-safety checks. "
              "Analysis inputs and project data roots are never deleted."),
    )
    p.add_argument("--figure-name", default="Figure_ARC_ME_multipanel")
    p.add_argument("--um-per-px", type=float, default=0.755)
    p.add_argument("--scalebar-um", type=int, choices=[20, 100], default=100)
    p.add_argument("--inset-field-um", type=float, default=50.0)
    p.add_argument("--inset-scalebar-um", type=float, default=20.0)
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--spatial-dir", type=Path, default=None,
                   help=("Output directory of 05_analyze_spatial_distributions.py, holding "
                         "the two spatial occurrence panels rendered as panels I and J. "
                         "Defaults to <analysis-dir>/spatial."))
    p.add_argument("--panel-i", type=Path, default=None,
                   help="Explicit path to the c-FOS spatial occurrence panel PNG.")
    p.add_argument("--panel-j", type=Path, default=None,
                   help="Explicit path to the c-FOS+POMC spatial occurrence panel PNG.")
    p.add_argument("--seed", type=int, default=20260812)
    p.add_argument(
        "--representative-sample", default=DEFAULT_REPRESENTATIVE_SAMPLE,
        help="Original reviewed representative (default: FR722).",
    )
    p.add_argument(
        "--representative-section", type=int,
        default=DEFAULT_REPRESENTATIVE_SECTION,
        help="Representative section index (default: 2).",
    )
    return p


def as_bool(values: pd.Series) -> pd.Series:
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def clean_tile(value: object) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return ""
    return text[:-2] if text.endswith(".0") and text[:-2].isdigit() else text


def consistent_group_metadata(
    frame: pd.DataFrame,
    group_columns: Sequence[str],
    metadata_columns: Sequence[str],
    *,
    context: str,
) -> pd.DataFrame:
    """Resolve one optional text value per group while rejecting conflicts."""
    rows: List[Dict[str, object]] = []
    for keys, group in frame.groupby(
        list(group_columns), dropna=False, sort=False,
    ):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        record: Dict[str, object] = dict(zip(group_columns, key_values))
        for column in metadata_columns:
            values = []
            for value in group[column]:
                if pd.isna(value):
                    continue
                text = str(value).strip()
                if text and text not in values:
                    values.append(text)
            if len(values) > 1:
                raise ValueError(
                    f"Conflicting {column} values for {context} {keys}: {values}"
                )
            record[column] = values[0] if values else pd.NA
        rows.append(record)
    return pd.DataFrame(rows, columns=[*group_columns, *metadata_columns])


def section_directory(root: Path, row: Mapping[str, object]) -> Path:
    tile = clean_tile(row.get("tile", ""))
    suffix = f"_T{tile}" if tile else ""
    return root / "reconstructed_sections" / str(row["animal_id"]) / f"{row['sample']}_S{int(row['section_index']):02d}{suffix}"

def microscopy_sources(folder: Path) -> Dict[str, Path]:
    return {
        "dapi": folder / "reconstructed_DAPI.tif",
        "cfos": folder / "reconstructed_cFOS.tif",
        "pomc": folder / "reconstructed_POMC.tif",
        "npy": folder / "reconstructed_NPY.tif",
        "dapi_labels": folder / "reconstructed_DAPI_labels.tif",
        "cfos_labels": folder / "reconstructed_cFOS_labels.tif",
        "pomc_labels": folder / "reconstructed_POMC_labels.tif",
        "npy_labels": folder / "reconstructed_NPY_labels.tif",
        "native_proposal": folder / "reconstructed_ARC_ME_region_mask.tif",
        "accepted_regions": folder / "reconstructed_ARC_ME_accepted_mask.tif",
        "walls": folder / "reconstructed_projected_ventricle_walls.csv",
        "floor": folder / "reconstructed_ventricle_floor_curve.csv",
        "registration": folder / "stitch_registration.json",
    }



def load_size_filtered_pomc(sources):
    original = read_gray(sources["pomc_labels"]).astype(np.int32)
    dapi = read_gray(sources["dapi_labels"]).astype(np.int32)
    regions = read_gray(sources["accepted_regions"]).astype(np.uint8)
    expected, audit = filter_pomc_by_local_dapi(original, dapi, regions)
    path = sources["pomc_labels"].with_name("reconstructed_POMC_size_filtered_labels.tif")
    receipt_path = path.with_name("pomc_size_qc.json")
    if not path.is_file() or not receipt_path.is_file():
        raise ValueError("POMC analysis predates local DAPI size QC; rerun the analyzer")
    receipt = json.loads(receipt_path.read_text())
    actual = read_gray(path).astype(np.int32)
    if not np.array_equal(actual, expected) or receipt.get("filtered_mask_sha256") != sha256_file(path):
        raise ValueError("POMC size-filtered mask differs from the deterministic regional rule")
    if receipt.get("original_mask_sha256") != sha256_file(sources["pomc_labels"]):
        raise ValueError("POMC size-QC receipt is not bound to the original masks")
    return actual, audit, path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _valid_sha256(value: object) -> bool:
    text = str(value).strip().lower()
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def _timezone_aware(value: object) -> bool:
    try:
        parsed = dt.datetime.fromisoformat(str(value).strip().replace("Z", "+00:00"))
        return parsed.tzinfo is not None
    except Exception:
        return False


def _strict_bool(value: object, *, field: str) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"true", "1"}:
        return True
    if text in {"false", "0"}:
        return False
    raise ValueError(f"{field} must be an explicit boolean, got {value!r}")


def _finite_int(value: object, *, field: str) -> int:
    number = float(value)
    if not np.isfinite(number) or not number.is_integer():
        raise ValueError(f"{field} must be a finite integer, got {value!r}")
    return int(number)


def _optional_float(value: object, *, field: str) -> float:
    """Parse an optional numeric provenance value, mapping JSON null to NaN.

    Registration methods do not all emit every NCC diagnostic. The analyzer
    serializes unavailable optional metrics as JSON ``null`` and writes them as
    blank CSV cells. Treating a present-but-null JSON key with ``float(None)``
    made the renderer reject otherwise valid registered sections. Missing/null
    values remain auditable as NaN; malformed or infinite values still fail.
    """
    if value is None:
        return float("nan")
    if isinstance(value, str):
        stripped = value.strip()
        if stripped.lower() in {"", "na", "nan", "none", "null"}:
            return float("nan")
        value = stripped
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be numeric or null, got {value!r}") from exc
    if np.isnan(number):
        return float("nan")
    if not np.isfinite(number):
        raise ValueError(f"{field} must be finite or null, got {value!r}")
    return number


def _lexical_absolute_path(path: Path, *, base: Optional[Path] = None) -> Path:
    """Normalize provenance text without touching an historical filesystem."""
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = Path(base) / candidate if base is not None else Path.cwd() / candidate
    return Path(os.path.abspath(os.fspath(candidate)))


def _resolve_provenance_path(value: object, *, folder: Path, root: Path) -> Path:
    raw = Path(str(value).strip()).expanduser()
    if raw.is_absolute():
        normalized = _lexical_absolute_path(raw)
        local_roots = tuple(
            _lexical_absolute_path(path) for path in (folder, root, PAPER_ROOT)
        )
        if any(
            normalized == local or normalized.is_relative_to(local)
            for local in local_roots
        ):
            return normalized.resolve()
        # Accepted receipts retain their original absolute stage path.  It is
        # provenance text only; callers bind it to the current section by its
        # path tail plus cryptographic hashes, so never probe that old tree.
        return normalized
    candidates = (folder / raw, root / raw, PAPER_ROOT / raw)
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return _lexical_absolute_path(candidates[0])


# <animal_id>/<sample>_S##/<leaf file>: the tail that identifies one reconstructed
# section inside any preparation or analysis stage.
RECEIPT_SECTION_PATH_DEPTH = 3


def _same_reconstructed_section(receipt_path: Path, expected_path: Path) -> bool:
    """True when both paths name the same reconstructed section, in any stage.

    A receipt records the directory the review images were served from, which is
    the immutable prepared review input; the analysis row records the directory
    this run rebuilt them into. Those are the same composite in two places, and
    a downloaded tree makes them differ by construction. The caller has already
    required the reviewed DAPI SHA-256 to equal the analyzed source DAPI's, so
    the stage prefix carries no identity that the hash does not.
    """
    if receipt_path == expected_path:
        return True
    return (
        receipt_path.parts[-RECEIPT_SECTION_PATH_DEPTH:]
        == expected_path.parts[-RECEIPT_SECTION_PATH_DEPTH:]
    )


def _registration_output_shape(record: Mapping[str, object]) -> Tuple[int, int]:
    raw = record.get("output_shape_px")
    if isinstance(raw, Mapping):
        height = raw.get("height", raw.get("height_px"))
        width = raw.get("width", raw.get("width_px"))
        return (
            _finite_int(height, field="output_shape_px.height"),
            _finite_int(width, field="output_shape_px.width"),
        )
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)) and len(raw) == 2:
        return (
            _finite_int(raw[0], field="output_shape_px[0]"),
            _finite_int(raw[1], field="output_shape_px[1]"),
        )
    raise ValueError("output_shape_px must contain registered height and width")


def _same_flat_value(observed: object, expected: object, *, field: str) -> None:
    if isinstance(expected, bool):
        if _strict_bool(observed, field=field) != expected:
            raise ValueError(f"{field} differs between summary and registration JSON")
        return
    if isinstance(expected, (int, float, np.integer, np.floating)):
        left = float(observed); right = float(expected)
        if not (np.isfinite(left) and np.isfinite(right) and np.isclose(left, right, atol=1e-9, rtol=0)):
            raise ValueError(
                f"{field} differs between summary ({observed!r}) and registration JSON ({expected!r})"
            )
        return
    if str(observed).strip() != str(expected).strip():
        raise ValueError(
            f"{field} differs between summary ({observed!r}) and registration JSON ({expected!r})"
        )


def validate_registered_reconstructions(
    quantified: pd.DataFrame, root: Path,
) -> Tuple[Dict[Tuple[str, str, int, str], Dict[str, object]], int]:
    """Verify every reconstructed section uses the current DAPI overlap contract.

    The current JSON is authoritative for the analyzed reconstruction and its
    exact bytes are linked from the section summary. A historical HIL receipt
    may carry an earlier JSON byte hash only when it still binds the identical
    registered DAPI TIFF, dimensions, section identity, and confirmation. Thus
    neither an old vertical-only reconstruction nor a replaced DAPI field can
    feed panels A/B or the regional statistics.
    """
    required = {
        "animal_id", "sample", "section_index", "tile", "has_left", "has_right",
        "has_pomc_channel", "has_npy_channel",
        "stitch_alignment_algorithm_version", "stitch_alignment_status",
        "stitch_alignment_requires_manual_review", "stitch_alignment_method",
        "stitch_alignment_right_origin_x_px", "stitch_alignment_shift_right_px",
        "stitch_alignment_right_medial_crop_px", "stitch_alignment_overlap_width_px",
        "stitch_alignment_raw_overlap_ncc", "stitch_alignment_gradient_overlap_ncc",
        "stitch_alignment_ncc_best_score", "stitch_alignment_ncc_second_score",
        "stitch_alignment_ncc_uniqueness_margin", "stitch_registration_path",
        "stitch_registration_sha256", "manual_region_review_path",
    }
    missing = sorted(required.difference(quantified.columns))
    if missing:
        raise ValueError(
            "Figure 4 analysis predates corrected DAPI overlap registration; "
            "rerun Fig4/02 and HIL before rendering. Missing fields: "
            + ", ".join(missing)
        )

    registrations: Dict[Tuple[str, str, int, str], Dict[str, object]] = {}
    newest_registration_mtime_ns = 0
    problems: List[str] = []
    for row in quantified.to_dict("records"):
        identity = (
            str(row["animal_id"]), str(row["sample"]),
            int(row["section_index"]), clean_tile(row.get("tile", "")),
        )
        identity_text = f"{identity[0]}/{identity[1]}_S{identity[2]:02d}"
        folder = section_directory(root, row)
        registration_path = folder / "stitch_registration.json"
        try:
            if not registration_path.is_file():
                raise FileNotFoundError(registration_path)
            row_registration_path = _resolve_provenance_path(
                row["stitch_registration_path"], folder=folder, root=root,
            )
            if row_registration_path != registration_path.resolve():
                raise ValueError(
                    "section summary points to a different stitch_registration.json"
                )
            registration_sha = sha256_file(registration_path)
            if (
                not _valid_sha256(row["stitch_registration_sha256"])
                or str(row["stitch_registration_sha256"]).strip().lower() != registration_sha
            ):
                raise ValueError("stitch registration JSON SHA-256 mismatch")
            record = json.loads(registration_path.read_text(encoding="utf-8"))
            if not isinstance(record, dict):
                raise ValueError("registration JSON root must be an object")
            schema = str(record.get("schema_version", ""))
            algorithm = str(record.get("algorithm", ""))
            status = str(record.get("status", ""))
            method = str(record.get("method", ""))
            requires_review = _strict_bool(
                record.get("requires_manual_review"), field="requires_manual_review",
            )
            if schema != STITCH_REGISTRATION_SCHEMA:
                raise ValueError(f"unsupported registration schema {schema!r}")
            if algorithm != STITCH_REGISTRATION_ALGORITHM:
                raise ValueError(f"legacy/unsupported registration algorithm {algorithm!r}")
            if status not in STITCH_REGISTRATION_STATUSES:
                raise ValueError(f"unsupported registration status: {status!r}")

            # Every new ROI receipt binds the exact registration JSON.  A
            # low-confidence pair additionally needs explicit image-by-image
            # seam confirmation. This deliberately cannot be satisfied by an
            # older ROI receipt created against the former y-only stitch.
            review_path = Path(
                str(row.get("manual_region_review_path", ""))
            ).expanduser().resolve()
            if not review_path.is_file():
                raise ValueError("registered composite lacks its accepted HIL receipt")
            review = json.loads(review_path.read_text(encoding="utf-8"))
            expected_confirmation = {
                "verified_overlap": "automated_overlap_qc_verified",
                "manual_review_required": "reviewed_registered_composite_image_by_image",
                "single_half_native": "single_half_native_no_registration_needed",
            }[status]
            registration_output_for_receipt = record.get("output")
            receipt_registration_sha = str(
                review.get("stitch_registration_sha256", "")
            ).strip().lower() if isinstance(review, dict) else ""
            exact_receipt_registration_hash = (
                _valid_sha256(receipt_registration_sha)
                and receipt_registration_sha == registration_sha
            )
            registered_dapi_sha = (
                str(registration_output_for_receipt.get("dapi_sha256", "")).strip().lower()
                if isinstance(registration_output_for_receipt, dict) else ""
            )
            receipt_dapi_sha = (
                str(review.get("source_dapi_sha256", "")).strip().lower()
                if isinstance(review, dict) else ""
            )
            registered_dapi_path = folder / "reconstructed_DAPI.tif"
            if not registered_dapi_path.is_file():
                raise ValueError("registered DAPI TIFF is missing")
            registered_dapi_array = np.asarray(tifffile.imread(registered_dapi_path))
            registered_dapi_shape = tuple(
                int(v) for v in registered_dapi_array.shape[:2]
            )
            registered_dapi_file_sha = sha256_file(registered_dapi_path)
            reviewed_composite_identical = bool(
                _valid_sha256(receipt_registration_sha)
                and _valid_sha256(receipt_dapi_sha)
                and receipt_dapi_sha == registered_dapi_sha
                and receipt_dapi_sha == registered_dapi_file_sha
                and _registration_output_shape(record) == registered_dapi_shape
            )
            registration_receipt_bound = bool(
                isinstance(review, dict)
                and review.get("accepted") is True
                and str(review.get("stitch_registration_schema_version", "")) == schema
                and (exact_receipt_registration_hash or reviewed_composite_identical)
                and str(review.get("stitch_registration_status", "")) == status
                and str(review.get("stitch_registration_method", "")) == method
                and review.get("stitch_registration_requires_manual_review") is requires_review
                and review.get("stitch_registration_confirmation") == expected_confirmation
            )
            if not registration_receipt_bound:
                raise ValueError(
                    "HIL receipt does not bind the current registration or the exact "
                    "registered DAPI composite"
                )
            if requires_review != (status == "manual_review_required"):
                raise ValueError(
                    "registration status/requires_manual_review fields are contradictory"
                )
            manual_registration_confirmed = False
            if status == "manual_review_required":
                manual_registration_confirmed = bool(
                    review.get("registered_composite_reviewed") is True
                    and review.get("stitch_registration_confirmation")
                        == "reviewed_registered_composite_image_by_image"
                    and (exact_receipt_registration_hash or reviewed_composite_identical)
                )
                if not manual_registration_confirmed:
                    raise ValueError(
                        "uncertain registration was not manually confirmed against the "
                        "current registered DAPI composite"
                    )

            has_left = _strict_bool(row["has_left"], field="has_left")
            has_right = _strict_bool(row["has_right"], field="has_right")
            if has_left and has_right:
                if status not in {"verified_overlap", "manual_review_required"}:
                    raise ValueError("paired section lacks an overlap-registration status")
                if method not in STITCH_REGISTRATION_PAIR_METHODS:
                    raise ValueError(f"unsupported paired registration method {method!r}")
            elif has_left != has_right:
                if status != "single_half_native" or method != "single_half_native":
                    raise ValueError("single-half section lacks single_half_native provenance")
            else:
                raise ValueError("section has neither left nor right DAPI half")

            transform = record.get("transform")
            qc = record.get("qc")
            output = record.get("output")
            if not isinstance(transform, dict) or not isinstance(qc, dict) or not isinstance(output, dict):
                raise ValueError("registration JSON needs transform, qc, and output objects")
            right_origin = _finite_int(
                transform.get("right_origin_x_px"), field="transform.right_origin_x_px",
            )
            shift_y = _finite_int(
                transform.get("right_shift_y_px"), field="transform.right_shift_y_px",
            )
            right_crop = _finite_int(
                transform.get("right_medial_crop_px"), field="transform.right_medial_crop_px",
            )
            overlap = _finite_int(
                transform.get("overlap_width_px"), field="transform.overlap_width_px",
            )
            gap = _finite_int(transform.get("gap_px"), field="transform.gap_px")
            if right_origin < 0 or right_crop < 0 or overlap < 0 or gap < 0:
                raise ValueError("registration crop/origin/overlap/gap values must be nonnegative")
            if _strict_bool(transform.get("resampling"), field="transform.resampling"):
                raise ValueError("registered integer stitch unexpectedly resampled microscopy")
            orientation = str(transform.get("orientation", ""))
            if orientation not in {"left-right", "right-left"}:
                raise ValueError(f"unsupported stitch orientation {orientation!r}")

            registration_input = record.get("input")
            if not isinstance(registration_input, dict):
                raise ValueError("registration JSON lacks input DAPI provenance")
            registered_input_hashes: Dict[str, str] = {}
            for side, present in (("left", has_left), ("right", has_right)):
                input_path_text = str(registration_input.get(f"{side}_dapi_path", "")).strip()
                input_sha = str(registration_input.get(f"{side}_dapi_sha256", "")).strip().lower()
                input_shape = registration_input.get(f"{side}_dapi_shape_px")
                if not present:
                    if input_path_text or input_sha or input_shape not in ([], None):
                        raise ValueError(f"absent {side} half has contradictory DAPI provenance")
                    continue
                input_path = _resolve_provenance_path(
                    input_path_text, folder=folder, root=root,
                )
                if not input_path.is_file():
                    raise ValueError(f"registered {side} input DAPI is missing")
                if not _valid_sha256(input_sha) or sha256_file(input_path) != input_sha:
                    raise ValueError(f"registered {side} input DAPI SHA-256 mismatch")
                if not isinstance(input_shape, list) or len(input_shape) != 2:
                    raise ValueError(f"registered {side} input DAPI shape is missing")
                if list(read_gray(input_path).shape) != [int(value) for value in input_shape]:
                    raise ValueError(f"registered {side} input DAPI shape/provenance mismatch")
                registered_input_hashes[side] = input_sha

            dapi_path = (folder / "reconstructed_DAPI.tif").resolve()
            output_dapi_path = _resolve_provenance_path(
                output.get("dapi_path", ""), folder=folder, root=root,
            )
            if output_dapi_path != dapi_path or not dapi_path.is_file():
                raise ValueError("registration output DAPI path is not this reconstructed section")
            dapi_sha = sha256_file(dapi_path)
            output_dapi_sha = str(output.get("dapi_sha256", "")).strip().lower()
            if not _valid_sha256(output_dapi_sha) or output_dapi_sha != dapi_sha:
                raise ValueError("registered output DAPI TIFF SHA-256 mismatch")
            output_shape = _registration_output_shape(record)
            dapi = read_gray(dapi_path)
            if tuple(dapi.shape) != tuple(output_shape):
                raise ValueError(
                    f"registered DAPI shape {dapi.shape} != provenance {output_shape}"
                )

            # Prove that Panel D channels and every segmentation/count raster
            # were emitted by the same immutable registration transaction.
            required_artifacts = ["dapi_labels", "cfos", "cfos_labels"]
            if _strict_bool(row.get("has_pomc_channel", False), field="has_pomc_channel"):
                required_artifacts.extend(["pomc", "pomc_labels"])
            if _strict_bool(row.get("has_npy_channel", False), field="has_npy_channel"):
                required_artifacts.extend(["npy", "npy_labels"])
            registered_artifact_hashes = {"dapi": dapi_sha}
            for prefix in required_artifacts:
                marker, is_label = (
                    prefix.removesuffix("_labels"), prefix.endswith("_labels")
                )
                marker_name = {"dapi": "DAPI", "cfos": "cFOS", "pomc": "POMC", "npy": "NPY"}[marker]
                expected_path = (
                    folder / f"reconstructed_{marker_name}{'_labels' if is_label else ''}.tif"
                ).resolve()
                artifact_path = _resolve_provenance_path(
                    output.get(f"{prefix}_path", ""), folder=folder, root=root,
                )
                if artifact_path != expected_path or not artifact_path.is_file():
                    raise ValueError(f"registered output {prefix} path is missing/different")
                artifact_sha = sha256_file(artifact_path)
                recorded_sha = str(output.get(f"{prefix}_sha256", "")).strip().lower()
                if not _valid_sha256(recorded_sha) or recorded_sha != artifact_sha:
                    raise ValueError(f"registered output {prefix} TIFF SHA-256 mismatch")
                recorded_shape = output.get(f"{prefix}_shape_px")
                if not isinstance(recorded_shape, list) or len(recorded_shape) < 2:
                    raise ValueError(f"registered output {prefix} shape is missing")
                artifact_shape = tuple(
                    int(value) for value in np.asarray(tifffile.imread(str(artifact_path))).shape
                )
                if list(artifact_shape) != [int(value) for value in recorded_shape]:
                    raise ValueError(f"registered output {prefix} TIFF shape/provenance mismatch")
                if tuple(artifact_shape[-2:]) != tuple(output_shape):
                    raise ValueError(
                        f"registered output {prefix} does not share DAPI pixel bounds"
                    )
                registered_artifact_hashes[prefix] = artifact_sha

            qc_values = {
                key: _optional_float(qc.get(key), field=f"qc.{key}")
                for key in (
                    "raw_overlap_ncc",
                    "gradient_overlap_ncc",
                    "ncc_best_score",
                    "ncc_second_score",
                    "ncc_uniqueness_margin",
                )
            }
            if status == "verified_overlap":
                if not (
                    np.isfinite(qc_values["raw_overlap_ncc"])
                    or np.isfinite(qc_values["gradient_overlap_ncc"])
                ):
                    raise ValueError("verified overlap lacks raw/gradient overlap NCC")
                if method == "direct_ncc" and not np.isfinite(qc_values["ncc_best_score"]):
                    raise ValueError("direct-NCC overlap lacks a finite best NCC score")

            flat_expected = {
                "stitch_alignment_algorithm_version": algorithm,
                "stitch_alignment_status": status,
                "stitch_alignment_requires_manual_review": requires_review,
                "stitch_alignment_method": method,
                "stitch_alignment_right_origin_x_px": right_origin,
                "stitch_alignment_shift_right_px": shift_y,
                "stitch_alignment_right_medial_crop_px": right_crop,
                "stitch_alignment_overlap_width_px": overlap,
                "stitch_alignment_raw_overlap_ncc": qc_values["raw_overlap_ncc"],
                "stitch_alignment_gradient_overlap_ncc": qc_values["gradient_overlap_ncc"],
                "stitch_alignment_ncc_best_score": qc_values["ncc_best_score"],
                "stitch_alignment_ncc_second_score": qc_values["ncc_second_score"],
                "stitch_alignment_ncc_uniqueness_margin": qc_values["ncc_uniqueness_margin"],
            }
            for field, expected in flat_expected.items():
                # Single-half rows may correctly report undefined NCC values.
                if isinstance(expected, float) and not np.isfinite(expected):
                    if pd.notna(row[field]):
                        raise ValueError(f"{field} must be blank for a single-half section")
                else:
                    _same_flat_value(row[field], expected, field=field)

            record = dict(record)
            record.update({
                "_registration_path": str(registration_path.resolve()),
                "_registration_sha256": registration_sha,
                "_dapi_path": str(dapi_path),
                "_dapi_sha256": dapi_sha,
                "_artifact_sha256s": registered_artifact_hashes,
                "_input_dapi_sha256s": registered_input_hashes,
                "_output_shape_px": output_shape,
                "_right_origin_x_px": right_origin,
                "_right_shift_y_px": shift_y,
                "_right_medial_crop_px": right_crop,
                "_overlap_width_px": overlap,
                "_gap_px": gap,
                "_registration_receipt_bound": registration_receipt_bound,
                "_manual_registration_confirmed": manual_registration_confirmed,
                **{f"_{key}": value for key, value in qc_values.items()},
            })
            registrations[identity] = record
            newest_registration_mtime_ns = max(
                newest_registration_mtime_ns,
                registration_path.stat().st_mtime_ns,
                dapi_path.stat().st_mtime_ns,
            )
        except Exception as exc:
            problems.append(f"{identity_text}: {exc}")

    if problems:
        detail = "\n".join(f"  - {item}" for item in problems[:8])
        if len(problems) > 8:
            detail += f"\n  - ... and {len(problems) - 8} more"
        raise ValueError(
            f"Corrected Figure 4 registration audit failed for {len(problems)} of "
            f"{len(quantified)} quantifiable sections.\n{detail}"
        )
    if len(registrations) != len(quantified):
        raise ValueError("Registration identities are not one-to-one with quantifiable sections")
    return registrations, newest_registration_mtime_ns


def recompute_registered_statistical_inputs(
    *, root: Path, quantified: pd.DataFrame,
    newest_registration_mtime_ns: int, required_fresh_paths: Sequence[Path],
) -> Tuple[pd.DataFrame, pd.DataFrame, Path]:
    """Reaggregate Figure 4 endpoints from current registered section counts.

    No pre-existing plot-values or definitive-endpoint CSV is trusted for the
    publication render.  The persisted per-animal table is audited against this
    independent reaggregation, then panels C/D and all exact tests consume the
    newly calculated in-memory tables.
    """
    section_region_path = root / "per_reconstructed_section_region_summary.csv"
    per_animal_path = root / "per_animal_region_summary.csv"
    for path in [*required_fresh_paths, section_region_path, per_animal_path]:
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_mtime_ns < newest_registration_mtime_ns:
            raise ValueError(
                "Statistics/inset input predates the corrected registered DAPI output: "
                f"{path}. Rerun the complete Figure 4 analysis after registration/HIL."
            )

    section_region = pd.read_csv(section_region_path, low_memory=False)
    required = {
        "sample", "animal_id", "animal_key", "section_index", "tile", "region",
        "dapi_nuclei", "cfos_cells", "pomc_cells", "cfos_pomc_all",
        "pomc_assay_dapi_nuclei", "cond", "genotype", "sex", "cage",
    }
    missing = sorted(required.difference(section_region.columns))
    if missing:
        raise ValueError(
            "Registered section-region summary lacks statistical inputs: "
            + ", ".join(missing)
        )
    section_region = section_region[
        section_region["region"].astype(str).str.upper().isin({"ARC", "ME"})
    ].copy()
    section_region["region"] = section_region["region"].astype(str).str.upper()
    if section_region.empty:
        raise ValueError("Registered section-region summary has no ARC/ME rows")
    identity_columns = ["animal_id", "sample", "section_index", "tile", "region"]
    normalized_identity = section_region[identity_columns].copy()
    normalized_identity["tile"] = normalized_identity["tile"].map(clean_tile)
    if normalized_identity.duplicated(identity_columns, keep=False).any():
        raise ValueError("Registered section-region summary has duplicate section/region rows")
    quantified_identities = {
        (str(row.animal_id), str(row.sample), int(row.section_index), clean_tile(row.tile))
        for row in quantified[["animal_id", "sample", "section_index", "tile"]].itertuples(index=False)
    }
    observed_identities = {
        (str(row.animal_id), str(row.sample), int(row.section_index), clean_tile(row.tile))
        for row in section_region[["animal_id", "sample", "section_index", "tile"]].itertuples(index=False)
    }
    unexpected = sorted(observed_identities.difference(quantified_identities))
    if unexpected:
        raise ValueError(
            "Section-region statistics contain rows outside the registered/HIL set: "
            + ", ".join(map(str, unexpected[:5]))
        )

    count_columns = [
        "dapi_nuclei", "cfos_cells", "pomc_cells", "cfos_pomc_all",
        "pomc_assay_dapi_nuclei",
    ]
    for column in count_columns:
        section_region[column] = pd.to_numeric(section_region[column], errors="coerce")
        values = section_region[column].to_numpy(float)
        if (
            not np.all(np.isfinite(values))
            or np.any(values < 0)
            or not np.allclose(values, np.rint(values), atol=1e-9, rtol=0)
        ):
            raise ValueError(f"Registered section-region {column} must be finite nonnegative counts")
    if np.any(section_region["cfos_cells"] > section_region["dapi_nuclei"]):
        raise ValueError("Registered regional c-FOS counts exceed DAPI nuclei")
    if np.any(section_region["cfos_pomc_all"] > section_region["pomc_cells"]):
        raise ValueError("Registered activated-POMC counts exceed total POMC counts")

    group_columns = ["animal_id", "animal_key", "region"]
    metadata_columns = ["cond", "genotype", "sex", "cage"]
    metadata = consistent_group_metadata(
        section_region,
        group_columns,
        metadata_columns,
        context="registered animal-region",
    )
    values = (
        section_region.groupby(group_columns, dropna=False, sort=False)[count_columns]
        .sum().reset_index().merge(metadata, on=group_columns, how="left", validate="one_to_one")
    )
    values["cfos_over_dapi"] = np.where(
        values["dapi_nuclei"] > 0,
        values["cfos_cells"] / values["dapi_nuclei"], np.nan,
    )
    values["cfos_pomc_over_pomc"] = np.where(
        values["pomc_cells"] > 0,
        values["cfos_pomc_all"] / values["pomc_cells"], np.nan,
    )

    # Audit the analyzer's persisted aggregation but never use its ratio columns
    # as the displayed/statistical values.
    persisted = pd.read_csv(per_animal_path, low_memory=False)
    needed_persisted = set(group_columns + metadata_columns + count_columns)
    missing_persisted = sorted(needed_persisted.difference(persisted.columns))
    if missing_persisted:
        raise ValueError(
            "Persisted animal-region summary lacks reaggregation audit fields: "
            + ", ".join(missing_persisted)
        )
    persisted = persisted[
        persisted["region"].astype(str).str.upper().isin({"ARC", "ME"})
    ].copy()
    persisted["region"] = persisted["region"].astype(str).str.upper()
    if persisted.duplicated(group_columns, keep=False).any():
        raise ValueError("Persisted animal-region summary has duplicate identities")
    comparison = values.merge(
        persisted[group_columns + metadata_columns + count_columns],
        on=group_columns, how="outer", suffixes=("_recomputed", "_persisted"),
        indicator=True, validate="one_to_one",
    )
    if not comparison["_merge"].eq("both").all():
        raise ValueError("Persisted animal-region identities differ from registered reaggregation")
    for column in metadata_columns:
        if not comparison[f"{column}_recomputed"].fillna("").astype(str).eq(
            comparison[f"{column}_persisted"].fillna("").astype(str)
        ).all():
            raise ValueError(f"Persisted {column} differs from registered reaggregation")
    for column in count_columns:
        left = pd.to_numeric(comparison[f"{column}_recomputed"], errors="coerce").to_numpy(float)
        right = pd.to_numeric(comparison[f"{column}_persisted"], errors="coerce").to_numpy(float)
        if not np.allclose(left, right, atol=1e-9, rtol=0, equal_nan=True):
            raise ValueError(
                f"Persisted {column} differs from current registered-section reaggregation"
            )

    # Negative-control immunostaining validates specificity but is not one of
    # the three biological treatment groups in Figure 4.
    values = values[values["cond"].astype(str).isin(CONDITION_ORDER)].copy()
    observed_conditions = set(values["cond"].astype(str))
    if observed_conditions != set(CONDITION_ORDER):
        raise ValueError(
            "Registered Figure 4 values do not contain exactly the Water/Sucrose/Allulose "
            f"treatment set; found {sorted(observed_conditions)}"
        )

    endpoint_rows: List[Dict[str, object]] = []
    endpoint_specs = {
        "total_cfos_activation_fraction": ("cfos_cells", "dapi_nuclei", "cfos_over_dapi"),
        "primary_pomc_activation_fraction": (
            "cfos_pomc_all", "pomc_cells", "cfos_pomc_over_pomc",
        ),
    }
    for row in values.to_dict("records"):
        for endpoint, (numerator_column, denominator_column, value_column) in endpoint_specs.items():
            numerator = float(row[numerator_column])
            denominator = float(row[denominator_column])
            raw_value = float(row[value_column])
            valid = bool(denominator > 0 and np.isfinite(raw_value))
            endpoint_rows.append({
                **{key: row[key] for key in [
                    "animal_id", "animal_key", "region", "cond", "genotype", "sex", "cage",
                ]},
                "endpoint": endpoint,
                "raw_value": raw_value,
                "numerator": numerator,
                "denominator": denominator,
                "valid_for_inference": valid,
                "exclusion_reason": "" if valid else f"nonpositive_{denominator_column}",
            })
    definitive_values = pd.DataFrame(endpoint_rows)
    values = values.sort_values(group_columns, kind="mergesort").reset_index(drop=True)
    definitive_values = definitive_values.sort_values(
        ["endpoint", "region", "cond", "cage", "animal_id"], kind="mergesort",
    ).reset_index(drop=True)
    return values, definitive_values, section_region_path


def _accepted_review_paths(
    manual_root: Path, row: Mapping[str, object], folder: Path,
) -> Tuple[Path, Path]:
    animal = str(row["animal_id"])
    section = folder.name
    directory = manual_root / animal / section
    stem = f"{animal}_{section}_ARC_ME_labels"
    return directory / f"{stem}.json", directory / f"{stem}.tif"


def _read_label_mask(path: Path, shape: Tuple[int, int]) -> np.ndarray:
    raw = np.squeeze(np.asarray(tifffile.imread(str(path))))
    if raw.ndim != 2 or tuple(raw.shape) != tuple(shape):
        raise ValueError(f"ARC/ME mask shape {raw.shape} != microscopy shape {shape}: {path}")
    if not np.issubdtype(raw.dtype, np.number) or not np.all(np.isfinite(raw)):
        raise ValueError(f"ARC/ME mask must be finite numeric labels: {path}")
    if not np.all(raw == np.rint(raw)):
        raise ValueError(f"ARC/ME mask contains non-integer labels: {path}")
    labels = np.asarray(raw, dtype=np.uint8)
    codes = {int(value) for value in np.unique(labels)}
    if not codes.issubset({0, 1, 2, 3}) or not codes.intersection({1, 2, 3}):
        raise ValueError(f"ARC/ME/VMN mask needs any positive ROI with codes 0/1/2/3; got {sorted(codes)}: {path}")
    return labels


def validate_hil_analysis_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Refuse any plotted run that is not wholly HIL accepted."""
    required = {
        "status", "sample", "animal_id", "section_index", "tile", "region_source",
        "manual_region_review_status", "manual_region_reviewer",
        "manual_region_accepted_at", "manual_region_accepted_source",
        "manual_region_path", "manual_region_review_path",
        "manual_region_sha256", "manual_region_review_sha256",
        "manual_region_source_dapi_path", "manual_region_source_dapi_sha256",
        "manual_region_outside_policy",
        "dapi_nuclei", "dapi_nuclei_whole_image", "dapi_nuclei_code0_excluded",
    }
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError(
            "Final POMC run predates the mandatory all-HIL/code-0 policy; "
            f"missing columns: {missing}"
        )
    quantified = summary[summary["status"].astype(str).eq("ok")].copy()
    identity_columns = ["animal_id", "sample", "section_index", "tile"]
    if summary.duplicated(identity_columns, keep=False).any():
        raise ValueError("POMC section summary contains duplicate repetition-safe identities")
    if len(quantified) != EXPECTED_POMC_QUANTIFIABLE_SECTIONS:
        raise ValueError(
            "Refusing incomplete POMC figure input: expected exactly "
            f"{EXPECTED_POMC_QUANTIFIABLE_SECTIONS} quantifiable repetition-safe "
            f"sections, found {len(quantified)}"
        )
    if quantified.duplicated(identity_columns, keep=False).any():
        raise ValueError("Quantifiable POMC set contains duplicate repetition-safe identities")
    if quantified.empty:
        raise ValueError("Final POMC run has no quantifiable status=ok sections")
    problems: List[str] = []
    allowed_sources = {"automated_mask_unchanged", "manual_replacement_polygons"}
    for row in quantified.itertuples(index=False):
        identity = f"{getattr(row, 'animal_id')}/{getattr(row, 'sample')}_S{int(getattr(row, 'section_index')):02d}"
        tile = clean_tile(getattr(row, "tile"))
        expected_section = (
            f"{getattr(row, 'sample')}_S{int(getattr(row, 'section_index')):02d}"
            f"{('_T' + tile) if tile else ''}"
        )
        if str(getattr(row, "region_source")) != "manual_region_tiff":
            problems.append(f"{identity}: region_source is not manual_region_tiff")
            continue
        if str(getattr(row, "manual_region_review_status")) != "accepted":
            problems.append(f"{identity}: HIL review status is not accepted")
            continue
        if not str(getattr(row, "manual_region_reviewer")).strip():
            problems.append(f"{identity}: reviewer is blank")
            continue
        if not _timezone_aware(getattr(row, "manual_region_accepted_at")):
            problems.append(f"{identity}: accepted_at is not timezone-aware")
            continue
        if str(getattr(row, "manual_region_accepted_source")) not in allowed_sources:
            problems.append(f"{identity}: accepted source is not native POMC HIL")
            continue
        outside_policy = str(getattr(row, "manual_region_outside_policy"))
        if "code0 excluded from all inferential" not in outside_policy:
            problems.append(f"{identity}: code-0 exclusion policy is absent")
            continue
        total = float(getattr(row, "dapi_nuclei"))
        whole = float(getattr(row, "dapi_nuclei_whole_image"))
        excluded = float(getattr(row, "dapi_nuclei_code0_excluded"))
        if total < 0 or whole < total or excluded < 0 or not np.isclose(total + excluded, whole):
            problems.append(f"{identity}: ARC+ME and code-0 DAPI totals are inconsistent")
            continue
        mask_path = Path(str(getattr(row, "manual_region_path"))).expanduser().resolve()
        receipt_path = Path(str(getattr(row, "manual_region_review_path"))).expanduser().resolve()
        if not mask_path.is_file() or not receipt_path.is_file():
            problems.append(f"{identity}: accepted TIFF/JSON pair is missing")
            continue
        mask_sha = str(getattr(row, "manual_region_sha256")).lower().strip()
        receipt_sha = str(getattr(row, "manual_region_review_sha256")).lower().strip()
        if not _valid_sha256(mask_sha) or sha256_file(mask_path) != mask_sha:
            problems.append(f"{identity}: accepted mask hash no longer matches")
            continue
        if not _valid_sha256(receipt_sha) or sha256_file(receipt_path) != receipt_sha:
            problems.append(f"{identity}: review receipt hash no longer matches")
            continue
        source_dapi_path = Path(
            str(getattr(row, "manual_region_source_dapi_path"))
        ).expanduser().resolve()
        source_dapi_sha = str(
            getattr(row, "manual_region_source_dapi_sha256")
        ).lower().strip()
        if not source_dapi_path.is_file():
            problems.append(f"{identity}: reconstructed source DAPI is missing")
            continue
        if not _valid_sha256(source_dapi_sha) or sha256_file(source_dapi_path) != source_dapi_sha:
            problems.append(f"{identity}: reconstructed source DAPI hash no longer matches")
            continue
        try:
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        except Exception as exc:
            problems.append(f"{identity}: review receipt is not valid JSON ({exc})")
            continue
        source_paths = receipt.get("source_paths") if isinstance(receipt, dict) else None
        if not isinstance(source_paths, dict):
            problems.append(f"{identity}: receipt source_paths is invalid")
            continue
        receipt_dapi_path = _resolve_provenance_path(
            source_paths.get("dapi", ""),
            folder=source_dapi_path.parent,
            root=PAPER_ROOT,
        )
        # Named individually so a failure says which binding broke rather than
        # only that "provenance differs".
        receipt_checks = {
            "receipt is not marked accepted": receipt.get("accepted") is True,
            "animal identity": str(receipt.get("animal", "")) == str(getattr(row, "animal_id")),
            "source sample identity": (
                str(receipt.get("source_sample", "")) == str(getattr(row, "sample"))
            ),
            "section identity": str(receipt.get("section", "")) == expected_section,
            "analyzed section folder name": source_dapi_path.parent.name == expected_section,
            "dataset": receipt.get("dataset") == "pomc",
            "reviewer": (
                str(receipt.get("reviewer", "")).strip()
                == str(getattr(row, "manual_region_reviewer")).strip()
            ),
            "acceptance timestamp": (
                str(receipt.get("accepted_at", ""))
                == str(getattr(row, "manual_region_accepted_at"))
            ),
            "accepted source": (
                str(receipt.get("accepted_source", ""))
                == str(getattr(row, "manual_region_accepted_source"))
            ),
            "accepted mask SHA-256": str(receipt.get("mask_sha256", "")).lower() == mask_sha,
            "reviewed section": _same_reconstructed_section(
                receipt_dapi_path, source_dapi_path,
            ),
            "reviewed source DAPI SHA-256": (
                str(receipt.get("source_dapi_sha256", "")).lower() == source_dapi_sha
            ),
        }
        failed = [name for name, ok in receipt_checks.items() if not ok]
        if failed:
            problems.append(
                f"{identity}: receipt provenance differs from the analyzed row "
                f"({', '.join(failed)})"
            )
            continue
    if problems:
        detail = "\n".join(f"  - {item}" for item in problems[:8])
        if len(problems) > 8:
            detail += f"\n  - ... and {len(problems) - 8} more"
        raise ValueError(
            f"Refusing stale/mixed POMC figure input: {len(problems)} of "
            f"{len(quantified)} quantifiable sections fail mandatory HIL validation.\n{detail}"
        )
    return quantified


def validate_analysis_summary(
    summary: pd.DataFrame, *, require_manual_regions: bool,
) -> pd.DataFrame:
    """Validate automated/optional-HIL runs, or delegate to strict HIL checks."""
    if require_manual_regions:
        return validate_hil_analysis_summary(summary)
    required = {
        "status", "sample", "animal_id", "section_index", "tile",
        "dapi_nuclei", "arc_area_px", "me_area_px",
    }
    missing = sorted(required.difference(summary.columns))
    if missing:
        raise ValueError(
            "POMC analysis summary lacks required automated/HIL provenance fields: "
            + ", ".join(missing)
        )
    identity_columns = ["animal_id", "sample", "section_index", "tile"]
    if summary.duplicated(identity_columns, keep=False).any():
        raise ValueError("POMC section summary contains duplicate repetition-safe identities")
    quantified = summary[summary["status"].astype(str).eq("ok")].copy()
    if quantified.empty:
        raise ValueError("POMC analysis has no quantifiable status=ok sections")
    if len(quantified) != EXPECTED_POMC_QUANTIFIABLE_SECTIONS:
        raise ValueError(
            "Refusing incomplete POMC figure input: expected exactly "
            f"{EXPECTED_POMC_QUANTIFIABLE_SECTIONS} quantifiable sections, "
            f"found {len(quantified)}"
        )
    if "region_source" in quantified.columns:
        sources = quantified["region_source"].fillna("").astype(str).str.lower()
        forbidden = sources.str.contains("transfer|register|warp", regex=True)
        if forbidden.any():
            found = sorted(set(sources[forbidden]))
            raise ValueError(
                "Transferred/registered ARC/ME provenance is forbidden: "
                + ", ".join(found)
            )
        manual = sources.eq("manual_region_tiff")
        if manual.any() and not manual.all():
            raise ValueError(
                "Mixed manual/automated POMC ARC/ME inference is forbidden; "
                "use complete 76/76 review or all native automated masks."
            )
        if manual.all():
            return validate_hil_analysis_summary(summary)
        # Old v17 summaries predate region_source and therefore contain a
        # missing/blank value. New runs state native_automated_v17 explicitly.
        # Reject every other nonblank source even when its model label happens
        # to resemble v17, so transferred/unknown provenance cannot slip into
        # a publication render.
        allowed_automated_sources = {"", "native_automated_v17"}
        unexpected_sources = sorted(
            value for value in set(sources) if value not in allowed_automated_sources
        )
        if unexpected_sources:
            raise ValueError(
                "Automated POMC figure input has unsupported region_source: "
                + ", ".join(unexpected_sources)
            )
    expected_model = "oval_aware_projected_walls_floor_only_ARC_ME_v17"
    model = quantified["region_ai_model"].fillna("").astype(str) if "region_ai_model" in quantified else pd.Series([], dtype=str)
    mode = quantified["region_mode"].fillna("").astype(str) if "region_mode" in quantified else pd.Series([], dtype=str)
    if len(model) != len(quantified) or not model.eq(expected_model).all():
        found = sorted(set(model)) if len(model) else ["missing"]
        raise ValueError(
            "Automated POMC figure input must use the native v17 projected-wall "
            f"ARC/ME model; found {found}"
        )
    if len(mode) != len(quantified) or not mode.eq("arc-me-walls").all():
        found = sorted(set(mode)) if len(mode) else ["missing"]
        raise ValueError(
            "Automated POMC figure input must use region_mode=arc-me-walls; "
            f"found {found}"
        )
    return quantified


def validate_native_automated_masks(
    quantified: pd.DataFrame, root: Path,
) -> None:
    """Verify all 76 saved native masks against analyzed ARC/ME pixel areas."""
    problems: List[str] = []
    for row in quantified.to_dict("records"):
        folder = section_directory(root, row)
        sources = microscopy_sources(folder)
        identity = f"{row['animal_id']}/{folder.name}"
        try:
            dapi = read_gray(sources["dapi"])
            labels = _read_label_mask(sources["native_proposal"], dapi.shape)
            arc_px = int(np.count_nonzero(labels == 1))
            me_px = int(np.count_nonzero(labels == 2))
            expected_arc = int(round(float(row["arc_area_px"])))
            expected_me = int(round(float(row["me_area_px"])))
            if (arc_px, me_px) != (expected_arc, expected_me):
                raise ValueError(
                    f"mask areas {(arc_px, me_px)} != summary "
                    f"{(expected_arc, expected_me)}"
                )
        except Exception as exc:
            problems.append(f"{identity}: {exc}")
    if problems:
        detail = "\n".join(f"  - {item}" for item in problems[:8])
        if len(problems) > 8:
            detail += f"\n  - ... and {len(problems) - 8} more"
        raise ValueError(
            f"Native automated POMC mask audit failed for {len(problems)} "
            f"sections.\n{detail}"
        )


def load_accepted_pomc_review(
    row: Mapping[str, object], folder: Path, manual_root: Path,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """Validate the representative receipt and return its accepted native mask."""
    sources = microscopy_sources(folder)
    dapi = read_gray(sources["dapi"])
    review_path, mask_path = _accepted_review_paths(manual_root, row, folder)
    if not review_path.is_file() or not mask_path.is_file():
        raise FileNotFoundError(
            f"Representative HIL-accepted TIFF+JSON pair is missing: {review_path}, {mask_path}"
        )
    review = json.loads(review_path.read_text(encoding="utf-8"))
    if not isinstance(review, dict) or review.get("accepted") is not True:
        raise ValueError("Representative receipt must be one object with accepted=true")
    if review.get("dataset") != "pomc":
        raise ValueError("Representative receipt dataset must equal pomc")
    if str(review.get("animal", "")) != str(row["animal_id"]):
        raise ValueError("Representative receipt animal does not match analysis")
    if str(review.get("source_sample", "")) != str(row["sample"]):
        raise ValueError("Representative receipt source_sample does not match analysis")
    if str(review.get("section", "")) != folder.name:
        raise ValueError("Representative receipt section does not match reconstructed field")
    reviewer = str(review.get("reviewer", "")).strip()
    accepted_at = str(review.get("accepted_at", "")).strip()
    if not reviewer or not _timezone_aware(accepted_at):
        raise ValueError("Representative receipt needs reviewer and timezone-aware accepted_at")
    accepted_source = str(review.get("accepted_source", ""))
    if accepted_source not in {"automated_mask_unchanged", "manual_replacement_polygons"}:
        raise ValueError(f"Non-native/transfer accepted source is forbidden: {accepted_source!r}")
    for forbidden in (
        "transferred_proposal_sha256", "source_manifest_sha256",
        "proposal_confirmation", "proposal_mask_sha256",
    ):
        if str(review.get(forbidden, "")).strip():
            raise ValueError(f"Transfer/legacy proposal field is forbidden: {forbidden}")
    receipt_codes = review.get("label_codes")
    if receipt_codes not in ({"background": 0, "ARC": 1, "ME": 2},
                              {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}):
        raise ValueError("Representative receipt has the wrong ARC/ME/VMN label-code contract")
    if (int(review.get("image_height_px", -1)), int(review.get("image_width_px", -1))) != dapi.shape:
        raise ValueError("Representative receipt dimensions do not match native DAPI")
    source_paths = review.get("source_paths")
    if not isinstance(source_paths, dict):
        raise ValueError("Representative receipt source_paths must be an object")
    source_dapi_path = _resolve_provenance_path(
        source_paths.get("dapi", ""), folder=folder, root=PAPER_ROOT,
    )
    if not _same_reconstructed_section(source_dapi_path, sources["dapi"].resolve()):
        raise ValueError("Representative receipt points to a different native DAPI field")
    source_dapi_sha = sha256_file(sources["dapi"])
    if str(review.get("source_dapi_sha256", "")).lower() != source_dapi_sha:
        raise ValueError("Representative native DAPI TIFF-container SHA-256 mismatch")
    mask = _read_label_mask(mask_path, dapi.shape)
    mask_sha = sha256_file(mask_path)
    if str(review.get("mask_sha256", "")).lower() != mask_sha:
        raise ValueError("Representative accepted mask SHA-256 mismatch")
    if receipt_codes == {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}:
        raw_status = review.get("region_status")
        if not isinstance(raw_status, dict):
            raise ValueError("Representative schema-v3 receipt lacks region_status")
        region_status = {
            name: str(raw_status.get(name, ""))
            for name in ("ARC", "ME", "VMN")
        }
        if any(value not in {"drawn", "explicitly_absent", "not_delineated"} for value in region_status.values()):
            raise ValueError("Representative schema-v3 receipt has invalid region_status")
        reasons = review.get("region_absence_reasons")
        if not isinstance(reasons, dict):
            raise ValueError("Representative schema-v3 receipt lacks absence-reason map")
    else:
        me_absent_legacy = review.get("me_absent") is True
        region_status = {
            "ARC": "drawn" if np.any(mask == 1) else "not_delineated",
            "ME": "drawn" if np.any(mask == 2) else "explicitly_absent" if me_absent_legacy else "not_delineated",
            "VMN": "not_delineated",
        }
        reasons = {"ARC": "", "ME": str(review.get("me_absent_reason", "")).strip(), "VMN": ""}
    for name, code in (("ARC", 1), ("ME", 2), ("VMN", 3)):
        area = int(np.count_nonzero(mask == code))
        state = region_status[name]
        reason = str(reasons.get(name, "")).strip()
        if (state == "drawn") != (area > 0):
            raise ValueError(f"Representative {name} status conflicts with mask area")
        if state == "explicitly_absent" and not reason:
            raise ValueError(f"Representative explicit {name} absence lacks a reason")
        if state != "explicitly_absent" and reason:
            raise ValueError(f"Representative {name} absence reason conflicts with {state}")
    if not any(value == "drawn" for value in region_status.values()):
        raise ValueError("Representative mask has no delineated ARC, ME, or VMN ROI")
    me_absent = region_status["ME"] == "explicitly_absent"
    me_reason = str(reasons.get("ME", "")).strip()

    native_sha = ""
    native_current_sha = ""
    native_current_match = False
    native_validation = "not_applicable"
    receipt_native_path_text = ""
    if accepted_source == "automated_mask_unchanged":
        if review.get("automated_mask_confirmation") != "reviewed_image_by_image":
            raise ValueError("Native automated acceptance lacks image-by-image confirmation")
        receipt_native_path_text = str(source_paths.get("automated_mask", "")).strip()
        native_sha = str(review.get("automated_mask_sha256", "")).strip().lower()
        if not receipt_native_path_text or not _valid_sha256(native_sha):
            raise ValueError("Representative historical native-proposal provenance is incomplete")
        provenance = review.get("automated_mask_provenance")
        if (
            not isinstance(provenance, dict)
            or provenance.get("dataset") != "pomc"
            or provenance.get("kind") != "native_pomc_automated_mask"
        ):
            raise ValueError("Representative native proposal provenance is invalid")

        # The accepted HIL TIFF is authoritative after the image-by-image
        # decision.  A later preparation/full-analysis pass may rebuild the
        # automated starting proposal, so current-proposal disagreement is
        # recorded rather than used to discard a valid accepted mask.
        native_validation = "accepted_tiff_receipt"
        current_path = sources["native_proposal"]
        if current_path.is_file():
            native_current_sha = sha256_file(current_path)
            try:
                historical_path = _resolve_provenance_path(
                    receipt_native_path_text, folder=folder, root=PAPER_ROOT,
                )
                current_mask = _read_label_mask(current_path, dapi.shape)
                native_current_match = bool(
                    historical_path == current_path.resolve()
                    and native_current_sha == native_sha
                    and np.array_equal(mask, current_mask)
                )
            except Exception:
                native_current_match = False
        if native_current_match:
            native_validation = "exact_current_proposal"

    accepted_output = sources["accepted_regions"]
    if not accepted_output.is_file():
        raise FileNotFoundError(
            f"Analyzed final accepted mask is missing: {accepted_output}"
        )
    analyzed_mask = _read_label_mask(accepted_output, dapi.shape)
    if not np.array_equal(mask, analyzed_mask):
        raise ValueError("Analyzed final mask differs from the accepted HIL mask")

    row_checks = {
        "manual_region_path": mask_path.resolve(),
        "manual_region_review_path": review_path.resolve(),
    }
    for column, expected in row_checks.items():
        if Path(str(row[column])).expanduser().resolve() != expected:
            raise ValueError(f"Analysis row {column} does not match representative receipt")
    if str(row["manual_region_sha256"]).lower() != mask_sha:
        raise ValueError("Analysis row mask hash does not match representative receipt")
    if str(row["manual_region_review_sha256"]).lower() != sha256_file(review_path):
        raise ValueError("Analysis row receipt hash does not match representative receipt")
    if str(row["manual_region_source_dapi_sha256"]).lower() != source_dapi_sha:
        raise ValueError("Analysis row source DAPI hash does not match representative field")
    if str(row["manual_region_reviewer"]).strip() != reviewer:
        raise ValueError("Analysis row reviewer does not match representative receipt")
    if str(row["manual_region_accepted_source"]) != accepted_source:
        raise ValueError("Analysis row accepted source does not match representative receipt")

    provenance: Dict[str, object] = {
        "human_reviewed": True,
        "review_status": "accepted",
        "reviewer": reviewer,
        "accepted_at": accepted_at,
        "accepted_source": accepted_source,
        "mask_path": str(mask_path.resolve()),
        "mask_sha256": mask_sha,
        "review_path": str(review_path.resolve()),
        "review_sha256": sha256_file(review_path),
        "source_dapi_path": str(sources["dapi"].resolve()),
        "source_dapi_sha256": source_dapi_sha,
        "native_proposal_path": receipt_native_path_text,
        "native_proposal_sha256": native_sha,
        "current_native_proposal_path": str(sources["native_proposal"].resolve()),
        "current_native_proposal_sha256": native_current_sha,
        "current_native_proposal_matches_historical_review": native_current_match,
        "native_proposal_validation": native_validation,
        "analyzed_accepted_mask_path": str(accepted_output.resolve()),
        "analyzed_accepted_mask_sha256": sha256_file(accepted_output),
        "me_absent": me_absent,
        "me_absent_reason": me_reason,
        "region_status": region_status,
    }
    return mask, provenance


def load_current_pomc_regions(
    row: Mapping[str, object], folder: Path, manual_root: Path,
    *, require_manual_regions: bool,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """Load a valid accepted mask, else the section's native automated mask."""
    row_source = str(row.get("region_source", ""))
    row_review = str(row.get("manual_region_review_status", ""))
    if require_manual_regions or (
        row_source == "manual_region_tiff" and row_review == "accepted"
    ):
        return load_accepted_pomc_review(row, folder, manual_root)

    sources = microscopy_sources(folder)
    dapi = read_gray(sources["dapi"])
    mask_path = sources["native_proposal"]
    if not mask_path.is_file():
        raise FileNotFoundError(
            f"Current native automated ARC/ME mask is missing: {mask_path}"
        )
    mask = _read_label_mask(mask_path, dapi.shape)
    mask_sha = sha256_file(mask_path)
    dapi_sha = sha256_file(sources["dapi"])
    has_me = bool(np.any(mask == 2))
    provenance: Dict[str, object] = {
        "human_reviewed": False,
        "review_status": "not_requested_optional",
        "reviewer": "not performed",
        "accepted_at": "not applicable",
        "accepted_source": "native_automated_mask_unreviewed",
        "mask_path": str(mask_path.resolve()),
        "mask_sha256": mask_sha,
        "review_path": "",
        "review_sha256": "",
        "source_dapi_path": str(sources["dapi"].resolve()),
        "source_dapi_sha256": dapi_sha,
        "native_proposal_path": str(mask_path.resolve()),
        "native_proposal_sha256": mask_sha,
        "analyzed_accepted_mask_path": "",
        "analyzed_accepted_mask_sha256": "",
        "me_absent": not has_me,
        "me_absent_reason": "native automated mask contains no ME" if not has_me else "",
    }
    return mask, provenance

def select_representative(summary: pd.DataFrame, root: Path, sample: str, section: int) -> Tuple[pd.Series, Path, pd.DataFrame]:
    required = {
        "sample", "animal_id", "section_index", "status", "cond",
        "has_pomc_channel", "pomc_channel_complete", "has_left", "has_right",
        "cfos_over_dapi", "dapi_nuclei", "pomc_cells", "arc_area_px", "me_area_px",
    }
    missing = sorted(required - set(summary.columns))
    if missing:
        raise ValueError(f"Missing columns in per_reconstructed_section_summary.csv: {missing}")
    arc_area = pd.to_numeric(summary["arc_area_px"], errors="coerce").fillna(0)
    me_area = pd.to_numeric(summary["me_area_px"], errors="coerce").fillna(0)
    vmn_area = (
        pd.to_numeric(summary["vmn_area_px"], errors="coerce").fillna(0)
        if "vmn_area_px" in summary.columns
        else pd.Series(0.0, index=summary.index)
    )
    # A schema-v3 review may legitimately delineate any nonempty subset of
    # ARC/ME/VMN. For this ARC/ME publication panel, retain ARC-only or ME-only
    # reviewed sections and prefer (but never fabricate) a section containing
    # both ARC and ME.
    q = summary[
        summary["status"].astype(str).eq("ok")
        & (arc_area.gt(0) | me_area.gt(0) | vmn_area.gt(0))
        & as_bool(summary["has_pomc_channel"])
        & as_bool(summary["pomc_channel_complete"])
        & as_bool(summary["has_left"])
        & as_bool(summary["has_right"])
        & ~summary["cond"].astype(str).str.lower().eq("control")
    ].copy()
    # Prefer an automatically verified overlap, then a low-confidence overlap
    # explicitly confirmed during image-by-image HIL, and only then a native
    # single half. Legacy selected_corr/y-only registration is never consulted.
    if "stitch_alignment_status" not in q.columns:
        raise ValueError("Representative selection requires corrected registration status")
    q["registration_status_priority"] = q["stitch_alignment_status"].map({
        "verified_overlap": 0,
        "manual_review_required": 1,
        "single_half_native": 2,
    }).fillna(9).astype(int)
    ncc_best = pd.to_numeric(
        q.get("stitch_alignment_ncc_best_score"), errors="coerce",
    )
    ncc_raw = pd.to_numeric(
        q.get("stitch_alignment_raw_overlap_ncc"), errors="coerce",
    )
    q["registration_ncc_priority"] = -ncc_best.fillna(ncc_raw).fillna(-np.inf)
    if sample:
        q = q[q["sample"].astype(str).eq(sample) | q["animal_id"].astype(str).eq(sample)].copy()
        if section > 0:
            q = q[pd.to_numeric(q["section_index"], errors="coerce").eq(section)]
        if q.empty:
            raise ValueError(f"No complete representative candidate for {sample!r}, section {section or 'any'}")
    q_arc = pd.to_numeric(q["arc_area_px"], errors="coerce").fillna(0)
    q_me = pd.to_numeric(q["me_area_px"], errors="coerce").fillna(0)
    q_vmn = (
        pd.to_numeric(q["vmn_area_px"], errors="coerce").fillna(0)
        if "vmn_area_px" in q.columns else pd.Series(0.0, index=q.index)
    )
    q["arc_me_pair_priority"] = np.select(
        [q_arc.gt(0) & q_me.gt(0), q_arc.gt(0) | q_me.gt(0), q_vmn.gt(0)],
        [0, 1, 2],
        default=3,
    )
    q["me_fraction"] = q_me / (q_me + q_arc).replace(0, np.nan)
    metrics = ["cfos_over_dapi", "dapi_nuclei", "pomc_cells", "me_fraction"]
    q["representativeness_score"] = 0.0
    for metric in metrics:
        values = pd.to_numeric(q[metric], errors="coerce")
        median = float(values.median())
        mad = float((values - median).abs().median())
        if not np.isfinite(mad) or mad <= 0:
            mad = float(values.std(ddof=1))
        if not np.isfinite(mad) or mad <= 0:
            mad = 1.0
        q[f"selection_median_{metric}"] = median
        q[f"selection_mad_{metric}"] = mad
        q["representativeness_score"] += (values - median).abs().fillna(5.0 * mad) / mad
    valid: List[int] = []
    folders: Dict[int, Path] = {}
    for idx, row in q.iterrows():
        folder = section_directory(root, row)
        src = microscopy_sources(folder)
        if all(src[k].exists() for k in ["dapi", "cfos", "pomc", "dapi_labels"]):
            valid.append(idx); folders[idx] = folder
    q = q.loc[valid].sort_values(
        [
            "registration_status_priority", "arc_me_pair_priority",
            "representativeness_score", "registration_ncc_priority",
            "animal_id", "sample", "section_index",
        ],
        kind="mergesort",
    )
    if q.empty:
        raise FileNotFoundError("No representative candidate has all reconstructed TIFFs.")
    chosen = q.iloc[0]
    return chosen, folders[int(chosen.name)], q


def read_gray(path: Path) -> np.ndarray:
    arr = np.asarray(tifffile.imread(str(path)))
    while arr.ndim > 3:
        arr = arr.max(axis=0)
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        arr = arr[..., :3].astype(np.float32).max(axis=-1)
    elif arr.ndim == 3:
        arr = arr.max(axis=0)
    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2-D microscopy plane: {path} -> {arr.shape}")
    return arr



def load_raw_cellpose_labels(path: Path) -> np.ndarray:
    """Load one local Cellpose ``*_seg.npy`` label plane without using its preview."""
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    payload = np.load(str(source), allow_pickle=True)
    if isinstance(payload, np.ndarray) and payload.dtype == object and payload.shape == ():
        payload = payload.item()
    if not isinstance(payload, Mapping) or "masks" not in payload:
        raise ValueError(f"Cellpose raw NPY lacks a masks array: {source}")
    labels = np.squeeze(np.asarray(payload["masks"]))
    if labels.ndim != 2:
        raise ValueError(f"Cellpose masks must be two-dimensional: {source} -> {labels.shape}")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError(f"Cellpose masks must contain integer instance labels: {source}")
    if np.any(labels < 0):
        raise ValueError(f"Cellpose masks contain negative labels: {source}")
    return labels.astype(np.int32, copy=False)


def _integer_overlap_views(
    left: np.ndarray, right: np.ndarray, dx: int, dy: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the two label views that occupy identical registered coordinates."""
    x0 = max(0, int(dx)); x1 = min(int(left.shape[1]), int(dx) + int(right.shape[1]))
    y0 = max(0, int(dy)); y1 = min(int(left.shape[0]), int(dy) + int(right.shape[0]))
    if x1 <= x0 or y1 <= y0:
        return (
            np.empty((0, 0), dtype=left.dtype),
            np.empty((0, 0), dtype=right.dtype),
        )
    return (
        left[y0:y1, x0:x1],
        right[y0 - int(dy):y1 - int(dy), x0 - int(dx):x1 - int(dx)],
    )


def _place_vertical_pair(
    left: np.ndarray, right: np.ndarray, shift_right_px: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Pad, never resample or crop, while vertically registering two masks."""
    shift = int(shift_right_px)
    left_top = max(0, -shift); right_top = max(0, shift)
    height = max(left_top + left.shape[0], right_top + right.shape[0])
    left_out = np.zeros((height, left.shape[1]), dtype=left.dtype)
    right_out = np.zeros((height, right.shape[1]), dtype=right.dtype)
    left_out[left_top:left_top + left.shape[0], :] = left
    right_out[right_top:right_top + right.shape[0], :] = right
    return left_out, right_out


def _reconcile_overlapping_cellpose_ids(
    left: np.ndarray, right: np.ndarray, dx: int, dy: int,
) -> np.ndarray:
    """Reproduce the analyzer's one-to-one ID mapping in the shared DAPI field."""
    left_overlap, right_overlap = _integer_overlap_views(left, right, dx, dy)
    matched: Dict[int, int] = {}
    if left_overlap.size and right_overlap.size:
        ll = left_overlap.astype(np.int64, copy=False).reshape(-1)
        rr = right_overlap.astype(np.int64, copy=False).reshape(-1)
        both = (ll > 0) & (rr > 0)
        if np.any(both):
            max_right = int(max(1, right.max()))
            encoded = ll[both] * (max_right + 1) + rr[both]
            keys, intersections = np.unique(encoded, return_counts=True)
            left_ids, left_areas = np.unique(ll[ll > 0], return_counts=True)
            right_ids, right_areas = np.unique(rr[rr > 0], return_counts=True)
            left_area = {int(key): int(value) for key, value in zip(left_ids, left_areas)}
            right_area = {int(key): int(value) for key, value in zip(right_ids, right_areas)}
            candidates: List[Tuple[float, float, int, int, int]] = []
            for key, intersection in zip(keys, intersections):
                left_id = int(key // (max_right + 1)); right_id = int(key % (max_right + 1))
                inter = int(intersection)
                union = left_area.get(left_id, 0) + right_area.get(right_id, 0) - inter
                iou = float(inter / union) if union > 0 else 0.0
                fraction = float(
                    inter / max(1, min(left_area.get(left_id, 0), right_area.get(right_id, 0)))
                )
                if inter >= 3 and (iou >= 0.20 or fraction >= 0.50):
                    candidates.append((iou, fraction, inter, left_id, right_id))
            used_left = set(); used_right = set()
            for _iou, _fraction, _intersection, left_id, right_id in sorted(
                candidates, reverse=True,
            ):
                if left_id in used_left or right_id in used_right:
                    continue
                matched[right_id] = left_id
                used_left.add(left_id); used_right.add(right_id)

    max_right = int(right.max()) if right.size else 0
    lut = np.zeros(max_right + 1, dtype=np.int32)
    next_id = int(left.max()) + 1 if left.size else 1
    for right_id in range(1, max_right + 1):
        if right_id in matched:
            lut[right_id] = matched[right_id]
        else:
            lut[right_id] = next_id
            next_id += 1
    return lut[right] if max_right > 0 else right.astype(np.int32, copy=True)


def stitch_raw_cellpose_labels(
    left: Optional[np.ndarray], right: Optional[np.ndarray],
    transform: Mapping[str, object],
) -> np.ndarray:
    """Apply the audited hard seam to raw Cellpose labels without resampling."""
    if left is None and right is None:
        raise ValueError("At least one raw Cellpose half is required")
    if left is None:
        return np.asarray(right, dtype=np.int32)
    if right is None:
        return np.asarray(left, dtype=np.int32)
    left = np.asarray(left, dtype=np.int32); right = np.asarray(right, dtype=np.int32)
    if left.ndim != 2 or right.ndim != 2:
        raise ValueError("Raw Cellpose halves must be two-dimensional")
    if str(transform.get("model", "")) != "integer_translation_hard_seam":
        raise ValueError("Figure 4 requires the audited integer hard-seam registration")
    if _strict_bool(transform.get("resampling"), field="transform.resampling"):
        raise ValueError("Figure 4 refuses a resampled Cellpose registration")
    orientation = str(transform.get("orientation", ""))
    if orientation != "left-right":
        raise ValueError(f"Figure 4 requires left-right registration, got {orientation!r}")
    dx = _finite_int(transform.get("right_origin_x_px"), field="transform.right_origin_x_px")
    dy = _finite_int(transform.get("right_shift_y_px"), field="transform.right_shift_y_px")
    seam = _finite_int(transform.get("left_retained_width_px"), field="transform.left_retained_width_px")
    right_crop = _finite_int(transform.get("right_medial_crop_px"), field="transform.right_medial_crop_px")
    overlap = _finite_int(transform.get("overlap_width_px"), field="transform.overlap_width_px")
    gap = _finite_int(transform.get("gap_px"), field="transform.gap_px")
    if not (0 < seam <= left.shape[1]) or not (0 <= right_crop < right.shape[1]):
        raise ValueError("Figure 4 Cellpose seam/crop falls outside a raw mask")
    if min(dx, overlap, gap) < 0:
        raise ValueError("Figure 4 Cellpose transform has a negative horizontal parameter")
    if overlap > 0:
        right = _reconcile_overlapping_cellpose_ids(left, right, dx, dy)
    else:
        right = right.copy()
        positive = right > 0
        right[positive] += int(left.max())
    left = left[:, :seam]
    right = right[:, right_crop:]
    left, right = _place_vertical_pair(left, right, dy)
    pieces = [left]
    if gap:
        pieces.append(np.zeros((left.shape[0], gap), dtype=np.int32))
    pieces.append(right)
    return np.concatenate(pieces, axis=1).astype(np.int32, copy=False)


def _filter_labels_by_area(labels: np.ndarray, min_area: int) -> np.ndarray:
    out = np.zeros(labels.shape, dtype=np.int32)
    current = 0
    for cell in regionprops(np.asarray(labels, dtype=np.int32)):
        if int(cell.area) < int(min_area):
            continue
        current += 1
        coords = cell.coords
        out[coords[:, 0], coords[:, 1]] = current
    return out


def _remove_pomc_npy_overlaps(
    pomc_labels: np.ndarray, npy_labels: np.ndarray, threshold: float,
) -> Tuple[np.ndarray, int]:
    out = np.asarray(pomc_labels, dtype=np.int32).copy()
    npy_positive = np.asarray(npy_labels) > 0
    removed = 0
    for cell in regionprops(np.asarray(pomc_labels, dtype=np.int32)):
        coords = cell.coords
        fraction = float(np.count_nonzero(npy_positive[coords[:, 0], coords[:, 1]])) / float(
            max(1, cell.area)
        )
        if fraction >= float(threshold):
            out[coords[:, 0], coords[:, 1]] = 0
            removed += 1
    return component_label(out > 0, connectivity=1).astype(np.int32), removed


def load_cellpose_cartoon_masks(
    registration: Mapping[str, object], folder: Path, root: Path,
    sources: Mapping[str, Path],
) -> Tuple[Dict[str, np.ndarray], Dict[str, Dict[str, object]]]:
    """Rebuild analyzed DAPI/POMC/c-FOS masks from raw NPYs and prove identity."""
    registration_input = registration.get("input")
    transform = registration.get("transform")
    if not isinstance(registration_input, Mapping) or not isinstance(transform, Mapping):
        raise ValueError("Panel F requires registration input and transform provenance")

    stitched: Dict[str, np.ndarray] = {}
    provenance: Dict[str, Dict[str, object]] = {}
    for display_name, marker_key in CELLPOSE_MARKERS.items():
        halves: Dict[str, Optional[np.ndarray]] = {"left": None, "right": None}
        record: Dict[str, object] = {}
        for side in ("left", "right"):
            path_text = str(registration_input.get(f"{side}_{marker_key}_seg_path", "")).strip()
            digest = str(registration_input.get(f"{side}_{marker_key}_seg_sha256", "")).strip().lower()
            if not path_text:
                record[f"{side}_path"] = ""
                record[f"{side}_sha256"] = ""
                record[f"{side}_raw_instances"] = 0
                continue
            path = _resolve_provenance_path(path_text, folder=folder, root=root)
            if not path.is_file() or not _valid_sha256(digest) or sha256_file(path) != digest:
                raise ValueError(f"Panel F raw {display_name} {side} NPY path/hash mismatch")
            labels = load_raw_cellpose_labels(path)
            expected_shape = registration_input.get(f"{side}_dapi_shape_px")
            if not isinstance(expected_shape, list) or list(labels.shape) != [int(v) for v in expected_shape]:
                raise ValueError(f"Panel F raw {display_name} {side} mask shape mismatch")
            halves[side] = labels
            record[f"{side}_path"] = str(path)
            record[f"{side}_sha256"] = digest
            record[f"{side}_raw_instances"] = int(np.count_nonzero(np.unique(labels) > 0))
        stitched[display_name] = stitch_raw_cellpose_labels(
            halves["left"], halves["right"], transform,
        )
        if tuple(stitched[display_name].shape) != _registration_output_shape(registration):
            raise ValueError(f"Panel F stitched {display_name} mask shape differs from Panel D")
        provenance[display_name] = record

    processed = {
        "DAPI": stitched["DAPI"],
        "c-FOS": _filter_labels_by_area(stitched["c-FOS"], min_area=8),
        "NPY": stitched["NPY"],
    }
    processed["POMC"], pomc_removed = _remove_pomc_npy_overlaps(
        stitched["POMC"], processed["NPY"], threshold=0.20,
    )
    processing = {
        "DAPI": "raw Cellpose instances; audited hard-seam stitch",
        "c-FOS": "raw Cellpose instances; audited hard-seam stitch; minimum area 8 px",
        "POMC": "raw Cellpose instances; audited hard-seam stitch; remove NPY overlap >= 0.20",
        "NPY": "raw Cellpose instances; audited hard-seam stitch; used only for POMC cleanup",
    }
    source_keys = {"DAPI": "dapi_labels", "c-FOS": "cfos_labels", "POMC": "pomc_labels", "NPY": "npy_labels"}
    for name, labels in processed.items():
        expected = read_gray(sources[source_keys[name]]).astype(np.int32)
        if not np.array_equal(labels, expected):
            differing = int(np.count_nonzero(labels != expected))
            raise ValueError(
                f"Panel F raw-NPY reconstruction differs from analyzed {name} labels at {differing} pixels"
            )
        positive = labels > 0
        provenance[name].update({
            "processing": processing[name],
            "processed_instances": int(np.count_nonzero(np.unique(labels) > 0)),
            "positive_pixels": int(np.count_nonzero(positive)),
            "registered_binary_sha256": hashlib.sha256(
                np.ascontiguousarray(positive, dtype=np.uint8).tobytes()
            ).hexdigest(),
            "matches_registered_analysis_labels": True,
        })
    provenance["POMC"]["raw_rois_removed_by_npy_overlap"] = int(pomc_removed)
    processed["POMC"], size_qc, size_path = load_size_filtered_pomc(sources)
    provenance["POMC"]["processed_instances_before_size_qc"] = provenance["POMC"]["processed_instances"]
    provenance["POMC"]["processed_instances"] = size_qc["objects_after"]
    provenance["POMC"]["post_segmentation_size_qc"] = size_qc
    provenance["POMC"]["size_filtered_label_source"] = str(size_path)
    return {name: processed[name] > 0 for name in ("DAPI", "POMC", "c-FOS")}, provenance


def hil_arc_ventricle_lumen(
    regions: np.ndarray,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """Derive the ventricular corridor from the two accepted HIL ARC lobes.

    The medial edge of the left ARC and medial edge of the right ARC define the
    lumen row by row. Their first valid positions are extended to the cropped
    image top, correcting the unsupported upper DAPI-density closure. The
    corridor stops at the first accepted HIL ME row. No automated wall curve or
    marker intensity participates.
    """
    regions = np.asarray(regions, dtype=np.uint8)
    arc_components, count = ndi_label(regions == 1)
    ranked = sorted(
        range(1, count + 1),
        key=lambda code: int(np.count_nonzero(arc_components == code)),
        reverse=True,
    )
    if len(ranked) < 2:
        raise ValueError("Panel D HIL ventricle correction requires two accepted ARC lobes")
    two = ranked[:2]
    two.sort(key=lambda code: float(np.mean(np.where(arc_components == code)[1])))
    left_arc = arc_components == two[0]
    right_arc = arc_components == two[1]
    height, width = regions.shape
    medial_rows: List[Tuple[int, float, float]] = []
    for y in range(height):
        left_x = np.flatnonzero(left_arc[y])
        right_x = np.flatnonzero(right_arc[y])
        if left_x.size and right_x.size and int(left_x.max()) < int(right_x.min()):
            medial_rows.append((y, float(left_x.max()), float(right_x.min())))
    if len(medial_rows) < 2:
        raise ValueError("Accepted HIL ARC lobes do not define a ventricular corridor")
    valid_y = np.asarray([row[0] for row in medial_rows], dtype=float)
    left_medial = np.asarray([row[1] for row in medial_rows], dtype=float)
    right_medial = np.asarray([row[2] for row in medial_rows], dtype=float)
    me_rows = np.flatnonzero(np.any(regions == 2, axis=1))
    stop_y = int(me_rows.min()) if me_rows.size else int(valid_y.max()) + 1
    target_y = np.arange(0, min(stop_y, height), dtype=float)
    left_interp = np.interp(target_y, valid_y, left_medial)
    right_interp = np.interp(target_y, valid_y, right_medial)
    lumen = np.zeros((height, width), dtype=bool)
    for y, left_x, right_x in zip(target_y.astype(int), left_interp, right_interp):
        x0 = max(0, int(math.floor(left_x)) + 1)
        x1 = min(width, int(math.ceil(right_x)))
        if x1 > x0:
            lumen[y, x0:x1] = True
    # Accepted HIL regions always win; this guard makes an accidental overlap
    # impossible even if future HIL polygons touch across the corridor.
    lumen[regions > 0] = False
    info: Dict[str, object] = {
        "source": "medial boundaries of the two accepted HIL ARC lobes",
        "uses_automated_wall_csv": False,
        "uses_marker_intensity": False,
        "first_bilateral_arc_row_px": int(valid_y.min()),
        "top_extension_rows_px": int(valid_y.min()),
        "stop_before_first_hil_me_row_px": stop_y,
        "bilateral_hil_arc_rows": int(len(valid_y)),
        "lumen_pixels": int(np.count_nonzero(lumen)),
    }
    return lumen, info


def cellpose_tissue_envelope(
    dapi_mask: np.ndarray, regions: np.ndarray, hil_ventricle_lumen: np.ndarray,
) -> np.ndarray:
    """Create Panel D DAPI tissue support, then enforce the HIL ARC lumen."""
    dapi_positive = np.asarray(dapi_mask, dtype=bool)
    if not (
        dapi_positive.shape == np.asarray(regions).shape
        == np.asarray(hil_ventricle_lumen).shape
    ):
        raise ValueError("Panel D DAPI and accepted-region masks must share coordinates")
    sigma = max(6.0, 0.025 * float(dapi_positive.shape[0]))
    density = gaussian_filter(dapi_positive.astype(np.float32), sigma=sigma, mode="constant")
    tissue = remove_small_objects(
        density > 0.02, min_size=max(256, int(round(dapi_positive.size * 0.002))),
    )
    tissue |= np.asarray(regions) > 0
    tissue[np.asarray(hil_ventricle_lumen, dtype=bool)] = False
    if not np.any(tissue):
        raise ValueError("Panel D could not derive a DAPI-supported tissue envelope")
    return tissue


def region_label_points(regions: np.ndarray) -> List[Tuple[str, float, float]]:
    """Return stable interior points: both ARC lobes and the single accepted ME."""
    points: List[Tuple[str, float, float]] = []
    for region_name, code, keep_count in (("ARC", 1, 2), ("ME", 2, 1)):
        components, count = ndi_label(np.asarray(regions) == code)
        ranked = sorted(
            range(1, count + 1),
            key=lambda component: int(np.count_nonzero(components == component)),
            reverse=True,
        )[:keep_count]
        for component in ranked:
            inside = components == component
            distance = distance_transform_edt(inside)
            y, x = np.unravel_index(int(np.argmax(distance)), distance.shape)
            points.append((region_name, float(x), float(y)))
    return points


def draw_cellpose_cartoon(
    ax: plt.Axes, marker: str, marker_mask: np.ndarray,
    tissue: np.ndarray, regions: np.ndarray, args: argparse.Namespace,
) -> None:
    """Draw one exact-coordinate, vector-like channel cartoon for Panel D."""
    shape = np.asarray(marker_mask).shape
    if not (shape == np.asarray(tissue).shape == np.asarray(regions).shape):
        raise ValueError("Panel D cartoon layers do not share one registered field")
    rgb = np.ones((*shape, 3), dtype=np.float32)
    tissue_rgb = np.asarray(matplotlib.colors.to_rgb(TISSUE_FILL), dtype=np.float32)
    rgb[np.asarray(tissue, dtype=bool)] = tissue_rgb
    for region_name, code in (("ARC", 1), ("ME", 2)):
        inside = np.asarray(regions) == code
        color = np.asarray(matplotlib.colors.to_rgb(REGION_COLORS[region_name]), dtype=np.float32)
        rgb[inside] = 0.88 * rgb[inside] + 0.12 * color
    ax.imshow(np.clip(rgb, 0, 1), interpolation="nearest")

    cell_layer = np.zeros((*shape, 4), dtype=np.float32)
    cell_layer[..., :3] = matplotlib.colors.to_rgb(CHANNEL_COLORS[marker])
    cell_layer[..., 3] = np.asarray(marker_mask, dtype=np.float32) * (0.78 if marker == "DAPI" else 0.90)
    ax.imshow(cell_layer, interpolation="nearest")
    ax.contour(tissue, [0.5], colors=["white"], linewidths=3.5, alpha=0.95)
    ax.contour(tissue, [0.5], colors=[TISSUE_EDGE], linewidths=1.7)
    for region_name, code in (("ARC", 1), ("ME", 2)):
        inside = np.asarray(regions) == code
        if not np.any(inside):
            continue
        ax.contour(inside, [0.5], colors=["#172033"], linewidths=3.2, alpha=0.92)
        ax.contour(inside, [0.5], colors=[REGION_COLORS[region_name]], linewidths=1.65)
    for region_name, x, y in region_label_points(regions):
        ax.text(
            x, y, region_name, ha="center", va="center",
            fontsize=7.2, fontweight="bold", color=REGION_COLORS[region_name], zorder=20,
            bbox={"boxstyle": "round,pad=0.20", "facecolor": "white",
                  "edgecolor": REGION_COLORS[region_name], "linewidth": 0.8, "alpha": 0.92},
        )
    if marker == "DAPI":
        ax.text(
            0.016, 0.035, "Tissue boundary", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=6.3, fontweight="bold", color=TISSUE_EDGE,
            bbox={"boxstyle": "round,pad=0.20", "facecolor": "white",
                  "edgecolor": TISSUE_EDGE, "linewidth": 0.7, "alpha": 0.90},
            zorder=22,
        )
    add_scalebar(ax, shape, args.um_per_px, args.scalebar_um, font_size=6.2)
    ax.set_title(marker, color=CHANNEL_COLORS[marker], fontsize=10.5, fontweight="bold", pad=3)
    ax.set_xlim(-0.5, float(shape[1]) - 0.5)
    ax.set_ylim(float(shape[0]) - 0.5, -0.5)
    ax.set_axis_off()


def robust_normalize(arr: np.ndarray, gamma: float = 0.82) -> np.ndarray:
    x = np.asarray(arr, dtype=np.float32)
    finite = x[np.isfinite(x)]
    positive = finite[finite > 0]
    basis = positive if positive.size >= 100 else finite
    if basis.size == 0:
        return np.zeros(x.shape, dtype=np.float32)
    lo, hi = np.percentile(basis, [0.5, 99.7])
    if not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    return np.power(np.clip((x - lo) / (hi - lo), 0, 1), gamma).astype(np.float32)


def single_channel_view(channel: np.ndarray, marker: str) -> np.ndarray:
    """One isolated constituent of the merged field, in its composite tint.

    raw_composite is the clipped sum of these three views over the same
    registered pixels, so Panels A-C are exactly what gives rise to Panel E.
    """
    if marker not in COMPOSITE_CHANNEL_TINTS:
        raise ValueError(f"Unknown microscopy channel: {marker!r}")
    tint = np.asarray(COMPOSITE_CHANNEL_TINTS[marker], dtype=np.float64)
    return np.clip(robust_normalize(channel)[..., None] * tint, 0, 1)


def raw_composite(dapi: np.ndarray, cfos: np.ndarray, pomc: np.ndarray) -> np.ndarray:
    if not (dapi.shape == cfos.shape == pomc.shape):
        raise ValueError(f"Raw channel shapes differ: {dapi.shape}, {cfos.shape}, {pomc.shape}")
    channels = {"DAPI": dapi, "c-FOS": cfos, "POMC": pomc}
    rgb = np.zeros((*np.asarray(dapi).shape, 3), dtype=np.float64)
    for marker in COMPOSITE_CHANNEL_ORDER:
        tint = np.asarray(COMPOSITE_CHANNEL_TINTS[marker], dtype=np.float64)
        rgb += robust_normalize(channels[marker])[..., None] * tint
    return np.clip(rgb, 0, 1)


def npy_gfp_display_normalize(arr: np.ndarray, gamma: float = 0.72) -> np.ndarray:
    """Background-window NPY-GFP while retaining native signal ordering.

    This channel has a nonzero whole-field camera background (representative
    field median 40 on uint8) whereas POMC contains true zero-background pixels.
    Using the generic positive-pixel 0.5th percentile therefore paints ordinary
    tissue green. A deterministic median/99.7th display window suppresses that
    background but applies no denoising, mask, spatial filter, or resampling.
    """
    x = np.asarray(arr, dtype=np.float32)
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return np.zeros(x.shape, dtype=np.float32)
    lo, hi = np.percentile(finite, [50.0, 99.7])
    if not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    return np.power(np.clip((x - lo) / (hi - lo), 0, 1), gamma).astype(np.float32)


def roi_microscopy_view(
    channel: np.ndarray, roi_labels: np.ndarray, marker: str,
) -> np.ndarray:
    """Retain a full-field microscopy display only at exact marker-ROI pixels.

    Normalize before applying the mask, so changing ROI membership cannot
    stretch the intensity of another ROI. Label values provide support only;
    all displayed color and within-cell variation come from the source channel.
    """
    channel = np.asarray(channel)
    roi_labels = np.asarray(roi_labels)
    if channel.ndim != 2 or channel.shape != roi_labels.shape:
        raise ValueError("Microscopy and ROI labels must share one 2-D registered field")
    if marker == "POMC":
        rgb = single_channel_view(channel, marker)
    elif marker == "NPY-GFP":
        rgb = (npy_gfp_display_normalize(channel)[..., None]
               * np.asarray(NPY_GFP_COMPOSITE_TINT, dtype=np.float64))
    else:
        raise ValueError(f"Unsupported ROI microscopy channel: {marker!r}")
    return np.where((roi_labels > 0)[..., None], rgb, 0.0)


def composite_with_roi_masked_pomc_signal_and_npy_gfp(
    dapi: np.ndarray, cfos: np.ndarray, pomc: np.ndarray, npy_gfp: np.ndarray,
    pomc_roi_labels: np.ndarray,
) -> np.ndarray:
    """Panel F microscopy with real POMC intensity gated by accepted ROIs.

    DAPI, c-FOS and POMC use the same full-field channel normalization shown in
    A-C/E. The exact Panel C POMC display is retained only where the registered
    processed POMC label TIFF is > 0, then screen-composited over DAPI/c-FOS.
    NPY-GFP receives its audited display window over that composite. No extra
    POMC gain, ROI-specific normalization, flat ROI paint, outline, dilation,
    spatial filter, interpolation, or resampling is used.
    """
    pomc_roi_labels = np.asarray(pomc_roi_labels)
    if not (
        dapi.shape == cfos.shape == pomc.shape == npy_gfp.shape
        == pomc_roi_labels.shape
    ):
        raise ValueError(
            "Registered DAPI/c-FOS/POMC/NPY-GFP/POMC-ROI shapes differ: "
            f"{dapi.shape}, {cfos.shape}, {pomc.shape}, {npy_gfp.shape}, "
            f"{pomc_roi_labels.shape}"
        )
    pomc_roi_mask = pomc_roi_labels > 0
    if not pomc_roi_mask.any():
        raise ValueError("Panel F POMC ROI label TIFF contains no selected ROIs")
    dapi_cfos = np.clip(
        single_channel_view(dapi, "DAPI")
        + single_channel_view(cfos, "c-FOS"),
        0, 1,
    )
    pomc_signal = roi_microscopy_view(pomc, pomc_roi_labels, "POMC")
    base = 1.0 - (1.0 - dapi_cfos) * (1.0 - pomc_signal)
    base[~pomc_roi_mask] = dapi_cfos[~pomc_roi_mask]
    npy_rgb = (
        npy_gfp_display_normalize(npy_gfp)[..., None]
        * np.asarray(NPY_GFP_COMPOSITE_TINT, dtype=np.float64)
    )
    screened = np.clip(1.0 - (1.0 - base) * (1.0 - npy_rgb), 0, 1)
    result = np.maximum(base, screened)
    result[..., 0] = base[..., 0]
    return result


def draw_single_channel_panel(
    ax: plt.Axes, channel: np.ndarray, marker: str, args: argparse.Namespace,
) -> None:
    """Draw one isolated composite constituent in the exact registered field."""
    shape = np.asarray(channel).shape
    ax.imshow(single_channel_view(channel, marker), interpolation="nearest")
    add_scalebar(ax, shape, args.um_per_px, args.scalebar_um, font_size=6.2)
    ax.set_title(marker, color=CHANNEL_COLORS[marker], fontsize=10.5, fontweight="bold", pad=3)
    ax.set_xlim(-0.5, float(shape[1]) - 0.5)
    ax.set_ylim(float(shape[0]) - 0.5, -0.5)
    ax.set_axis_off()

def accepted_segmentation_overlay(dapi: np.ndarray, regions: np.ndarray) -> np.ndarray:
    """Show native DAPI beneath the accepted ARC/ME labels in the same pixels."""
    if dapi.shape != regions.shape:
        raise ValueError(f"Native DAPI and accepted ARC/ME mask differ: {dapi.shape} vs {regions.shape}")
    base = robust_normalize(dapi, gamma=0.72)
    rgb = base[..., None] * np.asarray([0.18, 0.34, 0.92], dtype=np.float32)
    for code, color in (
        (1, np.asarray([0.00, 0.74, 0.83], dtype=np.float32)),
        (2, np.asarray([1.00, 0.60, 0.00], dtype=np.float32)),
        (3, np.asarray([0.75, 0.24, 0.68], dtype=np.float32)),
    ):
        inside = regions == code
        rgb[inside] = 0.70 * rgb[inside] + 0.30 * color
        boundary = find_boundaries(inside, connectivity=2, mode="thick")
        boundary = binary_dilation(boundary, disk(2))
        rgb[boundary] = color
    return np.clip(rgb, 0, 1)



def dapi_pomc_composite(dapi: np.ndarray, pomc: np.ndarray) -> np.ndarray:
    d, p = robust_normalize(dapi, 0.76), robust_normalize(pomc, 0.76)
    return np.clip(d[..., None] * [0.10, 0.25, 1.0] + p[..., None] * [1.00, 0.65, 0.12], 0, 1)


def classification_overlay(labels: np.ndarray, region_mask: np.ndarray, cells: pd.DataFrame) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.int32); regions = np.asarray(region_mask, dtype=np.uint8)
    if labels.shape != regions.shape:
        raise ValueError(f"DAPI-label and region-mask shapes differ: {labels.shape}, {regions.shape}")
    lut = np.zeros(int(labels.max()) + 1, dtype=np.uint8)
    for row in cells.itertuples(index=False):
        lab = int(getattr(row, "nucleus_label"))
        if not (0 < lab < lut.size):
            continue
        cfos = str(getattr(row, "is_cfos")).lower() in {"true", "1"}
        pomc = str(getattr(row, "is_pomc")).lower() in {"true", "1"}
        lut[lab] = 3 if cfos and pomc else 1 if cfos else 2 if pomc else 0
    cls = lut[labels]
    rgb = np.zeros((*labels.shape, 3), dtype=np.float32)
    rgb[labels > 0] = [0.42, 0.42, 0.47]
    rgb[cls == 1] = [1.00, 0.02, 0.82]
    # Amber POMC, matching COMPOSITE_CHANNEL_TINTS. The double-positive class
    # moves to white because yellow sat too close to amber to be told apart.
    rgb[cls == 2] = [1.00, 0.65, 0.12]
    rgb[cls == 3] = [1.00, 1.00, 1.00]
    rgb[find_boundaries(labels, mode="inner")] = np.maximum(rgb[find_boundaries(labels, mode="inner")], 0.92)
    for code, color in [(1, np.asarray([0.0, 0.55, 0.70])), (2, np.asarray([0.76, 0.25, 0.05])), (3, np.asarray([0.75, 0.24, 0.68]))]:
        mask = (regions == code) & (labels == 0)
        rgb[mask] = 0.07 * color
    return np.clip(rgb, 0, 1)


def add_scalebar(ax: plt.Axes, shape: Tuple[int, int], um_per_px: float, length_um: float, *, font_size: float = 12.0, inset: bool = False) -> None:
    """High-contrast white scale bar on a wide opaque black backing."""
    height, width = shape
    length_px = float(length_um) / float(um_per_px)
    if length_px <= 0 or length_px >= 0.48 * width:
        raise ValueError(f"Invalid scale bar: {length_um} µm -> {length_px:.1f}px for width {width}")
    margin_x = max(8.0 if inset else 30.0, (0.04 if inset else 0.026) * width)
    margin_y = max(8.0 if inset else 25.0, (0.06 if inset else 0.055) * height)
    bar_h = max(4.0 if inset else 11.0, (0.028 if inset else 0.017) * height)
    x0 = width - margin_x - length_px
    y0 = height - margin_y
    pad_x = max(5.0, 0.025 * length_px)
    text_h = max(16.0, 0.060 * height)
    if not inset:
        ax.add_patch(Rectangle(
            (x0 - pad_x, y0 - bar_h - text_h), length_px + 2 * pad_x, bar_h + text_h + 7,
            facecolor="black", edgecolor="black", linewidth=1.2, alpha=0.92, zorder=30,
        ))
    # The inset scale bar deliberately has no filled backing: the image remains visible
    # through and around the 20 µm annotation.
    ax.add_patch(Rectangle(
        (x0, y0 - bar_h), length_px, bar_h,
        facecolor="white", edgecolor="black" if inset else "white",
        linewidth=0.8 if inset else 1.0, zorder=31,
    ))
    txt = ax.text(x0 + length_px / 2, y0 - bar_h - 4, f"{length_um:g} µm", color="white",
                  ha="center", va="bottom", fontsize=font_size, fontweight="bold", zorder=32)
    txt.set_path_effects([pe.withStroke(linewidth=3.2 if not inset else 2.0, foreground="black")])


def draw_wall_geometry(ax: plt.Axes, sources: Mapping[str, Path]) -> None:
    if sources["walls"].exists():
        wall = pd.read_csv(sources["walls"])
        if {"y_px", "left_wall_x_px", "right_wall_x_px"}.issubset(wall.columns):
            y = pd.to_numeric(wall["y_px"], errors="coerce").to_numpy(float)
            for col in ["left_wall_x_px", "right_wall_x_px"]:
                x = pd.to_numeric(wall[col], errors="coerce").to_numpy(float)
                good = np.isfinite(x) & np.isfinite(y)
                ax.plot(x[good], y[good], color="#061b5c", linewidth=4.2, zorder=12)
                ax.plot(x[good], y[good], color="white", linewidth=1.8, zorder=13)
    if sources["floor"].exists():
        floor = pd.read_csv(sources["floor"])
        if {"x_px", "floor_y_px"}.issubset(floor.columns):
            x = pd.to_numeric(floor["x_px"], errors="coerce").to_numpy(float)
            y = pd.to_numeric(floor["floor_y_px"], errors="coerce").to_numpy(float)
            good = np.isfinite(x) & np.isfinite(y)
            ax.plot(x[good], y[good], color=REGION_COLORS["ME"], linewidth=2.2, linestyle="--", zorder=14)


def add_panel_label(ax: plt.Axes, label: str, size: float = 20.0) -> None:
    ax.text(-0.015, 1.025, label, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=size, fontweight="bold", clip_on=False)


def add_field_identity(ax: plt.Axes, chosen: pd.Series, font_size: float = 11.0) -> None:
    """Name the animal, section and condition the displayed field comes from."""
    sample = str(chosen.get("sample", "")).strip()
    section = chosen.get("section_index", "")
    condition = str(chosen.get("cond", "")).strip()
    try:
        section_text = f"S{int(section):02d}"
    except (TypeError, ValueError):
        section_text = str(section)
    label = "  ·  ".join(part for part in (f"{sample} {section_text}".strip(), condition) if part)
    txt = ax.text(
        0.018, 0.980, label, transform=ax.transAxes, ha="left", va="top",
        color="white", fontsize=font_size, fontweight="bold", zorder=43,
    )
    txt.set_path_effects([pe.withStroke(linewidth=2.8, foreground="black")])


def add_marker_distribution_inset(
    ax: plt.Axes, pomc: np.ndarray, npy_gfp: np.ndarray,
    pomc_roi_labels: np.ndarray, npy_gfp_roi_labels: np.ndarray,
    tissue: np.ndarray, ventricle_lumen: np.ndarray,
    um_per_px: float, scale_bar_um: float,
) -> Dict[str, object]:
    """Show microscopy within exact ROIs, with tissue and ventricular guides.

    These are the processed registered marker-label TIFFs, not DAPI-nucleus
    surrogates. Each marker uses its full-field microscopy display before masking,
    identical to its contribution to the main field. The black background carries
    no tissue or ARC/ME fill. Marker ROIs have no contour or uniform fill.
    Thin white dashed lines delimit the DAPI-supported
    tissue exterior and the HIL-ARC-derived ventricular lumen. There is no halo,
    dilation, smoothing, or resampling.
    """
    pomc_roi_labels = np.asarray(pomc_roi_labels, dtype=np.int32)
    npy_gfp_roi_labels = np.asarray(npy_gfp_roi_labels, dtype=np.int32)
    tissue = np.asarray(tissue, dtype=bool)
    ventricle_lumen = np.asarray(ventricle_lumen, dtype=bool)
    if not (
        pomc_roi_labels.shape == npy_gfp_roi_labels.shape == tissue.shape
        == ventricle_lumen.shape
    ):
        raise ValueError(
            "Panel F ROI/tissue/ventricle shapes differ: "
            f"{pomc_roi_labels.shape}, {npy_gfp_roi_labels.shape}, "
            f"{tissue.shape}, {ventricle_lumen.shape}"
        )
    pomc_mask = pomc_roi_labels > 0
    npy_gfp_mask = npy_gfp_roi_labels > 0
    counts = {
        "POMC": int(np.count_nonzero(np.unique(pomc_roi_labels))),
        "NPY-GFP": int(np.count_nonzero(np.unique(npy_gfp_roi_labels))),
    }
    if not pomc_mask.any() or not npy_gfp_mask.any():
        raise ValueError("Panel F requires nonempty POMC and NPY-GFP ROI label TIFFs")
    pomc_rgb = roi_microscopy_view(pomc, pomc_roi_labels, "POMC")
    npy_rgb = roi_microscopy_view(npy_gfp, npy_gfp_roi_labels, "NPY-GFP")
    # Screen both measured signals at shared pixels without replacing either
    # channel by a label color or clipping their sum to an opaque patch.
    rgb = 1.0 - (1.0 - pomc_rgb) * (1.0 - npy_rgb)

    # Forty percent is modestly smaller than the preceding 46% inset while still
    # rendering the native ROI interiors clearly at the fixed 600 dpi export.
    iax = inset_axes(ax, width="40.0%", height="40.0%", loc="upper right", borderpad=0.55)
    iax.imshow(np.clip(rgb, 0, 1), interpolation="nearest")
    # Fill the carved ventricular corridor back into the tissue support only
    # for the outer silhouette. The lumen itself is then delimited separately.
    tissue_outer = tissue | ventricle_lumen
    iax.contour(
        tissue_outer.astype(np.uint8), levels=[0.5], colors=["white"],
        linewidths=0.85, linestyles="--", alpha=0.96, zorder=5,
    )
    iax.contour(
        ventricle_lumen.astype(np.uint8), levels=[0.5], colors=["white"],
        linewidths=1.05, linestyles="--", alpha=0.98, zorder=6,
    )
    for spine in iax.spines.values():
        spine.set_visible(False)
    iax.set_xticks([]); iax.set_yticks([])
    add_scalebar(
        iax, pomc_roi_labels.shape, um_per_px, scale_bar_um, font_size=8.0, inset=True,
    )
    entries = [
        ("POMC ROIs", POMC_MERGE_LABEL_COLOR, counts["POMC"]),
        ("NPY-GFP ROIs", NPY_GFP_INSET_COLOR, counts["NPY-GFP"]),
    ]
    for index, (name, color, count) in enumerate(entries):
        iax.text(
            0.028, 0.965 - 0.125 * index, f"{name} {count}", transform=iax.transAxes,
            ha="left", va="top", color=color, fontsize=8.0, fontweight="bold",
            path_effects=[pe.withStroke(linewidth=1.6, foreground="black")], zorder=8,
        )
    return {
        "pomc_rois": counts["POMC"],
        "npy_gfp_rois": counts["NPY-GFP"],
        "pomc_roi_pixels": int(np.count_nonzero(pomc_mask)),
        "npy_gfp_roi_pixels": int(np.count_nonzero(npy_gfp_mask)),
        "overlap_pixels": int(np.count_nonzero(pomc_mask & npy_gfp_mask)),
        "morphological_dilation_px": 0,
        "roi_contour_width_pt": 0.0,
        "anatomical_contour_width_pt": 1.05,
        "tissue_background_drawn": False,
        "tissue_boundary_drawn": True,
        "tissue_boundary_color": "#ffffff",
        "tissue_boundary_width_pt": 0.85,
        "tissue_boundary_linestyle": "dashed",
        "ventricle_boundary_drawn": True,
        "ventricle_boundary_color": "#ffffff",
        "ventricle_boundary_width_pt": 1.05,
        "ventricle_boundary_linestyle": "dashed",
        "arc_me_region_contours_drawn": False,
        "display": "registered microscopy intensities within exact marker ROIs",
        "uniform_roi_fill": False,
        "channel_display_settings_shared_with_main_field": True,
        "marker_blend": "screen",
        "npy_gfp_tint_rgb": list(NPY_GFP_COMPOSITE_TINT),
        "pomc_tint_rgb": list(COMPOSITE_CHANNEL_TINTS["POMC"]),
        "width_percent": 40.0,
        "height_percent": 40.0,
        "source": (
            "registered POMC/NPY-GFP microscopy TIFFs gated by their exact "
            "processed label TIFF interiors; "
            "DAPI tissue envelope; HIL-ARC-derived ventricular lumen"
        ),
    }


def add_merge_channel_key(
    ax: plt.Axes, font_size: float = 12.0, *, include_npy_gfp: bool = False,
    pomc_rois_only: bool = False,
) -> None:
    """Large vertical channel key over the upper-left native field."""
    labels = []
    if include_npy_gfp:
        labels.append(("NPY-GFP", NPY_GFP_INSET_COLOR))
    labels.extend([
        ("POMC (within ROIs)" if pomc_rois_only else "POMC",
         POMC_MERGE_LABEL_COLOR),
        ("c-FOS", "#ff37d4"),
        ("DAPI", "#4f86ff"),
    ])
    entries = [
        (0.018, 0.885 - 0.085 * index, label, color)
        for index, (label, color) in enumerate(labels)
    ]
    for x, y, label, color in entries:
        txt = ax.text(
            x, y, label, transform=ax.transAxes,
            ha="left", va="top", color=color, fontsize=font_size,
            fontweight="bold", zorder=42,
        )
        txt.set_path_effects([pe.withStroke(linewidth=2.8, foreground="black")])

def add_region_color_key(
    ax: plt.Axes, displayed_regions: Sequence[str], font_size: float = 7.4,
) -> None:
    """Key only the ROI subset actually present in the displayed mask."""
    handles = [
        Line2D(
            [0], [0], color=REGION_COLORS[region], linewidth=3.4,
            path_effects=[pe.Stroke(linewidth=5.8, foreground="black"), pe.Normal()],
            label=f"{region} segmentation",
        )
        for region in displayed_regions
    ]
    if not handles:
        raise ValueError("Panel E cannot key an empty accepted ROI subset")
    legend = ax.legend(
        handles=handles, loc="lower left", bbox_to_anchor=(0.008, 0.018),
        ncol=3, frameon=True, facecolor="black", edgecolor="white",
        framealpha=0.82, fontsize=font_size, labelcolor="white",
        handlelength=1.45, handletextpad=0.42, columnspacing=0.8, borderpad=0.38,
    )
    legend.set_zorder(40)


def add_review_provenance(ax: plt.Axes, provenance: Mapping[str, object]) -> None:
    if not bool(provenance.get("human_reviewed", False)):
        text = (
            "Current automated ARC/ME\n"
            "HIL review: not requested\n"
            "Source: native POMC segmentation\n"
            f"Mask SHA-256: {str(provenance['mask_sha256'])[:16]}…"
        )
        ax.text(
            0.992, 0.988, text, transform=ax.transAxes,
            ha="right", va="top", color="white", fontsize=6.7,
            fontweight="bold", linespacing=1.18, zorder=45, clip_on=False,
            bbox={"boxstyle": "round,pad=0.38", "facecolor": "black",
                  "edgecolor": "white", "linewidth": 0.8, "alpha": 0.78},
        )
        return
    drawn = tuple(
        region for region in ("ARC", "ME", "VMN")
        if provenance.get("region_status", {}).get(region) == "drawn"
    )
    if not drawn:
        raise ValueError("HIL-reviewed panel provenance has no drawn ROI")
    text = f"HIL-accepted ROIs: {', '.join(drawn)}"
    ax.text(
        0.992, 0.988, text, transform=ax.transAxes,
        ha="right", va="top", color="white", fontsize=8.2,
        fontweight="bold", zorder=45, clip_on=False,
        bbox={"boxstyle": "round,pad=0.30", "facecolor": "black",
              "edgecolor": "white", "linewidth": 0.8, "alpha": 0.78},
    )


def draw_bold_region_contours(
    ax: plt.Axes,
    labels: np.ndarray,
    displayed_regions: Sequence[str],
) -> None:
    """Draw publication-visible ROI boundaries with a dark contrast halo."""
    code_by_region = {"ARC": 1, "ME": 2, "VMN": 3}
    for region in displayed_regions:
        inside = np.asarray(labels) == code_by_region[region]
        if not np.any(inside):
            raise ValueError(f"Displayed {region} status conflicts with an empty mask")
        # The underlay prevents cyan/orange/magenta from disappearing over bright
        # DAPI; the colored stroke is the anatomical boundary reported in the key.
        ax.contour(inside, [0.5], colors=["black"], linewidths=6.4, alpha=0.92)
        ax.contour(
            inside, [0.5], colors=[REGION_COLORS[region]], linewidths=3.5,
        )


def select_pomc_inset(cells: pd.DataFrame, labels: np.ndarray, cfos: np.ndarray, pomc: np.ndarray, field_px: int) -> Tuple[Dict[str, object], Tuple[int, int, int, int]]:
    """Select an isolated, representative POMC+/c-FOS+ nucleus for the inset."""
    candidates = cells[as_bool(cells["is_pomc"]) & as_bool(cells["is_cfos"])].copy()
    if candidates.empty:
        raise ValueError("No POMC+/c-FOS+ nuclei are available for the requested inset.")
    if "region" in candidates.columns and candidates["region"].astype(str).eq("ARC").any():
        candidates = candidates[candidates["region"].astype(str).eq("ARC")].copy()
    half = max(16, int(field_px // 2)); h, w = labels.shape
    coords = candidates[["centroid_y", "centroid_x"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    scored: List[Dict[str, object]] = []
    for row in candidates.itertuples(index=False):
        lab = int(getattr(row, "nucleus_label")); y = int(round(float(getattr(row, "centroid_y")))); x = int(round(float(getattr(row, "centroid_x"))))
        if x - half < 0 or x + half >= w or y - half < 0 or y + half >= h or abs(x - w / 2) < half + 12:
            continue
        pomc_patch = np.asarray(pomc[y-half:y+half, x-half:x+half], dtype=float)
        cfos_patch = np.asarray(cfos[y-half:y+half, x-half:x+half], dtype=float)
        pomc_positive = pomc_patch[pomc_patch > 0]
        cfos_positive = cfos_patch[cfos_patch > 0]
        pomc_score = float(np.percentile(pomc_positive, 95)) if pomc_positive.size else 0.0
        cfos_score = float(np.percentile(cfos_positive, 95)) if cfos_positive.size else 0.0
        distances = np.sqrt((coords[:, 0] - y) ** 2 + (coords[:, 1] - x) ** 2)
        neighbor_count = int(np.count_nonzero((distances > 1.0) & (distances < 0.62 * field_px)))
        nearest_other = float(np.min(distances[distances > 1.0])) if np.any(distances > 1.0) else float("inf")
        scored.append({"nucleus_label": lab, "centroid_y": y, "centroid_x": x,
                       "region": getattr(row, "region", ""), "is_cfos": getattr(row, "is_cfos", False),
                       "pomc_signal_score": pomc_score, "cfos_signal_score": cfos_score,
                       "pomc_neighbors_in_inset": neighbor_count,
                       "nearest_pomc_centroid_px": nearest_other})
    if not scored:
        raise ValueError("No POMC+/c-FOS+ nucleus is far enough from image boundaries for the inset.")
    table = pd.DataFrame(scored).sort_values(["pomc_signal_score", "nucleus_label"])
    table["pomc_signal_rank"] = table["pomc_signal_score"].rank(method="average", pct=True)
    table["cfos_signal_rank"] = table["cfos_signal_score"].rank(method="average", pct=True)
    table["distance_to_representative_signal"] = (
        (table["pomc_signal_rank"] - 0.75).abs() + (table["cfos_signal_rank"] - 0.75).abs()
    )
    chosen = table.sort_values(["pomc_neighbors_in_inset", "distance_to_representative_signal", "nucleus_label"]).iloc[0].to_dict()
    y, x = int(chosen["centroid_y"]), int(chosen["centroid_x"])
    return chosen, (y-half, y+half, x-half, x+half)


def add_pomc_inset(ax: plt.Axes, dapi: np.ndarray, cfos: np.ndarray, pomc: np.ndarray, bounds: Tuple[int, int, int, int], args: argparse.Namespace) -> None:
    y0, y1, x0, x1 = bounds
    crop = raw_composite(dapi[y0:y1, x0:x1], cfos[y0:y1, x0:x1], pomc[y0:y1, x0:x1])
    # Mark the field the inset was cropped from. Without it the inset is a
    # magnification of nowhere in particular: the crop is a few tens of
    # micrometres inside a field hundreds of micrometres wide, so a reader
    # cannot locate it by eye. Drawn from the same bounds the crop uses, so the
    # box cannot drift from what is displayed. The pixel centres of the crop are
    # x0..x1-1, and imshow puts pixel centres on integers, so the outline sits
    # half a pixel outside them on every side.
    source_box = Rectangle(
        (x0 - 0.5, y0 - 0.5), (x1 - x0), (y1 - y0),
        fill=False, edgecolor="white", linewidth=1.5,
        linestyle=(0, (3.5, 2.0)), zorder=7,
    )
    source_box.set_path_effects(
        [pe.withStroke(linewidth=3.2, foreground="#0b0b0b", alpha=0.9)]
    )
    ax.add_patch(source_box)
    iax = inset_axes(ax, width="15.5%", height="39.0%", loc="lower left", borderpad=0.55)
    iax.imshow(crop, interpolation="nearest")
    for spine in iax.spines.values():
        spine.set_visible(True); spine.set_linewidth(2.4); spine.set_color("white")
    iax.set_xticks([]); iax.set_yticks([])
    # Leaders from the marked field to the inset that carries it.
    for corner_xy, inset_xy in (((x0 - 0.5, y0 - 0.5), (1.0, 1.0)),
                                ((x0 - 0.5, y1 - 0.5), (1.0, 0.0))):
        leader = ConnectionPatch(
            xyA=corner_xy, coordsA=ax.transData,
            xyB=inset_xy, coordsB=iax.transAxes,
            color="white", linewidth=1.1, linestyle=(0, (3.5, 2.0)),
            zorder=6, clip_on=False,
        )
        leader.set_path_effects(
            [pe.withStroke(linewidth=2.6, foreground="#0b0b0b", alpha=0.85)]
        )
        ax.add_artist(leader)
    add_scalebar(iax, crop.shape[:2], args.um_per_px, args.inset_scalebar_um, font_size=7.5, inset=True)


def iqr_outliers(values: np.ndarray, k: float = 1.5) -> np.ndarray:
    x = np.asarray(values, dtype=float); out = np.zeros(x.shape, dtype=bool); finite = np.isfinite(x)
    if finite.sum() < 4:
        return out
    q1, q3 = np.percentile(x[finite], [25, 75]); iqr = q3 - q1
    if np.isfinite(iqr) and iqr > 0:
        out[finite] = (x[finite] < q1 - k * iqr) | (x[finite] > q3 + k * iqr)
    return out


def safe_anova(groups: Iterable[np.ndarray]) -> float:
    """Ordinary equal-variance one-way ANOVA p value on independent units."""
    usable = [np.asarray(group, dtype=float) for group in groups]
    usable = [group[np.isfinite(group)] for group in usable]
    usable = [group for group in usable if group.size >= 2]
    if len(usable) < 2:
        return float("nan")
    try:
        return float(f_oneway(*usable).pvalue)
    except Exception:
        return 1.0 if np.unique(np.concatenate(usable)).size == 1 else float("nan")


def safe_mwu(a: np.ndarray, b: np.ndarray, max_partitions: int = 1_000_000) -> float:
    """Exact two-sided label-permutation p for the Mann-Whitney rank statistic."""
    x = np.asarray(a, dtype=float); y = np.asarray(b, dtype=float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return float("nan")
    pooled = np.concatenate((x, y))
    n_first = int(x.size)
    total_partitions = math.comb(int(pooled.size), n_first)
    if total_partitions > max_partitions:
        return float(mannwhitneyu(x, y, alternative="two-sided").pvalue)
    ranks = rankdata(pooled, method="average")
    expected_rank_sum = n_first * (pooled.size + 1.0) / 2.0
    observed = abs(float(np.sum(ranks[:n_first])) - expected_rank_sum)
    extreme = sum(
        abs(float(np.sum(ranks[list(selection)])) - expected_rank_sum) >= observed - 1e-12
        for selection in combinations(range(int(pooled.size)), n_first)
    )
    return float(extreme / total_partitions)


def cliffs_delta(second: np.ndarray, first: np.ndarray) -> float:
    """Cliff's delta for the displayed direction second minus first."""
    left = np.asarray(second, dtype=float); right = np.asarray(first, dtype=float)
    left = left[np.isfinite(left)]; right = right[np.isfinite(right)]
    if not left.size or not right.size:
        return float("nan")
    differences = left[:, None] - right[None, :]
    return float((np.count_nonzero(differences > 0) - np.count_nonzero(differences < 0)) / differences.size)


def fmt_p(value: float) -> str:
    if not np.isfinite(value): return "NA"
    return f"{value:.1e}" if value < 0.001 else f"{value:.3g}"


def fmt_metric(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.3f}".rstrip("0").rstrip(".")


SENSITIVITY_ENDPOINTS = {
    "cfos_over_dapi": {
        "long_endpoint": "total_cfos_activation_fraction",
        "label": "c-FOS/DAPI",
        "numerator": "c-FOS-positive nuclei",
        "denominator": "DAPI nuclei",
    },
    "cfos_pomc_over_pomc": {
        "long_endpoint": "primary_pomc_activation_fraction",
        "label": "c-FOS∧POMC/POMC",
        "numerator": "c-FOS∧POMC-positive cells",
        "denominator": "POMC-positive cells",
    },
}


def anova_details(groups: Mapping[str, np.ndarray]) -> Dict[str, object]:
    """Ordinary equal-variance one-way ANOVA details."""
    arrays = [np.asarray(groups[condition], dtype=float) for condition in CONDITION_ORDER]
    arrays = [values[np.isfinite(values)] for values in arrays]
    if any(values.size < 2 for values in arrays):
        return {"statistic": np.nan, "p_value_raw": np.nan,
                "df_between": np.nan, "df_within": np.nan,
                "extreme_labelings": np.nan, "enumerated_labelings": np.nan}
    result = f_oneway(*arrays)
    return {"statistic": float(result.statistic), "p_value_raw": float(result.pvalue),
            "df_between": int(len(arrays) - 1),
            "df_within": int(sum(len(values) for values in arrays) - len(arrays)),
            "extreme_labelings": np.nan, "enumerated_labelings": np.nan}


def exact_mwu_details(first: np.ndarray, second: np.ndarray) -> Dict[str, object]:
    """Exhaustive unadjusted two-sided rank-sum/Mann–Whitney label permutation."""
    left = np.asarray(first, dtype=float); right = np.asarray(second, dtype=float)
    left = left[np.isfinite(left)]; right = right[np.isfinite(right)]
    if left.size < 2 or right.size < 2:
        return {"statistic": np.nan, "p_value_raw": np.nan,
                "extreme_labelings": 0, "enumerated_labelings": 0}
    pooled = np.concatenate((left, right)); ranks = rankdata(pooled, method="average")
    n_first = int(left.size)
    expected_rank_sum = n_first * (pooled.size + 1.0) / 2.0
    observed_distance = abs(float(np.sum(ranks[:n_first])) - expected_rank_sum)
    total_partitions = math.comb(int(pooled.size), n_first)
    extreme = sum(
        abs(float(np.sum(ranks[list(selection)])) - expected_rank_sum)
        >= observed_distance - 1e-12
        for selection in combinations(range(int(pooled.size)), n_first)
    )
    statistic = float(mannwhitneyu(left, right, alternative="two-sided").statistic)
    return {"statistic": statistic, "p_value_raw": float(extreme / total_partitions),
            "extreme_labelings": int(extreme),
            "enumerated_labelings": int(total_partitions)}


def binomial_deviance_statistic(
    successes: np.ndarray, totals: np.ndarray, labels: np.ndarray,
) -> float:
    """Likelihood-ratio deviance for condition-specific versus common binomial rates."""
    successes = np.asarray(successes, dtype=float)
    totals = np.asarray(totals, dtype=float)
    labels = np.asarray(labels, dtype=int)
    total_successes = float(np.sum(successes)); total_trials = float(np.sum(totals))
    if total_trials <= 0:
        return float("nan")
    common_rate = total_successes / total_trials; statistic = 0.0
    for label in np.unique(labels):
        selected = labels == label
        group_successes = float(np.sum(successes[selected]))
        group_trials = float(np.sum(totals[selected]))
        if group_trials <= 0:
            return float("nan")
        group_rate = group_successes / group_trials
        failures = group_trials - group_successes
        if group_successes > 0:
            statistic += 2.0 * group_successes * math.log(group_rate / common_rate)
        if failures > 0:
            statistic += 2.0 * failures * math.log((1.0 - group_rate) / (1.0 - common_rate))
    return float(statistic)


def exact_binomial_deviance_details(
    successes: np.ndarray, totals: np.ndarray, group_sizes: Sequence[int],
) -> Dict[str, object]:
    """Exhaustively permute whole-cage labels; no Monte-Carlo plus-one is used."""
    successes = np.asarray(successes, dtype=float); totals = np.asarray(totals, dtype=float)
    sizes = [int(value) for value in group_sizes]
    if len(sizes) not in {2, 3} or any(value < 2 for value in sizes):
        return {"statistic": np.nan, "p_value_raw": np.nan,
                "extreme_labelings": 0, "enumerated_labelings": 0}
    observed_labels = np.concatenate([
        np.full(size, index, dtype=int) for index, size in enumerate(sizes)
    ])
    observed = binomial_deviance_statistic(successes, totals, observed_labels)
    indices = tuple(range(int(successes.size)))
    total_partitions = math.comb(len(indices), sizes[0])
    if len(sizes) == 3:
        total_partitions *= math.comb(len(indices) - sizes[0], sizes[1])
    extreme = 0; enumerated = 0
    for first in combinations(indices, sizes[0]):
        first_set = set(first)
        remainder = tuple(index for index in indices if index not in first_set)
        second_selections = combinations(remainder, sizes[1]) if len(sizes) == 3 else [remainder]
        for second in second_selections:
            labels = np.full(len(indices), 2 if len(sizes) == 3 else 1, dtype=int)
            labels[list(first)] = 0; labels[list(second)] = 1
            statistic = binomial_deviance_statistic(successes, totals, labels)
            extreme += int(statistic >= observed - 1e-12); enumerated += 1
    if enumerated != total_partitions:
        raise RuntimeError(f"Enumerated {enumerated} cage labelings; expected {total_partitions}")
    return {"statistic": observed, "p_value_raw": float(extreme / enumerated),
            "extreme_labelings": int(extreme), "enumerated_labelings": int(enumerated)}


def compute_sensitivity_statistics(
    definitive_values: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Cage-mean ANOVA/MWU and count-aware exact sensitivity analyses."""
    required = {
        "animal_id", "region", "cond", "cage", "endpoint", "raw_value",
        "numerator", "denominator", "valid_for_inference", "exclusion_reason",
    }
    missing = sorted(required.difference(definitive_values.columns))
    if missing:
        raise ValueError(f"Definitive endpoint table lacks sensitivity columns: {missing}")
    long_endpoints = [str(spec["long_endpoint"]) for spec in SENSITIVITY_ENDPOINTS.values()]
    selected = definitive_values[
        definitive_values["endpoint"].astype(str).isin(long_endpoints)
        & definitive_values["region"].astype(str).str.upper().isin({"ARC", "ME"})
    ].copy()
    for column in ["raw_value", "numerator", "denominator"]:
        selected[column] = pd.to_numeric(selected[column], errors="coerce")
    selected["valid_for_inference"] = as_bool(selected["valid_for_inference"])
    selected["positive_finite_denominator"] = (
        np.isfinite(selected["denominator"]) & (selected["denominator"] > 0)
    )
    selected["low_denominator_lt_10"] = (
        selected["positive_finite_denominator"] & (selected["denominator"] < 10)
    )
    selected["included_in_ratio_analyses"] = (
        selected["valid_for_inference"] & selected["positive_finite_denominator"]
        & np.isfinite(selected["raw_value"])
    )
    selected["included_in_count_analyses"] = (
        selected["valid_for_inference"] & selected["positive_finite_denominator"]
        & np.isfinite(selected["numerator"]) & (selected["numerator"] >= 0)
        & (selected["numerator"] <= selected["denominator"])
    )
    if selected.duplicated(["endpoint", "region", "animal_id"], keep=False).any():
        raise ValueError("Definitive endpoint table has duplicate endpoint/region/animal rows")
    reverse_endpoint = {
        str(spec["long_endpoint"]): endpoint for endpoint, spec in SENSITIVITY_ENDPOINTS.items()
    }
    selected["display_endpoint"] = selected["endpoint"].map(reverse_endpoint)
    denominator_audit = selected[[
        "display_endpoint", "endpoint", "region", "animal_id", "cond", "cage",
        "numerator", "denominator", "raw_value", "valid_for_inference",
        "positive_finite_denominator", "low_denominator_lt_10",
        "included_in_ratio_analyses", "included_in_count_analyses", "exclusion_reason",
    ]].sort_values(
        ["display_endpoint", "region", "cond", "cage", "animal_id"], kind="mergesort"
    ).reset_index(drop=True)

    cage_rows: List[Dict[str, object]] = []
    stats_rows: List[Dict[str, object]] = []
    for display_endpoint, endpoint_spec in SENSITIVITY_ENDPOINTS.items():
        long_endpoint = str(endpoint_spec["long_endpoint"])
        for region in ["ME", "ARC"]:
            subset = selected[
                selected["endpoint"].astype(str).eq(long_endpoint)
                & selected["region"].astype(str).str.upper().eq(region)
            ].copy()
            for cage_record in subset[["cage", "cond"]].drop_duplicates().itertuples(index=False):
                cage_subset = subset[
                    subset["cage"].astype(str).eq(str(cage_record.cage))
                    & subset["cond"].astype(str).eq(str(cage_record.cond))
                ]
                ratio_subset = cage_subset[cage_subset["included_in_ratio_analyses"]]
                count_subset = cage_subset[cage_subset["included_in_count_analyses"]]
                sum_numerator = float(count_subset["numerator"].sum()) if not count_subset.empty else 0.0
                sum_denominator = float(count_subset["denominator"].sum()) if not count_subset.empty else 0.0
                cage_rows.append({
                    "endpoint": display_endpoint, "source_endpoint": long_endpoint,
                    "region": region, "cage": str(cage_record.cage),
                    "condition": str(cage_record.cond),
                    "n_animals_total": int(len(cage_subset)),
                    "n_animals_ratio": int(len(ratio_subset)),
                    "n_animals_count": int(len(count_subset)),
                    "cage_mean_animal_ratio": (
                        float(ratio_subset["raw_value"].mean()) if not ratio_subset.empty else np.nan
                    ),
                    "cage_sum_numerator": sum_numerator,
                    "cage_sum_denominator": sum_denominator,
                    "cage_pooled_ratio": (
                        sum_numerator / sum_denominator if sum_denominator > 0 else np.nan
                    ),
                    "minimum_positive_animal_denominator": (
                        float(count_subset["denominator"].min()) if not count_subset.empty else np.nan
                    ),
                    "contains_low_denominator_lt_10": bool(cage_subset["low_denominator_lt_10"].any()),
                    "count_estimable": bool(sum_denominator > 0),
                })
            cage_table = pd.DataFrame([
                row for row in cage_rows
                if row["endpoint"] == display_endpoint and row["region"] == region
            ])

            rank_table = cage_table[np.isfinite(cage_table["cage_mean_animal_ratio"])].copy()
            rank_groups = {
                condition: rank_table.loc[
                    rank_table["condition"].eq(condition), "cage_mean_animal_ratio"
                ].to_numpy(float)
                for condition in CONDITION_ORDER
            }
            global_rank = anova_details(rank_groups)
            stats_rows.append({
                "analysis_family": "cage_mean_anova", "endpoint": display_endpoint,
                "source_endpoint": long_endpoint, "region": region,
                "test": "One-way ANOVA", "group1": "all", "group2": "",
                **global_rank,
                "effect_size": np.nan,
                "effect_size_definition": "ordinary one-way ANOVA of arithmetic cage means",
                "experimental_unit": "cage", "no_plus_one": np.nan,
            })
            rank_pairs: List[Tuple[str, str, Dict[str, object]]] = []
            for first_index, first in enumerate(CONDITION_ORDER):
                for second in CONDITION_ORDER[first_index + 1:]:
                    rank_pairs.append((first, second, exact_mwu_details(
                        rank_groups[first], rank_groups[second]
                    )))
            for first, second, details in rank_pairs:
                stats_rows.append({
                    "analysis_family": "cage_mean_anova", "endpoint": display_endpoint,
                    "source_endpoint": long_endpoint, "region": region,
                    "test": "Mann-Whitney U", "group1": first, "group2": second,
                    **details,
                    "effect_size": cliffs_delta(rank_groups[second], rank_groups[first]),
                    "effect_size_definition": f"Cliff's delta ({second} minus {first})",
                    "permutation_unit": "cage", "no_plus_one": True,
                })

            count_table = cage_table[cage_table["count_estimable"]].copy()
            count_table["condition_order"] = count_table["condition"].map(
                {condition: index for index, condition in enumerate(CONDITION_ORDER)}
            )
            count_table = count_table.sort_values(["condition_order", "cage"], kind="mergesort")
            count_sizes = [int(count_table["condition"].eq(condition).sum()) for condition in CONDITION_ORDER]
            global_count = exact_binomial_deviance_details(
                count_table["cage_sum_numerator"].to_numpy(float),
                count_table["cage_sum_denominator"].to_numpy(float), count_sizes,
            )
            stats_rows.append({
                "analysis_family": "cage_count_binomial_deviance",
                "endpoint": display_endpoint, "source_endpoint": long_endpoint,
                "region": region, "test": "binomial deviance", "group1": "all",
                "group2": "", **global_count,
                "effect_size": np.nan,
                "effect_size_definition": "global likelihood-ratio deviance",
                "permutation_unit": "cage", "no_plus_one": True,
            })
            count_pairs: List[Tuple[str, str, Dict[str, object]]] = []
            for first_index, first in enumerate(CONDITION_ORDER):
                for second in CONDITION_ORDER[first_index + 1:]:
                    pair_table = count_table[count_table["condition"].isin([first, second])].copy()
                    pair_table["condition_order"] = pair_table["condition"].map({first: 0, second: 1})
                    pair_table = pair_table.sort_values(["condition_order", "cage"], kind="mergesort")
                    pair_sizes = [int(pair_table["condition"].eq(condition).sum()) for condition in [first, second]]
                    details = exact_binomial_deviance_details(
                        pair_table["cage_sum_numerator"].to_numpy(float),
                        pair_table["cage_sum_denominator"].to_numpy(float), pair_sizes,
                    )
                    count_pairs.append((first, second, details))
            for first, second, details in count_pairs:
                first_rows = count_table[count_table["condition"].eq(first)]
                second_rows = count_table[count_table["condition"].eq(second)]
                first_rate = float(first_rows["cage_sum_numerator"].sum() / first_rows["cage_sum_denominator"].sum())
                second_rate = float(second_rows["cage_sum_numerator"].sum() / second_rows["cage_sum_denominator"].sum())
                stats_rows.append({
                    "analysis_family": "cage_count_binomial_deviance",
                    "endpoint": display_endpoint, "source_endpoint": long_endpoint,
                    "region": region, "test": "binomial deviance", "group1": first,
                    "group2": second, **details,
                    "effect_size": second_rate - first_rate,
                    "effect_size_definition": f"pooled proportion difference ({second} minus {first})",
                    "permutation_unit": "cage", "no_plus_one": True,
                })

    cage_values = pd.DataFrame(cage_rows).sort_values(
        ["endpoint", "region", "condition", "cage"], kind="mergesort"
    ).reset_index(drop=True)
    sensitivity = pd.DataFrame(stats_rows)
    return sensitivity, cage_values, denominator_audit


def endpoint_column(values: pd.DataFrame, base: str) -> str:
    return f"plot_{base}" if f"plot_{base}" in values.columns else base


def plot_region_endpoint(ax: plt.Axes, data: pd.DataFrame, region: str, column: str, ylabel: str,
                         seed: int, stats: List[Dict[str, object]], endpoint: str,
                         font_size: float = 10.0) -> None:
    sub = data[data["region"].astype(str).str.upper().eq(region)].copy()
    groups: Dict[str, np.ndarray] = {}
    for cond in CONDITION_ORDER:
        vals = pd.to_numeric(sub.loc[sub["cond"].astype(str).eq(cond), column], errors="coerce").to_numpy(float)
        groups[cond] = vals[np.isfinite(vals)]
    x = np.arange(3, dtype=float)
    means = [float(np.mean(groups[c])) if groups[c].size else np.nan for c in CONDITION_ORDER]
    sds = [float(np.std(groups[c], ddof=1)) if groups[c].size >= 2 else 0.0 for c in CONDITION_ORDER]
    ax.bar(x, means, yerr=sds, width=0.52, color=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
           edgecolor="black", linewidth=1.7, capsize=4, error_kw={"elinewidth": 1.6, "capthick": 1.6}, zorder=2)
    rng = np.random.default_rng(seed); tops: List[float] = []
    for i, cond in enumerate(CONDITION_ORDER):
        vals = groups[cond]; outs = iqr_outliers(vals); jitter = rng.uniform(-0.065, 0.065, vals.size)
        ax.scatter(i + jitter[~outs], vals[~outs], facecolor=CONDITION_COLORS[cond], edgecolor="black", linewidth=1.0, s=32, zorder=4)
        if outs.any():
            ax.scatter(i + jitter[outs], vals[outs], marker="x", color="black", linewidth=1.8, s=42, zorder=5)
        tops.append(max(float(means[i] + sds[i]) if np.isfinite(means[i]) else 0.0, float(np.max(vals)) if vals.size else 0.0))
    # Reserve the upper ~45% for omnibus/pairwise annotations and group Ns so
    # that outlier symbols/error bars never compete with the statistical text.
    y_top = max(max(tops), 1e-4) * 1.85
    ax.set_ylim(0, y_top)
    for i, cond in enumerate(CONDITION_ORDER):
        if groups[cond].size:
            ax.text(i, tops[i] + 0.025*y_top, f"n={groups[cond].size}", ha="center", va="bottom", fontsize=font_size-1.2, fontweight="bold")
    p = safe_anova(groups.values())
    finite_groups = [values for values in groups.values() if values.size >= 2]
    anova_result = f_oneway(*finite_groups) if len(finite_groups) >= 2 else None
    stats.append({
        "endpoint": endpoint, "region": region, "test": "One-way ANOVA",
        "p_value_method": "ordinary equal-variance one-way ANOVA F",
        "group1": "all", "group2": "", "p_value_raw": p,
        "p_value": p,
        "statistic": float(anova_result.statistic) if anova_result is not None else np.nan,
        "df_between": len(finite_groups) - 1 if anova_result is not None else np.nan,
        "df_within": sum(len(values) for values in finite_groups) - len(finite_groups) if anova_result is not None else np.nan,
        "effect_size": np.nan, "effect_size_definition": "ordinary one-way ANOVA",
    })
    pair_records: List[Tuple[str, str, float]] = []
    for i, first in enumerate(CONDITION_ORDER):
        for second in CONDITION_ORDER[i+1:]:
            pair_records.append((first, second, safe_mwu(groups[first], groups[second])))
    pairwise_lines: List[str] = []
    for first, second, p_pair in pair_records:
        delta = cliffs_delta(groups[second], groups[first])
        stats.append({
            "endpoint": endpoint, "region": region, "test": "Mann-Whitney U",
            "p_value_method": "exhaustive independent-label permutation of rank-sum/U",
            "group1": first, "group2": second, "p_value_raw": p_pair,
            "p_value": p_pair,
            "n1": groups[first].size, "n2": groups[second].size,
            "effect_size": delta, "effect_size_definition": f"Cliff's delta ({second} minus {first})",
        })
        pairwise_lines.append(f"MWU {first[0]}–{second[0]} p = {fmt_p(p_pair)}")
    ax.text(
        0.025, 0.975, f"One-way ANOVA p = {fmt_p(p)}", transform=ax.transAxes,
        ha="left", va="top", fontsize=max(8.8, font_size + 0.6), fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
        zorder=11,
    )
    ax.text(
        0.025, 0.865, "\n".join(pairwise_lines), transform=ax.transAxes,
        ha="left", va="top", fontsize=max(6.4, font_size - 2.0), fontweight="bold",
        linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 1.2},
        zorder=10,
    )
    ax.set_title(region, fontsize=font_size+1.5, fontweight="bold", pad=4)
    ax.set_xticks(x); ax.set_xticklabels(CONDITION_ORDER, rotation=24, ha="right", fontsize=font_size, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=font_size+0.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=font_size-0.5, width=1.5, length=4)
    for spine in ax.spines.values(): spine.set_linewidth(1.6)
    ax.grid(False)


def plot_endpoint_pair(fig: plt.Figure, spec, values: pd.DataFrame, column: str, ylabel: str,
                       seed: int, stats: List[Dict[str, object]], endpoint: str, font_size: float) -> List[plt.Axes]:
    grid = spec.subgridspec(1, 2, wspace=0.34)
    axes = [fig.add_subplot(grid[0, i]) for i in range(2)]
    for i, region in enumerate(["ME", "ARC"]):
        plot_region_endpoint(axes[i], values, region, column, ylabel if i == 0 else "", seed+i, stats, endpoint, font_size)
    return axes


def save_image_panel(path: Path, image: np.ndarray, args: argparse.Namespace, *, label: str,
                     regions: Optional[np.ndarray] = None, sources: Optional[Mapping[str, Path]] = None,
                     inset_bounds: Optional[Tuple[int, int, int, int]] = None,
                     dapi: Optional[np.ndarray] = None, cfos: Optional[np.ndarray] = None,
                     pomc: Optional[np.ndarray] = None) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 3.75), facecolor="white")
    ax.imshow(image, interpolation="nearest")
    if regions is not None:
        if np.any(regions == 1): ax.contour(regions == 1, [0.5], colors=[REGION_COLORS["ARC"]], linewidths=3.2)
        if np.any(regions == 2): ax.contour(regions == 2, [0.5], colors=[REGION_COLORS["ME"]], linewidths=3.2)
        if np.any(regions == 3): ax.contour(regions == 3, [0.5], colors=[REGION_COLORS["VMN"]], linewidths=3.2)
    if sources is not None: draw_wall_geometry(ax, sources)
    if regions is not None: add_region_color_key(ax, font_size=7.2)
    add_scalebar(ax, image.shape[:2], args.um_per_px, args.scalebar_um)
    if inset_bounds is not None and dapi is not None and cfos is not None and pomc is not None: add_pomc_inset(ax, dapi, cfos, pomc, inset_bounds, args)
    if label == "D": add_merge_channel_key(ax, font_size=13.0)
    add_panel_label(ax, label, 22); ax.set_axis_off(); fig.subplots_adjust(0.015, 0.015, 0.99, 0.965)
    fig.savefig(path, bbox_inches="tight", metadata={"Title": f"ARC/ME panel {label}"}); plt.close(fig)


def save_plot_panel(path: Path, values: pd.DataFrame, column: str, ylabel: str, title: str,
                    endpoint: str, label: str, seed: int, stats: List[Dict[str, object]]) -> None:
    fig = plt.figure(figsize=(7.3, 3.65), facecolor="white")
    gs = fig.add_gridspec(1, 1, left=0.11, right=0.985, bottom=0.23, top=0.80)
    plot_endpoint_pair(fig, gs[0, 0], values, column, ylabel, seed, stats, endpoint, 11.0)
    fig.text(0.015, 0.94, label, fontsize=22, fontweight="bold", ha="left", va="top")
    fig.text(0.54, 0.94, title, fontsize=14, fontweight="bold", ha="center", va="top")
    fig.savefig(path, bbox_inches="tight", metadata={"Title": f"ARC/ME panel {label}: {title}"}); plt.close(fig)


def spatial_legend_descriptions(args: argparse.Namespace) -> Dict[str, str]:
    """Legend text for panels I and J, quoting the spatial analysis's own numbers.

    The statistics are read back from exact_spatial_permanova.csv rather than
    recomputed here, so the figure can never state a p-value the analysis did not
    produce.
    """
    panel_i, _panel_j = resolve_spatial_panels(args)
    summary_path = panel_i.parent / "exact_spatial_permanova.csv"
    if not summary_path.is_file():
        raise SystemExit(f"Spatial PERMANOVA summary not found beside the panels: {summary_path}")
    summary = pd.read_csv(summary_path).set_index("endpoint")
    shared = (
        "Within each reconstructed section the DAPI centroids define six "
        "covariance-normalized radial shells with approximately equal numbers of DAPI nuclei; an animal's shell "
        "counts are summed over its own sections, so each animal contributes one "
        "six-shell vector. Shell rates are arcsin(sqrt) transformed and compared by an "
        "exhaustive animal-label PERMANOVA on Euclidean distances, enumerating every "
        "allocation that preserves the observed group sizes. The biological animal is the "
        "experimental unit; cells and shells are feature-extraction units only. The same "
        "DAPI centroids define every shell and denominator, so both endpoints are "
        "exploratory and are reported with a Benjamini-Hochberg q across the two of them. "
        "Thin lines are animals, thick lines are condition means, and the heat map gives "
        "the condition mean per shell (W, Water; S, Sucrose; A, Allulose).\n"
    )
    described = {}
    for letter, endpoint, title in (
        ("I", "cfos_occurrence", "c-FOS-positive nuclei / DAPI nuclei"),
        ("J", "double_occurrence", "c-FOS and POMC double-positive nuclei / DAPI nuclei"),
    ):
        row = summary.loc[endpoint]
        excluded = int(row["animals_excluded_incomplete_pomc"])
        exclusion = (
            f" {excluded} animal(s) without a complete POMC channel are excluded from this "
            "endpoint only." if excluded else ""
        )
        described[letter] = (
            f"({letter}) Spatial occurrence of {title} across ARC and ME. "
            f"n = {int(row['n_animals'])} animals "
            f"({int(row['n_water'])} Water, {int(row['n_sucrose'])} Sucrose, "
            f"{int(row['n_allulose'])} Allulose); "
            f"{int(row['enumerated_labelings']):,} exhaustive labelings; "
            f"pseudo-F = {float(row['pseudo_f']):.4f}, R² = {float(row['r_squared']):.4f}, "
            f"exact p = {float(row['p_value_exact']):.4f}, "
            f"BH q = {float(row['bh_q_value_two_endpoint_family']):.4f}.{exclusion} "
            + shared
        )
    return described


def write_legends(outdir: Path, args: argparse.Namespace, chosen: pd.Series, inset: Mapping[str, object],
                  stats_df: pd.DataFrame, sensitivity_df: pd.DataFrame, values_path: Path,
                  review_provenance: Mapping[str, object],
                  registration: Mapping[str, object],
                  cellpose_provenance: Mapping[str, Mapping[str, object]],
                  panel_f_distribution: Mapping[str, object]) -> None:
    shift = int(registration["_right_shift_y_px"])
    right_origin = int(registration["_right_origin_x_px"])
    right_crop = int(registration["_right_medial_crop_px"])
    overlap = int(registration["_overlap_width_px"])
    gap = int(registration["_gap_px"])
    method_key = str(registration.get("method", "NA"))
    status = str(registration.get("status", "NA"))
    align_method = {
        "sift_ransac_raw_ncc": "DAPI SIFT/RANSAC followed by raw-intensity NCC",
        "direct_ncc": "direct DAPI normalized cross-correlation",
        "anatomical_y_fallback": "DAPI anatomical vertical fallback with image-by-image confirmation",
        "single_half_native": "single native half (no pair registration required)",
    }.get(method_key, method_key.replace("_", " "))
    ncc_text = (
        f"raw NCC {fmt_metric(float(registration['_raw_overlap_ncc']))}, "
        f"gradient NCC {fmt_metric(float(registration['_gradient_overlap_ncc']))}, "
        f"best score {fmt_metric(float(registration['_ncc_best_score']))}, "
        f"second score {fmt_metric(float(registration['_ncc_second_score']))}, "
        f"uniqueness margin {fmt_metric(float(registration['_ncc_uniqueness_margin']))}"
    )
    endpoint_labels = {endpoint: str(spec["label"]) for endpoint, spec in SENSITIVITY_ENDPOINTS.items()}

    def compact_line(table: pd.DataFrame, endpoint: str, region: str, family: str, prefix: str) -> str:
        rows = table[
            table["endpoint"].astype(str).eq(endpoint)
            & table["region"].astype(str).eq(region)
            & table["analysis_family"].astype(str).eq(family)
        ]
        global_rows = rows[rows["group1"].astype(str).eq("all")]
        if len(global_rows) != 1:
            raise ValueError(f"Expected one global {family} row for {endpoint}/{region}")
        global_row = global_rows.iloc[0]
        if family == "cage_count_binomial_deviance":
            pieces = [f"{prefix}: deviance exact p = {fmt_p(float(global_row['p_value_raw']))}"]
            pair_prefix = "pair"
        else:
            pieces = [f"{prefix}: one-way ANOVA p = {fmt_p(float(global_row['p_value_raw']))}"]
            pair_prefix = "MWU"
        for first, second in [("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose")]:
            pair = rows[
                rows["group1"].astype(str).eq(first)
                & rows["group2"].astype(str).eq(second)
            ]
            if len(pair) != 1:
                raise ValueError(f"Expected one {first}/{second} {family} row for {endpoint}/{region}")
            pieces.append(f"{pair_prefix} {first[0]}–{second[0]} p = {fmt_p(float(pair.iloc[0]['p_value_raw']))}")
        return "; ".join(pieces) + "."

    stat_lines: List[str] = []
    stat_blocks: Dict[str, List[str]] = {}
    for endpoint in ["cfos_over_dapi", "cfos_pomc_over_pomc"]:
        endpoint_lines: List[str] = []
        for region in ["ME", "ARC"]:
            label = endpoint_labels[endpoint]
            endpoint_lines.append(f"{label}, {region}:")
            endpoint_lines.append(compact_line(stats_df, endpoint, region, "animal_anova", "animal"))
            endpoint_lines.append(compact_line(sensitivity_df, endpoint, region, "cage_mean_anova", "cage mean"))
            endpoint_lines.append(compact_line(
                sensitivity_df, endpoint, region, "cage_count_binomial_deviance", "cage counts"
            ))
        stat_blocks[endpoint] = endpoint_lines
        stat_lines.extend(endpoint_lines)
    common_alignment = (
        f"Representative section: {chosen['animal_id']} ({chosen['cond']}), source {chosen['sample']}, section {int(chosen['section_index'])}.\n"
        f"Corrected registration {registration['schema_version']}/{registration['algorithm']}: "
        f"status {status}; method {align_method}; right origin x={right_origin} px, "
        f"right y shift={shift} px, right medial crop={right_crop} px, overlap={overlap} px, "
        f"gap={gap} px; {ncc_text}. Registration used DAPI only; the identical integer "
        "placement/crop was applied without resampling to every channel and mask. Panels A-F "
        f"share the exact registered field; DAPI TIFF SHA-256 {registration['_dapi_sha256']}; "
        f"NPY-GFP TIFF SHA-256 {registration['_artifact_sha256s']['npy']}.\n"
    )
    scale_text = f"Main microscopy and cartoon scale bars: {args.scalebar_um} µm. Inset scale bar: {args.inset_scalebar_um:g} µm.\n"
    cfos_yes = "yes" if str(inset["is_cfos"]).strip().lower() in {"true", "1", "yes"} else "no"
    panel_d_description = (
        "(D) Channel-separated segmentation cartoons of the exact Panel E field, in the same "
        "order as Panels A-C: DAPI is blue, c-FOS is magenta, and POMC is amber. The hash-verified left/right raw Cellpose NPY-GFP "
        "instance masks were rebuilt with the audited integer hard seam and no resampling. "
        "POMC instances smaller than the mean DAPI nuclear area of their section and anatomical region were excluded by a fixed size heuristic, applied to all conditions. Original fluorescence, segmentation inputs and accepted anatomy were retained. The c-FOS minimum-area filter (8 px) and POMC exclusion at NPY-GFP overlap >= 0.20 were "
        "reapplied and proved pixel-identical to the registered analysis labels "
        f"(DAPI n={int(cellpose_provenance['DAPI']['processed_instances'])}, "
        f"POMC n={int(cellpose_provenance['POMC']['processed_instances'])}, "
        f"c-FOS n={int(cellpose_provenance['c-FOS']['processed_instances'])}). "
        "The gray tissue boundary is DAPI-derived away from the ventricular corridor. The "
        "upper ventricular lumen is corrected row by row from the medial boundaries of the "
        "two accepted HIL ARC lobes, extended to the cropped image top and stopped before "
        "the first accepted HIL ME row; no automated wall CSV or marker intensity is used. "
        "Cyan ARC and orange ME fills, contours, and labels are the "
        "exact HIL-accepted mask used in D itself; no region was inferred from POMC, "
        "NPY-GFP, or c-FOS.\n"
    )
    channel_descriptions = {
        "A": (
            "(A) Isolated DAPI channel of the registered field, shown in the blue tint it "
            "contributes to the Panel E merge. The field is labeled FR722 S02, Allulose.\n"
        ),
        "B": (
            "(B) Isolated c-FOS channel of the same registered pixels, shown in the magenta "
            "tint it contributes to the Panel E merge. The field is labeled FR722 S02, Allulose.\n"
        ),
        "C": (
            "(C) Isolated POMC channel of the same registered pixels, shown in the amber tint "
            "it contributes to the Panel E merge. Panels A-C are single channels of one "
            "acquisition; their clipped sum over the identical pixels is exactly Panel E, and "
            "no channel was rescaled, shifted, or resampled relative to the others. The field "
            "is labeled FR722 S02, Allulose.\n"
        ),
    }
    panel_e_description = (
        f"(E) Aligned 20× composite of Panels A-C: DAPI is blue, c-FOS is magenta, and POMC is amber. "
        f"Green is reserved for NPY-GFP throughout the paper, which is why POMC is amber here. "
        f"The clipped sum of the exact displayed A/B/C pixels is Panel E. "
        f"The field is named at the upper left. The white-bordered inset at the lower left shows a "
        f"representative isolated POMC-positive/c-FOS-positive {str(inset.get('region', '') or 'accepted-ROI')} cell "
        f"(nucleus label {int(inset['nucleus_label'])}, centroid x/y {int(inset['centroid_x'])}/{int(inset['centroid_y'])} px; "
        f"c-FOS-positive: {cfos_yes}); the dashed box on the field marks exactly the region it magnifies "
        f"and the dashed leaders join the two. The inset contains DAPI, c-FOS, and POMC. "
        f"No NPY-GFP intensity or segmentation layer is present in Panel E.\n"
    )
    panel_f_description = (
        "(F) Registered DAPI/c-FOS/POMC/NPY-GFP microscopy. DAPI remains blue and c-FOS "
        "magenta. NPY-GFP uses the real "
        "registered TIFF in fluorescent green with a deterministic whole-field "
        "median/99.7th-percentile background display window and gamma 0.72, followed by a "
        "monotone screen blend. This true transgenic signal is not gated by the NPY-GFP ROI "
        "mask; that ROI mask is used only in the inset. The real registered POMC "
        "intensity uses exactly the same full-field normalization and amber tint as C/E, "
        "with no additional gain. It is retained strictly inside the exact processed "
        f"POMC ROIs (n={int(panel_f_distribution['pomc_rois'])}; "
        f"{int(panel_f_distribution['pomc_roi_pixels']):,} pixels). No flat ROI paint is used, "
        "and POMC contributes zero signal outside those ROIs. POMC and NPY-GFP are "
        "screen-composited over DAPI/c-FOS to avoid additional clipping of summed channels. "
        "The upper-right inset shows those same POMC microscopy pixels and real NPY-GFP "
        f"intensities within {int(panel_f_distribution['npy_gfp_rois'])} NPY-GFP ROIs "
        f"({int(panel_f_distribution['npy_gfp_roi_pixels']):,} pixels), using the same full-field "
        "display settings and screen blend on black. ROI labels define where signal is "
        "shown; they never replace the measured within-cell intensity variation with a "
        "uniform fill. The inset is 40% of the panel width. Thin white dashed lines delimit the DAPI-supported "
        "tissue exterior and HIL-ARC-derived ventricular lumen. No tissue or ARC/ME fill, marker "
        "ROI outline, ARC/ME contour, halo, dot dilation, smoothing, or resampling is present. "
        f"The {int(panel_f_distribution['overlap_pixels'])} pixels shared by the processed marker "
        "masks retain both color components. These display choices change no source TIFF, ROI "
        "membership, count, overlap call, or quantitative endpoint.\n"
    )
    spatial_descriptions = spatial_legend_descriptions(args)
    master = (
        "Figure ARC/ME.\n\n"
        + common_alignment + scale_text
        + channel_descriptions["A"] + channel_descriptions["B"] + channel_descriptions["C"]
        + panel_d_description
        + panel_e_description
        + panel_f_description
        + "(G) Animal-level c-FOS/DAPI ratios in ME and ARC.\n"
        + "(H) Animal-level c-FOS∧POMC/total-POMC ratios in ME and ARC.\n"
        + spatial_descriptions["I"]
        + spatial_descriptions["J"]
        + "Ring construction is illustrated once in Figure 3G. Corresponding-ring positive and DAPI counts are summed across each animal's sections before calculating percentages.\n"
        + "Bars are mean ± SD; dots are biological animals; crosses mark Tukey-IQR values retained in inference. Animal omnibus values use ordinary one-way ANOVA; W–S/W–A/S–A values use unadjusted exact two-sided MWU. Cage-mean sensitivity uses the same ANOVA/MWU family, while cage-summed numerator/denominator sensitivity uses exact binomial deviance; low positive denominators are retained (W, Water; S, Sucrose; A, Allulose).\n"
        + f"Animal-level values source: {values_path}.\n\nStatistics:\n"
        + "\n".join(stat_lines) + "\n"
    )
    (outdir / f"{args.figure_name}_LEGEND.txt").write_text(master, encoding="utf-8")
    # Compatibility alias used by earlier versions of this analysis package.
    (outdir / f"{args.figure_name}_caption.txt").write_text(master, encoding="utf-8")
    channel_scale_text = f"Single-channel scale bars: {args.scalebar_um} µm.\n"
    panel_text = {
        "A_dapi": common_alignment + channel_scale_text + channel_descriptions["A"],
        "B_cfos": common_alignment + channel_scale_text + channel_descriptions["B"],
        "C_pomc": common_alignment + channel_scale_text + channel_descriptions["C"],
        "D_cellpose_cartoon": common_alignment + scale_text + panel_d_description,
        "E_microscopy": common_alignment + scale_text + panel_e_description,
        "F_microscopy_npy_gfp": common_alignment + scale_text + panel_f_description,
        "G_cfos_dapi": "(G) Animal-level c-FOS/DAPI ratios in ME and ARC. Bars are mean ± SD; dots are animals; crosses mark IQR values retained in inference. Omnibus p values use ordinary one-way ANOVA; pairwise p values use unadjusted exact two-sided MWU.\n" + "\n".join(stat_blocks["cfos_over_dapi"]) + "\n",
        "H_pomc_activation": "(H) Animal-level c-FOS∧POMC/total POMC ratios in ME and ARC. Bars are mean ± SD; dots are animals; crosses mark IQR values retained in inference. This is the prespecified primary POMC-activation endpoint; omnibus p values use ordinary one-way ANOVA and pairwise p values use unadjusted exact two-sided MWU.\n" + "\n".join(stat_blocks["cfos_pomc_over_pomc"]) + "\n",
        "I_spatial_cfos": spatial_descriptions["I"],
        "J_spatial_cfos_pomc": spatial_descriptions["J"],
    }
    if set(panel_text) != {f"{letter}_{suffix}" for letter, suffix in PANEL_FILE_SUFFIXES.items()}:
        raise ValueError("Panel legend keys and panel file suffixes disagree")
    for suffix, text in panel_text.items():
        (outdir / f"Figure_ARC_ME_panel_{suffix}_LEGEND.txt").write_text(text, encoding="utf-8")
    spanish_panels = {
        "A_dapi": (
            "(A) Canal DAPI aislado del campo registrado, mostrado en azul como en la "
            "composición del panel E. Campo representativo FR722 S02, Alulosa."
        ),
        "B_cfos": (
            "(B) Canal c-FOS aislado de los mismos píxeles registrados, mostrado en "
            "magenta como en la composición del panel E."
        ),
        "C_pomc": (
            "(C) Canal POMC aislado de los mismos píxeles, mostrado en ámbar. La suma "
            "recortada de A, B y C sobre coordenadas idénticas produce exactamente E."
        ),
        "D_cellpose_cartoon": (
            "(D) Esquemas de segmentación Cellpose separados por canal para el campo "
            "de E. DAPI, c-FOS y POMC usan las máscaras procesadas verificadas; ARC y "
            "ME provienen únicamente de la máscara HIL aceptada."
        ),
        "E_microscopy": (
            "(E) Composición alineada 20× de DAPI, c-FOS y POMC. El recuadro muestra "
            "una célula POMC-positiva representativa y la caja discontinua identifica "
            "exactamente el campo ampliado."
        ),
        "F_microscopy_npy_gfp": (
            "(F) Microscopía registrada DAPI/c-FOS/POMC/NPY-GFP. NPY-GFP conserva la "
            "intensidad transgénica completa; POMC conserva su intensidad real solo "
            "dentro de las ROIs aceptadas, con la misma escala y color que C, sin ganancia "
            "adicional. POMC y NPY-GFP usan mezcla de pantalla sobre DAPI/c-FOS. El recuadro "
            "muestra las intensidades microscópicas reales de ambos marcadores dentro de "
            "sus ROIs, con los mismos ajustes del campo completo; las máscaras no aportan "
            "rellenos uniformes de color."
        ),
        "G_cfos_dapi": (
            "(G) Cociente c-FOS/DAPI por animal en ME y ARC."
        ),
        "H_pomc_activation": (
            "(H) Cociente c-FOS∧POMC/POMC total por animal en ME y ARC; es la variable "
            "primaria preespecificada de activación de POMC."
        ),
        "I_spatial_cfos": (
            "(I) Ocurrencia espacial de núcleos c-FOS positivos en seis capas radiales "
            "con aproximadamente igual número de núcleos DAPI, resumida una vez por animal."
        ),
        "J_spatial_cfos_pomc": (
            "(J) Ocurrencia espacial de núcleos doble positivos c-FOS/POMC usando las "
            "mismas capas y denominadores DAPI del panel I. La definición de los anillos se "
            "ilustra en la Figura 3G; los recuentos correspondientes se suman entre "
            "las secciones de cada animal antes de calcular los porcentajes."
        ),
    }
    spanish_shared = (
        "La unidad biológica es el animal. Las barras muestran media ± DE, los puntos "
        "son animales y las cruces señalan valores fuera de las cercas de Tukey, que "
        "se conservaron. Los valores globales usan ANOVA ordinario de una vía y las "
        "comparaciones Agua–Sacarosa, Agua–Alulosa y Sacarosa–Alulosa usan Mann–Whitney "
        "exacta bilateral sin ajuste. Las sensibilidades usan medias por jaula y una "
        "desviancia binomial exacta sobre recuentos sumados por jaula. Los paneles I/J "
        "usan PERMANOVA exacta por permutación de etiquetas animales y corrección BH "
        "sobre las dos variables. Todas las regiones ARC/ME/VMN cuantificadas proceden "
        "de la anotación HIL nativa completa; células y secciones no son réplicas "
        "biológicas. Se excluyeron máscaras POMC menores que el área media de DAPI en la misma sección y región, mediante la misma regla en todos los tratamientos. Se conservaron las imágenes, máscaras originales y anatomía aceptada."
    )
    spanish_master = (
        "Figura 4. POMC/c-FOS en ARC y ME con anatomía HIL aceptada.\n\n"
        + "\n".join(spanish_panels.values())
        + "\n\n" + spanish_shared + "\n"
    )
    (outdir / f"{args.figure_name}_LEGEND_spanish.txt").write_text(
        spanish_master, encoding="utf-8"
    )
    for suffix, text in spanish_panels.items():
        (outdir / f"Figure_ARC_ME_panel_{suffix}_LEGEND_spanish.txt").write_text(
            "Figura 4. POMC/c-FOS en ARC y ME.\n\n"
            + text + "\n\n" + spanish_shared + "\n",
            encoding="utf-8",
        )


def make_figure(args: argparse.Namespace) -> Tuple[Path, ...]:
    root = args.analysis_dir.expanduser().resolve()
    manual_root = args.manual_region_dir.expanduser().resolve()
    outdir_candidate = args.outdir.expanduser().resolve()
    summary_path = root / "per_reconstructed_section_summary.csv"
    cell_path = root / "per_cell_reconstructed_measurements.csv"
    for path in [summary_path, cell_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    # Every scientific/provenance check precedes output replacement, so a stale
    # automated run cannot erase the last valid publication figure under --force.
    if not bool(args.require_manual_regions):
        raise ValueError("Figure 4 publication rendering requires complete 76/76 native HIL")
    summary = pd.read_csv(summary_path, low_memory=False)
    quantified = validate_analysis_summary(
        summary, require_manual_regions=bool(args.require_manual_regions),
    )
    all_manual = bool(
        "region_source" in quantified.columns
        and quantified["region_source"].fillna("").astype(str)
            .eq("manual_region_tiff").all()
    )
    if not all_manual:
        raise ValueError(
            "Figure 4 publication rendering refuses automated or mixed ARC/ME/VMN masks"
        )
    if bool(args.require_manual_regions) and not manual_root.is_dir():
        raise FileNotFoundError(
            f"Mandatory POMC HIL root is missing: {manual_root}"
        )
    registrations, newest_registration_mtime_ns = validate_registered_reconstructions(
        quantified, root,
    )
    values, definitive_values, section_region_path = recompute_registered_statistical_inputs(
        root=root, quantified=quantified,
        newest_registration_mtime_ns=newest_registration_mtime_ns,
        required_fresh_paths=[summary_path, cell_path],
    )
    chosen, folder, candidates = select_representative(
        quantified, root, args.representative_sample, args.representative_section,
    )
    chosen_registration_key = (
        str(chosen["animal_id"]), str(chosen["sample"]),
        int(chosen["section_index"]), clean_tile(chosen.get("tile", "")),
    )
    registration = registrations[chosen_registration_key]
    sources = microscopy_sources(folder)
    dapi = read_gray(sources["dapi"])
    cfos = read_gray(sources["cfos"])
    pomc = read_gray(sources["pomc"])
    npy_gfp = read_gray(sources["npy"])
    labels = read_gray(sources["dapi_labels"]).astype(np.int32)
    pomc_roi_labels = read_gray(sources["pomc_labels"]).astype(np.int32)
    npy_gfp_roi_labels = read_gray(sources["npy_labels"]).astype(np.int32)
    if not (
        dapi.shape == cfos.shape == pomc.shape == npy_gfp.shape == labels.shape
        == pomc_roi_labels.shape == npy_gfp_roi_labels.shape
    ):
        raise ValueError(
            "Corrected registered DAPI/c-FOS/POMC/NPY-GFP/label and marker-ROI "
            "planes do not "
            "share one geometry: "
            f"{dapi.shape}, {cfos.shape}, {pomc.shape}, {npy_gfp.shape}, "
            f"{labels.shape}, {pomc_roi_labels.shape}, {npy_gfp_roi_labels.shape}"
        )
    if sha256_file(sources["dapi"]) != registration["_dapi_sha256"]:
        raise ValueError("Panel D DAPI differs from the audited corrected registration")
    if sha256_file(sources["npy"]) != registration["_artifact_sha256s"].get("npy"):
        raise ValueError("Panel F NPY-GFP differs from the audited corrected registration")
    for key, path in (("pomc_labels", sources["pomc_labels"]),
                      ("npy_labels", sources["npy_labels"])):
        if sha256_file(path) != registration["_artifact_sha256s"].get(key):
            raise ValueError(f"Panel F {key} differ from the audited corrected registration")
    regions, review_provenance = load_current_pomc_regions(
        chosen, folder, manual_root,
        require_manual_regions=bool(args.require_manual_regions),
    )
    if regions.shape != dapi.shape:
        raise ValueError("Panel E accepted ROI mask does not share Panel D registered bounds")
    pomc_roi_labels, pomc_size_qc, pomc_size_path = load_size_filtered_pomc(sources)
    if str(review_provenance["source_dapi_sha256"]) != str(registration["_dapi_sha256"]):
        raise ValueError("Panel E HIL receipt does not bind Panel D registered DAPI hash")
    displayed_regions = tuple(
        region for region, code in (("ARC", 1), ("ME", 2), ("VMN", 3))
        if np.any(regions == code)
    )
    if bool(review_provenance.get("human_reviewed", False)):
        receipt_regions = tuple(
            region for region in ("ARC", "ME", "VMN")
            if review_provenance.get("region_status", {}).get(region) == "drawn"
        )
        if displayed_regions != receipt_regions:
            raise ValueError(
                "Representative accepted-mask labels differ from receipt region_status: "
                f"mask={displayed_regions}, receipt={receipt_regions}"
            )
    cells = pd.read_csv(cell_path, low_memory=False)
    selected = cells[
        cells["sample"].astype(str).eq(str(chosen["sample"]))
        & pd.to_numeric(cells["section_index"], errors="coerce").eq(
            int(chosen["section_index"])
        )
    ].copy()
    chosen_tile = clean_tile(chosen.get("tile", ""))
    if chosen_tile and "tile" in selected.columns:
        selected = selected[selected["tile"].map(clean_tile).eq(chosen_tile)].copy()
    if selected.empty:
        raise ValueError("No per-cell rows match the representative section.")
    sensitivity_df, sensitivity_cages, denominator_audit = compute_sensitivity_statistics(
        definitive_values
    )
    composite = raw_composite(dapi, cfos, pomc)
    accepted_overlay = accepted_segmentation_overlay(dapi, regions)
    inset_field_px = int(round(args.inset_field_um / args.um_per_px))
    inset_info, inset_bounds = select_pomc_inset(
        selected, labels, cfos, pomc, inset_field_px,
    )
    cellpose_masks, cellpose_provenance = load_cellpose_cartoon_masks(
        registration, folder, root, sources,
    )
    if not np.array_equal(pomc_roi_labels > 0, cellpose_masks["POMC"]):
        raise ValueError("Panel F POMC ROI TIFF differs from raw-mask reconstruction")
    if not bool(cellpose_provenance["NPY"].get("matches_registered_analysis_labels", False)):
        raise ValueError("Panel F NPY-GFP ROI TIFF lacks raw-mask identity proof")
    composite_npy_gfp = composite_with_roi_masked_pomc_signal_and_npy_gfp(
        dapi, cfos, pomc, npy_gfp, pomc_roi_labels,
    )
    hil_ventricle_lumen, hil_ventricle_info = hil_arc_ventricle_lumen(regions)
    tissue_envelope = cellpose_tissue_envelope(
        cellpose_masks["DAPI"], regions, hil_ventricle_lumen,
    )

    outdir = _prepare_generated_output_dir(
        outdir_candidate,
        force=bool(args.force),
        protected_paths=[
            root, summary_path, cell_path, section_region_path, manual_root,
            PROJECT_ROOT, PROJECT_ROOT / "scripts",
            PROJECT_ROOT / "20x", PROJECT_ROOT / "20x_sanitized",
        ],
        preserved_names=("README.txt",),
    )
    panels_dir = outdir / "panels"
    legends_dir = outdir / "legends"
    source_data_dir = outdir / "source_data"
    provenance_dir = outdir / "provenance"
    for directory in (panels_dir, legends_dir, source_data_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=True)
    values_path = source_data_dir / "registered_recomputed_animal_values.csv"
    definitive_path = source_data_dir / "registered_recomputed_endpoint_values.csv"
    values.to_csv(values_path, index=False)
    definitive_values.to_csv(definitive_path, index=False)
    registration_audit_path = provenance_dir / "all_section_stitch_registration_audit.csv"
    registration_audit_rows: List[Dict[str, object]] = []
    for identity in sorted(registrations):
        record = registrations[identity]
        row: Dict[str, object] = {
            "animal_id": identity[0], "sample": identity[1],
            "section_index": identity[2], "tile": identity[3],
            "schema_version": record["schema_version"],
            "algorithm": record["algorithm"], "status": record["status"],
            "method": record["method"],
            "requires_manual_review": record["requires_manual_review"],
            "hil_receipt_bound_to_registration": record["_registration_receipt_bound"],
            "manual_registration_confirmed": record["_manual_registration_confirmed"],
            "right_origin_x_px": record["_right_origin_x_px"],
            "right_shift_y_px": record["_right_shift_y_px"],
            "right_medial_crop_px": record["_right_medial_crop_px"],
            "overlap_width_px": record["_overlap_width_px"],
            "gap_px": record["_gap_px"],
            "raw_overlap_ncc": record["_raw_overlap_ncc"],
            "gradient_overlap_ncc": record["_gradient_overlap_ncc"],
            "ncc_best_score": record["_ncc_best_score"],
            "ncc_second_score": record["_ncc_second_score"],
            "ncc_uniqueness_margin": record["_ncc_uniqueness_margin"],
            "registration_json": record["_registration_path"],
            "registration_json_sha256": record["_registration_sha256"],
            "registered_dapi_sha256": record["_dapi_sha256"],
        }
        row.update({
            f"registered_{key}_sha256": value
            for key, value in record["_artifact_sha256s"].items()
        })
        row.update({
            f"input_{side}_dapi_sha256": value
            for side, value in record["_input_dapi_sha256s"].items()
        })
        registration_audit_rows.append(row)
    pd.DataFrame(registration_audit_rows).to_csv(registration_audit_path, index=False)
    # savefig.dpi defaults to "figure", which is 100, and only the master PNG
    # passed dpi= explicitly. Every PDF here therefore embedded its microscopy at
    # 100 ppi while the PNG got 600: pdfimages showed 588x237 rasters for a
    # 720x1788 source. Setting it once covers the master and every isolated panel,
    # and an explicit dpi= at a call site still wins.
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans",
                         "axes.linewidth": 1.6, "savefig.dpi": int(args.dpi)})

    panel_paths = {
        letter: panels_dir / f"Figure_ARC_ME_panel_{letter}_{suffix}.pdf"
        for letter, suffix in PANEL_FILE_SUFFIXES.items()
    }
    spanish_panel_paths = {
        letter: panels_dir / f"Figure_ARC_ME_panel_{letter}_{suffix}_spanish.pdf"
        for letter, suffix in PANEL_FILE_SUFFIXES.items()
    }
    d_col = endpoint_column(values, "cfos_over_dapi")
    e_col = endpoint_column(values, "cfos_pomc_over_pomc")

    spatial_i, spatial_j = resolve_spatial_panels(args)
    spanish_spatial_i, spanish_spatial_j = spanish_spatial_variants(
        (spatial_i, spatial_j)
    )
    fig = plt.figure(figsize=(13.2, 15.1), facecolor="white")
    outer = fig.add_gridspec(
        # Row heights track each row's native 2.48:1 registered aspect, so the
        # added single-channel row does not open a letterboxed band above it.
        # The cartoon row now sits second, directly under the channels it
        # depicts, and the two spatial occurrence panels close the figure.
        # The spatial row's ratio is set so its two 1.97:1 panels fill the column
        # width rather than being height-limited and shrinking away from their
        # panel letters.
        5, 2, height_ratios=[0.54, 0.62, 0.82, 0.76, 1.02],
        left=0.035, right=0.992, bottom=0.036, top=0.984,
        wspace=0.15, hspace=0.30,
    )
    ax_a, ax_b = fig.add_subplot(outer[2, 0]), fig.add_subplot(outer[2, 1])
    registered_xlim = (-0.5, float(dapi.shape[1]) - 0.5)
    registered_ylim = (float(dapi.shape[0]) - 0.5, -0.5)
    channel_grid = outer[0, :].subgridspec(1, 3, wspace=0.055)
    channel_axes = [fig.add_subplot(channel_grid[0, index]) for index in range(3)]
    for ax, label, marker, channel in zip(
        channel_axes, "ABC", COMPOSITE_CHANNEL_ORDER, (dapi, cfos, pomc),
    ):
        draw_single_channel_panel(ax, channel, marker, args)
        add_field_identity(ax, chosen, font_size=7.2)
        add_panel_label(ax, label)
        ax.set_xlim(*registered_xlim)
        ax.set_ylim(*registered_ylim)

    ax_a.imshow(composite, interpolation="nearest")
    add_scalebar(ax_a, dapi.shape, args.um_per_px, args.scalebar_um)
    add_pomc_inset(ax_a, dapi, cfos, pomc, inset_bounds, args)
    add_merge_channel_key(ax_a, font_size=11.6)
    add_field_identity(ax_a, chosen)
    add_panel_label(ax_a, "E")
    ax_a.set_xlim(*registered_xlim); ax_a.set_ylim(*registered_ylim); ax_a.set_axis_off()

    ax_b.imshow(composite_npy_gfp, interpolation="nearest")
    add_scalebar(ax_b, dapi.shape, args.um_per_px, args.scalebar_um)
    panel_f_distribution_info = add_marker_distribution_inset(
        ax_b, pomc, npy_gfp, pomc_roi_labels, npy_gfp_roi_labels,
        tissue_envelope, hil_ventricle_lumen,
        args.um_per_px, args.scalebar_um,
    )
    add_merge_channel_key(
        ax_b, font_size=11.6, include_npy_gfp=True, pomc_rois_only=True,
    )
    add_field_identity(ax_b, chosen)
    add_panel_label(ax_b, "F")
    ax_b.set_xlim(*registered_xlim)
    ax_b.set_ylim(*registered_ylim)
    ax_b.set_axis_off()

    cartoon_grid = outer[1, :].subgridspec(1, 3, wspace=0.055)
    c_axes = [fig.add_subplot(cartoon_grid[0, index]) for index in range(3)]
    for ax, marker in zip(c_axes, CARTOON_MARKER_ORDER):
        draw_cellpose_cartoon(
            ax, marker, cellpose_masks[marker], tissue_envelope, regions, args,
        )
    c_pos = outer[1, :].get_position(fig)
    c_label = fig.text(c_pos.x0 - 0.012, c_pos.y1 + 0.012, "D", fontsize=20, fontweight="bold", ha="left", va="bottom")
    c_title = fig.text((c_pos.x0 + c_pos.x1) / 2, c_pos.y1 + 0.014, "Cellpose segmentation cartoon of the representative field", fontsize=12.5, fontweight="bold", ha="center", va="bottom")

    stats_multi: List[Dict[str, object]] = []
    d_axes = plot_endpoint_pair(fig, outer[3, 0], values, d_col, "c-FOS/DAPI ratio", args.seed, stats_multi, "cfos_over_dapi", 9.8)
    e_axes = plot_endpoint_pair(fig, outer[3, 1], values, e_col, "c-FOS∧POMC /\nPOMC ratio", args.seed+100, stats_multi, "cfos_pomc_over_pomc", 9.8)
    panel_artists = {
        "A": [channel_axes[0]], "B": [channel_axes[1]], "C": [channel_axes[2]],
        "D": [*c_axes, c_label, c_title], "E": [ax_a], "F": [ax_b],
    }
    # Panels I and J are the exact spatial occurrence analysis Figure 3 runs,
    # applied to these reconstructed ARC/ME sections by
    # Fig4/05_analyze_spatial_distributions.py and rendered here as images so the
    # statistics in the figure are the ones that script enumerated.
    spatial_grid = outer[4, :].subgridspec(1, 2, wspace=0.10)
    s_axes = []
    spatial_image_artists = []
    for column_index, panel_path in enumerate((spatial_i, spatial_j)):
        axis = fig.add_subplot(spatial_grid[0, column_index])
        spatial_image_artists.append(
            axis.imshow(plt.imread(str(panel_path)), interpolation="nearest")
        )
        axis.set_axis_off()
        s_axes.append(axis)
    for letter, axis in zip(("I", "J"), s_axes):
        pos = axis.get_position(fig)
        panel_label = fig.text(pos.x0 - 0.012, pos.y1 - 0.004, letter, fontsize=20,
                               fontweight="bold", ha="left", va="bottom")
        panel_artists[letter] = [axis, panel_label]
    for spec, label, title, axes in [
        (outer[3, 0], "G", "c-FOS/DAPI per animal", d_axes),
        (outer[3, 1], "H", "POMC activation per animal", e_axes),
    ]:
        pos = spec.get_position(fig)
        panel_label = fig.text(pos.x0-0.012, pos.y1+0.012, label, fontsize=20, fontweight="bold", ha="left", va="bottom")
        panel_title = fig.text((pos.x0+pos.x1)/2, pos.y1+0.014, title, fontsize=12.5, fontweight="bold", ha="center", va="bottom")
        panel_artists[label] = [*axes, panel_label, panel_title]

    pdf = outdir / f"{args.figure_name}.pdf"; png = outdir / f"{args.figure_name}.png"
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    panel_bboxes = {}
    for label, artists in panel_artists.items():
        boxes = [artist.get_tightbbox(renderer) for artist in artists]
        panel_bboxes[label] = Bbox.union(
            [box for box in boxes if box is not None]
        ).transformed(fig.dpi_scale_trans.inverted()).padded(0.04)

    fig.savefig(pdf, bbox_inches="tight", metadata={"Title": "Aligned ARC/ME c-FOS and POMC multipanel figure"})
    fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    for label, panel_path in panel_paths.items():
        fig.savefig(
            panel_path,
            bbox_inches=panel_bboxes[label],
            metadata={"Title": f"ARC/ME panel {label}"},
        )
        fig.savefig(
            panel_path.with_suffix(".png"),
            dpi=args.dpi,
            bbox_inches=panel_bboxes[label],
        )

    translation_receipt = translate_figure_texts_to_spanish(
        fig, extra=SPANISH_PANEL_TEXT
    )
    for artist, source in zip(
        spatial_image_artists, (spanish_spatial_i, spanish_spatial_j)
    ):
        artist.set_data(plt.imread(str(source)))
    fig.canvas.draw()
    for label, panel_path in spanish_panel_paths.items():
        fig.savefig(
            panel_path,
            bbox_inches=panel_bboxes[label],
            metadata={"Title": f"Panel ARC/ME {label} (español)"},
        )
        fig.savefig(
            panel_path.with_suffix(".png"),
            dpi=args.dpi,
            bbox_inches=panel_bboxes[label],
        )
    translation_receipt_path = (
        provenance_dir / "multipanel_spanish_translation_receipt.json"
    )
    translation_receipt_path.write_text(
        json.dumps(
            {
                "schema": "fig4_bilingual_isolated_panels_v1",
                "master_language": "en",
                "isolated_panel_locales": ["en", "es"],
                "panel_letters": list(PANEL_FILE_SUFFIXES),
                "raster_dpi": int(args.dpi),
                "text_replacements": translation_receipt,
                "spanish_spatial_sources": [
                    str(spanish_spatial_i), str(spanish_spatial_j),
                ],
            },
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)

    stats_df = pd.DataFrame(stats_multi)
    stats_df["analysis_family"] = "animal_anova"
    stats_path = source_data_dir / "multipanel_statistics.csv"
    sensitivity_path = source_data_dir / "multipanel_sensitivity_statistics.csv"
    cage_values_path = source_data_dir / "multipanel_sensitivity_cage_values.csv"
    denominator_audit_path = source_data_dir / "multipanel_denominator_audit.csv"
    methods_path = legends_dir / "multipanel_sensitivity_methods.txt"
    stats_df.to_csv(stats_path, index=False)
    sensitivity_df.to_csv(sensitivity_path, index=False)
    sensitivity_cages.to_csv(cage_values_path, index=False)
    denominator_audit.to_csv(denominator_audit_path, index=False)
    methods_path.write_text(
        "Figure 4 statistical methods\n"
        f"Source: {definitive_path}\n"
        f"Source SHA-256: {sha256_file(definitive_path)}\n"
        f"Registered section-count source: {section_region_path}\n"
        f"Registered section-count SHA-256: {sha256_file(section_region_path)}\n"
        "Animal endpoints were reaggregated by this figure maker from the current "
        "registered/HIL section-region counts; legacy plot-value and definitive-statistics "
        "CSVs were not used.\n"
        "Displayed animal analysis: ordinary equal-variance one-way ANOVA F; "
        "unadjusted exact exhaustive two-sided rank-sum/Mann-Whitney U pair tests. "
        "The displayed pairwise order is W–S, W–A, S–A; reported p values are raw.\n"
        "Cage-mean sensitivity: animal ratios are averaged arithmetically within each cage; "
        "the same ordinary one-way ANOVA and unadjusted exact two-sided MWU tests use cages as independent units.\n"
        "Count-aware sensitivity: numerator and denominator counts are summed within cage; "
        "the statistic is 2*sum_g[y_g*log(p_g/p_0)+(n_g-y_g)*log((1-p_g)/(1-p_0))], "
        "with zero terms omitted. Whole-cage labels are exhaustively permuted while preserving "
        "condition cage counts; reported p values are raw.\n"
        "Every exact MWU or deviance p value is extreme_labelings/enumerated_labelings; no Monte-Carlo plus-one "
        "correction is used. Positive denominators below 10 are flagged and retained. A zero or "
        "missing denominator is intrinsically non-estimable, is retained in the audit, and does "
        "not contribute a ratio or binomial likelihood. No finite positive-denominator value is filtered.\n",
        encoding="utf-8",
    )
    manifest = chosen.to_dict(); manifest.update({
        "selection_method": (
            "explicit original reviewed sample/section; registration, "
            "channel-completeness, and HIL gates still mandatory"
            if args.representative_sample else
            "verified overlap first, then confirmed manual overlap, then minimum "
            "summed robust biological deviation"
        ),
        "selection_candidate_count": int(len(candidates)), "selected_section_directory": str(folder),
        "microscopy_um_per_px": args.um_per_px, "microscopy_scalebar_um": args.scalebar_um,
        "inset_field_um": args.inset_field_um, "inset_scalebar_um": args.inset_scalebar_um,
        "raster_dpi": int(args.dpi),
        "master_language": "en",
        "isolated_subpanels_bilingual": True,
        "isolated_subpanel_locales": "en,es",
        "isolated_subpanel_pair_count": int(len(PANEL_FILE_SUFFIXES)),
        "spanish_translation_receipt": str(translation_receipt_path),
        "spanish_translation_receipt_sha256": sha256_file(translation_receipt_path),
        **{f"inset_{k}": v for k, v in inset_info.items()},
        "panel_E_content": "exact A+B+C merged microscopy with registered POMC-cell magnification inset",
        "panel_F_content": "registered DAPI/c-FOS/NPY-GFP microscopy plus real POMC intensity masked to accepted POMC ROIs; 40% microscopy-intensity ROI inset with dashed tissue/ventricle guides",
        **{
            f"panel_F_distribution_{key}": value
            for key, value in panel_f_distribution_info.items()
        },
        "panel_F_source": "hash-verified registered DAPI/c-FOS/POMC/NPY-GFP TIFFs plus processed POMC label TIFF; POMC intensity displayed only within accepted POMC ROIs",
        "panel_F_raw_pomc_intensity_displayed": True,
        "panel_F_raw_pomc_intensity_display_scope": "processed POMC label>0 pixels only",
        "panel_F_pomc_intensity_source": str(sources["pomc"]),
        "panel_F_pomc_intensity_sha256": sha256_file(sources["pomc"]),
        "panel_F_pomc_roi_label_source": str(pomc_size_path),
        "panel_F_pomc_size_qc": pomc_size_qc,
        "panel_F_pomc_roi_label_sha256": sha256_file(pomc_size_path),
        "panel_F_pomc_signal_gain": 1.0,
        "panel_F_pomc_signal_normalization": "same full-field positive-pixel 0.5/99.7 percentiles, gamma 0.82 and amber tint as panels C/E; no additional gain or ROI-specific normalization",
        "panel_F_pomc_signal_mask_method": "real normalized POMC intensity multiplied by label>0; zero contribution outside; no flat paint/outline/dilation/filter/resampling",
        "panel_F_pomc_signal_pixels_outside_roi": 0,
        "panel_F_pomc_blend": "screen over DAPI/c-FOS; no additive clipping",
        "panel_F_inset_uses_microscopy_intensities": True,
        "panel_F_inset_uniform_roi_fill": False,
        "panel_F_npy_gfp_intensity_source": str(sources["npy"]),
        "panel_F_npy_gfp_intensity_sha256": sha256_file(sources["npy"]),
        "panel_F_npy_gfp_display_normalization": "whole-field percentiles 50.0/99.7; clip 0/1; gamma 0.72; no mask/filter/resampling",
        "panel_F_npy_gfp_blend": "screen over exact DAPI+c-FOS plus real ROI-masked POMC intensity",
        "panel_F_npy_gfp_tint_rgb": ",".join(f"{value:g}" for value in NPY_GFP_COMPOSITE_TINT),
        "panel_F_npy_gfp_microscopy_roi_gated": False,
        "panel_F_npy_gfp_roi_labels_used_in_main_microscopy": False,
        "panel_F_npy_gfp_roi_labels_used_in_microscopy_inset": True,
        "panel_F_screen_blend_cannot_reduce_dapi_cfos_base_rgb_components": True,
        "panel_F_tissue_outline": "white dashed DAPI-supported exterior in ROI inset",
        "panel_F_ventricle_outline": "white dashed HIL-ARC-derived lumen in ROI inset",
        "panel_F_region_source": "accepted HIL ARC lobes used only to derive the ventricular guide; no ARC/ME fill or contour",
        "panel_F_region_mask_sha256": review_provenance["mask_sha256"],
        "panel_G_source": str(values_path), "panel_G_endpoint": d_col,
        "panel_H_source": str(values_path), "panel_H_endpoint": e_col,
        "sensitivity_source": str(definitive_path),
        "sensitivity_source_sha256": sha256_file(definitive_path),
        "statistics_recomputed_from_registered_sections": True,
        "registered_section_count_source": str(section_region_path),
        "registered_section_count_source_sha256": sha256_file(section_region_path),
        "all_section_registration_audit": str(registration_audit_path),
        "all_section_registration_audit_sha256": sha256_file(registration_audit_path),
        "sensitivity_statistics": str(sensitivity_path),
        "sensitivity_cage_values": str(cage_values_path),
        "denominator_audit": str(denominator_audit_path),
        "per_cell_source": str(cell_path), "summary_source": str(summary_path),
        "arc_me_review_policy": (
            "strict_complete_native_hil" if args.require_manual_regions
            else "optional_hil_all_automated_unless_complete_71_of_71"
        ),
        "human_review_gate_status": (
            "passed_strict" if args.require_manual_regions
            else ("complete_optional" if all_manual else "not_required_automated")
        ),
        "quantifiable_sections": int(len(quantified)),
        "stitch_registration_schema_version": registration["schema_version"],
        "stitch_registration_algorithm": registration["algorithm"],
        "stitch_registration_status": registration["status"],
        "stitch_registration_method": registration["method"],
        "stitch_registration_json": registration["_registration_path"],
        "stitch_registration_json_sha256": registration["_registration_sha256"],
        "stitch_registration_right_origin_x_px": registration["_right_origin_x_px"],
        "stitch_registration_right_shift_y_px": registration["_right_shift_y_px"],
        "stitch_registration_right_medial_crop_px": registration["_right_medial_crop_px"],
        "stitch_registration_overlap_width_px": registration["_overlap_width_px"],
        "stitch_registration_gap_px": registration["_gap_px"],
        "stitch_registration_raw_overlap_ncc": registration["_raw_overlap_ncc"],
        "stitch_registration_gradient_overlap_ncc": registration["_gradient_overlap_ncc"],
        "stitch_registration_ncc_best_score": registration["_ncc_best_score"],
        "stitch_registration_ncc_second_score": registration["_ncc_second_score"],
        "stitch_registration_ncc_uniqueness_margin": registration["_ncc_uniqueness_margin"],
        "stitch_registration_hil_receipt_bound": registration["_registration_receipt_bound"],
        "stitch_registration_manual_confirmation": registration["_manual_registration_confirmed"],
        **{
            f"stitch_registration_output_{key}_sha256": value
            for key, value in registration["_artifact_sha256s"].items()
        },
        **{
            f"stitch_registration_input_{side}_dapi_sha256": value
            for side, value in registration["_input_dapi_sha256s"].items()
        },
        **{
            f"panel_{letter}_single_channel": marker
            for letter, marker in zip("ABC", COMPOSITE_CHANNEL_ORDER)
        },
        **{
            f"panel_{letter}_composite_tint_rgb": ",".join(
                f"{value:g}" for value in COMPOSITE_CHANNEL_TINTS[marker]
            )
            for letter, marker in zip("ABC", COMPOSITE_CHANNEL_ORDER)
        },
        "panels_ABC_clipped_sum_equals_panel_E": True,
        "panels_ABC_registered_field_shape_px": f"{dapi.shape[0]}x{dapi.shape[1]}",
        "panels_ABC_registered_bounds_px": f"x[-0.5,{dapi.shape[1]-0.5}],y[{dapi.shape[0]-0.5},-0.5]",
        "panel_D_registered_field_shape_px": f"{dapi.shape[0]}x{dapi.shape[1]}",
        "panel_E_registered_field_shape_px": f"{regions.shape[0]}x{regions.shape[1]}",
        "panel_D_registered_dapi_sha256": registration["_dapi_sha256"],
        "panel_E_registered_dapi_sha256": review_provenance["source_dapi_sha256"],
        "panel_D_registered_bounds_px": f"x[-0.5,{dapi.shape[1]-0.5}],y[{dapi.shape[0]-0.5},-0.5]",
        **{
            f"panel_D_hil_arc_ventricle_{key}": value
            for key, value in hil_ventricle_info.items()
        },
        "panel_E_registered_bounds_px": f"x[-0.5,{dapi.shape[1]-0.5}],y[{dapi.shape[0]-0.5},-0.5]",
        "panel_E_same_registered_coordinates_hash_and_bounds_as_D": True,
        "panel_F_registered_field_shape_px": f"{tissue_envelope.shape[0]}x{tissue_envelope.shape[1]}",
        "panel_F_registered_bounds_px": f"x[-0.5,{dapi.shape[1]-0.5}],y[{dapi.shape[0]-0.5},-0.5]",
        "panel_F_same_registered_coordinates_as_E": True,
        "panel_F_dapi_cfos_contributions_match_E": True,
        "panel_F_pomc_contribution_matches_E": False,
        "panel_F_pomc_is_registered_roi_paint_only": False,
        "panel_F_pomc_is_registered_intensity_masked_to_rois": True,
        "panel_F_raw_npy_gfp_reconstruction_matches_analysis_labels": True,
        "panel_F_microscopy_uses_segmentation_masks": True,
        "panel_F_microscopy_segmentation_mask_scope": "POMC intensity gate only; DAPI/c-FOS/NPY-GFP are ungated microscopy intensities",
        "panel_F_microscopy_inset_uses_segmentation_masks": True,
        "panel_E_displayed_regions": "",
        "panel_F_displayed_regions": "",
        "panel_F_distribution_inset_region_contour_width_pt": 0.0,
        "panel_F_distribution_inset_roi_contour_width_pt": 0.0,
        "panel_F_distribution_inset_marker_intensities_are_masked_without_outlines": True,
        "panel_F_distribution_inset_has_white_dashed_tissue_boundary": True,
        "panel_F_distribution_inset_has_white_dashed_ventricle_boundary": True,
        "panel_E_displayed_mask_is_exact_accepted_HIL": False,
        "panel_F_displayed_mask_is_exact_accepted_HIL": False,
        "panel_E_wall_geometry_drawn": False,
        "panel_F_wall_geometry_drawn": False,
        **{f"human_review_{key}": value for key, value in review_provenance.items()},
    })
    for marker_name, marker_record in cellpose_provenance.items():
        safe_marker = marker_name.lower().replace("-", "")
        for key, value in marker_record.items():
            manifest[f"panel_F_{safe_marker}_{key}"] = value
    for name, path in sources.items():
        exists = path.is_file()
        manifest[f"source_{name}"] = str(path) if exists else ""
        manifest[f"source_{name}_exists"] = bool(exists)
    manifest_path = provenance_dir / "multipanel_source_manifest.csv"
    pd.DataFrame([manifest]).to_csv(manifest_path, index=False)
    size_records = []
    for qc_path in sorted((root / "reconstructed_sections").glob("*/*/pomc_size_qc.json")):
        qc = json.loads(qc_path.read_text())
        for obj in qc["objects"]:
            size_records.append({"sample": qc["sample"], "section": qc["section_index"], **obj})
    pd.DataFrame(size_records).to_csv(values_path.parent / "pomc_size_qc_all_objects.csv", index=False)
    write_legends(
        legends_dir, args, chosen, inset_info, stats_df, sensitivity_df,
        values_path, review_provenance, registration, cellpose_provenance,
        panel_f_distribution_info,
    )
    print(
        f"[SELECTED] {chosen['sample']} S{int(chosen['section_index']):02d} "
        f"({chosen['cond']}); method={registration['method']}; "
        f"right_origin_x={registration['_right_origin_x_px']} px; "
        f"shift_y={registration['_right_shift_y_px']} px; "
        f"overlap={registration['_overlap_width_px']} px; "
        f"status={registration['status']}"
    )
    for path in [*panel_paths.values(), *spanish_panel_paths.values(), pdf, png]:
        print(f"[FIGURE] {path}")
    print(f"[LEGEND] {legends_dir / (args.figure_name + '_LEGEND.txt')}")
    return (
        *panel_paths.values(), *spanish_panel_paths.values(), pdf, png,
        manifest_path, translation_receipt_path, stats_path,
        sensitivity_path, cage_values_path, denominator_audit_path, methods_path,
    )


def main() -> int:
    args = build_parser().parse_args()
    if int(args.dpi) != 600:
        raise SystemExit(f"Figure 4 publication graphics are fixed at 600 dpi, got {args.dpi}")
    if args.um_per_px <= 0 or args.inset_field_um <= 0: raise SystemExit("Pixel size and inset field must be > 0")
    if args.inset_scalebar_um >= args.inset_field_um * 0.48: raise SystemExit("Inset scale bar must be < 48% of inset field width")
    make_figure(args); return 0


if __name__ == "__main__":
    raise SystemExit(main())
