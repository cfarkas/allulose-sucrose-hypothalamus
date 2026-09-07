#!/usr/bin/env python3
"""Build one publication plate for ARC/ME microglia and GFAP analyses."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import csv
import datetime as dt
import json
import math
import os
import re
import sys
from itertools import combinations
from pathlib import Path
from typing import Any, Sequence

# Stabilize PDF CreationDate so command-level reruns are byte-reproducible.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1761264000")

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 1.8,
    "axes.titlesize": 14.0,
    "axes.titleweight": "bold",
    "axes.labelsize": 12.0,
    "axes.labelweight": "bold",
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 11.0,
})
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import tifffile
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, FancyArrowPatch, Patch, Rectangle
from matplotlib.ticker import ScalarFormatter
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import (
    convolve, distance_transform_edt, gaussian_filter,
    label as ndi_label, uniform_filter,
)
from scipy.stats import f_oneway, mannwhitneyu, rankdata
from skimage.morphology import remove_small_objects

HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")
PROJECT_ROOT = PAPER_ROOT
ROOT = PAPER_ROOT
FIG5_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig5"
CONDITIONS = ("Water", "Sucrose", "Allulose")
# Quantitative panels follow the canonical annotation-code order and include
# every analyzed hypothalamic ROI. VMN must never be silently dropped.
REGIONS = ("ARC", "ME", "VMN")
PANEL_A_INSET_SIZE = 112
COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
REGION_COLORS = {"ARC": "#00c7d9", "ME": "#ff9d2e", "VMN": "#be3eae"}
REGION_CODES = {"ARC": 1, "ME": 2, "VMN": 3}
SHARED_SCRIPTS = PAPER_ROOT / "scripts" / "shared"
if str(SHARED_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SHARED_SCRIPTS))
from spanish_panel_text import translate_figure_texts_to_spanish

SPANISH_PANEL_TEXT = {
    "HIL-accepted ARC/ME/VMN": "ARC/ME/VMN aceptados mediante HIL",
    "Accepted:": "Aceptado:",
    "manual polygons": "polígonos manuales",
    "native auto mask, HIL-confirmed": "máscara automática nativa, confirmada mediante HIL",
    "Ramified cell": "Célula ramificada",
    "GFAP processes": "Procesos GFAP+",
    "selected cell": "célula elegida",
    "Image-patch CNN": "CNN de parches de imagen",
    "Iba1 | mask | skeleton": "Iba1 | máscara | esqueleto",
    "four-class softmax": "softmax de cuatro clases",
    "pseudo-labels": "pseudoetiquetas",
    "animal-blocked": "partición por animal",
    "best epoch": "mejor época",
    "pseudo-label agreement": "concordancia de pseudoetiquetas",
    "Model-predicted cell examples": "Ejemplos de células predichas por el modelo",
    "No QC-valid": "Sin predicción válida por CC",
    " prediction": " predicha",
    "Representative-section microglial states": "Estados microgliales en sección representativa",
    "CNN-predicted cells (%)": "Células predichas por CNN (%)",
    "GFAP quantification": "Cuantificación de GFAP",
    "background correction": "corrección de fondo",
    "+ threshold": "+ umbral",
    "Raw GFAP + DAPI": "GFAP + DAPI crudos",
    "GFAP-positive mask": "Máscara GFAP positiva",
    "Corrected signal / ROI": "Señal corregida / ROI",
    "raw intensity": "intensidad cruda",
    "GFAP+ pixels / ROI": "píxeles GFAP+ / ROI",
    "positive fraction": "fracción positiva",
    "Channel analysis masks": "Máscaras de análisis por canal",
    "Tissue boundary": "Límite tisular",
    "Iba1 corrected signal per animal": "Señal corregida de Iba1 por animal",
    "Iba1 signal": "Señal Iba1",
    "Activated fraction": "Fracción activada",
    "GFAP corrected signal per animal": "Señal corregida de GFAP por animal",
    "GFAP signal": "Señal GFAP",
    "(raw a.u./µm²)": "(u.a. crudas/µm²)",
    "MWU W–S": "MWU Ag–Sac",
    "MWU W–A": "MWU Ag–Alu",
    "MWU S–A": "MWU Sac–Alu",
    "Accepted cells (%)": "Células aceptadas (%)",
    "Microglial state composition": "Composición de estados microgliales",
    "HIL anatomy": "Anatomía HIL",
    "Cellpose and HIL anatomy": "Cellpose y anatomía HIL",
    "Within-region microglial states": "Estados microgliales dentro de la región",
    "Microglia morphology classifier": "Clasificador de morfología microglial",
    "GFAP quantification workflow": "Flujo de cuantificación de GFAP",
    "Raw GFAP / DAPI": "GFAP / DAPI crudos",
    "Final GFAP-positive mask": "Máscara final GFAP positiva",
    "Channel-separated analysis masks": "Máscaras de análisis separadas por canal",
    "Raw-scale Iba1 signal / area": "Señal Iba1 en escala cruda / área",
    "Activated microglia fraction": "Fracción de microglía activada",
    "Raw-scale GFAP signal / area": "Señal GFAP en escala cruda / área",
    "Iba1 signal per animal": "Señal Iba1 por animal",
    "Activated microglia per animal": "Microglía activada por animal",
    "GFAP signal per animal": "Señal GFAP por animal",
    "Representative section": "Sección representativa",
    "Predicted state": "Estado predicho",
    "Ramified": "Ramificada",
    "Rod-like": "En bastón",
    "Activated": "Activada",
    "Amoeboid": "Ameboide",
    "Animal mean ± SD": "Media animal ± DE",
    "animal mean ± SD": "media animal ± DE",
    "individual animals": "animales individuales",
    "Raw signal / analyzed area": "Señal cruda / área analizada",
    "Fraction activated": "Fracción activada",
    "Cell count": "Recuento celular",
    "Analysis mask": "Máscara de análisis",
    "Source microscopy": "Microscopía de origen",
}
FIG5_FIGURE_ROOT = PAPER_ROOT / "Fig5"


def _latest_reviewed_root(parent, pattern):
    """Newest directory matching ``pattern`` under ``parent`` that holds receipts."""
    if not parent.is_dir():
        return None
    roots = sorted(
        path for path in parent.glob(pattern)
        if path.is_dir() and any(path.glob("*/*/*.json"))
    )
    return roots[-1] if roots else None


def _first_existing_dir(*candidates):
    """First candidate that exists, else the first named one.

    Directories under analyses/ are generated working output and are absent in a
    freshly unpacked copy, so each default falls back to what this tree ships.
    """
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    return next((candidate for candidate in candidates if candidate is not None), None)


# The analyzer's own default output comes first, so the two documented commands
# chain with no arguments; final_hybrid_run is the older working copy and is used
# only when no human-final run has been produced.
DEFAULT_ANALYSIS_ROOT = _first_existing_dir(
    FIG5_ANALYSIS_ROOT / "results" / "human_final_run",
    FIG5_ANALYSIS_ROOT / "results" / "final_hybrid_run",
)
DEFAULT_INPUT_ROOT = _first_existing_dir(
    FIG5_ANALYSIS_ROOT / "raw" / "input", FIG5_FIGURE_ROOT / "raw_data",
)
DEFAULT_MANUAL_REGION_ROOT = _first_existing_dir(
    _latest_reviewed_root(FIG5_ANALYSIS_ROOT / "hil", "accepted_annotations"),
    _latest_reviewed_root(FIG5_FIGURE_ROOT / "hil_review", "human_final_*"),
    FIG5_ANALYSIS_ROOT / "hil" / "accepted_annotations",
)


def fully_delineated_sections(manual_root) -> list:
    """Accepted receipts under ``manual_root`` that carry all three drawn regions.

    A representative panel has to show ARC, ME and VMN, so a review that left one
    of them not delineated on the chosen section blocks the render with nothing
    the reviewer can act on. Naming the sections that do qualify turns the
    refusal into a choice: draw the missing region, or pick one of these.
    """
    found = []
    for path in sorted(Path(manual_root).rglob("*.json")):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not record.get("accepted"):
            continue
        status = record.get("region_status", {})
        if all(status.get(region) == "drawn" for region in REGIONS):
            found.append((str(record.get("animal", "")), str(record.get("section", ""))))
    return found
EXPECTED_HUMAN_REVIEW_SECTIONS = 61
ALLOWED_NATIVE_ACCEPTED_SOURCES = {"manual_replacement_polygons", "automated_mask_unchanged"}
ALL_HIL_REGIONS = ("ARC", "ME", "VMN")
MICROGLIA_DENSITY_SCALE_UM2 = 100000.0
STATE_COLORS = {
    "Ramified": "#17627a",
    "Rod-like": "#08a7a4",
    "Activated": "#f1b83b",
    "Amoeboid": "#f06455",
}
CHANNEL_CARTOON_COLORS = {
    "DAPI": "#4f86ff",
    "Iba1": "#ff334d",
    "GFAP": "#008f45",
}
CARTOON_TISSUE_FILL = "#e8edf2"
CARTOON_TISSUE_EDGE = "#354052"
PANEL_ENDPOINTS = (
    ("G", "iba1_corrected_intensity_per_um2", "Iba1 corrected signal per animal", "Iba1 signal\n(raw a.u./µm²)"),
    ("H", "activated_microglia_fraction", "Activated microglia per animal", "Activated fraction"),
    ("I", "gfap_corrected_intensity_per_um2", "GFAP corrected signal per animal", "GFAP signal\n(raw a.u./µm²)"),
)


def _path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _portable_provenance_path(
    path: Path | str, analysis_dir: Path, output_dir: Path,
) -> str:
    """Describe a source without binding publication metadata to a build root."""
    resolved = Path(path).expanduser().resolve()
    for root, prefix in (
        (Path(analysis_dir).expanduser().resolve(), "analysis"),
        (Path(output_dir).expanduser().resolve(), "publication"),
        (PAPER_ROOT.resolve(), ""),
    ):
        try:
            relative = resolved.relative_to(root)
        except ValueError:
            continue
        return (Path(prefix) / relative).as_posix() if prefix else relative.as_posix()
    return resolved.as_posix()


def _lexical_absolute_path(path: Path, *, base: Path | None = None) -> Path:
    """Normalize historical receipt paths without filesystem probes."""
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = Path(base) / candidate if base is not None else Path.cwd() / candidate
    return Path(os.path.abspath(os.fspath(candidate)))


def _review_provenance_path_matches(
    value: Any,
    expected_path: Path,
    *,
    repository_root: Path = PAPER_ROOT,
) -> bool:
    """Match absolute, section-relative, or Paper-relative HIL source paths."""
    text = str(value or "").strip()
    if not text:
        return False
    expected = _lexical_absolute_path(Path(expected_path))
    candidate = Path(text).expanduser()
    root = _lexical_absolute_path(Path(repository_root))
    if candidate.is_absolute():
        candidates = [candidate]
        # Accepted receipts predate the relocatable public bundle and may
        # retain an absolute path beneath the original Paper/Fig5 tree. Only
        # rebase only its unique Fig5 suffix under the current repository; the equality
        # check below still requires the exact expected source.
        fig5_indices = [
            index for index, part in enumerate(candidate.parts) if part == "Fig5"
        ]
        if len(fig5_indices) == 1:
            rebased = root.joinpath(*candidate.parts[fig5_indices[0] :])
            if _path_is_within(rebased, root):
                candidates.append(rebased)
    else:
        candidates = [expected.parent / candidate, root / candidate]
    for path in candidates:
        try:
            if _lexical_absolute_path(path, base=root) == expected:
                return True
        except (OSError, RuntimeError):
            continue
    return False


def _receipt_base_animal(value: str) -> str:
    """Mirror the analyzer's conservative repetition-suffix normalization."""
    text = str(value).strip()
    match = re.search(r"\s*\(([^)]*)\)\s*$", text)
    if match:
        return re.sub(r"\s*\([^)]*\)\s*$", "", text).strip()
    match = re.match(
        r"^(.*?)(?:[_-]?(?:rep|repeat|repetition))([0-9]+)$",
        text,
        flags=re.IGNORECASE,
    )
    if match:
        return match.group(1).strip("_- ")
    match = re.match(
        r"^(FR\s*-?\d+(?:-\d+)?|FR\d+)(R)$",
        text,
        flags=re.IGNORECASE,
    )
    return match.group(1) if match else text


def _native_receipt_candidates(
    manual_root: Path,
    animal: str,
    section: str,
) -> list[tuple[str, Path, Path]]:
    animals = [str(animal)]
    base = _receipt_base_animal(str(animal))
    if base != str(animal):
        animals.append(base)
    candidates: list[tuple[str, Path, Path]] = []
    for receipt_animal in animals:
        pair_dir = Path(manual_root) / receipt_animal / section
        stem = f"{receipt_animal}_{section}_ARC_ME_labels"
        candidates.append(
            (
                receipt_animal,
                pair_dir / f"{stem}.tif",
                pair_dir / f"{stem}.json",
            )
        )
    return candidates


def _absolute_lexical(path: Path) -> Path:
    path = path.expanduser()
    return path if path.is_absolute() else Path.cwd() / path


def _force_clean_generated_dir(
    output_dir: Path,
    *,
    forbidden_roots: Sequence[Path],
    preserved_roots: Sequence[Path] = (),
    preserved_names: Sequence[str] = (),
) -> Path:
    """Empty one generated directory without following links or deleting annotations."""
    output_lex = _absolute_lexical(Path(output_dir))
    if output_lex.is_symlink():
        raise ValueError(f"Refusing --force on symlink output directory: {output_lex}")
    output = output_lex.resolve(strict=False)
    filesystem_root = Path(output.anchor).resolve()
    if output == filesystem_root:
        raise ValueError(f"Refusing --force on filesystem root: {output}")

    forbidden_pairs = [
        (_absolute_lexical(Path(path)), _absolute_lexical(Path(path)).resolve(strict=False))
        for path in forbidden_roots
        if path is not None
    ]
    for protected_lex, protected in forbidden_pairs:
        if (
            output_lex == protected_lex
            or output == protected
            or _path_is_within(protected_lex, output_lex)
            or _path_is_within(protected, output)
        ):
            raise ValueError(
                f"Refusing --force because output is or contains protected path "
                f"{protected}: {output}"
            )

    preserved_pairs = [
        (_absolute_lexical(Path(path)), _absolute_lexical(Path(path)).resolve(strict=False))
        for path in preserved_roots
        if path is not None
    ]
    for preserve_lex, preserve in preserved_pairs:
        if output_lex == preserve_lex or output == preserve:
            raise ValueError(f"Refusing --force on manual annotation root: {output}")
        if _path_is_within(output_lex, preserve_lex) or _path_is_within(output, preserve):
            raise ValueError(
                f"Refusing --force because output lies inside manual annotation root "
                f"{preserve}: {output}"
            )

    if output_lex.exists() and not output_lex.is_dir():
        raise ValueError(f"Refusing --force because output is not a directory: {output_lex}")
    if not output_lex.exists():
        return output

    def clean(directory_lex: Path) -> None:
        for child_lex in directory_lex.iterdir():
            if directory_lex == output_lex and child_lex.name in set(preserved_names):
                continue
            child_abs = _absolute_lexical(child_lex)
            child = child_abs.resolve(strict=False)
            exact_preserve = any(
                child_abs == preserve_lex or child == preserve
                for preserve_lex, preserve in preserved_pairs
            )
            contains_preserve = any(
                _path_is_within(preserve_lex, child_abs)
                or _path_is_within(preserve, child)
                for preserve_lex, preserve in preserved_pairs
            )
            if exact_preserve:
                continue
            if contains_preserve:
                if child_abs.is_symlink() or not child_abs.is_dir():
                    raise ValueError(
                        f"Cannot preserve nested manual annotation path through "
                        f"non-directory or symlink: {child_abs}"
                    )
                clean(child_abs)
                continue
            if child_abs.is_symlink() or not child_abs.is_dir():
                child_abs.unlink()
            else:
                shutil.rmtree(child_abs)

    clean(output_lex)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--manual-region-dir", type=Path, default=DEFAULT_MANUAL_REGION_ROOT, help="Native HIL-accepted ARC/ME/VMN TIFF+JSON root.")
    parser.add_argument(
        "--output-dir", type=Path,
        default=FIG5_ANALYSIS_ROOT / "figure" / "human_final_render",
        help=("Generated render directory. It must never default to the live Fig5 "
              "directory: promotion into it is a separate, deliberate copy."),
    )
    parser.add_argument("--iba1-animal", default="FR722", help="Representative animal for panel A (Iba1/DAPI).")
    parser.add_argument("--iba1-section", default="S03", help="All-three-region representative section for panel A (Iba1/DAPI).")
    parser.add_argument("--iba1-cell-uid", default="fr722__s03__iba1_00113", help="Audited large branched QC-valid cell to center in the panel A source box.")
    # Panel B needs a representative carrying drawn ARC, ME and VMN. FR8-3 was the
    # original choice, but its HIL review delineated ARC and VMN only, on all
    # four of its sections, so it contributes no ME at all and cannot be drawn.
    # FR8-1 is the same Allulose/WT/F group and all four of its sections carry
    # the three regions. Superseded default: --animal FR8-3 --section S04.
    parser.add_argument("--animal", default="FR8-1")
    parser.add_argument("--section", default="S04")
    parser.add_argument("--um-per-px", type=float, default=1.51)
    parser.add_argument("--seed", type=int, default=20260812, help="Deterministic animal-point jitter seed for panels G-I.")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--formats", default="pdf,png", help="Combined-figure formats; publication layout permits PDF and PNG only.")
    parser.add_argument("--force", action="store_true", help="Safely clean only the designated publication output before rendering; protected inputs and HIL annotations are never removed.")
    return parser.parse_args()


def marker_signal(array: np.ndarray) -> np.ndarray:
    array = np.asarray(array)
    if array.ndim == 2:
        return array.astype(float)
    if array.ndim == 3:
        return np.max(array[..., :3].astype(float), axis=2)
    return np.squeeze(array).astype(float)


def normalize(array: np.ndarray, low: float = 1.0, high: float = 99.7, gamma: float = 0.82) -> np.ndarray:
    array = np.asarray(array, dtype=float)
    finite = array[np.isfinite(array)]
    if not finite.size:
        return np.zeros_like(array, dtype=float)
    lo, hi = np.percentile(finite, (low, high))
    if hi <= lo:
        return np.zeros_like(array, dtype=float)
    scaled = np.clip((array - lo) / (hi - lo), 0.0, 1.0)
    return np.power(scaled, gamma)


def composite(dapi: np.ndarray, marker: np.ndarray, color: tuple[float, float, float]) -> np.ndarray:
    d = normalize(marker_signal(dapi), high=99.6, gamma=0.78)
    m = normalize(marker_signal(marker), high=99.7, gamma=0.72)
    rgb = np.zeros((*d.shape, 3), dtype=float)
    rgb[..., 2] += 0.82 * d
    rgb[..., 0] += 0.08 * d
    for channel, coefficient in enumerate(color):
        rgb[..., channel] += coefficient * m
    return np.clip(rgb, 0.0, 1.0)


def crop_square(array: np.ndarray, center: tuple[float, float], size: int) -> tuple[np.ndarray, tuple[int, int, int, int]]:
    height, width = array.shape[:2]
    half = size // 2
    cx, cy = int(round(center[0])), int(round(center[1]))
    x0, x1 = max(0, cx - half), min(width, cx + half)
    y0, y1 = max(0, cy - half), min(height, cy + half)
    result = np.zeros((size, size, *array.shape[2:]), dtype=array.dtype) if array.ndim == 3 else np.zeros((size, size), dtype=array.dtype)
    dx, dy = half - (cx - x0), half - (cy - y0)
    result[dy:dy + (y1 - y0), dx:dx + (x1 - x0)] = array[y0:y1, x0:x1]
    return result, (x0, y0, x1, y1)


def add_scalebar(ax: plt.Axes, width_px: int, height_px: int, bar_um: float, um_per_px: float, inset: bool = False) -> None:
    bar_px = bar_um / um_per_px
    # Small 96/118 px insets need enough left padding for the full "25 µm"
    # label and its dark backing; 5.5% clipped the leading digit in Panel C.
    x0 = width_px * (0.22 if inset else 0.055)
    y = height_px * 0.925
    line_width = 3.0 if inset else 4.5
    font_size = 6.5 if inset else 9.5
    ax.add_patch(Rectangle((x0 - 8, y - 18), bar_px + 16, 29, facecolor="black", edgecolor="none", alpha=0.72, zorder=20))
    ax.plot([x0, x0 + bar_px], [y, y], color="white", linewidth=line_width, solid_capstyle="butt", zorder=21)
    ax.text(x0 + bar_px / 2, y - 7, f"{bar_um:g} µm", color="white", ha="center", va="bottom", fontsize=font_size, fontweight="bold", zorder=21)


def channel_cartoon_tissue_envelope(
    dapi_mask: np.ndarray, regions: np.ndarray,
) -> np.ndarray:
    """Create one shared DAPI-supported tissue envelope for all Panel F channels."""
    dapi_positive = np.asarray(dapi_mask, dtype=bool)
    region_codes = np.asarray(regions)
    if dapi_positive.ndim != 2 or dapi_positive.shape != region_codes.shape:
        raise RuntimeError(
            "Panel F DAPI and accepted ARC/ME/VMN masks must be one registered 2D field"
        )
    sigma = max(6.0, 0.025 * float(dapi_positive.shape[0]))
    density = gaussian_filter(
        dapi_positive.astype(np.float32), sigma=sigma, mode="constant",
    )
    tissue = remove_small_objects(
        density > 0.02,
        min_size=max(256, int(round(dapi_positive.size * 0.002))),
    )
    tissue |= region_codes > 0
    if not np.any(tissue):
        raise RuntimeError("Panel F could not derive a DAPI-supported tissue envelope")
    return np.asarray(tissue, dtype=bool)


def channel_cartoon_region_label_points(
    regions: np.ndarray,
) -> list[tuple[str, float, float]]:
    """Return stable interior labels for both ARC and VMN lobes and the single ME."""
    region_codes = np.asarray(regions)
    if region_codes.ndim != 2:
        raise RuntimeError("Panel F accepted-region mask must be 2D")
    expected_components = {"ARC": 2, "ME": 1, "VMN": 2}
    points: list[tuple[str, float, float]] = []
    for region_name in REGIONS:
        components, count = ndi_label(region_codes == REGION_CODES[region_name])
        expected = expected_components[region_name]
        if int(count) != expected:
            raise RuntimeError(
                f"Panel F {region_name} geometry has {int(count)} connected components; "
                f"expected {expected} for the audited {REGIONS} representative"
            )
        ranked = sorted(
            range(1, int(count) + 1),
            key=lambda component: int(np.count_nonzero(components == component)),
            reverse=True,
        )
        for component in ranked:
            inside = components == component
            distance = distance_transform_edt(inside)
            y, x = np.unravel_index(int(np.argmax(distance)), distance.shape)
            points.append((region_name, float(x), float(y)))
    return points


def add_channel_cartoon_scalebar(
    ax: plt.Axes, shape: tuple[int, int], um_per_px: float,
    length_um: float = 100.0,
) -> None:
    """Use the high-contrast, bottom-right scale-bar design from Figure 4C."""
    height, width = shape
    length_px = float(length_um) / float(um_per_px)
    if length_px <= 0 or length_px >= 0.48 * width:
        raise RuntimeError(
            f"Invalid Panel F scale bar: {length_um:g} um -> {length_px:.1f}px"
        )
    margin_x = max(30.0, 0.026 * width)
    margin_y = max(25.0, 0.055 * height)
    bar_h = max(11.0, 0.017 * height)
    x0 = width - margin_x - length_px
    y0 = height - margin_y
    pad_x = max(5.0, 0.025 * length_px)
    text_h = max(16.0, 0.060 * height)
    ax.add_patch(Rectangle(
        (x0 - pad_x, y0 - bar_h - text_h),
        length_px + 2 * pad_x, bar_h + text_h + 7,
        facecolor="black", edgecolor="black", linewidth=1.2,
        alpha=0.92, zorder=30,
    ))
    ax.add_patch(Rectangle(
        (x0, y0 - bar_h), length_px, bar_h,
        facecolor="white", edgecolor="white", linewidth=1.0, zorder=31,
    ))
    label = ax.text(
        x0 + length_px / 2, y0 - bar_h - 4, f"{length_um:g} µm",
        color="white", ha="center", va="bottom", fontsize=6.2,
        fontweight="bold", zorder=32,
    )
    label.set_path_effects([pe.withStroke(linewidth=3.2, foreground="black")])


def draw_channel_analysis_cartoon(
    ax: plt.Axes, marker: str, marker_mask: np.ndarray,
    tissue: np.ndarray, regions: np.ndarray, um_per_px: float,
) -> None:
    """Draw one exact-coordinate channel analysis-mask cartoon for Panel F."""
    if marker not in CHANNEL_CARTOON_COLORS:
        raise RuntimeError(f"Unknown Panel F channel: {marker}")
    mask = np.asarray(marker_mask, dtype=bool)
    tissue_mask = np.asarray(tissue, dtype=bool)
    region_codes = np.asarray(regions)
    if mask.ndim != 2 or not (mask.shape == tissue_mask.shape == region_codes.shape):
        raise RuntimeError("Panel F cartoon layers do not share one registered 2D field")

    rgb = np.ones((*mask.shape, 3), dtype=np.float32)
    rgb[tissue_mask] = np.asarray(to_rgb(CARTOON_TISSUE_FILL), dtype=np.float32)
    for region_name in REGIONS:
        inside = region_codes == REGION_CODES[region_name]
        if not np.any(inside):
            raise RuntimeError(f"Panel F accepted mask is missing {region_name}")
        color = np.asarray(to_rgb(REGION_COLORS[region_name]), dtype=np.float32)
        rgb[inside] = 0.88 * rgb[inside] + 0.12 * color
    ax.imshow(np.clip(rgb, 0.0, 1.0), interpolation="nearest")

    marker_layer = np.zeros((*mask.shape, 4), dtype=np.float32)
    marker_layer[..., :3] = to_rgb(CHANNEL_CARTOON_COLORS[marker])
    marker_layer[..., 3] = mask.astype(np.float32) * (0.78 if marker == "DAPI" else 0.90)
    ax.imshow(marker_layer, interpolation="nearest")
    ax.contour(tissue_mask, [0.5], colors=["white"], linewidths=3.5, alpha=0.95)
    ax.contour(tissue_mask, [0.5], colors=[CARTOON_TISSUE_EDGE], linewidths=1.7)
    for region_name in REGIONS:
        inside = region_codes == REGION_CODES[region_name]
        ax.contour(inside, [0.5], colors=["#172033"], linewidths=3.2, alpha=0.92)
        ax.contour(
            inside, [0.5], colors=[REGION_COLORS[region_name]], linewidths=1.65,
        )
    for region_name, x, y in channel_cartoon_region_label_points(region_codes):
        ax.text(
            x, y, region_name, ha="center", va="center",
            fontsize=7.2, fontweight="bold", color=REGION_COLORS[region_name],
            zorder=20,
            bbox={
                "boxstyle": "round,pad=0.20", "facecolor": "white",
                "edgecolor": REGION_COLORS[region_name], "linewidth": 0.8,
                "alpha": 0.92,
            },
        )
    if marker == "DAPI":
        ax.text(
            0.016, 0.035, "Tissue boundary", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=6.3, fontweight="bold",
            color=CARTOON_TISSUE_EDGE,
            bbox={
                "boxstyle": "round,pad=0.20", "facecolor": "white",
                "edgecolor": CARTOON_TISSUE_EDGE, "linewidth": 0.7,
                "alpha": 0.90,
            },
            zorder=22,
        )
    add_channel_cartoon_scalebar(ax, mask.shape, um_per_px, 100.0)
    ax.set_title(
        marker, color=CHANNEL_CARTOON_COLORS[marker], fontsize=10.5,
        fontweight="bold", pad=3,
    )
    ax.set_xlim(-0.5, float(mask.shape[1]) - 0.5)
    ax.set_ylim(float(mask.shape[0]) - 0.5, -0.5)
    ax.set_axis_off()


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.045, 1.035, letter, transform=ax.transAxes, fontsize=22, fontweight="bold", ha="right", va="bottom", color="black", clip_on=False)


def select_microglia_inset(
    cells: pd.DataFrame,
    animal: str,
    section: str,
    bounds: tuple[int, int, int, int],
    preferred_uid: str = "",
) -> tuple[float, float, str]:
    x0, y0, x1, y1 = bounds
    subset = cells.loc[cells["animal_id"].astype(str).eq(animal) & cells["section"].astype(str).eq(section)].copy()
    qc_pass = subset.get("cell_qc_pass", pd.Series(False, index=subset.index)).astype(str).str.lower().isin(("true", "1", "yes"))
    subset = subset.loc[
        qc_pass
        & pd.to_numeric(subset["area_px"], errors="coerce").between(50, 450)
        & pd.to_numeric(subset["centroid_x"], errors="coerce").between(x0 + 60, x1 - 60)
        & pd.to_numeric(subset["centroid_y"], errors="coerce").between(y0 + 60, y1 - 60)
    ].copy()
    if subset.empty:
        raise RuntimeError(f"No visible QC-valid inset cell for {animal}/{section}")
    if str(preferred_uid).strip():
        preferred = subset.loc[
            subset["cell_uid"].astype(str).eq(str(preferred_uid).strip())
        ]
        if len(preferred) != 1:
            raise RuntimeError(
                f"Requested panel A cell {preferred_uid!r} is not one visible QC-valid "
                f"cell in {animal}/{section}"
            )
        row = preferred.iloc[0]
        return float(row["centroid_x"]), float(row["centroid_y"]), str(row["cell_uid"])
    coordinates = subset[["centroid_x", "centroid_y"]].apply(
        pd.to_numeric, errors="coerce",
    ).to_numpy(float)
    distances = np.sqrt(np.sum((coordinates[:, None, :] - coordinates[None, :, :]) ** 2, axis=2))
    np.fill_diagonal(distances, np.inf)
    nearest = np.min(distances, axis=1) if len(subset) > 1 else np.full(len(subset), 100.0)
    def percentile_rank(column: str) -> np.ndarray:
        values = pd.to_numeric(subset[column], errors="coerce")
        values = values.fillna(values.median() if values.notna().any() else 0.0)
        return values.rank(method="average", pct=True).to_numpy(float)
    subset["_visibility_score"] = (
        percentile_rank("iba1_mean_intensity")
        + 0.55 * percentile_rank("iba1_max_intensity")
        + 0.45 * pd.Series(nearest).rank(method="average", pct=True).to_numpy(float)
        + 0.20 * percentile_rank("classifier_confidence")
    )
    row = subset.sort_values(
        ["_visibility_score", "cell_uid"], ascending=[False, True],
    ).iloc[0]
    return float(row["centroid_x"]), float(row["centroid_y"]), str(row["cell_uid"])


def select_gfap_inset(gfap: np.ndarray, positive: np.ndarray, labels: np.ndarray, bounds: tuple[int, int, int, int], size: int = 76) -> tuple[float, float]:
    x0, y0, x1, y1 = bounds
    signal = normalize(marker_signal(gfap), high=99.5, gamma=1.0)
    foreground = signal * (positive > 0)
    score = uniform_filter(foreground.astype(float), size=max(9, size // 2), mode="constant")
    valid = np.zeros_like(score, dtype=bool)
    margin = size // 2 + 3
    valid[max(y0, margin):min(y1, score.shape[0] - margin), max(x0, margin):min(x1, score.shape[1] - margin)] = True
    valid &= labels > 0
    masked = np.where(valid, score, -np.inf)
    if not np.isfinite(masked).any():
        return (0.35 * x0 + 0.65 * x1, 0.50 * y0 + 0.50 * y1)
    cy, cx = np.unravel_index(np.argmax(masked), masked.shape)
    return float(cx), float(cy)


def draw_microscopy(
    fig: plt.Figure,
    spec,
    letter: str,
    title: str,
    dapi: np.ndarray,
    marker: np.ndarray,
    marker_color: tuple[float, float, float],
    bounds: tuple[int, int, int, int],
    inset_center: tuple[float, float],
    inset_label: str,
    um_per_px: float,
    inset_size: int = 76,
    source_dapi_sha256: str = "",
    inset_outline_color: str = "white",
    inset_width: str = "23%",
    inset_height: str = "30%",
    source_box_label: str = "",
) -> plt.Axes:
    ax = fig.add_subplot(spec)
    x0, y0, x1, y1 = bounds

    # Preserve the original whole-image normalization while allocating RGB only
    # for the displayed field and inset rather than for the complete source image.
    dapi_signal = marker_signal(dapi)
    marker_data = marker_signal(marker)

    def limits(signal: np.ndarray, low: float, high: float) -> tuple[float, float]:
        finite = np.asarray(signal, dtype=float)
        finite = finite[np.isfinite(finite)]
        if not finite.size:
            return 0.0, 0.0
        lo, hi = np.percentile(finite, (low, high))
        return float(lo), float(hi)

    dapi_lo, dapi_hi = limits(dapi_signal, 1.0, 99.6)
    marker_lo, marker_hi = limits(marker_data, 1.0, 99.7)

    def compose_window(dapi_window: np.ndarray, marker_window: np.ndarray) -> np.ndarray:
        if dapi_hi > dapi_lo:
            d = np.power(np.clip((dapi_window - dapi_lo) / (dapi_hi - dapi_lo), 0.0, 1.0), 0.78)
        else:
            d = np.zeros_like(dapi_window, dtype=float)
        if marker_hi > marker_lo:
            m = np.power(np.clip((marker_window - marker_lo) / (marker_hi - marker_lo), 0.0, 1.0), 0.72)
        else:
            m = np.zeros_like(marker_window, dtype=float)
        rgb = np.zeros((*d.shape, 3), dtype=float)
        rgb[..., 2] += 0.82 * d
        rgb[..., 0] += 0.08 * d
        for channel, coefficient in enumerate(marker_color):
            rgb[..., channel] += coefficient * m
        return np.clip(rgb, 0.0, 1.0)

    view = compose_window(dapi_signal[y0:y1, x0:x1], marker_data[y0:y1, x0:x1])
    ax.imshow(view, interpolation="nearest")
    ax.text(0.025, 0.96, title, transform=ax.transAxes, color="white", fontsize=15, fontweight="bold", va="top", bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.55, "pad": 3})
    add_scalebar(ax, view.shape[1], view.shape[0], 200.0, um_per_px)

    dapi_detail, detail_bounds = crop_square(dapi_signal, inset_center, inset_size)
    marker_detail, _ = crop_square(marker_data, inset_center, inset_size)
    detail = compose_window(dapi_detail, marker_detail)
    iax = inset_axes(ax, width=inset_width, height=inset_height, loc="upper right", borderpad=0.75)
    iax.imshow(detail, interpolation="nearest")
    iax.set_xticks([])
    iax.set_yticks([])
    for spine in iax.spines.values():
        spine.set_color(inset_outline_color)
        spine.set_linewidth(2.8)
    iax.text(0.04, 0.94, inset_label, transform=iax.transAxes, color="white", fontsize=6.7, fontweight="bold", va="top", bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.62, "pad": 1.5})
    add_scalebar(iax, inset_size, inset_size, 25.0, um_per_px, inset=True)
    sx0, sy0, sx1, sy1 = detail_bounds
    source_rect = (sx0 - x0, sy0 - y0, sx1 - sx0, sy1 - sy0)
    ax.add_patch(Rectangle(source_rect[:2], source_rect[2], source_rect[3], fill=False, edgecolor="#101820", linewidth=5.4, zorder=12))
    ax.add_patch(Rectangle(source_rect[:2], source_rect[2], source_rect[3], fill=False, edgecolor=inset_outline_color, linewidth=2.8, zorder=13))
    if source_box_label:
        ax.text(
            source_rect[0] + source_rect[2] / 2.0, source_rect[1] - 8.0,
            source_box_label, color=inset_outline_color, fontsize=7.4,
            fontweight="bold", ha="center", va="bottom", zorder=14,
            bbox={"facecolor": "#101820", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
        )
    ax.set_axis_off()
    ax.set_anchor("C")
    panel_label(ax, letter)
    ax._source_bounds_px = tuple(int(value) for value in bounds)
    ax._display_shape_px = tuple(int(value) for value in view.shape[:2])
    ax._source_dapi_sha256 = str(source_dapi_sha256)
    return ax


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str = "#24414b", width: float = 1.5) -> None:
    ax.add_patch(FancyArrowPatch(start, end, transform=ax.transAxes, arrowstyle="-|>", mutation_scale=11, linewidth=width, color=color))


def draw_feature_stack(ax: plt.Axes, x: float, height: float, color: str, channels: str, spatial: str) -> None:
    width = 0.055
    base_y = 0.52 - height / 2
    for offset in (0.018, 0.012, 0.006, 0.0):
        ax.add_patch(Rectangle((x + offset, base_y + offset), width, height, transform=ax.transAxes, facecolor=color, edgecolor="#173844", linewidth=1.0, alpha=0.92))
    ax.text(x + width / 2 + 0.009, base_y - 0.08, channels, transform=ax.transAxes, ha="center", fontsize=10.2, fontweight="bold", color="#173844")
    ax.text(x + width / 2 + 0.009, base_y - 0.145, spatial, transform=ax.transAxes, ha="center", fontsize=8.6, color="#38545d")



def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def validate_all_hil_analysis(
    analysis_dir: Path,
    manual_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Refuse statistics produced without complete native 61/61 HIL provenance."""
    coverage_path = analysis_dir / "human_review_coverage.csv"
    discovered_path = analysis_dir / "discovered_sections.csv"
    summary_path = analysis_dir / "per_section_region_summary.csv"
    if not coverage_path.is_file():
        raise RuntimeError(f"Missing native HIL coverage audit: {coverage_path}")
    if not discovered_path.is_file():
        raise RuntimeError(f"Missing current GFAP/Iba1 section inventory: {discovered_path}")
    if not summary_path.is_file():
        raise RuntimeError(f"Missing per-section analysis provenance: {summary_path}")
    coverage = pd.read_csv(coverage_path, dtype=str, keep_default_na=False)
    coverage_required = {
        "sample", "section", "review_status", "accepted_source",
        "manual_region_path", "manual_region_sha256",
        "manual_region_review_path", "manual_region_review_sha256",
        "reviewer", "accepted_at", "source_dapi_path", "source_dapi_sha256",
        "source_dapi_sha256_verified",
    }
    missing = sorted(coverage_required.difference(coverage.columns))
    if missing:
        raise RuntimeError("HIL coverage lacks provenance columns: " + ", ".join(missing))
    keys = list(zip(coverage["sample"].astype(str), coverage["section"].astype(str)))
    invalid_coverage = coverage.loc[
        coverage["review_status"].ne("accepted")
        | ~coverage["accepted_source"].isin(ALLOWED_NATIVE_ACCEPTED_SOURCES)
        | ~coverage["source_dapi_sha256_verified"].map(truthy)
        | coverage["manual_region_path"].eq("")
        | coverage["manual_region_sha256"].eq("")
        | coverage["manual_region_review_path"].eq("")
        | coverage["manual_region_review_sha256"].eq("")
        | coverage["reviewer"].str.strip().eq("")
        | coverage["accepted_at"].str.strip().eq("")
        | coverage["source_dapi_path"].eq("")
        | coverage["source_dapi_sha256"].eq("")
    ]
    if len(coverage) != EXPECTED_HUMAN_REVIEW_SECTIONS or len(set(keys)) != EXPECTED_HUMAN_REVIEW_SECTIONS or not invalid_coverage.empty:
        raise RuntimeError(
            "Publication requires exactly 61 unique, valid native HIL reviews; "
            f"found rows={len(coverage)}, unique={len(set(keys))}, invalid={len(invalid_coverage)}"
        )

    discovered = pd.read_csv(discovered_path, dtype=str, keep_default_na=False)
    discovered_required = {"sample", "section", "dapi_image"}
    missing = sorted(discovered_required.difference(discovered.columns))
    if missing:
        raise RuntimeError("Current section inventory lacks columns: " + ", ".join(missing))
    discovered_keys = list(zip(discovered["sample"], discovered["section"]))
    if (
        len(discovered) != EXPECTED_HUMAN_REVIEW_SECTIONS
        or len(set(discovered_keys)) != EXPECTED_HUMAN_REVIEW_SECTIONS
        or set(discovered_keys) != set(keys)
    ):
        raise RuntimeError(
            "HIL coverage does not match the exact current 61-section inventory: "
            f"inventory rows={len(discovered)}, inventory unique={len(set(discovered_keys))}, "
            f"coverage unique={len(set(keys))}"
        )
    dapi_by_key = {
        (str(row["sample"]), str(row["section"])): Path(str(row["dapi_image"])).expanduser().resolve()
        for row in discovered.to_dict("records")
    }

    # Re-open every accepted native receipt and hash its TIFF, JSON, and DAPI.
    # This makes the publication gate independent of a stale/tampered coverage CSV.
    manual_root = manual_root.expanduser().resolve()
    verified_coverage: dict[tuple[str, str], dict[str, Any]] = {}
    for coverage_index, row in coverage.iterrows():
        row = row.to_dict()
        sample = str(row["sample"])
        section = str(row["section"])
        dapi_path = Path(str(row["source_dapi_path"])).expanduser().resolve()
        if dapi_path != dapi_by_key[(sample, section)]:
            raise RuntimeError(
                f"Coverage DAPI is not the current inventory DAPI for {sample}/{section}: "
                f"{dapi_path} != {dapi_by_key[(sample, section)]}"
            )
        labels, provenance = load_native_accepted_arc_me(
            manual_root, sample, section, dapi_path
        )
        if not _path_is_within(provenance["mask_path"], manual_root):
            raise RuntimeError(
                f"Accepted mask lies outside native manual root: {provenance['mask_path']}"
            )
        expected = {
            "manual_region_path": str(provenance["mask_path"]),
            "manual_region_sha256": provenance["mask_sha256"],
            "manual_region_review_path": str(provenance["review_path"]),
            "manual_region_review_sha256": provenance["review_sha256"],
            "accepted_source": provenance["accepted_source"],
            "reviewer": provenance["reviewer"],
            "accepted_at": provenance["accepted_at"],
            "source_dapi_sha256": provenance["source_dapi_sha256"],
        }
        mismatches = {
            field: (str(row[field]), str(value))
            for field, value in expected.items()
            if str(row[field]) != str(value)
        }
        if mismatches:
            raise RuntimeError(
                f"Coverage/receipt provenance mismatch for {sample}/{section}: {mismatches}"
            )
        final_mask_path, final_geometry = verify_final_analysis_mask_identity(
            analysis_dir, sample, section, labels, provenance
        )
        coverage.loc[coverage_index, "final_analysis_mask_path"] = str(final_mask_path)
        coverage.loc[coverage_index, "final_analysis_mask_sha256"] = sha256_file(final_mask_path)
        coverage.loc[coverage_index, "final_analysis_mask_pixel_identical"] = "true"
        geometry_path = final_mask_path.parent / "region_geometry.json"
        coverage.loc[coverage_index, "final_analysis_geometry_path"] = str(
            geometry_path.resolve()
        )
        coverage.loc[coverage_index, "final_analysis_geometry_sha256"] = sha256_file(
            geometry_path
        )
        coverage.loc[coverage_index, "final_analysis_region_source"] = str(
            final_geometry.get("region_source", "")
        )
        verified_coverage[(sample, section)] = provenance

    summary = pd.read_csv(summary_path, dtype=str, keep_default_na=False)
    summary_required = {
        "sample", "section", "region", "region_source",
        "manual_region_review_status", "manual_region_accepted_source",
        "manual_region_path", "manual_region_sha256",
        "manual_region_source_dapi_sha256_verified", "manual_region_human_in_the_loop",
    }
    missing = sorted(summary_required.difference(summary.columns))
    if missing:
        raise RuntimeError(
            "Analysis is stale or automated; per-section table lacks native HIL columns: "
            + ", ".join(missing)
        )
    summary_keys = list(zip(summary["sample"].astype(str), summary["section"].astype(str)))
    unique_summary_keys = set(summary_keys)
    invalid_summary = summary.loc[
        summary["region_source"].ne("manual_region_tiff")
        | summary["manual_region_review_status"].ne("accepted")
        | ~summary["manual_region_accepted_source"].isin(ALLOWED_NATIVE_ACCEPTED_SOURCES)
        | ~summary["manual_region_source_dapi_sha256_verified"].map(truthy)
        | ~summary["manual_region_human_in_the_loop"].map(truthy)
        | summary["manual_region_path"].eq("")
        | summary["manual_region_sha256"].eq("")
    ]
    bad_region_sets = [
        f"{sample}/{section}"
        for (sample, section), group in summary.groupby(["sample", "section"], sort=False)
        if set(group["region"].astype(str)) != {
            name for name, status in verified_coverage[(str(sample), str(section))]
            ["region_status"].items() if status == "drawn"
        }
    ]
    if (
        len(unique_summary_keys) != EXPECTED_HUMAN_REVIEW_SECTIONS
        or set(unique_summary_keys) != set(verified_coverage)
        or not invalid_summary.empty
        or bad_region_sets
    ):
        raise RuntimeError(
            "Publication refuses non-HIL/stale analysis: "
            f"rows={len(summary)}, unique sections={len(unique_summary_keys)}, "
            f"invalid rows={len(invalid_summary)}, "
            f"bad drawn-region section sets={len(bad_region_sets)}"
        )
    for (sample, section), group in summary.groupby(["sample", "section"], sort=False):
        provenance = verified_coverage[(str(sample), str(section))]
        expected = {
            "manual_region_path": str(provenance["mask_path"]),
            "manual_region_sha256": provenance["mask_sha256"],
            "manual_region_accepted_source": provenance["accepted_source"],
        }
        for field, value in expected.items():
            actual_values = set(group[field].astype(str))
            if actual_values != {str(value)}:
                raise RuntimeError(
                    f"Analysis/receipt provenance mismatch for {sample}/{section} "
                    f"field {field}: {sorted(actual_values)} != {value!r}"
                )
    return coverage, summary


def _state_token(value: str) -> str:
    return "".join(character for character in str(value).lower() if character.isalnum())


def consistent_group_metadata(
    frame: pd.DataFrame,
    group_columns: Sequence[str],
    metadata_columns: Sequence[str],
    *,
    context: str,
) -> pd.DataFrame:
    """Resolve one optional text value per group while rejecting conflicts."""
    rows: list[dict[str, Any]] = []
    for keys, group in frame.groupby(
        list(group_columns), dropna=False, sort=False,
    ):
        key_values = keys if isinstance(keys, tuple) else (keys,)
        record: dict[str, Any] = dict(zip(group_columns, key_values))
        for column in metadata_columns:
            values: list[str] = []
            for value in group[column]:
                if pd.isna(value):
                    continue
                text = str(value).strip()
                if text and text not in values:
                    values.append(text)
            if len(values) > 1:
                raise RuntimeError(
                    f"Conflicting {column} values for {context} {keys}: {values}"
                )
            record[column] = values[0] if values else pd.NA
        rows.append(record)
    return pd.DataFrame(rows, columns=[*group_columns, *metadata_columns])


def reaggregate_animal_regions(
    section_summary: pd.DataFrame,
    cells: pd.DataFrame,
) -> pd.DataFrame:
    """Rebuild every animal/region endpoint used by Figure 5 from detail tables."""
    group_columns = ["animal_id", "animal_key", "region"]
    sum_columns = [
        "region_area_px", "region_area_um2", "dapi_nuclei",
        "gfap_positive_area_px", "gfap_positive_area_um2",
        "gfap_corrected_integrated_intensity", "iba1_positive_area_px",
        "iba1_corrected_integrated_intensity", "iba1_assay_area_um2",
        "iba1_objects_raw", "microglia_cells_qc_accepted",
        "iba1_objects_qc_rejected",
    ]
    metadata_columns = ["cond", "genotype", "sex", "cage"]
    required_section = set(group_columns + sum_columns + metadata_columns + [
        "sample", "section", "iba1_segmentation_valid",
    ])
    missing = sorted(required_section.difference(section_summary.columns))
    if missing:
        raise RuntimeError(
            "Cannot reaggregate Figure 5: section table lacks " + ", ".join(missing)
        )
    sections = section_summary.copy()
    sections = sections.loc[
        sections["region"].astype(str).isin(ALL_HIL_REGIONS)
    ].copy()
    if sections.empty:
        raise RuntimeError("Cannot reaggregate Figure 5: no drawn ARC/ME/VMN section rows")
    for column in group_columns:
        sections[column] = sections[column].astype(str)
    for column in sum_columns:
        numeric = pd.to_numeric(sections[column], errors="coerce")
        if numeric.isna().any():
            raise RuntimeError(f"Cannot reaggregate Figure 5: nonnumeric {column}")
        sections[column] = numeric
    duplicated = sections.duplicated(["sample", "section", "region"], keep=False)
    if duplicated.any():
        identities = sorted({
            "/".join(map(str, values))
            for values in sections.loc[
                duplicated, ["sample", "section", "region"]
            ].itertuples(index=False, name=None)
        })
        raise RuntimeError(
            "Section table has duplicate sample/section/region rows: "
            + ", ".join(identities[:8])
        )

    metadata = consistent_group_metadata(
        sections,
        ["animal_id", "animal_key"],
        metadata_columns,
        context="Figure 5 animal",
    )
    aggregated = (
        sections.groupby(group_columns, dropna=False, sort=True)[sum_columns]
        .sum().reset_index()
    )
    section_counts = (
        sections.groupby(group_columns, dropna=False, sort=True)
        .size().rename("n_sections").reset_index()
    )
    valid_sections = (
        sections.assign(
            _iba1_valid=sections["iba1_segmentation_valid"].map(truthy).astype(int)
        )
        .groupby(group_columns, dropna=False, sort=True)["_iba1_valid"]
        .sum().rename("n_iba1_valid_sections").reset_index()
    )
    aggregated = aggregated.merge(
        section_counts, on=group_columns, validate="one_to_one"
    ).merge(valid_sections, on=group_columns, validate="one_to_one")

    required_cells = {
        "animal_id", "animal_key", "region", "microglia_state",
        "classifier_uncertain", "cell_qc_pass",
    }
    missing = sorted(required_cells.difference(cells.columns))
    if missing:
        raise RuntimeError(
            "Cannot reaggregate Figure 5: cell table lacks " + ", ".join(missing)
        )
    accepted_cells = cells.copy()
    if not accepted_cells.empty and not accepted_cells["cell_qc_pass"].map(truthy).all():
        raise RuntimeError(
            "Publication cell table contains objects that failed biological-cell QC"
        )
    accepted_cells = accepted_cells.loc[
        accepted_cells["region"].astype(str).isin(ALL_HIL_REGIONS)
    ].copy()
    for column in group_columns:
        accepted_cells[column] = accepted_cells[column].astype(str)

    if accepted_cells.empty:
        aggregated["microglia_cells"] = 0
        for state in STATE_COLORS:
            aggregated[f"microglia_{_state_token(state)}_cells"] = 0
        aggregated["microglia_uncertain_cells"] = 0
    else:
        total = (
            accepted_cells.groupby(group_columns, dropna=False, sort=True)
            .size().rename("microglia_cells").reset_index()
        )
        aggregated = aggregated.merge(
            total, on=group_columns, how="left", validate="one_to_one"
        )
        for state in STATE_COLORS:
            column = f"microglia_{_state_token(state)}_cells"
            state_counts = (
                accepted_cells.loc[
                    accepted_cells["microglia_state"].astype(str).eq(state)
                ]
                .groupby(group_columns, dropna=False, sort=True)
                .size().rename(column).reset_index()
            )
            aggregated = aggregated.merge(
                state_counts, on=group_columns, how="left", validate="one_to_one"
            )
        uncertain = (
            accepted_cells.loc[accepted_cells["classifier_uncertain"].map(truthy)]
            .groupby(group_columns, dropna=False, sort=True)
            .size().rename("microglia_uncertain_cells").reset_index()
        )
        aggregated = aggregated.merge(
            uncertain, on=group_columns, how="left", validate="one_to_one"
        )
    count_columns = [
        "microglia_cells", "microglia_uncertain_cells",
        *(f"microglia_{_state_token(state)}_cells" for state in STATE_COLORS),
    ]
    for column in count_columns:
        if column not in aggregated:
            aggregated[column] = 0
        aggregated[column] = pd.to_numeric(
            aggregated[column], errors="coerce"
        ).fillna(0).astype(int)

    if not np.array_equal(
        aggregated["microglia_cells"].to_numpy(int),
        aggregated["microglia_cells_qc_accepted"].to_numpy(int),
    ):
        raise RuntimeError(
            "Cell table counts do not match section-level QC-accepted microglia counts"
        )
    area = aggregated["region_area_um2"].replace(0, np.nan)
    area_px = aggregated["region_area_px"].replace(0, np.nan)
    assay_area = aggregated["iba1_assay_area_um2"].replace(0, np.nan)
    aggregated["dapi_density_per_um2"] = aggregated["dapi_nuclei"] / area
    aggregated["gfap_area_fraction"] = aggregated["gfap_positive_area_px"] / area_px
    aggregated["gfap_corrected_intensity_per_um2"] = (
        aggregated["gfap_corrected_integrated_intensity"] / area
    )
    aggregated["gfap_corrected_mean_intensity"] = (
        aggregated["gfap_corrected_integrated_intensity"] / area_px
    )
    aggregated["gfap_reactivity_index"] = (
        aggregated["gfap_area_fraction"] * aggregated["gfap_corrected_mean_intensity"]
    )
    aggregated["iba1_area_fraction"] = aggregated["iba1_positive_area_px"] / area_px
    aggregated["iba1_corrected_intensity_per_um2"] = (
        aggregated["iba1_corrected_integrated_intensity"] / area
    )
    aggregated["microglia_density_per_um2"] = (
        aggregated["microglia_cells"] / assay_area
    )
    aggregated["microglia_density_per_100k_um2"] = (
        aggregated["microglia_density_per_um2"] * MICROGLIA_DENSITY_SCALE_UM2
    )
    for state in STATE_COLORS:
        token = _state_token(state)
        cell_column = f"microglia_{token}_cells"
        aggregated[f"microglia_{token}_fraction"] = (
            aggregated[cell_column] / aggregated["microglia_cells"].replace(0, np.nan)
        )
        aggregated[f"microglia_{token}_density_per_um2"] = (
            aggregated[cell_column] / assay_area
        )
        aggregated[f"microglia_{token}_density_per_100k_um2"] = (
            aggregated[f"microglia_{token}_density_per_um2"]
            * MICROGLIA_DENSITY_SCALE_UM2
        )
    aggregated["activated_microglia_cells"] = (
        aggregated["microglia_activated_cells"] + aggregated["microglia_amoeboid_cells"]
    )
    aggregated["resting_microglia_cells"] = aggregated["microglia_ramified_cells"]
    aggregated["activated_microglia_fraction"] = (
        aggregated["activated_microglia_cells"]
        / aggregated["microglia_cells"].replace(0, np.nan)
    )
    aggregated["resting_microglia_fraction"] = (
        aggregated["resting_microglia_cells"]
        / aggregated["microglia_cells"].replace(0, np.nan)
    )
    aggregated["activated_microglia_density_per_um2"] = (
        aggregated["activated_microglia_cells"] / assay_area
    )
    aggregated["activated_microglia_density_per_100k_um2"] = (
        aggregated["activated_microglia_density_per_um2"]
        * MICROGLIA_DENSITY_SCALE_UM2
    )
    return (
        aggregated.merge(
            metadata, on=["animal_id", "animal_key"], how="left", validate="many_to_one"
        )
        .sort_values(group_columns).reset_index(drop=True)
    )


def verify_canonical_animal_summary(
    canonical: pd.DataFrame,
    recomputed: pd.DataFrame,
) -> None:
    """Fail unless the canonical animal table equals detail-table reaggregation."""
    keys = ["animal_id", "animal_key", "region"]
    for name, frame in (("canonical", canonical), ("recomputed", recomputed)):
        missing = sorted(set(keys).difference(frame.columns))
        if missing:
            raise RuntimeError(f"{name} animal summary lacks keys: {', '.join(missing)}")
        if frame.duplicated(keys, keep=False).any():
            raise RuntimeError(f"{name} animal summary has duplicate animal/region keys")
    left = canonical.copy()
    right = recomputed.copy()
    for column in keys:
        left[column] = left[column].astype(str)
        right[column] = right[column].astype(str)
    left = left.sort_values(keys).reset_index(drop=True)
    right = right.sort_values(keys).reset_index(drop=True)
    if list(map(tuple, left[keys].to_numpy())) != list(map(tuple, right[keys].to_numpy())):
        raise RuntimeError(
            "Canonical animal summary keys do not match section/cell reaggregation"
        )
    missing = sorted(set(right.columns).difference(left.columns))
    if missing:
        raise RuntimeError(
            "Canonical animal summary lacks independently reconstructed columns: "
            + ", ".join(missing)
        )
    text_columns = {"cond", "genotype", "sex", "cage"}
    disagreements: list[str] = []
    for column in (name for name in right.columns if name not in keys):
        if column in text_columns:
            mismatch = left[column].fillna("").astype(str).ne(
                right[column].fillna("").astype(str)
            ).to_numpy(bool)
        else:
            observed = pd.to_numeric(left[column], errors="coerce").to_numpy(float)
            expected = pd.to_numeric(right[column], errors="coerce").to_numpy(float)
            mismatch = ~np.isclose(
                observed, expected, rtol=1e-10, atol=1e-8, equal_nan=True
            )
        if bool(np.any(mismatch)):
            disagreements.append(column)
    if disagreements:
        raise RuntimeError(
            "Canonical animal summary is stale or inconsistent with current section/cell data: "
            + ", ".join(disagreements)
        )


def load_native_accepted_arc_me(
    manual_root: Path,
    animal: str,
    section: str,
    dapi_path: Path,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Load one native accepted mask/receipt; registration-derived receipts are invalid."""
    for value, field in ((animal, "animal"), (section, "section")):
        if not str(value).strip() or any(token in str(value) for token in ("/", "\\", "\x00")):
            raise RuntimeError(f"Unsafe representative {field}: {value!r}")
    candidates = _native_receipt_candidates(manual_root, animal, section)
    _current_animal, current_mask, current_review = candidates[0]
    if current_mask.is_file() != current_review.is_file():
        missing = current_review if current_mask.is_file() else current_mask
        raise RuntimeError(
            f"Incomplete exact-source accepted pair for {animal}/{section}; "
            f"missing {missing}"
        )
    selected = next(
        (
            (receipt_animal, mask, review)
            for receipt_animal, mask, review in candidates
            if mask.is_file() and review.is_file()
        ),
        None,
    )
    if selected is None:
        expected = "; ".join(
            f"{mask} + {review}" for _receipt_animal, mask, review in candidates
        )
        raise RuntimeError(
            f"Accepted native ARC/ME/VMN TIFF+JSON pair is required for "
            f"{animal}/{section}; expected {expected}"
        )
    receipt_animal_path, mask_path, review_path = selected
    try:
        review = json.loads(review_path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise RuntimeError(f"Cannot parse accepted review receipt {review_path}: {exc}") from exc
    if not isinstance(review, dict) or review.get("accepted") is not True:
        raise RuntimeError(f"Receipt is not explicitly HIL accepted: {review_path}")
    accepted_source = str(review.get("accepted_source", "")).strip()
    if accepted_source not in ALLOWED_NATIVE_ACCEPTED_SOURCES:
        raise RuntimeError(
            f"Representative receipt is not native HIL (accepted_source={accepted_source!r}): {review_path}"
        )
    if str(review.get("dataset", "")).strip().lower() != "gfap":
        raise RuntimeError(f"Representative receipt requires dataset=gfap: {review_path}")
    receipt_animal = str(review.get("animal", "")).strip()
    allowed_animals = {str(animal), _receipt_base_animal(str(animal)), receipt_animal_path}
    if receipt_animal not in allowed_animals or str(review.get("section", "")).strip() != section:
        raise RuntimeError(f"Representative receipt identity does not match {animal}/{section}: {review_path}")
    reviewer = str(review.get("reviewer", "")).strip()
    accepted_at = str(review.get("accepted_at", "")).strip()
    if not reviewer or not accepted_at:
        raise RuntimeError(f"Representative receipt lacks reviewer/accepted_at: {review_path}")
    try:
        dt.datetime.fromisoformat(accepted_at.replace("Z", "+00:00"))
    except ValueError as exc:
        raise RuntimeError(f"Representative accepted_at is not ISO-8601: {review_path}") from exc
    receipt_codes = review.get("label_codes")
    if receipt_codes not in ({"background": 0, "ARC": 1, "ME": 2},
                              {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}):
        raise RuntimeError(f"Representative receipt has invalid ARC/ME/VMN codes: {review_path}")
    if not dapi_path.is_file():
        raise RuntimeError(f"Representative DAPI is missing: {dapi_path}")
    dapi_shape = marker_signal(tifffile.imread(dapi_path)).shape
    if (int(review.get("image_height_px", -1)), int(review.get("image_width_px", -1))) != tuple(dapi_shape):
        raise RuntimeError(f"Representative receipt dimensions do not match DAPI: {review_path}")
    dapi_sha256 = sha256_file(dapi_path)
    if str(review.get("source_dapi_sha256", "")).strip().lower() != dapi_sha256:
        raise RuntimeError(f"Representative receipt DAPI SHA-256 mismatch: {review_path}")
    mask_sha256 = sha256_file(mask_path)
    if str(review.get("mask_sha256", "")).strip().lower() != mask_sha256:
        raise RuntimeError(f"Representative accepted mask SHA-256 mismatch: {review_path}")
    labels = np.squeeze(np.asarray(tifffile.imread(mask_path)))
    if labels.ndim != 2 or labels.shape != dapi_shape:
        raise RuntimeError(f"Representative accepted mask shape {labels.shape} != DAPI {dapi_shape}: {mask_path}")
    if not np.issubdtype(labels.dtype, np.number) or not np.all(np.isfinite(labels)) or not np.all(labels == np.rint(labels)):
        raise RuntimeError(f"Representative accepted mask is not finite integer-valued: {mask_path}")
    labels = labels.astype(np.uint8, copy=False)
    codes = {int(value) for value in np.unique(labels)}
    if not codes.issubset({0, 1, 2, 3}) or not codes.intersection({1, 2, 3}):
        raise RuntimeError(f"Representative accepted mask violates 0/1/2/3 nonempty contract: {mask_path}")
    if receipt_codes == {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}:
        raw_status = review.get("region_status")
        reasons = review.get("region_absence_reasons")
        if not isinstance(raw_status, dict) or not isinstance(reasons, dict):
            raise RuntimeError(f"Schema-v3 receipt lacks region status/absence maps: {review_path}")
        region_status = {
            name: str(raw_status.get(name, ""))
            for name in ("ARC", "ME", "VMN")
        }
    else:
        legacy_absent = review.get("me_absent") is True
        region_status = {
            "ARC": "drawn" if 1 in codes else "not_delineated",
            "ME": "drawn" if 2 in codes else "explicitly_absent" if legacy_absent else "not_delineated",
            "VMN": "not_delineated",
        }
        reasons = {
            "ARC": "",
            "ME": str(review.get("me_absent_reason", "")).strip(),
            "VMN": "",
        }
    for name, code in REGION_CODES.items():
        area = int(np.count_nonzero(labels == code))
        state = region_status[name]
        reason = str(reasons.get(name, "")).strip()
        if state not in {"drawn", "explicitly_absent", "not_delineated"}:
            raise RuntimeError(f"Representative {name} status is invalid: {state!r}")
        if (state == "drawn") != (area > 0):
            raise RuntimeError(f"Representative {name} status conflicts with mask area")
        if state == "explicitly_absent" and not reason:
            raise RuntimeError(f"Representative explicit {name} absence lacks a reason")
        if state != "explicitly_absent" and reason:
            raise RuntimeError(f"Representative {name} reason conflicts with status {state}")
    me_absent = region_status["ME"] == "explicitly_absent"
    me_absent_reason = str(reasons.get("ME", "")).strip()
    source_paths = review.get("source_paths") if isinstance(review.get("source_paths"), dict) else {}
    receipt_dapi_text = str(source_paths.get("dapi", "") or "").strip()
    if not _review_provenance_path_matches(receipt_dapi_text, dapi_path):
        raise RuntimeError(f"Representative receipt source_paths.dapi mismatch: {review_path}")
    if (
        source_paths.get("manifest")
        or review.get("source_manifest_sha256")
        or review.get("registration") is not None
    ):
        raise RuntimeError(f"Registration/manifest-derived representative receipt is forbidden: {review_path}")
    proposal_verified = False
    if accepted_source == "automated_mask_unchanged":
        proposal_text = str(source_paths.get("automated_mask", "") or "").strip()
        proposal_candidate = Path(proposal_text).expanduser() if proposal_text else None
        proposal_path = (
            proposal_candidate.resolve()
            if proposal_candidate is not None and proposal_candidate.is_absolute()
            else (PAPER_ROOT / proposal_candidate).resolve()
            if proposal_candidate is not None
            else None
        )
        proposal_sha256 = str(
            review.get("automated_mask_sha256")
            or review.get("proposal_mask_sha256")
            or ""
        ).strip().lower()
        confirmation = (
            review.get("automated_mask_confirmation")
            or review.get("proposal_confirmation")
        )
        automated_provenance = (
            review.get("automated_mask_provenance")
            if isinstance(review.get("automated_mask_provenance"), dict)
            else {}
        )
        if (
            confirmation != "reviewed_image_by_image"
            or proposal_path is None or not proposal_path.is_file()
            or sha256_file(proposal_path) != proposal_sha256
            or str(automated_provenance.get("dataset", "")).strip().lower() != "gfap"
            or str(automated_provenance.get("kind", "")).strip()
            != "native_gfap_automated_mask"
        ):
            raise RuntimeError(f"Native automated proposal provenance is invalid: {review_path}")
        proposal_labels = np.squeeze(np.asarray(tifffile.imread(proposal_path)))
        if (
            proposal_labels.shape != labels.shape
            or not np.issubdtype(proposal_labels.dtype, np.number)
            or not np.all(np.isfinite(proposal_labels))
            or not np.all(proposal_labels == np.rint(proposal_labels))
            or not np.array_equal(proposal_labels.astype(np.uint8, copy=False), labels)
        ):
            raise RuntimeError(
                f"Accepted TIFF is not pixel-identical to unchanged native automated mask: {review_path}"
            )
        proposal_verified = True
    provenance = {
        "animal": animal,
        "section": section,
        "reviewer": reviewer,
        "accepted_at": accepted_at,
        "accepted_source": accepted_source,
        "me_absent": me_absent,
        "me_absent_reason": me_absent_reason,
        "region_status": region_status,
        "mask_path": mask_path.resolve(),
        "review_path": review_path.resolve(),
        "mask_sha256": mask_sha256,
        "review_sha256": sha256_file(review_path),
        "source_dapi_sha256": dapi_sha256,
        "proposal_verified": proposal_verified,
        "label_content_sha256": hashlib.sha256(np.ascontiguousarray(labels).tobytes()).hexdigest(),
    }
    return labels, provenance


def verify_final_analysis_mask_identity(
    analysis_dir: Path,
    animal: str,
    section: str,
    accepted_labels: np.ndarray,
    provenance: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    processed_dir = analysis_dir / "processed_sections" / animal / section
    final_mask_path = processed_dir / f"{animal}_{section}_ARC_ME_labels.tif"
    geometry_path = processed_dir / "region_geometry.json"
    if not final_mask_path.is_file() or not geometry_path.is_file():
        raise RuntimeError(f"Final analysis mask/provenance missing for {animal}/{section}: {processed_dir}")
    final_labels = np.squeeze(np.asarray(tifffile.imread(final_mask_path)))
    if final_labels.shape != accepted_labels.shape or not np.array_equal(final_labels, accepted_labels):
        raise RuntimeError(
            f"Final analysis ARC/ME mask is not pixel-identical to accepted HIL mask for {animal}/{section}"
        )
    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    expected = {
        "region_source": "manual_region_tiff",
        "manual_region_review_status": "accepted",
        "manual_region_accepted_source": provenance["accepted_source"],
        "manual_region_path": str(provenance["mask_path"]),
        "manual_region_sha256": provenance["mask_sha256"],
        "manual_region_review_path": str(provenance["review_path"]),
        "manual_region_review_sha256": provenance["review_sha256"],
        "manual_region_source_dapi_sha256": provenance["source_dapi_sha256"],
    }
    expected.update({
        f"manual_region_{region.lower()}_status": status
        for region, status in provenance["region_status"].items()
    })
    mismatches = {key: (geometry.get(key), value) for key, value in expected.items() if geometry.get(key) != value}
    if (
        mismatches
        or not truthy(geometry.get("manual_region_source_dapi_sha256_verified"))
        or not truthy(geometry.get("manual_region_human_in_the_loop"))
    ):
        raise RuntimeError(f"Final analysis provenance does not match accepted mask for {animal}/{section}: {mismatches}")
    return final_mask_path.resolve(), geometry


def _arc_me_codes(labels: np.ndarray, *, allow_me_absent: bool = False) -> dict[str, int]:
    del allow_me_absent
    positive_codes = {int(value) for value in np.unique(np.asarray(labels)) if int(value) > 0}
    unexpected = sorted(positive_codes.difference(set(REGION_CODES.values())))
    if not positive_codes or unexpected:
        raise RuntimeError(
            "ARC/ME/VMN label TIFF violates the canonical code contract; "
            f"positive={sorted(positive_codes)}, unexpected={unexpected}"
        )
    return dict(REGION_CODES)


def load_region_geometry_arrays(processed_dir: Path) -> dict[str, np.ndarray]:
    path = processed_dir / "region_geometry_arrays.npz"
    if not path.is_file():
        return {}
    with np.load(path, allow_pickle=False) as archive:
        return {name: np.asarray(archive[name]) for name in archive.files}


def label_lumen_center_x(labels: np.ndarray) -> float:
    """Infer the medial tissue-free gap from the accepted label mask."""
    accepted = np.asarray(labels)
    support = accepted > 0
    height, width = support.shape
    mid = width // 2
    candidate_midpoints: list[float] = []
    positive_rows = np.flatnonzero(np.any(support, axis=1))
    for y in positive_rows.tolist():
        left = np.flatnonzero(support[y, :mid])
        right = np.flatnonzero(support[y, mid:])
        if not left.size or not right.size:
            continue
        left_edge = int(left[-1])
        right_edge = int(mid + right[0])
        if right_edge - left_edge >= 2:
            candidate_midpoints.append(0.5 * (left_edge + right_edge))
    if candidate_midpoints:
        return float(np.median(candidate_midpoints))
    support_x = np.flatnonzero(np.any(support, axis=0))
    return float(np.median(support_x)) if support_x.size else float(mid)


def draw_anatomy_qc(
    fig: plt.Figure,
    spec,
    dapi: np.ndarray,
    labels: np.ndarray,
    bounds: tuple[int, int, int, int],
    provenance: dict[str, Any],
    um_per_px: float,
) -> plt.Axes:
    """Show the accepted ARC/ME mask over the exact native microscopy field."""
    ax = fig.add_subplot(spec)
    accepted = np.asarray(labels)
    codes = _arc_me_codes(accepted, allow_me_absent=bool(provenance.get("me_absent")))
    x0, y0, x1, y1 = (int(value) for value in bounds)
    dapi_signal = marker_signal(dapi)
    if accepted.shape != dapi_signal.shape:
        raise RuntimeError(
            f"Accepted mask/DAPI shape mismatch in publication panel: {accepted.shape} != {dapi_signal.shape}"
        )
    dapi_view = dapi_signal[y0:y1, x0:x1]
    accepted_view = accepted[y0:y1, x0:x1]
    if dapi_view.shape != accepted_view.shape or not dapi_view.size:
        raise RuntimeError(f"Native/accepted panel crop mismatch or empty bounds: {bounds}")
    signal = normalize(dapi_view, high=99.8, gamma=0.85)
    rgb = np.zeros((*signal.shape, 3), dtype=float)
    rgb[..., 0] = 0.05 * signal
    rgb[..., 1] = 0.12 * signal
    rgb[..., 2] = 0.92 * signal
    display_regions = tuple(
        region for region in ("ARC", "ME", "VMN")
        if provenance.get("region_status", {}).get(region) == "drawn"
    )
    for region in display_regions:
        color = np.asarray(to_rgb(REGION_COLORS[region]), dtype=float)
        region_mask = accepted_view == codes[region]
        rgb[region_mask] = 0.62 * rgb[region_mask] + 0.38 * color
    ax.imshow(rgb, interpolation="nearest")
    add_scalebar(
        ax, accepted_view.shape[1], accepted_view.shape[0], 200.0, um_per_px,
    )
    for region in display_regions:
        region_mask = accepted_view == codes[region]
        if np.any(region_mask):
            # A dark underlay keeps the colored anatomical boundary visible over
            # bright DAPI while preserving the exact accepted mask geometry.
            ax.contour(
                region_mask, [0.5], colors=["black"], linewidths=6.4, alpha=0.92,
            )
            ax.contour(
                region_mask, [0.5], colors=[REGION_COLORS[region]], linewidths=3.5,
            )
    ax.text(
        0.025, 0.97, "HIL-accepted ARC/ME/VMN", transform=ax.transAxes,
        ha="left", va="top", fontsize=10.0, fontweight="bold", color="white",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.70, "pad": 2.5},
    )
    accepted_date = str(provenance["accepted_at"]).split("T", 1)[0]
    source_label = (
        "manual polygons" if provenance["accepted_source"] == "manual_replacement_polygons"
        else "native auto mask, HIL-confirmed"
    )
    ax.text(
        0.025, 0.865,
        f"Accepted: {accepted_date}\n{source_label}",
        transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color="white",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.64, "pad": 2.0},
    )
    handles = [
        Rectangle((0, 0), 1, 1, facecolor=REGION_COLORS[region], edgecolor="none", label=region)
        for region in display_regions
    ]
    ax.legend(
        handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.015), ncol=max(1, len(handles)),
        frameon=True, facecolor="black", edgecolor="none", framealpha=0.62,
        fontsize=7.0, handlelength=0.9, labelcolor="white",
    )
    ax.set_axis_off()
    ax.set_anchor("C")
    ax._source_bounds_px = tuple(int(value) for value in bounds)
    ax._display_shape_px = tuple(int(value) for value in accepted_view.shape)
    ax._accepted_mask_sha256 = str(provenance["mask_sha256"])
    ax._source_dapi_sha256 = str(provenance["source_dapi_sha256"])
    ax._region_legend_labels = display_regions
    return ax


def validate_aligned_representative_pair(native_ax: plt.Axes, accepted_ax: plt.Axes) -> None:
    if getattr(native_ax, "_source_bounds_px", None) != getattr(accepted_ax, "_source_bounds_px", None):
        raise RuntimeError("Native microscopy and accepted ARC/ME panels do not use identical source coordinates")
    if getattr(native_ax, "_display_shape_px", None) != getattr(accepted_ax, "_display_shape_px", None):
        raise RuntimeError("Native microscopy and accepted ARC/ME panels do not have identical crop dimensions")
    native_dapi_sha = str(getattr(native_ax, "_source_dapi_sha256", ""))
    accepted_dapi_sha = str(getattr(accepted_ax, "_source_dapi_sha256", ""))
    if not native_dapi_sha or native_dapi_sha != accepted_dapi_sha:
        raise RuntimeError(
            "Native microscopy and accepted ROI panel do not resolve to the same source DAPI SHA-256"
        )
    labels = tuple(getattr(accepted_ax, "_region_legend_labels", ()))
    if not labels or not set(labels).issubset({"ARC", "ME", "VMN"}):
        raise RuntimeError(
            "Accepted ARC/ME/VMN panel has no valid delineated-region key")


def section_state_composition(cells: pd.DataFrame, animal: str, section: str,
                              regions: Sequence[str]) -> dict[str, dict[str, int]]:
    subset = cells.loc[cells["animal_id"].astype(str).eq(animal) & cells["section"].astype(str).eq(section)].copy()
    if "cell_qc_pass" in subset:
        accepted = subset["cell_qc_pass"].astype(str).str.lower().isin(("true", "1", "yes"))
        subset = subset.loc[accepted]
    state_column = "microglia_state" if "microglia_state" in subset else "predicted_state"
    return {
        region: {state: int(((subset["region"] == region) & (subset[state_column] == state)).sum()) for state in STATE_COLORS}
        for region in regions
    }


def _draw_section_state_qc_legacy(fig: plt.Figure, spec, counts: dict[str, dict[str, int]]) -> None:
    """Classic 100% composition plot for the representative section."""
    ax = fig.add_subplot(spec)
    states = list(STATE_COLORS)
    regions = tuple(counts)
    if not regions:
        raise RuntimeError("Representative section has no delineated ARC/ME/VMN ROI")
    y_positions = dict(zip(regions, np.linspace(1.25, 0.48, len(regions))))
    totals = {
        region: int(sum(int(counts.get(region, {}).get(state, 0)) for state in states))
        for region in regions
    }
    left = {region: 0.0 for region in regions}
    for state in states:
        for region in regions:
            total = totals[region]
            count = int(counts.get(region, {}).get(state, 0))
            percent = 100.0 * count / total if total else 0.0
            if percent > 0:
                ax.barh(
                    y_positions[region], percent, left=left[region], height=0.31,
                    color=STATE_COLORS[state], edgecolor="black", linewidth=0.42,
                    align="center", zorder=2,
                )
                if percent >= 13.0:
                    text_color = "#173844" if state == "Activated" else "white"
                    ax.text(
                        left[region] + percent / 2.0, y_positions[region], f"{percent:.0f}%",
                        ha="center", va="center", fontsize=5.2, fontweight="bold",
                        color=text_color, clip_on=True, zorder=3,
                    )
            left[region] += percent
        
    for region in regions:
        ax.text(102.5, y_positions[region], f"n={totals[region]}", ha="left", va="center", fontsize=5.8, color="black")

    ax.set_xlim(0, 116)
    ax.set_ylim(0.0, 1.58)
    ax.set_yticks([y_positions[region] for region in regions], list(regions))
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("Accepted cells (%)", fontsize=6.2, labelpad=2.0)
    ax.tick_params(axis="x", labelsize=5.8, width=0.7, length=2.5, pad=1.5)
    ax.tick_params(axis="y", labelsize=6.3, width=0.0, length=0, pad=3.0)
    for tick, region in zip(ax.get_yticklabels(), regions):
        tick.set_color(REGION_COLORS[region])
        tick.set_fontweight("bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("black")
    ax.spines["bottom"].set_color("black")
    ax.spines["left"].set_linewidth(0.75)
    ax.spines["bottom"].set_linewidth(0.75)
    ax.set_title("Microglial state composition", fontsize=8.2, fontweight="bold", pad=15.0, color="#173844")
    ax.text(0.5, 1.035, "Representative section", transform=ax.transAxes, ha="center", va="bottom", fontsize=5.8, color="#536970")
    legend_labels = {"Ramified": "Ramified", "Rod-like": "Rod-like", "Activated": "Activated", "Amoeboid": "Amoeboid"}
    handles = [Patch(facecolor=STATE_COLORS[state], edgecolor="black", linewidth=0.35, label=legend_labels[state]) for state in states]
    ax.legend(
        handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.855),
        ncol=4, frameon=False, fontsize=4.9, handlelength=0.9,
        handletextpad=0.35, columnspacing=0.8, borderaxespad=0.0,
    )

    # companion-heading-overlap repair: reserve internal headroom below anatomy.
    ax.set_ylim(0.20, 1.80)
    for _artist in [ax.title, ax._left_title, ax._right_title, *ax.texts]:
        if _artist.get_text() in {
            "Microglial state composition",
            "Representative section",
        }:
            _artist.set_visible(False)
    ax.text(
        0.0,
        0.965,
        "Microglial state composition",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9.0,
        fontweight="bold",
        color="#111111",
        zorder=20,
    )
    ax.text(
        0.0,
        0.805,
        "Representative section",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.1,
        color="#555555",
        zorder=20,
    )


def _apply_ggplot_theme(ax: plt.Axes, *, grid_axis: str = "both") -> None:
    """Apply a restrained ggplot2-like publication theme."""
    ax.set_facecolor("#EBEBEB")
    ax.set_axisbelow(True)
    ax.grid(True, axis=grid_axis, color="white", linewidth=1.15)
    ax.tick_params(colors="#333333", width=0.8, length=3.0)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_section_state_qc(
    fig: plt.Figure,
    spec,
    counts: dict[str, dict[str, int]],
    *,
    header_spec=None,
    panel_letter: str | None = None,
) -> None:
    """Compact unstacked state-by-region bars for the representative section."""
    ax = fig.add_subplot(spec)
    regions = tuple(region for region in REGIONS if region in counts)
    if set(regions) != set(REGIONS):
        raise RuntimeError(
            "Panel D representative must contain ARC, ME, and VMN composition"
        )
    states = tuple(STATE_COLORS)
    totals = {
        region: int(sum(int(counts[region].get(state, 0)) for state in states))
        for region in regions
    }
    x = np.arange(len(regions), dtype=float)
    width = 0.18
    offsets = (np.arange(len(states), dtype=float) - 1.5) * 0.20
    maximum = 0.0
    for index, state in enumerate(states):
        values = np.asarray([
            100.0 * int(counts[region].get(state, 0)) / totals[region]
            if totals[region] else 0.0
            for region in regions
        ])
        maximum = max(maximum, float(np.max(values)))
        ax.bar(
            x + offsets[index], values, width=width,
            color=STATE_COLORS[state], edgecolor="white", linewidth=0.7,
            label=state, zorder=3,
        )
    upper = max(20.0, 10.0 * math.ceil((maximum * 1.25) / 10.0))
    ax.set_ylim(0, upper)
    ax.set_xlim(-0.55, len(regions) - 0.45)
    ax.set_yticks(np.arange(0.0, upper + 0.1, 20.0))
    ax.set_ylabel("CNN-predicted cells (%)", fontsize=7.5, fontweight="bold", labelpad=2.0)
    ax.set_xticks(x, [f"{region}\nn={totals[region]}" for region in regions])
    for tick, region in zip(ax.get_xticklabels(), regions):
        tick.set_color(REGION_COLORS[region])
        tick.set_fontweight("bold")
        tick.set_fontsize(7.2)
    ax.tick_params(axis="y", labelsize=6.8)
    _apply_ggplot_theme(ax, grid_axis="y")
    handles, labels = ax.get_legend_handles_labels()
    display_labels = ["Rod-like" if label == "Rod-like" else label for label in labels]
    if header_spec is None:
        if panel_letter:
            panel_label(ax, panel_letter)
        ax.text(
            0.012, 0.975, "Representative-section microglial states",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=8.6, fontweight="bold", color="#111111", zorder=5,
        )
        ax.legend(
            handles, display_labels,
            loc="upper right", bbox_to_anchor=(0.995, 0.995), ncol=4,
            frameon=False, fontsize=5.8, handlelength=0.9,
            handletextpad=0.25, columnspacing=0.55, borderaxespad=0.0,
        )
    else:
        header_ax = fig.add_subplot(header_spec)
        header_ax.set_axis_off()
        if panel_letter:
            header_ax.text(
                -0.018, 1.03, panel_letter, transform=header_ax.transAxes,
                fontsize=22, fontweight="bold", ha="right", va="bottom",
                color="black", clip_on=False,
            )
        header_ax.text(
            0.01, 0.58, "Representative-section microglial states",
            transform=header_ax.transAxes, ha="left", va="center",
            fontsize=8.6, fontweight="bold", color="#111111",
        )
        header_ax.legend(
            handles, display_labels,
            loc="center right", bbox_to_anchor=(0.995, 0.58), ncol=4,
            frameon=False, fontsize=5.8, handlelength=0.9,
            handletextpad=0.25, columnspacing=0.55, borderaxespad=0.0,
        )


def gfap_network_metrics(positive: np.ndarray, labels: np.ndarray, um_per_px: float, *, me_absent: bool = False) -> dict[str, dict[str, float]]:
    from skimage.morphology import skeletonize
    positive = np.asarray(positive, dtype=bool)
    labels = np.asarray(labels)
    codes = _arc_me_codes(labels, allow_me_absent=me_absent)
    pixel_area_um2 = float(um_per_px) ** 2
    metrics: dict[str, dict[str, float]] = {}
    kernel = np.ones((3, 3), dtype=np.uint8)
    kernel[1, 1] = 0
    drawn_regions = tuple(
        region for region in ("ARC", "ME", "VMN") if np.any(labels == codes[region])
    )
    for region in drawn_regions:
        region_mask = labels == codes[region]
        regional_skeleton = skeletonize(positive & region_mask)
        neighbors = convolve(regional_skeleton.astype(np.uint8), kernel, mode="constant", cval=0)
        branchpoints = int(np.count_nonzero(regional_skeleton & (neighbors >= 3)))
        endpoints = int(np.count_nonzero(regional_skeleton & (neighbors == 1)))
        horizontal = int(np.count_nonzero(regional_skeleton[:, :-1] & regional_skeleton[:, 1:]))
        vertical = int(np.count_nonzero(regional_skeleton[:-1, :] & regional_skeleton[1:, :]))
        diagonal = int(np.count_nonzero(regional_skeleton[:-1, :-1] & regional_skeleton[1:, 1:])) + int(np.count_nonzero(regional_skeleton[1:, :-1] & regional_skeleton[:-1, 1:]))
        skeleton_length_um = float((horizontal + vertical + math.sqrt(2.0) * diagonal) * float(um_per_px))
        region_area_um2 = float(np.count_nonzero(region_mask) * pixel_area_um2)
        scale = 1000.0 / region_area_um2 if region_area_um2 > 0 else float("nan")
        metrics[region] = {
            "region_area_um2": region_area_um2,
            "skeleton_length_um": skeleton_length_um,
            "branchpoints": branchpoints,
            "endpoints": endpoints,
            "skeleton_length_um_per_1000um2": skeleton_length_um * scale,
            "branchpoints_per_1000um2": branchpoints * scale,
            "endpoints_per_1000um2": endpoints * scale,
        }
    return metrics


def draw_gfap_network_qc(fig: plt.Figure, spec, metrics: dict[str, dict[str, float]]) -> None:
    """Three-region Cleveland dot plot in a restrained ggplot-like style."""
    ax = fig.add_subplot(spec)
    rows = (
        ("Skeleton length", "skeleton_length_um_per_1000um2"),
        ("Branch points", "branchpoints_per_1000um2"),
        ("Endpoints", "endpoints_per_1000um2"),
    )
    display_regions = tuple(region for region in REGIONS if region in metrics)
    if set(display_regions) != set(REGIONS):
        raise RuntimeError("Panel B representative must contain ARC, ME, and VMN")
    y = np.arange(len(rows), dtype=float)[::-1]
    values_by_region = {
        region: np.asarray(
            [float(metrics[region][metric]) for _, metric in rows], dtype=float
        )
        for region in display_regions
    }
    finite_positive = np.concatenate(list(values_by_region.values()))
    finite_positive = finite_positive[np.isfinite(finite_positive) & (finite_positive > 0)]
    if not finite_positive.size:
        raise RuntimeError("GFAP network morphology requires positive finite densities for logarithmic display")

    marker_by_region = {"ARC": "o", "ME": "^", "VMN": "s"}
    offsets = {"ARC": 0.16, "ME": 0.0, "VMN": -0.16}
    for region, values in values_by_region.items():
        valid = np.isfinite(values) & (values > 0)
        ax.scatter(
            values[valid], y[valid] + offsets[region], s=58,
            marker=marker_by_region[region], facecolor=REGION_COLORS[region],
            edgecolor="white", linewidth=0.8, label=region, zorder=3,
        )

    lower = 10.0 ** math.floor(math.log10(float(np.min(finite_positive)) * 0.82))
    upper = 10.0 ** math.ceil(math.log10(float(np.max(finite_positive)) * 1.18))
    if upper <= lower:
        upper = lower * 10.0
    ax.set_xscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(-0.45, 2.45)
    ax.set_yticks(y, [label_text for label_text, _ in rows])
    ax.set_xlabel("Density per 1,000 µm² (log scale)", fontsize=8.4, fontweight="bold")
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(axis="both", labelsize=7.6)
    ax.set_title(
        "GFAP network morphology · representative section",
        loc="left", fontsize=10.0, fontweight="bold", pad=9.0,
    )
    _apply_ggplot_theme(ax, grid_axis="x")
    ax.legend(
        loc="lower right", bbox_to_anchor=(1.0, 1.015), ncol=3,
        frameon=False, fontsize=7.0, handletextpad=0.35,
        columnspacing=1.0, borderaxespad=0.0,
    )
    ax.set_xlabel(
        r"Density per 1,000 $\mu m^2$ (log scale)",
        fontsize=7.1,
        labelpad=2.5,
    )


def draw_cnn(
    ax: plt.Axes,
    input_patch: np.ndarray,
    metadata: dict,
    history: pd.DataFrame,
    *,
    um_per_px: float,
    inside_label: bool = False,
    panel_letter: str = "C",
) -> None:
    ax.set_axis_off()
    if inside_label:
        ax.text(
            0.01, 1.03, panel_letter, transform=ax.transAxes,
            fontsize=22, fontweight="bold", ha="left", va="bottom",
            color="black", clip_on=False,
        )
    else:
        panel_label(ax, panel_letter)
    ax.text(0.005, 0.98, "Image-patch CNN", transform=ax.transAxes, fontsize=16.0, fontweight="bold", va="top", color="#173844")
    patch_ax = ax.inset_axes([0.015, 0.37, 0.105, 0.45])
    patch_ax.imshow(input_patch)
    add_scalebar(
        patch_ax, input_patch.shape[1], input_patch.shape[0], 25.0, um_per_px, inset=True,
    )
    patch_ax.set_xticks([])
    patch_ax.set_yticks([])
    for spine in patch_ax.spines.values():
        spine.set_color("#173844")
        spine.set_linewidth(1.3)
    ax.text(0.068, 0.29, "96 x 96 x 3", transform=ax.transAxes, ha="center", fontsize=10.0, fontweight="bold")
    ax.text(0.068, 0.20, "Iba1 | mask | skeleton", transform=ax.transAxes, ha="center", fontsize=8.8)
    arrow(ax, (0.125, 0.59), (0.18, 0.59))
    stacks = ((0.19, 0.40, "#d9e8ed", "16", "48 x 48"), (0.31, 0.34, "#bcdde0", "32", "24 x 24"), (0.43, 0.28, "#95ceca", "64", "12 x 12"), (0.55, 0.24, "#6dbbb4", "96", "12 x 12"))
    for index, (x, height, color, channels, spatial) in enumerate(stacks):
        draw_feature_stack(ax, x, height, color, channels, spatial)
        if index < len(stacks) - 1:
            arrow(ax, (x + 0.085, 0.59), (stacks[index + 1][0] - 0.01, 0.59))
    arrow(ax, (0.635, 0.59), (0.69, 0.59))
    ax.add_patch(Circle((0.72, 0.59), 0.036, transform=ax.transAxes, facecolor="#f2cf78", edgecolor="#173844", linewidth=1.3))
    ax.text(0.72, 0.59, "GAP", transform=ax.transAxes, ha="center", va="center", fontsize=8.2, fontweight="bold")
    for index, state in enumerate(STATE_COLORS):
        y = 0.79 - index * 0.135
        arrow(ax, (0.758, 0.59), (0.81, y), color="#536970", width=0.9)
        ax.add_patch(Circle((0.835, y), 0.021, transform=ax.transAxes, facecolor=STATE_COLORS[state], edgecolor="#173844", linewidth=1.0))
        ax.text(0.866, y, state, transform=ax.transAxes, ha="left", va="center", fontsize=9.5, fontweight="bold")
    candidates = int(metadata.get("training_candidates_n", 10848) or 10848)
    if not history.empty and "val_loss" in history:
        best = history.loc[pd.to_numeric(history["val_loss"], errors="coerce").idxmin()]
        best_epoch = int(best["epoch"])
        balanced = float(best.get("val_balanced_accuracy", np.nan))
    else:
        best_epoch, balanced = 25, np.nan
    metric = f" | pseudo-label agreement = {balanced:.2f}" if math.isfinite(balanced) else ""
    ax.text(0.18, 0.105, "Conv 3x3: 16 -> 32 -> 64 -> 96 | GAP | four-class softmax", transform=ax.transAxes, fontsize=9.5, fontweight="bold", color="#38545d")
    ax.text(0.18, 0.025, f"{candidates:,} pseudo-labels | animal-blocked 80/20 | best epoch {best_epoch}{metric}", transform=ax.transAxes, fontsize=9.0, fontweight="bold", color="#173844")


def resolve_patch(row: pd.Series, analysis_dir: Path | None = None) -> Path:
    path = Path(str(row["patch_path"]))
    if path.is_file():
        return path
    root = (
        Path(analysis_dir).expanduser().resolve()
        if analysis_dir is not None
        else FIG5_ANALYSIS_ROOT
    )
    candidates = [
        root / "microglia_patches" / str(row["animal_id"])
        / str(row["section"]) / f"{row['cell_uid']}.png",
    ]
    if "microglia_patches" in path.parts:
        offset = path.parts.index("microglia_patches")
        candidates.append(root.joinpath(*path.parts[offset:]))
    candidates.append(
        FIG5_ANALYSIS_ROOT / "archive" / "patches"
        / "microglia_patch_snapshots" / "final_hybrid_run"
        / str(row["animal_id"]) / str(row["section"])
        / f"{row['cell_uid']}.png"
    )
    return next((candidate for candidate in candidates if candidate.is_file()), candidates[0])


def _cell_example_qc_reasons(row: pd.Series, state: str) -> list[str]:
    reasons: list[str] = []
    maximum_area = {"Ramified": 500.0, "Rod-like": 800.0, "Activated": 500.0, "Amoeboid": 600.0}

    def value(name: str) -> float:
        return float(pd.to_numeric(pd.Series([row.get(name, np.nan)]), errors="coerce").iloc[0])

    source = str(row.get("classifier_source", "")).lower()
    uncertain = str(row.get("classifier_uncertain", False)).lower() in ("true", "1", "yes")
    confidence, area = value("classifier_confidence"), value("area_px")
    circularity, solidity = value("circularity"), value("solidity")
    eccentricity, extent = value("eccentricity"), value("extent")
    if "cnn" not in source:
        reasons.append("not a CNN prediction")
    if uncertain:
        reasons.append("classifier uncertain")
    if not math.isfinite(confidence) or confidence < 0.75:
        reasons.append("confidence <0.75")
    if not math.isfinite(area) or not 20.0 <= area <= maximum_area[state]:
        reasons.append(f"area outside 20-{maximum_area[state]:g} px")
    if not math.isfinite(circularity) or not 0.0 < circularity <= 1.05:
        reasons.append("invalid circularity")
    if not math.isfinite(solidity) or not 0.0 < solidity <= 1.01:
        reasons.append("invalid solidity")
    if not math.isfinite(eccentricity) or not 0.0 <= eccentricity < 1.0:
        reasons.append("invalid eccentricity")
    if not math.isfinite(extent) or not 0.05 <= extent <= 1.0:
        reasons.append("invalid extent")
    if state == "Amoeboid":
        if math.isfinite(circularity) and circularity < 0.35:
            reasons.append("amoeboid circularity <0.35")
        if math.isfinite(solidity) and solidity < 0.70:
            reasons.append("amoeboid solidity <0.70")
        if math.isfinite(eccentricity) and eccentricity > 0.85:
            reasons.append("amoeboid eccentricity >0.85")
        if math.isfinite(extent) and extent < 0.35:
            reasons.append("amoeboid extent <0.35")
    return reasons


def _patch_example_qc(
    row: pd.Series, state: str, analysis_dir: Path,
) -> tuple[bool, str]:
    path = resolve_patch(row, analysis_dir)
    if not path.is_file():
        return False, "patch file missing"
    try:
        patch = np.asarray(mpimg.imread(path))
    except Exception:
        return False, "patch unreadable"
    if patch.ndim != 3 or patch.shape[2] < 2 or min(patch.shape[:2]) < 48:
        return False, "invalid patch dimensions"
    green = np.asarray(patch[..., 1], dtype=float)
    if float(np.nanmax(green)) > 1.5:
        green = green / 255.0
    mask = np.isfinite(green) & (green >= 0.50)
    _component_labels, component_count = ndi_label(mask)
    if int(component_count) != 1:
        return False, f"selected-cell mask has {int(component_count)} components"
    mask_pixels = int(mask.sum())
    expected_area = float(pd.to_numeric(pd.Series([row.get("area_px", np.nan)]), errors="coerce").iloc[0])
    if math.isfinite(expected_area) and abs(mask_pixels - expected_area) > max(3.0, 0.15 * expected_area):
        return False, "patch mask disagrees with measured area"
    clearance = 3 if state == "Amoeboid" else 2
    if mask[:clearance, :].any() or mask[-clearance:, :].any() or mask[:, :clearance].any() or mask[:, -clearance:].any():
        return False, f"cell mask enters {clearance}-px patch border"
    yy, xx = np.nonzero(mask)
    center_y, center_x = 0.5 * (patch.shape[0] - 1), 0.5 * (patch.shape[1] - 1)
    if math.hypot(float(xx.mean()) - center_x, float(yy.mean()) - center_y) > 12.0:
        return False, "cell mask is not centered"
    return True, "passed cell and patch QC"


def select_examples(
    cells: pd.DataFrame, analysis_dir: Path,
) -> dict[str, dict[str, object]]:
    features = ["area_px", "circularity", "solidity", "soma_fraction", "skeleton_length_px", "skeleton_branchpoints"]
    selections: dict[str, dict[str, object]] = {}
    for state in STATE_COLORS:
        state_rows = cells.loc[cells["microglia_state"].eq(state)].copy()
        rejection_counts: dict[str, int] = {}
        passed_indices: list[object] = []
        for index, row in state_rows.iterrows():
            reasons = _cell_example_qc_reasons(row, state)
            if reasons:
                for reason in reasons:
                    rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
            else:
                passed_indices.append(index)
        subset = state_rows.loc[passed_indices].copy()
        selected: pd.Series | None = None
        patches_checked = 0
        if not subset.empty:
            numeric = subset[features].apply(pd.to_numeric, errors="coerce")
            med = numeric.median()
            mad = numeric.sub(med).abs().median().replace(0, 1).fillna(1)
            subset["_representative_score"] = numeric.sub(med).div(mad).abs().sum(axis=1, min_count=1).fillna(np.inf) - 0.4 * pd.to_numeric(subset["classifier_confidence"], errors="coerce").fillna(0)
            subset = subset.sort_values(["_representative_score", "classifier_confidence", "cell_uid"], ascending=[True, False, True])
            for _, row in subset.iterrows():
                patches_checked += 1
                passed, reason = _patch_example_qc(row, state, analysis_dir)
                if passed:
                    selected = row
                    break
                rejection_counts[reason] = rejection_counts.get(reason, 0) + 1
        if selected is None:
            reason = "; ".join(f"{count} {label}" for label, count in sorted(rejection_counts.items(), key=lambda item: (-item[1], item[0])))
            reason = reason or "no predictions available"
        else:
            reason = f"selected after checking {patches_checked} patch{'es' if patches_checked != 1 else ''}"
        selections[state] = {
            "row": selected,
            "total_predictions": int(len(state_rows)),
            "cell_qc_candidates": int(len(subset)),
            "patches_checked": int(patches_checked),
            "reason": reason,
        }
    return selections


def draw_examples(
    ax: plt.Axes, selections: dict[str, dict[str, object]], um_per_px: float,
    analysis_dir: Path,
) -> None:
    ax.set_axis_off()
    ax.text(0.015, 1.08, "Model-predicted cell examples", transform=ax.transAxes, fontsize=12.0, fontweight="bold", va="bottom", color="#173844", clip_on=False)
    for index, state in enumerate(STATE_COLORS):
        entry = selections[state]
        row = entry["row"]
        left = 0.015 + index * 0.245
        image_ax = ax.inset_axes([left, 0.08, 0.225, 0.86])
        if row is None:
            image_ax.set_facecolor("#f1f2ef")
            image_ax.text(0.5, 0.55, f"No QC-valid\n{state} prediction", transform=image_ax.transAxes, ha="center", va="center", fontsize=9.0, fontweight="bold", color="#59656a")
            image_ax.text(0.5, 0.16, f"0 / {entry['total_predictions']}", transform=image_ax.transAxes, ha="center", va="center", fontsize=8.2, color="#768287")
        else:
            patch = mpimg.imread(resolve_patch(row, analysis_dir))
            image_ax.imshow(patch)
            add_scalebar(
                image_ax, patch.shape[1], patch.shape[0], 25.0, um_per_px, inset=True,
            )
        image_ax.set_xticks([])
        image_ax.set_yticks([])
        for spine in image_ax.spines.values():
            spine.set_color(STATE_COLORS[state] if row is not None else "#9aa3a6")
            spine.set_linewidth(2.0)
            spine.set_linestyle("-")
        ax.text(left + 0.1125, -0.005, state, transform=ax.transAxes, ha="center", va="top", fontsize=11.5, fontweight="bold", color=STATE_COLORS[state] if row is not None else "#59656a")


def draw_gfap_workflow(
    ax: plt.Axes,
    dapi: np.ndarray,
    gfap: np.ndarray,
    positive: np.ndarray,
    labels: np.ndarray,
    center: tuple[float, float],
    *,
    um_per_px: float,
    panel_letter: str | None = None,
) -> None:
    ax.set_axis_off()
    if panel_letter:
        ax.text(
            -0.018, 0.965, panel_letter, transform=ax.transAxes,
            fontsize=22, fontweight="bold", ha="right", va="top",
            color="black", clip_on=False,
        )
    ax.text(0.01, 0.965, "GFAP quantification", transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="top", color="#173844")
    size = 118
    raw_crop, _ = crop_square(composite(dapi, gfap, (0.08, 1.0, 0.32)), center, size)
    positive_crop, _ = crop_square(positive.astype(float), center, size)
    images = (raw_crop, positive_crop)
    titles = ("Raw GFAP + DAPI", "GFAP-positive mask")
    for index, (image_data, title) in enumerate(zip(images, titles)):
        left = 0.05 + index * 0.49
        iax = ax.inset_axes([left, 0.47, 0.40, 0.37])
        iax.imshow(image_data, cmap="gray" if index == 1 else None, interpolation="nearest")
        add_scalebar(
            iax, image_data.shape[1], image_data.shape[0], 25.0, um_per_px, inset=True,
        )
        iax.set_xticks([])
        iax.set_yticks([])
        for spine in iax.spines.values():
            spine.set_color("#173844")
            spine.set_linewidth(1.3)
        ax.text(
            left + 0.20, 0.405, title, transform=ax.transAxes,
            ha="center", va="top", fontsize=9.6, fontweight="bold",
        )
        if index == 0:
            arrow(ax, (0.455, 0.64), (0.535, 0.64))
            ax.text(0.495, 0.73, "background correction\n+ threshold", transform=ax.transAxes, ha="center", va="center", fontsize=7.2, color="#536970")
    outputs = (
        (0.07, "Corrected signal / ROI", "raw intensity"),
        (0.53, "GFAP+ pixels / ROI", "positive fraction"),
    )
    for x, heading, detail in outputs:
        ax.add_patch(Rectangle((x, 0.10), 0.40, 0.20, transform=ax.transAxes, facecolor="#EBEBEB", edgecolor="white", linewidth=1.5))
        ax.text(x + 0.20, 0.215, heading, transform=ax.transAxes, ha="center", va="center", fontsize=9.4, fontweight="bold", color="#173844")
        ax.text(x + 0.20, 0.145, detail, transform=ax.transAxes, ha="center", va="center", fontsize=7.8, color="#536970")


def format_p(value: float) -> str:
    if not math.isfinite(value):
        return "p = NA"
    if value < 0.0001:
        return "p < 0.0001"
    if value < 0.001:
        return f"p = {value:.4f}"
    return f"p = {value:.3f}"


def format_p_token(value: float) -> str:
    """Match the compact p-value token used inside Figure 4D."""
    if not math.isfinite(value):
        return "NA"
    return f"{value:.1e}" if value < 0.001 else f"{value:.3g}"


def tukey_outliers(values: np.ndarray) -> np.ndarray:
    """Use Figure 4D's finite-value Tukey 1.5-IQR display rule exactly."""
    array = np.asarray(values, dtype=float)
    outliers = np.zeros(array.shape, dtype=bool)
    finite = np.isfinite(array)
    if int(finite.sum()) < 4:
        return outliers
    q1, q3 = np.percentile(array[finite], (25.0, 75.0))
    iqr = q3 - q1
    if math.isfinite(float(iqr)) and iqr > 0:
        outliers[finite] = (
            (array[finite] < q1 - 1.5 * iqr)
            | (array[finite] > q3 + 1.5 * iqr)
        )
    return outliers


def cliffs_delta(second: np.ndarray, first: np.ndarray) -> float:
    left = np.asarray(second, dtype=float); right = np.asarray(first, dtype=float)
    left = left[np.isfinite(left)]; right = right[np.isfinite(right)]
    if not left.size or not right.size:
        return float("nan")
    differences = left[:, None] - right[None, :]
    return float((np.count_nonzero(differences > 0) - np.count_nonzero(differences < 0)) / differences.size)


def _multinomial_partition_count(sizes: Sequence[int]) -> int:
    total = 1
    remaining = int(sum(sizes))
    for size in list(sizes)[:-1]:
        total *= math.comb(remaining, int(size))
        remaining -= int(size)
    return int(total)


def one_way_anova_receipt(
    groups: list[np.ndarray],
) -> dict[str, float | int | bool]:
    """Ordinary equal-variance one-way ANOVA result and degrees of freedom."""
    usable = [np.asarray(group, dtype=float) for group in groups]
    usable = [group[np.isfinite(group)] for group in usable if np.count_nonzero(np.isfinite(group)) >= 2]
    if len(usable) < 2:
        return {
            "statistic": float("nan"), "p_value": float("nan"),
            "df_between": float("nan"), "df_within": float("nan"),
            "extreme_labelings": float("nan"), "enumerated_labelings": float("nan"),
            "permutation_exhaustive": float("nan"),
        }
    try:
        result = f_oneway(*usable)
    except Exception:
        return {
            "statistic": float("nan"), "p_value": float("nan"),
            "df_between": float("nan"), "df_within": float("nan"),
            "extreme_labelings": float("nan"), "enumerated_labelings": float("nan"),
            "permutation_exhaustive": float("nan"),
        }
    return {
        "statistic": float(result.statistic), "p_value": float(result.pvalue),
        "df_between": int(len(usable) - 1),
        "df_within": int(sum(group.size for group in usable) - len(usable)),
        "extreme_labelings": float("nan"),
        "enumerated_labelings": float("nan"),
        "permutation_exhaustive": float("nan"),
    }


def exact_permutation_mwu_receipt(
    first: np.ndarray, second: np.ndarray, max_partitions: int = 1_000_000,
) -> dict[str, float | int | bool]:
    """Exact two-sided rank-sum permutation result and enumeration receipt."""
    left = np.asarray(first, dtype=float); right = np.asarray(second, dtype=float)
    left = left[np.isfinite(left)]; right = right[np.isfinite(right)]
    if left.size < 2 or right.size < 2:
        return {
            "statistic": float("nan"), "p_value": float("nan"),
            "extreme_labelings": 0, "enumerated_labelings": 0,
            "permutation_exhaustive": False,
        }
    pooled = np.concatenate((left, right)); n_first = int(left.size)
    total_partitions = math.comb(int(pooled.size), n_first)
    if total_partitions > max_partitions:
        raise RuntimeError(
            f"Exact rank-sum enumeration requires {total_partitions} labelings; "
            f"cap is {max_partitions}"
        )
    ranks = rankdata(pooled, method="average")
    expected_rank_sum = n_first * (pooled.size + 1.0) / 2.0
    observed = abs(float(np.sum(ranks[:n_first])) - expected_rank_sum)
    extreme = sum(
        abs(float(np.sum(ranks[list(selection)])) - expected_rank_sum) >= observed - 1e-12
        for selection in combinations(range(int(pooled.size)), n_first)
    )
    observed_u = float(mannwhitneyu(left, right, alternative="two-sided").statistic)
    return {
        "statistic": observed_u, "p_value": float(extreme / total_partitions),
        "extreme_labelings": int(extreme),
        "enumerated_labelings": int(total_partitions),
        "permutation_exhaustive": True,
    }


def exact_permutation_mwu(first: np.ndarray, second: np.ndarray, max_partitions: int = 1_000_000) -> float:
    """Exact two-sided independent-label permutation p for rank-sum/U."""
    return float(exact_permutation_mwu_receipt(first, second, max_partitions)["p_value"])


def endpoint_statistics(data: pd.DataFrame, endpoint: str, region: str) -> tuple[dict, list[dict]]:
    groups: dict[str, np.ndarray] = {}
    rows: list[dict] = []
    for condition in CONDITIONS:
        values = pd.to_numeric(data.loc[data["region"].eq(region) & data["cond"].eq(condition), endpoint], errors="coerce").dropna().to_numpy(float)
        groups[condition] = values
        rows.append({
            "endpoint": endpoint, "region": region, "test": "descriptive",
            "comparison": condition, "n": int(values.size),
            "mean": float(np.mean(values)) if values.size else np.nan,
            "sd": float(np.std(values, ddof=1)) if values.size > 1 else np.nan,
            "statistic": np.nan, "p_value_raw": np.nan,
            "effect_size": np.nan, "experimental_unit": "animal",
            "extreme_labelings": np.nan, "enumerated_labelings": np.nan,
            "permutation_exhaustive": np.nan,
        })
    usable = [groups[name] for name in CONDITIONS if groups[name].size]
    global_result = one_way_anova_receipt(usable)
    global_p = float(global_result["p_value"])
    rows.append({
        "endpoint": endpoint, "region": region, "test": "One_way_ANOVA",
        "p_value_method": "ordinary equal-variance one-way ANOVA F",
        "comparison": "Water|Sucrose|Allulose", "n": int(sum(v.size for v in usable)),
        "mean": np.nan, "sd": np.nan,
        "statistic": float(global_result["statistic"]), "p_value_raw": global_p,
        "effect_size": np.nan,
        "experimental_unit": "animal",
        "df_between": global_result["df_between"],
        "df_within": global_result["df_within"],
        "extreme_labelings": global_result["extreme_labelings"],
        "enumerated_labelings": global_result["enumerated_labelings"],
        "permutation_exhaustive": global_result["permutation_exhaustive"],
    })
    pairs_raw: dict[tuple[str, str], float] = {}
    pair_receipts: dict[tuple[str, str], dict[str, float | int | bool]] = {}
    pair_order = (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose"))
    for left, right in pair_order:
        receipt = exact_permutation_mwu_receipt(groups[left], groups[right])
        pair_receipts[(left, right)] = receipt
        pairs_raw[(left, right)] = float(receipt["p_value"])
    for left, right in pair_order:
        receipt = pair_receipts[(left, right)]
        rows.append({
            "endpoint": endpoint, "region": region, "test": "Mann_Whitney_U_two_sided",
            "p_value_method": "exhaustive independent-label permutation of rank-sum/U",
            "comparison": f"{left}|{right}", "n": int(groups[left].size + groups[right].size),
            "n_left": int(groups[left].size), "n_right": int(groups[right].size),
            "mean": np.nan, "sd": np.nan,
            "statistic": float(receipt["statistic"]),
            "p_value_raw": pairs_raw[(left, right)],
            "effect_size": cliffs_delta(groups[right], groups[left]),
            "effect_size_definition": f"Cliff's delta ({right} minus {left})",
            "experimental_unit": "animal",
            "extreme_labelings": int(receipt["extreme_labelings"]),
            "enumerated_labelings": int(receipt["enumerated_labelings"]),
            "permutation_exhaustive": bool(receipt["permutation_exhaustive"]),
        })
    return {"global_p": global_p, "pairs_raw": pairs_raw, "groups": groups}, rows


def finalize_animal_statistics(stats: pd.DataFrame) -> pd.DataFrame:
    """Validate the animal-level ANOVA/exact-MWU result family."""
    result = stats.copy()
    global_mask = result["test"].eq("One_way_ANOVA")
    global_rows = result.loc[global_mask]
    expected = {(endpoint, region) for _, endpoint, _, _ in PANEL_ENDPOINTS for region in REGIONS}
    observed = set(zip(global_rows["endpoint"], global_rows["region"]))
    if len(global_rows) != len(expected) or observed != expected:
        raise RuntimeError(
            "Animal-level global ANOVA family must contain exactly "
            f"{len(PANEL_ENDPOINTS)} endpoints x {len(REGIONS)} regions"
        )
    anova_rows = result.loc[global_mask]
    if not np.isfinite(pd.to_numeric(anova_rows["p_value_raw"], errors="coerce")).all():
        raise RuntimeError("Every animal-level one-way ANOVA must be finite")
    inferential = result["test"].eq("Mann_Whitney_U_two_sided")
    if not result.loc[inferential, "permutation_exhaustive"].astype(bool).all():
        raise RuntimeError("Every animal-level MWU result must be exhaustive")
    if not (
        result.loc[inferential, "extreme_labelings"].astype(int)
        <= result.loc[inferential, "enumerated_labelings"].astype(int)
    ).all():
        raise RuntimeError("Invalid animal-level exact-permutation receipt")
    return result


def _validate_cage_assignments(data: pd.DataFrame) -> None:
    required = {"animal_id", "region", "cond", "cage"}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError("Cage sensitivity columns are missing: " + ", ".join(missing))
    if data["cage"].isna().any() or data["cage"].astype(str).str.strip().eq("").any():
        raise ValueError("Every Figure 5 animal requires a nonempty cage identifier")
    cage_conditions = data[["cage", "cond"]].drop_duplicates()
    if cage_conditions["cage"].duplicated().any():
        conflicts = sorted(cage_conditions.loc[cage_conditions["cage"].duplicated(False), "cage"].astype(str).unique())
        raise ValueError("Cages map to multiple conditions: " + ", ".join(conflicts))
    if set(cage_conditions["cond"].astype(str)) != set(CONDITIONS):
        raise ValueError("Cage assignments must contain Water, Sucrose, and Allulose")
    duplicate = data.duplicated(["animal_id", "region"], keep=False)
    if duplicate.any():
        raise ValueError("Expected one animal-level row per animal and region")


def activated_denominator_audit(data: pd.DataFrame) -> pd.DataFrame:
    """Audit every activated/total-microglia pair without outcome filtering."""
    required = {"activated_microglia_cells", "microglia_cells", "activated_microglia_fraction"}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError("Activated-count columns are missing: " + ", ".join(missing))
    rows: list[dict[str, object]] = []
    for record in data.sort_values(["region", "cond", "cage", "animal_id"]).itertuples(index=False):
        successes = float(record.activated_microglia_cells)
        trials = float(record.microglia_cells)
        if (
            not math.isfinite(successes) or not math.isfinite(trials)
            or successes < 0 or trials < 0 or successes > trials
            or not math.isclose(successes, round(successes), rel_tol=0, abs_tol=1e-9)
            or not math.isclose(trials, round(trials), rel_tol=0, abs_tol=1e-9)
        ):
            raise ValueError(
                f"Invalid activated/total microglia counts for {record.animal_id}/{record.region}"
            )
        reported = float(record.activated_microglia_fraction)
        ratio_estimable = trials > 0
        recomputed = successes / trials if ratio_estimable else float("nan")
        if ratio_estimable and (
            not math.isfinite(reported)
            or not math.isclose(reported, recomputed, rel_tol=1e-10, abs_tol=1e-12)
        ):
            raise ValueError(
                f"Activated fraction/count mismatch for {record.animal_id}/{record.region}"
            )
        if not ratio_estimable and math.isfinite(reported):
            raise ValueError(
                f"Zero denominator has a finite fraction for {record.animal_id}/{record.region}"
            )
        rows.append({
            "animal_id": str(record.animal_id), "region": str(record.region),
            "condition": str(record.cond), "cage": str(record.cage),
            "activated_microglia_cells": int(round(successes)),
            "total_microglia_cells": int(round(trials)),
            "activated_fraction_reported": reported,
            "activated_fraction_recomputed": recomputed,
            "counts_valid": True, "ratio_estimable": ratio_estimable,
            "included_in_cage_count_aggregation": True,
            "independent_permutation_unit": False,
            "aggregation_status": (
                "retained_valid_counts"
                if ratio_estimable else
                "retained_zero_over_zero_counts; ratio intrinsically non-estimable"
            ),
            "outcome_dependent_filtering": False,
        })
    audit = pd.DataFrame(rows)
    if len(audit) != len(data) or not audit["included_in_cage_count_aggregation"].all():
        raise RuntimeError("Activated denominator audit lost an animal-region row")
    return audit


def multipanel_sensitivity_cage_values(
    data: pd.DataFrame, denominator_audit: pd.DataFrame,
) -> pd.DataFrame:
    """Create unweighted cage means and activated-count cage aggregates."""
    _validate_cage_assignments(data)
    rows: list[dict[str, object]] = []
    endpoints = [endpoint for _, endpoint, _, _ in PANEL_ENDPOINTS]
    for (region, cage, condition), group in data.groupby(
        ["region", "cage", "cond"], sort=True, observed=True,
    ):
        n_animals = int(group["animal_id"].nunique())
        for endpoint in endpoints:
            numeric = pd.to_numeric(group[endpoint], errors="coerce")
            finite = numeric[np.isfinite(numeric)]
            if finite.empty:
                raise RuntimeError(
                    f"Cage mean is intrinsically non-estimable for {endpoint}/{region}/{cage}"
                )
            rows.append({
                "sensitivity_family": "cage_mean_anova", "endpoint": endpoint,
                "region": str(region), "cage": str(cage), "condition": str(condition),
                "value": float(finite.mean()), "numerator": np.nan,
                "denominator": np.nan, "n_animals_total": n_animals,
                "n_animals_estimable": int(finite.size),
                "n_zero_denominator_animals": (
                    int((pd.to_numeric(group["microglia_cells"], errors="coerce") == 0).sum())
                    if endpoint == "activated_microglia_fraction" else 0
                ),
                "aggregation": "unweighted mean of finite animal-level values",
                "no_outcome_dependent_filtering": True,
            })
        activated = pd.to_numeric(group["activated_microglia_cells"], errors="raise").astype(float)
        total = pd.to_numeric(group["microglia_cells"], errors="raise").astype(float)
        successes = int(round(float(activated.sum())))
        trials = int(round(float(total.sum())))
        if trials <= 0 or successes < 0 or successes > trials:
            raise RuntimeError(
                f"Cage count denominator is not estimable for {region}/{cage}"
            )
        audit_group = denominator_audit.loc[
            denominator_audit["region"].eq(str(region))
            & denominator_audit["cage"].eq(str(cage))
        ]
        if len(audit_group) != len(group):
            raise RuntimeError(f"Denominator audit/cage aggregation mismatch for {region}/{cage}")
        rows.append({
            "sensitivity_family": "cage_count_binomial_deviance",
            "endpoint": "activated_microglia_fraction", "region": str(region),
            "cage": str(cage), "condition": str(condition),
            "value": float(successes / trials), "numerator": successes,
            "denominator": trials, "n_animals_total": n_animals,
            "n_animals_estimable": int((total > 0).sum()),
            "n_zero_denominator_animals": int((total == 0).sum()),
            "aggregation": "sum activated microglia / sum total microglia within cage",
            "no_outcome_dependent_filtering": True,
        })
    cage_values = pd.DataFrame(rows).sort_values(
        ["sensitivity_family", "endpoint", "region", "condition", "cage"]
    ).reset_index(drop=True)
    expected_cages = data["cage"].nunique()
    for family, expected_endpoints in (
        ("cage_mean_anova", len(PANEL_ENDPOINTS)),
        ("cage_count_binomial_deviance", 1),
    ):
        subset = cage_values.loc[cage_values["sensitivity_family"].eq(family)]
        expected_rows = expected_cages * len(REGIONS) * expected_endpoints
        if len(subset) != expected_rows:
            raise RuntimeError(f"{family} cage-value table has {len(subset)} rows; expected {expected_rows}")
    return cage_values


def _fixed_count_labelings(counts: Sequence[int]) -> list[np.ndarray]:
    """Enumerate unique label vectors with fixed condition counts."""
    sizes = tuple(int(value) for value in counts)
    n_total = int(sum(sizes))
    indices = tuple(range(n_total))
    output: list[np.ndarray] = []

    def recurse(label_index: int, remaining: tuple[int, ...], labels: np.ndarray) -> None:
        if label_index == len(sizes) - 1:
            labels[list(remaining)] = label_index
            output.append(labels.copy())
            return
        for selected in combinations(remaining, sizes[label_index]):
            updated = labels.copy()
            updated[list(selected)] = label_index
            selected_set = set(selected)
            recurse(
                label_index + 1,
                tuple(index for index in remaining if index not in selected_set),
                updated,
            )

    recurse(0, indices, np.full(n_total, -1, dtype=np.int8))
    expected = _multinomial_partition_count(sizes)
    if len(output) != expected or any(np.any(labels < 0) for labels in output):
        raise RuntimeError(f"Fixed-count labeling enumeration failed: {len(output)} != {expected}")
    return output


def _xlog_ratio(observed: float, expected: float) -> float:
    if observed == 0:
        return 0.0
    if expected <= 0:
        return float("inf")
    return float(observed * math.log(observed / expected))


def cage_binomial_deviance(
    successes: np.ndarray, trials: np.ndarray, labels: np.ndarray,
) -> float:
    """Binomial deviance for condition-specific versus pooled cage counts."""
    successes = np.asarray(successes, dtype=float)
    trials = np.asarray(trials, dtype=float)
    labels = np.asarray(labels, dtype=int)
    if (
        successes.shape != trials.shape or successes.shape != labels.shape
        or np.any(successes < 0) or np.any(trials <= 0) or np.any(successes > trials)
        or not np.allclose(successes, np.round(successes), rtol=0, atol=1e-9)
        or not np.allclose(trials, np.round(trials), rtol=0, atol=1e-9)
    ):
        raise ValueError("Count sensitivity requires positive integer cage denominators")
    pooled_rate = float(successes.sum() / trials.sum())
    if pooled_rate <= 0 or pooled_rate >= 1:
        return 0.0
    deviance = 0.0
    for label_value in sorted(np.unique(labels)):
        selected = labels == label_value
        group_successes = float(successes[selected].sum())
        group_trials = float(trials[selected].sum())
        group_failures = group_trials - group_successes
        deviance += 2.0 * (
            _xlog_ratio(group_successes, group_trials * pooled_rate)
            + _xlog_ratio(group_failures, group_trials * (1.0 - pooled_rate))
        )
    return float(max(0.0, deviance))


def exact_cage_label_binomial_deviance(
    successes: np.ndarray, trials: np.ndarray, observed_labels: np.ndarray,
) -> dict[str, float | int | bool]:
    """Exhaustively permute whole-cage labels for the binomial deviance."""
    labels = np.asarray(observed_labels, dtype=int)
    unique = sorted(np.unique(labels))
    if len(unique) < 2 or unique != list(range(len(unique))):
        raise ValueError("Observed cage labels must contain consecutive condition codes")
    counts = [int(np.count_nonzero(labels == label)) for label in unique]
    labelings = _fixed_count_labelings(counts)
    observed = cage_binomial_deviance(successes, trials, labels)
    statistics = np.asarray(
        [cage_binomial_deviance(successes, trials, candidate) for candidate in labelings],
        dtype=float,
    )
    extreme = int(np.count_nonzero(statistics >= observed - 1e-12))
    return {
        "statistic": observed, "p_value": float(extreme / len(labelings)),
        "extreme_labelings": extreme, "enumerated_labelings": len(labelings),
        "permutation_exhaustive": True,
    }


def multipanel_sensitivity_statistics(cage_values: pd.DataFrame) -> pd.DataFrame:
    """Cage-unit ANOVA/MWU and exact activated-count sensitivity analyses."""
    pair_order = (("Water", "Sucrose"), ("Water", "Allulose"), ("Sucrose", "Allulose"))
    rows: list[dict[str, object]] = []
    rank_values = cage_values.loc[cage_values["sensitivity_family"].eq("cage_mean_anova")]
    for endpoint in [endpoint for _, endpoint, _, _ in PANEL_ENDPOINTS]:
        for region in REGIONS:
            subset = rank_values.loc[
                rank_values["endpoint"].eq(endpoint) & rank_values["region"].eq(region)
            ]
            groups = {
                condition: subset.loc[subset["condition"].eq(condition), "value"].to_numpy(float)
                for condition in CONDITIONS
            }
            if any(values.size < 2 or not np.isfinite(values).all() for values in groups.values()):
                raise RuntimeError(f"Cage-mean ANOVA is incomplete for {endpoint}/{region}")
            global_result = one_way_anova_receipt([groups[c] for c in CONDITIONS])
            rows.append({
                "endpoint": endpoint, "region": region,
                "sensitivity_family": "cage_mean_anova", "test": "One_way_ANOVA",
                "comparison": "Water|Sucrose|Allulose",
                "statistic": global_result["statistic"], "n_cages": int(len(subset)),
                "n_left": np.nan, "n_right": np.nan,
                "p_value_raw": global_result["p_value"],
                "effect_size": np.nan, "effect_size_definition": "",
                "permutation_unit": "whole cage",
                "p_value_method": "ordinary equal-variance one-way ANOVA F on cage means",
                "df_between": global_result["df_between"],
                "df_within": global_result["df_within"],
                "extreme_labelings": global_result["extreme_labelings"],
                "enumerated_labelings": global_result["enumerated_labelings"],
                "permutation_exhaustive": global_result["permutation_exhaustive"],
                "formula": "unweighted animal mean within cage; ordinary one-way ANOVA F",
            })
            pair_results = {
                pair: exact_permutation_mwu_receipt(groups[pair[0]], groups[pair[1]])
                for pair in pair_order
            }
            for left, right in pair_order:
                receipt = pair_results[(left, right)]
                rows.append({
                    "endpoint": endpoint, "region": region,
                    "sensitivity_family": "cage_mean_anova",
                    "test": "Mann_Whitney_U_two_sided", "comparison": f"{left}|{right}",
                    "statistic": receipt["statistic"],
                    "n_cages": int(groups[left].size + groups[right].size),
                    "n_left": int(groups[left].size), "n_right": int(groups[right].size),
                    "p_value_raw": receipt["p_value"],
                    "effect_size": cliffs_delta(groups[right], groups[left]),
                    "effect_size_definition": f"Cliff's delta ({right} minus {left})",
                    "permutation_unit": "whole cage",
                    "p_value_method": "exhaustive two-sided cage-label permutation of rank-sum/U",
                    "extreme_labelings": receipt["extreme_labelings"],
                    "enumerated_labelings": receipt["enumerated_labelings"],
                    "permutation_exhaustive": receipt["permutation_exhaustive"],
                    "formula": "unweighted animal mean within cage; two-sided rank-sum/U",
                })

    count_values = cage_values.loc[
        cage_values["sensitivity_family"].eq("cage_count_binomial_deviance")
    ]
    condition_code = {condition: index for index, condition in enumerate(CONDITIONS)}
    for region in REGIONS:
        subset = count_values.loc[count_values["region"].eq(region)].sort_values("cage")
        if len(subset) != count_values["cage"].nunique():
            raise RuntimeError(f"Activated-count cage table is incomplete for {region}")
        successes = subset["numerator"].to_numpy(float)
        trials = subset["denominator"].to_numpy(float)
        observed_labels = np.asarray([condition_code[str(value)] for value in subset["condition"]], dtype=int)
        global_result = exact_cage_label_binomial_deviance(successes, trials, observed_labels)
        rows.append({
            "endpoint": "activated_microglia_fraction", "region": region,
            "sensitivity_family": "cage_count_binomial_deviance",
            "test": "binomial_deviance", "comparison": "Water|Sucrose|Allulose",
            "statistic": global_result["statistic"], "n_cages": int(len(subset)),
            "n_left": np.nan, "n_right": np.nan,
            "p_value_raw": global_result["p_value"],
            "effect_size": np.nan, "effect_size_definition": "",
            "permutation_unit": "whole cage (numerator and denominator move together)",
            "p_value_method": "exhaustive cage-label permutation of binomial deviance",
            "extreme_labelings": global_result["extreme_labelings"],
            "enumerated_labelings": global_result["enumerated_labelings"],
            "permutation_exhaustive": global_result["permutation_exhaustive"],
            "formula": "D=2*sum_g[y_g*log(y_g/(n_g*p))+(n_g-y_g)*log((n_g-y_g)/(n_g*(1-p)))]; p=sum(y)/sum(n)",
        })
        pair_receipts: dict[tuple[str, str], dict[str, float | int | bool]] = {}
        pair_effects: dict[tuple[str, str], float] = {}
        for left, right in pair_order:
            pair = subset.loc[subset["condition"].isin((left, right))].sort_values("cage")
            labels = np.asarray([0 if value == left else 1 for value in pair["condition"]], dtype=int)
            pair_receipts[(left, right)] = exact_cage_label_binomial_deviance(
                pair["numerator"].to_numpy(float), pair["denominator"].to_numpy(float), labels,
            )
            rates = {}
            for condition in (left, right):
                selected = pair["condition"].eq(condition)
                rates[condition] = float(
                    pair.loc[selected, "numerator"].sum() / pair.loc[selected, "denominator"].sum()
                )
            pair_effects[(left, right)] = rates[right] - rates[left]
        for left, right in pair_order:
            pair = subset.loc[subset["condition"].isin((left, right))]
            receipt = pair_receipts[(left, right)]
            rows.append({
                "endpoint": "activated_microglia_fraction", "region": region,
                "sensitivity_family": "cage_count_binomial_deviance",
                "test": "binomial_deviance_two_condition", "comparison": f"{left}|{right}",
                "statistic": receipt["statistic"], "n_cages": int(len(pair)),
                "n_left": int(pair["condition"].eq(left).sum()),
                "n_right": int(pair["condition"].eq(right).sum()),
                "p_value_raw": receipt["p_value"],
                "effect_size": pair_effects[(left, right)],
                "effect_size_definition": f"aggregate activated fraction ({right} minus {left})",
                "permutation_unit": "whole cage (numerator and denominator move together)",
                "p_value_method": "exhaustive two-condition cage-label permutation of binomial deviance",
                "extreme_labelings": receipt["extreme_labelings"],
                "enumerated_labelings": receipt["enumerated_labelings"],
                "permutation_exhaustive": receipt["permutation_exhaustive"],
                "formula": "two-condition version of D; cage count vectors are never split",
            })

    result = pd.DataFrame(rows)
    for family, expected_global in (
        ("cage_mean_anova", len(PANEL_ENDPOINTS) * len(REGIONS)),
        ("cage_count_binomial_deviance", len(REGIONS)),
    ):
        selected = result["sensitivity_family"].eq(family) & result["comparison"].eq("Water|Sucrose|Allulose")
        if int(selected.sum()) != expected_global:
            raise RuntimeError(f"{family} global family has {int(selected.sum())} rows; expected {expected_global}")
    exact_rows = ~result["test"].eq("One_way_ANOVA")
    if not result.loc[exact_rows, "permutation_exhaustive"].astype(bool).all():
        raise RuntimeError("Every exact cage sensitivity must use exhaustive label enumeration")
    if not (
        result.loc[exact_rows, "extreme_labelings"].astype(int)
        <= result.loc[exact_rows, "enumerated_labelings"].astype(int)
    ).all():
        raise RuntimeError("Invalid cage sensitivity enumeration receipt")
    return result.sort_values(
        ["sensitivity_family", "endpoint", "region", "test", "comparison"]
    ).reset_index(drop=True)


def draw_region_plot(
    ax: plt.Axes,
    data: pd.DataFrame,
    endpoint: str,
    region: str,
    ylabel: str,
    show_ylabel: bool,
    seed: int,
    font_size: float = 9.8,
) -> list[dict]:
    """Render one region using the Figure 4D bar/dot/statistics grammar."""
    result, rows = endpoint_statistics(data, endpoint, region)
    positions = np.arange(3, dtype=float)
    groups = result["groups"]
    means = [float(np.mean(groups[condition])) if groups[condition].size else np.nan for condition in CONDITIONS]
    sds = [float(np.std(groups[condition], ddof=1)) if groups[condition].size >= 2 else 0.0 for condition in CONDITIONS]
    ax.bar(
        positions, means, yerr=sds, width=0.52,
        color=[COLORS[condition] for condition in CONDITIONS],
        edgecolor="black", linewidth=1.7, capsize=4,
        error_kw={"elinewidth": 1.6, "capthick": 1.6}, zorder=2,
    )
    rng = np.random.default_rng(int(seed))
    tops: list[float] = []
    for index, condition in enumerate(CONDITIONS):
        values = groups[condition]
        outliers = tukey_outliers(values)
        jitter = rng.uniform(-0.065, 0.065, values.size)
        ax.scatter(
            index + jitter[~outliers], values[~outliers],
            facecolor=COLORS[condition], edgecolor="black",
            linewidth=1.0, s=32, zorder=4,
        )
        if outliers.any():
            ax.scatter(
                index + jitter[outliers], values[outliers], marker="x",
                color="black", linewidth=1.8, s=42, zorder=5,
            )
        tops.append(max(
            float(means[index] + sds[index]) if math.isfinite(means[index]) else 0.0,
            float(np.max(values)) if values.size else 0.0,
        ))
    y_top = max(max(tops), 1e-4) * 1.85
    ax.set_ylim(0, y_top)
    for index, condition in enumerate(CONDITIONS):
        if groups[condition].size:
            ax.text(
                index, tops[index] + 0.025 * y_top,
                f"n={groups[condition].size}", ha="center", va="bottom",
                fontsize=font_size - 1.2, fontweight="bold",
            )
    pairwise_lines = (
        "MWU W–S p = " + format_p_token(result["pairs_raw"][("Water", "Sucrose")]),
        "MWU W–A p = " + format_p_token(result["pairs_raw"][("Water", "Allulose")]),
        "MWU S–A p = " + format_p_token(result["pairs_raw"][("Sucrose", "Allulose")]),
    )
    ax.text(
        0.025, 0.975, "One-way ANOVA p = " + format_p_token(result["global_p"]),
        transform=ax.transAxes, ha="left", va="top",
        fontsize=max(8.8, font_size + 0.6), fontweight="bold",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
        zorder=11,
    )
    ax.text(
        0.025, 0.865, "\n".join(pairwise_lines),
        transform=ax.transAxes, ha="left", va="top",
        fontsize=max(6.4, font_size - 2.0), fontweight="bold",
        linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.76, "pad": 1.2},
        zorder=10,
    )
    ax.set_title(region, fontsize=font_size + 1.5, fontweight="bold", pad=4)
    ax.set_xticks(positions)
    ax.set_xticklabels(CONDITIONS, rotation=24, ha="right", fontsize=font_size, fontweight="bold")
    ax.set_ylabel(ylabel if show_ylabel else "", fontsize=font_size + 0.5, fontweight="bold")
    ax.tick_params(axis="both", labelsize=font_size - 0.5, width=1.5, length=4)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.6)
    ax.grid(False)
    return rows


def draw_endpoint_pair(
    fig: plt.Figure,
    spec,
    panel: tuple[str, str, str, str],
    data: pd.DataFrame,
    seed: int,
    show_condition_labels: bool = True,
) -> list[dict]:
    axes_before = {id(axis) for axis in fig.axes}
    letter, endpoint, title, ylabel = panel
    inner = spec.subgridspec(
        2, len(REGIONS),
        height_ratios=(0.18, 0.82), hspace=0.02, wspace=0.34,
    )
    heading_ax = fig.add_subplot(inner[0, :])
    heading_ax.set_axis_off()
    heading_ax.text(
        0.0, 0.95, letter, transform=heading_ax.transAxes,
        fontsize=20, fontweight="bold", ha="left", va="top",
    )
    heading_ax.text(
        0.5, 0.95, title, transform=heading_ax.transAxes,
        fontsize=12.5, fontweight="bold", ha="center", va="top",
    )
    axes = [fig.add_subplot(inner[1, index]) for index in range(len(REGIONS))]
    rows: list[dict] = []
    for index, (ax, region) in enumerate(zip(axes, REGIONS)):
        rows.extend(draw_region_plot(
            ax, data, endpoint, region, ylabel, index == 0,
            seed=int(seed) + index, font_size=9.8,
        ))
    if not show_condition_labels:
        for axis in fig.axes:
            if id(axis) not in axes_before:
                axis.tick_params(axis="x", which="both", labelbottom=False)
    return rows


def statistics_lines(stats: pd.DataFrame, endpoint: str, label: str) -> list[str]:
    lines = []
    for region in REGIONS:
        subset = stats.loc[(stats["endpoint"] == endpoint) & (stats["region"] == region)]
        raw_lookup = {(row.test, row.comparison): row.p_value_raw for row in subset.itertuples()}
        lines.append(
            f"{label}, {region}: one-way ANOVA {format_p(float(raw_lookup.get(('One_way_ANOVA', 'Water|Sucrose|Allulose'), np.nan)))}, "
            f"MWU W-S {format_p(float(raw_lookup.get(('Mann_Whitney_U_two_sided', 'Water|Sucrose'), np.nan)))}, "
            f"MWU W-A {format_p(float(raw_lookup.get(('Mann_Whitney_U_two_sided', 'Water|Allulose'), np.nan)))}, "
            f"MWU S-A {format_p(float(raw_lookup.get(('Mann_Whitney_U_two_sided', 'Sucrose|Allulose'), np.nan)))}."
        )
    return lines


def sensitivity_statistics_lines(
    sensitivity: pd.DataFrame, endpoint: str, label: str,
) -> list[str]:
    lines: list[str] = []
    for region in REGIONS:
        subset = sensitivity.loc[
            sensitivity["endpoint"].eq(endpoint)
            & sensitivity["region"].eq(region)
            & sensitivity["sensitivity_family"].eq("cage_mean_anova")
        ]
        lookup = {(row.test, row.comparison): row for row in subset.itertuples()}
        global_row = lookup.get(("One_way_ANOVA", "Water|Sucrose|Allulose"))
        if global_row is None:
            raise RuntimeError(f"Missing cage-mean global sensitivity for {endpoint}/{region}")
        lines.append(
            f"{label}, {region}: one-way ANOVA {format_p(float(global_row.p_value_raw))}, "
            f"MWU W-S {format_p(float(lookup[('Mann_Whitney_U_two_sided', 'Water|Sucrose')].p_value_raw))}, "
            f"MWU W-A {format_p(float(lookup[('Mann_Whitney_U_two_sided', 'Water|Allulose')].p_value_raw))}, "
            f"MWU S-A {format_p(float(lookup[('Mann_Whitney_U_two_sided', 'Sucrose|Allulose')].p_value_raw))}."
        )
    return lines


def activated_count_statistics_lines(sensitivity: pd.DataFrame) -> list[str]:
    lines: list[str] = []
    for region in REGIONS:
        subset = sensitivity.loc[
            sensitivity["region"].eq(region)
            & sensitivity["sensitivity_family"].eq("cage_count_binomial_deviance")
        ]
        lookup = {row.comparison: row for row in subset.itertuples()}
        global_row = lookup.get("Water|Sucrose|Allulose")
        if global_row is None:
            raise RuntimeError(f"Missing activated-count sensitivity for {region}")
        lines.append(
            f"Count H, {region}: D exact {format_p(float(global_row.p_value_raw))}, "
            f"D W-S {format_p(float(lookup['Water|Sucrose'].p_value_raw))}, "
            f"D W-A {format_p(float(lookup['Water|Allulose'].p_value_raw))}, "
            f"D S-A {format_p(float(lookup['Sucrose|Allulose'].p_value_raw))}."
        )
    return lines


def write_multipanel_sensitivity_outputs(
    source_data_dir: Path,
    legends_dir: Path,
    provenance_dir: Path,
    values_path: Path,
    cage_values: pd.DataFrame,
    denominator_audit: pd.DataFrame,
    sensitivity_stats: pd.DataFrame,
) -> dict[str, Path]:
    """Write deterministic sensitivity tables, methods, and provenance."""
    stem = "Figure_GFAP_Iba1_microglia_ARC_ME"
    paths = {
        "sensitivity_statistics": source_data_dir / f"{stem}_multipanel_sensitivity_statistics.csv",
        "sensitivity_cage_values": source_data_dir / f"{stem}_multipanel_sensitivity_cage_values.csv",
        "denominator_audit": source_data_dir / f"{stem}_multipanel_denominator_audit.csv",
        "sensitivity_methods": legends_dir / f"{stem}_multipanel_sensitivity_methods.txt",
        "sensitivity_provenance": provenance_dir / f"{stem}_multipanel_sensitivity_provenance.json",
    }
    sensitivity_stats.to_csv(paths["sensitivity_statistics"], index=False)
    cage_values.to_csv(paths["sensitivity_cage_values"], index=False)
    denominator_audit.to_csv(paths["denominator_audit"], index=False)
    methods = "\n".join((
        "Figure 5 multipanel sensitivity methods",
        "",
        "Primary visible inference uses animals: ordinary one-way ANOVA globally and raw exact two-sided MWU for W-S, W-A, and S-A.",
        "Cage-mean sensitivity uses an unweighted mean of finite animal values within each cage, then ordinary one-way ANOVA and raw exact two-sided MWU.",
        "Activated-count sensitivity is restricted to panel H. Activated (Activated + Amoeboid) and total QC-accepted microglia counts are summed within cage. Raw exact whole-cage-label binomial-deviance tests are global and pairwise.",
        "No outcome-dependent filtering is used. Positive low denominators are retained. Zero denominators make an animal ratio intrinsically non-estimable but its valid 0/0 count pair remains in the cage sum; every cage must have a positive aggregate denominator.",
        "All exact tests enumerate every fixed-count labeling; numerator and denominator move together with the cage.",
        "",
        "Cage count formula:",
        "D = 2*sum_g[y_g*log(y_g/(n_g*p)) + (n_g-y_g)*log((n_g-y_g)/(n_g*(1-p)))], p = sum(y)/sum(n).",
        "",
    ))
    paths["sensitivity_methods"].write_text(methods, encoding="utf-8")
    source_bytes = values_path.read_bytes()
    provenance = {
        "schema": "figure6_multipanel_sensitivity_v2_raw_exact",
        "animal_values_path": f"analysis/{values_path.name}",
        "animal_values_sha256": hashlib.sha256(source_bytes).hexdigest(),
        "renderer_path": "Fig5/02_make_figure_5_gfap_iba1_microglia.py",
        "renderer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "conditions": list(CONDITIONS), "regions": list(REGIONS),
        "cages_by_condition": (
            cage_values.loc[cage_values["sensitivity_family"].eq("cage_count_binomial_deviance")]
            .drop_duplicates(["cage", "condition"])
            .groupby("condition", observed=True)["cage"].nunique().sort_index().to_dict()
        ),
        "animal_region_rows": int(len(denominator_audit)),
        "zero_denominator_animal_region_rows": int((~denominator_audit["ratio_estimable"]).sum()),
        "outcome_dependent_filtering": False,
        "primary_inference": "ordinary one-way ANOVA global and raw exact two-sided MWU W-S/W-A/S-A",
        "cage_mean_inference": "ordinary one-way ANOVA on cage means and raw exact two-sided MWU W-S/W-A/S-A",
        "activated_count_inference": "raw exact cage-label binomial deviance global and W-S/W-A/S-A",
        "pairwise_display": "raw exhaustive two-sided p values; W-S, W-A, S-A",
        "multiplicity_adjustment": "none",
        "canonical_generation_gate": "exactly 61/61 native GFAP/Iba1 HIL receipts",
    }
    paths["sensitivity_provenance"].write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return paths


def write_legends(
    output_dir: Path,
    values_path: Path,
    metadata_path: Path,
    iba1_animal: str,
    iba1_section: str,
    animal: str,
    section: str,
    micro_uid: str,
    micro_bbox: tuple[int, int, int, int],
    gfap_center: tuple[float, float],
    examples: dict[str, dict[str, object]],
    stats: pd.DataFrame,
    sensitivity_stats: pd.DataFrame,
    metadata: dict,
    data: pd.DataFrame,
    iba1_provenance: dict[str, Any],
    gfap_provenance: dict[str, Any],
    panel_f_audit: dict[str, Any],
) -> None:
    counts = data[["animal_id", "cond"]].drop_duplicates().groupby("cond", observed=True)["animal_id"].nunique().to_dict()
    training_n = int(metadata.get("training_candidates_n", 10848) or 10848)
    example_qc = " ".join(
        f"{state}: {'1 selected' if entry['row'] is not None else '0 selected'} from {entry['total_predictions']} predictions ({entry['reason']})."
        for state, entry in examples.items()
    )
    def accepted_subset_phrase(provenance: dict[str, Any]) -> str:
        colors = {"ARC": "cyan", "ME": "orange", "VMN": "magenta"}
        drawn = tuple(
            region for region in ("ARC", "ME", "VMN")
            if provenance.get("region_status", {}).get(region) == "drawn"
        )
        if not drawn:
            raise RuntimeError("Publication representative has no HIL-delineated ROI")
        keyed = ", ".join(f"{region} {colors[region]}" for region in drawn)
        return (
            f"the exact accepted subset {', '.join(drawn)} ({keyed}) with bold colored "
            "contours and a dark contrast halo; omitted ROIs remain missing"
        )

    iba1_subset = accepted_subset_phrase(iba1_provenance)
    gfap_subset = accepted_subset_phrase(gfap_provenance)
    panel_text = {
        "A_Iba1_microscopy": f"(A) Representative Iba1/DAPI section from {iba1_animal}, {iba1_section}. Iba1 is magenta and DAPI is blue; the high-contrast yellow {PANEL_A_INSET_SIZE} x {PANEL_A_INSET_SIZE}-px source box (x={micro_bbox[0]}-{micro_bbox[2]}, y={micro_bbox[1]}-{micro_bbox[3]}) is centered on QC-valid cell {micro_uid} and matches the enlarged inset. The immediately adjacent right-hand ROI view uses the identical source DAPI and source-pixel bounds and shows {iba1_subset}, HIL-accepted at {iba1_provenance['accepted_at']} from {iba1_provenance['accepted_source']}; mask SHA-256 {iba1_provenance['mask_sha256']}.",
        "B_GFAP_microscopy": f"(B) Representative GFAP/DAPI section from {animal}, {section}. GFAP is cyan-green and DAPI is blue; the inset is centered at x={gfap_center[0]:.1f}, y={gfap_center[1]:.1f} px. The immediately adjacent right-hand ROI view uses the identical source DAPI and source-pixel bounds and shows {gfap_subset}, HIL-accepted at {gfap_provenance['accepted_at']} from {gfap_provenance['accepted_source']}; mask SHA-256 {gfap_provenance['mask_sha256']}.",
        "C_microglia_CNN": f"(C) Implemented TinyMorphCNN. Three-channel 96 x 96 inputs encode background-corrected Iba1 intensity, the cell mask, and skeleton/distance morphology. Four convolutional blocks, global-average pooling, dropout, and softmax predict Ramified, Rod-like, Activated, or Amoeboid morphology. The hybrid run trained on {training_n:,} high-confidence morphometric pseudo-labels using a whole-animal 80/20 split. Displayed patches are deterministic within-class representatives that passed classifier, object-size, morphometry, centering, connected-component, and border-integrity QC; a neutral placeholder is shown when no prediction passes, with no out-of-QC fallback. {example_qc}",
        "D_microglia_state_composition": f"(D) Within-region CNN-predicted microglial-state percentages for the representative Iba1 section from {iba1_animal}, {iba1_section}. Grouped bars show Ramified, Rod-like, Activated, and Amoeboid assignments separately for ARC, ME, and VMN; n denotes QC-accepted cells. This representative-section summary is descriptive and no inferential test is applied.",
        "E_GFAP_workflow": "(E) GFAP quantification workflow for the representative raw GFAP/DAPI field and final GFAP-positive mask in panel B. P0.5-P99.8 normalization is used only to create the positive mask; a 35-um grey-opening background is subtracted and the section threshold is max(0.8 x Otsu, median + 3 x normalized MAD). Objects smaller than 8 um2 are removed. GFAP-positive area is normalized by ROI pixels. Quantitative integrated intensity is background-corrected on the unscaled raw exported channel, then normalized by ROI area before animal-level aggregation.",
        "F_DAPI_Iba1_GFAP_mask_cartoon": f"(F) Vertical channel-separated analysis-mask cartoons from the exact panel A field ({iba1_animal}, {iba1_section}; source bounds {panel_f_audit['bounds_px']}). DAPI nuclei are blue ({panel_f_audit['instance_counts']['DAPI_instances']:,} labeled instances; {panel_f_audit['positive_pixels']['DAPI']:,} positive pixels), final QC-accepted Iba1-positive cells are red-magenta ({panel_f_audit['instance_counts']['Iba1_QC_accepted_cells']:,} labeled cells; {panel_f_audit['positive_pixels']['Iba1']:,} positive pixels), and the final threshold-derived GFAP-positive signal is green ({panel_f_audit['positive_pixels']['GFAP']:,} positive pixels; pixels, not cells). Every channel uses the same DAPI-derived tissue visualization envelope and the same pixel-identical HIL-accepted/final ARC, ME, and VMN mask; the tissue envelope is visualization support, not an anatomical ROI or quantitative denominator.",
        "G_Iba1_signal": "(G) Raw-scale background-corrected Iba1 signal per analyzed regional area.",
        "H_activated_microglia": "(H) Exploratory fraction of segmented Iba1-positive cells assigned to the activated morphology grouping (Activated plus Amoeboid) by the pseudo-label-trained classifier.",
        "I_GFAP_signal": "(I) Raw-scale background-corrected GFAP signal per analyzed regional area.",
    }
    shared = (
        f"ARC, ME, and VMN are shown in canonical annotation-code order in all quantitative panels. Bars are animal means and error bars are SD; circles are individual biological animals and crosses identify values outside within-condition Tukey 1.5-IQR fences. All finite animal values, including marked outliers, were retained. Primary global values are ordinary one-way ANOVA p values. Pairwise W-S, W-A, and S-A values are raw exact two-sided Mann–Whitney U p values. Cage sensitivity uses unweighted cage means with one-way ANOVA and exact MWU. Panel H also uses raw exact cage-label binomial deviance on cage-summed activated/total microglia counts. Positive low denominators are retained; a zero denominator is an intrinsically undefined animal ratio but its 0/0 counts remain in the cage sum. No outcome-dependent filtering is used (W, Water; S, Sucrose; A, Allulose). Animals available anywhere in the dataset were Water n={counts.get('Water', 0)}, Sucrose n={counts.get('Sucrose', 0)}, and Allulose n={counts.get('Allulose', 0)}; actual finite animal Ns can be smaller for an endpoint/region."
    )
    boundary = (
        "ARC/ME/VMN provenance: all statistics and representative fields are gated on exactly 61/61 native image-by-image HIL reviews; registration-derived or unreviewed automated masks are forbidden. The displayed accepted segmentations and the final analysis masks are pixel-identical and use the same source-pixel bounds as their adjacent microscopy fields. Classifier provenance: this was a hybrid pseudo-label CNN, not an expert-annotated supervised classifier. Its validation metric measures agreement with K-means morphometric pseudo-labels and must not be interpreted as biological accuracy. Panels C, D, and H are classifier-dependent and exploratory: panel C documents the classifier and examples, panel D is a descriptive representative-section composition, and panel H is the exploratory quantitative activated fraction. Microglial cells, DAPI nuclei, fluorescence signals, and area denominators are assigned with the same final ARC/ME/VMN mask; the projected-wall classification is retained as an audit vote. Raw exported intensity remains vulnerable to acquisition/batch differences. The Water cohort is smaller and confounded with sex/genotype in the available material."
    )
    lines = ["Figure GFAP/Iba1 microglia ARC/ME/VMN.", "", *panel_text.values(), "", shared, "", boundary, "", f"Animal-level data source: analysis/{values_path.name}", f"Classifier metadata: analysis/{metadata_path.name}", "", "Statistics (all finite animal-level values):"]
    for letter, endpoint, title, _ in PANEL_ENDPOINTS:
        lines.extend(statistics_lines(stats, endpoint, f"Panel {letter}, {title}"))
    lines.extend(("", "Cage-mean sensitivity:"))
    for letter, endpoint, title, _ in PANEL_ENDPOINTS:
        lines.extend(sensitivity_statistics_lines(sensitivity_stats, endpoint, f"Cage {letter}, {title}"))
    lines.extend(("", "Activated-count sensitivity:"))
    lines.extend(activated_count_statistics_lines(sensitivity_stats))
    master = "\n".join(lines) + "\n"
    stem = "Figure_GFAP_Iba1_microglia_ARC_ME"
    (output_dir / f"{stem}_LEGEND.txt").write_text(master, encoding="utf-8")
    (output_dir / f"{stem}_caption.txt").write_text(master, encoding="utf-8")
    for slug, text in panel_text.items():
        (output_dir / f"{stem}_panel_{slug}_LEGEND.txt").write_text("Figure GFAP/Iba1 microglia ARC/ME/VMN.\n\n" + text + "\n\n" + shared + "\n", encoding="utf-8")
    spanish_panel_text = {
        "A_Iba1_microscopy": (
            "(A) Sección representativa Iba1/DAPI con el campo de origen y la máscara "
            "ARC/ME/VMN aceptada mediante HIL en coordenadas idénticas. Iba1 se "
            "muestra en magenta, DAPI en azul y el recuadro amplía una célula válida."
        ),
        "B_GFAP_microscopy": (
            "(B) Sección representativa GFAP/DAPI con su máscara anatómica aceptada "
            "mediante HIL en los mismos píxeles. GFAP se muestra en verde cian y DAPI en azul."
        ),
        "C_microglia_CNN": (
            "(C) TinyMorphCNN implementada con entradas de intensidad Iba1, máscara "
            "celular y morfología de esqueleto/distancia. Las clases Ramificada, En "
            "bastón, Activada y Ameboide proceden de pseudoetiquetas morfométricas, no "
            "de etiquetas biológicas expertas."
        ),
        "D_microglia_state_composition": (
            "(D) Porcentajes de estados microgliales predichos por la CNN dentro de ARC, "
            "ME y VMN en la sección Iba1 representativa. Es un resumen descriptivo."
        ),
        "E_GFAP_workflow": (
            "(E) Flujo de cuantificación GFAP desde el campo GFAP/DAPI crudo hasta la "
            "máscara positiva final. La intensidad integrada se corrige por fondo y se "
            "normaliza por el área regional antes de agregar por animal."
        ),
        "F_DAPI_Iba1_GFAP_mask_cartoon": (
            "(F) Esquemas verticales de las máscaras DAPI, Iba1 y GFAP del mismo campo "
            "del panel A. Todas usan la misma máscara final ARC/ME/VMN aceptada."
        ),
        "G_Iba1_signal": (
            "(G) Señal Iba1 corregida por fondo y normalizada por área regional analizada."
        ),
        "H_activated_microglia": (
            "(H) Fracción exploratoria de células Iba1 positivas asignadas al grupo "
            "Activada más Ameboide por el clasificador entrenado con pseudoetiquetas."
        ),
        "I_GFAP_signal": (
            "(I) Señal GFAP corregida por fondo y normalizada por área regional analizada."
        ),
    }
    spanish_shared = (
        "ARC, ME y VMN se muestran en el orden canónico del código de anotación. "
        "Las barras son medias animales con DE, los círculos son animales biológicos y "
        "las cruces señalan valores fuera de las cercas de Tukey, conservados en el "
        "análisis. Los valores globales usan ANOVA ordinario de una vía y las "
        "comparaciones Agua–Sacarosa, Agua–Alulosa y Sacarosa–Alulosa usan "
        "Mann–Whitney exacta bilateral sin ajuste. La sensibilidad emplea medias por "
        "jaula y, para el recuento activado, desviancia binomial exacta con etiquetas "
        "de jaula. Las 61 de 61 secciones usan máscaras ARC/ME/VMN nativas revisadas "
        "imagen por imagen. Sin embargo, los paneles C, D y H dependen de un clasificador "
        "entrenado con pseudoetiquetas y siguen siendo exploratorios; su exactitud "
        "biológica no ha sido validada independientemente."
    )
    spanish_master = (
        "Figura 5. Microglía Iba1 y señal GFAP en ARC/ME/VMN.\n\n"
        + "\n".join(spanish_panel_text.values())
        + "\n\n" + spanish_shared + "\n"
    )
    (output_dir / f"{stem}_LEGEND_spanish.txt").write_text(
        spanish_master, encoding="utf-8"
    )
    for slug, text in spanish_panel_text.items():
        (output_dir / f"{stem}_panel_{slug}_LEGEND_spanish.txt").write_text(
            "Figura 5. Microglía Iba1 y señal GFAP en ARC/ME/VMN.\n\n"
            + text + "\n\n" + spanish_shared + "\n",
            encoding="utf-8",
        )



def export_panel_crops(
    fig: plt.Figure,
    panel_axes: dict[str, list[plt.Axes]],
    output_dir: Path,
    suffix: str = "",
) -> None:
    """Crop exact live panels in one language from the same figure artists."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    all_axes = list(fig.axes)
    all_figure_texts = list(fig.texts)
    panel_texts: dict[str, list[Text]] = {}
    panel_bboxes: dict[str, Bbox] = {}

    panel_letters = tuple("ABCDEF") + tuple(panel[0] for panel in PANEL_ENDPOINTS)
    quantitative_letters = {panel[0] for panel in PANEL_ENDPOINTS}
    for letter in panel_letters:
        axes = panel_axes.get(letter, [])
        explicit_texts: list[Text] = []
        for axis in axes:
            explicit_texts.extend(getattr(axis, "_panel_export_texts", ()))
        panel_texts[letter] = explicit_texts
        boxes = []
        for axis in axes:
            if axis.get_visible():
                box = axis.get_tightbbox(renderer)
                if box is not None and np.isfinite(box.extents).all():
                    boxes.append(box)
        for artist in explicit_texts:
            box = artist.get_window_extent(renderer)
            if box is not None and np.isfinite(box.extents).all():
                boxes.append(box)
        if not boxes:
            raise RuntimeError(f"No drawable artists recorded for panel {letter}")
        pixel_box = Bbox.union(boxes)
        inch_box = pixel_box.transformed(fig.dpi_scale_trans.inverted())
        # Matplotlib's 100-dpi tight-bbox estimate can be a little narrower
        # than the final 600-dpi text geometry. Quantitative panels need a
        # deliberately larger left allowance for the vertical endpoint label.
        left_pad_inches = 0.28 if letter in quantitative_letters else 0.08
        top_pad_inches = 0.10 if letter in quantitative_letters else 0.08
        panel_bboxes[letter] = Bbox.from_extents(
            inch_box.x0 - left_pad_inches,
            inch_box.y0 - 0.08,
            inch_box.x1 + 0.08,
            inch_box.y1 + top_pad_inches,
        )

    for letter in panel_letters:
        target_axes = set(panel_axes[letter])
        target_texts = set(panel_texts[letter])
        axis_visibility = [(axis, axis.get_visible()) for axis in all_axes]
        text_visibility = [(artist, artist.get_visible()) for artist in all_figure_texts]
        try:
            for axis in all_axes:
                axis.set_visible(axis in target_axes)
            for artist in all_figure_texts:
                artist.set_visible(artist in target_texts)
            stem = output_dir / f"Figure_GFAP_Iba1_microglia_ARC_ME_Panel_{letter}{suffix}"
            fig.savefig(stem.with_suffix(".pdf"), format="pdf", dpi=600, bbox_inches=panel_bboxes[letter], facecolor="white")
            fig.savefig(stem.with_suffix(".png"), format="png", dpi=600, bbox_inches=panel_bboxes[letter], facecolor="white")
        finally:
            for axis, visible in axis_visibility:
                axis.set_visible(visible)
            for artist, visible in text_visibility:
                artist.set_visible(visible)


def main() -> int:
    args = parse_args()
    args.analysis_dir = args.analysis_dir.expanduser().resolve()
    args.input_root = args.input_root.expanduser().resolve()
    manual_root = args.manual_region_dir.expanduser().resolve()
    output_arg = args.output_dir.expanduser()
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    invalid_formats = sorted(set(formats).difference({"pdf", "png"}))
    if not formats or invalid_formats:
        raise RuntimeError(
            f"Publication formats must be a nonempty subset of pdf,png; invalid={invalid_formats}"
        )
    if int(args.dpi) <= 0 or float(args.um_per_px) <= 0:
        raise RuntimeError("--dpi and --um-per-px must both be positive")

    # Scientific/provenance preflight is deliberately first. A stale automated
    # analysis or incomplete 61-section review must not create, clean, or modify
    # the publication output directory, including when --force was requested.
    coverage, section_summary = validate_all_hil_analysis(args.analysis_dir, manual_root)
    values_path = args.analysis_dir / "per_animal_region_summary.csv"
    if not values_path.is_file():
        raise RuntimeError(f"Missing canonical animal/region summary: {values_path}")
    canonical_data = pd.read_csv(values_path)
    required = {"animal_id", "region", "cond", *(panel[1] for panel in PANEL_ENDPOINTS)}
    missing = sorted(required.difference(canonical_data.columns))
    if missing:
        raise ValueError("Missing animal-level columns: " + ", ".join(missing))
    missing_regions = sorted(
        set(REGIONS).difference(canonical_data["region"].astype(str).unique())
    )
    if missing_regions:
        raise ValueError(
            "Canonical animal summary is missing quantitative regions: "
            + ", ".join(missing_regions)
        )

    cells_path = args.analysis_dir / "per_microglia_cell_measurements.csv"
    if not cells_path.is_file():
        raise RuntimeError(f"Missing current per-cell analysis table: {cells_path}")
    cells = pd.read_csv(cells_path)
    recomputed_animal_data = reaggregate_animal_regions(section_summary, cells)
    verify_canonical_animal_summary(canonical_data, recomputed_animal_data)
    data = canonical_data.loc[
        canonical_data["region"].isin(REGIONS)
        & canonical_data["cond"].isin(CONDITIONS)
    ].copy()
    section_dir = args.input_root / args.animal / args.section
    processed_dir = args.analysis_dir / "processed_sections" / args.animal / args.section
    iba1_section_dir = args.input_root / args.iba1_animal / args.iba1_section
    iba1_processed_dir = (
        args.analysis_dir / "processed_sections" / args.iba1_animal / args.iba1_section
    )
    iba1_dapi_path = iba1_section_dir / f"Image_{args.iba1_section}_DAPI.tif"
    iba1_path = iba1_section_dir / f"Image_{args.iba1_section}_Iba1.tif"
    panel_f_raw_gfap_path = iba1_section_dir / f"Image_{args.iba1_section}_GFAP.tif"
    panel_f_dapi_labels_path = (
        iba1_processed_dir
        / f"{args.iba1_animal}_{args.iba1_section}_DAPI_labels.tif"
    )
    panel_f_iba1_cells_path = (
        iba1_processed_dir
        / f"{args.iba1_animal}_{args.iba1_section}_Iba1_cells.tif"
    )
    panel_f_gfap_positive_path = (
        iba1_processed_dir
        / f"{args.iba1_animal}_{args.iba1_section}_GFAP_positive.tif"
    )
    gfap_dapi_path = section_dir / f"Image_{args.section}_DAPI.tif"
    gfap_path = section_dir / f"Image_{args.section}_GFAP.tif"
    positive_path = processed_dir / f"{args.animal}_{args.section}_GFAP_positive.tif"

    # Reload the two representative native pairs for display. The preflight above
    # already verified every one of the 61 final-analysis masks pixel-for-pixel.
    iba1_labels, iba1_provenance = load_native_accepted_arc_me(
        manual_root, args.iba1_animal, args.iba1_section, iba1_dapi_path
    )
    labels, gfap_provenance = load_native_accepted_arc_me(
        manual_root, args.animal, args.section, gfap_dapi_path
    )
    for panel_name, provenance, rep_animal, rep_section, rep_flags in (
        ("A", iba1_provenance, args.iba1_animal, args.iba1_section,
         "--iba1-animal and --iba1-section"),
        ("B", gfap_provenance, args.animal, args.section,
         "--animal and --section"),
    ):
        missing_representative_regions = [
            region for region in REGIONS
            if provenance.get("region_status", {}).get(region) != "drawn"
        ]
        if missing_representative_regions:
            candidates = fully_delineated_sections(manual_root)
            available = (
                f"{len(candidates)} accepted sections carry all three: "
                + ", ".join(f"{animal} {section}" for animal, section in candidates[:12])
                + (" ..." if len(candidates) > 12 else "")
                if candidates else
                "No accepted section in this review carries all three regions."
            )
            raise RuntimeError(
                f"Panel {panel_name} representative {rep_animal} {rep_section} must have "
                f"reviewed ARC, ME, and VMN; missing "
                f"{', '.join(missing_representative_regions)}. Either delineate the missing "
                f"region on that section in the reviewer, or name a different representative "
                f"with {rep_flags}. {available}"
            )
    iba1_final_mask_path, iba1_geometry = verify_final_analysis_mask_identity(
        args.analysis_dir,
        args.iba1_animal,
        args.iba1_section,
        iba1_labels,
        iba1_provenance,
    )
    gfap_final_mask_path, gfap_geometry = verify_final_analysis_mask_identity(
        args.analysis_dir,
        args.animal,
        args.section,
        labels,
        gfap_provenance,
    )
    panel_f_required_paths = (
        panel_f_raw_gfap_path,
        panel_f_dapi_labels_path,
        panel_f_iba1_cells_path,
        panel_f_gfap_positive_path,
        iba1_final_mask_path,
    )
    missing_panel_f_paths = [path for path in panel_f_required_paths if not path.is_file()]
    if missing_panel_f_paths:
        raise RuntimeError(
            "Panel F requires current FR722/S03 raw channels and final analysis masks; "
            f"missing={[str(path) for path in missing_panel_f_paths]}"
        )
    iba1_dapi = tifffile.imread(iba1_dapi_path)
    iba1 = tifffile.imread(iba1_path)
    panel_f_raw_gfap = tifffile.imread(panel_f_raw_gfap_path)
    gfap_dapi = tifffile.imread(gfap_dapi_path)
    gfap = tifffile.imread(gfap_path)
    positive = tifffile.imread(positive_path) > 0
    for name, image, dapi in (
        ("Iba1", iba1, iba1_dapi),
        ("GFAP", gfap, gfap_dapi),
        ("GFAP positive mask", positive, gfap_dapi),
    ):
        if marker_signal(image).shape != marker_signal(dapi).shape:
            raise RuntimeError(
                f"Representative {name}/DAPI shape mismatch: "
                f"{marker_signal(image).shape} != {marker_signal(dapi).shape}"
            )
    panel_f_dapi_labels = np.squeeze(tifffile.imread(panel_f_dapi_labels_path))
    panel_f_iba1_cells = np.squeeze(tifffile.imread(panel_f_iba1_cells_path))
    panel_f_gfap_positive = np.squeeze(
        tifffile.imread(panel_f_gfap_positive_path)
    ) > 0
    panel_f_regions = np.squeeze(tifffile.imread(iba1_final_mask_path))
    panel_f_layers = {
        "DAPI labels": panel_f_dapi_labels,
        "Iba1 cells": panel_f_iba1_cells,
        "GFAP positive": panel_f_gfap_positive,
        "ARC/ME/VMN regions": panel_f_regions,
    }
    if any(np.asarray(layer).ndim != 2 for layer in panel_f_layers.values()):
        raise RuntimeError(
            "Panel F final analysis masks must each squeeze to exactly two dimensions: "
            + ", ".join(
                f"{name}={np.asarray(layer).shape}"
                for name, layer in panel_f_layers.items()
            )
        )
    panel_f_shape = tuple(int(value) for value in panel_f_dapi_labels.shape)
    mismatched_panel_f_shapes = {
        name: tuple(int(value) for value in np.asarray(layer).shape)
        for name, layer in panel_f_layers.items()
        if np.asarray(layer).shape != panel_f_dapi_labels.shape
    }
    raw_panel_f_shapes = {
        "raw DAPI": marker_signal(iba1_dapi).shape,
        "raw Iba1": marker_signal(iba1).shape,
        "raw GFAP": marker_signal(panel_f_raw_gfap).shape,
    }
    mismatched_panel_f_shapes.update({
        name: tuple(int(value) for value in shape)
        for name, shape in raw_panel_f_shapes.items()
        if tuple(shape) != panel_f_shape
    })
    if mismatched_panel_f_shapes:
        raise RuntimeError(
            f"Panel F layers must share the exact Panel A field {panel_f_shape}; "
            f"mismatches={mismatched_panel_f_shapes}"
        )
    for name, labels_array in (
        ("DAPI labels", panel_f_dapi_labels),
        ("Iba1 cells", panel_f_iba1_cells),
    ):
        if not np.issubdtype(labels_array.dtype, np.integer):
            raise RuntimeError(f"Panel F {name} must retain integer instance labels")
        if np.any(labels_array < 0) or not np.any(labels_array > 0):
            raise RuntimeError(f"Panel F {name} must contain nonnegative, nonempty labels")
    region_values = set(int(value) for value in np.unique(panel_f_regions))
    if region_values != {0, *REGION_CODES.values()}:
        raise RuntimeError(
            "Panel F final region mask must contain exactly background plus ARC/ME/VMN; "
            f"found={sorted(region_values)}"
        )
    if not np.array_equal(panel_f_regions, iba1_labels):
        raise RuntimeError(
            "Panel F final ARC/ME/VMN mask is not pixel-identical to the accepted Panel A mask"
        )
    if not np.any(panel_f_gfap_positive):
        raise RuntimeError("Panel F final GFAP-positive mask is empty")
    panel_f_tissue = channel_cartoon_tissue_envelope(
        panel_f_dapi_labels > 0, panel_f_regions,
    )
    # This also gates the audited two-ARC/one-ME/two-VMN topology before output.
    panel_f_region_label_points = channel_cartoon_region_label_points(panel_f_regions)
    panel_f_instance_counts = {
        "DAPI_instances": int(np.count_nonzero(np.unique(panel_f_dapi_labels) > 0)),
        "Iba1_QC_accepted_cells": int(np.count_nonzero(np.unique(panel_f_iba1_cells) > 0)),
        "GFAP_positive_pixels": int(np.count_nonzero(panel_f_gfap_positive)),
    }
    panel_f_positive_pixels = {
        "DAPI": int(np.count_nonzero(panel_f_dapi_labels > 0)),
        "Iba1": int(np.count_nonzero(panel_f_iba1_cells > 0)),
        "GFAP": int(np.count_nonzero(panel_f_gfap_positive)),
    }
    panel_f_binary_hashes = {
        channel: hashlib.sha256(
            np.ascontiguousarray(mask, dtype=np.uint8).tobytes()
        ).hexdigest()
        for channel, mask in (
            ("DAPI", panel_f_dapi_labels > 0),
            ("Iba1", panel_f_iba1_cells > 0),
            ("GFAP", panel_f_gfap_positive),
        )
    }
    iba1_bounds = (0, 0, iba1_dapi.shape[1], min(720, iba1_dapi.shape[0]))
    gfap_bounds = (0, 0, gfap_dapi.shape[1], min(720, gfap_dapi.shape[0]))
    panel_f_bounds = (0, 0, panel_f_shape[1], panel_f_shape[0])
    if panel_f_bounds != iba1_bounds:
        raise RuntimeError(
            f"Panel F bounds {panel_f_bounds} do not match Panel A bounds {iba1_bounds}"
        )
    panel_f_audit = {
        "bounds_px": ",".join(str(value) for value in panel_f_bounds),
        "shape_yx": panel_f_shape,
        "instance_counts": panel_f_instance_counts,
        "positive_pixels": panel_f_positive_pixels,
        "binary_sha256": panel_f_binary_hashes,
        "region_label_points": panel_f_region_label_points,
    }
    micro_x, micro_y, micro_uid = select_microglia_inset(
        cells, args.iba1_animal, args.iba1_section, iba1_bounds,
        preferred_uid=args.iba1_cell_uid,
    )
    _, micro_bbox = crop_square(
        marker_signal(iba1), (micro_x, micro_y), PANEL_A_INSET_SIZE,
    )
    inset_state_rows = cells.loc[cells["cell_uid"].astype(str).eq(str(micro_uid)), "microglia_state"]
    inset_state = str(inset_state_rows.iloc[0]) if not inset_state_rows.empty else ""
    inset_label = f"{inset_state} cell" if inset_state in STATE_COLORS else "Microglial cell"
    gfap_center = select_gfap_inset(gfap, positive, labels, gfap_bounds)
    examples = select_examples(cells, args.analysis_dir)
    iba1_drawn_regions = tuple(
        region for region in ("ARC", "ME", "VMN")
        if iba1_provenance["region_status"].get(region) == "drawn"
    )
    section_state_counts = section_state_composition(cells, args.iba1_animal, args.iba1_section, iba1_drawn_regions)
    if any(
        sum(section_state_counts.get(region, {}).values()) <= 0
        for region in REGIONS
    ):
        raise RuntimeError(
            "Panel D representative-section state plot has no accepted "
            "microglia in one or more regions"
        )
    first_example = next((entry["row"] for entry in examples.values() if entry["row"] is not None), None)
    input_patch = mpimg.imread(resolve_patch(first_example, args.analysis_dir)) if first_example is not None else np.zeros((96, 96, 3), dtype=float)
    metadata_path = args.analysis_dir / "microglia_classifier_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.is_file() else {}
    history_path = args.analysis_dir / "microglia_cnn_training_history.csv"
    history = pd.read_csv(history_path) if history_path.is_file() else pd.DataFrame()

    # Only after every analysis, receipt, source-image, and representative-mask
    # gate has passed may publication output be cleaned or created.
    if bool(args.force):
        try:
            args.output_dir = _force_clean_generated_dir(
                output_arg,
                forbidden_roots=(
                    args.analysis_dir,
                    args.input_root,
                    PROJECT_ROOT,
                    PROJECT_ROOT / "scripts",
                    PROJECT_ROOT / "Iba1_GFAP_final_10x",
                    PROJECT_ROOT / "20x",
                ),
                preserved_roots=(manual_root,),
                preserved_names=("README.txt",),
            )
        except ValueError as exc:
            raise SystemExit(f"ERROR: {exc}") from exc
    else:
        args.output_dir = output_arg.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    panels_dir = args.output_dir / "panels"
    legends_dir = args.output_dir / "legends"
    source_data_dir = args.output_dir / "source_data"
    provenance_dir = args.output_dir / "provenance"
    for directory in (panels_dir, legends_dir, source_data_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=True)
    mask_audit_path = (
        source_data_dir
        / "Figure_GFAP_Iba1_microglia_ARC_ME_all_HIL_final_mask_audit.csv"
    )
    recomputed_values_path = (
        source_data_dir
        / "Figure_GFAP_Iba1_microglia_ARC_ME_animal_values_recomputed.csv"
    )
    # Publication audit retains receipt hashes/timestamps but omits reviewer PII.
    coverage.drop(columns=["reviewer"], errors="ignore").to_csv(
        mask_audit_path, index=False,
    )
    recomputed_animal_data.to_csv(recomputed_values_path, index=False)

    # Keep the A/B microscopy pairs large, compact the explanatory row, and
    # preserve the established physical height of the quantitative panels.
    fig = plt.figure(figsize=(16.0, 17.0755), facecolor="white")
    panel_axes: dict[str, list[plt.Axes]] = {}

    def record_panel(letter: str, axes_before: set[int]) -> None:
        new_axes = [axis for axis in fig.axes if id(axis) not in axes_before]
        panel_axes.setdefault(letter, []).extend(new_axes)

    outer = fig.add_gridspec(
        3, 12,
        left=0.035, right=0.99, bottom=0.035, top=0.99,
        wspace=0.20, hspace=0.07,
        height_ratios=(2.90, 3.85, 8.83),
    )
    microscopy_row = outer[0, 0:12].subgridspec(
        1, 2, width_ratios=(1.0, 1.0), wspace=0.08,
    )
    panel_a_microscopy = microscopy_row[0, 0].subgridspec(
        1, 2, width_ratios=(1.0, 1.0), wspace=0.025,
    )
    axes_before = {id(axis) for axis in fig.axes}
    iba1_native_ax = draw_microscopy(
        fig, panel_a_microscopy[0, 0], "A", "Iba1 / DAPI", iba1_dapi, iba1,
        (1.0, 0.14, 0.20), iba1_bounds, (micro_x, micro_y), inset_label,
        args.um_per_px, inset_size=PANEL_A_INSET_SIZE,
        source_dapi_sha256=str(iba1_provenance["source_dapi_sha256"]),
        inset_outline_color="#ffe84a", inset_width="28%", inset_height="36%",
        source_box_label="selected cell",
    )
    iba1_accepted_ax = draw_anatomy_qc(
        fig, panel_a_microscopy[0, 1], iba1_dapi, iba1_labels,
        iba1_bounds, iba1_provenance, args.um_per_px,
    )
    validate_aligned_representative_pair(iba1_native_ax, iba1_accepted_ax)
    record_panel("A", axes_before)

    panel_b_microscopy = microscopy_row[0, 1].subgridspec(
        1, 2, width_ratios=(1.0, 1.0), wspace=0.025,
    )
    axes_before = {id(axis) for axis in fig.axes}
    gfap_native_ax = draw_microscopy(
        fig, panel_b_microscopy[0, 0], "B", "GFAP / DAPI", gfap_dapi, gfap,
        (0.08, 1.0, 0.32), gfap_bounds, gfap_center, "GFAP processes",
        args.um_per_px,
        source_dapi_sha256=str(gfap_provenance["source_dapi_sha256"]),
    )
    gfap_accepted_ax = draw_anatomy_qc(
        fig, panel_b_microscopy[0, 1], gfap_dapi, labels,
        gfap_bounds, gfap_provenance, args.um_per_px,
    )
    validate_aligned_representative_pair(gfap_native_ax, gfap_accepted_ax)
    record_panel("B", axes_before)

    explanation_row = outer[1, 0:12].subgridspec(
        1, 2, width_ratios=(8.5, 6.5), wspace=0.08,
    )

    axes_before = {id(axis) for axis in fig.axes}
    c_band = explanation_row[0, 0].subgridspec(
        2, 1, height_ratios=(3.30, 0.55), hspace=0.0,
    )
    micro_method = c_band[0, 0].subgridspec(
        2, 1, height_ratios=(1.65, 1.25), hspace=0.28,
    )
    cnn_ax = fig.add_subplot(micro_method[0, 0])
    draw_cnn(
        cnn_ax, input_patch, metadata, history,
        um_per_px=args.um_per_px,
        inside_label=True, panel_letter="C",
    )
    examples_ax = fig.add_subplot(micro_method[1, 0])
    draw_examples(examples_ax, examples, args.um_per_px, args.analysis_dir)
    record_panel("C", axes_before)

    explanatory_right = explanation_row[0, 1].subgridspec(
        3, 1, height_ratios=(1.35, 0.30, 2.20), hspace=0.0,
    )
    d_grid = explanatory_right[0, 0].subgridspec(
        2, 1, height_ratios=(0.30, 1.05), hspace=0.0,
    )
    axes_before = {id(axis) for axis in fig.axes}
    draw_section_state_qc(
        fig, d_grid[1, 0], section_state_counts,
        header_spec=d_grid[0, 0], panel_letter="D",
    )
    record_panel("D", axes_before)

    axes_before = {id(axis) for axis in fig.axes}
    gfap_ax = fig.add_subplot(explanatory_right[2, 0])
    draw_gfap_workflow(
        gfap_ax, gfap_dapi, gfap, positive, labels, gfap_center,
        um_per_px=args.um_per_px,
        panel_letter="E",
    )
    record_panel("E", axes_before)

    # Preserve the established quantitative width while using the former left
    # whitespace for the exact-coordinate, vertical channel-mask Panel F.
    lower_band = outer[2, 0:12].subgridspec(
        1, 2, width_ratios=(3.0, 9.0), wspace=0.0,
    )
    panel_f_slot = lower_band[0, 0].subgridspec(
        1, 2, width_ratios=(2.4, 0.6), wspace=0.0,
    )
    panel_f_grid = panel_f_slot[0, 0].subgridspec(
        4, 1, height_ratios=(0.15, 1.0, 1.0, 1.0), hspace=0.08,
    )
    axes_before = {id(axis) for axis in fig.axes}
    panel_f_header_ax = fig.add_subplot(panel_f_grid[0, 0])
    panel_f_header_ax.set_axis_off()
    panel_f_header_ax.text(
        -0.018, 0.50, "F", transform=panel_f_header_ax.transAxes,
        fontsize=22, fontweight="bold", ha="right", va="center",
        color="black", clip_on=False,
    )
    panel_f_header_ax.text(
        0.01, 0.50, "Channel analysis masks",
        transform=panel_f_header_ax.transAxes, fontsize=8.6,
        fontweight="bold", ha="left", va="center", color="black",
    )
    for panel_f_row, (marker, marker_mask) in enumerate((
        ("DAPI", panel_f_dapi_labels > 0),
        ("Iba1", panel_f_iba1_cells > 0),
        ("GFAP", panel_f_gfap_positive),
    ), start=1):
        channel_ax = fig.add_subplot(panel_f_grid[panel_f_row, 0])
        draw_channel_analysis_cartoon(
            channel_ax, marker, marker_mask,
            panel_f_tissue, panel_f_regions, args.um_per_px,
        )
    record_panel("F", axes_before)

    graph_grid = lower_band[0, 1].subgridspec(
        len(PANEL_ENDPOINTS), 1, hspace=0.26,
    )
    all_stats: list[dict] = []
    for row, endpoint in enumerate(PANEL_ENDPOINTS):
        letter = endpoint[0]
        spec = graph_grid[row, 0]
        axes_before = {id(axis) for axis in fig.axes}
        all_stats.extend(
            draw_endpoint_pair(
                fig, spec, endpoint, data,
                seed=int(args.seed) + row * 100,
                show_condition_labels=True,
            )
        )
        record_panel(letter, axes_before)

    expected_panel_letters = set("ABCDEFGHI")
    if set(panel_axes) != expected_panel_letters:
        raise RuntimeError(
            "Figure 5 panel registry must contain exactly A-I; "
            f"found={sorted(panel_axes)}"
        )
    empty_panels = sorted(letter for letter, axes in panel_axes.items() if not axes)
    if empty_panels:
        raise RuntimeError(f"Figure 5 panel registry contains empty panels: {empty_panels}")

    stem = args.output_dir / "Figure_GFAP_Iba1_microglia_ARC_ME"
    for extension in formats:
        # The PNG used to be capped at 300 while the PDF got the full --dpi, so the
        # two deliverables disagreed. Both now honour --dpi.
        save_dpi = int(args.dpi)
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=save_dpi, facecolor="white", bbox_inches="tight")
    export_panel_crops(fig, panel_axes, panels_dir)
    translation_receipt = translate_figure_texts_to_spanish(
        fig, extra=SPANISH_PANEL_TEXT
    )
    # The longer Spanish y-axis label reaches the header band. Move only the
    # Spanish Panel D letter left, preserving the wording and all data geometry.
    spanish_panel_d_letters = [
        artist
        for axis in panel_axes["D"]
        for artist in axis.texts
        if artist.get_text() == "D" and math.isclose(artist.get_fontsize(), 22.0)
    ]
    if len(spanish_panel_d_letters) != 1:
        raise RuntimeError("Could not identify the unique Spanish Panel D letter")
    spanish_panel_d_letters[0].set_x(-0.11)
    export_panel_crops(fig, panel_axes, panels_dir, suffix="_spanish")
    translation_receipt_path = (
        provenance_dir
        / "Figure_GFAP_Iba1_microglia_ARC_ME_spanish_translation_receipt.json"
    )
    translation_receipt_path.write_text(
        json.dumps(
            {
                "schema": "fig5_bilingual_isolated_panels_v1",
                "master_language": "en",
                "isolated_panel_locales": ["en", "es"],
                "panel_letters": list("ABCDEFGHI"),
                "raster_dpi": int(args.dpi),
                "text_replacements": translation_receipt,
            },
            indent=2,
            sort_keys=True,
            ensure_ascii=False,
        ) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)
    stats = finalize_animal_statistics(pd.DataFrame(all_stats))
    denominator_audit = activated_denominator_audit(data)
    cage_values = multipanel_sensitivity_cage_values(data, denominator_audit)
    sensitivity_stats = multipanel_sensitivity_statistics(cage_values)
    stats.to_csv(source_data_dir / "Figure_GFAP_Iba1_microglia_ARC_ME_statistics.csv", index=False)
    sensitivity_paths = write_multipanel_sensitivity_outputs(
        source_data_dir, legends_dir, provenance_dir, values_path, cage_values, denominator_audit, sensitivity_stats,
    )
    write_legends(
        legends_dir, values_path, metadata_path,
        args.iba1_animal, args.iba1_section, args.animal, args.section,
        micro_uid, micro_bbox, gfap_center,
        examples, stats, sensitivity_stats, metadata, data,
        iba1_provenance, gfap_provenance, panel_f_audit,
    )

    source_rows: list[tuple[str, str]] = []
    for role, path in (
        ("animal_values", values_path),
        ("animal_values_independently_recomputed", recomputed_values_path),
        ("per_cell_values", cells_path),
        ("all_HIL_coverage", args.analysis_dir / "human_review_coverage.csv"),
        ("all_HIL_section_summary", args.analysis_dir / "per_section_region_summary.csv"),
        ("all_HIL_final_mask_audit", mask_audit_path),
        ("panel_A_DAPI", iba1_dapi_path),
        ("panel_A_Iba1", iba1_path),
        ("panel_A_native_accepted_ARC_ME_labels", iba1_provenance["mask_path"]),
        ("panel_A_native_review_receipt", iba1_provenance["review_path"]),
        ("panel_A_final_analysis_ARC_ME_labels", iba1_final_mask_path),
        ("panel_A_final_analysis_geometry", iba1_final_mask_path.parent / "region_geometry.json"),
        ("panel_B_DAPI", gfap_dapi_path),
        ("panel_B_GFAP", gfap_path),
        ("panel_E_GFAP_positive_mask", positive_path),
        ("panel_B_native_accepted_ARC_ME_labels", gfap_provenance["mask_path"]),
        ("panel_B_native_review_receipt", gfap_provenance["review_path"]),
        ("panel_B_final_analysis_ARC_ME_labels", gfap_final_mask_path),
        ("panel_B_final_analysis_geometry", gfap_final_mask_path.parent / "region_geometry.json"),
        ("panel_F_raw_DAPI", iba1_dapi_path),
        ("panel_F_raw_Iba1", iba1_path),
        ("panel_F_raw_GFAP", panel_f_raw_gfap_path),
        ("panel_F_DAPI_labels", panel_f_dapi_labels_path),
        ("panel_F_Iba1_QC_accepted_cells", panel_f_iba1_cells_path),
        ("panel_F_GFAP_positive_mask", panel_f_gfap_positive_path),
        ("panel_F_native_accepted_ARC_ME_VMN_labels", iba1_provenance["mask_path"]),
        ("panel_F_final_analysis_ARC_ME_VMN_labels", iba1_final_mask_path),
        ("classifier_metadata", metadata_path),
        ("panel_C_classifier_metadata", metadata_path),
        ("spanish_translation_receipt", translation_receipt_path),
        ("multipanel_sensitivity_statistics", sensitivity_paths["sensitivity_statistics"]),
        ("multipanel_sensitivity_cage_values", sensitivity_paths["sensitivity_cage_values"]),
        ("multipanel_denominator_audit", sensitivity_paths["denominator_audit"]),
        ("multipanel_sensitivity_methods", sensitivity_paths["sensitivity_methods"]),
        ("multipanel_sensitivity_provenance", sensitivity_paths["sensitivity_provenance"]),
    ):
        source_rows.append((
            role,
            _portable_provenance_path(
                path, args.analysis_dir, args.output_dir,
            ),
        ))
    source_rows.extend((
        ("all_HIL_required_sections", str(EXPECTED_HUMAN_REVIEW_SECTIONS)),
        ("all_HIL_verified_coverage_rows", str(len(coverage))),
        ("all_HIL_verified_section_region_rows", str(len(section_summary))),
        ("all_HIL_final_masks_pixel_identical", "61/61"),
        ("animal_values_equal_independent_reaggregation", "true"),
        ("all_HIL_accepted_source_counts", json.dumps(coverage["accepted_source"].value_counts().sort_index().to_dict(), sort_keys=True)),
        ("panel_A_representative_section", f"{args.iba1_animal}/{args.iba1_section}"),
        ("panel_A_source_bounds_px", ",".join(str(value) for value in iba1_bounds)),
        ("panel_A_native_and_accepted_bounds_identical", "true"),
        ("panel_A_HIL_layout", "accepted ROI mask immediately right of corresponding Iba1/DAPI field"),
        ("panel_A_bold_contour_widths_pt", "colored=3.5;dark_halo=6.4"),
        ("panel_A_accepted_mask_sha256", str(iba1_provenance["mask_sha256"])),
        ("panel_A_accepted_receipt_sha256", str(iba1_provenance["review_sha256"])),
        ("panel_A_source_dapi_sha256", str(iba1_provenance["source_dapi_sha256"])),
        ("panel_A_accepted_source", str(iba1_provenance["accepted_source"])),
        ("panel_A_accepted_at", str(iba1_provenance["accepted_at"])),
        ("panel_A_final_mask_pixel_identical_to_accepted", "true"),
        ("panel_A_final_geometry_region_source", str(iba1_geometry.get("region_source", ""))),
        ("panel_A_inset_bbox_px", ",".join(str(value) for value in micro_bbox)),
        ("panel_B_representative_section", f"{args.animal}/{args.section}"),
        ("panel_B_source_bounds_px", ",".join(str(value) for value in gfap_bounds)),
        ("panel_B_native_and_accepted_bounds_identical", "true"),
        ("panel_B_HIL_layout", "accepted ROI mask immediately right of corresponding GFAP/DAPI field"),
        ("panel_B_bold_contour_widths_pt", "colored=3.5;dark_halo=6.4"),
        ("panel_B_accepted_mask_sha256", str(gfap_provenance["mask_sha256"])),
        ("panel_B_accepted_receipt_sha256", str(gfap_provenance["review_sha256"])),
        ("panel_B_source_dapi_sha256", str(gfap_provenance["source_dapi_sha256"])),
        ("panel_B_accepted_source", str(gfap_provenance["accepted_source"])),
        ("panel_B_accepted_at", str(gfap_provenance["accepted_at"])),
        ("panel_B_final_mask_pixel_identical_to_accepted", "true"),
        ("panel_B_final_geometry_region_source", str(gfap_geometry.get("region_source", ""))),
        ("panel_D_section_state_QC_representative_only", json.dumps(section_state_counts, sort_keys=True)),
        ("panel_E_GFAP_workflow", "GFAP quantification workflow using the representative raw GFAP/DAPI field and final GFAP-positive mask"),
        ("panel_F_representative_section", f"{args.iba1_animal}/{args.iba1_section}"),
        ("panel_F_source_bounds_px", panel_f_audit["bounds_px"]),
        ("panel_F_same_source_field_as_panel_A", "true"),
        ("panel_F_all_layers_share_exact_coordinates", "true"),
        ("panel_F_final_region_mask_pixel_identical_to_panel_A_accepted_mask", "true"),
        ("panel_F_channel_instance_counts", json.dumps(panel_f_instance_counts, sort_keys=True)),
        ("panel_F_channel_positive_pixels", json.dumps(panel_f_positive_pixels, sort_keys=True)),
        ("panel_F_channel_binary_sha256", json.dumps(panel_f_binary_hashes, sort_keys=True)),
        ("panel_F_region_label_points", json.dumps(panel_f_region_label_points)),
        ("panel_F_tissue_envelope", "shared DAPI-density visualization support union final ARC/ME/VMN mask; not an anatomical ROI or quantitative denominator"),
        ("panel_F_channel_colors", json.dumps(CHANNEL_CARTOON_COLORS, sort_keys=True)),
        ("panel_F_scalebar_um", "100"),
        ("panel_G_I_animal_point_jitter_seed", str(int(args.seed))),
        ("primary_statistics", "animal one-way ANOVA global; raw exact two-sided MWU W-S/W-A/S-A; no multiplicity adjustment"),
        ("cage_mean_sensitivity", "unweighted animal means per cage; one-way ANOVA and raw exact two-sided MWU; no multiplicity adjustment"),
        ("activated_count_sensitivity", "cage-summed activated/total microglia; raw exact cage-label binomial deviance global and W-S/W-A/S-A; no multiplicity adjustment"),
        ("sensitivity_filtering", "none; positive low denominators retained; zero-denominator animal count pairs retained in cage sums and audited as ratio-non-estimable"),
        ("microglia_inset_cell_uid", str(micro_uid)),
        ("microglia_inset_predicted_state", inset_state or "NA"),
        ("microglia_inset_center_px", f"{micro_x:.3f},{micro_y:.3f}"),
        ("gfap_inset_center_px", f"{gfap_center[0]:.3f},{gfap_center[1]:.3f}"),
    ))
    for state, entry in examples.items():
        row = entry["row"]
        if row is None:
            value = f"NONE|selected=0/{entry['total_predictions']}|cell_qc_candidates={entry['cell_qc_candidates']}|patches_checked={entry['patches_checked']}|reason={entry['reason']}"
        else:
            patch_source = _portable_provenance_path(
                resolve_patch(row, args.analysis_dir),
                args.analysis_dir,
                args.output_dir,
            )
            value = f"{row['cell_uid']}|{patch_source}|selected=1/{entry['total_predictions']}|cell_qc_candidates={entry['cell_qc_candidates']}|patches_checked={entry['patches_checked']}|reason={entry['reason']}"
        source_rows.append((f"panel_C_example_{state}", value))

    manifest_path = provenance_dir / "Figure_GFAP_Iba1_microglia_ARC_ME_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("role", "path_or_value"))
        writer.writerows(source_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
