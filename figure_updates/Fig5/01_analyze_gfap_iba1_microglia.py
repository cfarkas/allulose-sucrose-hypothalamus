#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reproducible Keyence 10x DAPI/GFAP/Iba1 ARC/ME and microglial morphology pipeline.

Numbered copy created 2026-08-12 from
``../../keyence10x_gfap_iba1_microglia_arc_me.py``
(source SHA256 e466b2baa129bb01dcad16220b5164510bf24eabd1082adad06a5e7e19dbf830).

CH1=DAPI, CH3=GFAP, CH4=Iba1. Each XY folder is one full bilateral section.
The validated v17 projected-third-ventricle-wall ARC/ME anatomy and the exact v18/v19
ME-left/ARC-right figure style are retained. Existing segmentations are preserved.

When ``--region-mode arc-me-dl`` is selected, the tiny U-Net is trained from DAPI-only
anatomical features and a projected-wall prior. Its ME prediction remains constrained to the closed ventricular-floor prior or to
a small explicitly configured boundary band. GFAP, Iba1 and condition labels never
enter ARC/ME placement. Microglial state classification remains a separate morphology +
appearance analysis and is summarized only after region assignment.
"""

from __future__ import annotations

# stdlib only here; scientific imports happen after optional conda bootstrap
import argparse
import csv
import datetime as dt
import hashlib
import logging
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import warnings
from dataclasses import dataclass, asdict
from itertools import combinations

warnings.filterwarnings("ignore", category=FutureWarning, module="skimage")
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


# -----------------------------------------------------------------------------
# Defaults for this experiment
# -----------------------------------------------------------------------------

HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")
PROJECT_ROOT = PAPER_ROOT
FIG5_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig5"



def _first_existing(*candidates: Path) -> Path:
    """First candidate that exists, else the first, so defaults name real paths.

    The GFAP/Iba1 acquisition ships as Fig5/raw_data. analyses/Fig5/raw/input is
    the optional private working copy and wins only when it is actually present,
    so the documented command also works in a freshly unpacked tree.
    """
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


FIG5_FIGURE_ROOT = PAPER_ROOT / "Fig5"
DEFAULT_RAW_ROOT = _first_existing(
    FIG5_ANALYSIS_ROOT / "raw" / "input", FIG5_FIGURE_ROOT / "raw_data",
)
DEFAULT_SANITIZED_ROOT = DEFAULT_RAW_ROOT
DEFAULT_UM_PER_PX = 0.755
FIG5_HIL_REVIEW_ROOT = FIG5_FIGURE_ROOT / "hil_review"

DEFAULT_REVIEWED_DAPI_ROOT = FIG5_HIL_REVIEW_ROOT / "dapi_final_20260830_v1"
DEFAULT_REVIEWED_IBA1_ROOT = FIG5_HIL_REVIEW_ROOT / "iba1_final_20260830_v1"

def _latest_accepted_hil_root():
    """Newest shipped review under Fig5/hil_review that actually holds receipts."""
    if not FIG5_HIL_REVIEW_ROOT.is_dir():
        return None
    roots = sorted(
        path for path in FIG5_HIL_REVIEW_ROOT.glob("human_final_*")
        if path.is_dir() and any(path.glob("*/*/*.json"))
    )
    return roots[-1] if roots else None


def _first_reviewed(*candidates):
    """First candidate holding accepted receipts, else the first candidate.

    Existence is not enough here: an earlier run leaves an empty
    analyses/Fig5/hil/accepted_annotations behind, and choosing that empty
    directory reads as a review of zero sections rather than as the completed
    61/61 review this tree ships.
    """
    for candidate in candidates:
        if candidate is not None and candidate.is_dir() and any(candidate.glob("*/*/*.json")):
            return candidate
    return next((candidate for candidate in candidates if candidate is not None), None)


# A fresh review writes analyses/Fig5/hil/accepted_annotations; the accepted
# 61/61 review ships under Fig5/hil_review, and is used when no working copy
# holds receipts, so the documented command works in a freshly unpacked tree.
DEFAULT_MANUAL_REGION_ROOT = _first_reviewed(
    FIG5_ANALYSIS_ROOT / "hil" / "accepted_annotations", _latest_accepted_hil_root(),
)
# Generated analysis output. It is a working directory, never shipped input, so
# the renderer's default --analysis-dir resolves to it and the two documented
# commands chain with no arguments on a freshly unpacked tree. It must never
# default inside the raw acquisition, which is read-only.
DEFAULT_ANALYSIS_OUTPUT = FIG5_ANALYSIS_ROOT / "results" / "human_final_run"

DEFAULT_MODEL_POMC = FIG5_ANALYSIS_ROOT / "raw" / "models" / "unused_pomc_model"
DEFAULT_MODEL_CFOS = FIG5_ANALYSIS_ROOT / "raw" / "models" / "unused_cfos_model"
DEFAULT_MODEL_DAPI = FIG5_ANALYSIS_ROOT / "raw" / "models" / "unused_dapi_model"

CHANNEL_TO_MARKER = {
    "CH1": "DAPI",
    "CH2": "NPY",
    "CH3": "POMC",
    "CH4": "cFOS",
}
MARKER_TO_CH = {v: k for k, v in CHANNEL_TO_MARKER.items()}
MARKERS = ["DAPI", "NPY", "POMC", "cFOS"]

DEFAULT_SAMPLESHEET_CSV = """sample,cond,genotype,sex,cage
FR-721,Allulose,WT,M,JAULA1
FR-722,Allulose,Tg,M,JAULA1
FR7-23,Sucrose,WT,M,JAULA2
FR-724,Sucrose,WT,M,JAULA2
FR8-1,Allulose,WT,F,JAULA3
FR6-1,Allulose,Tg,F,JAULA3
FR8-2,Sucrose,WT,F,JAULA4
FR6-2,Sucrose,Tg,F,JAULA4
FR8-3,Allulose,WT,F,JAULA5
FR6-3,Allulose,WT,F,JAULA5
FR6-4,Water,Tg,F,JAULA6
FR7-1,Water,WT,F,JAULA6
FR6-5,Water,WT,F,JAULA7
FR7-2,Water,WT,F,JAULA7
FR7-3,Sucrose,WT,F,JAULA8
FR7-4,Sucrose,WT,F,JAULA8
Negative_ctr,Control,NA,NA,NA
"""

# regexes
RE_XY_DIR = re.compile(r"^XY(?P<num>\d+)$", re.IGNORECASE)
RE_S_PLANE_DIR = re.compile(r"^S(?P<num>\d{1,3})_(?P<side>left|right)(?:_T(?P<tile>\d+))?$", re.IGNORECASE)
RE_CH_TIF = re.compile(
    r"^(?P<prefix>.*?)(?:_(?P<tile>\d{3,6}))?_CH(?P<ch>[1-4])\.(?P<ext>tif|tiff)$",
    re.IGNORECASE,
)
RE_CH_SEG = re.compile(
    r"^(?P<prefix>.*?)(?:_(?P<tile>\d{3,6}))?_CH(?P<ch>[1-4])_seg\.npy$",
    re.IGNORECASE,
)
RE_MARKER_TIF = re.compile(r"_(?P<marker>DAPI|NPY|POMC|cFOS)\.(tif|tiff)$", re.IGNORECASE)
RE_MARKER_SEG = re.compile(r"_(?P<marker>DAPI|NPY|POMC|cFOS)_seg\.npy$", re.IGNORECASE)


# Scientific module globals populated by import_science_stack().
np = None
pd = None
plt = None
matplotlib = None
regionprops = None
label = None
binary_dilation = None
remove_small_objects = None
remove_small_holes = None
binary_opening = None
binary_closing = None
disk = None
expand_labels = None
convex_hull_image = None
tifffile = None
f_oneway = None
mannwhitneyu = None
wilcoxon = None
ttest_rel = None


# -----------------------------------------------------------------------------
# CLI / bootstrap
# -----------------------------------------------------------------------------

def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Full Keyence 20x manual-segmentation-aware ARC/ME-only pipeline using validated v17 projected-wall anatomy and final compact plots plus definitive animal-level statistics."
    )

    # Paths
    ap.add_argument("--raw-input", "--raw-root", dest="raw_root", default=str(DEFAULT_RAW_ROOT),
                    help="Raw 20x folder. Used only if sanitized root is missing or --force-migrate is used.")
    ap.add_argument("--input", "--sanitized-root", "--root", dest="sanitized_root", default=str(DEFAULT_SANITIZED_ROOT),
                    help="Input/working folder with sanitized images and manual *_seg.npy masks. Existing manual segmentations here are kept.")
    ap.add_argument("--output", "--outdir", dest="outdir", default=None,
                    help=("Output directory for CSVs/plots/logs. Default: "
                          "analyses/Fig5/results/human_final_run, which is also the "
                          "renderer's default --analysis-dir."))
    ap.add_argument("--samplesheet", default=None,
                    help="Optional CSV with columns sample,cond,genotype,sex,cage. If omitted, embedded default is used.")
    ap.add_argument("--include-samples", nargs="*", default=None,
                    help="Optional sample names to process, e.g. --include-samples FR6-1 FR6-2. Matching ignores hyphens/parentheses.")
    ap.add_argument("--exclude-samples", nargs="*", default=None,
                    help="Optional sample names to exclude.")
    ap.add_argument("--limit-planes", type=int, default=0,
                    help="Debug option: process only the first N discovered planes after sample filtering. 0 = no limit.")

    # Stage controls
    ap.add_argument("--dry-run", action="store_true", help="Print what would be done; do not modify files or run Cellpose.")
    ap.add_argument("--full-run", action="store_true", help="Run the normal complete pipeline using defaults. This flag is optional; it is provided for a minimal, readable command.")
    ap.add_argument("--migrate-only", action="store_true", help="Only create/update sanitized folder, then stop.")
    ap.add_argument("--segment-only", action="store_true", help="Only check/run missing segmentations, then stop.")
    ap.add_argument("--quantify-only", action="store_true", help="Only quantify/plot existing segmentations.")
    ap.add_argument("--force-migrate", action="store_true", help="Run migration even if sanitized root already exists. Existing files are not overwritten unless --overwrite-sanitized-files.")
    ap.add_argument("--overwrite-sanitized-files", action="store_true", help="During migration, overwrite existing sanitized images/masks. Usually keep this OFF.")
    ap.add_argument("--migration-tile-mode", choices=["all", "first", "last", "collapse"], default="all",
                    help="If raw XY folders contain Image_XY##_00001...00009_CH*.tif, preserve all tile IDs, or keep only first/last/collapse.")

    seg_group = ap.add_mutually_exclusive_group()
    seg_group.add_argument("--segment-missing", dest="segment_missing", action="store_true", default=False,
                           help="Maintenance mode: run Cellpose only for missing raw-local *_seg.npy files.")
    seg_group.add_argument("--no-segment-missing", dest="segment_missing", action="store_false",
                           help="Do not run Cellpose. Quantify only available manual segmentations.")
    ap.add_argument("--overwrite-seg", action="store_true", help="Overwrite existing segmentation masks. Default is to preserve manual masks.")

    # Conda / environment
    ap.add_argument("--env-name", default="paper_apotome_repro", help="Unified Paper conda environment.")
    ap.add_argument("--allow-existing-env-update", action="store_true", default=False,
                    help="Safety override: allow installing packages into base/cellpose_env if --env-name points there. Not recommended.")
    ap.add_argument("--auto-conda", dest="auto_conda", action="store_true", default=True,
                    help="If not in --env-name, create/use conda env and re-run automatically.")
    ap.add_argument("--no-auto-conda", dest="auto_conda", action="store_false", help="Disable automatic conda bootstrap/re-exec.")
    ap.add_argument("--setup-only", action="store_true", help="Create/install conda environment, then stop.")
    ap.add_argument("--install-cellpose", action="store_true", default=False,
                    help="When creating a new env, also pip install cellpose. Enabled automatically if --segment-missing needs missing masks.")

    # Cellpose models and parameters
    ap.add_argument("--model-pomc", default=str(DEFAULT_MODEL_POMC), help="POMC Cellpose model path.")
    ap.add_argument("--model-cfos", default=str(DEFAULT_MODEL_CFOS), help="c-FOS Cellpose model path.")
    ap.add_argument("--model-dapi", default=str(DEFAULT_MODEL_DAPI), help="DAPI Cellpose model path.")
    ap.add_argument("--model-npy", default=None, help="Optional NPY Cellpose model path. If omitted, uses DAPI model as fallback.")
    ap.add_argument("--gpu", "--use-gpu", action="store_true", help="Use GPU for missing Cellpose segmentations.")
    ap.add_argument("--gpu-device", default="0", help="GPU device id for Cellpose.")
    ap.add_argument("--diameter-dapi", type=float, default=None)
    ap.add_argument("--diameter-pomc", type=float, default=None)
    ap.add_argument("--diameter-npy", type=float, default=None)
    ap.add_argument("--diameter-cfos", type=float, default=None)
    ap.add_argument("--flow-threshold", type=float, default=None)
    ap.add_argument("--cellprob-threshold", type=float, default=None)
    ap.add_argument("--min-size-seg", type=int, default=None)

    # Quantification parameters
    ap.add_argument("--um-per-px", "--pixel-size-um", type=float, default=DEFAULT_UM_PER_PX,
                    help="Microns per pixel.")
    ap.add_argument("--cfos-association-mode", choices=["strict", "blur-aware", "centroid-fallback"], default="strict",
                    help="How c-FOS objects are assigned to DAPI nuclei. strict = original DAPI-mask overlap only; blur-aware = also allows expanded DAPI labels; centroid-fallback = last-resort nearest centroid.")
    ap.add_argument("--cfos-dapi-overlap", "--overlap", type=float, default=0.40,
                    help="Strict DAPI/c-FOS overlap threshold: accept a c-FOS ROI only when best DAPI overlap / c-FOS ROI area is >= this fraction.")
    ap.add_argument("--cfos-blur-um", type=float, default=2.0,
                    help="Only used with --cfos-association-mode blur-aware or centroid-fallback. Expands DAPI labels by this radius to tolerate c-FOS Z-focus blur.")
    ap.add_argument("--cfos-min-direct-overlap", type=float, default=None,
                    help="Deprecated alias/override for --cfos-dapi-overlap. If provided, it overrides --cfos-dapi-overlap.")
    ap.add_argument("--cfos-min-expanded-overlap", type=float, default=0.05,
                    help="Only used with blur-aware/centroid-fallback modes. Minimum c-FOS ROI fraction inside expanded DAPI labels.")
    ap.add_argument("--cfos-max-centroid-distance-um", type=float, default=4.0,
                    help="Only used with centroid-fallback mode. Accept c-FOS ROI if centroid is within this distance of a DAPI centroid.")
    ap.add_argument("--cfos-min-area-px", type=int, default=8, help="Minimum c-FOS ROI area in pixels.")
    ap.add_argument("--cfos-max-area-px", type=int, default=0,
                    help="Optional maximum c-FOS ROI area in pixels; 0 disables.")
    ap.add_argument("--cfos-max-expanded-nuclei-touched", type=int, default=12,
                    help="Reject very diffuse c-FOS ROIs touching more than this many expanded nuclei. 0 disables.")
    ap.add_argument("--pomc-expand-um", type=float, default=7.5,
                    help="Expand DAPI labels for POMC cell-body association.")
    ap.add_argument("--npy-expand-um", type=float, default=7.5,
                    help="Expand DAPI labels for NPY cell-body association.")
    ap.add_argument("--marker-min-pixels", type=int, default=5,
                    help="Minimum marker pixels assigned to a nucleus to call POMC/NPY positive.")
    ap.add_argument("--pct-pomc-npy-overlap", type=float, default=0.20,
                    help="Remove POMC ROIs with POMC area fraction overlapping NPY >= this value.")
    ap.add_argument("--pomc-in-npy", type=float, default=0.60,
                    help="For NPY+POMC ambiguous nuclei, keep POMC only when POMC object is inside NPY by this fraction.")
    ap.add_argument("--allow-missing-channel", action="store_true", default=True,
                    help="Treat missing NPY/POMC as empty rather than failing. DAPI and c-FOS remain required for quantification.")

    # QC / plots
    qc_group = ap.add_mutually_exclusive_group()
    qc_group.add_argument("--save-qc", dest="save_qc", action="store_true", default=True)
    qc_group.add_argument("--no-qc", dest="save_qc", action="store_false")
    plot_group = ap.add_mutually_exclusive_group()
    plot_group.add_argument("--make-plots", dest="make_plots", action="store_true", default=True)
    plot_group.add_argument("--no-plots", dest="make_plots", action="store_false")
    ap.add_argument("--exclude-control-plots", action="store_true", default=True,
                    help="Exclude Control/Negative controls from condition plots and stats tables.")
    ap.add_argument("--qc-mask-dirname", default="qc_masks_strict_dapi_fos")
    ap.add_argument("--scalebar-um", type=float, default=200.0)

    # Final ARC/ME publication-plot controls. The default suite produces only compact
    # two-panel ME/ARC figures in the style requested by the user.
    ap.add_argument("--plot-suite", choices=["final", "all"], default="final",
                    help="final writes only the compact ME/ARC endpoint figures; all also writes the legacy reconstructed plots.")
    ap.add_argument("--plot-error", choices=["sd", "sem"], default="sd",
                    help="Error bars for final ARC/ME plots.")
    ap.add_argument("--plot-density-area-um2", type=float, default=100000.0,
                    help="Display density as cells per this many square microns. Exact cells/um^2 columns are still written to CSV.")
    ap.add_argument("--plot-bar-width", type=float, default=0.40,
                    help="Bar width for the compact final plots. Smaller values make the bars thinner.")
    ap.add_argument("--plot-linewidth", type=float, default=1.8,
                    help="Bold outline/error-bar/spine width for final plots.")
    ap.add_argument("--plot-legend-fontsize", type=float, default=13.0,
                    help="Legend font size for final plots.")
    ap.add_argument("--plot-title-fontsize", type=float, default=17.0,
                    help="Main title font size for final plots.")
    ap.add_argument("--plot-axis-fontsize", type=float, default=9.5,
                    help="Axis/tick font size for final plots.")
    ap.add_argument("--plot-stats-fontsize", type=float, default=7.8,
                    help="ANOVA/Mann-Whitney text size for final plots.")
    ap.add_argument("--plot-figure-width", type=float, default=8.0,
                    help="Width in inches of each compact two-panel plot.")
    ap.add_argument("--plot-figure-height", type=float, default=3.8,
                    help="Height in inches of each compact two-panel plot.")
    ap.add_argument("--plot-formats", default="pdf,png",
                    help="Comma-separated individual plot formats, e.g. pdf,png,svg. A combined multipage PDF is always written.")
    ap.add_argument("--plot-dpi", type=int, default=600,
                    help="Raster DPI for PNG output.")

    # Reconstruction / anatomical region controls
    ap.add_argument("--n-bins", type=int, default=12,
                    help="Debug/output bin count. Bins are kept only for optional QC and backwards-compatible CSVs; anatomical statistics use ARC/ME CV masks by default.")
    ap.add_argument("--region-mode", choices=["arc-me-walls", "arc-me-dl", "arc-me-cv", "bins"], default="arc-me-walls",
                    help=("arc-me-walls (recommended) fits the two third-ventricle walls and projects them ventrally; "
                          "ME is only the closed floor below the curved lumen-bottom envelope, while ventricular sidewalls and all remaining tissue are ARC. "
                          "arc-me-dl applies an optional tiny U-Net refinement but keeps ME strictly constrained to that floor-only prior; "
                          "arc-me-cv is retained as a backwards-compatible alias of arc-me-walls; bins is only a fallback."))
    ap.add_argument("--arc-me-only", dest="arc_me_only", action="store_true", default=True,
                    help="Restrict regional outputs to ARC and ME. Cells outside inferred ARC/ME remain in per-cell CSV as OUTSIDE_ARC_ME but are excluded from regional summaries.")
    ap.add_argument("--include-outside-region", dest="arc_me_only", action="store_false",
                    help="Also keep OUTSIDE_ARC_ME as a reported region. Usually avoid this for the current ARC/ME-only dataset.")
    ap.add_argument("--outside-arc-me-action", choices=["arc", "assign-by-y", "exclude"], default="arc",
                    help=("What to do with nuclei whose centroids fall outside the inferred tissue contour. "
                          "Default arc implements the requested rule: only nuclei below the closed floor curve and between the projected walls are ME; "
                          "ventricular sidewalls and every other nucleus are ARC. assign-by-y is a deprecated alias of arc. "
                          "Use exclude only for a conservative sensitivity analysis."))
    ap.add_argument("--region-bins", default="default",
                    help=("Fallback bin-to-region map for --region-mode bins. Default is ARC/ME only: "
                          "ARC for dorsal/middle bins and ME for the most ventral bins. Example: 'ARC:0-9,ME:10-11'."))
    ap.add_argument("--region-order", default="ARC,ME,VMN",
                    help="Comma-separated plotting/order preference for delineated anatomical regions. Default: ARC,ME,VMN.")
    ap.add_argument("--arc-top-frac", type=float, default=0.06,
                    help="Top of ARC as fraction of detected tissue height, measured from dorsal/top to ventral/bottom. Used by the CV ARC wedge model.")
    ap.add_argument("--me-top-frac", type=float, default=0.72,
                    help="Top of ME as fraction of detected tissue height. ME is modeled as the ventral band below this boundary.")
    ap.add_argument("--arc-top-width-um", type=float, default=180.0,
                    help="ARC triangular wedge half-side width at its dorsal/top boundary, in µm from the third-ventricle wall.")
    ap.add_argument("--arc-bottom-width-um", type=float, default=650.0,
                    help="ARC triangular wedge half-side width near the ARC/ME boundary, in µm from the third-ventricle wall.")
    ap.add_argument("--arc-ventricle-margin-um", type=float, default=8.0,
                    help="Small margin from the detected third-ventricle edge before ARC starts. Use 0-15 µm depending on segmentation alignment.")
    ap.add_argument("--me-lateral-width-um", type=float, default=2000.0,
                    help="ME lateral extent on each side from the inferred midline, in µm. Increase if ME is wider in your field.")
    ap.add_argument("--region-tissue-close-um", type=float, default=7.5,
                    help="Morphological closing radius used to convert DAPI nuclei into a smooth tissue mask for third-ventricle/ARC/ME inference.")
    ap.add_argument("--region-dl-epochs", type=int, default=8,
                    help="Epochs for the optional per-section tiny U-Net used by --region-mode arc-me-dl. The final ME label remains strictly constrained to the closed floor-only anatomical prior between the projected walls. Use 0 or --region-mode arc-me-walls for deterministic anatomy-only inference.")
    ap.add_argument("--region-dl-size", type=int, default=384,
                    help="Longest side, in pixels, used for ARC/ME U-Net self-training/inference. Higher is slower but more detailed.")
    ap.add_argument("--region-dl-lr", type=float, default=1e-3,
                    help="Learning rate for the ARC/ME U-Net self-training step.")
    ap.add_argument("--region-dl-device", choices=["auto", "cpu", "cuda"], default="cpu",
                    help="Device for ARC/ME deep-learning inference. Default is cpu to avoid CUDA-driver incompatibilities; auto uses CUDA only when safely available.")
    ap.add_argument("--region-dl-min-confidence", type=float, default=0.45,
                    help="Minimum softmax confidence for replacing the anatomical prior with the U-Net prediction. Low-confidence pixels keep the anatomical prior.")
    ap.add_argument("--region-dl-base-channels", type=int, default=16,
                    help="Base channels for the small ARC/ME U-Net. Increase only if you want a larger model.")
    ap.add_argument("--region-dl-seed", type=int, default=20260620,
                    help="Deterministic initialization seed for the DAPI-only ARC/ME U-Net.")
    region_dl_verbose_group = ap.add_mutually_exclusive_group()
    region_dl_verbose_group.add_argument("--region-dl-verbose", dest="region_dl_verbose", action="store_true", default=True,
                                         help="Log per-section ventricular-wall fitting, ME-corridor diagnostics, and optional U-Net epoch losses.")
    region_dl_verbose_group.add_argument("--no-region-dl-verbose", dest="region_dl_verbose", action="store_false",
                                         help="Silence per-section ARC/ME self-training details; the CSV/QC metadata are still written.")
    ap.add_argument("--region-dl-log-every", type=int, default=1,
                    help="Log ARC/ME U-Net training loss every N epochs when --region-dl-verbose is active.")
    ap.add_argument("--region-dl-me-class-weight", type=float, default=2.50,
                    help="Extra supervised weight for the ME class during self-training. Helps the tiny U-Net respect sparse floor tissue while keeping ARC separate.")
    ap.add_argument("--region-dl-arc-class-weight", type=float, default=1.00,
                    help="Extra supervised weight for the ARC class during self-training.")
    ap.add_argument("--region-dl-background-class-weight", type=float, default=0.20,
                    help="Background class weight during ARC/ME self-training.")
    ap.add_argument("--region-dl-cache", dest="region_dl_cache", action="store_true", default=True,
                    help="Keep/save ARC/ME region masks; reruns overwrite only when section is recomputed.")
    ap.add_argument("--no-region-dl-cache", dest="region_dl_cache", action="store_false",
                    help="Disable ARC/ME region mask cache flag; masks may still be written if --save-region-mask-tifs is enabled.")
    ap.add_argument("--ventricle-search-half-width-frac", type=float, default=0.24,
                    help="Search window around image center for the central third-ventricle gap, as fraction of image width on each side.")
    ap.add_argument("--ventricle-min-gap-um", type=float, default=30.0,
                    help="Minimum central empty gap width in µm to accept a row as containing the third ventricle.")

    # Ventricle-wall / floor-only model (v16). The projected walls provide the lateral
    # limits; ME is only the closed floor below the curved lumen-bottom envelope.
    ap.add_argument("--wall-fit-lower-fraction", type=float, default=0.55,
                    help="Fraction of the lower/main third-ventricle wall run used for robust straight-line fitting. Default 0.55 follows the ventral wall tangent while ignoring the dorsal bulb.")
    ap.add_argument("--wall-oval-mode", choices=["auto", "always", "never"], default="auto",
                    help=("How to handle an oval/closed third-ventricle lumen. auto detects an interior width maximum "
                          "and projects near-vertical tangents from the oval's lateral sides, matching the drawn schematic."))
    ap.add_argument("--wall-oval-min-bulge-ratio", type=float, default=1.16,
                    help="Minimum peak-width/end-width ratio required to call an oval lumen in auto mode.")
    ap.add_argument("--wall-oval-min-height-um", type=float, default=70.0,
                    help="Minimum supported lumen height required for automatic oval detection.")
    ap.add_argument("--wall-oval-tangent-window-frac", type=float, default=0.16,
                    help="Fraction of the supported lumen height used around the widest row to estimate oval-side tangents.")
    ap.add_argument("--wall-oval-max-tangent-slope", type=float, default=0.30,
                    help="Maximum absolute x/y slope for an oval-side tangent; keeps the extension approximately vertical.")
    oval_direction_group = ap.add_mutually_exclusive_group()
    oval_direction_group.add_argument("--wall-oval-enforce-opposite-tangents", dest="wall_oval_enforce_opposite_tangents", action="store_true", default=True,
                                      help=("For an oval lumen, force the left projected tangent to open ventrolaterally to the left "
                                            "and the right tangent to open ventrolaterally to the right. This prevents both tangents "
                                            "from drifting in the same direction and artificially shrinking the ME corridor."))
    oval_direction_group.add_argument("--no-wall-oval-enforce-opposite-tangents", dest="wall_oval_enforce_opposite_tangents", action="store_false",
                                      help="Keep the two independently fitted oval tangent signs. Intended only for sensitivity checks.")
    ap.add_argument("--wall-oval-min-outward-slope", type=float, default=0.02,
                    help="Minimum absolute outward x/y slope imposed on each oval tangent when opposite-direction correction is enabled.")
    oval_divergence_group = ap.add_mutually_exclusive_group()
    oval_divergence_group.add_argument("--wall-oval-monotonic-divergence", dest="wall_oval_monotonic_divergence", action="store_true", default=True,
                                       help="Make the oval-wall corridor stay the same width or widen ventrally; it may not reconverge below the tangent points.")
    oval_divergence_group.add_argument("--no-wall-oval-monotonic-divergence", dest="wall_oval_monotonic_divergence", action="store_false",
                                       help="Allow the fitted oval-wall corridor to reconverge ventrally.")
    ap.add_argument("--wall-fit-min-rows", type=int, default=35,
                    help="Minimum number of supported wall rows used to fit each projected ventricular wall.")
    ap.add_argument("--wall-fit-max-gap-rows", type=int, default=18,
                    help="Fill short missing runs in row-wise ventricle detection before selecting the main wall segment.")
    ap.add_argument("--wall-fit-tip-exclude-um", type=float, default=6.0,
                    help="Exclude this many micrometres from the very ventral end of the detected lumen when fitting the wall tangent, avoiding floor curvature.")
    ap.add_argument("--wall-support-probe-um", type=float, default=24.0,
                    help="Width of the DAPI-tissue probe on each side of a candidate ventricular wall. Rows lacking tissue support on either wall are rejected.")
    ap.add_argument("--wall-support-min-fraction", type=float, default=0.06,
                    help="Minimum tissue fraction in each wall-support probe for a row to be accepted as true third-ventricle wall.")
    ap.add_argument("--wall-max-gap-frac", type=float, default=0.42,
                    help="Reject central empty runs wider than this fraction of image width; such runs usually represent external background rather than the third ventricle.")
    ap.add_argument("--wall-center-tolerance-frac", type=float, default=0.08,
                    help="Maximum row-wise displacement of the detected lumen center from the global midline, as a fraction of image width.")
    ap.add_argument("--wall-max-asymmetry", type=float, default=0.40,
                    help="Maximum left/right half-gap asymmetry around the fitted midline. Rejects one-sided background gaps after the true ventricular wall ends.")
    ap.add_argument("--wall-max-abs-slope", type=float, default=1.35,
                    help="Maximum absolute fitted wall slope in pixels of x per pixel of y. Larger values are clipped for stability.")
    ap.add_argument("--wall-corridor-margin-um", type=float, default=0.0,
                    help="Signed margin applied to the projected corridor. Positive values expand ME laterally; negative values shrink it. Default 0 follows the wall lines exactly.")
    ap.add_argument("--wall-corridor-min-width-um", type=float, default=20.0,
                    help="Minimum allowed distance between projected wall lines after clamping.")
    ap.add_argument("--wall-corridor-max-width-frac", type=float, default=0.70,
                    help="Maximum projected corridor width as a fraction of image width. Prevents unstable wall fits from swallowing the whole section.")
    ap.add_argument("--wall-projection-top-margin-um", type=float, default=12.0,
                    help="Activate the ME corridor this far dorsal to the last supported ventricle-lumen row, including floor cells touching the ventricular tip.")
    ap.add_argument("--wall-projection-bottom-margin-um", type=float, default=15.0,
                    help="Allow the projected wall lines to extend this far beyond the lowest DAPI centroid/tissue row for complete floor coverage.")
    ap.add_argument("--wall-floor-open-bottom-margin-um", type=float, default=12.0,
                    help="If the central lumen component reaches this close to the projected/image bottom, the floor is open and ME is absent.")
    ap.add_argument("--wall-floor-curve-smooth-um", type=float, default=12.0,
                    help="Horizontal smoothing scale for the detected curved ventricular-floor envelope.")
    ap.add_argument("--wall-floor-max-rise-um", type=float, default=90.0,
                    help="Maximum dorsal rise retained in the curved floor envelope relative to its ventral-most portion; suppresses sidewall contamination.")
    ap.add_argument("--wall-floor-centroid-margin-um", type=float, default=0.0,
                    help="Additional ventral offset a nucleus centroid must clear below the local floor curve to be ME.")
    ap.add_argument("--wall-sidewall-exclusion-um", type=float, default=18.0,
                    help="Nuclei this close to either projected ventricular wall are assigned to ARC rather than ME near the floor transition.")
    ap.add_argument("--wall-sidewall-release-depth-um", type=float, default=35.0,
                    help="At this depth below the local floor curve, the sidewall exclusion is relaxed so deep floor tissue is retained as ME.")
    ap.add_argument("--wall-floor-component-close-um", type=float, default=10.0,
                    help="Small closing radius used to seal DAPI-derived tissue gaps before testing whether the ventricle has a closed floor.")
    sparse_floor_group = ap.add_mutually_exclusive_group()
    sparse_floor_group.add_argument("--wall-floor-sparse-rescue", dest="wall_floor_sparse_rescue", action="store_true", default=True,
                                    help=("When the empty-lumen component leaks through a genuinely present but sparsely nucleated floor, "
                                          "recover the ME floor from a bilateral, low-density DAPI bridge inside the projected walls."))
    sparse_floor_group.add_argument("--no-wall-floor-sparse-rescue", dest="wall_floor_sparse_rescue", action="store_false",
                                    help="Disable the sparse-floor rescue and require a morphologically closed lumen component.")
    ap.add_argument("--wall-floor-sparse-horizontal-link-um", type=float, default=55.0,
                    help="Horizontal linking radius used only to connect sparse floor nuclei into a candidate ME bridge.")
    ap.add_argument("--wall-floor-sparse-vertical-link-um", type=float, default=10.0,
                    help="Vertical linking radius used only for sparse-floor rescue; kept small so the open lumen is not artificially filled.")
    ap.add_argument("--wall-floor-sparse-top-margin-um", type=float, default=25.0,
                    help="Allow sparse-floor evidence to begin this far dorsal to the detected ventral lumen tip.")
    ap.add_argument("--wall-floor-sparse-curve-margin-um", type=float, default=35.0,
                    help="Lateral safety margin added around sparse floor nuclei when interpolating the floor curve, bounded by the projected walls.")
    ap.add_argument("--wall-floor-sparse-min-nuclei", type=int, default=3,
                    help="Minimum DAPI nuclei needed for sparse-floor rescue. The usual closed-floor path remains stricter.")
    ap.add_argument("--wall-floor-sparse-min-bilateral-nuclei", type=int, default=1,
                    help="Minimum original DAPI nuclei required on each side of the fitted midline for sparse-floor rescue.")
    ap.add_argument("--wall-me-min-nuclei", type=int, default=5,
                    help="Minimum DAPI nuclei with centroids inside the projected wall corridor before ME is considered present in automatic mode.")
    ap.add_argument("--wall-me-min-midline-nuclei", type=int, default=1,
                    help="Minimum enclosed nuclei near the ventricular midline. Set 0 to allow very asymmetric/partial ME sections.")
    ap.add_argument("--wall-me-midline-half-width-um", type=float, default=140.0,
                    help="Half-width around the fitted midline used for --wall-me-min-midline-nuclei.")
    ap.add_argument("--wall-me-min-depth-um", type=float, default=18.0,
                    help="Minimum ventral depth spanned by enclosed ME nuclei in automatic mode.")
    ap.add_argument("--wall-me-tissue-dilate-um", type=float, default=10.0,
                    help="Small dilation around DAPI nuclei used only to draw a continuous ME tissue contour; nucleus classification remains centroid-based.")
    ap.add_argument("--wall-me-tissue-close-um", type=float, default=12.0,
                    help="Closing radius used only to smooth the displayed/saved ME tissue contour inside the wall corridor.")
    ap.add_argument("--save-wall-projection-csv", dest="save_wall_projection_csv", action="store_true", default=True,
                    help="Write per-row projected ventricular-wall coordinates for every reconstructed section.")
    ap.add_argument("--no-wall-projection-csv", dest="save_wall_projection_csv", action="store_false",
                    help="Do not write per-row wall-projection CSV files.")
    ap.add_argument("--arc-me-boundary-mode", choices=["auto", "fraction"], default="auto",
                    help="ARC/ME boundary mode. auto uses ventricle/tissue plus weak marker/c-FOS guidance; fraction uses only --me-top-frac.")
    ap.add_argument("--auto-me-boundary", dest="auto_me_boundary", action="store_true", default=True,
                    help="Automatically refine the ARC/ME boundary from DAPI tissue density plus c-FOS/NPY/POMC density. The fixed --me-top-frac remains the fallback.")
    ap.add_argument("--no-auto-me-boundary", dest="auto_me_boundary", action="store_false",
                    help="Disable automatic ARC/ME boundary refinement and use --me-top-frac only.")
    ap.add_argument("--me-boundary-min-frac", type=float, default=0.55,
                    help="Lowest allowed ARC/ME boundary as fraction of tissue height for automatic boundary search.")
    ap.add_argument("--me-boundary-max-frac", type=float, default=0.84,
                    help="Highest allowed ARC/ME boundary as fraction of tissue height for automatic boundary search.")
    ap.add_argument("--me-boundary-fos-weight", "--region-fos-guide-weight", dest="me_boundary_fos_weight", type=float, default=0.35,
                    help="Weight of c-FOS/marker density when refining the ARC/ME boundary. 0 uses DAPI/ventricle geometry only; higher follows activation clusters more.")
    ap.add_argument("--me-detection-mode", choices=["auto", "always", "never"], default="auto",
                    help=("Ventricle-floor ME detection. auto calls ME only when a closed curved floor is found and enough DAPI nuclei lie below it; "
                          "always labels every enclosed tissue nucleus as ME; never disables ME and treats all nuclei as ARC."))
    ap.add_argument("--arc-fill-mode", choices=["tissue-minus-me", "shape-prior"], default="tissue-minus-me",
                    help=("How to define ARC after ME detection. tissue-minus-me treats all non-ME imaged tissue as ARC, "
                          "which matches these ARC/ME-only fields and avoids missing cells. shape-prior keeps the older triangle/oval ARC prior."))
    ap.add_argument("--me-floor-search-min-frac", type=float, default=0.58,
                    help="Earliest fraction of tissue height to start looking for the ventral ME floor. Lower values detect sections with a large ME.")
    ap.add_argument("--me-floor-search-max-frac", type=float, default=0.98,
                    help="Latest fraction of tissue height for the ME floor search window.")
    ap.add_argument("--me-floor-core-half-width-um", type=float, default=160.0,
                    help="Half-width around the midline used to decide whether a true ME floor is present under the third ventricle.")
    ap.add_argument("--me-floor-lateral-width-um", type=float, default=850.0,
                    help="Maximum half-width around the midline used to build the dynamic oval/contour ME mask once floor tissue is detected.")
    ap.add_argument("--me-floor-min-core-nuclei", type=int, default=10,
                    help="Minimum DAPI nuclei in the central ventral floor core required for automatic ME detection.")
    ap.add_argument("--me-floor-bridge-half-width-um", type=float, default=95.0,
                    help="Very central half-width around the third-ventricle midline used to detect a real ME floor bridge. Helps avoid false ME when only lateral ARC walls are present.")
    ap.add_argument("--me-floor-min-bridge-nuclei", type=int, default=3,
                    help="Minimum DAPI nuclei in the very central floor bridge required for automatic ME detection.")
    ap.add_argument("--me-floor-min-total-nuclei", type=int, default=45,
                    help="Minimum DAPI nuclei in the wider ventral floor window required for automatic ME detection.")
    ap.add_argument("--me-floor-min-core-density", type=float, default=0.003,
                    help="Minimum DAPI pixel density in the central floor core. This prevents drawing ME when the third-ventricle floor is absent.")
    ap.add_argument("--me-floor-threshold-rel", type=float, default=0.22,
                    help="Relative threshold applied to the smoothed central floor density profile to find the top of the ME floor.")
    ap.add_argument("--me-floor-min-run-um", type=float, default=18.0,
                    help="Minimum vertical run of floor-density-positive rows needed to accept a ME floor.")
    ap.add_argument("--me-floor-top-margin-um", type=float, default=25.0,
                    help="Small dorsal margin above the detected floor top included when drawing the dynamic ME contour.")
    ap.add_argument("--me-floor-dilate-um", type=float, default=18.0,
                    help="Dilation radius for the floor-aware ME seed before contour/oval smoothing.")
    ap.add_argument("--me-floor-close-um", type=float, default=28.0,
                    help="Closing radius for smoothing the dynamic ME contour.")
    ap.add_argument("--me-floor-oval-margin-um", type=float, default=45.0,
                    help="Extra margin around the floor nuclei when creating the dynamic ME oval. Increase if ME is undercalled.")
    ap.add_argument("--me-floor-use-convex-hull", dest="me_floor_use_convex_hull", action="store_true", default=False,
                    help="Use a convex-hull component to smooth the ME floor contour. v16 keeps this OFF by default because convex hulls can overcall dense ARC shelves as ME.")
    ap.add_argument("--no-me-floor-convex-hull", dest="me_floor_use_convex_hull", action="store_false",
                    help="Disable convex-hull smoothing of the ME floor mask.")
    hanging_group = ap.add_mutually_exclusive_group()
    hanging_group.add_argument("--me-floor-hanging-only", dest="me_floor_hanging_only", action="store_true", default=True,
                               help="Require ME to be the sparse tissue hanging from the third-ventricle floor. Dense ARC shelves below the old horizontal line are kept as ARC.")
    hanging_group.add_argument("--no-me-floor-hanging-only", dest="me_floor_hanging_only", action="store_false",
                               help="Disable hanging-floor restriction and allow a broader floor mask. Use only for sensitivity checks.")
    conn_group = ap.add_mutually_exclusive_group()
    conn_group.add_argument("--me-floor-require-connected-bridge", dest="me_floor_require_connected_bridge", action="store_true", default=True,
                            help="Require the ME floor component to be connected across/around the third-ventricle midline. If tissue hangs from only one side, ME is absent.")
    conn_group.add_argument("--no-me-floor-require-connected-bridge", dest="me_floor_require_connected_bridge", action="store_false",
                            help="Allow one-sided floor components to be called ME. Usually not recommended for these reconstructed images.")
    sparse_group = ap.add_mutually_exclusive_group()
    sparse_group.add_argument("--me-floor-use-sparse-gating", dest="me_floor_use_sparse_gating", action="store_true", default=True,
                              help="Use local DAPI density to separate dense ARC from sparse ME floor tissue.")
    sparse_group.add_argument("--no-me-floor-sparse-gating", dest="me_floor_use_sparse_gating", action="store_false",
                              help="Disable sparse-vs-dense DAPI gating for ME. Use only for QC sensitivity.")
    ap.add_argument("--me-floor-density-sigma-um", type=float, default=24.0,
                    help="Gaussian sigma in µm for DAPI centroid density used to distinguish dense ARC from sparse ME.")
    ap.add_argument("--me-floor-max-dense-arc-ratio", type=float, default=0.62,
                    help="Maximum ME/ARC local DAPI-density ratio. Lower is stricter; 0.55-0.70 usually separates sparse ME from dense ARC.")
    ap.add_argument("--me-floor-min-sparse-fraction", type=float, default=0.55,
                    help="Minimum fraction of a candidate hanging component that must be sparse relative to ARC density.")
    ap.add_argument("--me-floor-min-cross-midline-span-um", type=float, default=120.0,
                    help="Minimum horizontal span across the midline required for a connected ME floor bridge.")
    ap.add_argument("--me-floor-anchor-depth-um", type=float, default=95.0,
                    help="Depth below the detected floor top used to find the central hanging-floor anchor.")
    ap.add_argument("--me-floor-min-hanging-depth-um", type=float, default=45.0,
                    help="Minimum ventral depth of the connected hanging component before ME is accepted.")
    ap.add_argument("--me-floor-arc-reference-band-um", type=float, default=150.0,
                    help="Dorsal/juxtafloor ARC band used as high-density reference for sparse ME gating.")
    ap.add_argument("--outside-me-near-radius-um", type=float, default=18.0,
                    help="If a DAPI centroid falls outside the ME mask but very near it, assign to ME; otherwise outside centroids default to ARC, avoiding horizontal ME overcalls.")
    ap.add_argument("--arc-extend-if-no-me", dest="arc_extend_if_no_me", action="store_true", default=True,
                    help="If no ME floor is detected, classify the full visible tissue/oval field as ARC instead of drawing a fake ME line.")
    ap.add_argument("--no-arc-extend-if-no-me", dest="arc_extend_if_no_me", action="store_false",
                    help="If no ME floor is detected, keep only the ARC shape prior rather than all visible tissue.")
    ap.add_argument("--region-marker-guide-weight", type=float, default=1.0,
                    help="Weight for NPY/POMC marker centroids when expanding the ARC triangle. DAPI remains the anatomical base.")
    ap.add_argument("--arc-auto-width", dest="arc_auto_width", action="store_true", default=True,
                    help="Automatically widen the ARC triangle from DAPI + c-FOS/NPY/POMC centroids so active cells are not missed.")
    ap.add_argument("--no-arc-auto-width", dest="arc_auto_width", action="store_false",
                    help="Disable automatic ARC width estimation and use --arc-bottom-width-um / --arc-top-width-um only.")
    ap.add_argument("--arc-cell-quantile", "--arc-lateral-quantile", dest="arc_cell_quantile", type=float, default=0.985,
                    help="Quantile of cell/marker distance from the ventricle wall used to set ARC lateral width.")
    ap.add_argument("--arc-cell-margin-um", "--arc-safety-margin-um", dest="arc_cell_margin_um", type=float, default=45.0,
                    help="Extra margin added to the ARC auto-width, in µm.")
    ap.add_argument("--arc-shape-mode", choices=["adaptive", "triangle", "oval", "hybrid"], default="adaptive",
                    help=("ARC shape model. adaptive uses ventricle-following triangle when adequate, "
                          "but switches/expands to c-FOS/marker-guided oval or hybrid support in ventral/caudal cuts. "
                          "triangle reproduces the previous dark-blue triangle; oval uses only the ventricle-referenced ellipse; "
                          "hybrid always uses triangle union oval."))
    ap.add_argument("--arc-shape-min-coverage", type=float, default=0.92,
                    help="Minimum weighted c-FOS/marker centroid coverage required for the triangle before adaptive mode expands to oval/hybrid.")
    ap.add_argument("--arc-shape-switch-min-gain", type=float, default=0.05,
                    help="Minimum weighted coverage gain of the oval over the triangle required to trigger hybrid mode in --arc-shape-mode adaptive.")
    ap.add_argument("--arc-fos-guide-weight", type=float, default=None,
                    help=("Weight of c-FOS centroids for ARC shape inference. Default derives from --me-boundary-fos-weight. "
                          "Increase this if c-FOS activation should pull the ARC contour wider/oval."))
    ap.add_argument("--arc-dapi-guide-weight", type=float, default=0.35,
                    help="DAPI centroid weight for ARC oval/shape inference. Lower values prevent all DAPI nuclei from over-expanding the ROI.")
    ap.add_argument("--arc-marker-guide-weight", type=float, default=1.25,
                    help="NPY/POMC centroid weight for ARC oval/shape inference.")
    ap.add_argument("--arc-oval-quantile", type=float, default=0.92,
                    help="Weighted quantile used to set the ventricle-referenced oval radius around DAPI/c-FOS/marker centroids.")
    ap.add_argument("--arc-oval-margin-um", type=float, default=55.0,
                    help="Extra margin added to the adaptive oval ARC ROI, in µm.")
    ap.add_argument("--arc-oval-min-rx-um", type=float, default=190.0,
                    help="Minimum lateral radius of the adaptive oval ARC ROI, in µm.")
    ap.add_argument("--arc-oval-min-ry-um", type=float, default=120.0,
                    help="Minimum vertical radius of the adaptive oval ARC ROI, in µm.")
    ap.add_argument("--arc-oval-max-rx-frac", type=float, default=0.48,
                    help="Maximum adaptive oval lateral radius as a fraction of image width.")
    ap.add_argument("--arc-oval-max-ry-frac", type=float, default=0.48,
                    help="Maximum adaptive oval vertical radius as a fraction of detected tissue height.")
    ap.add_argument("--arc-include-fos-cluster-mask", dest="arc_include_fos_cluster_mask", action="store_true", default=True,
                    help="Add smoothed c-FOS/marker activation clusters to the ARC support before the neural refinement. This helps avoid missing activated ARC cells.")
    ap.add_argument("--no-arc-include-fos-cluster-mask", dest="arc_include_fos_cluster_mask", action="store_false",
                    help="Disable direct inclusion of c-FOS/marker activation clusters in the ARC support.")
    ap.add_argument("--arc-fos-cluster-quantile", type=float, default=0.82,
                    help="Quantile threshold for the smoothed c-FOS/marker density cluster mask used to expand ARC.")
    ap.add_argument("--arc-fos-cluster-dilate-um", type=float, default=25.0,
                    help="Dilation radius for c-FOS/marker activation clusters added to ARC, in µm.")
    ap.add_argument("--region-dl-preserve-geometry", dest="region_dl_preserve_geometry", action="store_true", default=True,
                    help="Constrain the deep-learning ARC/ME output to the ventricle-interpolated triangular geometry. This avoids biologically impossible contours.")
    ap.add_argument("--no-region-dl-preserve-geometry", dest="region_dl_preserve_geometry", action="store_false",
                    help="Allow the U-Net refinement to override the triangular geometry more freely. Use only after manual QC.")
    me_prior_group = ap.add_mutually_exclusive_group()
    me_prior_group.add_argument("--region-dl-strict-me-prior", dest="region_dl_strict_me_prior", action="store_true", default=True,
                                help="Do not allow U-Net refinement to create ME outside the floor prior or the explicitly configured boundary band. Recommended default.")
    me_prior_group.add_argument("--no-region-dl-strict-me-prior", dest="region_dl_strict_me_prior", action="store_false",
                                help="Allow U-Net refinement to expand ME outside the prior. Use only as a sensitivity/QC experiment.")
    ap.add_argument(
        "--region-dl-boundary-refine-um", type=float, default=0.0,
        help=("Maximum anatomy-constrained ME expansion around the deterministic prior "
              "when strict mode is active. 0 makes the deep mask exactly conservative; "
              "a small value such as 20 permits a genuine pixel-level GFAP boundary "
              "refinement while cell-region assignments remain anatomy-authoritative."),
    )
    ap.add_argument("--save-region-mask-tifs", dest="save_region_mask_tifs", action="store_true", default=True,
                    help="Save inferred ARC/ME region masks as TIFF files next to reconstructed sections.")
    ap.add_argument("--no-region-mask-tifs", dest="save_region_mask_tifs", action="store_false",
                    help="Do not save inferred ARC/ME region mask TIFFs.")
    ap.add_argument("--stitch-gap-px", type=int, default=0,
                    help="Pixel gap inserted between left and right halves in reconstructed sections.")
    ap.add_argument("--stitch-orientation", choices=["left-right", "right-left"], default="left-right",
                    help="How to place S##_left and S##_right when reconstructing the section.")
    ap.add_argument("--keep-repetitions-as-same-animal", dest="keep_repetitions_as_same_animal", action="store_true", default=True,
                    help="Treat folders like FR6-2(R) or FR721(R) as repeated sections of the same animal, not as a new N.")
    ap.add_argument("--treat-repetitions-as-new-animals", dest="keep_repetitions_as_same_animal", action="store_false",
                    help="Override: count folders with trailing parenthetical repetition labels as independent animals.")
    recon_tif_group = ap.add_mutually_exclusive_group()
    recon_tif_group.add_argument("--save-reconstructed-tifs", dest="save_reconstructed_tifs", action="store_true", default=True,
                                 help="Save reconstructed stitched channel TIFFs and label masks.")
    recon_tif_group.add_argument("--no-reconstructed-tifs", dest="save_reconstructed_tifs", action="store_false",
                                 help="Do not save reconstructed channel TIFFs/label masks.")
    ap.add_argument("--reconstructed-dirname", default="reconstructed_sections",
                    help="Output subfolder for reconstructed channel images, stitched label masks, and ARC/ME masks.")
    ap.add_argument("--reconstructed-qc-dirname", default="reconstructed_qc_arc_me",
                    help="Output subfolder for reconstructed ARC/ME QC PDFs/PNGs. No left/right QC is generated in this mode.")

    # Outlier/statistics controls
    ap.add_argument("--outlier-method", choices=["none", "iqr"], default="iqr",
                    help="Outlier method used for reports/plot markers. IQR uses Tukey fences Q1-k*IQR and Q3+k*IQR within each condition/metric/side.")
    ap.add_argument("--iqr-k", type=float, default=1.5,
                    help="Multiplier for IQR outlier fences.")
    ap.add_argument("--exclude-outliers-stats", action="store_true", default=False,
                    help="Compute plot p-values from IQR-filtered values. Plotted bars still show all values; outlier-filtered and all-value CSV stat tables are always written unless --no-stat-tables.")
    outlier_mark_group = ap.add_mutually_exclusive_group()
    outlier_mark_group.add_argument("--mark-outliers", dest="mark_outliers", action="store_true", default=True,
                                    help="Mark IQR outliers with black x symbols in plots.")
    outlier_mark_group.add_argument("--no-mark-outliers", dest="mark_outliers", action="store_false",
                                    help="Do not mark IQR outliers in plots.")
    stat_group = ap.add_mutually_exclusive_group()
    stat_group.add_argument("--write-stat-tables", dest="write_stat_tables", action="store_true", default=True,
                            help="Write condition and paired statistical CSV tables, including IQR-filtered versions.")
    stat_group.add_argument("--no-stat-tables", dest="write_stat_tables", action="store_false",
                            help="Do not write statistical CSV tables.")

    definitive_group = ap.add_mutually_exclusive_group()
    definitive_group.add_argument("--write-definitive-stats", dest="write_definitive_stats", action="store_true", default=True,
                                  help="Write definitive animal-level endpoint and statistics CSVs for ARC/ME, including count denominators and area offsets. Does not alter figures.")
    definitive_group.add_argument("--no-definitive-stats", dest="write_definitive_stats", action="store_false",
                                  help="Disable the additional definitive statistics CSVs.")
    ap.add_argument("--definitive-density-scale-um2", type=float, default=100000.0,
                    help="Scale exact density endpoints for definitive statistics, e.g. cells per 100000 µm². Raw cells/µm² are also retained.")
    ap.add_argument("--definitive-stats-primary", choices=["animal_anova_mwu"], default="animal_anova_mwu",
                    help="Primary inferential strategy: animal-level one-way ANOVA globally with exact MWU pairwise tests.")
    ap.add_argument("--verbose", action="store_true")

    return ap


def run_cmd(cmd: Sequence[str], *, check: bool = True, logger: Optional[logging.Logger] = None) -> subprocess.CompletedProcess:
    if logger:
        logger.info("RUN: %s", " ".join(map(str, cmd)))
    proc = subprocess.run(list(map(str, cmd)), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if logger and proc.stdout.strip():
        logger.info("STDOUT:\n%s", proc.stdout.rstrip())
    if logger and proc.stderr.strip():
        logger.warning("STDERR:\n%s", proc.stderr.rstrip())
    if check and proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(map(str, cmd))}\nSTDERR:\n{proc.stderr}")
    return proc


def conda_env_exists(conda_exe: str, env_name: str) -> bool:
    proc = subprocess.run([conda_exe, "env", "list"], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        return False
    for line in proc.stdout.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        first = line.split()[0]
        if first == env_name:
            return True
        # Sometimes `conda env list` shows path only; check basename too.
        if Path(first).name == env_name:
            return True
        parts = line.split()
        if parts and Path(parts[-1]).name == env_name:
            return True
    return False


def currently_in_env(env_name: str) -> bool:
    if os.environ.get("CONDA_DEFAULT_ENV") == env_name:
        return True
    if Path(sys.prefix).name == env_name:
        return True
    return False


def _conda_run_python(conda_exe: str, env_name: str, code: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [conda_exe, "run", "-n", env_name, "python", "-c", code],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )


def _env_has_required_packages(conda_exe: str, env_name: str, need_cellpose: bool, need_region_dl: bool = False) -> bool:
    modules = ["numpy", "pandas", "scipy", "skimage", "matplotlib", "tifffile", "PIL", "imagecodecs"]
    if need_cellpose:
        modules.append("cellpose")
    if need_region_dl:
        modules.append("torch")
    code = "\n".join([f"import {m}" for m in modules])
    proc = _conda_run_python(conda_exe, env_name, code)
    return proc.returncode == 0


def _install_pipeline_packages(conda_exe: str, env_name: str, *, need_cellpose: bool, need_region_dl: bool = False) -> None:
    """Install packages only inside the dedicated pipeline env."""
    core = ["numpy", "pandas", "scipy", "scikit-image", "matplotlib", "tifffile", "imagecodecs", "pillow"]
    pkgs = core + (["cellpose"] if need_cellpose else []) + (["torch"] if need_region_dl else [])
    subprocess.check_call([conda_exe, "run", "-n", env_name, "python", "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([conda_exe, "run", "-n", env_name, "python", "-m", "pip", "install", "--upgrade"] + pkgs)


def bootstrap_conda_if_needed(args: argparse.Namespace) -> None:
    """Create/use a dedicated conda environment before importing scientific modules.

    Important: the default env is paper_apotome_repro, so this script does not modify
    base or your existing cellpose_env. Package installation happens only inside
    that dedicated environment.
    """
    if not args.auto_conda:
        return

    protected_envs = {"base", "cellpose_env"}
    if str(args.env_name) in protected_envs and not getattr(args, "allow_existing_env_update", False):
        raise RuntimeError(
            f"Refusing to install/update packages in existing environment '{args.env_name}'. "
            "Use the default --env-name paper_apotome_repro, or pass --allow-existing-env-update if you really intend this."
        )

    conda_exe = shutil.which("conda")
    if conda_exe is None:
        print("[WARN] conda was not found in PATH; continuing in current Python.", file=sys.stderr)
        return

    need_cellpose = bool(args.segment_missing or args.install_cellpose or args.setup_only)
    need_region_dl = str(getattr(args, "region_mode", "arc-me-dl")) == "arc-me-dl"

    if not currently_in_env(args.env_name):
        created = False
        if not conda_env_exists(conda_exe, args.env_name):
            print(f"[SETUP] Creating dedicated conda environment: {args.env_name}")
            subprocess.check_call([conda_exe, "create", "-n", args.env_name, "python=3.10", "-y"])
            created = True

        if created or args.setup_only or not _env_has_required_packages(conda_exe, args.env_name, need_cellpose, need_region_dl):
            print(f"[SETUP] Installing pipeline packages inside {args.env_name} only")
            _install_pipeline_packages(conda_exe, args.env_name, need_cellpose=need_cellpose, need_region_dl=need_region_dl)

        if args.setup_only:
            print(f"[OK] Dedicated environment {args.env_name} is ready.")
            raise SystemExit(0)

        # Re-run this same script inside the requested env, disabling auto-conda to avoid recursion.
        argv = [a for a in sys.argv[1:] if a != "--auto-conda"]
        if "--no-auto-conda" not in argv:
            argv.append("--no-auto-conda")
        cmd = [conda_exe, "run", "-n", args.env_name, "python", str(Path(__file__).resolve())] + argv
        print("[INFO] Re-running inside dedicated conda env:")
        print("       " + " ".join(cmd))
        os.execvp(conda_exe, cmd)

    # We are already inside the dedicated env. It is safe to install missing
    # packages here because this is not base/cellpose_env unless the user
    # explicitly requested that env name.
    if not _env_has_required_packages(conda_exe, args.env_name, need_cellpose, need_region_dl):
        print(f"[SETUP] Installing missing pipeline packages inside current dedicated env: {args.env_name}")
        _install_pipeline_packages(conda_exe, args.env_name, need_cellpose=need_cellpose, need_region_dl=need_region_dl)

    if args.setup_only:
        print("[OK] Setup complete.")
        raise SystemExit(0)

def import_science_stack() -> None:
    global np, pd, plt, matplotlib, regionprops, label, binary_dilation, remove_small_objects, remove_small_holes
    global binary_opening, binary_closing, disk, expand_labels, convex_hull_image, tifffile
    global f_oneway, mannwhitneyu, wilcoxon, ttest_rel, rankdata

    try:
        import numpy as _np
        import pandas as _pd
        import tifffile as _tifffile
        import matplotlib as _matplotlib
        _matplotlib.use("Agg")
        import matplotlib.pyplot as _plt
        from skimage.measure import regionprops as _regionprops, label as _label
        from skimage.morphology import (
            binary_dilation as _binary_dilation,
            remove_small_objects as _remove_small_objects,
            remove_small_holes as _remove_small_holes,
            binary_opening as _binary_opening,
            binary_closing as _binary_closing,
            disk as _disk,
            convex_hull_image as _convex_hull_image,
        )
        from skimage.segmentation import expand_labels as _expand_labels
        from scipy.stats import f_oneway as _f_oneway, mannwhitneyu as _mannwhitneyu, wilcoxon as _wilcoxon, ttest_rel as _ttest_rel, rankdata as _rankdata
    except Exception as e:
        raise RuntimeError(
            "Missing Python dependencies. Run:\n"
            "  python keyence20x_arc_me_definitive_stats_v19.py --setup-only\n"
            "or install: numpy pandas scipy scikit-image matplotlib tifffile imagecodecs pillow torch\n"
            f"Original error: {e}"
        )

    np = _np
    pd = _pd
    tifffile = _tifffile
    matplotlib = _matplotlib
    plt = _plt
    regionprops = _regionprops
    label = _label
    binary_dilation = _binary_dilation
    remove_small_objects = _remove_small_objects
    remove_small_holes = _remove_small_holes
    binary_opening = _binary_opening
    binary_closing = _binary_closing
    disk = _disk
    expand_labels = _expand_labels
    convex_hull_image = _convex_hull_image
    f_oneway = _f_oneway
    mannwhitneyu = _mannwhitneyu
    wilcoxon = _wilcoxon
    ttest_rel = _ttest_rel
    rankdata = _rankdata


# -----------------------------------------------------------------------------
# Logging / small utilities
# -----------------------------------------------------------------------------

def setup_logger(outdir: Path, verbose: bool) -> logging.Logger:
    outdir.mkdir(parents=True, exist_ok=True)
    log_path = outdir / "keyence20x_pipeline.log"
    logger = logging.getLogger("keyence20x_pipeline")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_path, mode="a")
    fh.setFormatter(fmt)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    ch.setLevel(logging.INFO if verbose else logging.WARNING)
    logger.addHandler(ch)

    logger.info("Log file: %s", log_path)
    return logger


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _path_is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _absolute_lexical(path: Path) -> Path:
    path = path.expanduser()
    return path if path.is_absolute() else Path.cwd() / path


def _force_clean_generated_dir(
    output_dir: Path,
    *,
    forbidden_roots: Sequence[Path],
    preserved_roots: Sequence[Path] = (),
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
        directory = directory_lex.resolve(strict=False)
        for child_lex in directory_lex.iterdir():
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

def norm_token(s: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(s).lower())


def canonical_sample_id(s: Any) -> str:
    if s is None:
        return ""
    ss = str(s).strip()
    ss = re.sub(r"\(.*?\)", "", ss)
    return norm_token(ss)


def xy_to_sname(xy_num: int) -> str:
    side = "left" if xy_num % 2 == 1 else "right"
    s_idx = (xy_num + 1) // 2
    return f"S{s_idx:02d}_{side}"


def parse_xy_num(text: str) -> Optional[int]:
    m = re.search(r"XY(\d+)", text, flags=re.IGNORECASE)
    return int(m.group(1)) if m else None


def sha256_file(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def safe_copy(src: Path, dst: Path, *, overwrite: bool, dry_run: bool, logger: logging.Logger) -> str:
    if dst.exists() and not overwrite:
        return "exists_kept"
    if dry_run:
        logger.info("[DRY] copy %s -> %s", src, dst)
        return "dry_run"
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    return "copied" if not dst.exists() else "copied_or_overwritten"


def _normalize_image_array(arr: Any) -> Any:
    """Normalize microscopy image arrays to 2D grayscale or RGB(A)."""
    arr = np.asarray(arr)
    while arr.ndim > 3:
        arr = arr.max(axis=0)
    if arr.ndim == 3 and arr.shape[-1] not in (3, 4):
        # Treat (Z/Y-like, Y, X) or (C, Y, X) arrays as stacks and project.
        arr = arr.max(axis=0)
    return arr


def _read_tiff_with_pillow_fallback(path: Path) -> Any:
    """Read a TIFF robustly, including LZW-compressed Keyence TIFFs.

    tifffile needs the optional imagecodecs package for LZW-compressed TIFFs.
    The dedicated conda setup now installs imagecodecs, but this Pillow fallback
    keeps the pipeline usable even if a user disables auto-conda or runs in an
    older environment.
    """
    try:
        return tifffile.imread(str(path))
    except Exception as e:
        msg = str(e)
        needs_codecs = ("imagecodecs" in msg.lower()) or ("COMPRESSION.LZW" in msg) or ("<COMPRESSION.LZW" in msg)
        if not needs_codecs:
            raise
        try:
            from PIL import Image, ImageSequence
            with Image.open(str(path)) as im:
                frames = [np.asarray(frame.copy()) for frame in ImageSequence.Iterator(im)]
            if not frames:
                raise RuntimeError(f"Pillow opened {path}, but found no frames")
            if len(frames) == 1:
                return frames[0]
            shapes = {tuple(fr.shape) for fr in frames}
            if len(shapes) == 1:
                return np.stack(frames, axis=0)
            # Keyence TIFFs often contain a full-resolution page plus a small
            # thumbnail page. If shapes differ, keep the largest page.
            return max(frames, key=lambda a: int(np.prod(a.shape[:2])))
        except Exception as pil_e:
            raise RuntimeError(
                f"Could not read LZW-compressed TIFF: {path}\n"
                "Install the missing decoder in the dedicated environment with:\n"
                f"  conda run -n paper_apotome_repro python -m pip install --upgrade imagecodecs\n"
                "or rerun this script with --setup-only.\n"
                f"Original tifffile error: {e}\n"
                f"Pillow fallback error: {pil_e}"
            ) from e


def read_image(path: Path) -> Any:
    arr = _read_tiff_with_pillow_fallback(path)
    return _normalize_image_array(arr)


def image_to_grayscale_for_cellpose(path: Path) -> Any:
    arr = read_image(path)
    if arr.ndim == 2:
        return arr
    if arr.ndim == 3 and arr.shape[-1] in (3, 4):
        # Use max projection across RGB channels. The images are pseudo-colored display TIFFs.
        return arr[..., :3].max(axis=-1)
    return np.squeeze(arr)


# -----------------------------------------------------------------------------
# Samplesheet
# -----------------------------------------------------------------------------

def load_samplesheet(samplesheet_path: Optional[Path], outdir: Path, logger: logging.Logger) -> Any:
    if samplesheet_path is not None:
        df = pd.read_csv(samplesheet_path)
        logger.info("Loaded samplesheet: %s", samplesheet_path)
    else:
        from io import StringIO
        df = pd.read_csv(StringIO(DEFAULT_SAMPLESHEET_CSV))
        logger.info("Using embedded default samplesheet")

    required = {"sample", "cond"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Samplesheet is missing required columns: {sorted(missing)}")
    df = df.copy()
    df["sample"] = df["sample"].astype(str)
    df["sample_key"] = df["sample"].map(canonical_sample_id)
    df.to_csv(outdir / "samplesheet_used.csv", index=False)
    if samplesheet_path is None:
        (outdir / "samplesheet_embedded_default.csv").write_text(DEFAULT_SAMPLESHEET_CSV)
    return df


def attach_sample_metadata(df: Any, samplesheet: Any) -> Any:
    if df.empty:
        return df
    ss = samplesheet.copy()
    ss["sample_key"] = ss["sample"].map(canonical_sample_id)
    ss = ss.drop_duplicates("sample_key", keep="first")
    out = df.copy()
    out["sample_key"] = out["sample"].map(canonical_sample_id)
    meta_cols = [c for c in ss.columns if c != "sample"]
    out = out.merge(ss[meta_cols], on="sample_key", how="left", suffixes=("", "_samplesheet"))
    return out


# -----------------------------------------------------------------------------
# Migration / sanitization
# -----------------------------------------------------------------------------

@dataclass
class MigrationRow:
    sample: str
    src_position: str
    dst_plane: str
    src_file: str
    dst_file: str
    file_type: str
    src_channel: str
    marker: str
    action: str


def raw_position_to_plane_name(pos_dir_name: str, tile: Optional[str], mode: str) -> Optional[str]:
    xy = parse_xy_num(pos_dir_name)
    if xy is None:
        sm = RE_S_PLANE_DIR.match(pos_dir_name)
        if sm:
            base = f"S{int(sm.group('num')):02d}_{sm.group('side').lower()}"
        else:
            return None
    else:
        base = xy_to_sname(xy)

    if mode == "all" and tile:
        return f"{base}_T{tile}"
    return base


def marker_dest_filename(plane_name: str, marker: str, ext: str = ".tif") -> str:
    return f"Image_{plane_name}_{marker}{ext}"


def overlay_dest_filename(plane_name: str, ext: str = ".tif") -> str:
    return f"Image_{plane_name}_Overlay{ext}"


def seg_dest_filename(plane_name: str, marker: str) -> str:
    return f"Image_{plane_name}_{marker}_seg.npy"


def detect_tile_from_name(name: str) -> Optional[str]:
    m = RE_CH_TIF.match(name) or RE_CH_SEG.match(name)
    if m and m.group("tile"):
        return m.group("tile")
    m2 = re.search(r"_(\d{3,6})_Overlay\.(tif|tiff)$", name, flags=re.IGNORECASE)
    if m2:
        return m2.group(1)
    return None


def choose_tile_files(files: List[Path], mode: str) -> List[Path]:
    if mode == "all":
        return files
    # first/last/collapse reduce multiple tile ids. collapse currently equals last because
    # old script would overwrite repeatedly; last is a better approximation of the existing sanitized folder.
    tiles = []
    for f in files:
        t = detect_tile_from_name(f.name)
        if t:
            tiles.append(t)
    if not tiles:
        return files
    chosen = sorted(set(tiles))[0] if mode == "first" else sorted(set(tiles))[-1]
    return [f for f in files if detect_tile_from_name(f.name) in (None, chosen)]


def migrate_raw_to_sanitized(raw_root: Path, sanitized_root: Path, args: argparse.Namespace, logger: logging.Logger) -> List[MigrationRow]:
    rows: List[MigrationRow] = []
    if not raw_root.exists():
        logger.warning("Raw root does not exist: %s", raw_root)
        return rows

    if sanitized_root.exists() and not args.force_migrate:
        logger.info("Sanitized root already exists; migration skipped: %s", sanitized_root)
        return rows

    logger.info("Migrating raw root -> sanitized root")
    logger.info("Raw: %s", raw_root)
    logger.info("Sanitized: %s", sanitized_root)

    for sample_dir in sorted([p for p in raw_root.iterdir() if p.is_dir() and not p.name.startswith(".")]):
        sample_name = sample_dir.name
        # Position directories can be XY## or already S##_left/right.
        pos_dirs = [p for p in sorted(sample_dir.iterdir()) if p.is_dir() and (RE_XY_DIR.match(p.name) or RE_S_PLANE_DIR.match(p.name))]
        if not pos_dirs:
            continue

        for pos_dir in pos_dirs:
            files_all = [p for p in sorted(pos_dir.iterdir()) if p.is_file()]
            files = choose_tile_files(files_all, args.migration_tile_mode)
            for src in files:
                name = src.name
                tile = detect_tile_from_name(name)
                plane = raw_position_to_plane_name(pos_dir.name, tile, args.migration_tile_mode)
                if plane is None:
                    continue
                dst_plane_dir = sanitized_root / sample_name / plane

                # Channel TIFF
                m = RE_CH_TIF.match(name)
                if m:
                    ch = f"CH{m.group('ch')}"
                    marker = CHANNEL_TO_MARKER.get(ch, ch)
                    ext = "." + m.group("ext").lower()
                    dst = dst_plane_dir / marker_dest_filename(plane, marker, ext=ext)
                    action = safe_copy(src, dst, overwrite=args.overwrite_sanitized_files, dry_run=args.dry_run, logger=logger)
                    rows.append(MigrationRow(sample_name, pos_dir.name, plane, str(src), str(dst), "channel_tif", ch, marker, action))
                    continue

                # Channel segmentation
                m = RE_CH_SEG.match(name)
                if m:
                    ch = f"CH{m.group('ch')}"
                    marker = CHANNEL_TO_MARKER.get(ch, ch)
                    dst = dst_plane_dir / seg_dest_filename(plane, marker)
                    action = safe_copy(src, dst, overwrite=args.overwrite_sanitized_files, dry_run=args.dry_run, logger=logger)
                    rows.append(MigrationRow(sample_name, pos_dir.name, plane, str(src), str(dst), "seg_npy", ch, marker, action))
                    continue

                # Overlay
                if re.search(r"Overlay\.(tif|tiff)$", name, flags=re.IGNORECASE):
                    ext = src.suffix.lower()
                    dst = dst_plane_dir / overlay_dest_filename(plane, ext=ext)
                    action = safe_copy(src, dst, overwrite=args.overwrite_sanitized_files, dry_run=args.dry_run, logger=logger)
                    rows.append(MigrationRow(sample_name, pos_dir.name, plane, str(src), str(dst), "overlay", "", "", action))
                    continue

                # GCI and other metadata: keep one copy per plane/tile.
                if src.suffix.lower() in {".gci", ".txt", ".xml"}:
                    dst = dst_plane_dir / src.name
                    action = safe_copy(src, dst, overwrite=args.overwrite_sanitized_files, dry_run=args.dry_run, logger=logger)
                    rows.append(MigrationRow(sample_name, pos_dir.name, plane, str(src), str(dst), "metadata", "", "", action))
                    continue

        # Copy models/training dirs if present, without overwriting.
        for item in sorted(sample_dir.iterdir()):
            if item.is_dir() and not (RE_XY_DIR.match(item.name) or RE_S_PLANE_DIR.match(item.name)):
                dst_dir = sanitized_root / sample_name / item.name
                if dst_dir.exists() and not args.overwrite_sanitized_files:
                    rows.append(MigrationRow(sample_name, "", "", str(item), str(dst_dir), "dir", "", "", "exists_kept"))
                    continue
                if args.dry_run:
                    logger.info("[DRY] copytree %s -> %s", item, dst_dir)
                    rows.append(MigrationRow(sample_name, "", "", str(item), str(dst_dir), "dir", "", "", "dry_run"))
                else:
                    if dst_dir.exists():
                        shutil.rmtree(dst_dir)
                    ensure_dir(dst_dir.parent)
                    shutil.copytree(item, dst_dir, copy_function=shutil.copy2)
                    rows.append(MigrationRow(sample_name, "", "", str(item), str(dst_dir), "dir", "", "", "copied"))

    return rows


# -----------------------------------------------------------------------------
# Discovery / segmentation
# -----------------------------------------------------------------------------

@dataclass
class PlaneFiles:
    sample: str
    sample_dir: Path
    plane: str
    plane_dir: Path
    side: str
    plane_index: int
    tile: str
    images: Dict[str, Optional[Path]]
    segs: Dict[str, Optional[Path]]


def find_plane_dirs(sample_dir: Path) -> List[Path]:
    dirs = []
    for p in sorted(sample_dir.iterdir()):
        if p.is_dir() and RE_S_PLANE_DIR.match(p.name):
            dirs.append(p)
    def key(p: Path) -> Tuple[int, str, str]:
        m = RE_S_PLANE_DIR.match(p.name)
        if not m:
            return (9999, p.name, "")
        side = m.group("side").lower()
        tile = m.group("tile") or ""
        return (int(m.group("num")), side, tile)
    return sorted(dirs, key=key)


def find_marker_file(folder: Path, marker: str, is_seg: bool) -> Optional[Path]:
    marker_l = marker.lower()
    suffix = f"_{marker_l}_seg.npy" if is_seg else f"_{marker_l}.tif"
    suffix2 = f"_{marker_l}.tiff"
    hits = []
    for f in folder.iterdir():
        if not f.is_file():
            continue
        fn = f.name.lower()
        if is_seg:
            if fn.endswith(suffix):
                hits.append(f)
        else:
            if fn.endswith(suffix) or fn.endswith(suffix2):
                hits.append(f)
    if hits:
        return sorted(hits)[0]

    # Fallback to CH names if non-sanitized files remain in the plane folder.
    ch = MARKER_TO_CH.get(marker)
    if ch:
        ch_l = ch.lower()
        for f in folder.iterdir():
            if not f.is_file():
                continue
            fn = f.name.lower()
            if is_seg and fn.endswith(f"_{ch_l}_seg.npy"):
                return f
            if not is_seg and (fn.endswith(f"_{ch_l}.tif") or fn.endswith(f"_{ch_l}.tiff")):
                return f
    return None


def discover_planes(sanitized_root: Path, logger: logging.Logger) -> List[PlaneFiles]:
    planes: List[PlaneFiles] = []
    if not sanitized_root.exists():
        logger.warning("Sanitized root not found: %s", sanitized_root)
        return planes
    for sample_dir in sorted([p for p in sanitized_root.iterdir() if p.is_dir() and not p.name.startswith(".")]):
        for plane_dir in find_plane_dirs(sample_dir):
            m = RE_S_PLANE_DIR.match(plane_dir.name)
            if not m:
                continue
            images = {mk: find_marker_file(plane_dir, mk, is_seg=False) for mk in MARKERS}
            segs = {mk: find_marker_file(plane_dir, mk, is_seg=True) for mk in MARKERS}
            planes.append(
                PlaneFiles(
                    sample=sample_dir.name,
                    sample_dir=sample_dir,
                    plane=plane_dir.name,
                    plane_dir=plane_dir,
                    side=m.group("side").lower(),
                    plane_index=int(m.group("num")),
                    tile=m.group("tile") or "",
                    images=images,
                    segs=segs,
                )
            )
    logger.info("Discovered %d plane directories", len(planes))
    return planes


def plane_discovery_dataframe(planes: List[PlaneFiles]) -> Any:
    rows = []
    for p in planes:
        row = {
            "sample": p.sample,
            "plane": p.plane,
            "side": p.side,
            "plane_index": p.plane_index,
            "tile": p.tile,
            "plane_dir": str(p.plane_dir),
        }
        for mk in MARKERS:
            row[f"image_{mk}"] = str(p.images.get(mk)) if p.images.get(mk) else ""
            row[f"seg_{mk}"] = str(p.segs.get(mk)) if p.segs.get(mk) else ""
            row[f"has_image_{mk}"] = bool(p.images.get(mk))
            row[f"has_seg_{mk}"] = bool(p.segs.get(mk))
        rows.append(row)
    return pd.DataFrame(rows)


def filter_planes_by_args(planes: List[PlaneFiles], args: argparse.Namespace) -> List[PlaneFiles]:
    out = list(planes)
    if args.include_samples:
        include = {canonical_sample_id(x) for x in args.include_samples}
        out = [p for p in out if canonical_sample_id(p.sample) in include]
    if args.exclude_samples:
        exclude = {canonical_sample_id(x) for x in args.exclude_samples}
        out = [p for p in out if canonical_sample_id(p.sample) not in exclude]
    if int(args.limit_planes or 0) > 0:
        out = out[:int(args.limit_planes)]
    return out


def expected_seg_for_image(image_path: Path) -> Path:
    return image_path.with_name(image_path.stem + "_seg.npy")


def resolve_model_paths(args: argparse.Namespace, sanitized_root: Path) -> Dict[str, Optional[Path]]:
    candidates: Dict[str, List[Optional[Path]]] = {
        "POMC": [Path(args.model_pomc) if args.model_pomc else None,
                 sanitized_root / "FR6-1" / "models" / "POMC_2",
                 DEFAULT_RAW_ROOT / "FR6-1" / "models" / "POMC_2"],
        "cFOS": [Path(args.model_cfos) if args.model_cfos else None,
                 sanitized_root / "FR6-1" / "c-FOS_training" / "models" / "cfos_Keyence1",
                 DEFAULT_RAW_ROOT / "FR6-1" / "c-FOS_training" / "models" / "cfos_Keyence1"],
        "DAPI": [Path(args.model_dapi) if args.model_dapi else None, DEFAULT_MODEL_DAPI],
        "NPY":  [Path(args.model_npy) if args.model_npy else None,
                 Path(args.model_dapi) if args.model_dapi else None,
                 DEFAULT_MODEL_DAPI],
    }
    out: Dict[str, Optional[Path]] = {}
    for mk, cands in candidates.items():
        out[mk] = None
        for c in cands:
            if c and c.exists():
                out[mk] = c
                break
    return out


def get_diameter_for_marker(args: argparse.Namespace, marker: str) -> Optional[float]:
    return {
        "DAPI": args.diameter_dapi,
        "POMC": args.diameter_pomc,
        "NPY": args.diameter_npy,
        "cFOS": args.diameter_cfos,
    }.get(marker)


def save_cellpose_seg(seg_path: Path, masks: Any, flows: Any = None, styles: Any = None, diams: Any = None) -> None:
    obj = {"masks": np.asarray(masks), "flows": flows, "styles": styles, "diams": diams}
    np.save(str(seg_path), obj, allow_pickle=True)


def run_cellpose_one_image(image_path: Path, seg_path: Path, model_path: Path, marker: str, args: argparse.Namespace, logger: logging.Logger) -> bool:
    if args.dry_run:
        logger.info("[DRY] would segment missing %s mask: %s -> %s", marker, image_path, seg_path)
        return True

    try:
        from cellpose import models
    except Exception as e:
        raise RuntimeError(
            "Cellpose is required to segment missing masks. Activate/install cellpose_env, or run with --no-segment-missing.\n"
            f"Original import error: {e}"
        )

    img = image_to_grayscale_for_cellpose(image_path)
    gpu_arg = bool(args.gpu)
    try:
        model = models.CellposeModel(gpu=gpu_arg, pretrained_model=str(model_path), device=None)
    except TypeError:
        model = models.CellposeModel(gpu=gpu_arg, pretrained_model=str(model_path))

    eval_kwargs: Dict[str, Any] = {}
    diameter = get_diameter_for_marker(args, marker)
    if diameter is not None:
        eval_kwargs["diameter"] = diameter
    if args.flow_threshold is not None:
        eval_kwargs["flow_threshold"] = args.flow_threshold
    if args.cellprob_threshold is not None:
        eval_kwargs["cellprob_threshold"] = args.cellprob_threshold
    if args.min_size_seg is not None:
        eval_kwargs["min_size"] = args.min_size_seg

    # channels=[0,0] is safe for single-channel 2D input. Cellpose v4 may ignore it.
    try:
        result = model.eval(img, channels=[0, 0], **eval_kwargs)
    except TypeError:
        result = model.eval(img, **eval_kwargs)

    if isinstance(result, tuple):
        masks = result[0]
        flows = result[1] if len(result) > 1 else None
        styles = result[2] if len(result) > 2 else None
        diams = result[3] if len(result) > 3 else None
    else:
        masks, flows, styles, diams = result, None, None, None

    save_cellpose_seg(seg_path, masks, flows=flows, styles=styles, diams=diams)
    logger.info("Segmented %s: %s", marker, seg_path)
    return True


def segment_missing_masks(planes: List[PlaneFiles], args: argparse.Namespace, logger: logging.Logger, outdir: Path) -> Any:
    model_paths = resolve_model_paths(args, Path(args.sanitized_root))
    rows = []
    for mk, mp in model_paths.items():
        logger.info("Model %s: %s", mk, mp if mp else "NOT FOUND")

    for p in planes:
        for marker in MARKERS:
            img = p.images.get(marker)
            seg = p.segs.get(marker)
            if img is None:
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": "no_image", "image": "", "seg": ""})
                continue
            expected = expected_seg_for_image(img)
            if seg is not None and seg.exists() and not args.overwrite_seg:
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": "existing_seg_kept", "image": str(img), "seg": str(seg)})
                continue
            if expected.exists() and not args.overwrite_seg:
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": "expected_seg_exists_kept", "image": str(img), "seg": str(expected)})
                continue
            mp = model_paths.get(marker)
            if mp is None or not mp.exists():
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": "missing_model_skip", "image": str(img), "seg": str(expected)})
                logger.warning("Missing model for %s; cannot segment %s", marker, img)
                continue
            try:
                ok = run_cellpose_one_image(img, expected, mp, marker, args, logger)
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": "segmented" if ok else "segment_failed", "image": str(img), "seg": str(expected)})
            except Exception as e:
                rows.append({"sample": p.sample, "plane": p.plane, "marker": marker, "action": f"error:{e}", "image": str(img), "seg": str(expected)})
                logger.exception("Cellpose failed for %s", img)

    df = pd.DataFrame(rows)
    df.to_csv(outdir / "segmentation_report.csv", index=False)
    return df


# -----------------------------------------------------------------------------
# Segmentation mask loading and analysis helpers
# -----------------------------------------------------------------------------

def load_seg_mask(path: Path) -> Any:
    obj = np.load(str(path), allow_pickle=True)
    if isinstance(obj, np.ndarray) and obj.dtype == object and obj.shape == ():
        obj = obj.item()
    if isinstance(obj, dict) and "masks" in obj:
        return np.asarray(obj["masks"])
    if isinstance(obj, tuple) and len(obj) > 0:
        return np.asarray(obj[0])
    if isinstance(obj, list) and len(obj) > 0:
        return np.asarray(obj[0])
    return np.asarray(obj)


def ensure_label_image(mask: Any) -> Any:
    arr = np.asarray(mask)
    if arr.ndim > 2:
        # For Cellpose masks saved as H,W or stacks; use first plane if needed.
        if arr.shape[-1] in (3, 4):
            arr = arr[..., 0]
        else:
            arr = arr[0]
    if np.issubdtype(arr.dtype, np.integer):
        return arr.astype(np.int32, copy=False)
    return label(arr > 0, connectivity=1).astype(np.int32)


def count_labels(lbl: Any) -> int:
    if lbl is None:
        return 0
    return int(np.max(lbl)) if np.size(lbl) else 0


def filter_rois_by_area(lbl: Any, min_area: int = 0, max_area: int = 0) -> Any:
    if lbl is None or int(lbl.max()) <= 0:
        return lbl
    out = np.zeros_like(lbl, dtype=np.int32)
    current = 0
    for r in regionprops(lbl):
        area = int(r.area)
        if min_area and area < min_area:
            continue
        if max_area and area > max_area:
            continue
        current += 1
        coords = r.coords
        out[coords[:, 0], coords[:, 1]] = current
    return out


def filter_pomc_rois_by_npy_overlap(pomc_lbl: Any, npy_mask: Any, overlap_thr: float) -> Tuple[Any, int]:
    if pomc_lbl is None or int(pomc_lbl.max()) <= 0 or overlap_thr <= 0:
        return pomc_lbl, 0
    out = pomc_lbl.copy().astype(np.int32)
    npy_bin = np.asarray(npy_mask > 0, dtype=bool)
    removed = 0
    for r in regionprops(pomc_lbl):
        coords = r.coords
        frac = float(np.count_nonzero(npy_bin[coords[:, 0], coords[:, 1]])) / float(max(1, r.area))
        if frac >= overlap_thr:
            out[coords[:, 0], coords[:, 1]] = 0
            removed += 1
    out = label(out > 0, connectivity=1).astype(np.int32)
    return out, removed


def assigned_pixel_counts(expanded_dapi_lbl: Any, marker_mask: Any) -> Dict[int, int]:
    vals = expanded_dapi_lbl[np.asarray(marker_mask, dtype=bool)]
    vals = vals[vals > 0]
    if vals.size == 0:
        return {}
    counts = np.bincount(vals.astype(np.int64))
    return {int(i): int(c) for i, c in enumerate(counts) if i > 0 and c > 0}


def nuclei_positive_by_expanded_label(dapi_lbl: Any, marker_lbl: Any, expand_px: int, min_pixels: int) -> Set[int]:
    if dapi_lbl is None or marker_lbl is None or int(dapi_lbl.max()) <= 0 or int(marker_lbl.max()) <= 0:
        return set()
    expanded = expand_labels(dapi_lbl.astype(np.int32), distance=int(max(0, expand_px)))
    counts = assigned_pixel_counts(expanded, marker_lbl > 0)
    return {lab for lab, cnt in counts.items() if cnt >= int(min_pixels)}


def build_centroid_table(lbl: Any) -> Dict[int, Tuple[float, float]]:
    return {int(r.label): (float(r.centroid[0]), float(r.centroid[1])) for r in regionprops(lbl)}


def nearest_nucleus_by_centroid(y: float, x: float, centroids: Dict[int, Tuple[float, float]], max_dist_px: float) -> Tuple[int, float]:
    best_lab = 0
    best_d = float("inf")
    for lab, (cy, cx) in centroids.items():
        d = math.hypot(float(y) - cy, float(x) - cx)
        if d < best_d:
            best_lab = lab
            best_d = d
    if best_d <= max_dist_px:
        return best_lab, best_d
    return 0, best_d


@dataclass
class FosAssignment:
    fos_roi: int
    nucleus: int
    accepted: bool
    reason: str
    roi_area_px: int
    direct_overlap_px: int
    expanded_overlap_px: int
    direct_overlap_frac: float
    expanded_overlap_frac: float
    n_expanded_nuclei_touched: int
    centroid_distance_px: float


def assign_cfos_rois_bluraware(dapi_lbl: Any, fos_lbl: Any, args: argparse.Namespace) -> Tuple[Set[int], List[FosAssignment]]:
    """
    Assign each c-FOS ROI to at most one DAPI nucleus.

    Default mode is STRICT because the biological readout should be nuclear:
      --cfos-association-mode strict
      --cfos-dapi-overlap 0.40

    Modes:
      strict            : original DAPI-mask overlap only.
      blur-aware        : strict first, then expanded DAPI labels using --cfos-blur-um.
      centroid-fallback : blur-aware, then nearest DAPI centroid as last resort.

    Expanded labels from skimage.segmentation.expand_labels are non-overlapping, so a blurred
    c-FOS object cannot be counted for multiple nuclei.
    """
    fos_positive_nuclei: Set[int] = set()
    assignments: List[FosAssignment] = []

    if dapi_lbl is None or fos_lbl is None or int(dapi_lbl.max()) <= 0 or int(fos_lbl.max()) <= 0:
        return fos_positive_nuclei, assignments

    mode = str(getattr(args, "cfos_association_mode", "strict")).lower()
    if mode not in {"strict", "blur-aware", "centroid-fallback"}:
        mode = "strict"

    # Backward-compatible override: if --cfos-min-direct-overlap was passed explicitly, use it.
    direct_thr = getattr(args, "cfos_dapi_overlap", None)
    if direct_thr is None:
        direct_thr = 0.40
    if getattr(args, "cfos_min_direct_overlap", None) is not None:
        direct_thr = float(args.cfos_min_direct_overlap)
    direct_thr = float(max(0.0, min(1.0, direct_thr)))

    um_per_px = float(args.um_per_px)
    blur_px = int(round(float(args.cfos_blur_um) / um_per_px)) if um_per_px > 0 else int(round(args.cfos_blur_um))
    blur_px = max(0, blur_px)
    expanded_dapi = expand_labels(dapi_lbl.astype(np.int32), distance=blur_px) if mode in {"blur-aware", "centroid-fallback"} and blur_px > 0 else dapi_lbl
    centroids = build_centroid_table(dapi_lbl) if mode == "centroid-fallback" else {}
    max_centroid_px = float(args.cfos_max_centroid_distance_um) / um_per_px if um_per_px > 0 else float(args.cfos_max_centroid_distance_um)
    max_nuc_touched = int(args.cfos_max_expanded_nuclei_touched or 0)
    expanded_thr = float(max(0.0, min(1.0, float(args.cfos_min_expanded_overlap))))

    for r in regionprops(fos_lbl):
        roi_lab = int(r.label)
        roi_area = int(r.area)
        coords = r.coords
        ys = coords[:, 0]
        xs = coords[:, 1]
        reason = "rejected"
        accepted = False
        assigned_nuc = 0
        centroid_distance = float("nan")

        if roi_area < int(args.cfos_min_area_px):
            assignments.append(FosAssignment(roi_lab, 0, False, "too_small", roi_area, 0, 0, 0.0, 0.0, 0, float("nan")))
            continue
        if int(args.cfos_max_area_px or 0) > 0 and roi_area > int(args.cfos_max_area_px):
            assignments.append(FosAssignment(roi_lab, 0, False, "too_large", roi_area, 0, 0, 0.0, 0.0, 0, float("nan")))
            continue

        # Direct DAPI overlap: this is the strict nuclear c-FOS rule.
        direct_vals = dapi_lbl[ys, xs]
        direct_counts = np.bincount(direct_vals.astype(np.int64), minlength=int(dapi_lbl.max()) + 1)
        if direct_counts.size > 0:
            direct_counts[0] = 0
        direct_best = int(np.argmax(direct_counts)) if direct_counts.size else 0
        direct_overlap = int(direct_counts[direct_best]) if direct_best > 0 else 0
        direct_frac = direct_overlap / float(max(1, roi_area))

        # Expanded DAPI overlap is computed for reporting even in strict mode, but not accepted in strict mode.
        exp_vals = expanded_dapi[ys, xs]
        exp_labels_unique = set(int(v) for v in np.unique(exp_vals) if int(v) > 0)
        n_exp_touch = len(exp_labels_unique)
        exp_counts = np.bincount(exp_vals.astype(np.int64), minlength=int(dapi_lbl.max()) + 1)
        if exp_counts.size > 0:
            exp_counts[0] = 0
        exp_best = int(np.argmax(exp_counts)) if exp_counts.size else 0
        exp_overlap = int(exp_counts[exp_best]) if exp_best > 0 else 0
        exp_frac = exp_overlap / float(max(1, roi_area))

        if max_nuc_touched > 0 and n_exp_touch > max_nuc_touched:
            assignments.append(FosAssignment(roi_lab, 0, False, "too_diffuse_many_nuclei", roi_area, direct_overlap, exp_overlap, direct_frac, exp_frac, n_exp_touch, float("nan")))
            continue

        if direct_best > 0 and direct_frac >= direct_thr:
            assigned_nuc = direct_best
            accepted = True
            reason = "strict_direct_dapi_overlap"
        elif mode in {"blur-aware", "centroid-fallback"} and exp_best > 0 and exp_frac >= expanded_thr:
            assigned_nuc = exp_best
            accepted = True
            reason = "expanded_dapi_blur_overlap"
        elif mode == "centroid-fallback":
            cy, cx = r.centroid
            nearest, dist = nearest_nucleus_by_centroid(float(cy), float(cx), centroids, max_centroid_px)
            centroid_distance = float(dist)
            if nearest > 0:
                assigned_nuc = nearest
                accepted = True
                reason = "nearest_centroid_blur_fallback"
            else:
                reason = "no_dapi_association"
        else:
            reason = "strict_overlap_below_threshold"

        if accepted and assigned_nuc > 0:
            fos_positive_nuclei.add(int(assigned_nuc))

        assignments.append(
            FosAssignment(
                fos_roi=roi_lab,
                nucleus=int(assigned_nuc),
                accepted=bool(accepted),
                reason=reason,
                roi_area_px=roi_area,
                direct_overlap_px=direct_overlap,
                expanded_overlap_px=exp_overlap,
                direct_overlap_frac=float(direct_frac),
                expanded_overlap_frac=float(exp_frac),
                n_expanded_nuclei_touched=int(n_exp_touch),
                centroid_distance_px=float(centroid_distance) if np.isfinite(centroid_distance) else float("nan"),
            )
        )

    return fos_positive_nuclei, assignments

def best_marker_label_in_nucleus(nuc_prop: Any, marker_lbl: Any) -> int:
    sl = nuc_prop.slice
    nuc_mask = nuc_prop.image
    vals = marker_lbl[sl][nuc_mask]
    vals = vals[vals > 0]
    if vals.size == 0:
        return 0
    counts = np.bincount(vals.astype(np.int64))
    counts[0] = 0
    return int(np.argmax(counts)) if counts.size else 0


def marker_overlap_fraction(marker_a_lbl: Any, marker_b_lbl: Any, label_a: int, label_b: int) -> float:
    if label_a <= 0 or label_b <= 0:
        return 0.0
    mask_a = marker_a_lbl == int(label_a)
    area = int(np.count_nonzero(mask_a))
    if area <= 0:
        return 0.0
    inter = int(np.count_nonzero(mask_a & (marker_b_lbl == int(label_b))))
    return float(inter) / float(area)


def apply_strict_pomc_vs_npy(dapi_lbl: Any, npy_set: Set[int], pomc_set: Set[int], npy_lbl: Any, pomc_lbl: Any, thr: float) -> Tuple[Set[int], int]:
    if thr <= 0:
        return set(pomc_set), 0
    both = set(npy_set) & set(pomc_set)
    if not both:
        return set(pomc_set), 0

    props = {int(p.label): p for p in regionprops(dapi_lbl) if int(p.label) in both}
    out = set(pomc_set)
    dropped = 0
    for nuc in both:
        prop = props.get(int(nuc))
        if prop is None:
            out.discard(nuc)
            dropped += 1
            continue
        best_npy = best_marker_label_in_nucleus(prop, npy_lbl)
        best_pomc = best_marker_label_in_nucleus(prop, pomc_lbl)
        if best_npy <= 0 or best_pomc <= 0:
            out.discard(nuc)
            dropped += 1
            continue
        frac = marker_overlap_fraction(pomc_lbl, npy_lbl, best_pomc, best_npy)
        if frac < thr:
            out.discard(nuc)
            dropped += 1
    return out, dropped


# -----------------------------------------------------------------------------
# Quantification
# -----------------------------------------------------------------------------

@dataclass
class PlaneQuantResult:
    sample: str
    plane: str
    side: str
    plane_index: int
    tile: str
    plane_dir: str
    dapi_nuclei: int
    cfos_rois_total: int
    cfos_rois_accepted: int
    cfos_cells: int
    npy_cells: int
    pomc_cells_raw: int
    pomc_cells: int
    pomc_rois_removed_by_npy: int
    pomc_nuclei_dropped_strict_npy: int
    cfos_only: int
    cfos_npy: int
    cfos_pomc: int
    cfos_npy_pomc: int
    missing_required: str
    status: str


def empty_label_like(shape: Tuple[int, int]) -> Any:
    return np.zeros(shape, dtype=np.int32)


def load_plane_labels(p: PlaneFiles) -> Tuple[Optional[Any], Optional[Any], Optional[Any], Optional[Any], List[str]]:
    missing: List[str] = []
    dapi_path = p.segs.get("DAPI")
    cfos_path = p.segs.get("cFOS")
    npy_path = p.segs.get("NPY")
    pomc_path = p.segs.get("POMC")

    dapi_lbl = None
    cfos_lbl = None
    npy_lbl = None
    pomc_lbl = None
    if dapi_path is None or not dapi_path.exists():
        missing.append("DAPI")
    else:
        dapi_lbl = ensure_label_image(load_seg_mask(dapi_path))
    if cfos_path is None or not cfos_path.exists():
        missing.append("cFOS")
    else:
        cfos_lbl = ensure_label_image(load_seg_mask(cfos_path))

    # Optional markers become empty masks with DAPI shape.
    if dapi_lbl is not None:
        shape = dapi_lbl.shape
        if npy_path is not None and npy_path.exists():
            npy_lbl = ensure_label_image(load_seg_mask(npy_path))
        else:
            npy_lbl = empty_label_like(shape)
        if pomc_path is not None and pomc_path.exists():
            pomc_lbl = ensure_label_image(load_seg_mask(pomc_path))
        else:
            pomc_lbl = empty_label_like(shape)

    return dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl, missing


def validate_same_shape(arrays: Sequence[Any]) -> bool:
    shapes = [tuple(a.shape) for a in arrays if a is not None]
    return len(set(shapes)) <= 1


def process_plane_quant(p: PlaneFiles, args: argparse.Namespace, logger: logging.Logger, qc_dir: Path) -> Tuple[PlaneQuantResult, List[Dict[str, Any]]]:
    assignment_rows: List[Dict[str, Any]] = []
    dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl, missing = load_plane_labels(p)

    if "DAPI" in missing or "cFOS" in missing:
        res = PlaneQuantResult(p.sample, p.plane, p.side, p.plane_index, p.tile, str(p.plane_dir),
                               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ";".join(missing), "missing_required")
        return res, assignment_rows

    if not validate_same_shape([dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl]):
        res = PlaneQuantResult(p.sample, p.plane, p.side, p.plane_index, p.tile, str(p.plane_dir),
                               count_labels(dapi_lbl), count_labels(cfos_lbl), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                               "shape_mismatch", "shape_mismatch")
        return res, assignment_rows

    # Ensure labels and remove tiny c-FOS noise.
    dapi_lbl = ensure_label_image(dapi_lbl)
    cfos_lbl = ensure_label_image(cfos_lbl)
    cfos_lbl = filter_rois_by_area(cfos_lbl, min_area=int(args.cfos_min_area_px), max_area=int(args.cfos_max_area_px or 0))
    npy_lbl = ensure_label_image(npy_lbl)
    pomc_lbl = ensure_label_image(pomc_lbl)

    dapi_n = count_labels(dapi_lbl)
    cfos_rois_total = count_labels(cfos_lbl)
    npy_rois_total = count_labels(npy_lbl)
    pomc_rois_raw = count_labels(pomc_lbl)

    # POMC/NPY de-overlap before nucleus assignment.
    pomc_lbl_clean, pomc_removed = filter_pomc_rois_by_npy_overlap(pomc_lbl, npy_lbl > 0, float(args.pct_pomc_npy_overlap))
    pomc_rois_clean = count_labels(pomc_lbl_clean)

    # c-FOS blur-aware assignment.
    cfos_nuclei, fos_assignments = assign_cfos_rois_bluraware(dapi_lbl, cfos_lbl, args)
    for fa in fos_assignments:
        row = asdict(fa)
        row.update({"sample": p.sample, "plane": p.plane, "side": p.side, "plane_index": p.plane_index, "tile": p.tile})
        assignment_rows.append(row)

    # POMC/NPY DAPI-associated assignment.
    um_per_px = float(args.um_per_px)
    pomc_expand_px = int(round(float(args.pomc_expand_um) / um_per_px)) if um_per_px > 0 else int(round(args.pomc_expand_um))
    npy_expand_px = int(round(float(args.npy_expand_um) / um_per_px)) if um_per_px > 0 else int(round(args.npy_expand_um))
    pomc_nuclei = nuclei_positive_by_expanded_label(dapi_lbl, pomc_lbl_clean, pomc_expand_px, int(args.marker_min_pixels))
    npy_nuclei = nuclei_positive_by_expanded_label(dapi_lbl, npy_lbl, npy_expand_px, int(args.marker_min_pixels))

    pomc_nuclei_strict, n_pomc_dropped = apply_strict_pomc_vs_npy(
        dapi_lbl, npy_nuclei, pomc_nuclei, npy_lbl, pomc_lbl_clean, float(args.pomc_in_npy)
    )
    pomc_nuclei = pomc_nuclei_strict

    triple = cfos_nuclei & npy_nuclei & pomc_nuclei
    cfos_npy_only = (cfos_nuclei & npy_nuclei) - pomc_nuclei
    cfos_pomc_only = (cfos_nuclei & pomc_nuclei) - npy_nuclei
    cfos_only = cfos_nuclei - npy_nuclei - pomc_nuclei

    res = PlaneQuantResult(
        sample=p.sample,
        plane=p.plane,
        side=p.side,
        plane_index=p.plane_index,
        tile=p.tile,
        plane_dir=str(p.plane_dir),
        dapi_nuclei=int(dapi_n),
        cfos_rois_total=int(cfos_rois_total),
        cfos_rois_accepted=int(sum(1 for x in fos_assignments if x.accepted)),
        cfos_cells=int(len(cfos_nuclei)),
        npy_cells=int(len(npy_nuclei)),
        pomc_cells_raw=int(len(nuclei_positive_by_expanded_label(dapi_lbl, pomc_lbl, pomc_expand_px, int(args.marker_min_pixels)))),
        pomc_cells=int(len(pomc_nuclei)),
        pomc_rois_removed_by_npy=int(pomc_removed),
        pomc_nuclei_dropped_strict_npy=int(n_pomc_dropped),
        cfos_only=int(len(cfos_only)),
        cfos_npy=int(len(cfos_npy_only)),
        cfos_pomc=int(len(cfos_pomc_only)),
        cfos_npy_pomc=int(len(triple)),
        missing_required="",
        status="ok",
    )

    if args.save_qc:
        try:
            sample_qc = qc_dir / p.sample
            ensure_dir(sample_qc)
            render_qc_pdf(
                sample_qc / f"{p.plane}_qc_strict_dapi_fos.pdf",
                dapi_lbl=dapi_lbl,
                sample=p.sample,
                plane=p.plane,
                cfos_nuclei=cfos_nuclei,
                npy_nuclei=npy_nuclei,
                pomc_nuclei=pomc_nuclei,
                cfos_only=cfos_only,
                cfos_npy=cfos_npy_only,
                cfos_pomc=cfos_pomc_only,
                triple=triple,
                res=res,
                args=args,
            )
        except Exception:
            logger.exception("QC rendering failed for %s", p.plane_dir)

    return res, assignment_rows


def render_qc_pdf(
    out_path: Path,
    *,
    dapi_lbl: Any,
    sample: str,
    plane: str,
    cfos_nuclei: Set[int],
    npy_nuclei: Set[int],
    pomc_nuclei: Set[int],
    cfos_only: Set[int],
    cfos_npy: Set[int],
    cfos_pomc: Set[int],
    triple: Set[int],
    res: PlaneQuantResult,
    args: argparse.Namespace,
) -> None:
    H, W = dapi_lbl.shape
    max_lab = int(dapi_lbl.max())
    cls = np.zeros(max_lab + 1, dtype=np.uint8)

    # 1 NPY, 2 POMC, 3 cFOS-only, 4 cFOS+NPY, 5 cFOS+POMC, 6 triple
    for n in npy_nuclei:
        if 0 < n <= max_lab:
            cls[n] = 1
    for n in pomc_nuclei:
        if 0 < n <= max_lab:
            cls[n] = 2
    for n in cfos_only:
        if 0 < n <= max_lab:
            cls[n] = 3
    for n in cfos_npy:
        if 0 < n <= max_lab:
            cls[n] = 4
    for n in cfos_pomc:
        if 0 < n <= max_lab:
            cls[n] = 5
    for n in triple:
        if 0 < n <= max_lab:
            cls[n] = 6

    class_img = cls[dapi_lbl]
    rgb = np.zeros((H, W, 3), dtype=np.float32)
    rgb[dapi_lbl > 0] = (1, 1, 1)  # all DAPI nuclei white if negative

    colors = {
        1: (0.00, 0.85, 0.00),   # NPY green
        2: (1.00, 0.00, 0.00),   # POMC red
        3: (0.56, 0.00, 1.00),   # cFOS violet
        4: (1.00, 1.00, 0.00),   # cFOS NPY yellow
        5: (244/255, 187/255, 1.0),  # cFOS POMC pink
        6: (0.00, 1.00, 1.00),   # triple cyan
    }
    for k, col in colors.items():
        rgb[class_img == k] = col

    img_height = 8.0
    img_width = img_height * (W / max(1, H))
    panel_width = 5.2
    fig = plt.figure(figsize=(img_width + panel_width, img_height), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[img_width, panel_width])
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(rgb, interpolation="nearest")
    ax.set_axis_off()

    # scale bar
    if args.um_per_px and args.um_per_px > 0:
        bar_px = int(round(float(args.scalebar_um) / float(args.um_per_px)))
        if bar_px > 0 and bar_px < W:
            x0 = W - bar_px - 50
            y0 = H - 60
            ax.add_patch(plt.Rectangle((x0, y0), bar_px, 12, color="white", ec="black", lw=1.5))
            ax.text(x0 + bar_px / 2, y0 - 10, f"{int(args.scalebar_um)} µm", color="white", ha="center", va="bottom", fontsize=16, weight="bold",
                    path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")])

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_axis_off()
    lines = [
        f"{sample}",
        f"{plane}",
        "",
        f"DAPI nuclei: {res.dapi_nuclei}",
        f"c-FOS ROIs: {res.cfos_rois_total}",
        f"c-FOS accepted: {res.cfos_rois_accepted}",
        f"c-FOS+ nuclei: {res.cfos_cells}",
        f"POMC+ nuclei: {res.pomc_cells}",
        f"NPY+ nuclei: {res.npy_cells}",
        "",
        f"c-FOS only: {res.cfos_only}",
        f"c-FOS∧POMC: {res.cfos_pomc}",
        f"c-FOS∧NPY: {res.cfos_npy}",
        f"triple: {res.cfos_npy_pomc}",
        "",
        f"c-FOS DAPI expansion: {args.cfos_blur_um:.2f} µm",
        f"POMC DAPI expansion: {args.pomc_expand_um:.2f} µm",
        f"POMC ROIs removed by NPY: {res.pomc_rois_removed_by_npy}",
    ]
    ax2.text(0.02, 0.98, "\n".join(lines), ha="left", va="top", fontsize=15, weight="bold")
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=colors[2], label="POMC"),
        Patch(facecolor=colors[1], label="NPY"),
        Patch(facecolor=colors[3], label="c-FOS only"),
        Patch(facecolor=colors[5], label="c-FOS∧POMC"),
        Patch(facecolor=colors[4], label="c-FOS∧NPY"),
        Patch(facecolor=colors[6], label="Triple"),
    ]
    ax2.legend(handles=handles, loc="lower left", fontsize=13, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def aggregate_summary(per_plane_df: Any) -> Any:
    if per_plane_df.empty:
        return per_plane_df
    ok = per_plane_df[per_plane_df["status"] == "ok"].copy()
    if ok.empty:
        return pd.DataFrame()

    metrics = [
        "dapi_nuclei", "cfos_rois_total", "cfos_rois_accepted", "cfos_cells",
        "npy_cells", "pomc_cells_raw", "pomc_cells", "pomc_rois_removed_by_npy",
        "pomc_nuclei_dropped_strict_npy", "cfos_only", "cfos_npy", "cfos_pomc", "cfos_npy_pomc"
    ]
    rows = []
    for sample, sdf in ok.groupby("sample", sort=True):
        row: Dict[str, Any] = {"sample": sample, "n_planes_total": int(sdf.shape[0])}
        for side in ["left", "right"]:
            ss = sdf[sdf["side"] == side]
            row[f"n_planes_{side}"] = int(ss.shape[0])
            for m in metrics:
                row[f"{m}_{side}"] = int(pd.to_numeric(ss[m], errors="coerce").fillna(0).sum()) if not ss.empty else 0
        for m in metrics:
            row[f"{m}_total"] = int(pd.to_numeric(sdf[m], errors="coerce").fillna(0).sum())

        # Old-compatible aliases
        row["cfos_pomc_left"] = row.get("cfos_pomc_left", 0)
        row["cfos_pomc_right"] = row.get("cfos_pomc_right", 0)
        row["cfos_pomc_total"] = row.get("cfos_pomc_total", 0)
        row["cfos_npy_left"] = row.get("cfos_npy_left", 0)
        row["cfos_npy_right"] = row.get("cfos_npy_right", 0)
        row["cfos_npy_total"] = row.get("cfos_npy_total", 0)
        rows.append(row)
    return pd.DataFrame(rows)


def run_quantification(planes: List[PlaneFiles], args: argparse.Namespace, logger: logging.Logger, outdir: Path, samplesheet: Any) -> Tuple[Any, Any, Any]:
    # Refresh segmentation paths after segmentation stage, then apply any debug/sample filters.
    planes = filter_planes_by_args(discover_planes(Path(args.sanitized_root), logger), args)
    per_plane: List[Dict[str, Any]] = []
    assignment_rows: List[Dict[str, Any]] = []
    qc_dir = outdir / args.qc_mask_dirname
    ensure_dir(qc_dir)

    for p in planes:
        try:
            res, rows = process_plane_quant(p, args, logger, qc_dir)
            per_plane.append(asdict(res))
            assignment_rows.extend(rows)
        except Exception as e:
            logger.exception("Quantification failed: %s", p.plane_dir)
            per_plane.append({
                "sample": p.sample, "plane": p.plane, "side": p.side, "plane_index": p.plane_index,
                "tile": p.tile, "plane_dir": str(p.plane_dir), "status": f"error:{e}", "missing_required": "error"
            })

    per_plane_df = pd.DataFrame(per_plane)
    if not per_plane_df.empty:
        per_plane_df = attach_sample_metadata(per_plane_df, samplesheet)
    per_plane_df.to_csv(outdir / "cfos_npy_pomc_per_plane_strict_dapi_fos.csv", index=False)

    assign_df = pd.DataFrame(assignment_rows)
    if not assign_df.empty:
        assign_df = attach_sample_metadata(assign_df, samplesheet)
    assign_df.to_csv(outdir / "cfos_object_assignments_strict_dapi_fos.csv", index=False)

    summary_df = aggregate_summary(per_plane_df)
    if not summary_df.empty:
        summary_df = attach_sample_metadata(summary_df, samplesheet)
    summary_df.to_csv(outdir / "cfos_npy_pomc_summary_strict_dapi_fos.csv", index=False)

    # Also write legacy name for old plotting/Excel compatibility.
    summary_df.to_csv(outdir / "cfos_npy_pomc_summary.csv", index=False)

    # Compact c-FOS-only table for quick review.
    if not summary_df.empty:
        fos_cols = [c for c in [
            "sample", "cond", "genotype", "sex", "cage",
            "cfos_cells_left", "cfos_cells_right", "cfos_cells_total",
            "cfos_only_left", "cfos_only_right", "cfos_only_total",
            "cfos_pomc_left", "cfos_pomc_right", "cfos_pomc_total",
            "cfos_npy_left", "cfos_npy_right", "cfos_npy_total",
            "cfos_npy_pomc_left", "cfos_npy_pomc_right", "cfos_npy_pomc_total",
        ] if c in summary_df.columns]
        if fos_cols:
            summary_df[fos_cols].to_csv(outdir / "fos_only_and_fos_subtypes_summary.csv", index=False)
    return per_plane_df, summary_df, assign_df


# -----------------------------------------------------------------------------
# Plots, IQR outlier handling, and statistical tables
# -----------------------------------------------------------------------------

COND_COLORS = {
    "Water": "#b9e3f2",
    "Sucrose": "#e31a1c",
    "Allulose": "#2ecc71",
    "Control": "#999999",
}
COND_ALIASES = {
    "water": "Water", "agua": "Water",
    "sucrose": "Sucrose", "sacarosa": "Sucrose",
    "allulose": "Allulose", "alulosa": "Allulose",
    "control": "Control", "negativectr": "Control", "negativectr": "Control",
}

# metric, y-axis label, title.  These are used consistently by plots,
# outlier reports, and stats tables.
PLOT_SPECS = [
    ("cfos_cells", "c-FOS+ nuclei", "c-FOS+ nuclei by condition"),
    ("cfos_only_cells", "c-FOS-only nuclei", "c-FOS-only nuclei by condition"),
    ("cfos_pomc_cells", "c-FOS∧POMC nuclei", "Activated POMC nuclei by condition"),
    ("cfos_only_over_cfos", "c-FOS-only / c-FOS", "Fraction of c-FOS+ nuclei that are c-FOS-only"),
    ("cfos_pomc_over_cfos", "(c-FOS∧POMC) / c-FOS", "Fraction of c-FOS+ nuclei that are POMC"),
    ("cfos_pomc_over_pomc", "(c-FOS∧POMC) / POMC", "POMC activation ratio"),
]


def get_plot_specs() -> List[Tuple[str, str, str]]:
    return list(PLOT_SPECS)


def clean_condition(c: Any) -> str:
    if c is None or (isinstance(c, float) and math.isnan(c)):
        return "NA"
    s = str(c).strip()
    return COND_ALIASES.get(norm_token(s), s)


def condition_order(vals: Iterable[Any]) -> List[str]:
    conds = [clean_condition(v) for v in vals if str(v) != "nan"]
    seen: List[str] = []
    for c in conds:
        if c not in seen:
            seen.append(c)
    preferred = ["Water", "Sucrose", "Allulose", "Control"]
    out = [c for c in preferred if c in seen]
    out += [c for c in sorted(seen) if c not in out]
    return out


def metric_values(df: Any, metric: str, side: str) -> Any:
    """Return a numeric pandas Series for metric/side, aligned to df rows."""
    if metric == "cfos_cells":
        return pd.to_numeric(df.get(f"cfos_cells_{side}", np.nan), errors="coerce")
    if metric == "cfos_only_cells":
        return pd.to_numeric(df.get(f"cfos_only_{side}", np.nan), errors="coerce")
    if metric == "cfos_pomc_cells":
        return pd.to_numeric(df.get(f"cfos_pomc_{side}", np.nan), errors="coerce")
    if metric == "pomc_cells":
        return pd.to_numeric(df.get(f"pomc_cells_{side}", np.nan), errors="coerce")
    if metric == "cfos_pomc_over_cfos":
        num = pd.to_numeric(df.get(f"cfos_pomc_{side}", np.nan), errors="coerce")
        den = pd.to_numeric(df.get(f"cfos_cells_{side}", np.nan), errors="coerce")
        return num / den.replace(0, np.nan)
    if metric == "cfos_only_over_cfos":
        num = pd.to_numeric(df.get(f"cfos_only_{side}", np.nan), errors="coerce")
        den = pd.to_numeric(df.get(f"cfos_cells_{side}", np.nan), errors="coerce")
        return num / den.replace(0, np.nan)
    if metric == "cfos_pomc_over_pomc":
        num = pd.to_numeric(df.get(f"cfos_pomc_{side}", np.nan), errors="coerce")
        den = pd.to_numeric(df.get(f"pomc_cells_{side}", np.nan), errors="coerce")
        return num / den.replace(0, np.nan)
    raise ValueError(metric)


def safe_anova(groups: List[Any]) -> float:
    clean = [np.asarray(g, dtype=float)[np.isfinite(np.asarray(g, dtype=float))] for g in groups]
    clean = [g for g in clean if g.size >= 2]
    if len(clean) < 2:
        return float("nan")
    try:
        return float(f_oneway(*clean).pvalue)
    except Exception:
        return float("nan")


def cliffs_delta(a: Any, b: Any) -> float:
    """Cliff's delta for b vs a: positive means b tends to be larger than a."""
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size == 0 or y.size == 0:
        return float("nan")
    gt = 0
    lt = 0
    for yy in y:
        gt += int(np.sum(yy > x))
        lt += int(np.sum(yy < x))
    return float((gt - lt) / float(x.size * y.size))


def safe_mwu(a: Any, b: Any, max_partitions: int = 1_000_000) -> float:
    """Exact two-sided label-permutation p for the rank-sum/U statistic."""
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return float("nan")
    pooled = np.concatenate((x, y)); n_first = int(x.size)
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


def safe_wilcoxon_arr(a: Any, b: Any) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return float("nan")
    try:
        if np.allclose(x, y):
            return 1.0
        return float(wilcoxon(x, y, zero_method="wilcox", correction=False).pvalue)
    except Exception:
        return float("nan")


def safe_ttest_rel_arr(a: Any, b: Any) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return float("nan")
    try:
        return float(ttest_rel(x, y, nan_policy="omit").pvalue)
    except Exception:
        return float("nan")


def iqr_bounds(values: Any, k: float = 1.5) -> Tuple[float, float, float, float, float]:
    """Return Q1, Q3, IQR, low fence, high fence. NaNs are ignored."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 4:
        return float("nan"), float("nan"), float("nan"), float("nan"), float("nan")
    q1 = float(np.nanpercentile(v, 25))
    q3 = float(np.nanpercentile(v, 75))
    iqr = q3 - q1
    return q1, q3, iqr, q1 - float(k) * iqr, q3 + float(k) * iqr


def iqr_outlier_mask(values: Any, k: float = 1.5, method: str = "iqr") -> Any:
    """Boolean mask, same length as values. Uses Tukey 1.5*IQR by default."""
    arr = np.asarray(values, dtype=float)
    mask = np.zeros(arr.shape, dtype=bool)
    if str(method).lower() != "iqr":
        return mask
    _, _, _, lo, hi = iqr_bounds(arr, k=float(k))
    if not (np.isfinite(lo) and np.isfinite(hi)):
        return mask
    return (arr < lo) | (arr > hi)


def _prepare_plot_df(summary_df: Any, exclude_control: bool) -> Any:
    df = summary_df.copy()
    if df.empty or "cond" not in df.columns:
        return pd.DataFrame()
    df["cond_clean"] = df["cond"].map(clean_condition)
    if exclude_control:
        df = df[df["cond_clean"] != "Control"].copy()
    return df


def _finite_values_for_group(df: Any, metric: str, side: str) -> Tuple[Any, Any]:
    vals = metric_values(df, metric, side).astype(float).replace([np.inf, -np.inf], np.nan)
    m = vals.notna().to_numpy()
    return vals.to_numpy(dtype=float)[m], df.loc[m].copy()


def _filter_values_for_stats(vals: Any, args: argparse.Namespace) -> Tuple[Any, Any]:
    vals = np.asarray(vals, dtype=float)
    out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
    if bool(getattr(args, "exclude_outliers_stats", False)):
        return vals[~out], out
    return vals, out


def _fmt_p(p: float) -> str:
    return f"{p:.3g}" if np.isfinite(p) else "NA"


def write_outlier_report(summary_df: Any, outdir: Path, args: argparse.Namespace, specs: List[Tuple[str, str, str]], logger: logging.Logger) -> None:
    if summary_df is None or summary_df.empty:
        return
    df = _prepare_plot_df(summary_df, bool(args.exclude_control_plots))
    if df.empty:
        return
    rows: List[Dict[str, Any]] = []
    sides = ["left", "right", "total"]
    for metric, ylabel, title in specs:
        for side in sides:
            for cond in condition_order(df["cond_clean"]):
                sub = df[df["cond_clean"] == cond].copy()
                if sub.empty:
                    continue
                vals = metric_values(sub, metric, side).astype(float).replace([np.inf, -np.inf], np.nan).to_numpy()
                q1, q3, iqr, lo, hi = iqr_bounds(vals, k=float(args.iqr_k))
                out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
                for idx, (_, r) in enumerate(sub.iterrows()):
                    v = vals[idx] if idx < len(vals) else np.nan
                    if not np.isfinite(v):
                        continue
                    rows.append({
                        "metric": metric,
                        "metric_label": ylabel,
                        "side": side,
                        "condition": cond,
                        "sample": r.get("sample", ""),
                        "genotype": r.get("genotype", ""),
                        "sex": r.get("sex", ""),
                        "cage": r.get("cage", ""),
                        "value": float(v),
                        "outlier_method": args.outlier_method,
                        "iqr_k": float(args.iqr_k),
                        "q1": q1,
                        "q3": q3,
                        "iqr": iqr,
                        "lower_fence": lo,
                        "upper_fence": hi,
                        "is_outlier": bool(out[idx]) if idx < len(out) else False,
                    })
    pd.DataFrame(rows).to_csv(outdir / "outlier_report_iqr.csv", index=False)
    logger.warning("IQR outlier report written: %s", outdir / "outlier_report_iqr.csv")


def _group_summary_payload(vals: Any, outmask: Optional[Any] = None) -> Dict[str, Any]:
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    outmask = np.asarray(outmask, dtype=bool) if outmask is not None else np.zeros(vals.shape, dtype=bool)
    return {
        "n": int(vals.size),
        "n_outliers": int(np.count_nonzero(outmask)) if outmask.size == vals.size else 0,
        "mean": float(np.nanmean(vals)) if vals.size else float("nan"),
        "sd": float(np.nanstd(vals, ddof=1)) if vals.size >= 2 else float("nan"),
        "median": float(np.nanmedian(vals)) if vals.size else float("nan"),
    }


def write_stat_tables(summary_df: Any, outdir: Path, args: argparse.Namespace, specs: List[Tuple[str, str, str]], logger: logging.Logger) -> None:
    """Write all-value and IQR-filtered stats tables for condition and paired analyses."""
    if summary_df is None or summary_df.empty:
        return
    df = _prepare_plot_df(summary_df, bool(args.exclude_control_plots))
    if df.empty:
        return

    def one_condition_table(filtered: bool) -> Any:
        rows: List[Dict[str, Any]] = []
        sides = ["left", "right", "total"]
        for metric, ylabel, title in specs:
            for side in sides:
                conds = condition_order(df["cond_clean"])
                group_vals: Dict[str, Any] = {}
                group_out: Dict[str, Any] = {}
                for cond in conds:
                    sub = df[df["cond_clean"] == cond]
                    vals = metric_values(sub, metric, side).astype(float).replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
                    out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
                    group_out[cond] = out
                    group_vals[cond] = vals[~out] if filtered else vals

                p_anova = safe_anova([group_vals[c] for c in conds])
                base = {
                    "metric": metric,
                    "metric_label": ylabel,
                    "side": side,
                    "stats_values": "iqr_filtered" if filtered else "all_values",
                    "outlier_method": args.outlier_method,
                    "iqr_k": float(args.iqr_k),
                    "test": "one_way_anova",
                    "group1": "all_conditions",
                    "group2": "",
                    "p_value": p_anova,
                }
                for cond in conds:
                    payload_all = _group_summary_payload(metric_values(df[df["cond_clean"] == cond], metric, side).astype(float).replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float), group_out[cond])
                    payload_use = _group_summary_payload(group_vals[cond])
                    base[f"n_{cond}"] = payload_use["n"]
                    base[f"mean_{cond}"] = payload_use["mean"]
                    base[f"sd_{cond}"] = payload_use["sd"]
                    base[f"median_{cond}"] = payload_use["median"]
                    base[f"n_outliers_{cond}"] = payload_all["n_outliers"]
                rows.append(dict(base))

                for i in range(len(conds)):
                    for j in range(i + 1, len(conds)):
                        a, b = conds[i], conds[j]
                        va, vb = group_vals[a], group_vals[b]
                        p_mwu = safe_mwu(va, vb)
                        rows.append({
                            "metric": metric,
                            "metric_label": ylabel,
                            "side": side,
                            "stats_values": "iqr_filtered" if filtered else "all_values",
                            "outlier_method": args.outlier_method,
                            "iqr_k": float(args.iqr_k),
                            "test": "mann_whitney_u",
                            "group1": a,
                            "group2": b,
                            "p_value": p_mwu,
                            "n_group1": int(np.asarray(va).size),
                            "n_group2": int(np.asarray(vb).size),
                            "mean_group1": float(np.nanmean(va)) if np.asarray(va).size else float("nan"),
                            "mean_group2": float(np.nanmean(vb)) if np.asarray(vb).size else float("nan"),
                            "sd_group1": float(np.nanstd(va, ddof=1)) if np.asarray(va).size >= 2 else float("nan"),
                            "sd_group2": float(np.nanstd(vb, ddof=1)) if np.asarray(vb).size >= 2 else float("nan"),
                        })
        return pd.DataFrame(rows)

    def one_paired_table(filtered: bool) -> Any:
        rows: List[Dict[str, Any]] = []
        for metric, ylabel, title in specs:
            for cond in condition_order(df["cond_clean"]):
                sub = df[df["cond_clean"] == cond].copy()
                left = metric_values(sub, metric, "left").astype(float).replace([np.inf, -np.inf], np.nan).to_numpy(dtype=float)
                right = metric_values(sub, metric, "right").astype(float).replace([np.inf, -np.inf], np.nan).to_numpy(dtype=float)
                finite = np.isfinite(left) & np.isfinite(right)
                left = left[finite]
                right = right[finite]
                out_left = iqr_outlier_mask(left, k=float(args.iqr_k), method=str(args.outlier_method))
                out_right = iqr_outlier_mask(right, k=float(args.iqr_k), method=str(args.outlier_method))
                keep = ~(out_left | out_right) if filtered else np.ones(left.shape, dtype=bool)
                l_use = left[keep]
                r_use = right[keep]
                rows.append({
                    "metric": metric,
                    "metric_label": ylabel,
                    "condition": cond,
                    "stats_values": "iqr_filtered" if filtered else "all_values",
                    "outlier_method": args.outlier_method,
                    "iqr_k": float(args.iqr_k),
                    "test": "paired_left_right",
                    "n_pairs": int(l_use.size),
                    "n_pairs_removed_as_outlier": int(np.count_nonzero(~keep)) if keep.size else 0,
                    "wilcoxon_p": safe_wilcoxon_arr(l_use, r_use),
                    "paired_ttest_p": safe_ttest_rel_arr(l_use, r_use),
                    "mean_left": float(np.nanmean(l_use)) if l_use.size else float("nan"),
                    "mean_right": float(np.nanmean(r_use)) if r_use.size else float("nan"),
                    "sd_left": float(np.nanstd(l_use, ddof=1)) if l_use.size >= 2 else float("nan"),
                    "sd_right": float(np.nanstd(r_use, ddof=1)) if r_use.size >= 2 else float("nan"),
                })
        return pd.DataFrame(rows)

    cond_all = one_condition_table(filtered=False)
    cond_iqr = one_condition_table(filtered=True)
    pair_all = one_paired_table(filtered=False)
    pair_iqr = one_paired_table(filtered=True)
    cond_all.to_csv(outdir / "stats_condition_tests_all_values.csv", index=False)
    cond_iqr.to_csv(outdir / "stats_condition_tests_iqr_filtered.csv", index=False)
    pair_all.to_csv(outdir / "stats_paired_left_right_all_values.csv", index=False)
    pair_iqr.to_csv(outdir / "stats_paired_left_right_iqr_filtered.csv", index=False)
    logger.warning("Stat tables written: %s", outdir)


def write_outlier_and_stat_tables(summary_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    specs = get_plot_specs()
    if getattr(args, "write_stat_tables", True):
        write_outlier_report(summary_df, outdir, args, specs, logger)
        write_stat_tables(summary_df, outdir, args, specs, logger)


def _plot_legend_handles(conds: List[str], mark_outliers: bool) -> List[Any]:
    from matplotlib.lines import Line2D
    from matplotlib import patheffects as pe
    handles: List[Any] = []
    for c in conds:
        handles.append(Line2D([0], [0], marker="o", linestyle="none",
                              markerfacecolor=COND_COLORS.get(c, "#cccccc"),
                              markeredgecolor="black", markersize=8, label=c))
    if mark_outliers:
        handles.append(Line2D([0], [0], marker="x", linestyle="none", color="black", markersize=9,
                              label="IQR outlier"))
    return handles


def plot_metric_by_condition(summary_df: Any, out_pdf: Path, metric: str, ylabel: str, title: str, args: argparse.Namespace) -> None:
    df = _prepare_plot_df(summary_df, bool(args.exclude_control_plots))
    if df.empty:
        return
    conds = condition_order(df["cond_clean"])
    sides = ["left", "right", "total"]

    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 12})
    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.8), sharey=True)
    fig.suptitle(title, fontsize=18, weight="bold", y=1.05)
    handles = _plot_legend_handles(conds, bool(args.mark_outliers) and args.outlier_method == "iqr")
    if handles:
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.17),
                   ncol=min(len(handles), 6), frameon=False)

    rng = np.random.default_rng(123)
    stats_filtered = bool(args.exclude_outliers_stats)
    for ax, side in zip(axes, sides):
        groups_all: List[Any] = []
        groups_stats: List[Any] = []
        group_outs: List[Any] = []
        means: List[float] = []
        errs: List[float] = []
        for c in conds:
            vals = metric_values(df[df["cond_clean"] == c], metric, side).astype(float).replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
            out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
            vals_stats = vals[~out] if stats_filtered else vals
            groups_all.append(vals)
            groups_stats.append(vals_stats)
            group_outs.append(out)
            means.append(float(np.nanmean(vals)) if vals.size else float("nan"))
            errs.append(float(np.nanstd(vals, ddof=1)) if vals.size >= 2 else 0.0)

        x = np.arange(len(conds))
        colors = [COND_COLORS.get(c, "#cccccc") for c in conds]
        ax.bar(x, means, yerr=errs, capsize=5, color=colors, edgecolor="black", linewidth=1.8, alpha=0.88)
        for i, vals in enumerate(groups_all):
            if vals.size:
                jitter = rng.normal(0, 0.06, vals.size)
                xs = np.full(vals.size, i) + jitter
                ax.scatter(xs, vals, s=60, color=colors[i], edgecolor="black", linewidth=1.1, zorder=3)
                out = group_outs[i]
                if bool(args.mark_outliers) and str(args.outlier_method) == "iqr" and out.size == vals.size and out.any():
                    ax.scatter(xs[out], vals[out], marker="x", s=120, linewidth=2.6, color="black", zorder=4)
                span = ax.get_ylim()[1] - ax.get_ylim()[0]
                label = f"n={vals.size}"
                if out.size == vals.size and np.count_nonzero(out) > 0:
                    label += f"\nout={int(np.count_nonzero(out))}"
                ax.text(i, means[i] + errs[i] + 0.03 * span, label, ha="center", va="bottom", fontsize=9)

        ax.set_xticks(x)
        ax.set_xticklabels(conds, rotation=0)
        ax.set_title(side.upper(), weight="bold")
        ax.grid(axis="y", alpha=0.22)
        ax.set_axisbelow(True)
        if side == "left":
            ax.set_ylabel(ylabel)

        # Add headroom after plotting points and n labels.
        ymin, ymax = ax.get_ylim()
        span = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 1.0
        ax.set_ylim(ymin, ymax + 0.20 * span)

        p_anova = safe_anova(groups_stats)
        stat_lines = ["stats: IQR-filtered" if stats_filtered else "stats: all values"]
        stat_lines.append(f"ANOVA p={_fmt_p(p_anova)}")
        for i in range(len(conds)):
            for j in range(i + 1, len(conds)):
                p = safe_mwu(groups_stats[i], groups_stats[j])
                stat_lines.append(f"{conds[i]} vs {conds[j]} p={_fmt_p(p)}")
        ax.text(0.98, 0.98, "\n".join(stat_lines), transform=ax.transAxes, ha="right", va="top", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_paired_left_right(summary_df: Any, out_pdf: Path, metric: str, ylabel: str, title: str, args: argparse.Namespace) -> None:
    df = _prepare_plot_df(summary_df, bool(args.exclude_control_plots))
    if df.empty:
        return
    conds = condition_order(df["cond_clean"])
    if not conds:
        return
    ncols = min(3, len(conds))
    nrows = int(math.ceil(len(conds) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 5.1 * nrows), sharey=True)
    axes = np.asarray(axes).reshape(-1)
    rng = np.random.default_rng(123)
    fig.suptitle(title, fontsize=18, weight="bold", y=1.04)
    stats_filtered = bool(args.exclude_outliers_stats)

    for ax, c in zip(axes, conds):
        sub = df[df["cond_clean"] == c]
        left = metric_values(sub, metric, "left").astype(float).replace([np.inf, -np.inf], np.nan).to_numpy(dtype=float)
        right = metric_values(sub, metric, "right").astype(float).replace([np.inf, -np.inf], np.nan).to_numpy(dtype=float)
        m = np.isfinite(left) & np.isfinite(right)
        left, right = left[m], right[m]
        out_left = iqr_outlier_mask(left, k=float(args.iqr_k), method=str(args.outlier_method))
        out_right = iqr_outlier_mask(right, k=float(args.iqr_k), method=str(args.outlier_method))
        keep_stats = ~(out_left | out_right) if stats_filtered else np.ones(left.shape, dtype=bool)
        left_stats, right_stats = left[keep_stats], right[keep_stats]

        color = COND_COLORS.get(c, "#cccccc")
        means = [float(np.nanmean(left)) if left.size else float("nan"), float(np.nanmean(right)) if right.size else float("nan")]
        errs = [float(np.nanstd(left, ddof=1)) if left.size >= 2 else 0.0,
                float(np.nanstd(right, ddof=1)) if right.size >= 2 else 0.0]
        ax.bar([0, 1], means, yerr=errs, capsize=5, color=[color, color], edgecolor="black", linewidth=1.8, alpha=0.85)
        for k, (lv, rv) in enumerate(zip(left, right)):
            j = float(rng.normal(0, 0.03))
            x0, x1 = 0 + j, 1 + j
            ax.plot([x0, x1], [lv, rv], color="black", alpha=0.35, lw=1)
            ax.scatter([x0, x1], [lv, rv], s=55, color=color, edgecolor="black", zorder=3)
            if bool(args.mark_outliers) and str(args.outlier_method) == "iqr":
                if k < len(out_left) and bool(out_left[k]):
                    ax.scatter([x0], [lv], marker="x", s=120, linewidth=2.6, color="black", zorder=4)
                if k < len(out_right) and bool(out_right[k]):
                    ax.scatter([x1], [rv], marker="x", s=120, linewidth=2.6, color="black", zorder=4)
        p_w = safe_wilcoxon_arr(left_stats, right_stats)
        p_t = safe_ttest_rel_arr(left_stats, right_stats)
        suffix = "IQR-filtered" if stats_filtered else "all values"
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Left", "Right"])
        ax.set_title(f"{c}\nn={left_stats.size} pairs; {suffix}\nWilcoxon p={_fmt_p(p_w)}; t p={_fmt_p(p_t)}", weight="bold", fontsize=11)
        ax.grid(axis="y", alpha=0.22)
        ax.set_ylabel(ylabel)
    for ax in axes[len(conds):]:
        ax.set_axis_off()
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def make_plots(summary_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    if summary_df is None or summary_df.empty:
        logger.warning("No summary rows; skipping plots.")
        return
    plot_dir = outdir / "plots"
    ensure_dir(plot_dir)
    for metric, ylabel, title in get_plot_specs():
        try:
            plot_metric_by_condition(summary_df, plot_dir / f"{metric}_by_condition.pdf", metric, ylabel, title, args)
            plot_paired_left_right(summary_df, plot_dir / f"{metric}_left_vs_right_paired.pdf", metric, ylabel, title + " — left vs right", args)
        except Exception:
            logger.exception("Plot failed for metric %s", metric)

    # Group summaries.
    try:
        if "cond" in summary_df.columns:
            numeric_cols = [c for c in summary_df.columns if c.endswith("_total") or c.endswith("_left") or c.endswith("_right")]
            group_cols = [c for c in ["cond", "genotype", "sex"] if c in summary_df.columns]
            if group_cols and numeric_cols:
                summary_df.groupby(group_cols, dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_group_summary_cond_genotype_sex.csv")
            group_cols2 = [c for c in ["cond", "genotype"] if c in summary_df.columns]
            if group_cols2 and numeric_cols:
                summary_df.groupby(group_cols2, dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_group_summary_cond_genotype.csv")
            if "cage" in summary_df.columns:
                summary_df.groupby(["cage"], dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_cage_summary.csv")
    except Exception:
        logger.exception("Could not write group summaries")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
# Reconstructed full-section / regional binning analysis
# -----------------------------------------------------------------------------

@dataclass
class SectionPair:
    sample: str
    sample_key: str
    animal_key: str
    animal_id: str
    is_repetition: bool
    repetition_label: str
    section_index: int
    tile: str
    left: Optional[PlaneFiles]
    right: Optional[PlaneFiles]


def strip_repetition_suffix(sample: Any) -> Tuple[str, bool, str]:
    """Return base sample, whether it looked like a repetition, and the suffix label.

    Examples:
      FR6-2(R)  -> FR6-2, True, R
      FR721(R)  -> FR721, True, R
      FR8-3R    -> FR8-3, True, R   (last-resort common typo pattern)
    """
    s = str(sample).strip()
    m = re.search(r"\s*\(([^)]*)\)\s*$", s)
    if m:
        lab = m.group(1).strip()
        base = re.sub(r"\s*\([^)]*\)\s*$", "", s).strip()
        return base, True, lab
    m2 = re.match(r"^(.*?)(?:[_-]?(?:rep|repeat|repetition))([0-9]+)$", s, flags=re.IGNORECASE)
    if m2:
        return m2.group(1).strip("_- "), True, "rep" + m2.group(2)
    # Handle names like FR8-3R conservatively, but avoid changing Negative_ctr.
    m3 = re.match(r"^(FR\s*-?\d+(?:-\d+)?|FR\d+)(R)$", s, flags=re.IGNORECASE)
    if m3:
        return m3.group(1), True, m3.group(2)
    return s, False, ""


def animal_key_from_sample(sample: Any, args: argparse.Namespace) -> str:
    if bool(getattr(args, "keep_repetitions_as_same_animal", True)):
        base, _is_rep, _lab = strip_repetition_suffix(sample)
        return canonical_sample_id(base)
    return canonical_sample_id(sample)


def animal_id_from_sample(sample: Any, samplesheet: Any, args: argparse.Namespace) -> str:
    key = animal_key_from_sample(sample, args)
    if samplesheet is not None and not samplesheet.empty and "sample_key" in samplesheet.columns:
        ss = samplesheet.drop_duplicates("sample_key", keep="first")
        hit = ss[ss["sample_key"] == key]
        if not hit.empty:
            return str(hit.iloc[0]["sample"])
    base, _is_rep, _lab = strip_repetition_suffix(sample) if bool(getattr(args, "keep_repetitions_as_same_animal", True)) else (str(sample), False, "")
    return str(base)


def build_section_pairs(planes: List[PlaneFiles], samplesheet: Any, args: argparse.Namespace) -> List[SectionPair]:
    grouped: Dict[Tuple[str, int, str], Dict[str, Any]] = {}
    for p in planes:
        # Stitch left/right within each source sample folder. Repetitions are only merged later at animal level.
        key = (p.sample, int(p.plane_index), str(p.tile or ""))
        if key not in grouped:
            base, is_rep, rep_label = strip_repetition_suffix(p.sample)
            grouped[key] = {
                "sample": p.sample,
                "sample_key": canonical_sample_id(p.sample),
                "animal_key": animal_key_from_sample(p.sample, args),
                "animal_id": animal_id_from_sample(p.sample, samplesheet, args),
                "is_repetition": is_rep,
                "repetition_label": rep_label,
                "section_index": int(p.plane_index),
                "tile": str(p.tile or ""),
                "left": None,
                "right": None,
            }
        if p.side == "left":
            grouped[key]["left"] = p
        elif p.side == "right":
            grouped[key]["right"] = p
    pairs = [SectionPair(**v) for v in grouped.values()]
    pairs.sort(key=lambda x: (x.animal_key, x.sample, x.section_index, x.tile))
    return pairs


def parse_region_bins(region_bins: str, n_bins: int) -> Dict[int, str]:
    n_bins = int(max(1, n_bins))
    spec = str(region_bins or "default").strip()
    out: Dict[int, str] = {}
    if spec.lower() in {"", "default", "auto"}:
        # This projected-wall fallback maps only ARC/ME; accepted HIL masks can add VMN.
        # y=0 is dorsal/top and y=max is ventral/bottom. ME is the most ventral band.
        for i in range(n_bins):
            f = (i + 0.5) / float(n_bins)
            out[i] = "ME" if f >= 0.84 else "ARC"
        return out

    for block in spec.split(","):
        block = block.strip()
        if not block:
            continue
        if ":" not in block:
            raise ValueError(f"Invalid --region-bins block '{block}'. Expected REGION:0-2 or REGION:0|1|2")
        region, ranges = block.split(":", 1)
        region = region.strip()
        for part in re.split(r"[|;+]", ranges.strip()):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                a, b = part.split("-", 1)
                lo, hi = int(a), int(b)
                if hi < lo:
                    lo, hi = hi, lo
                for k in range(lo, hi + 1):
                    if 0 <= k < n_bins:
                        out[k] = region
            else:
                k = int(part)
                if 0 <= k < n_bins:
                    out[k] = region
    for i in range(n_bins):
        out.setdefault(i, "ARC")
    return out


def region_order_from_args(args: argparse.Namespace, region_map: Dict[int, str]) -> List[str]:
    # Native HIL summaries may contain any delineated subset of ARC/ME/VMN;
    # OUTSIDE_ARC_ME remains excluded unless explicitly requested.
    requested = [x.strip() for x in str(args.region_order or "ARC,ME,VMN").split(",") if x.strip()]
    ordered: List[str] = []
    for r in requested:
        rr = r.upper() if r.lower() in {"arc", "me", "vmn"} else r
        if bool(getattr(args, "arc_me_only", True)) and rr not in {"ARC", "ME", "VMN"}:
            continue
        if rr not in ordered:
            ordered.append(rr)
    for r in ["ARC", "ME", "VMN"]:
        if r not in ordered:
            ordered.append(r)
    if not bool(getattr(args, "arc_me_only", True)) and "OUTSIDE_ARC_ME" not in ordered:
        ordered.append("OUTSIDE_ARC_ME")
    return ordered


def bin_index_from_y(y: float, height: int, n_bins: int) -> int:
    if height <= 0:
        return 0
    idx = int(math.floor(float(y) / (float(height) / float(max(1, n_bins)))))
    return int(max(0, min(n_bins - 1, idx)))


def _empty_label(shape: Tuple[int, int]) -> Any:
    return np.zeros(shape, dtype=np.int32)


def _pad_to_height(arr: Any, height: int) -> Any:
    arr = np.asarray(arr)
    if arr.shape[0] == height:
        return arr
    out_shape = (height, arr.shape[1]) + tuple(arr.shape[2:])
    out = np.zeros(out_shape, dtype=arr.dtype)
    out[:arr.shape[0], :arr.shape[1], ...] = arr
    return out


def _offset_labels(lbl: Any, offset: int) -> Any:
    arr = np.asarray(lbl).astype(np.int32, copy=True)
    if offset > 0:
        m = arr > 0
        arr[m] += int(offset)
    return arr


def stitch_label_pair(left: Optional[Any], right: Optional[Any], gap_px: int = 0, orientation: str = "left-right") -> Any:
    if left is None and right is None:
        return None
    if left is None:
        left = _empty_label(tuple(right.shape))
    if right is None:
        right = _empty_label(tuple(left.shape))
    left = ensure_label_image(left)
    right = ensure_label_image(right)
    h = int(max(left.shape[0], right.shape[0]))
    left = _pad_to_height(left, h)
    right = _pad_to_height(right, h)
    gap_px = int(max(0, gap_px))
    gap = np.zeros((h, gap_px), dtype=np.int32) if gap_px > 0 else None
    if orientation == "right-left":
        first = _offset_labels(right, 0)
        second = _offset_labels(left, int(first.max()))
    else:
        first = _offset_labels(left, 0)
        second = _offset_labels(right, int(first.max()))
    if gap is None:
        return np.concatenate([first, second], axis=1).astype(np.int32)
    return np.concatenate([first, gap, second], axis=1).astype(np.int32)


def _as_gray_image(arr: Any) -> Any:
    if arr is None:
        return None
    a = np.asarray(arr)
    while a.ndim > 3:
        a = a.max(axis=0)
    if a.ndim == 3 and a.shape[-1] in (3, 4):
        a = a[..., :3].max(axis=-1)
    elif a.ndim == 3:
        a = a.max(axis=0)
    return np.asarray(a)


def load_marker_image_from_plane(p: Optional[PlaneFiles], marker: str) -> Optional[Any]:
    if p is None:
        return None
    path = p.images.get(marker)
    if path is None or not path.exists():
        return None
    try:
        return _as_gray_image(read_image(path))
    except Exception as e:
        # The reconstructed raw channel TIFFs are only a QC convenience.
        # Quantification uses label masks already loaded above, so one old/LZW
        # TIFF decoder problem should not stop the full biological analysis.
        print(
            f"[WARN] Could not read raw {marker} image for reconstructed TIFF export: {path} ({e})",
            file=sys.stderr,
        )
        return None


def stitch_image_pair(left: Optional[Any], right: Optional[Any], gap_px: int = 0, orientation: str = "left-right") -> Optional[Any]:
    if left is None and right is None:
        return None
    if left is None:
        left = np.zeros_like(right)
    if right is None:
        right = np.zeros_like(left)
    left = np.asarray(left)
    right = np.asarray(right)
    h = int(max(left.shape[0], right.shape[0]))
    left = _pad_to_height(left, h)
    right = _pad_to_height(right, h)
    gap_px = int(max(0, gap_px))
    gap_shape = (h, gap_px) + tuple(left.shape[2:])
    gap = np.zeros(gap_shape, dtype=left.dtype) if gap_px > 0 else None
    if orientation == "right-left":
        pieces = [right, left]
    else:
        pieces = [left, right]
    if gap is not None:
        pieces = [pieces[0], gap, pieces[1]]
    return np.concatenate(pieces, axis=1)


def load_and_stitch_section_labels(pair: SectionPair, args: argparse.Namespace) -> Tuple[Optional[Any], Optional[Any], Optional[Any], Optional[Any], List[str]]:
    labels: Dict[str, Dict[str, Optional[Any]]] = {"left": {}, "right": {}}
    missing: List[str] = []
    for side_name, p in [("left", pair.left), ("right", pair.right)]:
        if p is None:
            for mk in MARKERS:
                labels[side_name][mk] = None
            continue
        dapi, cfos, npy, pomc, miss = load_plane_labels(p)
        labels[side_name]["DAPI"] = dapi
        labels[side_name]["cFOS"] = cfos
        labels[side_name]["NPY"] = npy
        labels[side_name]["POMC"] = pomc
        for m in miss:
            missing.append(f"{side_name}:{m}")

    # Fill missing markers per half with zero masks matching that half. This allows
    # quantifying the valid half when the other half lacks one segmentation mask.
    for side_name in ["left", "right"]:
        ref_shape = None
        for ref_mk in ["DAPI", "cFOS", "POMC", "NPY"]:
            ref = labels[side_name].get(ref_mk)
            if ref is not None:
                ref_shape = tuple(np.asarray(ref).shape[:2])
                break
        if ref_shape is not None:
            for mk in MARKERS:
                if labels[side_name].get(mk) is None:
                    labels[side_name][mk] = _empty_label(ref_shape)

    dapi_lbl = stitch_label_pair(labels["left"].get("DAPI"), labels["right"].get("DAPI"), int(args.stitch_gap_px), str(args.stitch_orientation))
    cfos_lbl = stitch_label_pair(labels["left"].get("cFOS"), labels["right"].get("cFOS"), int(args.stitch_gap_px), str(args.stitch_orientation))
    npy_lbl = stitch_label_pair(labels["left"].get("NPY"), labels["right"].get("NPY"), int(args.stitch_gap_px), str(args.stitch_orientation))
    pomc_lbl = stitch_label_pair(labels["left"].get("POMC"), labels["right"].get("POMC"), int(args.stitch_gap_px), str(args.stitch_orientation))

    if dapi_lbl is None or cfos_lbl is None:
        return dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl, missing
    # Missing optional channels become empty masks with the reconstructed DAPI shape.
    if npy_lbl is None:
        npy_lbl = _empty_label(tuple(dapi_lbl.shape))
    if pomc_lbl is None:
        pomc_lbl = _empty_label(tuple(dapi_lbl.shape))
    return dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl, missing


def load_and_stitch_section_images(pair: SectionPair, args: argparse.Namespace) -> Dict[str, Optional[Any]]:
    out: Dict[str, Optional[Any]] = {}
    for mk in MARKERS:
        li = load_marker_image_from_plane(pair.left, mk)
        ri = load_marker_image_from_plane(pair.right, mk)
        out[mk] = stitch_image_pair(li, ri, int(args.stitch_gap_px), str(args.stitch_orientation))
    return out


def save_reconstructed_section_tifs(pair: SectionPair, imgs: Dict[str, Optional[Any]], labels: Dict[str, Optional[Any]], outdir: Path, args: argparse.Namespace) -> None:
    if not bool(args.save_reconstructed_tifs) or bool(args.dry_run):
        return
    base = outdir / args.reconstructed_dirname / pair.animal_id / f"{pair.sample}_S{pair.section_index:02d}{('_T'+pair.tile) if pair.tile else ''}"
    ensure_dir(base)
    for mk, arr in imgs.items():
        if arr is not None:
            tifffile.imwrite(str(base / f"reconstructed_{mk}.tif"), np.asarray(arr))
    for mk, arr in labels.items():
        if arr is not None:
            tifffile.imwrite(str(base / f"reconstructed_{mk}_labels.tif"), np.asarray(arr).astype(np.int32))


def classify_phenotype(is_cfos: bool, is_npy: bool, is_pomc: bool) -> str:
    if is_cfos and is_npy and is_pomc:
        return "cFOS_NPY_POMC"
    if is_cfos and is_npy:
        return "cFOS_NPY"
    if is_cfos and is_pomc:
        return "cFOS_POMC"
    if is_cfos:
        return "cFOS_only"
    if is_npy and is_pomc:
        return "NPY_POMC"
    if is_npy:
        return "NPY_only"
    if is_pomc:
        return "POMC_only"
    return "DAPI_only"


def summarize_cell_rows(cell_df: Any, group_cols: List[str], all_groups: Optional[List[Dict[str, Any]]] = None) -> Any:
    """Summarize nuclei while retaining marker-assay availability.

    POMC-specific denominators must only use sections in which a POMC segmentation
    was actually available. Otherwise animals stained for NPY rather than POMC would
    be incorrectly treated as biological zeroes.
    """
    metric_cols = [
        "dapi_nuclei", "cfos_cells", "npy_cells", "pomc_cells",
        "cfos_only", "cfos_npy", "cfos_pomc", "cfos_npy_pomc",
        "pomc_assay_dapi_nuclei", "npy_assay_dapi_nuclei",
    ]
    if cell_df is None or cell_df.empty:
        base = pd.DataFrame(all_groups or [])
        for c in metric_cols:
            base[c] = 0
        add_ratio_columns(base)
        return base

    df = cell_df.copy()
    df["dapi_nuclei"] = 1
    df["cfos_cells"] = df["is_cfos"].astype(int)
    df["npy_cells"] = df["is_npy"].astype(int)
    df["pomc_cells"] = df["is_pomc"].astype(int)
    df["cfos_only"] = ((df["is_cfos"]) & (~df["is_npy"]) & (~df["is_pomc"])).astype(int)
    df["cfos_npy"] = ((df["is_cfos"]) & (df["is_npy"]) & (~df["is_pomc"])).astype(int)
    df["cfos_pomc"] = ((df["is_cfos"]) & (~df["is_npy"]) & (df["is_pomc"])).astype(int)
    df["cfos_npy_pomc"] = ((df["is_cfos"]) & (df["is_npy"]) & (df["is_pomc"])).astype(int)

    has_pomc = df["has_pomc_channel"].astype(bool) if "has_pomc_channel" in df.columns else pd.Series(False, index=df.index)
    has_npy = df["has_npy_channel"].astype(bool) if "has_npy_channel" in df.columns else pd.Series(False, index=df.index)
    df["pomc_assay_dapi_nuclei"] = has_pomc.astype(int)
    df["npy_assay_dapi_nuclei"] = has_npy.astype(int)

    g = df.groupby(group_cols, dropna=False)[metric_cols].sum().reset_index()
    if all_groups:
        full = pd.DataFrame(all_groups)
        g = full.merge(g, on=group_cols, how="left")
        for c in metric_cols:
            g[c] = pd.to_numeric(g[c], errors="coerce").fillna(0).astype(int)
    add_ratio_columns(g)
    return g


def add_ratio_columns(df: Any) -> Any:
    if df is None or df.empty:
        return df

    def ratio(num: str, den: str) -> Any:
        if num not in df.columns or den not in df.columns:
            return pd.Series(np.nan, index=df.index, dtype=float)
        n = pd.to_numeric(df[num], errors="coerce").astype(float)
        d = pd.to_numeric(df[den], errors="coerce").astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = n / d
        return r.replace([np.inf, -np.inf], np.nan)

    # Total activated POMC/NPY include the very rare triple-positive class as well.
    triple_src = df["cfos_npy_pomc"] if "cfos_npy_pomc" in df.columns else pd.Series(0, index=df.index)
    triple = pd.to_numeric(triple_src, errors="coerce").fillna(0)
    if "cfos_pomc" in df.columns:
        df["cfos_pomc_all"] = pd.to_numeric(df["cfos_pomc"], errors="coerce").fillna(0) + triple
    if "cfos_npy" in df.columns:
        df["cfos_npy_all"] = pd.to_numeric(df["cfos_npy"], errors="coerce").fillna(0) + triple

    df["cfos_over_dapi"] = ratio("cfos_cells", "dapi_nuclei")
    df["pomc_over_dapi"] = ratio("pomc_cells", "pomc_assay_dapi_nuclei")
    df["npy_over_dapi"] = ratio("npy_cells", "npy_assay_dapi_nuclei")
    df["cfos_pomc_over_dapi"] = ratio("cfos_pomc_all", "pomc_assay_dapi_nuclei")
    df["cfos_pomc_over_pomc"] = ratio("cfos_pomc_all", "pomc_cells")
    df["cfos_pomc_over_cfos"] = ratio("cfos_pomc_all", "cfos_cells")
    df["cfos_npy_over_dapi"] = ratio("cfos_npy_all", "npy_assay_dapi_nuclei")
    df["cfos_npy_over_npy"] = ratio("cfos_npy_all", "npy_cells")
    df["cfos_npy_over_cfos"] = ratio("cfos_npy_all", "cfos_cells")
    df["fos_only_over_cfos"] = ratio("cfos_only", "cfos_cells")
    return df

def add_area_normalized_columns(
    df: Any,
    *,
    overall_area_col: str = "region_area_um2",
    pomc_area_col: str = "pomc_assay_area_um2",
    npy_area_col: str = "npy_assay_area_um2",
) -> Any:
    """Add exact cells/um^2 columns.

    c-FOS and DAPI use the full analyzed ARC/ME area. POMC and NPY endpoints use
    only the area from sections where the corresponding marker was actually assayed.
    """
    if df is None or df.empty:
        return df

    def numeric_series(col: str) -> Any:
        src = df[col] if col in df.columns else pd.Series(np.nan, index=df.index, dtype=float)
        return pd.to_numeric(src, errors="coerce").astype(float)

    overall_area = numeric_series(overall_area_col)
    pomc_area = numeric_series(pomc_area_col)
    npy_area = numeric_series(npy_area_col)

    def density(num_col: str, area: Any) -> Any:
        num = numeric_series(num_col)
        with np.errstate(divide="ignore", invalid="ignore"):
            out = num / area
        return out.replace([np.inf, -np.inf], np.nan)

    df["dapi_per_um2"] = density("dapi_nuclei", overall_area)
    df["cfos_per_um2"] = density("cfos_cells", overall_area)
    df["cfos_only_per_um2"] = density("cfos_only", overall_area)
    df["pomc_per_um2"] = density("pomc_cells", pomc_area)
    df["cfos_pomc_per_um2"] = density("cfos_pomc_all", pomc_area)
    df["npy_per_um2"] = density("npy_cells", npy_area)
    df["cfos_npy_per_um2"] = density("cfos_npy_all", npy_area)
    return df


# -----------------------------------------------------------------------------
# ARC / ME computer-vision region inference
# -----------------------------------------------------------------------------

REGION_CODE_TO_NAME = {0: "OUTSIDE_ARC_ME", 1: "ARC", 2: "ME", 3: "VMN"}
REGION_NAME_TO_CODE = {v: k for k, v in REGION_CODE_TO_NAME.items()}


def um_to_px(value_um: float, args: argparse.Namespace, *, min_px: int = 0) -> int:
    try:
        um = float(value_um)
        px_size = float(args.um_per_px)
    except Exception:
        return int(max(min_px, round(float(value_um))))
    if px_size <= 0:
        return int(max(min_px, round(um)))
    return int(max(min_px, round(um / px_size)))


def smooth_tissue_mask_from_dapi(dapi_lbl: Any, args: argparse.Namespace) -> Any:
    """Convert sparse DAPI nuclei labels into a smooth tissue mask without filling the third ventricle."""
    mask = np.asarray(dapi_lbl) > 0
    if mask.size == 0:
        return mask
    close_px = um_to_px(float(getattr(args, "region_tissue_close_um", 7.5)), args, min_px=1)
    try:
        if close_px > 0:
            mask = binary_closing(mask, disk(close_px))
            mask = binary_dilation(mask, disk(max(1, close_px // 2)))
        # Remove tiny isolated islands while keeping fragmented ventral ME tissue.
        min_obj = max(64, int((close_px + 1) ** 2))
        mask = remove_small_objects(mask.astype(bool), min_size=min_obj)
        mask = remove_small_holes(mask.astype(bool), area_threshold=max(64, int((close_px + 1) ** 2)))
    except Exception:
        mask = np.asarray(dapi_lbl) > 0
    return mask.astype(bool)


def tissue_bounds(mask: Any) -> Tuple[int, int, int, int]:
    ys, xs = np.where(np.asarray(mask) > 0)
    if ys.size == 0:
        H, W = np.asarray(mask).shape[:2]
        return 0, max(0, H - 1), 0, max(0, W - 1)
    return int(ys.min()), int(ys.max()), int(xs.min()), int(xs.max())


def longest_false_run(row_bool: Any, start_x: int, end_x: int) -> Tuple[int, int, int]:
    """Find longest empty/background run in row_bool[start_x:end_x]. Returns x0,x1,length."""
    row = np.asarray(row_bool).astype(bool)
    start_x = int(max(0, start_x)); end_x = int(min(row.size, end_x))
    if end_x <= start_x:
        return start_x, start_x, 0
    empty = ~row[start_x:end_x]
    best_len = 0; best0 = start_x; best1 = start_x
    cur0 = None
    for i, v in enumerate(empty):
        if bool(v) and cur0 is None:
            cur0 = i
        elif (not bool(v)) and cur0 is not None:
            L = i - cur0
            if L > best_len:
                best_len = L; best0 = start_x + cur0; best1 = start_x + i
            cur0 = None
    if cur0 is not None:
        L = len(empty) - cur0
        if L > best_len:
            best_len = L; best0 = start_x + cur0; best1 = start_x + len(empty)
    return int(best0), int(best1), int(best_len)


def interpolate_nan_vector(v: Any, fallback: float) -> Any:
    arr = np.asarray(v, dtype=float)
    x = np.arange(arr.size)
    good = np.isfinite(arr)
    if good.sum() == 0:
        return np.full(arr.size, float(fallback), dtype=float)
    if good.sum() == 1:
        return np.full(arr.size, float(arr[good][0]), dtype=float)
    return np.interp(x, x[good], arr[good])


def detect_third_ventricle_edges(tissue_mask: Any, args: argparse.Namespace) -> Dict[str, Any]:
    """
    Detect the central third-ventricle gap row-by-row from the smoothed DAPI tissue mask.

    This is deliberately train-free: it learns the midline and ventricle walls from each
    reconstructed section, then ARC/ME are inferred relative to those local landmarks.
    """
    tissue = np.asarray(tissue_mask).astype(bool)
    H, W = tissue.shape
    y0, y1, _x0, _x1 = tissue_bounds(tissue)
    mid_guess = W / 2.0
    half = float(getattr(args, "ventricle_search_half_width_frac", 0.24)) * W
    sx0 = int(max(0, round(mid_guess - half)))
    sx1 = int(min(W, round(mid_guess + half)))
    min_gap_px = um_to_px(float(getattr(args, "ventricle_min_gap_um", 30.0)), args, min_px=5)
    left_edge = np.full(H, np.nan, dtype=float)
    right_edge = np.full(H, np.nan, dtype=float)
    gap_width = np.zeros(H, dtype=float)
    accepted = np.zeros(H, dtype=bool)
    for y in range(H):
        a, b, L = longest_false_run(tissue[y, :], sx0, sx1)
        if L >= min_gap_px:
            left_edge[y] = float(a)
            right_edge[y] = float(b)
            gap_width[y] = float(L)
            accepted[y] = True
    # Use robust central-gap rows within tissue vertical bounds only.
    tissue_rows = np.zeros(H, dtype=bool)
    tissue_rows[y0:y1 + 1] = True
    good = accepted & tissue_rows
    if good.any():
        mids = (left_edge[good] + right_edge[good]) / 2.0
        mid_x = float(np.nanmedian(mids))
        half_gap = float(np.nanmedian((right_edge[good] - left_edge[good]) / 2.0))
        coverage = float(good.sum() / max(1, (y1 - y0 + 1)))
    else:
        mid_x = float(mid_guess)
        half_gap = float(max(um_to_px(35.0, args, min_px=8), W * 0.025))
        coverage = 0.0
    left_edge_i = interpolate_nan_vector(left_edge, mid_x - half_gap)
    right_edge_i = interpolate_nan_vector(right_edge, mid_x + half_gap)
    return {
        "mid_x": mid_x,
        "left_edge": left_edge_i,
        "right_edge": right_edge_i,
        "gap_width": gap_width,
        "accepted_rows": accepted,
        "coverage": coverage,
        "tissue_y_min": int(y0),
        "tissue_y_max": int(y1),
        "search_x0": int(sx0),
        "search_x1": int(sx1),
        "min_gap_px": int(min_gap_px),
    }


def _polygon_to_mask(shape: Tuple[int, int], points: Sequence[Tuple[float, float]]) -> Any:
    """Rasterize a polygon defined as (x, y) points into a boolean mask."""
    H, W = int(shape[0]), int(shape[1])
    if len(points) < 3:
        return np.zeros((H, W), dtype=bool)
    try:
        from matplotlib.path import Path as MplPath
        yy, xx = np.mgrid[0:H, 0:W]
        pts = np.vstack((xx.ravel(), yy.ravel())).T
        mask = MplPath([(float(x), float(y)) for x, y in points]).contains_points(pts)
        return mask.reshape((H, W))
    except Exception:
        return np.zeros((H, W), dtype=bool)


def _safe_edge_at(edge: Any, y: int, fallback: float) -> float:
    arr = np.asarray(edge, dtype=float)
    if arr.size == 0:
        return float(fallback)
    yy = int(max(0, min(arr.size - 1, round(float(y)))))
    val = float(arr[yy])
    return val if np.isfinite(val) else float(fallback)



def _smooth1d_vector(v: Any, sigma: float) -> Any:
    arr = np.asarray(v, dtype=np.float32)
    if arr.size == 0:
        return arr
    try:
        from scipy.ndimage import gaussian_filter1d
        return gaussian_filter1d(arr, sigma=float(max(0.1, sigma)), mode="nearest")
    except Exception:
        return arr


def _normalize_1d(v: Any) -> Any:
    arr = np.asarray(v, dtype=np.float32)
    if arr.size == 0:
        return arr
    lo = float(np.nanpercentile(arr, 2))
    hi = float(np.nanpercentile(arr, 98))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        mx = float(np.nanmax(arr)) if np.isfinite(np.nanmax(arr)) else 1.0
        return np.clip(arr / max(mx, 1e-6), 0, 1).astype(np.float32)
    return np.clip((arr - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32)


def _row_density_from_label_centroids(lbl: Any, H: int, *, weight: float = 1.0, sigma: float = 10.0) -> Any:
    arr = np.asarray(lbl) if lbl is not None else np.zeros((H, 1), dtype=np.int32)
    out = np.zeros(int(H), dtype=np.float32)
    if arr.ndim != 2 or arr.shape[0] != H or int(np.max(arr)) <= 0:
        return out
    try:
        for r in regionprops(arr.astype(np.int32, copy=False)):
            y = int(max(0, min(H - 1, round(float(r.centroid[0])))))
            out[y] += float(weight)
    except Exception:
        # Fallback to pixel density if regionprops fails.
        out = (arr > 0).mean(axis=1).astype(np.float32) * float(weight)
    return _smooth1d_vector(out, sigma=sigma)


def _weighted_quantile(values: Sequence[float], weights: Sequence[float], q: float) -> float:
    v = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    ok = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if ok.sum() == 0:
        return float("nan")
    v = v[ok]
    w = w[ok]
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cw = np.cumsum(w)
    if cw[-1] <= 0:
        return float("nan")
    return float(np.interp(float(q) * cw[-1], cw, v))


def _auto_me_top_y(
    dapi_lbl: Any,
    cfos_lbl: Optional[Any],
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    tissue: Any,
    vent: Dict[str, Any],
    args: argparse.Namespace,
    y_min: int,
    y_max: int,
) -> Tuple[int, Dict[str, Any]]:
    """Estimate the horizontal ARC/ME base line from tissue geometry plus marker activation.

    Primary signal: row with low DAPI/tissue density in the ARC-to-ME transition window.
    Secondary signal: top edge of c-FOS/NPY/POMC activation density. This is deliberately
    weighted and clamped, because activation depends on treatment and should not alone define anatomy.
    """
    H, W = np.asarray(dapi_lbl).shape
    span = max(1, int(y_max) - int(y_min) + 1)
    fixed_y = int(round(int(y_min) + float(getattr(args, "me_top_frac", 0.72)) * span))
    if (str(getattr(args, "arc_me_boundary_mode", "auto")) == "fraction") or (not bool(getattr(args, "auto_me_boundary", True))):
        return int(max(y_min + 1, min(y_max, fixed_y))), {"me_boundary_method": "fixed_me_top_frac", "me_top_y_fixed": int(fixed_y)}

    lo = int(round(int(y_min) + float(getattr(args, "me_boundary_min_frac", 0.55)) * span))
    hi = int(round(int(y_min) + float(getattr(args, "me_boundary_max_frac", 0.84)) * span))
    lo = int(max(y_min + 1, min(y_max - 1, lo)))
    hi = int(max(lo + 3, min(y_max, hi)))

    # DAPI centroid density: the ARC/ME base often sits in a low-density horizontal trough.
    sigma_px = max(3.0, um_to_px(35.0, args, min_px=3))
    dapi_row = _row_density_from_label_centroids(dapi_lbl, H, weight=1.0, sigma=sigma_px)
    dapi_n = _normalize_1d(dapi_row)
    search = dapi_n[lo:hi + 1]
    if search.size > 0 and np.isfinite(search).any():
        dapi_min_y = int(lo + int(np.nanargmin(search)))
    else:
        dapi_min_y = int(fixed_y)

    # Marker/activation density: c-FOS is the strongest guide; NPY/POMC contribute smaller priors.
    cfos_row = _row_density_from_label_centroids(cfos_lbl, H, weight=1.0, sigma=sigma_px) if cfos_lbl is not None else np.zeros(H, dtype=np.float32)
    npy_row = _row_density_from_label_centroids(npy_lbl, H, weight=0.65, sigma=sigma_px) if npy_lbl is not None else np.zeros(H, dtype=np.float32)
    pomc_row = _row_density_from_label_centroids(pomc_lbl, H, weight=0.45, sigma=sigma_px) if pomc_lbl is not None else np.zeros(H, dtype=np.float32)
    marker = _normalize_1d(cfos_row + npy_row + pomc_row)
    marker_edge_y = int(fixed_y)
    marker_ok = False
    if marker[lo:hi + 1].size > 5 and float(np.nanmax(marker[lo:hi + 1])) > 0.05:
        # The top of the ventral activation/NPY cluster is approximated by the strongest
        # positive derivative in the search window. This helps put the ME line just above
        # the dense ventral c-FOS/NPY cluster instead of cutting through it.
        grad = np.gradient(_smooth1d_vector(marker, sigma=max(2.0, sigma_px / 2.0)))
        g = grad[lo:hi + 1]
        if g.size and np.isfinite(g).any():
            marker_edge_y = int(lo + int(np.nanargmax(g)))
            marker_ok = True

    w_marker = float(getattr(args, "region_fos_guide_weight", getattr(args, "me_boundary_fos_weight", 0.35))) if marker_ok else 0.0
    w_marker = max(0.0, min(0.90, w_marker))
    # Fixed anatomical fraction stabilizes sections where marker signal is weak or treatment-dependent.
    y_est = (1.0 - w_marker) * float(dapi_min_y) + w_marker * float(marker_edge_y)
    # Do not allow a wild jump from the anatomical fallback.
    max_jump = max(um_to_px(120.0, args, min_px=50), int(0.12 * span))
    y_est = max(float(fixed_y - max_jump), min(float(fixed_y + max_jump), y_est))
    y_est = int(max(lo, min(hi, round(y_est))))
    return y_est, {
        "me_boundary_method": "auto_dapi_trough_plus_fos_marker_edge",
        "me_top_y_fixed": int(fixed_y),
        "me_top_y_dapi_trough": int(dapi_min_y),
        "me_top_y_marker_edge": int(marker_edge_y),
        "me_boundary_marker_used": bool(marker_ok),
        "me_boundary_fos_weight": float(w_marker),
        "me_boundary_search_lo": int(lo),
        "me_boundary_search_hi": int(hi),
    }


def _auto_arc_widths_from_cells(
    dapi_lbl: Any,
    cfos_lbl: Optional[Any],
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    left_edge: Any,
    right_edge: Any,
    mid_x: float,
    arc_top_y: int,
    me_top_y: int,
    default_bottom_w: int,
    default_top_w: int,
    args: argparse.Namespace,
) -> Dict[str, int]:
    """Estimate lateral triangle widths so the ARC ROI does not miss active cells.

    Distances are measured from the local third-ventricle wall to each DAPI/marker centroid
    within the ARC vertical band. c-FOS centroids are intentionally up-weighted because the
    user wants activation to help guide ROI size, but DAPI remains the main anatomical mass.
    """
    if float(getattr(args, "arc_lateral_quantile", getattr(args, "arc_cell_quantile", 0.985))) <= 0:
        return {
            "left_bottom_w": int(default_bottom_w), "right_bottom_w": int(default_bottom_w),
            "left_top_w": int(default_top_w), "right_top_w": int(default_top_w),
            "arc_auto_width_used": 0,
        }
    H, W = np.asarray(dapi_lbl).shape
    q = float(getattr(args, "arc_lateral_quantile", getattr(args, "arc_cell_quantile", 0.985)))
    q = max(0.50, min(0.999, q))
    margin = um_to_px(float(getattr(args, "arc_safety_margin_um", getattr(args, "arc_cell_margin_um", 45.0))), args, min_px=0)
    le = np.asarray(left_edge, dtype=float)
    re = np.asarray(right_edge, dtype=float)

    left_d: List[float] = []
    left_w: List[float] = []
    right_d: List[float] = []
    right_w: List[float] = []

    def add_label(lbl: Any, wt: float) -> None:
        if lbl is None:
            return
        arr = np.asarray(lbl)
        if arr.ndim != 2 or arr.shape != (H, W) or int(np.max(arr)) <= 0:
            return
        try:
            props = regionprops(arr.astype(np.int32, copy=False))
        except Exception:
            return
        for r in props:
            y = float(r.centroid[0]); x = float(r.centroid[1])
            if y < arc_top_y or y > me_top_y:
                continue
            yy = int(max(0, min(H - 1, round(y))))
            if x < mid_x:
                wall = _safe_edge_at(le, yy, mid_x - max(8, W * 0.03))
                dist = float(wall - x)
                if dist > 0:
                    left_d.append(dist); left_w.append(float(wt))
            else:
                wall = _safe_edge_at(re, yy, mid_x + max(8, W * 0.03))
                dist = float(x - wall)
                if dist > 0:
                    right_d.append(dist); right_w.append(float(wt))

    marker_w = float(getattr(args, "region_marker_guide_weight", 1.0))
    fos_w = float(getattr(args, "region_fos_guide_weight", getattr(args, "me_boundary_fos_weight", 0.35)))
    add_label(dapi_lbl, 1.00)
    # c-FOS is a guide, not the sole anatomical definition. Increasing
    # --region-fos-guide-weight makes active clusters expand the ROI more.
    add_label(cfos_lbl, max(0.0, 1.0 + 4.0 * fos_w))
    add_label(npy_lbl, 1.35 * marker_w)
    add_label(pomc_lbl, 1.20 * marker_w)

    lq = _weighted_quantile(left_d, left_w, q)
    rq = _weighted_quantile(right_d, right_w, q)
    left_bottom = int(max(default_bottom_w, (0 if not np.isfinite(lq) else round(lq)) + margin))
    right_bottom = int(max(default_bottom_w, (0 if not np.isfinite(rq) else round(rq)) + margin))
    # Keep top wider than a point so dorsal ARC cells are not clipped, but still triangular.
    left_top = int(max(default_top_w, round(0.25 * left_bottom)))
    right_top = int(max(default_top_w, round(0.25 * right_bottom)))
    max_w = int(max(50, 0.48 * W))
    return {
        "left_bottom_w": int(min(max_w, left_bottom)),
        "right_bottom_w": int(min(max_w, right_bottom)),
        "left_top_w": int(min(max_w, left_top)),
        "right_top_w": int(min(max_w, right_top)),
        "arc_auto_width_used": 1,
        "arc_left_quantile_distance_px": float(lq) if np.isfinite(lq) else np.nan,
        "arc_right_quantile_distance_px": float(rq) if np.isfinite(rq) else np.nan,
    }


def _ventricle_wall_polygon_side(
    shape: Tuple[int, int],
    edge: Any,
    *,
    side: str,
    mid_x: float,
    y_top: int,
    y_base: int,
    top_width: int,
    bottom_width: int,
    margin: int,
) -> Tuple[Any, List[Tuple[float, float]]]:
    """Build ARC ROI as ventricle-wall-following triangular wedge for one hemisphere."""
    H, W = int(shape[0]), int(shape[1])
    ys = np.linspace(int(y_top), int(y_base), num=max(8, min(64, int(y_base - y_top + 1))))
    inner: List[Tuple[float, float]] = []
    for y in ys:
        yy = int(max(0, min(H - 1, round(float(y)))))
        if side == "left":
            wall = min(_safe_edge_at(edge, yy, mid_x - max(8, W * 0.03)), mid_x - max(8, W * 0.03))
            x = max(0.0, float(wall) - float(margin))
        else:
            wall = max(_safe_edge_at(edge, yy, mid_x + max(8, W * 0.03)), mid_x + max(8, W * 0.03))
            x = min(float(W - 1), float(wall) + float(margin))
        inner.append((float(x), float(yy)))
    inner_top = inner[0]
    inner_base = inner[-1]
    if side == "left":
        outer_top = (max(0.0, inner_top[0] - float(top_width)), float(y_top))
        outer_base = (max(0.0, inner_base[0] - float(bottom_width)), float(y_base))
        pts = inner + [outer_base, outer_top]
    else:
        outer_top = (min(float(W - 1), inner_top[0] + float(top_width)), float(y_top))
        outer_base = (min(float(W - 1), inner_base[0] + float(bottom_width)), float(y_base))
        pts = inner + [outer_base, outer_top]
    return _polygon_to_mask((H, W), pts), pts




def _arc_fos_weight(args: argparse.Namespace) -> float:
    """c-FOS weight for anatomical ARC shape inference, independent of c-FOS/DAPI positivity calling."""
    val = getattr(args, "arc_fos_guide_weight", None)
    if val is not None:
        try:
            return float(val)
        except Exception:
            pass
    bw = float(getattr(args, "me_boundary_fos_weight", 0.35))
    return float(max(0.0, 1.0 + 4.0 * bw))


def _collect_arc_shape_points(
    dapi_lbl: Any,
    cfos_lbl: Optional[Any],
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    left_edge: Any,
    right_edge: Any,
    mid_x: float,
    arc_top_y: int,
    me_top_y: int,
    args: argparse.Namespace,
) -> Dict[str, List[Tuple[float, float, float, float]]]:
    """Collect weighted centroids for ARC shape fitting.

    Returns dict side -> records of (y, x, distance_from_local_ventricle_wall, weight).
    c-FOS receives a user-controllable weight because activated neurons can indicate
    the relevant ARC territory in ventral/caudal sections, while a low DAPI weight
    prevents every nucleus in the field from over-expanding the ARC.
    """
    dapi = np.asarray(dapi_lbl)
    H, W = dapi.shape
    le = np.asarray(left_edge, dtype=float)
    re = np.asarray(right_edge, dtype=float)
    records: Dict[str, List[Tuple[float, float, float, float]]] = {"left": [], "right": []}
    ylo = int(max(0, min(H - 1, arc_top_y)))
    yhi = int(max(ylo + 1, min(H - 1, me_top_y)))
    max_dist = max(um_to_px(float(getattr(args, "arc_bottom_width_um", 650.0)) * 1.8, args, min_px=120), int(0.55 * W))

    def add(lbl: Any, wt: float, name: str) -> None:
        if lbl is None or wt <= 0:
            return
        arr = np.asarray(lbl)
        if arr.ndim != 2 or arr.shape != (H, W) or int(np.max(arr)) <= 0:
            return
        try:
            props = regionprops(arr.astype(np.int32, copy=False))
        except Exception:
            return
        for rp in props:
            y, x = float(rp.centroid[0]), float(rp.centroid[1])
            if y < ylo or y > yhi:
                continue
            yy = int(max(0, min(H - 1, round(y))))
            if x < mid_x:
                wall = _safe_edge_at(le, yy, mid_x - max(8, W * 0.03))
                dist = float(wall - x)
                if 0 < dist < max_dist:
                    aw = float(wt) * min(3.0, max(0.75, math.sqrt(max(1.0, float(rp.area))) / 7.0))
                    records["left"].append((float(y), float(x), float(dist), float(aw)))
            else:
                wall = _safe_edge_at(re, yy, mid_x + max(8, W * 0.03))
                dist = float(x - wall)
                if 0 < dist < max_dist:
                    aw = float(wt) * min(3.0, max(0.75, math.sqrt(max(1.0, float(rp.area))) / 7.0))
                    records["right"].append((float(y), float(x), float(dist), float(aw)))

    dapi_w = float(getattr(args, "arc_dapi_guide_weight", 0.35))
    marker_w = float(getattr(args, "arc_marker_guide_weight", getattr(args, "region_marker_guide_weight", 1.25)))
    fos_w = _arc_fos_weight(args)
    add(dapi_lbl, dapi_w, "DAPI")
    add(cfos_lbl, fos_w, "cFOS")
    add(npy_lbl, marker_w, "NPY")
    add(pomc_lbl, marker_w, "POMC")
    return records


def _weighted_median(values: Sequence[float], weights: Sequence[float]) -> float:
    return _weighted_quantile(values, weights, 0.50)


def _ventricle_wall_ellipse_side(
    shape: Tuple[int, int],
    edge: Any,
    *,
    side: str,
    mid_x: float,
    records: List[Tuple[float, float, float, float]],
    y_top: int,
    y_base: int,
    default_rx: int,
    default_ry: int,
    args: argparse.Namespace,
) -> Tuple[Any, Dict[str, float]]:
    """Build a ventricle-wall-referenced oval ARC support for ventral/caudal sections."""
    H, W = int(shape[0]), int(shape[1])
    if not records:
        return np.zeros((H, W), dtype=bool), {
            f"{side}_oval_used": 0,
            f"{side}_oval_center_y": float("nan"),
            f"{side}_oval_center_dist": float("nan"),
            f"{side}_oval_rx": float(default_rx),
            f"{side}_oval_ry": float(default_ry),
        }
    ys = np.asarray([r[0] for r in records], dtype=float)
    ds = np.asarray([r[2] for r in records], dtype=float)
    ws = np.asarray([r[3] for r in records], dtype=float)
    q = float(getattr(args, "arc_oval_quantile", 0.92))
    q = max(0.60, min(0.995, q))
    cy = _weighted_median(ys, ws)
    cd = _weighted_median(ds, ws)
    if not np.isfinite(cy):
        cy = (float(y_top) + float(y_base)) / 2.0
    if not np.isfinite(cd):
        cd = max(20.0, float(default_rx) * 0.5)
    margin = um_to_px(float(getattr(args, "arc_oval_margin_um", 55.0)), args, min_px=0)
    min_rx = um_to_px(float(getattr(args, "arc_oval_min_rx_um", 190.0)), args, min_px=35)
    min_ry = um_to_px(float(getattr(args, "arc_oval_min_ry_um", 120.0)), args, min_px=35)
    rx = _weighted_quantile(np.abs(ds - cd), ws, q)
    ry = _weighted_quantile(np.abs(ys - cy), ws, q)
    rx = int(max(min_rx, default_rx * 0.45, (0 if not np.isfinite(rx) else round(rx)) + margin))
    ry = int(max(min_ry, default_ry * 0.45, (0 if not np.isfinite(ry) else round(ry)) + margin))
    rx = int(min(rx, max(40, round(float(getattr(args, "arc_oval_max_rx_frac", 0.48)) * W))))
    tissue_span = max(1, int(y_base) - int(y_top) + 1)
    ry = int(min(ry, max(30, round(float(getattr(args, "arc_oval_max_ry_frac", 0.48)) * tissue_span))))
    cy = float(max(y_top, min(y_base, cy)))
    yy, xx = np.indices((H, W))
    edge_arr = np.asarray(edge, dtype=float)
    wall = np.zeros(H, dtype=float)
    for y in range(H):
        if side == "left":
            wall[y] = min(_safe_edge_at(edge_arr, y, mid_x - max(8, W * 0.03)), mid_x - max(8, W * 0.03))
        else:
            wall[y] = max(_safe_edge_at(edge_arr, y, mid_x + max(8, W * 0.03)), mid_x + max(8, W * 0.03))
    wall2 = wall[:, None]
    if side == "left":
        dist = wall2 - xx.astype(float)
        side_mask = xx < wall2
    else:
        dist = xx.astype(float) - wall2
        side_mask = xx > wall2
    oval = (((dist - cd) / max(1.0, float(rx))) ** 2 + ((yy.astype(float) - cy) / max(1.0, float(ry))) ** 2) <= 1.0
    oval &= side_mask
    oval &= (yy >= int(y_top)) & (yy <= int(y_base))
    oval &= dist > 0
    try:
        oval = binary_closing(oval.astype(bool), disk(max(1, um_to_px(10.0, args, min_px=1))))
        oval = remove_small_holes(oval.astype(bool), area_threshold=max(64, int(0.0008 * H * W)))
        oval = remove_small_objects(oval.astype(bool), min_size=max(64, int(0.0005 * H * W)))
    except Exception:
        pass
    return oval.astype(bool), {
        f"{side}_oval_used": 1,
        f"{side}_oval_center_y": float(cy),
        f"{side}_oval_center_dist": float(cd),
        f"{side}_oval_rx": float(rx),
        f"{side}_oval_ry": float(ry),
        f"{side}_oval_points": int(len(records)),
    }


def _mask_weighted_point_coverage(mask: Any, records_by_side: Dict[str, List[Tuple[float, float, float, float]]]) -> float:
    m = np.asarray(mask).astype(bool)
    H, W = m.shape
    total = 0.0
    hit = 0.0
    for records in records_by_side.values():
        for y, x, _d, w in records:
            yy = int(max(0, min(H - 1, round(y))))
            xx = int(max(0, min(W - 1, round(x))))
            ww = float(max(0.0, w))
            total += ww
            if m[yy, xx]:
                hit += ww
    if total <= 0:
        return float("nan")
    return float(hit / total)


def _activation_cluster_arc_mask(
    shape: Tuple[int, int],
    cfos_lbl: Optional[Any],
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    dapi_lbl: Any,
    left_edge: Any,
    right_edge: Any,
    mid_x: float,
    arc_top_y: int,
    me_top_y: int,
    args: argparse.Namespace,
) -> Any:
    """Dense c-FOS/marker cluster support used to prevent missing activated ARC cells."""
    H, W = int(shape[0]), int(shape[1])
    if not bool(getattr(args, "arc_include_fos_cluster_mask", True)):
        return np.zeros((H, W), dtype=bool)
    sigma = max(2.0, um_to_px(25.0, args, min_px=2))
    density = np.zeros((H, W), dtype=np.float32)
    try:
        from scipy.ndimage import gaussian_filter
        if cfos_lbl is not None:
            density += float(_arc_fos_weight(args)) * gaussian_filter((np.asarray(cfos_lbl) > 0).astype(np.float32), sigma=sigma)
        if npy_lbl is not None:
            density += float(getattr(args, "arc_marker_guide_weight", 1.25)) * gaussian_filter((np.asarray(npy_lbl) > 0).astype(np.float32), sigma=sigma)
        if pomc_lbl is not None:
            density += float(getattr(args, "arc_marker_guide_weight", 1.25)) * gaussian_filter((np.asarray(pomc_lbl) > 0).astype(np.float32), sigma=sigma)
        density += 0.10 * gaussian_filter((np.asarray(dapi_lbl) > 0).astype(np.float32), sigma=sigma)
    except Exception:
        if cfos_lbl is not None:
            density += (np.asarray(cfos_lbl) > 0).astype(np.float32)
        if npy_lbl is not None:
            density += (np.asarray(npy_lbl) > 0).astype(np.float32)
        if pomc_lbl is not None:
            density += (np.asarray(pomc_lbl) > 0).astype(np.float32)
    yy, xx = np.indices((H, W))
    le = np.asarray(left_edge, dtype=float)
    re = np.asarray(right_edge, dtype=float)
    dist_ok = np.zeros((H, W), dtype=bool)
    max_dist = max(um_to_px(float(getattr(args, "arc_bottom_width_um", 650.0)) * 1.35, args, min_px=160), int(0.35 * W))
    for y in range(max(0, arc_top_y), min(H, me_top_y + 1)):
        lw = _safe_edge_at(le, y, mid_x - max(8, W * 0.03))
        rw = _safe_edge_at(re, y, mid_x + max(8, W * 0.03))
        dist_ok[y, :] = (((lw - xx[y, :]) > 0) & ((lw - xx[y, :]) < max_dist)) | (((xx[y, :] - rw) > 0) & ((xx[y, :] - rw) < max_dist))
    support = dist_ok & (yy >= int(arc_top_y)) & (yy <= int(me_top_y))
    vals = density[support]
    if vals.size == 0 or float(np.nanmax(vals)) <= 0:
        return np.zeros((H, W), dtype=bool)
    q = float(getattr(args, "arc_fos_cluster_quantile", 0.82))
    q = max(0.50, min(0.99, q))
    pos = vals[vals > 0]
    thr = float(np.nanquantile(pos, q)) if pos.size else float(np.nanmax(vals))
    cluster = (density >= thr) & support
    dil = um_to_px(float(getattr(args, "arc_fos_cluster_dilate_um", 25.0)), args, min_px=1)
    try:
        cluster = binary_dilation(cluster.astype(bool), disk(max(1, dil)))
        cluster = binary_closing(cluster.astype(bool), disk(max(1, dil // 2)))
        cluster = remove_small_objects(cluster.astype(bool), min_size=max(20, int(0.0002 * H * W)))
        cluster = remove_small_holes(cluster.astype(bool), area_threshold=max(20, int(0.0002 * H * W)))
    except Exception:
        pass
    return (cluster & support).astype(bool)


def _first_persistent_true_run(flags: Any, offset: int, min_run: int) -> Optional[int]:
    """Return the first row index where a True run of at least min_run starts."""
    arr = np.asarray(flags).astype(bool)
    if arr.size == 0:
        return None
    min_run = int(max(1, min_run))
    run_start = None
    run_len = 0
    for i, v in enumerate(arr):
        if bool(v):
            if run_start is None:
                run_start = i
                run_len = 1
            else:
                run_len += 1
            if run_len >= min_run:
                return int(offset + run_start)
        else:
            run_start = None
            run_len = 0
    return None


def _centroid_records_for_floor(lbl: Any, H: int, W: int, *, weight: float = 1.0) -> List[Tuple[float, float, float, int]]:
    """Return (y, x, weight, area) for label centroids."""
    arr = np.asarray(lbl) if lbl is not None else np.zeros((H, W), dtype=np.int32)
    if arr.ndim != 2 or arr.shape != (H, W) or int(np.max(arr)) <= 0:
        return []
    out: List[Tuple[float, float, float, int]] = []
    try:
        for rp in regionprops(arr.astype(np.int32, copy=False)):
            out.append((float(rp.centroid[0]), float(rp.centroid[1]), float(weight), int(rp.area)))
    except Exception:
        pass
    return out


def _row_density_in_x_window_from_labels(lbl: Any, H: int, W: int, x0: int, x1: int, *, sigma: float) -> Any:
    """Smoothed per-row centroid count inside x0:x1."""
    out = np.zeros(H, dtype=np.float32)
    arr = np.asarray(lbl) if lbl is not None else np.zeros((H, W), dtype=np.int32)
    if arr.ndim != 2 or arr.shape != (H, W) or int(np.max(arr)) <= 0 or x1 <= x0:
        return out
    try:
        for rp in regionprops(arr.astype(np.int32, copy=False)):
            y, x = float(rp.centroid[0]), float(rp.centroid[1])
            if x0 <= x <= x1:
                yy = int(max(0, min(H - 1, round(y))))
                out[yy] += 1.0
    except Exception:
        # fallback to pixel occupancy if regionprops fails
        out = (arr[:, x0:x1] > 0).sum(axis=1).astype(np.float32)
    return _smooth1d_vector(out, sigma=max(1.0, float(sigma)))



def _region_debug(args: argparse.Namespace, msg: str, *vals: Any) -> None:
    """Verbose logger for anatomical ARC/ME inference."""
    try:
        if not bool(getattr(args, "region_dl_verbose", True)) and not bool(getattr(args, "verbose", False)):
            return
        logging.getLogger("keyence20x_pipeline").info("[ARC/ME] " + str(msg), *vals)
    except Exception:
        pass


def _connected_label_image(mask: Any) -> Any:
    """Connected components with a scipy fallback."""
    m = np.asarray(mask).astype(bool)
    try:
        return label(m, connectivity=1).astype(np.int32)
    except Exception:
        try:
            from scipy import ndimage as ndi
            lbl, _n = ndi.label(m)
            return lbl.astype(np.int32)
        except Exception:
            return m.astype(np.int32)


def _density_map_from_label_centroids(lbl: Any, sigma_px: float) -> Any:
    """Raw centroid-density map. Dense ARC has high values; sparse ME has low values."""
    arr = np.asarray(lbl)
    H, W = arr.shape
    imp = np.zeros((H, W), dtype=np.float32)
    try:
        for rp in regionprops(arr.astype(np.int32)):
            y, x = rp.centroid
            yy = int(max(0, min(H - 1, round(float(y)))))
            xx = int(max(0, min(W - 1, round(float(x)))))
            # One event per nucleus.  A weak size term stabilizes tiny/large manual masks.
            imp[yy, xx] += float(min(1.8, max(0.6, math.sqrt(max(1.0, float(rp.area))) / 9.0)))
    except Exception:
        imp[arr > 0] = 1.0
    try:
        from scipy.ndimage import gaussian_filter
        return gaussian_filter(imp, sigma=max(1.0, float(sigma_px))).astype(np.float32)
    except Exception:
        return imp.astype(np.float32)


def _component_bbox_and_stats(mask: Any, density: Any, mid_x: float, floor_top: int, args: argparse.Namespace) -> Dict[str, Any]:
    ys, xs = np.where(np.asarray(mask).astype(bool))
    if ys.size == 0:
        return {"area": 0, "x0": 0, "x1": 0, "y0": 0, "y1": 0, "span_px": 0, "depth_px": 0, "crosses_midline": False, "sparse_fraction": 0.0, "density_median": 0.0}
    min_span_px = um_to_px(float(getattr(args, "me_floor_min_cross_midline_span_um", 120.0)), args, min_px=8)
    left_ok = bool(np.any(xs <= float(mid_x) - min_span_px / 2.0))
    right_ok = bool(np.any(xs >= float(mid_x) + min_span_px / 2.0))
    den = np.asarray(density)
    dens_vals = den[ys, xs] if den.shape == mask.shape else np.zeros_like(ys, dtype=float)
    return {
        "area": int(ys.size),
        "x0": int(xs.min()),
        "x1": int(xs.max()),
        "y0": int(ys.min()),
        "y1": int(ys.max()),
        "span_px": int(xs.max() - xs.min() + 1),
        "depth_px": int(max(0, ys.max() - int(floor_top) + 1)),
        "crosses_midline": bool(left_ok and right_ok),
        "sparse_fraction": float(np.nan),
        "density_median": float(np.nanmedian(dens_vals)) if dens_vals.size else 0.0,
    }


def _dynamic_me_floor_mask(
    dapi_lbl: Any,
    cfos_lbl: Optional[Any],
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    tissue: Any,
    vent: Dict[str, Any],
    arc_shape_mask: Any,
    preliminary_me_top_y: int,
    args: argparse.Namespace,
    y_min: int,
    y_max: int,
    x_min: int,
    x_max: int,
    mid_x: float,
) -> Tuple[Any, Dict[str, Any]]:
    """Detect a true sparse, hanging ME floor and avoid horizontal ME overcalling.

    Biological rules implemented in v16
    -----------------------------------
    1. If there is no tissue hanging below the third-ventricle floor, ME is absent.
    2. If the ventral tissue is only on one side and is not connected across/around the
       midline floor, ME is absent.
    3. Dense lateral arcuate tissue remains ARC.  ME is the sparser hanging floor tissue.
    4. c-FOS/NPY/POMC can help shape support, but c-FOS alone cannot create ME.
    """
    dapi = np.asarray(dapi_lbl)
    H, W = dapi.shape
    tissue = np.asarray(tissue).astype(bool)
    arc_shape_mask = np.asarray(arc_shape_mask).astype(bool)
    yy, xx = np.indices((H, W))
    span = max(1, int(y_max) - int(y_min) + 1)
    mode = str(getattr(args, "me_detection_mode", "auto")).lower()

    def _absent(method: str, extra: Optional[Dict[str, Any]] = None) -> Tuple[Any, Dict[str, Any]]:
        info = {
            "me_floor_present": False,
            "me_floor_detection_mode": mode,
            "me_floor_method": method,
            "me_top_y": None,
            "me_floor_top_y": None,
            "me_floor_preliminary_y": int(preliminary_me_top_y),
            "me_floor_component_count": 0,
            "me_floor_good_component_count": 0,
            "me_floor_bridge_connected": False,
            "me_floor_cross_midline": False,
            "me_floor_sparse_density_ref": float("nan"),
            "me_floor_sparse_threshold": float("nan"),
            "me_floor_sparse_fraction_best": 0.0,
        }
        if extra:
            info.update(extra)
        _region_debug(args, "ME absent: %s | %s", method, {k: info.get(k) for k in ["me_floor_bridge_nuclei", "me_floor_core_nuclei", "me_floor_total_nuclei", "me_floor_component_count", "me_floor_good_component_count"]})
        return np.zeros((H, W), dtype=bool), info

    if mode == "never":
        return _absent("disabled_by_user")

    # Ventricle landmarks.  The floor must be related to the central third ventricle,
    # not to arbitrary ventral/lateral tissue.
    left_edge = np.asarray(vent.get("left_edge", np.full(H, float(mid_x) - 40.0)), dtype=float)
    right_edge = np.asarray(vent.get("right_edge", np.full(H, float(mid_x) + 40.0)), dtype=float)
    accepted_rows = np.asarray(vent.get("accepted_rows", np.zeros(H, dtype=bool))).astype(bool)
    gap_width = np.asarray(vent.get("gap_width", np.zeros(H, dtype=float)), dtype=float)
    vent_good = accepted_rows & (np.arange(H) >= int(y_min)) & (np.arange(H) <= int(y_max))
    if np.any(vent_good):
        vent_bottom_y = int(np.nanpercentile(np.where(vent_good)[0], 98))
    else:
        vent_bottom_y = int(preliminary_me_top_y)

    core_half = um_to_px(float(getattr(args, "me_floor_core_half_width_um", 160.0)), args, min_px=max(10, int(0.04 * W)))
    core_half = int(min(core_half, max(12, int(0.14 * W))))
    bridge_half = um_to_px(float(getattr(args, "me_floor_bridge_half_width_um", 95.0)), args, min_px=max(5, int(0.02 * W)))
    bridge_half = int(min(bridge_half, max(8, int(0.10 * W))))
    lat_half = um_to_px(float(getattr(args, "me_floor_lateral_width_um", 850.0)), args, min_px=max(30, int(0.20 * W)))
    lat_half = int(min(lat_half, max(40, int(0.49 * W))))
    core_x0 = int(max(0, round(float(mid_x) - core_half)))
    core_x1 = int(min(W, round(float(mid_x) + core_half)))
    bridge_x0 = int(max(0, round(float(mid_x) - bridge_half)))
    bridge_x1 = int(min(W, round(float(mid_x) + bridge_half)))
    lat_x0 = int(max(0, round(float(mid_x) - lat_half)))
    lat_x1 = int(min(W, round(float(mid_x) + lat_half)))

    lo_frac = float(getattr(args, "me_floor_search_min_frac", 0.58))
    hi_frac = float(getattr(args, "me_floor_search_max_frac", 0.98))
    lo = int(round(int(y_min) + lo_frac * span))
    hi = int(round(int(y_min) + hi_frac * span))
    top_margin = um_to_px(float(getattr(args, "me_floor_top_margin_um", 25.0)), args, min_px=4)
    # The search is allowed to start slightly above the ventricle bottom to capture large ME,
    # but a still-open ventricle at the image floor will fail the connected bridge test later.
    lo = int(max(y_min + 1, min(y_max - 2, min(lo, int(preliminary_me_top_y) + top_margin, vent_bottom_y + top_margin))))
    hi = int(max(lo + 3, min(y_max, hi)))

    sigma_px = max(2.0, um_to_px(28.0, args, min_px=3))
    dapi_bin = dapi > 0
    # Use the *narrow bridge* first.  A wide core can accidentally see dense ARC walls.
    bridge_rows = _row_density_in_x_window_from_labels(dapi, H, W, bridge_x0, bridge_x1, sigma=sigma_px)
    bridge_pix = dapi_bin[:, bridge_x0:bridge_x1].mean(axis=1).astype(np.float32) if bridge_x1 > bridge_x0 else np.zeros(H, dtype=np.float32)
    bridge_pix = _smooth1d_vector(bridge_pix, sigma=max(1.0, sigma_px / 2.0))
    bridge_score = _normalize_1d(bridge_rows) * 0.70 + _normalize_1d(bridge_pix) * 0.30
    search = bridge_score[lo:hi + 1]
    max_bridge_density = float(np.nanmax(bridge_pix[lo:hi + 1])) if search.size else 0.0
    if search.size == 0 or not np.isfinite(search).any():
        floor_top = None
    else:
        rel = float(getattr(args, "me_floor_threshold_rel", 0.22))
        thr = max(0.025, float(np.nanmax(search)) * max(0.02, min(0.90, rel)))
        min_run = um_to_px(float(getattr(args, "me_floor_min_run_um", 18.0)), args, min_px=3)
        floor_top = _first_persistent_true_run(search >= thr, lo, min_run)

    if floor_top is None:
        if mode == "always":
            floor_top = int(preliminary_me_top_y)
        else:
            return _absent("no_persistent_midline_floor_bridge_density", {
                "me_floor_search_lo": int(lo), "me_floor_search_hi": int(hi),
                "me_floor_max_bridge_density": float(max_bridge_density),
                "me_floor_vent_bottom_y": int(vent_bottom_y),
            })

    # Count actual DAPI nuclei below the putative floor top.
    y_seed0 = int(max(y_min, min(y_max, int(floor_top) - top_margin)))
    anchor_depth = um_to_px(float(getattr(args, "me_floor_anchor_depth_um", 95.0)), args, min_px=15)
    core_nuclei = bridge_nuclei = total_nuclei = 0
    for y, x, _w, _area in _centroid_records_for_floor(dapi, H, W, weight=1.0):
        if y < y_seed0 or y > y_max:
            continue
        if bridge_x0 <= x <= bridge_x1 and y <= y_seed0 + anchor_depth:
            bridge_nuclei += 1
        if core_x0 <= x <= core_x1:
            core_nuclei += 1
        if lat_x0 <= x <= lat_x1:
            total_nuclei += 1

    min_core = int(getattr(args, "me_floor_min_core_nuclei", 10))
    min_bridge = int(getattr(args, "me_floor_min_bridge_nuclei", 3))
    min_total = int(getattr(args, "me_floor_min_total_nuclei", 45))
    min_density = float(getattr(args, "me_floor_min_core_density", 0.003))
    if mode != "always" and not (bridge_nuclei >= min_bridge and core_nuclei >= min_core and total_nuclei >= min_total and max_bridge_density >= min_density):
        return _absent("not_enough_true_central_floor_nuclei", {
            "me_floor_top_y_candidate": int(y_seed0),
            "me_floor_search_lo": int(lo), "me_floor_search_hi": int(hi),
            "me_floor_bridge_nuclei": int(bridge_nuclei),
            "me_floor_core_nuclei": int(core_nuclei),
            "me_floor_total_nuclei": int(total_nuclei),
            "me_floor_max_bridge_density": float(max_bridge_density),
            "me_floor_min_bridge_nuclei": int(min_bridge),
            "me_floor_min_core_nuclei": int(min_core),
            "me_floor_min_total_nuclei": int(min_total),
            "me_floor_min_core_density": float(min_density),
        })

    # Sparse-vs-dense model: ARC nuclei form tight/high-density fields; ME floor is sparser.
    dens_sigma = um_to_px(float(getattr(args, "me_floor_density_sigma_um", 24.0)), args, min_px=3)
    density = _density_map_from_label_centroids(dapi, sigma_px=dens_sigma)
    ref_band = um_to_px(float(getattr(args, "me_floor_arc_reference_band_um", 150.0)), args, min_px=25)
    arc_ref = arc_shape_mask.astype(bool) & (yy >= max(y_min, y_seed0 - ref_band)) & (yy <= min(y_max, y_seed0 + max(4, ref_band // 3)))
    # Exclude the narrow midline bridge from the ARC reference; we want dense ARC shoulder/lobe density.
    arc_ref &= ~((xx >= bridge_x0) & (xx <= bridge_x1) & (yy >= y_seed0 - ref_band))
    ref_vals = density[arc_ref & (density > 0)]
    if ref_vals.size < 30:
        dorsal_ref = (dapi > 0) & (yy >= max(y_min, y_seed0 - ref_band)) & (yy <= min(y_max, y_seed0 + ref_band // 2)) & (xx >= x_min) & (xx <= x_max)
        ref_vals = density[dorsal_ref & (density > 0)]
    if ref_vals.size == 0:
        ref_density = float(np.nanpercentile(density[density > 0], 75)) if np.any(density > 0) else 1.0
    else:
        ref_density = float(np.nanpercentile(ref_vals, 75))
    ratio = float(getattr(args, "me_floor_max_dense_arc_ratio", 0.62))
    sparse_thr = max(1e-6, ref_density * max(0.05, min(1.25, ratio)))
    sparse_mask = density <= sparse_thr
    # Always allow the very narrow bridge/anchor itself to seed the floor; density gating is applied to the final component statistics.
    anchor_y0 = int(max(y_min, floor_top - top_margin))
    anchor_y1 = int(min(y_max, floor_top + anchor_depth))
    central_anchor_window = (yy >= anchor_y0) & (yy <= anchor_y1) & (xx >= bridge_x0) & (xx <= bridge_x1)

    hanging_window = (yy >= y_seed0) & (yy <= int(y_max)) & (xx >= lat_x0) & (xx <= lat_x1)
    if bool(getattr(args, "me_floor_use_sparse_gating", True)) and bool(getattr(args, "me_floor_hanging_only", True)):
        candidate = tissue & hanging_window & (sparse_mask | central_anchor_window)
    else:
        candidate = tissue & hanging_window

    # Mildly include marker pixels near the floor after DAPI floor is proven present, but only
    # as support around tissue; this cannot create ME without the DAPI bridge above.
    marker_support = np.zeros((H, W), dtype=bool)
    for extra in [cfos_lbl, npy_lbl, pomc_lbl]:
        if extra is None:
            continue
        arr = np.asarray(extra)
        if arr.shape == (H, W) and int(np.max(arr)) > 0:
            marker_support |= (arr > 0) & hanging_window & (sparse_mask | central_anchor_window)
    if np.any(marker_support):
        try:
            marker_support = binary_dilation(marker_support, disk(max(1, um_to_px(8.0, args, min_px=1))))
        except Exception:
            pass
        candidate |= marker_support & tissue

    dil = um_to_px(float(getattr(args, "me_floor_dilate_um", 18.0)), args, min_px=1)
    close = um_to_px(float(getattr(args, "me_floor_close_um", 28.0)), args, min_px=1)
    try:
        candidate = binary_closing(candidate.astype(bool), disk(max(1, close)))
        candidate = binary_dilation(candidate.astype(bool), disk(max(1, dil)))
        candidate &= hanging_window
        candidate &= binary_dilation(tissue.astype(bool), disk(max(1, um_to_px(10.0, args, min_px=1))))
        if bool(getattr(args, "me_floor_use_sparse_gating", True)) and bool(getattr(args, "me_floor_hanging_only", True)):
            sparse_support = binary_dilation((sparse_mask | central_anchor_window).astype(bool), disk(max(1, um_to_px(16.0, args, min_px=1))))
            candidate &= sparse_support
        candidate = remove_small_objects(candidate.astype(bool), min_size=max(64, int(0.00045 * H * W)))
        candidate = remove_small_holes(candidate.astype(bool), area_threshold=max(64, int(0.00025 * H * W)))
    except Exception:
        candidate = candidate.astype(bool)

    # Keep only the connected hanging component that actually touches the central floor bridge
    # and crosses the midline.  This rejects one-sided horizontal tissue and dense ARC shelves.
    anchor_mask = candidate & central_anchor_window
    cc = _connected_label_image(candidate)
    n_cc = int(cc.max()) if cc.size else 0
    good_labels: List[int] = []
    best_sparse_fraction = 0.0
    best_span = 0
    best_depth = 0
    min_hanging_depth_px = um_to_px(float(getattr(args, "me_floor_min_hanging_depth_um", 45.0)), args, min_px=8)
    min_span_px = um_to_px(float(getattr(args, "me_floor_min_cross_midline_span_um", 120.0)), args, min_px=8)
    min_sparse_frac = float(getattr(args, "me_floor_min_sparse_fraction", 0.55))
    require_bridge = bool(getattr(args, "me_floor_require_connected_bridge", True))
    for lab_id in range(1, n_cc + 1):
        comp = cc == lab_id
        if not np.any(comp):
            continue
        if np.count_nonzero(comp & anchor_mask) < max(3, int(0.00002 * H * W)) and mode != "always":
            continue
        ys, xs = np.where(comp)
        if ys.size == 0:
            continue
        crosses = (np.any(xs <= float(mid_x) - min_span_px / 2.0) and np.any(xs >= float(mid_x) + min_span_px / 2.0))
        depth = int(ys.max() - int(floor_top) + 1)
        span_px = int(xs.max() - xs.min() + 1)
        dvals = density[ys, xs]
        sparse_fraction = float(np.mean(dvals <= sparse_thr)) if dvals.size else 0.0
        med_ratio = float(np.nanmedian(dvals) / max(ref_density, 1e-6)) if dvals.size else 0.0
        best_sparse_fraction = max(best_sparse_fraction, sparse_fraction)
        best_span = max(best_span, span_px)
        best_depth = max(best_depth, depth)
        if bool(getattr(args, "me_floor_hanging_only", True)) and sparse_fraction < min_sparse_frac and mode != "always":
            continue
        if require_bridge and (not crosses) and mode != "always":
            continue
        if span_px < min_span_px and mode != "always":
            continue
        if depth < min_hanging_depth_px and mode != "always":
            continue
        # If the median density is dense-ARC-like, reject unless user forced always.
        if bool(getattr(args, "me_floor_use_sparse_gating", True)) and med_ratio > max(0.95, ratio * 1.35) and mode != "always":
            continue
        good_labels.append(int(lab_id))

    if not good_labels:
        return _absent("no_connected_sparse_hanging_floor_component", {
            "me_floor_top_y_candidate": int(y_seed0),
            "me_floor_search_lo": int(lo), "me_floor_search_hi": int(hi),
            "me_floor_vent_bottom_y": int(vent_bottom_y),
            "me_floor_bridge_nuclei": int(bridge_nuclei),
            "me_floor_core_nuclei": int(core_nuclei),
            "me_floor_total_nuclei": int(total_nuclei),
            "me_floor_component_count": int(n_cc),
            "me_floor_good_component_count": 0,
            "me_floor_bridge_connected": False,
            "me_floor_cross_midline": False,
            "me_floor_best_sparse_fraction": float(best_sparse_fraction),
            "me_floor_best_span_px": int(best_span),
            "me_floor_best_depth_px": int(best_depth),
            "me_floor_sparse_density_ref": float(ref_density),
            "me_floor_sparse_threshold": float(sparse_thr),
            "me_floor_max_dense_arc_ratio": float(ratio),
        })

    me_core = np.isin(cc, good_labels)
    # Use only gentle smoothing around the accepted component.  Avoid default convex hulls,
    # because they can convert the dense ARC shoulder into false ME.
    try:
        me_mask = binary_closing(me_core.astype(bool), disk(max(1, close)))
        me_mask = binary_dilation(me_mask.astype(bool), disk(max(1, dil)))
    except Exception:
        me_mask = me_core.astype(bool)

    if bool(getattr(args, "me_floor_use_convex_hull", False)):
        try:
            hull = convex_hull_image(me_core.astype(bool))
            # Strong guard: hull may only occupy sparse hanging support.
            sparse_support = binary_dilation((sparse_mask | central_anchor_window).astype(bool), disk(max(1, um_to_px(20.0, args, min_px=1))))
            me_mask |= (hull & sparse_support & hanging_window)
        except Exception:
            pass

    # Final guards: attached to floor tissue, not a horizontal band, not dense ARC.
    try:
        tissue_support = binary_dilation(tissue.astype(bool), disk(max(1, um_to_px(12.0, args, min_px=1))))
    except Exception:
        tissue_support = tissue.astype(bool)
    me_mask &= tissue_support & hanging_window
    if bool(getattr(args, "me_floor_hanging_only", True)) and bool(getattr(args, "me_floor_use_sparse_gating", True)):
        sparse_support = binary_dilation((sparse_mask | central_anchor_window).astype(bool), disk(max(1, um_to_px(18.0, args, min_px=1))))
        me_mask &= sparse_support
    me_mask &= (yy >= max(y_min, int(floor_top) - top_margin))

    if int(np.count_nonzero(me_mask)) < max(50, int(0.00035 * H * W)):
        return _absent("accepted_floor_component_erased_by_final_guards", {
            "me_floor_top_y_candidate": int(y_seed0),
            "me_floor_component_count": int(n_cc),
            "me_floor_good_component_count": int(len(good_labels)),
            "me_floor_sparse_density_ref": float(ref_density),
            "me_floor_sparse_threshold": float(sparse_thr),
        })

    # Recompute final stats.
    fs = _component_bbox_and_stats(me_mask, density, float(mid_x), int(floor_top), args)
    final_vals = density[me_mask]
    final_sparse_fraction = float(np.mean(final_vals <= sparse_thr)) if final_vals.size else 0.0
    fs["sparse_fraction"] = final_sparse_fraction
    me_top_final = int(np.min(np.where(me_mask)[0])) if np.any(me_mask) else int(y_seed0)
    _region_debug(args,
        "ME present: top=%s area=%s comps=%s good=%s sparse_frac=%.2f span_px=%s depth_px=%s ref=%.4g thr=%.4g",
        me_top_final, int(np.count_nonzero(me_mask)), int(n_cc), int(len(good_labels)), final_sparse_fraction, int(fs.get("span_px", 0)), int(fs.get("depth_px", 0)), float(ref_density), float(sparse_thr)
    )
    return me_mask.astype(bool), {
        "me_floor_present": True,
        "me_floor_detection_mode": mode,
        "me_floor_method": "ventricle_wall_projection_component_v16",
        "me_top_y": int(me_top_final),
        "me_floor_top_y": int(me_top_final),
        "me_floor_raw_floor_top_y": int(floor_top),
        "me_floor_preliminary_y": int(preliminary_me_top_y),
        "me_floor_vent_bottom_y": int(vent_bottom_y),
        "me_floor_search_lo": int(lo),
        "me_floor_search_hi": int(hi),
        "me_floor_core_x0": int(core_x0),
        "me_floor_core_x1": int(core_x1),
        "me_floor_bridge_x0": int(bridge_x0),
        "me_floor_bridge_x1": int(bridge_x1),
        "me_floor_lat_x0": int(lat_x0),
        "me_floor_lat_x1": int(lat_x1),
        "me_floor_bridge_nuclei": int(bridge_nuclei),
        "me_floor_core_nuclei": int(core_nuclei),
        "me_floor_total_nuclei": int(total_nuclei),
        "me_floor_max_bridge_density": float(max_bridge_density),
        "me_floor_component_count": int(n_cc),
        "me_floor_good_component_count": int(len(good_labels)),
        "me_floor_bridge_connected": True,
        "me_floor_cross_midline": bool(fs.get("crosses_midline", False)),
        "me_floor_min_cross_midline_span_px": int(min_span_px),
        "me_floor_hanging_depth_px": int(fs.get("depth_px", 0)),
        "me_floor_horizontal_span_px": int(fs.get("span_px", 0)),
        "me_floor_sparse_density_ref": float(ref_density),
        "me_floor_sparse_threshold": float(sparse_thr),
        "me_floor_max_dense_arc_ratio": float(ratio),
        "me_floor_sparse_fraction_best": float(best_sparse_fraction),
        "me_floor_sparse_fraction_final": float(final_sparse_fraction),
        "me_floor_mask_area_px": int(np.count_nonzero(me_mask)),
        "me_floor_seed_area_px": int(np.count_nonzero(candidate)),
        "me_floor_oval_area_px": 0,
        "me_floor_hull_area_px": 0,
    }

def _arc_tissue_fill_mask(dapi_lbl: Any, tissue: Any, args: argparse.Namespace, y_min: int, y_max: int, x_min: int, x_max: int) -> Any:
    """Build an ARC support over the visible tissue field for ARC/ME-only images."""
    dapi = np.asarray(dapi_lbl)
    H, W = dapi.shape
    tissue = np.asarray(tissue).astype(bool)
    yy, xx = np.indices((H, W))
    pad_y = um_to_px(35.0, args, min_px=8)
    pad_x = um_to_px(60.0, args, min_px=12)
    bbox = (yy >= max(0, int(y_min) - pad_y)) & (yy <= min(H - 1, int(y_max) + pad_y)) & (xx >= max(0, int(x_min) - pad_x)) & (xx <= min(W - 1, int(x_max) + pad_x))
    mask = tissue.astype(bool) & bbox
    try:
        # Convert sparse nuclei into a readable oval/tissue contour while preserving gaps.
        close = um_to_px(float(getattr(args, "region_tissue_close_um", 7.5)) * 1.5, args, min_px=2)
        mask = binary_dilation(mask, disk(max(1, close)))
        mask = binary_closing(mask, disk(max(1, close)))
        mask = remove_small_holes(mask.astype(bool), area_threshold=max(64, int(0.001 * H * W)))
        mask = remove_small_objects(mask.astype(bool), min_size=max(64, int(0.0005 * H * W)))
    except Exception:
        pass
    if int(np.count_nonzero(mask)) < 50:
        mask = (dapi > 0) & bbox
    return mask.astype(bool)

def _fill_short_false_runs_1d(mask: Any, max_gap: int) -> Any:
    """Fill short False gaps bounded by True values in a 1-D boolean vector."""
    a = np.asarray(mask, dtype=bool).copy()
    n = int(a.size)
    i = 0
    max_gap = int(max(0, max_gap))
    while i < n:
        if a[i]:
            i += 1
            continue
        j = i
        while j < n and not a[j]:
            j += 1
        if i > 0 and j < n and (j - i) <= max_gap:
            a[i:j] = True
        i = j
    return a


def _true_runs_1d(mask: Any) -> List[Tuple[int, int]]:
    """Return inclusive (start, end) runs of True values."""
    a = np.asarray(mask, dtype=bool)
    runs: List[Tuple[int, int]] = []
    i = 0
    while i < int(a.size):
        if not a[i]:
            i += 1
            continue
        j = i + 1
        while j < int(a.size) and a[j]:
            j += 1
        runs.append((int(i), int(j - 1)))
        i = j
    return runs


def _robust_wall_line_fit(y: Any, x: Any, min_rows: int, max_abs_slope: float) -> Dict[str, Any]:
    """Robustly fit x = slope*y + intercept with iterative MAD trimming."""
    yy = np.asarray(y, dtype=float)
    xx = np.asarray(x, dtype=float)
    ok = np.isfinite(yy) & np.isfinite(xx)
    yy = yy[ok]
    xx = xx[ok]
    n0 = int(yy.size)
    if n0 == 0:
        return {"slope": 0.0, "intercept": float("nan"), "rmse": float("nan"), "n": 0, "status": "no_points"}
    if n0 == 1:
        return {"slope": 0.0, "intercept": float(xx[0]), "rmse": 0.0, "n": 1, "status": "single_point"}

    keep = np.ones(n0, dtype=bool)
    slope = 0.0
    intercept = float(np.nanmedian(xx))
    min_keep = max(6, min(int(min_rows), n0) // 2)
    for _ in range(6):
        if int(np.count_nonzero(keep)) < 2:
            break
        try:
            slope, intercept = np.polyfit(yy[keep], xx[keep], 1)
        except Exception:
            slope = 0.0
            intercept = float(np.nanmedian(xx[keep]))
        pred = slope * yy + intercept
        resid = xx - pred
        med = float(np.nanmedian(resid[keep]))
        mad = float(np.nanmedian(np.abs(resid[keep] - med)))
        scale = max(1.5, 1.4826 * mad)
        new_keep = np.abs(resid - med) <= max(3.0, 3.0 * scale)
        if int(np.count_nonzero(new_keep)) < min_keep:
            break
        if np.array_equal(new_keep, keep):
            keep = new_keep
            break
        keep = new_keep

    max_abs_slope = float(max(0.05, max_abs_slope))
    slope = float(np.clip(float(slope), -max_abs_slope, max_abs_slope))
    intercept = float(np.nanmedian(xx[keep] - slope * yy[keep])) if np.any(keep) else float(np.nanmedian(xx - slope * yy))
    residual = xx[keep] - (slope * yy[keep] + intercept) if np.any(keep) else xx - (slope * yy + intercept)
    rmse = float(np.sqrt(np.nanmean(residual ** 2))) if residual.size else float("nan")
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "rmse": float(rmse),
        "n": int(np.count_nonzero(keep)),
        "n_input": int(n0),
        "status": "ok" if int(np.count_nonzero(keep)) >= 2 else "fallback",
    }


def _ventricle_wall_support_fraction(tissue: Any, y: int, x: float, side: str, probe_px: int) -> float:
    """Measure DAPI-derived tissue immediately outside a putative lumen wall."""
    t = np.asarray(tissue, dtype=bool)
    H, W = t.shape
    yy = int(max(0, min(H - 1, y)))
    xi = int(round(float(x)))
    probe_px = int(max(2, probe_px))
    if str(side).lower() == "left":
        x0 = max(0, xi - probe_px)
        x1 = min(W, xi + 3)
    else:
        x0 = max(0, xi - 2)
        x1 = min(W, xi + probe_px + 1)
    if x1 <= x0:
        return 0.0
    return float(np.mean(t[yy, x0:x1]))



def _smooth_edge_profile(values: Any, sigma_rows: float) -> Any:
    """Smooth a row-wise wall/width profile without changing its length."""
    arr = np.asarray(values, dtype=float)
    if arr.size <= 2:
        return arr.copy()
    try:
        from scipy.ndimage import gaussian_filter1d
        return gaussian_filter1d(arr, sigma=max(0.5, float(sigma_rows)), mode="nearest")
    except Exception:
        return arr.copy()


def _oval_wall_profile(
    rows: Any,
    left_values: Any,
    right_values: Any,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """Detect a closed/oval lumen and define tangent rows at its lateral bulge.

    For an oval lumen, extending the lower rounded cap would make projected lines
    converge incorrectly.  Instead, the algorithm finds the widest part of the oval,
    estimates local side tangents there, and projects those tangents ventrally.
    """
    rr = np.asarray(rows, dtype=int)
    ll = np.asarray(left_values, dtype=float)
    rx = np.asarray(right_values, dtype=float)
    ok = np.isfinite(rr) & np.isfinite(ll) & np.isfinite(rx) & (rx > ll + 1)
    rr = rr[ok]; ll = ll[ok]; rx = rx[ok]
    info: Dict[str, Any] = {
        "detected": False,
        "mode": str(getattr(args, "wall_oval_mode", "auto")).lower(),
        "peak_y": None,
        "peak_width_px": float("nan"),
        "top_width_px": float("nan"),
        "bottom_width_px": float("nan"),
        "bulge_ratio_top": float("nan"),
        "bulge_ratio_bottom": float("nan"),
        "peak_fraction": float("nan"),
        "tangent_rows": np.asarray([], dtype=int),
        "left_smooth": ll,
        "right_smooth": rx,
    }
    if rr.size < 8:
        return info
    # The rows passed here are contiguous/interpolated, so profile smoothing is stable.
    sigma = max(1.0, 0.025 * float(rr.size))
    lsm = _smooth_edge_profile(ll, sigma)
    rsm = _smooth_edge_profile(rx, sigma)
    widths = np.maximum(1.0, rsm - lsm)
    wsm = _smooth_edge_profile(widths, sigma)
    peak_i = int(np.nanargmax(wsm))
    peak_y = int(rr[peak_i])
    edge_n = int(max(3, round(0.15 * rr.size)))
    top_w = float(np.nanmedian(wsm[:edge_n]))
    bot_w = float(np.nanmedian(wsm[-edge_n:]))
    peak_w = float(wsm[peak_i])
    ratio_top = float(peak_w / max(1.0, top_w))
    ratio_bottom = float(peak_w / max(1.0, bot_w))
    peak_frac = float((peak_i) / max(1, rr.size - 1))
    min_h_px = um_to_px(float(getattr(args, "wall_oval_min_height_um", 70.0)), args, min_px=20)
    ratio_thr = float(max(1.01, getattr(args, "wall_oval_min_bulge_ratio", 1.16)))
    mode = str(getattr(args, "wall_oval_mode", "auto")).lower()
    auto_detected = (
        int(rr[-1] - rr[0] + 1) >= min_h_px
        and 0.10 <= peak_frac <= 0.76
        and ratio_top >= ratio_thr
        and ratio_bottom >= ratio_thr
    )
    detected = mode == "always" or (mode == "auto" and auto_detected)
    if mode == "never":
        detected = False
    win_frac = float(np.clip(float(getattr(args, "wall_oval_tangent_window_frac", 0.16)), 0.05, 0.50))
    half_n = int(max(4, round(0.5 * win_frac * rr.size)))
    a = max(0, peak_i - half_n)
    b = min(rr.size, peak_i + half_n + 1)
    tangent_rows = rr[a:b]
    info.update({
        "detected": bool(detected),
        "peak_y": int(peak_y),
        "peak_width_px": float(peak_w),
        "top_width_px": float(top_w),
        "bottom_width_px": float(bot_w),
        "bulge_ratio_top": float(ratio_top),
        "bulge_ratio_bottom": float(ratio_bottom),
        "peak_fraction": float(peak_frac),
        "tangent_rows": tangent_rows.astype(int),
        "left_smooth": lsm.astype(float),
        "right_smooth": rsm.astype(float),
        "rows": rr.astype(int),
    })
    return info


def _line_coordinates_at_y(wall: Dict[str, Any], y: float) -> Tuple[float, float]:
    ys = np.asarray(wall.get("line_y", []), dtype=float)
    lx = np.asarray(wall.get("left_x", []), dtype=float)
    rx = np.asarray(wall.get("right_x", []), dtype=float)
    good = np.isfinite(ys) & np.isfinite(lx) & np.isfinite(rx)
    if int(np.count_nonzero(good)) < 2:
        mid = float(wall.get("mid_x", 0.0))
        half = float(wall.get("median_half_gap_px", 10.0))
        return mid - half, mid + half
    ys = ys[good]; lx = lx[good]; rx = rx[good]
    order = np.argsort(ys); ys = ys[order]; lx = lx[order]; rx = rx[order]
    yy = float(np.clip(float(y), float(ys[0]), float(ys[-1])))
    left_here = float(np.interp(yy, ys, lx))
    right_here = float(np.interp(yy, ys, rx))
    if right_here < left_here:
        left_here, right_here = right_here, left_here
    return left_here, right_here


def _infer_closed_ventricle_floor_curve_primary(
    tissue: Any,
    wall: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """Find the curved ventral boundary of the central lumen.

    The floor is derived from the bottom envelope of the central empty component that
    overlaps the observed third-ventricle lumen.  If that component remains open to the
    ventral image boundary, no ME floor exists.
    """
    t = np.asarray(tissue, dtype=bool)
    H, W = t.shape
    corridor = np.asarray(wall.get("corridor", np.zeros((H, W), dtype=bool)), dtype=bool)
    obs_y = np.asarray(wall.get("observed_y", []), dtype=int)
    obs_l = np.asarray(wall.get("observed_left_x", []), dtype=float)
    obs_r = np.asarray(wall.get("observed_right_x", []), dtype=float)
    line_y = np.asarray(wall.get("line_y", []), dtype=int)
    line_l = np.asarray(wall.get("left_x", []), dtype=float)
    line_r = np.asarray(wall.get("right_x", []), dtype=float)

    # Seal only small inter-nuclear gaps. This does not fill a genuinely open ventricle.
    close_px = um_to_px(float(getattr(args, "wall_floor_component_close_um", 10.0)), args, min_px=0)
    sealed = t.copy()
    if close_px > 0:
        try:
            sealed = binary_closing(sealed.astype(bool), disk(max(1, close_px)))
        except Exception:
            sealed = t.copy()

    domain = corridor.copy()
    # Include the actually observed lumen above the projection start, especially for an oval.
    for y, lx, rx in zip(obs_y, obs_l, obs_r):
        if y < 0 or y >= H or not np.isfinite(lx) or not np.isfinite(rx):
            continue
        x0 = int(max(0, min(W - 1, math.ceil(min(lx, rx)))))
        x1 = int(max(0, min(W - 1, math.floor(max(lx, rx)))))
        if x1 >= x0:
            domain[int(y), x0:x1 + 1] = True

    empty = (~sealed) & domain
    seed = np.zeros((H, W), dtype=bool)
    # Seed from the observed lumen interior, not from arbitrary background below it.
    for y, lx, rx in zip(obs_y, obs_l, obs_r):
        if y < 0 or y >= H or not np.isfinite(lx) or not np.isfinite(rx):
            continue
        x0 = int(max(0, min(W - 1, math.ceil(min(lx, rx) + 2))))
        x1 = int(max(0, min(W - 1, math.floor(max(lx, rx) - 2))))
        if x1 >= x0:
            seed[int(y), x0:x1 + 1] = True
    seed &= empty

    try:
        cc = label(empty.astype(bool), connectivity=1).astype(np.int32)
    except Exception:
        cc = _connected_label_image(empty)
    ids = cc[seed]
    ids = ids[ids > 0]
    if ids.size == 0:
        return {
            "floor_found": False,
            "open_bottom": False,
            "status": "no_seeded_lumen_component",
            "lumen_mask": np.zeros((H, W), dtype=bool),
            "curve_x": np.asarray([], dtype=float),
            "curve_y": np.asarray([], dtype=float),
            "curve_y_raw": np.asarray([], dtype=float),
        }
    counts = np.bincount(ids)
    counts[0] = 0
    lumen_id = int(np.argmax(counts))
    lumen = cc == lumen_id
    ly, lx = np.where(lumen)
    if ly.size == 0:
        return {
            "floor_found": False,
            "open_bottom": False,
            "status": "empty_seeded_lumen_component",
            "lumen_mask": lumen,
            "curve_x": np.asarray([], dtype=float),
            "curve_y": np.asarray([], dtype=float),
            "curve_y_raw": np.asarray([], dtype=float),
        }

    projection_end = int(wall.get("projection_y_end", H - 1))
    open_margin = um_to_px(float(getattr(args, "wall_floor_open_bottom_margin_um", 12.0)), args, min_px=1)
    max_lumen_y = int(np.max(ly))
    open_bottom = bool(max_lumen_y >= min(H - 1, projection_end) - open_margin)

    floor_raw = np.full(W, np.nan, dtype=float)
    for x in np.unique(lx):
        yy = ly[lx == x]
        if yy.size:
            floor_raw[int(x)] = float(np.max(yy))
    valid = np.isfinite(floor_raw)
    runs = _true_runs_1d(valid)
    mid_x = float(wall.get("mid_x", W / 2.0))
    chosen: Optional[Tuple[int, int]] = None
    for a, b in runs:
        if a <= mid_x <= b:
            chosen = (a, b)
            break
    if chosen is None and runs:
        chosen = max(runs, key=lambda ab: (ab[1] - ab[0] + 1) - 0.25 * abs(((ab[0] + ab[1]) / 2.0) - mid_x))
    if chosen is None:
        return {
            "floor_found": False,
            "open_bottom": bool(open_bottom),
            "status": "no_contiguous_floor_columns",
            "lumen_mask": lumen,
            "curve_x": np.asarray([], dtype=float),
            "curve_y": np.asarray([], dtype=float),
            "curve_y_raw": np.asarray([], dtype=float),
            "max_lumen_y": int(max_lumen_y),
        }

    a, b = chosen
    curve_x = np.arange(a, b + 1, dtype=float)
    raw = floor_raw[a:b + 1]
    good = np.isfinite(raw)
    if int(np.count_nonzero(good)) < 3:
        return {
            "floor_found": False,
            "open_bottom": bool(open_bottom),
            "status": "too_few_floor_columns",
            "lumen_mask": lumen,
            "curve_x": np.asarray([], dtype=float),
            "curve_y": np.asarray([], dtype=float),
            "curve_y_raw": np.asarray([], dtype=float),
            "max_lumen_y": int(max_lumen_y),
        }
    raw_interp = np.interp(curve_x, curve_x[good], raw[good])
    sigma_px = max(1.0, um_to_px(float(getattr(args, "wall_floor_curve_smooth_um", 12.0)), args, min_px=1) / 3.0)
    smooth = _smooth_edge_profile(raw_interp, sigma_px)
    bottom_ref = float(np.nanpercentile(smooth, 90))
    max_rise = float(um_to_px(float(getattr(args, "wall_floor_max_rise_um", 90.0)), args, min_px=1))
    # Do not let sidewall portions of the empty component create a high false floor.
    gate = np.maximum(smooth, bottom_ref - max_rise)
    floor_found = bool((not open_bottom) and curve_x.size >= max(8, um_to_px(25.0, args, min_px=8)))
    return {
        "floor_found": bool(floor_found),
        "open_bottom": bool(open_bottom),
        "status": "closed_lumen_floor_curve" if floor_found else ("lumen_open_to_bottom" if open_bottom else "floor_curve_insufficient"),
        "lumen_mask": lumen.astype(bool),
        "curve_x": curve_x.astype(float),
        "curve_y": gate.astype(float),
        "curve_y_raw": raw_interp.astype(float),
        "floor_bottom_y": float(bottom_ref),
        "max_lumen_y": int(max_lumen_y),
        "curve_columns": int(curve_x.size),
        "curve_x0": int(a),
        "curve_x1": int(b),
    }


def _sparse_floor_bridge_rescue(
    dapi_lbl: Any,
    tissue: Any,
    wall: Dict[str, Any],
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """Recover a real but sparsely nucleated ME floor when the lumen mask leaks.

    The normal path remains the seeded empty-lumen component.  This rescue is invoked
    only after that path fails.  It looks for a low-density, bilateral DAPI bridge below
    the ventral lumen tip and inside the projected wall corridor.  The linking kernel is
    deliberately much wider horizontally than vertically, so sparse floor nuclei can be
    joined without filling the open third-ventricle lumen above them.
    """
    dapi = np.asarray(dapi_lbl).astype(np.int32, copy=False)
    t = np.asarray(tissue, dtype=bool)
    H, W = dapi.shape
    empty_result: Dict[str, Any] = {
        "floor_found": False,
        "open_bottom": True,
        "status": "sparse_floor_rescue_not_found",
        "sparse_floor_rescue": False,
        "sparse_floor_candidate_nuclei": 0,
        "sparse_floor_left_nuclei": 0,
        "sparse_floor_right_nuclei": 0,
        "sparse_floor_midline_nuclei": 0,
        "sparse_floor_span_px": 0.0,
        "sparse_floor_depth_px": 0.0,
        "sparse_floor_sparse_fraction": 0.0,
        "sparse_floor_density_ratio": float("nan"),
        "lumen_mask": np.zeros((H, W), dtype=bool),
        "curve_x": np.asarray([], dtype=float),
        "curve_y": np.asarray([], dtype=float),
        "curve_y_raw": np.asarray([], dtype=float),
    }
    if not bool(getattr(args, "wall_floor_sparse_rescue", True)):
        empty_result["status"] = "sparse_floor_rescue_disabled"
        return empty_result

    corridor = np.asarray(wall.get("corridor", np.zeros((H, W), dtype=bool)), dtype=bool)
    line_y = np.asarray(wall.get("line_y", []), dtype=float)
    if corridor.shape != (H, W) or line_y.size < 2 or not np.any(corridor):
        empty_result["status"] = "sparse_floor_no_projected_corridor"
        return empty_result

    mid_x = float(wall.get("mid_x", W / 2.0))
    projection_end = int(min(H - 1, wall.get("projection_y_end", H - 1)))
    ventral_tip = int(wall.get("ventral_tip_y", wall.get("main_run_end_y", int(np.nanmax(line_y)))))
    top_margin = um_to_px(float(getattr(args, "wall_floor_sparse_top_margin_um", 25.0)), args, min_px=2)
    search_start = int(max(int(np.nanmin(line_y)), ventral_tip - top_margin))
    search_end = int(min(H - 1, projection_end))
    if search_end <= search_start + 2:
        empty_result["status"] = "sparse_floor_search_window_too_short"
        return empty_result

    side_excl = float(um_to_px(float(getattr(args, "wall_sidewall_exclusion_um", 18.0)), args, min_px=0))
    records: List[Dict[str, Any]] = []
    all_props = list(regionprops(dapi))
    for rp in all_props:
        y, x = float(rp.centroid[0]), float(rp.centroid[1])
        yi = int(max(0, min(H - 1, round(y))))
        xi = int(max(0, min(W - 1, round(x))))
        if y < search_start or y > search_end or not bool(corridor[yi, xi]):
            continue
        left_here, right_here = _line_coordinates_at_y(wall, y)
        if x < left_here or x > right_here:
            continue
        # The rescue should be driven by floor evidence, not the ventricular sidewalls.
        if min(x - left_here, right_here - x) < max(2.0, 0.55 * side_excl):
            continue
        records.append({"label": int(rp.label), "y": y, "x": x, "yi": yi, "xi": xi})

    min_n = int(max(2, getattr(args, "wall_floor_sparse_min_nuclei", 3)))
    if len(records) < min_n:
        empty_result.update({
            "status": "sparse_floor_too_few_corridor_nuclei",
            "sparse_floor_candidate_nuclei": int(len(records)),
        })
        return empty_result

    point_mask = np.zeros((H, W), dtype=bool)
    for r in records:
        point_mask[int(r["yi"]), int(r["xi"])] = True
    hx = um_to_px(float(getattr(args, "wall_floor_sparse_horizontal_link_um", 55.0)), args, min_px=8)
    hy = um_to_px(float(getattr(args, "wall_floor_sparse_vertical_link_um", 10.0)), args, min_px=2)
    yy_s, xx_s = np.ogrid[-hy:hy + 1, -hx:hx + 1]
    structure = ((yy_s / max(1.0, float(hy))) ** 2 + (xx_s / max(1.0, float(hx))) ** 2) <= 1.0
    try:
        from scipy import ndimage as ndi
        linked = ndi.binary_dilation(point_mask, structure=structure)
        # A small close repairs one-pixel gaps but does not add a second large dilation.
        close_hx = max(2, hx // 3)
        close_hy = max(1, hy // 2)
        yy_c, xx_c = np.ogrid[-close_hy:close_hy + 1, -close_hx:close_hx + 1]
        close_structure = ((yy_c / max(1.0, float(close_hy))) ** 2 + (xx_c / max(1.0, float(close_hx))) ** 2) <= 1.0
        linked = ndi.binary_closing(linked, structure=close_structure)
    except Exception:
        linked = point_mask.copy()
    row_gate = np.zeros((H, 1), dtype=bool)
    row_gate[search_start:search_end + 1, 0] = True
    linked &= corridor & row_gate
    cc = _connected_label_image(linked)

    dens_sigma = um_to_px(float(getattr(args, "me_floor_density_sigma_um", 24.0)), args, min_px=3)
    density = _density_map_from_label_centroids(dapi, sigma_px=dens_sigma)
    ref_band = um_to_px(float(getattr(args, "me_floor_arc_reference_band_um", 150.0)), args, min_px=20)
    ref_vals: List[float] = []
    for rp in all_props:
        y, x = float(rp.centroid[0]), float(rp.centroid[1])
        yi = int(max(0, min(H - 1, round(y))))
        xi = int(max(0, min(W - 1, round(x))))
        if y < max(0, search_start - ref_band) or y > search_end:
            continue
        if not bool(corridor[yi, xi]):
            ref_vals.append(float(density[yi, xi]))
    if len(ref_vals) < 10:
        for rp in all_props:
            y, x = float(rp.centroid[0]), float(rp.centroid[1])
            yi = int(max(0, min(H - 1, round(y))))
            xi = int(max(0, min(W - 1, round(x))))
            if max(0, search_start - ref_band) <= y < search_start and abs(x - mid_x) <= 0.45 * W:
                ref_vals.append(float(density[yi, xi]))
    ref_density = float(np.nanpercentile(np.asarray(ref_vals, dtype=float), 75)) if ref_vals else float("nan")
    density_ratio_limit = float(getattr(args, "me_floor_max_dense_arc_ratio", 0.62))
    sparse_fraction_min = float(getattr(args, "me_floor_min_sparse_fraction", 0.55))
    sparse_thr = ref_density * density_ratio_limit if np.isfinite(ref_density) and ref_density > 0 else float("inf")

    min_span_px = um_to_px(float(getattr(args, "me_floor_min_cross_midline_span_um", 120.0)), args, min_px=12)
    min_depth_px = um_to_px(float(getattr(args, "wall_me_min_depth_um", 18.0)), args, min_px=3)
    min_bilateral = int(max(0, getattr(args, "wall_floor_sparse_min_bilateral_nuclei", 1)))
    mid_half = um_to_px(float(getattr(args, "wall_me_midline_half_width_um", 140.0)), args, min_px=10)

    best: Optional[Dict[str, Any]] = None
    max_cc = int(np.max(cc)) if cc.size else 0
    for cid in range(1, max_cc + 1):
        comp = cc == cid
        if not np.any(comp):
            continue
        comp_records = [r for r in records if bool(comp[int(r["yi"]), int(r["xi"])])]
        n = len(comp_records)
        if n < min_n:
            continue
        xs = np.asarray([float(r["x"]) for r in comp_records], dtype=float)
        ys = np.asarray([float(r["y"]) for r in comp_records], dtype=float)
        left_n = int(np.count_nonzero(xs < mid_x))
        right_n = int(np.count_nonzero(xs > mid_x))
        mid_n = int(np.count_nonzero(np.abs(xs - mid_x) <= mid_half))
        comp_y, comp_x = np.where(comp)
        comp_span = float(comp_x.max() - comp_x.min() + 1) if comp_x.size else 0.0
        actual_span = float(xs.max() - xs.min()) if xs.size > 1 else 0.0
        depth = float(comp_y.max() - comp_y.min() + 1) if comp_y.size else 0.0
        # Original nuclei, not only the dilated linking component, must be represented
        # on both sides of the fitted midline.  This preserves the established rule that
        # a one-sided hanging fragment is not ME.
        crosses = bool(left_n >= min_bilateral and right_n >= min_bilateral)
        dens_vals = np.asarray([float(density[int(r["yi"]), int(r["xi"])]) for r in comp_records], dtype=float)
        sparse_fraction = float(np.mean(dens_vals <= sparse_thr)) if np.isfinite(sparse_thr) and dens_vals.size else 1.0
        density_ratio = float(np.nanmedian(dens_vals) / max(ref_density, 1e-9)) if dens_vals.size and np.isfinite(ref_density) and ref_density > 0 else float("nan")
        if not crosses or comp_span < min_span_px or depth < min_depth_px:
            continue
        if np.isfinite(ref_density) and sparse_fraction < sparse_fraction_min:
            continue
        if np.isfinite(density_ratio) and density_ratio > max(0.95, density_ratio_limit * 1.35):
            continue
        center_penalty = abs(float(np.nanmedian(xs)) - mid_x)
        score = 30.0 * n + comp_span + 0.5 * actual_span + 0.5 * depth + 120.0 * sparse_fraction - 0.15 * center_penalty
        cand = {
            "component_id": int(cid), "records": comp_records, "n": int(n),
            "left_n": left_n, "right_n": right_n, "mid_n": mid_n,
            "span_px": comp_span, "actual_span_px": actual_span, "depth_px": depth,
            "sparse_fraction": sparse_fraction, "density_ratio": density_ratio,
            "score": float(score), "component": comp,
        }
        if best is None or float(cand["score"]) > float(best["score"]):
            best = cand

    if best is None:
        empty_result.update({
            "status": "sparse_floor_no_bilateral_low_density_bridge",
            "sparse_floor_candidate_nuclei": int(len(records)),
            "sparse_floor_reference_density": float(ref_density),
            "sparse_floor_density_threshold": float(sparse_thr),
        })
        return empty_result

    selected_labels = [int(r["label"]) for r in best["records"]]
    selected_mask = np.isin(dapi, np.asarray(selected_labels, dtype=np.int32))
    sy, sx = np.where(selected_mask)
    if sy.size == 0:
        empty_result["status"] = "sparse_floor_selected_labels_empty"
        return empty_result

    y_ref = float(np.nanmedian(sy))
    left_at_ref, right_at_ref = _line_coordinates_at_y(wall, y_ref)
    curve_margin = um_to_px(float(getattr(args, "wall_floor_sparse_curve_margin_um", 35.0)), args, min_px=0)
    x0 = int(max(0, math.ceil(left_at_ref), int(np.min(sx)) - curve_margin))
    x1 = int(min(W - 1, math.floor(right_at_ref), int(np.max(sx)) + curve_margin))
    if x1 <= x0 + 2:
        empty_result["status"] = "sparse_floor_curve_span_too_short"
        return empty_result

    top_raw = np.full(W, np.nan, dtype=float)
    for xcol in np.unique(sx):
        vals = sy[sx == xcol]
        if vals.size:
            top_raw[int(xcol)] = float(np.min(vals))
    valid_cols = np.isfinite(top_raw[x0:x1 + 1])
    if int(np.count_nonzero(valid_cols)) < 3:
        # Centroids are enough to define a low-count floor even when the masks occupy few columns.
        for r in best["records"]:
            top_raw[int(r["xi"])] = float(r["y"])
        valid_cols = np.isfinite(top_raw[x0:x1 + 1])
    if int(np.count_nonzero(valid_cols)) < 2:
        empty_result["status"] = "sparse_floor_too_few_curve_columns"
        return empty_result

    curve_x = np.arange(x0, x1 + 1, dtype=float)
    local_raw = top_raw[x0:x1 + 1]
    good = np.isfinite(local_raw)
    raw_interp = np.interp(curve_x, curve_x[good], local_raw[good])
    sigma_px = max(1.0, um_to_px(float(getattr(args, "wall_floor_curve_smooth_um", 12.0)), args, min_px=1) / 3.0)
    smooth = _smooth_edge_profile(raw_interp, sigma_px)
    bottom_ref = float(np.nanpercentile(smooth, 90))
    max_rise = float(um_to_px(float(getattr(args, "wall_floor_max_rise_um", 90.0)), args, min_px=1))
    curve = np.maximum(smooth, bottom_ref - max_rise)

    _region_debug(
        args,
        "Sparse floor rescue: n=%s left/right/mid=%s/%s/%s span=%.1f depth=%.1f sparse=%.2f density_ratio=%s search=%s-%s curve=%s-%s",
        int(best["n"]), int(best["left_n"]), int(best["right_n"]), int(best["mid_n"]),
        float(best["span_px"]), float(best["depth_px"]), float(best["sparse_fraction"]),
        f"{float(best['density_ratio']):.2f}" if np.isfinite(float(best["density_ratio"])) else "NA",
        int(search_start), int(search_end), int(x0), int(x1),
    )
    return {
        "floor_found": True,
        "open_bottom": False,
        "status": "sparse_floor_bridge_rescue",
        "sparse_floor_rescue": True,
        "sparse_floor_candidate_nuclei": int(best["n"]),
        "sparse_floor_left_nuclei": int(best["left_n"]),
        "sparse_floor_right_nuclei": int(best["right_n"]),
        "sparse_floor_midline_nuclei": int(best["mid_n"]),
        "sparse_floor_span_px": float(best["span_px"]),
        "sparse_floor_depth_px": float(best["depth_px"]),
        "sparse_floor_sparse_fraction": float(best["sparse_fraction"]),
        "sparse_floor_density_ratio": float(best["density_ratio"]),
        "sparse_floor_reference_density": float(ref_density),
        "sparse_floor_density_threshold": float(sparse_thr),
        "sparse_floor_search_start_y": int(search_start),
        "sparse_floor_search_end_y": int(search_end),
        "sparse_floor_selected_labels": selected_labels,
        "lumen_mask": np.zeros((H, W), dtype=bool),
        "curve_x": curve_x.astype(float),
        "curve_y": curve.astype(float),
        "curve_y_raw": raw_interp.astype(float),
        "floor_bottom_y": float(bottom_ref),
        "max_lumen_y": int(ventral_tip),
        "curve_columns": int(curve_x.size),
        "curve_x0": int(x0),
        "curve_x1": int(x1),
    }


def _infer_closed_ventricle_floor_curve(
    tissue: Any,
    wall: Dict[str, Any],
    args: argparse.Namespace,
    dapi_lbl: Optional[Any] = None,
) -> Dict[str, Any]:
    """Primary closed-lumen floor inference with a conservative sparse-floor fallback."""
    primary = _infer_closed_ventricle_floor_curve_primary(tissue, wall, args)
    if bool(primary.get("floor_found", False)):
        primary.setdefault("sparse_floor_rescue", False)
        return primary
    if dapi_lbl is None or not bool(getattr(args, "wall_floor_sparse_rescue", True)):
        primary.setdefault("sparse_floor_rescue", False)
        return primary
    rescue = _sparse_floor_bridge_rescue(dapi_lbl, tissue, wall, args)
    if bool(rescue.get("floor_found", False)):
        rescue["primary_floor_status"] = str(primary.get("status", "NA"))
        rescue["primary_open_bottom"] = bool(primary.get("open_bottom", False))
        # Preserve the observed lumen mask for QC even though the rescue curve is DAPI-based.
        rescue["lumen_mask"] = np.asarray(primary.get("lumen_mask", rescue.get("lumen_mask")), dtype=bool)
        return rescue
    primary.setdefault("sparse_floor_rescue", False)
    for key, value in rescue.items():
        if key.startswith("sparse_floor_"):
            primary[key] = value
    primary["sparse_floor_rescue_status"] = str(rescue.get("status", "NA"))
    return primary


def _point_is_me_floor(
    wall: Dict[str, Any],
    floor_info: Dict[str, Any],
    y: float,
    x: float,
    args: argparse.Namespace,
) -> bool:
    """Authoritative point rule: below the closed floor curve, not on a sidewall."""
    if not bool(floor_info.get("floor_found", False)):
        return False
    yy = float(y); xx = float(x)
    line_y = np.asarray(wall.get("line_y", []), dtype=float)
    if line_y.size < 2 or yy < float(np.nanmin(line_y)) or yy > float(np.nanmax(line_y)):
        return False
    left_here, right_here = _line_coordinates_at_y(wall, yy)
    if xx < left_here or xx > right_here:
        return False
    cx = np.asarray(floor_info.get("curve_x", []), dtype=float)
    cy = np.asarray(floor_info.get("curve_y", []), dtype=float)
    good = np.isfinite(cx) & np.isfinite(cy)
    if int(np.count_nonzero(good)) < 2:
        return False
    cx = cx[good]; cy = cy[good]
    order = np.argsort(cx); cx = cx[order]; cy = cy[order]
    if xx < float(cx[0]) or xx > float(cx[-1]):
        return False
    floor_y = float(np.interp(xx, cx, cy))
    margin = float(um_to_px(float(getattr(args, "wall_floor_centroid_margin_um", 0.0)), args, min_px=0))
    if yy < floor_y + margin:
        return False
    side_dist = min(xx - left_here, right_here - xx)
    side_excl = float(um_to_px(float(getattr(args, "wall_sidewall_exclusion_um", 18.0)), args, min_px=0))
    release = float(um_to_px(float(getattr(args, "wall_sidewall_release_depth_um", 35.0)), args, min_px=0))
    if side_dist < side_excl and yy < floor_y + release:
        return False
    return True

def fit_projected_ventricle_walls(tissue: Any, vent: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    """Fit and ventrally project the two third-ventricle walls.

    The original robust wall logic is retained.  The only geometric addition is an
    oval-aware branch: when the lumen has a clear interior width maximum and narrows at
    both ends, local tangents are fitted around the widest lateral points and projected
    ventrally.  This matches the user's oval-plus-straight-extension drawing and avoids
    following the rounded ventral cap of a closed oval.
    """
    t = np.asarray(tissue, dtype=bool)
    H, W = t.shape
    y_min = int(vent.get("tissue_y_min", 0))
    y_max = int(vent.get("tissue_y_max", H - 1))
    mid_x = float(vent.get("mid_x", W / 2.0))
    left = np.asarray(vent.get("left_edge", np.full(H, mid_x - W * 0.03)), dtype=float)
    right = np.asarray(vent.get("right_edge", np.full(H, mid_x + W * 0.03)), dtype=float)
    accepted = np.asarray(vent.get("accepted_rows", np.zeros(H, dtype=bool)), dtype=bool)

    probe_px = um_to_px(float(getattr(args, "wall_support_probe_um", 24.0)), args, min_px=4)
    min_support = float(getattr(args, "wall_support_min_fraction", 0.06))
    max_gap_px = float(max(20.0, float(getattr(args, "wall_max_gap_frac", 0.42)) * W))
    mid_tol = float(max(um_to_px(45.0, args, min_px=12), float(getattr(args, "wall_center_tolerance_frac", 0.08)) * W))
    max_asym = float(np.clip(float(getattr(args, "wall_max_asymmetry", 0.40)), 0.05, 0.95))

    valid = np.zeros(H, dtype=bool)
    left_support = np.zeros(H, dtype=np.float32)
    right_support = np.zeros(H, dtype=np.float32)
    for y in range(max(0, y_min), min(H, y_max + 1)):
        if not bool(accepted[y]) or not np.isfinite(left[y]) or not np.isfinite(right[y]):
            continue
        width = float(right[y] - left[y])
        center = float((right[y] + left[y]) / 2.0)
        dl = float(mid_x - left[y])
        dr = float(right[y] - mid_x)
        asym = abs(dl - dr) / max(1.0, dl + dr)
        if width <= 0 or width > max_gap_px or abs(center - mid_x) > mid_tol or dl <= 0 or dr <= 0 or asym > max_asym:
            continue
        ls = _ventricle_wall_support_fraction(t, y, left[y], "left", probe_px)
        rs = _ventricle_wall_support_fraction(t, y, right[y], "right", probe_px)
        left_support[y] = ls
        right_support[y] = rs
        if ls >= min_support and rs >= min_support:
            valid[y] = True

    filled = _fill_short_false_runs_1d(valid, int(getattr(args, "wall_fit_max_gap_rows", 18)))
    runs = _true_runs_1d(filled)
    min_fit_rows = int(max(8, getattr(args, "wall_fit_min_rows", 35)))
    candidate_runs = [r for r in runs if (r[1] - r[0] + 1) >= max(8, min_fit_rows // 2)]
    if not candidate_runs:
        candidate_runs = runs

    if candidate_runs:
        def run_score(run: Tuple[int, int]) -> float:
            a, b = run
            rows = np.arange(a, b + 1)
            support = float(np.nanmean((left_support[rows] + right_support[rows]) / 2.0)) if rows.size else 0.0
            length = float(b - a + 1)
            ventral_bonus = 1.0 + 0.35 * max(0.0, (b - (y_min + y_max) / 2.0) / max(1.0, y_max - y_min))
            return length * (0.5 + support) * ventral_bonus
        run_start, run_end = max(candidate_runs, key=run_score)
        fit_status = "supported_main_run"
    else:
        rows = np.where(accepted & (np.arange(H) >= y_min) & (np.arange(H) <= y_max))[0]
        if rows.size:
            run_start, run_end = int(rows.min()), int(rows.max())
            fit_status = "accepted_rows_fallback"
        else:
            run_start, run_end = int(y_min), int(y_max)
            fit_status = "synthetic_vertical_fallback"

    profile_rows = np.where(valid & (np.arange(H) >= run_start) & (np.arange(H) <= run_end))[0]
    if profile_rows.size < 2:
        profile_rows = np.arange(run_start, run_end + 1, dtype=int)
    observed_y = np.arange(run_start, run_end + 1, dtype=int)
    if profile_rows.size >= 2:
        observed_left = np.interp(observed_y, profile_rows, left[profile_rows])
        observed_right = np.interp(observed_y, profile_rows, right[profile_rows])
    else:
        observed_left = np.full(observed_y.size, mid_x - max(10.0, W * 0.03), dtype=float)
        observed_right = np.full(observed_y.size, mid_x + max(10.0, W * 0.03), dtype=float)

    oval = _oval_wall_profile(observed_y, observed_left, observed_right, args)

    # Default non-oval fit: preserve the successful v15/v16 lower-wall logic.
    fit_rows_all = profile_rows.copy()
    tip_exclude = um_to_px(float(getattr(args, "wall_fit_tip_exclude_um", 6.0)), args, min_px=0)
    if tip_exclude > 0:
        trimmed = fit_rows_all[fit_rows_all <= int(run_end - tip_exclude)]
        if trimmed.size >= max(8, min_fit_rows // 2):
            fit_rows_all = trimmed
    lower_frac = float(np.clip(float(getattr(args, "wall_fit_lower_fraction", 0.55)), 0.15, 1.0))
    take_n = int(max(min_fit_rows, round(fit_rows_all.size * lower_frac)))
    take_n = int(min(max(2, take_n), max(2, fit_rows_all.size)))
    fit_rows = fit_rows_all[-take_n:]
    left_fit_values = left[fit_rows]
    right_fit_values = right[fit_rows]
    max_slope = float(getattr(args, "wall_max_abs_slope", 1.35))
    line_y_start = int(max(y_min, run_start))

    if bool(oval.get("detected", False)):
        tang_rows = np.asarray(oval.get("tangent_rows", []), dtype=int)
        obs_rows = np.asarray(oval.get("rows", observed_y), dtype=int)
        obs_lsm = np.asarray(oval.get("left_smooth", observed_left), dtype=float)
        obs_rsm = np.asarray(oval.get("right_smooth", observed_right), dtype=float)
        if tang_rows.size >= 2 and obs_rows.size == obs_lsm.size == obs_rsm.size:
            fit_rows = tang_rows
            left_fit_values = np.interp(fit_rows, obs_rows, obs_lsm)
            right_fit_values = np.interp(fit_rows, obs_rows, obs_rsm)
            max_slope = float(min(max_slope, max(0.02, getattr(args, "wall_oval_max_tangent_slope", 0.30))))
            line_y_start = int(oval.get("peak_y", fit_rows[len(fit_rows) // 2]))
            fit_status += "_oval_side_tangent"

    left_fit = _robust_wall_line_fit(fit_rows, left_fit_values, min_fit_rows, max_slope)
    right_fit = _robust_wall_line_fit(fit_rows, right_fit_values, min_fit_rows, max_slope)
    oval_left_slope_raw = float(left_fit.get("slope", 0.0))
    oval_right_slope_raw = float(right_fit.get("slope", 0.0))
    oval_opposite_tangent_corrected = False
    oval_monotonic_divergence_applied = False

    # In oval mode force each tangent through the observed lateral point at the widest row.
    if bool(oval.get("detected", False)) and oval.get("peak_y") is not None:
        py = int(oval["peak_y"])
        obs_rows = np.asarray(oval.get("rows", observed_y), dtype=float)
        obs_lsm = np.asarray(oval.get("left_smooth", observed_left), dtype=float)
        obs_rsm = np.asarray(oval.get("right_smooth", observed_right), dtype=float)
        if obs_rows.size >= 2:
            la = float(np.interp(py, obs_rows, obs_lsm))
            ra = float(np.interp(py, obs_rows, obs_rsm))
            left_fit["slope"] = float(np.clip(left_fit["slope"], -max_slope, max_slope))
            right_fit["slope"] = float(np.clip(right_fit["slope"], -max_slope, max_slope))
            left_fit["intercept"] = float(la - left_fit["slope"] * py)
            right_fit["intercept"] = float(ra - right_fit["slope"] * py)

            # Independent fits can inherit a common tilt of an oval lumen and make both
            # projected tangents point in the same direction.  That shifts and narrows the
            # ME corridor.  Keep each side's fitted magnitude, but enforce ventrolateral
            # signs: the left tangent opens left and the right tangent opens right.
            if bool(getattr(args, "wall_oval_enforce_opposite_tangents", True)):
                min_out = float(max(0.0, getattr(args, "wall_oval_min_outward_slope", 0.02)))
                old_ls = float(left_fit.get("slope", 0.0))
                old_rs = float(right_fit.get("slope", 0.0))
                left_mag = float(min(max_slope, max(min_out, abs(old_ls))))
                right_mag = float(min(max_slope, max(min_out, abs(old_rs))))
                left_fit["slope"] = -left_mag
                right_fit["slope"] = right_mag
                left_fit["intercept"] = float(la - left_fit["slope"] * py)
                right_fit["intercept"] = float(ra - right_fit["slope"] * py)
                oval_opposite_tangent_corrected = bool(old_ls >= 0.0 or old_rs <= 0.0 or old_ls * old_rs >= 0.0)
                if oval_opposite_tangent_corrected:
                    fit_status += "_opposite_tangents_corrected"

    median_half_gap = float(np.nanmedian(np.maximum(1.0, (right[profile_rows] - left[profile_rows]) / 2.0))) if profile_rows.size else max(10.0, W * 0.03)
    if not np.isfinite(left_fit.get("intercept", np.nan)):
        left_fit["slope"] = 0.0
        left_fit["intercept"] = float(mid_x - median_half_gap)
    if not np.isfinite(right_fit.get("intercept", np.nan)):
        right_fit["slope"] = 0.0
        right_fit["intercept"] = float(mid_x + median_half_gap)

    top_margin = um_to_px(float(getattr(args, "wall_projection_top_margin_um", 12.0)), args, min_px=0)
    bottom_margin = um_to_px(float(getattr(args, "wall_projection_bottom_margin_um", 15.0)), args, min_px=0)
    ventral_tip_y = int(max(line_y_start, min(y_max, run_end)))
    presence_y_start = int(max(line_y_start, min(y_max, ventral_tip_y - top_margin)))
    projection_y_start = int(line_y_start)
    projection_y_end = int(max(projection_y_start, min(H - 1, y_max + bottom_margin)))
    line_y = np.arange(projection_y_start, projection_y_end + 1, dtype=int)

    lraw = float(left_fit["slope"]) * line_y.astype(float) + float(left_fit["intercept"])
    rraw = float(right_fit["slope"]) * line_y.astype(float) + float(right_fit["intercept"])
    signed_margin = float(getattr(args, "wall_corridor_margin_um", 0.0)) / max(float(getattr(args, "um_per_px", DEFAULT_UM_PER_PX)), 1e-6)
    lraw = lraw - signed_margin
    rraw = rraw + signed_margin
    if bool(oval.get("detected", False)) and bool(getattr(args, "wall_oval_monotonic_divergence", True)) and lraw.size:
        # Below an oval tangent point the corridor may stay vertical or widen, but it may
        # not reconverge.  This also protects against later center/width clipping.
        lraw = np.minimum.accumulate(lraw)
        rraw = np.maximum.accumulate(rraw)
        oval_monotonic_divergence_applied = True

    min_width = float(max(2, um_to_px(float(getattr(args, "wall_corridor_min_width_um", 20.0)), args, min_px=2)))
    max_width = float(max(min_width + 2.0, float(getattr(args, "wall_corridor_max_width_frac", 0.70)) * W))
    center_limit = float(max(um_to_px(100.0, args, min_px=20), 0.10 * W))
    left_x = np.zeros_like(lraw, dtype=float)
    right_x = np.zeros_like(rraw, dtype=float)
    for i, (lx, rx) in enumerate(zip(lraw, rraw)):
        if not np.isfinite(lx) or not np.isfinite(rx):
            lx, rx = mid_x - median_half_gap, mid_x + median_half_gap
        center = float(np.clip((lx + rx) / 2.0, mid_x - center_limit, mid_x + center_limit))
        width = float(np.clip(rx - lx, min_width, max_width))
        left_x[i] = float(np.clip(center - width / 2.0, 0, W - 2))
        right_x[i] = float(np.clip(center + width / 2.0, 1, W - 1))
        if right_x[i] <= left_x[i] + 1:
            right_x[i] = min(W - 1.0, left_x[i] + max(2.0, min_width))

    if bool(oval.get("detected", False)) and bool(getattr(args, "wall_oval_monotonic_divergence", True)) and left_x.size:
        left_x = np.minimum.accumulate(left_x)
        right_x = np.maximum.accumulate(right_x)
        left_x = np.clip(left_x, 0.0, W - 2.0)
        right_x = np.clip(right_x, 1.0, W - 1.0)
        for i in range(left_x.size):
            width = float(right_x[i] - left_x[i])
            if width < min_width:
                center = float((left_x[i] + right_x[i]) / 2.0)
                left_x[i] = float(np.clip(center - min_width / 2.0, 0.0, W - 2.0))
                right_x[i] = float(np.clip(center + min_width / 2.0, 1.0, W - 1.0))
            elif width > max_width:
                center = float((left_x[i] + right_x[i]) / 2.0)
                left_x[i] = float(np.clip(center - max_width / 2.0, 0.0, W - 2.0))
                right_x[i] = float(np.clip(center + max_width / 2.0, 1.0, W - 1.0))

    corridor = np.zeros((H, W), dtype=bool)
    for y, lx, rx in zip(line_y, left_x, right_x):
        x0 = int(max(0, min(W - 1, math.ceil(float(lx)))))
        x1 = int(max(0, min(W - 1, math.floor(float(rx)))))
        if x1 >= x0:
            corridor[int(y), x0:x1 + 1] = True

    _region_debug(
        args,
        "Wall fit: status=%s oval=%s peak_y=%s run=%s-%s fit_rows=%s-%s n=%s left(x=%.4fy%+.1f rmse=%.2f) right(x=%.4fy%+.1f rmse=%.2f) projection=%s-%s",
        fit_status, bool(oval.get("detected", False)), oval.get("peak_y"), int(run_start), int(run_end),
        int(fit_rows.min()) if fit_rows.size else -1, int(fit_rows.max()) if fit_rows.size else -1,
        int(fit_rows.size), float(left_fit["slope"]), float(left_fit["intercept"]),
        float(left_fit.get("rmse", np.nan)), float(right_fit["slope"]), float(right_fit["intercept"]),
        float(right_fit.get("rmse", np.nan)), int(projection_y_start), int(projection_y_end),
    )

    return {
        "mid_x": float(mid_x),
        "valid_rows": valid,
        "main_run_start_y": int(run_start),
        "main_run_end_y": int(run_end),
        "fit_row_start_y": int(fit_rows.min()) if fit_rows.size else int(run_start),
        "fit_row_end_y": int(fit_rows.max()) if fit_rows.size else int(run_end),
        "fit_rows_n": int(fit_rows.size),
        "line_y_start": int(line_y_start),
        "ventral_tip_y": int(ventral_tip_y),
        "presence_y_start": int(presence_y_start),
        "projection_y_start": int(projection_y_start),
        "projection_y_end": int(projection_y_end),
        "left_slope": float(left_fit["slope"]),
        "left_intercept": float(left_fit["intercept"]),
        "left_rmse": float(left_fit.get("rmse", np.nan)),
        "left_fit_n": int(left_fit.get("n", 0)),
        "right_slope": float(right_fit["slope"]),
        "right_intercept": float(right_fit["intercept"]),
        "right_rmse": float(right_fit.get("rmse", np.nan)),
        "right_fit_n": int(right_fit.get("n", 0)),
        "fit_status": str(fit_status),
        "line_y": line_y,
        "left_x": left_x,
        "right_x": right_x,
        "corridor": corridor,
        "median_half_gap_px": float(median_half_gap),
        "valid_row_count": int(np.count_nonzero(valid)),
        "observed_y": np.asarray(oval.get("rows", observed_y), dtype=int),
        "observed_left_x": np.asarray(oval.get("left_smooth", observed_left), dtype=float),
        "observed_right_x": np.asarray(oval.get("right_smooth", observed_right), dtype=float),
        "oval_detected": bool(oval.get("detected", False)),
        "oval_mode": str(oval.get("mode", "auto")),
        "oval_peak_y": oval.get("peak_y"),
        "oval_peak_width_px": float(oval.get("peak_width_px", np.nan)),
        "oval_top_width_px": float(oval.get("top_width_px", np.nan)),
        "oval_bottom_width_px": float(oval.get("bottom_width_px", np.nan)),
        "oval_bulge_ratio_top": float(oval.get("bulge_ratio_top", np.nan)),
        "oval_bulge_ratio_bottom": float(oval.get("bulge_ratio_bottom", np.nan)),
        "oval_peak_fraction": float(oval.get("peak_fraction", np.nan)),
        "oval_left_slope_raw": float(oval_left_slope_raw),
        "oval_right_slope_raw": float(oval_right_slope_raw),
        "oval_opposite_tangent_corrected": bool(oval_opposite_tangent_corrected),
        "oval_monotonic_divergence_applied": bool(oval_monotonic_divergence_applied),
        "oval_projection_width_start_px": float(right_x[0] - left_x[0]) if right_x.size else float("nan"),
        "oval_projection_width_end_px": float(right_x[-1] - left_x[-1]) if right_x.size else float("nan"),
    }


def _me_mask_from_projected_walls(
    dapi_lbl: Any,
    tissue: Any,
    wall: Dict[str, Any],
    args: argparse.Namespace,
) -> Tuple[Any, Dict[str, Any]]:
    """Classify ME as the closed ventricular floor, not the ventricular sidewalls.

    The projected wall corridor is retained, but it is only a lateral constraint.  The
    central empty lumen is followed to its ventral boundary and a curved floor envelope
    is estimated.  A DAPI nucleus is ME only when its centroid is below that local floor
    curve, between the projected walls, and not in the near-wall exclusion band.  Thus
    sidewall nuclei are heuristically reassigned to ARC.
    """
    dapi = np.asarray(dapi_lbl).astype(np.int32, copy=False)
    t = np.asarray(tissue, dtype=bool)
    H, W = dapi.shape
    corridor = np.asarray(wall.get("corridor", np.zeros((H, W), dtype=bool)), dtype=bool)
    mid_x = float(wall.get("mid_x", W / 2.0))
    floor_info = _infer_closed_ventricle_floor_curve(t, wall, args, dapi_lbl=dapi)

    all_labels: List[int] = []
    floor_labels: List[int] = []
    floor_centroids: List[Tuple[float, float, float]] = []
    sidewall_reassigned = 0
    cx = np.asarray(floor_info.get("curve_x", []), dtype=float)
    cy = np.asarray(floor_info.get("curve_y", []), dtype=float)
    for rp in regionprops(dapi):
        y, x = float(rp.centroid[0]), float(rp.centroid[1])
        yi = int(max(0, min(H - 1, round(y))))
        xi = int(max(0, min(W - 1, round(x))))
        if bool(corridor[yi, xi]):
            all_labels.append(int(rp.label))
            if _point_is_me_floor(wall, floor_info, y, x, args):
                fy = float(np.interp(x, cx, cy)) if cx.size >= 2 else float("nan")
                floor_labels.append(int(rp.label))
                floor_centroids.append((y, x, fy))
            else:
                sidewall_reassigned += 1

    n_all = int(len(all_labels))
    n_floor = int(len(floor_labels))
    half_mid_px = um_to_px(float(getattr(args, "wall_me_midline_half_width_um", 140.0)), args, min_px=1)
    midline_n = int(sum(1 for _y, x, _fy in floor_centroids if abs(float(x) - mid_x) <= half_mid_px))
    if floor_centroids:
        depths = np.asarray([max(0.0, y - fy) for y, _x, fy in floor_centroids if np.isfinite(fy)], dtype=float)
        ys_v = np.asarray([p[0] for p in floor_centroids], dtype=float)
        xs_v = np.asarray([p[1] for p in floor_centroids], dtype=float)
        depth_px = float(np.nanmax(depths)) if depths.size else 0.0
        span_px = float(np.nanmax(xs_v) - np.nanmin(xs_v)) if xs_v.size > 1 else 0.0
        centroid_top_y = int(math.floor(float(np.nanmin(ys_v))))
        centroid_bottom_y = int(math.ceil(float(np.nanmax(ys_v))))
    else:
        depth_px = 0.0
        span_px = 0.0
        centroid_top_y = None
        centroid_bottom_y = None

    mode = str(getattr(args, "me_detection_mode", "auto")).lower()
    min_n = int(max(0, getattr(args, "wall_me_min_nuclei", 5)))
    min_mid = int(max(0, getattr(args, "wall_me_min_midline_nuclei", 1)))
    min_depth_px = float(um_to_px(float(getattr(args, "wall_me_min_depth_um", 18.0)), args, min_px=0))
    floor_geometry_ok = bool(floor_info.get("floor_found", False)) and not bool(floor_info.get("open_bottom", False))
    sparse_rescue_used = bool(floor_info.get("sparse_floor_rescue", False))
    sparse_rescue_n = int(floor_info.get("sparse_floor_candidate_nuclei", 0) or 0)
    sparse_rescue_min_n = int(max(2, getattr(args, "wall_floor_sparse_min_nuclei", 3)))
    sparse_bilateral_ok = bool(
        int(floor_info.get("sparse_floor_left_nuclei", 0) or 0) >= int(max(0, getattr(args, "wall_floor_sparse_min_bilateral_nuclei", 1)))
        and int(floor_info.get("sparse_floor_right_nuclei", 0) or 0) >= int(max(0, getattr(args, "wall_floor_sparse_min_bilateral_nuclei", 1)))
    )
    if mode == "never":
        present = False
    elif mode == "always":
        present = bool(floor_geometry_ok and n_floor > 0)
    elif sparse_rescue_used:
        # The rescue has already passed bilateral, span, low-density and depth checks.
        # Do not re-impose the normal five-nucleus floor threshold that caused the
        # sparse-but-real ME example to be discarded.
        present = bool(
            floor_geometry_ok
            and n_floor > 0
            and sparse_rescue_n >= sparse_rescue_min_n
            and sparse_bilateral_ok
        )
    else:
        present = bool(
            floor_geometry_ok
            and n_floor >= min_n
            and midline_n >= min_mid
            and (depth_px >= min_depth_px or n_floor >= max(min_n * 2, min_n + 5))
        )

    floor_nuc_mask = np.isin(dapi, np.asarray(floor_labels, dtype=np.int32)) if floor_labels else np.zeros((H, W), dtype=bool)
    me_mask = np.zeros((H, W), dtype=bool)
    if present and cx.size >= 2:
        floor_by_x = np.full(W, np.nan, dtype=float)
        x0 = max(0, int(math.floor(float(np.nanmin(cx)))))
        x1 = min(W - 1, int(math.ceil(float(np.nanmax(cx)))))
        xxv = np.arange(x0, x1 + 1, dtype=float)
        floor_by_x[x0:x1 + 1] = np.interp(xxv, cx, cy)
        yy_grid = np.arange(H, dtype=float)[:, None]
        floor_gate = yy_grid >= (floor_by_x[None, :] + float(um_to_px(float(getattr(args, "wall_floor_centroid_margin_um", 0.0)), args, min_px=0)))
        floor_gate[:, ~np.isfinite(floor_by_x)] = False
        candidate = corridor & t & floor_gate

        # Exclude sidewall tissue near the projected lines unless it lies well below the
        # local floor curve. This is the explicit ARC reassignment requested by the user.
        line_ys = np.asarray(wall.get("line_y", []), dtype=float)
        line_l = np.asarray(wall.get("left_x", []), dtype=float)
        line_r = np.asarray(wall.get("right_x", []), dtype=float)
        good = np.isfinite(line_ys) & np.isfinite(line_l) & np.isfinite(line_r)
        if int(np.count_nonzero(good)) >= 2:
            line_ys = line_ys[good]; line_l = line_l[good]; line_r = line_r[good]
            order = np.argsort(line_ys); line_ys = line_ys[order]; line_l = line_l[order]; line_r = line_r[order]
            rows = np.arange(H, dtype=float)
            lrow = np.interp(rows, line_ys, line_l, left=np.nan, right=np.nan)
            rrow = np.interp(rows, line_ys, line_r, left=np.nan, right=np.nan)
            xx_grid = np.arange(W, dtype=float)[None, :]
            side_dist = np.minimum(xx_grid - lrow[:, None], rrow[:, None] - xx_grid)
            side_excl = float(um_to_px(float(getattr(args, "wall_sidewall_exclusion_um", 18.0)), args, min_px=0))
            release = float(um_to_px(float(getattr(args, "wall_sidewall_release_depth_um", 35.0)), args, min_px=0))
            shallow = yy_grid < (floor_by_x[None, :] + release)
            shallow[:, ~np.isfinite(floor_by_x)] = False
            candidate &= ~((side_dist < side_excl) & shallow)

        try:
            cc = label(candidate.astype(bool), connectivity=1)
            keep_ids = np.unique(cc[floor_nuc_mask])
            keep_ids = keep_ids[keep_ids > 0]
            me_mask = np.isin(cc, keep_ids) if keep_ids.size else candidate.copy()
        except Exception:
            me_mask = candidate.copy()
        dil = um_to_px(float(getattr(args, "wall_me_tissue_dilate_um", 10.0)), args, min_px=0)
        close = um_to_px(float(getattr(args, "wall_me_tissue_close_um", 12.0)), args, min_px=0)
        try:
            if dil > 0:
                me_mask = binary_dilation((me_mask | floor_nuc_mask).astype(bool), disk(max(1, dil)))
            else:
                me_mask |= floor_nuc_mask
            if close > 0:
                me_mask = binary_closing(me_mask.astype(bool), disk(max(1, close)))
            me_mask = remove_small_holes(me_mask.astype(bool), area_threshold=max(32, int(0.00015 * H * W)))
        except Exception:
            me_mask |= floor_nuc_mask
        # Never let display smoothing climb back into the lumen sidewalls.
        me_mask &= (floor_gate & corridor) | floor_nuc_mask
        me_mask |= floor_nuc_mask

    _region_debug(
        args,
        "ME floor-only: present=%s mode=%s floor_status=%s open_bottom=%s corridor_nuclei=%s floor_nuclei=%s sidewall_to_ARC=%s midline=%s depth_px=%.1f span_px=%.1f",
        bool(present), mode, floor_info.get("status"), bool(floor_info.get("open_bottom", False)),
        n_all, n_floor, sidewall_reassigned, midline_n, depth_px, span_px,
    )
    return me_mask.astype(bool), {
        "me_floor_present": bool(present),
        "me_present": bool(present),
        "me_floor_method": ("sparse_bilateral_floor_rescue_v17" if sparse_rescue_used else "closed_curved_ventricle_floor_excluding_sidewalls_v17"),
        "me_floor_geometry_status": str(floor_info.get("status", "NA")),
        "me_floor_open_bottom": bool(floor_info.get("open_bottom", False)),
        "me_floor_sparse_rescue": bool(sparse_rescue_used),
        "me_floor_primary_status": str(floor_info.get("primary_floor_status", floor_info.get("status", "NA"))),
        "me_floor_primary_open_bottom": bool(floor_info.get("primary_open_bottom", floor_info.get("open_bottom", False))),
        "me_floor_sparse_candidate_nuclei": int(sparse_rescue_n),
        "me_floor_sparse_left_nuclei": int(floor_info.get("sparse_floor_left_nuclei", 0) or 0),
        "me_floor_sparse_right_nuclei": int(floor_info.get("sparse_floor_right_nuclei", 0) or 0),
        "me_floor_sparse_midline_nuclei": int(floor_info.get("sparse_floor_midline_nuclei", 0) or 0),
        "me_floor_sparse_span_px": float(floor_info.get("sparse_floor_span_px", 0.0) or 0.0),
        "me_floor_sparse_depth_px": float(floor_info.get("sparse_floor_depth_px", 0.0) or 0.0),
        "me_floor_sparse_fraction": float(floor_info.get("sparse_floor_sparse_fraction", 0.0) or 0.0),
        "me_floor_sparse_density_ratio": float(floor_info.get("sparse_floor_density_ratio", np.nan)),
        "me_wall_enclosed_nuclei": int(n_all),
        "me_wall_floor_nuclei": int(n_floor),
        "me_wall_enclosed_ventral_nuclei": int(n_floor),
        "me_wall_sidewall_reassigned_nuclei": int(sidewall_reassigned),
        "me_wall_midline_nuclei": int(midline_n),
        "me_wall_depth_px": float(depth_px),
        "me_wall_span_px": float(span_px),
        "me_wall_centroid_top_y": centroid_top_y,
        "me_wall_centroid_bottom_y": centroid_bottom_y,
        "me_wall_line_y_start": int(wall.get("line_y_start", wall.get("projection_y_start", 0))),
        "me_wall_ventral_tip_y": int(wall.get("ventral_tip_y", wall.get("main_run_end_y", 0))),
        "me_wall_presence_y_start": int(wall.get("presence_y_start", wall.get("ventral_tip_y", 0))),
        "me_floor_top_y": int(math.floor(float(np.nanmin(cy)))) if cy.size else None,
        "me_top_y": int(math.floor(float(np.nanmin(cy)))) if cy.size else None,
        "me_wall_min_nuclei": int(min_n),
        "me_wall_min_midline_nuclei": int(min_mid),
        "me_wall_min_depth_px": float(min_depth_px),
        "me_wall_corridor_area_px": int(np.count_nonzero(corridor)),
        "me_wall_tissue_area_px": int(np.count_nonzero(me_mask)),
        "me_floor_bottom_y": float(floor_info.get("floor_bottom_y", np.nan)),
        "me_floor_curve_columns": int(floor_info.get("curve_columns", 0)),
        "me_sidewall_exclusion_px": int(um_to_px(float(getattr(args, "wall_sidewall_exclusion_um", 18.0)), args, min_px=0)),
        "me_sidewall_release_depth_px": int(um_to_px(float(getattr(args, "wall_sidewall_release_depth_um", 35.0)), args, min_px=0)),
        "me_floor_centroid_margin_px": int(um_to_px(float(getattr(args, "wall_floor_centroid_margin_um", 0.0)), args, min_px=0)),
        "_me_floor_curve_x": np.asarray(floor_info.get("curve_x", []), dtype=float),
        "_me_floor_curve_y": np.asarray(floor_info.get("curve_y", []), dtype=float),
        "_me_floor_curve_y_raw": np.asarray(floor_info.get("curve_y_raw", []), dtype=float),
    }


def infer_arc_me_region_mask(
    dapi_lbl: Any,
    npy_lbl: Optional[Any],
    pomc_lbl: Optional[Any],
    args: argparse.Namespace,
    cfos_lbl: Optional[Any] = None,
) -> Tuple[Any, Dict[str, Any]]:
    """Infer ARC/ME from projected third-ventricle walls.

    Anatomical rule implemented in v17:
      * robustly follow the left and right walls of the third ventricle;
      * for an oval lumen, project tangents from its widest lateral points;
      * detect the closed curved ventral floor of the lumen;
      * only nuclei below that floor and between the projected walls are ME;
      * ventricular sidewall nuclei and every remaining nucleus are ARC;
      * if the lumen is open ventrally, ME is absent.

    c-FOS, POMC and NPY are deliberately not used to place the ARC/ME border, avoiding
    treatment-dependent circularity. They are quantified only after the anatomical mask
    has been fixed.
    """
    dapi = np.asarray(dapi_lbl).astype(np.int32, copy=False)
    H, W = dapi.shape
    region = np.zeros((H, W), dtype=np.uint8)
    tissue = smooth_tissue_mask_from_dapi(dapi, args)
    y_min, y_max, x_min, x_max = tissue_bounds(tissue)
    vent = detect_third_ventricle_edges(tissue, args)
    wall = fit_projected_ventricle_walls(tissue, vent, args)
    me_mask, me_info = _me_mask_from_projected_walls(dapi, tissue, wall, args)
    me_present = bool(me_info.get("me_floor_present", False))

    arc_support = _arc_tissue_fill_mask(dapi, tissue, args, y_min, y_max, x_min, x_max)
    arc_support |= (dapi > 0)
    region[arc_support] = REGION_NAME_TO_CODE["ARC"]
    if me_present:
        region[me_mask] = REGION_NAME_TO_CODE["ME"]
        point_info = dict(wall)
        point_info.update(me_info)
        # region_name_from_projected_walls consumes the stable public geometry
        # keys used downstream. fit_projected_ventricle_walls exposes the same
        # arrays under its internal line_y/left_x/right_x names; bridge them
        # here before classifying DAPI centroids. Without this bridge every
        # nucleus silently falls back to ARC while later microglia calls can
        # still vote ME from the renamed output keys.
        point_info.update({
            "_wall_y": np.asarray(wall.get("line_y", []), dtype=float),
            "_wall_left_x": np.asarray(wall.get("left_x", []), dtype=float),
            "_wall_right_x": np.asarray(wall.get("right_x", []), dtype=float),
        })
        for rp in regionprops(dapi):
            y, x = float(rp.centroid[0]), float(rp.centroid[1])
            if region_name_from_projected_walls(point_info, y, x) == "ME":
                region[dapi == int(rp.label)] = REGION_NAME_TO_CODE["ME"]
            else:
                region[dapi == int(rp.label)] = REGION_NAME_TO_CODE["ARC"]

    pad_x = um_to_px(80.0, args, min_px=20)
    pad_y = um_to_px(40.0, args, min_px=10)
    bbox = np.zeros_like(region, dtype=bool)
    bbox[max(0, y_min - pad_y):min(H, y_max + pad_y + 1), max(0, x_min - pad_x):min(W, x_max + pad_x + 1)] = True
    region[(~bbox) & (dapi == 0)] = 0

    line_y = np.asarray(wall.get("line_y", []), dtype=int)
    left_x = np.asarray(wall.get("left_x", []), dtype=float)
    right_x = np.asarray(wall.get("right_x", []), dtype=float)
    info: Dict[str, Any] = {
        "region_mode": "arc-me-walls",
        "region_ai_model": "oval_aware_projected_walls_floor_only_ARC_ME_v17",
        "region_rule": "ME=closed_floor_below_curved_lumen_boundary_between_projected_walls; ventricular_sidewalls=ARC; all_remaining_tissue=ARC",
        "mid_x": float(wall.get("mid_x", vent.get("mid_x", W / 2.0))),
        "ventricle_coverage": float(vent.get("coverage", 0.0)),
        "tissue_y_min": int(y_min),
        "tissue_y_max": int(y_max),
        "tissue_x_min": int(x_min),
        "tissue_x_max": int(x_max),
        "arc_top_y": int(y_min),
        "me_top_y": int(wall.get("projection_y_start", y_max)),
        "arc_shape_selected": "all_tissue_except_closed_floor_ME; ventricular_sidewalls_forced_to_ARC",
        "arc_area_px": int(np.count_nonzero(region == REGION_NAME_TO_CODE["ARC"])),
        "me_area_px": int(np.count_nonzero(region == REGION_NAME_TO_CODE["ME"])),
        "wall_fit_status": str(wall.get("fit_status", "NA")),
        "wall_main_run_start_y": int(wall.get("main_run_start_y", y_min)),
        "wall_main_run_end_y": int(wall.get("main_run_end_y", y_max)),
        "wall_fit_row_start_y": int(wall.get("fit_row_start_y", y_min)),
        "wall_fit_row_end_y": int(wall.get("fit_row_end_y", y_max)),
        "wall_fit_rows_n": int(wall.get("fit_rows_n", 0)),
        "wall_line_y_start": int(wall.get("line_y_start", wall.get("projection_y_start", y_min))),
        "wall_ventral_tip_y": int(wall.get("ventral_tip_y", wall.get("main_run_end_y", y_max))),
        "wall_presence_y_start": int(wall.get("presence_y_start", wall.get("main_run_end_y", y_max))),
        "wall_projection_y_start": int(wall.get("projection_y_start", y_min)),
        "wall_projection_y_end": int(wall.get("projection_y_end", y_max)),
        "wall_left_slope": float(wall.get("left_slope", np.nan)),
        "wall_left_intercept": float(wall.get("left_intercept", np.nan)),
        "wall_left_rmse": float(wall.get("left_rmse", np.nan)),
        "wall_left_fit_n": int(wall.get("left_fit_n", 0)),
        "wall_right_slope": float(wall.get("right_slope", np.nan)),
        "wall_right_intercept": float(wall.get("right_intercept", np.nan)),
        "wall_right_rmse": float(wall.get("right_rmse", np.nan)),
        "wall_right_fit_n": int(wall.get("right_fit_n", 0)),
        "wall_valid_row_count": int(wall.get("valid_row_count", 0)),
        "wall_median_half_gap_px": float(wall.get("median_half_gap_px", np.nan)),
        "wall_oval_detected": bool(wall.get("oval_detected", False)),
        "wall_oval_mode": str(wall.get("oval_mode", "auto")),
        "wall_oval_peak_y": wall.get("oval_peak_y"),
        "wall_oval_peak_width_px": float(wall.get("oval_peak_width_px", np.nan)),
        "wall_oval_bulge_ratio_top": float(wall.get("oval_bulge_ratio_top", np.nan)),
        "wall_oval_bulge_ratio_bottom": float(wall.get("oval_bulge_ratio_bottom", np.nan)),
        "wall_oval_left_slope_raw": float(wall.get("oval_left_slope_raw", np.nan)),
        "wall_oval_right_slope_raw": float(wall.get("oval_right_slope_raw", np.nan)),
        "wall_oval_opposite_tangent_corrected": bool(wall.get("oval_opposite_tangent_corrected", False)),
        "wall_oval_monotonic_divergence_applied": bool(wall.get("oval_monotonic_divergence_applied", False)),
        "wall_oval_projection_width_start_px": float(wall.get("oval_projection_width_start_px", np.nan)),
        "wall_oval_projection_width_end_px": float(wall.get("oval_projection_width_end_px", np.nan)),
        "_wall_observed_y": np.asarray(wall.get("observed_y", []), dtype=int),
        "_wall_observed_left_x": np.asarray(wall.get("observed_left_x", []), dtype=float),
        "_wall_observed_right_x": np.asarray(wall.get("observed_right_x", []), dtype=float),
        "_wall_y": line_y,
        "_wall_left_x": left_x,
        "_wall_right_x": right_x,
    }
    info.update(me_info)
    return region.astype(np.uint8), info

def _normalize_feature(arr: Any) -> Any:
    x = np.asarray(arr, dtype=np.float32)
    if x.size == 0:
        return x
    lo = float(np.nanpercentile(x, 1))
    hi = float(np.nanpercentile(x, 99))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        m = float(np.nanmax(x)) if np.isfinite(np.nanmax(x)) else 1.0
        return np.clip(x / max(m, 1e-6), 0, 1).astype(np.float32)
    return np.clip((x - lo) / max(hi - lo, 1e-6), 0, 1).astype(np.float32)


def _density_feature(lbl: Any, sigma: float = 6.0) -> Any:
    try:
        from scipy.ndimage import gaussian_filter
        arr = gaussian_filter((np.asarray(lbl) > 0).astype(np.float32), sigma=float(sigma))
        return _normalize_feature(arr)
    except Exception:
        return (np.asarray(lbl) > 0).astype(np.float32)


def _resize_float(arr: Any, shape: Tuple[int, int], order: int = 1) -> Any:
    from skimage.transform import resize
    return resize(np.asarray(arr), shape, order=order, mode="reflect", preserve_range=True, anti_aliasing=(order != 0)).astype(np.float32)


def _postprocess_arc_me_prediction(pred: Any, prior: Any, info: Dict[str, Any], args: argparse.Namespace) -> Any:
    pred = np.asarray(pred, dtype=np.uint8)
    prior = np.asarray(prior, dtype=np.uint8)
    H, W = prior.shape
    out = np.zeros((H, W), dtype=np.uint8)
    min_area = max(64, int(0.002 * H * W))
    for name in ["ARC", "ME"]:
        code = REGION_NAME_TO_CODE[name]
        m = pred == code
        try:
            m = remove_small_objects(m.astype(bool), min_size=min_area)
            m = remove_small_holes(m.astype(bool), area_threshold=max(64, min_area // 2))
        except Exception:
            pass
        out[m] = code
    # Safety: if neural prediction erases a class, recover the anatomical prior.
    for name in ["ARC", "ME"]:
        code = REGION_NAME_TO_CODE[name]
        prior_area = int(np.count_nonzero(prior == code))
        out_area = int(np.count_nonzero(out == code))
        if prior_area > 0 and out_area < max(50, int(0.20 * prior_area)):
            out[prior == code] = code
    # Preserve the projected-wall geometry.  No horizontal ARC/ME line is created.
    # If the wall corridor contains no convincing tissue, the prior has no ME and the
    # neural refinement cannot invent one.
    me_present = bool(info.get("me_floor_present", True))
    if bool(getattr(args, "region_dl_preserve_geometry", True)):
        outside = (prior == 0)
        out[outside] = 0
        missing = (prior > 0) & (out == 0)
        out[missing] = prior[missing]
    # v17 provides the closed floor-only prior. In strict mode the optional U-Net may
    # expand ME only inside a small, explicit boundary band; it cannot escape that enclosure.
    if me_present:
        me_code = REGION_NAME_TO_CODE["ME"]
        arc_code = REGION_NAME_TO_CODE["ARC"]
        if bool(getattr(args, "region_dl_strict_me_prior", True)):
            band_um = float(max(0.0, getattr(args, "region_dl_boundary_refine_um", 0.0)))
            band_px = um_to_px(band_um, args, min_px=0)
            if band_px > 0:
                allowed_me = binary_dilation(prior == me_code, disk(max(1, band_px))) & (prior > 0)
                out[(out == me_code) & (~allowed_me)] = arc_code
            else:
                out[(out == me_code) & (prior != me_code)] = arc_code
            info["region_dl_boundary_refine_um"] = float(band_um)
            info["region_dl_boundary_refine_px"] = int(band_px)
        out[prior == me_code] = me_code
    else:
        out[out == REGION_NAME_TO_CODE["ME"]] = REGION_NAME_TO_CODE["ARC"]
        out[prior == REGION_NAME_TO_CODE["ARC"]] = REGION_NAME_TO_CODE["ARC"]
    return out.astype(np.uint8)


def infer_arc_me_region_mask_deeplearning(dapi_lbl: Any, npy_lbl: Optional[Any], pomc_lbl: Optional[Any], args: argparse.Namespace, cfos_lbl: Optional[Any] = None) -> Tuple[Any, Dict[str, Any]]:
    """
    Deep-learning assisted ARC/ME segmentation.

    This optionally trains a small U-Net on-the-fly from the deterministic projected-
    ventricular-wall prior.  Only DAPI-derived anatomy, coordinates, and prior channels
    are used; c-FOS/NPY/POMC are deliberately excluded so treatment-dependent signal can
    never move the ARC/ME border.  DAPI nuclei are ultimately classified by the exact wall
    enclosure rule even when this optional visual-mask refinement is enabled.
    """
    prior, info = infer_arc_me_region_mask(dapi_lbl, npy_lbl, pomc_lbl, args, cfos_lbl=cfos_lbl)
    info = dict(info)
    info["region_mode"] = "arc-me-dl"
    info["region_ai_model"] = "ventricle_wall_floor_only_constrained_tiny_unet_ARC_ME_v17"
    epochs = int(max(0, getattr(args, "region_dl_epochs", 5)))
    if epochs <= 0:
        info["region_ai_model"] = "projected_ventricle_walls_prior_no_neural_refinement_epochs0"
        return prior, info

    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
    except Exception as e:
        info["region_ai_model"] = f"projected_ventricle_walls_prior_torch_unavailable_fallback: {e}"
        return prior, info

    dapi = np.asarray(dapi_lbl)
    H, W = dapi.shape
    max_side = int(max(128, getattr(args, "region_dl_size", 384)))
    scale = min(1.0, max_side / float(max(H, W)))
    h2 = int(max(64, round(H * scale)))
    w2 = int(max(64, round(W * scale)))

    # Feature stack is anatomy-only.  Do not leak c-FOS, NPY, or POMC into region
    # placement: those markers are outcomes to quantify, not anatomical training labels.
    yy, xx = np.mgrid[0:H, 0:W]
    mid_x = float(info.get("mid_x", W / 2.0))
    y_norm = yy.astype(np.float32) / max(1.0, float(H - 1))
    x_norm = xx.astype(np.float32) / max(1.0, float(W - 1))
    x_mid_dist = np.abs(xx.astype(np.float32) - mid_x) / max(1.0, float(W) / 2.0)
    tissue = smooth_tissue_mask_from_dapi(dapi, args).astype(np.float32)
    dens_sigma = um_to_px(24.0, args, min_px=3)
    dapi_density = _normalize_feature(_density_map_from_label_centroids(dapi, sigma_px=dens_sigma))
    sparse_density = (1.0 - dapi_density).astype(np.float32)
    features = [
        dapi_density,
        sparse_density,
        tissue,
        y_norm,
        x_norm,
        x_mid_dist.astype(np.float32),
        (prior > 0).astype(np.float32),
        (prior == REGION_NAME_TO_CODE["ARC"]).astype(np.float32),
        (prior == REGION_NAME_TO_CODE["ME"]).astype(np.float32),
    ]
    X = np.stack([_resize_float(f, (h2, w2), order=1) for f in features], axis=0).astype(np.float32)
    y = _resize_float(prior, (h2, w2), order=0).round().astype(np.int64)

    class TinyUNet(nn.Module):
        def __init__(self, in_ch: int, base: int = 16, n_cls: int = 3):
            super().__init__()
            base = int(max(4, base))
            self.c1 = nn.Sequential(nn.Conv2d(in_ch, base, 3, padding=1), nn.ReLU(inplace=True), nn.Conv2d(base, base, 3, padding=1), nn.ReLU(inplace=True))
            self.c2 = nn.Sequential(nn.Conv2d(base, base * 2, 3, padding=1), nn.ReLU(inplace=True), nn.Conv2d(base * 2, base * 2, 3, padding=1), nn.ReLU(inplace=True))
            self.c3 = nn.Sequential(nn.Conv2d(base * 3, base, 3, padding=1), nn.ReLU(inplace=True), nn.Conv2d(base, base, 3, padding=1), nn.ReLU(inplace=True))
            self.out = nn.Conv2d(base, n_cls, 1)
        def forward(self, x):
            e1 = self.c1(x)
            p = F.max_pool2d(e1, 2)
            e2 = self.c2(p)
            u = F.interpolate(e2, size=e1.shape[-2:], mode="bilinear", align_corners=False)
            return self.out(self.c3(torch.cat([u, e1], dim=1)))

    dev_arg = str(getattr(args, "region_dl_device", "cpu")).strip().lower()
    cuda_ok = False
    if dev_arg in {"auto", "cuda"}:
        try:
            cuda_ok = bool(torch.cuda.is_available())
        except Exception:
            cuda_ok = False
    if dev_arg == "auto":
        dev_arg = "cuda" if cuda_ok else "cpu"
    elif dev_arg == "cuda" and not cuda_ok:
        dev_arg = "cpu"
    if dev_arg not in {"cpu", "cuda"}:
        dev_arg = "cpu"
    device = torch.device(dev_arg)
    region_seed = int(getattr(args, "region_dl_seed", 20260620))
    torch.manual_seed(region_seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(region_seed)
    model = TinyUNet(X.shape[0], base=int(getattr(args, "region_dl_base_channels", 16)), n_cls=3).to(device)
    x_t = torch.from_numpy(X[None, ...]).to(device)
    y_t = torch.from_numpy(y[None, ...]).to(device)

    # Robust class weights plus explicit ARC/ME weights. Background is intentionally low.
    counts = np.bincount(y.ravel(), minlength=3).astype(np.float32)
    inv = 1.0 / np.maximum(counts, 1.0)
    weights = inv / max(float(inv.mean()), 1e-6)
    weights[0] = min(weights[0], float(getattr(args, "region_dl_background_class_weight", 0.20)))
    weights[REGION_NAME_TO_CODE["ARC"]] *= float(getattr(args, "region_dl_arc_class_weight", 1.0))
    weights[REGION_NAME_TO_CODE["ME"]] *= float(getattr(args, "region_dl_me_class_weight", 2.5))
    weights = weights / max(float(np.mean(weights[weights > 0])), 1e-6)
    w_t = torch.tensor(weights, dtype=torch.float32, device=device)
    opt = torch.optim.AdamW(model.parameters(), lr=float(getattr(args, "region_dl_lr", 1e-3)), weight_decay=1e-4)
    model.train()
    final_loss = float("nan")
    log_every = int(max(1, getattr(args, "region_dl_log_every", 1)))
    _region_debug(args,
        "DL start: device=%s epochs=%s size=%sx%s prior_ARC=%s prior_ME=%s weights=%s floor_present=%s method=%s",
        str(device), int(epochs), int(h2), int(w2), int(np.count_nonzero(prior == REGION_NAME_TO_CODE["ARC"])),
        int(np.count_nonzero(prior == REGION_NAME_TO_CODE["ME"])), [round(float(x), 3) for x in weights.tolist()],
        bool(info.get("me_floor_present", False)), str(info.get("me_floor_method", "NA"))
    )
    try:
        for ep in range(1, epochs + 1):
            opt.zero_grad(set_to_none=True)
            logits = model(x_t)
            loss = F.cross_entropy(logits, y_t, weight=w_t)
            loss.backward()
            opt.step()
            final_loss = float(loss.detach().cpu())
            if bool(getattr(args, "region_dl_verbose", True)) and (ep == 1 or ep == epochs or ep % log_every == 0):
                _region_debug(args, "DL epoch %03d/%03d loss=%.6f", int(ep), int(epochs), float(final_loss))
    except Exception as e:
        info["region_ai_model"] = f"projected_ventricle_walls_prior_unet_training_failed_fallback: {e}"
        return prior, info

    model.eval()
    with torch.no_grad():
        logits = model(x_t)
        prob = F.softmax(logits, dim=1)[0].detach().cpu().numpy()
    pred_small = prob.argmax(axis=0).astype(np.uint8)
    conf_small = prob.max(axis=0).astype(np.float32)
    pred = _resize_float(pred_small, (H, W), order=0).round().astype(np.uint8)
    conf = _resize_float(conf_small, (H, W), order=1).astype(np.float32)
    min_conf = float(getattr(args, "region_dl_min_confidence", 0.45))
    pred = np.where(conf >= min_conf, pred, prior).astype(np.uint8)
    pred = _postprocess_arc_me_prediction(pred, prior, info, args)

    info.update({
        "region_dl_epochs": int(epochs),
        "region_dl_size": int(max(h2, w2)),
        "region_dl_device": str(device),
        "region_dl_seed": int(region_seed),
        "region_dl_final_loss": float(final_loss) if np.isfinite(final_loss) else np.nan,
        "region_dl_min_confidence": float(min_conf),
        "region_dl_boundary_refine_um": float(max(0.0, getattr(args, "region_dl_boundary_refine_um", 0.0))),
        "region_dl_boundary_refine_px": int(um_to_px(float(max(0.0, getattr(args, "region_dl_boundary_refine_um", 0.0))), args, min_px=0)),
        "arc_area_px": int(np.count_nonzero(pred == REGION_NAME_TO_CODE["ARC"])),
        "me_area_px": int(np.count_nonzero(pred == REGION_NAME_TO_CODE["ME"])),
    })
    return pred.astype(np.uint8), info


def fallback_bin_region_for_centroid(y: float, H: int, args: argparse.Namespace, region_map: Dict[int, str]) -> Tuple[int, str]:
    n_bins = int(max(1, args.n_bins))
    b = bin_index_from_y(y, H, n_bins)
    r = str(region_map.get(b, "ARC"))
    if bool(getattr(args, "arc_me_only", True)) and r not in {"ARC", "ME"}:
        r = "OUTSIDE_ARC_ME"
    return int(b), r


def region_name_from_mask(region_mask: Any, y: float, x: float) -> str:
    if region_mask is None:
        return "OUTSIDE_ARC_ME"
    H, W = np.asarray(region_mask).shape[:2]
    yy = int(max(0, min(H - 1, round(float(y)))))
    xx = int(max(0, min(W - 1, round(float(x)))))
    return REGION_CODE_TO_NAME.get(int(region_mask[yy, xx]), "OUTSIDE_ARC_ME")


def region_name_from_projected_walls(region_info: Dict[str, Any], y: float, x: float) -> str:
    """Return ME only for the closed ventricular floor; sidewalls remain ARC."""
    if not bool((region_info or {}).get("me_floor_present", False)):
        return "ARC"
    ys = np.asarray((region_info or {}).get("_wall_y", []), dtype=float)
    lx = np.asarray((region_info or {}).get("_wall_left_x", []), dtype=float)
    rx = np.asarray((region_info or {}).get("_wall_right_x", []), dtype=float)
    good = np.isfinite(ys) & np.isfinite(lx) & np.isfinite(rx)
    if int(np.count_nonzero(good)) < 2:
        return "ARC"
    ys = ys[good]; lx = lx[good]; rx = rx[good]
    order = np.argsort(ys); ys = ys[order]; lx = lx[order]; rx = rx[order]
    yy = float(y); xx = float(x)
    if yy < float(ys[0]) or yy > float(ys[-1]):
        return "ARC"
    left_here = float(np.interp(yy, ys, lx))
    right_here = float(np.interp(yy, ys, rx))
    if right_here < left_here:
        left_here, right_here = right_here, left_here
    if xx < left_here or xx > right_here:
        return "ARC"

    cx = np.asarray((region_info or {}).get("_me_floor_curve_x", []), dtype=float)
    cy = np.asarray((region_info or {}).get("_me_floor_curve_y", []), dtype=float)
    goodc = np.isfinite(cx) & np.isfinite(cy)
    if int(np.count_nonzero(goodc)) < 2:
        return "ARC"
    cx = cx[goodc]; cy = cy[goodc]
    order = np.argsort(cx); cx = cx[order]; cy = cy[order]
    if xx < float(cx[0]) or xx > float(cx[-1]):
        return "ARC"
    floor_here = float(np.interp(xx, cx, cy))
    margin = float((region_info or {}).get("me_floor_centroid_margin_px", 0.0))
    if yy < floor_here + margin:
        return "ARC"
    side_dist = min(xx - left_here, right_here - xx)
    side_excl = float((region_info or {}).get("me_sidewall_exclusion_px", 0.0))
    release = float((region_info or {}).get("me_sidewall_release_depth_px", 0.0))
    if side_dist < side_excl and yy < floor_here + release:
        return "ARC"
    return "ME"


def region_order_from_args_arc_me(args: argparse.Namespace) -> List[str]:
    requested = [x.strip() for x in str(getattr(args, "region_order", "ARC,ME,VMN") or "ARC,ME,VMN").split(",") if x.strip()]
    out = []
    for r in requested:
        rr = r.upper() if r.lower() in {"arc", "me", "vmn"} else r
        if bool(getattr(args, "arc_me_only", True)) and rr not in {"ARC", "ME", "VMN"}:
            continue
        if rr not in out:
            out.append(rr)
    for r in ["ARC", "ME", "VMN"]:
        if r not in out:
            out.append(r)
    if not bool(getattr(args, "arc_me_only", True)) and "OUTSIDE_ARC_ME" not in out:
        out.append("OUTSIDE_ARC_ME")
    return out


def save_region_mask_tif(pair: SectionPair, region_mask: Any, outdir: Path, args: argparse.Namespace) -> None:
    if bool(args.dry_run) or not bool(getattr(args, "save_region_mask_tifs", True)):
        return
    try:
        base = outdir / args.reconstructed_dirname / pair.animal_id / f"{pair.sample}_S{pair.section_index:02d}{('_T'+pair.tile) if pair.tile else ''}"
        ensure_dir(base)
        tifffile.imwrite(str(base / "reconstructed_ARC_ME_region_mask.tif"), np.asarray(region_mask).astype(np.uint8))
    except Exception:
        pass

def save_wall_projection_csv(pair: SectionPair, region_info: Dict[str, Any], outdir: Path, args: argparse.Namespace) -> None:
    """Save projected walls, observed oval walls, and the curved floor envelope."""
    if bool(args.dry_run) or not bool(getattr(args, "save_wall_projection_csv", True)):
        return
    ys = np.asarray((region_info or {}).get("_wall_y", []), dtype=float)
    left = np.asarray((region_info or {}).get("_wall_left_x", []), dtype=float)
    right = np.asarray((region_info or {}).get("_wall_right_x", []), dtype=float)
    good = np.isfinite(ys) & np.isfinite(left) & np.isfinite(right)
    base = outdir / args.reconstructed_dirname / pair.animal_id / f"{pair.sample}_S{pair.section_index:02d}{('_T'+pair.tile) if pair.tile else ''}"
    ensure_dir(base)
    if int(np.count_nonzero(good)) > 0:
        ys2 = ys[good]; left2 = left[good]; right2 = right[good]
        presence_y = float((region_info or {}).get("wall_presence_y_start", (region_info or {}).get("me_wall_presence_y_start", np.nan)))
        tip_y = float((region_info or {}).get("wall_ventral_tip_y", (region_info or {}).get("me_wall_ventral_tip_y", np.nan)))
        df = pd.DataFrame({
            "y_px": ys2,
            "left_wall_x_px": left2,
            "right_wall_x_px": right2,
            "corridor_width_px": right2 - left2,
            "corridor_center_x_px": (right2 + left2) / 2.0,
            "at_or_below_ventral_tip": ys2 >= tip_y if np.isfinite(tip_y) else False,
            "used_for_me_presence_evidence": ys2 >= presence_y if np.isfinite(presence_y) else False,
            "oval_detected": bool((region_info or {}).get("wall_oval_detected", False)),
            "oval_tangent_start_y": (region_info or {}).get("wall_oval_peak_y", np.nan),
            "me_present": bool((region_info or {}).get("me_floor_present", False)),
        })
        df.to_csv(base / "reconstructed_projected_ventricle_walls.csv", index=False)

    oy = np.asarray((region_info or {}).get("_wall_observed_y", []), dtype=float)
    ol = np.asarray((region_info or {}).get("_wall_observed_left_x", []), dtype=float)
    orr = np.asarray((region_info or {}).get("_wall_observed_right_x", []), dtype=float)
    goodo = np.isfinite(oy) & np.isfinite(ol) & np.isfinite(orr)
    if int(np.count_nonzero(goodo)) > 0:
        pd.DataFrame({
            "y_px": oy[goodo],
            "observed_left_wall_x_px": ol[goodo],
            "observed_right_wall_x_px": orr[goodo],
            "observed_lumen_width_px": orr[goodo] - ol[goodo],
            "oval_detected": bool((region_info or {}).get("wall_oval_detected", False)),
        }).to_csv(base / "reconstructed_observed_ventricle_walls.csv", index=False)

    fx = np.asarray((region_info or {}).get("_me_floor_curve_x", []), dtype=float)
    fy = np.asarray((region_info or {}).get("_me_floor_curve_y", []), dtype=float)
    fyr = np.asarray((region_info or {}).get("_me_floor_curve_y_raw", []), dtype=float)
    goodf = np.isfinite(fx) & np.isfinite(fy)
    if int(np.count_nonzero(goodf)) > 0:
        raw = fyr[goodf] if fyr.size == fx.size else np.full(int(np.count_nonzero(goodf)), np.nan)
        pd.DataFrame({
            "x_px": fx[goodf],
            "floor_y_px": fy[goodf],
            "raw_lumen_bottom_y_px": raw,
            "me_present": bool((region_info or {}).get("me_floor_present", False)),
            "lumen_open_to_bottom": bool((region_info or {}).get("me_floor_open_bottom", False)),
            "sidewall_exclusion_px": (region_info or {}).get("me_sidewall_exclusion_px", 0),
        }).to_csv(base / "reconstructed_ventricle_floor_curve.csv", index=False)


def render_reconstructed_qc_pdf(
    out_path: Path,
    *,
    dapi_lbl: Any,
    section_row: Dict[str, Any],
    cell_df: Any,
    region_map: Dict[int, str],
    args: argparse.Namespace,
    region_mask: Optional[Any] = None,
    region_info: Optional[Dict[str, Any]] = None,
) -> None:
    """Save full-section QC with the actual projected ventricular-wall separators.

    v17 does not draw a fixed horizontal ARC/ME boundary. The white/navy lines are
    the fitted third-ventricle walls projected ventrally. ME is only the closed floor
    below the curved lumen-bottom envelope; ventricular sidewall nuclei remain ARC.
    """
    if bool(args.dry_run):
        return
    H, W = dapi_lbl.shape
    cls = np.zeros(int(dapi_lbl.max()) + 1, dtype=np.uint8)
    if cell_df is not None and not cell_df.empty:
        for r in cell_df.itertuples(index=False):
            lab = int(getattr(r, "nucleus_label"))
            if lab <= 0 or lab >= cls.size:
                continue
            is_cfos = bool(getattr(r, "is_cfos"))
            is_npy = bool(getattr(r, "is_npy"))
            is_pomc = bool(getattr(r, "is_pomc"))
            if is_cfos and is_npy and is_pomc:
                cls[lab] = 6
            elif is_cfos and is_npy:
                cls[lab] = 4
            elif is_cfos and is_pomc:
                cls[lab] = 5
            elif is_cfos:
                cls[lab] = 3
            elif is_npy:
                cls[lab] = 1
            elif is_pomc:
                cls[lab] = 2

    class_img = cls[dapi_lbl]
    rgb = np.zeros((H, W, 3), dtype=np.float32)
    rgb[dapi_lbl > 0] = (0.82, 0.82, 0.82)
    colors = {
        1: (0.00, 0.75, 0.00),
        2: (1.00, 0.10, 0.05),
        3: (0.95, 0.00, 0.95),
        4: (1.00, 1.00, 0.00),
        5: (1.00, 0.55, 0.95),
        6: (0.00, 1.00, 1.00),
    }
    for k, col in colors.items():
        rgb[class_img == k] = col

    img_h = 10.5
    img_w = img_h * (W / max(1, H))
    panel_w = 6.8
    fig = plt.figure(figsize=(img_w + panel_w, img_h), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[img_w, panel_w])
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(rgb, interpolation="nearest")
    ax.set_axis_off()
    region_info = region_info or {}
    me_present = bool(region_info.get("me_floor_present", False))

    # Draw only the ME tissue contour (when present). ARC is everything else and does
    # not need a misleading large outline.
    if region_mask is not None and me_present:
        rm = np.asarray(region_mask)
        try:
            if np.any(rm == REGION_NAME_TO_CODE["ME"]):
                ax.contour(
                    (rm == REGION_NAME_TO_CODE["ME"]).astype(float),
                    levels=[0.5], colors=["orange"], linewidths=2.3,
                )
        except Exception:
            pass

    # For an oval lumen, show the observed upper oval walls as a thin trace, then the
    # straight projected tangents as the authoritative separators below the lateral bulge.
    if bool(region_info.get("wall_oval_detected", False)):
        oy = np.asarray(region_info.get("_wall_observed_y", []), dtype=float)
        ol = np.asarray(region_info.get("_wall_observed_left_x", []), dtype=float)
        orr = np.asarray(region_info.get("_wall_observed_right_x", []), dtype=float)
        tangent_y = float(region_info.get("wall_oval_peak_y", np.nan))
        goodo = np.isfinite(oy) & np.isfinite(ol) & np.isfinite(orr)
        if np.isfinite(tangent_y):
            goodo &= oy <= tangent_y
        if int(np.count_nonzero(goodo)) >= 2:
            oo = np.argsort(oy[goodo])
            for wx in (ol[goodo][oo], orr[goodo][oo]):
                ax.plot(wx, oy[goodo][oo], color="navy", lw=3.8, solid_capstyle="round", zorder=6)
                ax.plot(wx, oy[goodo][oo], color="white", lw=1.5, solid_capstyle="round", zorder=7)

    # Authoritative projected tangents. A dark navy underlay keeps the white lines visible.
    wall_y = np.asarray(region_info.get("_wall_y", []), dtype=float)
    wall_l = np.asarray(region_info.get("_wall_left_x", []), dtype=float)
    wall_r = np.asarray(region_info.get("_wall_right_x", []), dtype=float)
    good = np.isfinite(wall_y) & np.isfinite(wall_l) & np.isfinite(wall_r)
    if int(np.count_nonzero(good)) >= 2:
        wy = wall_y[good]; wl = wall_l[good]; wr = wall_r[good]
        order = np.argsort(wy); wy = wy[order]; wl = wl[order]; wr = wr[order]
        for wx in (wl, wr):
            ax.plot(wx, wy, color="navy", lw=6.2, solid_capstyle="round", zorder=7)
            ax.plot(wx, wy, color="white", lw=3.0, solid_capstyle="round", zorder=8)

    # Curved ventricular floor used to separate ME floor from ARC sidewalls.
    fx = np.asarray(region_info.get("_me_floor_curve_x", []), dtype=float)
    fy = np.asarray(region_info.get("_me_floor_curve_y", []), dtype=float)
    goodf = np.isfinite(fx) & np.isfinite(fy)
    if int(np.count_nonzero(goodf)) >= 2:
        order = np.argsort(fx[goodf])
        ax.plot(fx[goodf][order], fy[goodf][order], color="orange", lw=2.2, ls="--", zorder=9)

    # Region-specific counts in the image corner.
    if cell_df is not None and not cell_df.empty and "region" in cell_df.columns:
        arc_df = cell_df[cell_df["region"] == "ARC"]
        me_df = cell_df[cell_df["region"] == "ME"]
        arc_cfos = int(arc_df["is_cfos"].sum()) if not arc_df.empty else 0
        arc_dapi = int(arc_df.shape[0])
        ax.text(
            12, int(H * 0.07), f"ARC: c-FOS/DAPI {arc_cfos}/{arc_dapi}",
            color="white", ha="left", va="top", fontsize=12, weight="bold",
            path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="navy")],
        )
        if me_present:
            me_cfos = int(me_df["is_cfos"].sum()) if not me_df.empty else 0
            me_dapi = int(me_df.shape[0])
            ax.text(
                12, int(H * 0.14), f"ME: c-FOS/DAPI {me_cfos}/{me_dapi}",
                color="orange", ha="left", va="top", fontsize=12, weight="bold",
                path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")],
            )
        else:
            ax.text(
                12, int(H * 0.14), "ME: absent (no closed ventricular floor)",
                color="orange", ha="left", va="top", fontsize=11.5, weight="bold",
                path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")],
            )

    if int(args.stitch_gap_px) > 0:
        ax.axvline(W / 2.0, color="cyan", lw=0.8, alpha=0.45)
    if args.um_per_px and args.um_per_px > 0:
        bar_px = int(round(float(args.scalebar_um) / float(args.um_per_px)))
        if 0 < bar_px < W:
            x0 = W - bar_px - 50
            y0 = H - 55
            ax.add_patch(plt.Rectangle((x0, y0), bar_px, 10, color="white", ec="black", lw=1.3))
            ax.text(
                x0 + bar_px / 2, y0 - 8, f"{int(args.scalebar_um)} µm", color="white",
                ha="center", va="bottom", fontsize=13, weight="bold",
                path_effects=[__import__("matplotlib.patheffects", fromlist=["withStroke"]).withStroke(linewidth=3, foreground="black")],
            )

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_axis_off()
    if cell_df is not None and not cell_df.empty and "region" in cell_df.columns:
        arc_df = cell_df[cell_df["region"] == "ARC"]
        me_df = cell_df[cell_df["region"] == "ME"]
        arc_line = f"ARC c-FOS/DAPI: {int(arc_df['is_cfos'].sum())}/{int(arc_df.shape[0])}"
        me_line = (
            f"ME c-FOS/DAPI: {int(me_df['is_cfos'].sum())}/{int(me_df.shape[0])}"
            if me_present else
            "ME: absent (lumen open or no closed floor tissue)"
        )
    else:
        arc_line = "ARC c-FOS/DAPI: NA"
        me_line = "ME c-FOS/DAPI: NA"

    def _fmt(v: Any, digits: int = 2) -> str:
        try:
            f = float(v)
            return f"{f:.{digits}f}" if np.isfinite(f) else "NA"
        except Exception:
            return "NA"

    lines = [
        f"Animal: {section_row.get('animal_id', '')}",
        f"Source folder: {section_row.get('sample', '')}",
        f"Section: S{int(section_row.get('section_index', 0)):02d}{' T'+str(section_row.get('tile')) if section_row.get('tile') else ''}",
        f"Repetition folder: {section_row.get('is_repetition', False)} {section_row.get('repetition_label', '')}",
        "",
        f"DAPI nuclei, whole image: {int(section_row.get('dapi_nuclei', 0))}",
        f"c-FOS+ nuclei, whole image: {int(section_row.get('cfos_cells', 0))}",
        f"POMC+ nuclei, whole image: {int(section_row.get('pomc_cells', 0))}",
        f"NPY+ nuclei, whole image: {int(section_row.get('npy_cells', 0))}",
        "",
        arc_line,
        me_line,
        "",
        f"c-FOS only: {int(section_row.get('cfos_only', 0))}",
        f"c-FOS∧POMC: {int(section_row.get('cfos_pomc', 0))}",
        f"c-FOS∧NPY: {int(section_row.get('cfos_npy', 0))}",
        f"Triple: {int(section_row.get('cfos_npy_pomc', 0))}",
        "",
        f"Region model: {region_info.get('region_ai_model', args.region_mode)}",
        "Rule: below closed curved ventricle floor + between walls = ME",
        "Rule: ventricular sidewall nuclei + all remaining nuclei = ARC",
        "Fixed horizontal ARC/ME boundary: none",
        f"Midline x: {_fmt(region_info.get('mid_x'), 1)} px",
        f"Ventricle row coverage: {_fmt(region_info.get('ventricle_coverage'), 2)}",
        f"Wall fit status: {region_info.get('wall_fit_status', 'NA')}",
        f"Oval lumen detected: {region_info.get('wall_oval_detected', False)}; tangent y: {region_info.get('wall_oval_peak_y', 'NA')}",
        f"Oval bulge ratios top/bottom: {_fmt(region_info.get('wall_oval_bulge_ratio_top'), 2)}/{_fmt(region_info.get('wall_oval_bulge_ratio_bottom'), 2)}",
        f"Oval raw slopes L/R: {_fmt(region_info.get('wall_oval_left_slope_raw'), 4)}/{_fmt(region_info.get('wall_oval_right_slope_raw'), 4)}; opposite corrected: {region_info.get('wall_oval_opposite_tangent_corrected', False)}",
        f"Oval corridor width start/end: {_fmt(region_info.get('wall_oval_projection_width_start_px'), 1)}/{_fmt(region_info.get('wall_oval_projection_width_end_px'), 1)} px; monotonic: {region_info.get('wall_oval_monotonic_divergence_applied', False)}",
        f"Supported wall run y: {region_info.get('wall_main_run_start_y', 'NA')}–{region_info.get('wall_main_run_end_y', 'NA')}",
        f"Fitted wall rows y: {region_info.get('wall_fit_row_start_y', 'NA')}–{region_info.get('wall_fit_row_end_y', 'NA')} (n={region_info.get('wall_fit_rows_n', 'NA')})",
        f"Left wall slope/RMSE: {_fmt(region_info.get('wall_left_slope'), 4)} / {_fmt(region_info.get('wall_left_rmse'), 2)} px",
        f"Right wall slope/RMSE: {_fmt(region_info.get('wall_right_slope'), 4)} / {_fmt(region_info.get('wall_right_rmse'), 2)} px",
        f"Wall lines y: {region_info.get('wall_line_y_start', 'NA')}–{region_info.get('wall_projection_y_end', 'NA')}",
        f"Ventral lumen tip y: {region_info.get('wall_ventral_tip_y', 'NA')}",
        f"ME detected: {me_present}; floor status: {region_info.get('me_floor_geometry_status', 'NA')}",
        f"Lumen open to bottom: {region_info.get('me_floor_open_bottom', 'NA')}; sparse rescue: {region_info.get('me_floor_sparse_rescue', False)}",
        f"Sparse floor evidence n L/R/mid: {region_info.get('me_floor_sparse_candidate_nuclei', 0)} {region_info.get('me_floor_sparse_left_nuclei', 0)}/{region_info.get('me_floor_sparse_right_nuclei', 0)}/{region_info.get('me_floor_sparse_midline_nuclei', 0)}",
        f"Sparse floor span/depth/fraction/ratio: {_fmt(region_info.get('me_floor_sparse_span_px'), 1)}/{_fmt(region_info.get('me_floor_sparse_depth_px'), 1)}/{_fmt(region_info.get('me_floor_sparse_fraction'), 2)}/{_fmt(region_info.get('me_floor_sparse_density_ratio'), 2)}",
        f"Corridor/floor nuclei: {region_info.get('me_wall_enclosed_nuclei', 'NA')}/{region_info.get('me_wall_floor_nuclei', 'NA')}",
        f"Sidewall nuclei reassigned to ARC: {region_info.get('me_wall_sidewall_reassigned_nuclei', 'NA')}",
        f"Floor midline nuclei: {region_info.get('me_wall_midline_nuclei', 'NA')}",
        f"Enclosed depth/span: {_fmt(region_info.get('me_wall_depth_px'), 1)}/{_fmt(region_info.get('me_wall_span_px'), 1)} px",
        f"Regions quantified: ARC and ME; outside action: {getattr(args, 'outside_arc_me_action', 'arc')}",
        "",
        f"c-FOS association: {args.cfos_association_mode}",
        f"Strict c-FOS/DAPI overlap: {float(args.cfos_dapi_overlap):.2f}",
        f"Blur expansion: {float(args.cfos_blur_um):.2f} µm",
    ]
    ax2.text(0.02, 0.985, "\n".join(lines), ha="left", va="top", fontsize=10.4, weight="bold")

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], color="navy", lw=6, label="Projected ventricle wall"),
        Line2D([0], [0], color="white", lw=2.5, label="Wall centerline"),
    ]
    if bool(region_info.get("wall_oval_detected", False)):
        handles.append(Line2D([0], [0], color="navy", lw=2.5, label="Observed oval wall"))
    if np.asarray(region_info.get("_me_floor_curve_x", [])).size >= 2:
        handles.append(Line2D([0], [0], color="orange", lw=2, ls="--", label="Curved ventricle floor"))
    if me_present:
        handles.append(Patch(facecolor="none", edgecolor="orange", linewidth=2, label="ME floor tissue"))
    handles += [
        Patch(facecolor=colors[2], label="POMC"),
        Patch(facecolor=colors[1], label="NPY"),
        Patch(facecolor=colors[3], label="c-FOS only"),
        Patch(facecolor=colors[5], label="c-FOS∧POMC"),
        Patch(facecolor=colors[4], label="c-FOS∧NPY"),
        Patch(facecolor=colors[6], label="Triple"),
    ]
    ax2.legend(handles=handles, loc="lower left", fontsize=9.4, framealpha=0.95)
    ensure_dir(out_path.parent)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

def process_reconstructed_section(pair: SectionPair, args: argparse.Namespace, outdir: Path, region_map: Dict[int, str], logger: logging.Logger) -> Tuple[List[Dict[str, Any]], Dict[str, Any], List[Dict[str, Any]]]:
    assignment_rows: List[Dict[str, Any]] = []
    cell_rows: List[Dict[str, Any]] = []

    def _seg_available(plane: Optional[PlaneFiles], marker: str) -> bool:
        if plane is None:
            return False
        p = plane.segs.get(marker)
        return bool(p is not None and Path(p).exists())

    pomc_left_available = _seg_available(pair.left, "POMC")
    pomc_right_available = _seg_available(pair.right, "POMC")
    npy_left_available = _seg_available(pair.left, "NPY")
    npy_right_available = _seg_available(pair.right, "NPY")
    existing_halves = int(pair.left is not None) + int(pair.right is not None)
    pomc_halves_available = int(pomc_left_available) + int(pomc_right_available)
    npy_halves_available = int(npy_left_available) + int(npy_right_available)
    has_pomc_channel = pomc_halves_available > 0
    has_npy_channel = npy_halves_available > 0
    pomc_channel_complete = bool(existing_halves > 0 and pomc_halves_available == existing_halves)
    npy_channel_complete = bool(existing_halves > 0 and npy_halves_available == existing_halves)

    dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl, missing = load_and_stitch_section_labels(pair, args)

    base_section_row: Dict[str, Any] = {
        "sample": pair.sample,
        "sample_key": pair.sample_key,
        "animal_id": pair.animal_id,
        "animal_key": pair.animal_key,
        "is_repetition": bool(pair.is_repetition),
        "repetition_label": pair.repetition_label,
        "section_index": int(pair.section_index),
        "tile": pair.tile,
        "has_left": pair.left is not None,
        "has_right": pair.right is not None,
        "left_dir": str(pair.left.plane_dir) if pair.left else "",
        "right_dir": str(pair.right.plane_dir) if pair.right else "",
        "missing": ";".join(missing),
        "has_pomc_channel": bool(has_pomc_channel),
        "pomc_channel_complete": bool(pomc_channel_complete),
        "pomc_halves_available": int(pomc_halves_available),
        "has_npy_channel": bool(has_npy_channel),
        "npy_channel_complete": bool(npy_channel_complete),
        "npy_halves_available": int(npy_halves_available),
        "status": "ok",
    }

    # DAPI and c-FOS are required for quantification, but a section is still usable
    # when one half is missing and the other half is valid. Missing sides are filled
    # with zeros during stitching. Therefore we fail only if the reconstructed masks
    # themselves are absent or incompatible.
    if dapi_lbl is None or cfos_lbl is None:
        base_section_row.update({"status": "missing_required", "dapi_nuclei": 0, "cfos_cells": 0, "npy_cells": 0, "pomc_cells": 0})
        return cell_rows, base_section_row, assignment_rows
    if not validate_same_shape([dapi_lbl, cfos_lbl, npy_lbl, pomc_lbl]):
        base_section_row.update({"status": "shape_mismatch", "dapi_nuclei": 0, "cfos_cells": 0, "npy_cells": 0, "pomc_cells": 0})
        return cell_rows, base_section_row, assignment_rows

    dapi_lbl = ensure_label_image(dapi_lbl)
    cfos_lbl = filter_rois_by_area(ensure_label_image(cfos_lbl), min_area=int(args.cfos_min_area_px), max_area=int(args.cfos_max_area_px or 0))
    npy_lbl = ensure_label_image(npy_lbl)
    pomc_lbl = ensure_label_image(pomc_lbl)

    pomc_lbl_clean, pomc_removed = filter_pomc_rois_by_npy_overlap(pomc_lbl, npy_lbl > 0, float(args.pct_pomc_npy_overlap))
    cfos_nuclei, fos_assignments = assign_cfos_rois_bluraware(dapi_lbl, cfos_lbl, args)
    for fa in fos_assignments:
        row = asdict(fa)
        row.update({
            "sample": pair.sample, "animal_id": pair.animal_id, "animal_key": pair.animal_key,
            "section_index": pair.section_index, "tile": pair.tile,
            "is_repetition": pair.is_repetition, "repetition_label": pair.repetition_label,
        })
        assignment_rows.append(row)

    um_per_px = float(args.um_per_px)
    pomc_expand_px = int(round(float(args.pomc_expand_um) / um_per_px)) if um_per_px > 0 else int(round(args.pomc_expand_um))
    npy_expand_px = int(round(float(args.npy_expand_um) / um_per_px)) if um_per_px > 0 else int(round(args.npy_expand_um))
    pomc_nuclei_raw = nuclei_positive_by_expanded_label(dapi_lbl, pomc_lbl, pomc_expand_px, int(args.marker_min_pixels))
    pomc_nuclei = nuclei_positive_by_expanded_label(dapi_lbl, pomc_lbl_clean, pomc_expand_px, int(args.marker_min_pixels))
    npy_nuclei = nuclei_positive_by_expanded_label(dapi_lbl, npy_lbl, npy_expand_px, int(args.marker_min_pixels))
    pomc_nuclei, n_pomc_dropped = apply_strict_pomc_vs_npy(dapi_lbl, npy_nuclei, pomc_nuclei, npy_lbl, pomc_lbl_clean, float(args.pomc_in_npy))

    H, W = dapi_lbl.shape
    n_bins = int(max(1, args.n_bins))
    region_mask = None
    region_info: Dict[str, Any] = {"region_mode": str(args.region_mode)}
    if str(args.region_mode) == "arc-me-dl":
        region_mask, region_info = infer_arc_me_region_mask_deeplearning(dapi_lbl, npy_lbl, pomc_lbl_clean, args, cfos_lbl=cfos_lbl)
    elif str(args.region_mode) in {"arc-me-walls", "arc-me-cv", "arc-me-ventricle"}:
        # arc-me-cv / arc-me-ventricle are retained as aliases for compatibility.
        region_mask, region_info = infer_arc_me_region_mask(dapi_lbl, npy_lbl, pomc_lbl_clean, args, cfos_lbl=cfos_lbl)
    if region_mask is not None:
        try:
            save_wall_projection_csv(pair, region_info, outdir, args)
        except Exception:
            logger.exception("Failed saving ventricle-wall coordinates for %s S%02d", pair.sample, pair.section_index)
    for r in regionprops(dapi_lbl):
        lab = int(r.label)
        y, x = float(r.centroid[0]), float(r.centroid[1])
        b, fallback_region = fallback_bin_region_for_centroid(y, H, args, region_map)
        if region_mask is not None and str(args.region_mode) in {"arc-me-walls", "arc-me-dl", "arc-me-cv", "arc-me-ventricle"}:
            # Authoritative v17 rule: ME requires a closed curved ventricular floor,
            # centroid below that local floor, and separation from the sidewall band.
            region = region_name_from_projected_walls(region_info, y, x)
            if str(getattr(args, "outside_arc_me_action", "arc")) == "exclude":
                raw_region = region_name_from_mask(region_mask, y, x)
                if raw_region not in {"ARC", "ME"}:
                    region = "OUTSIDE_ARC_ME"
        elif region_mask is not None:
            raw_region = region_name_from_mask(region_mask, y, x)
            if raw_region in {"ARC", "ME"}:
                region = raw_region
            elif str(getattr(args, "outside_arc_me_action", "arc")) == "exclude":
                region = "OUTSIDE_ARC_ME"
            else:
                region = "ARC"
        else:
            region = fallback_region if fallback_region in {"ARC", "ME"} else "ARC"
        if bool(getattr(args, "arc_me_only", True)) and region not in {"ARC", "ME"}:
            if str(getattr(args, "outside_arc_me_action", "arc")) == "exclude":
                region = "OUTSIDE_ARC_ME"
            else:
                region = "ARC"
        is_cfos = lab in cfos_nuclei
        is_npy = lab in npy_nuclei
        is_pomc = lab in pomc_nuclei
        cell_rows.append({
            "sample": pair.sample,
            "sample_key": pair.sample_key,
            "animal_id": pair.animal_id,
            "animal_key": pair.animal_key,
            "is_repetition": bool(pair.is_repetition),
            "repetition_label": pair.repetition_label,
            "section_index": int(pair.section_index),
            "tile": pair.tile,
            "nucleus_label": lab,
            "centroid_y": y,
            "centroid_x": x,
            "bin": int(b),
            "region": region,
            "area_px": int(r.area),
            "has_pomc_channel": bool(has_pomc_channel),
            "pomc_channel_complete": bool(pomc_channel_complete),
            "has_npy_channel": bool(has_npy_channel),
            "npy_channel_complete": bool(npy_channel_complete),
            "is_cfos": bool(is_cfos),
            "is_npy": bool(is_npy),
            "is_pomc": bool(is_pomc),
            "phenotype": classify_phenotype(is_cfos, is_npy, is_pomc),
        })

    section_df = pd.DataFrame(cell_rows)
    if section_df.empty:
        section_counts = {"dapi_nuclei": 0, "cfos_cells": 0, "npy_cells": 0, "pomc_cells": 0, "cfos_only": 0, "cfos_npy": 0, "cfos_pomc": 0, "cfos_npy_pomc": 0}
    else:
        section_summary = summarize_cell_rows(section_df, ["sample", "animal_id", "animal_key", "section_index", "tile"])
        section_counts = section_summary.iloc[0].to_dict() if not section_summary.empty else {}

    base_section_row.update(section_counts)
    base_section_row.update({
        "dapi_nuclei": int(count_labels(dapi_lbl)),
        "cfos_rois_total": int(count_labels(cfos_lbl)),
        "cfos_rois_accepted": int(sum(1 for x in fos_assignments if x.accepted)),
        "pomc_cells_raw": int(len(pomc_nuclei_raw)),
        "pomc_rois_removed_by_npy": int(pomc_removed),
        "pomc_nuclei_dropped_strict_npy": int(n_pomc_dropped),
        "image_height_px": int(H),
        "image_width_px": int(W),
        "n_bins": int(n_bins),
    })
    # Add reproducible ARC/ME CV geometry columns to the section summary.
    for _k, _v in (region_info or {}).items():
        if isinstance(_v, (str, int, float, bool)) or _v is None:
            base_section_row[str(_k)] = _v

    pixel_area_um2 = float(args.um_per_px) ** 2
    base_section_row["pixel_area_um2"] = pixel_area_um2
    base_section_row["arc_area_um2"] = float(base_section_row.get("arc_area_px", 0) or 0) * pixel_area_um2
    base_section_row["me_area_um2"] = float(base_section_row.get("me_area_px", 0) or 0) * pixel_area_um2
    base_section_row["analyzed_area_px"] = int(base_section_row.get("arc_area_px", 0) or 0) + int(base_section_row.get("me_area_px", 0) or 0)
    base_section_row["analyzed_area_um2"] = float(base_section_row["analyzed_area_px"]) * pixel_area_um2
    base_section_row["pomc_assay_area_px"] = int(base_section_row["analyzed_area_px"]) if has_pomc_channel else 0
    base_section_row["pomc_assay_area_um2"] = float(base_section_row["pomc_assay_area_px"]) * pixel_area_um2
    base_section_row["npy_assay_area_px"] = int(base_section_row["analyzed_area_px"]) if has_npy_channel else 0
    base_section_row["npy_assay_area_um2"] = float(base_section_row["npy_assay_area_px"]) * pixel_area_um2

    if bool(args.save_reconstructed_tifs):
        try:
            # Quantification is mask-based. Raw channel TIFF export is only for QC.
            # If a raw TIFF cannot be decoded, load_marker_image_from_plane returns
            # None and we still save the reconstructed label masks.
            imgs = load_and_stitch_section_images(pair, args)
            save_reconstructed_section_tifs(
                pair,
                imgs,
                {"DAPI": dapi_lbl, "cFOS": cfos_lbl, "NPY": npy_lbl, "POMC": pomc_lbl_clean},
                outdir,
                args,
            )
            if region_mask is not None:
                save_region_mask_tif(pair, region_mask, outdir, args)
        except Exception:
            logger.exception("Failed saving reconstructed TIFFs/labels for %s S%02d", pair.sample, pair.section_index)

    if bool(args.save_qc):
        try:
            qc_path = outdir / args.reconstructed_qc_dirname / pair.animal_id / f"{pair.sample}_S{pair.section_index:02d}{('_T'+pair.tile) if pair.tile else ''}_reconstructed_qc.pdf"
            render_reconstructed_qc_pdf(qc_path, dapi_lbl=dapi_lbl, section_row=base_section_row, cell_df=section_df, region_map=region_map, args=args, region_mask=region_mask, region_info=region_info)
        except Exception:
            logger.exception("Failed reconstructed QC for %s S%02d", pair.sample, pair.section_index)

    return cell_rows, base_section_row, assignment_rows


def attach_metadata_by_animal(df: Any, samplesheet: Any) -> Any:
    if df is None or df.empty:
        return df
    ss = samplesheet.copy()
    ss = ss.drop_duplicates("sample_key", keep="first")
    out = df.copy()
    if "animal_key" not in out.columns:
        out["animal_key"] = out["animal_id"].map(canonical_sample_id)
    meta_cols = [c for c in ss.columns if c != "sample"]
    out = out.merge(ss[meta_cols], on="sample_key" if "sample_key" in out.columns else "animal_key", how="left", suffixes=("", "_samplesheet")) if False else out
    # Explicit merge on animal_key -> sample_key avoids confusion with source sample keys from repetition folders.
    ss2 = ss.rename(columns={"sample_key": "animal_key"})
    keep_cols = [c for c in ss2.columns if c != "sample"]
    out = out.merge(ss2[keep_cols], on="animal_key", how="left", suffixes=("", "_samplesheet"))
    return out


def build_section_region_area_table(section_df: Any, region_order: List[str], args: argparse.Namespace) -> Any:
    """Long-form ARC/ME area table, including marker-assayed area only."""
    cols = [
        "sample", "animal_id", "animal_key", "section_index", "tile", "region",
        "region_area_px", "region_area_um2",
        "pomc_assay_area_px", "pomc_assay_area_um2",
        "npy_assay_area_px", "npy_assay_area_um2",
        "has_pomc_channel", "pomc_channel_complete",
        "has_npy_channel", "npy_channel_complete",
    ]
    if section_df is None or section_df.empty:
        return pd.DataFrame(columns=cols)

    pixel_area_um2 = float(args.um_per_px) ** 2
    rows: List[Dict[str, Any]] = []
    for rec in section_df.to_dict("records"):
        if str(rec.get("status", "ok")) != "ok":
            continue
        has_pomc = bool(rec.get("has_pomc_channel", False))
        has_npy = bool(rec.get("has_npy_channel", False))
        for region in region_order:
            area_px = float(rec.get(f"{str(region).lower()}_area_px", 0) or 0)
            rows.append({
                "sample": rec.get("sample", ""),
                "animal_id": rec.get("animal_id", ""),
                "animal_key": rec.get("animal_key", ""),
                "section_index": rec.get("section_index", 0),
                "tile": rec.get("tile", ""),
                "region": str(region),
                "region_area_px": area_px,
                "region_area_um2": area_px * pixel_area_um2,
                "pomc_assay_area_px": area_px if has_pomc else 0.0,
                "pomc_assay_area_um2": area_px * pixel_area_um2 if has_pomc else 0.0,
                "npy_assay_area_px": area_px if has_npy else 0.0,
                "npy_assay_area_um2": area_px * pixel_area_um2 if has_npy else 0.0,
                "has_pomc_channel": has_pomc,
                "pomc_channel_complete": bool(rec.get("pomc_channel_complete", False)),
                "has_npy_channel": has_npy,
                "npy_channel_complete": bool(rec.get("npy_channel_complete", False)),
            })
    return pd.DataFrame(rows, columns=cols)

def run_reconstructed_region_quantification(pairs: List[SectionPair], args: argparse.Namespace, logger: logging.Logger, outdir: Path, samplesheet: Any) -> Tuple[Any, Any, Any, Any, Any]:
    region_map = parse_region_bins(str(args.region_bins), int(args.n_bins))
    region_order = region_order_from_args(args, region_map)
    # Save bin map for reproducibility.
    pd.DataFrame([{"bin": b, "region": region_map[b]} for b in sorted(region_map)]).to_csv(outdir / "bin_region_map.csv", index=False)

    all_cell_rows: List[Dict[str, Any]] = []
    section_rows: List[Dict[str, Any]] = []
    assignment_rows: List[Dict[str, Any]] = []
    for pair in pairs:
        try:
            cells, sec, assigns = process_reconstructed_section(pair, args, outdir, region_map, logger)
            all_cell_rows.extend(cells)
            section_rows.append(sec)
            assignment_rows.extend(assigns)
        except Exception as e:
            logger.exception("Reconstructed quantification failed for %s S%02d", pair.sample, pair.section_index)
            section_rows.append({
                "sample": pair.sample, "animal_id": pair.animal_id, "animal_key": pair.animal_key,
                "is_repetition": pair.is_repetition, "repetition_label": pair.repetition_label,
                "section_index": pair.section_index, "tile": pair.tile, "status": f"error:{e}"
            })

    cells_df = pd.DataFrame(all_cell_rows)
    if not cells_df.empty:
        cells_df = attach_metadata_by_animal(cells_df, samplesheet)
    cells_df.to_csv(outdir / "per_cell_reconstructed_measurements.csv", index=False)

    section_df = pd.DataFrame(section_rows)
    if not section_df.empty:
        add_ratio_columns(section_df)
        section_df = attach_metadata_by_animal(section_df, samplesheet)
    section_df.to_csv(outdir / "per_reconstructed_section_summary.csv", index=False)

    assign_df = pd.DataFrame(assignment_rows)
    if not assign_df.empty:
        assign_df = attach_metadata_by_animal(assign_df, samplesheet)
    assign_df.to_csv(outdir / "cfos_object_assignments_reconstructed.csv", index=False)

    # Per-section / per-bin and per-region counts.
    section_bin_groups: List[Dict[str, Any]] = []
    section_region_groups: List[Dict[str, Any]] = []
    for sec in section_rows:
        if str(sec.get("status", "ok")) != "ok":
            continue
        for b in range(int(args.n_bins)):
            section_bin_groups.append({
                "sample": sec["sample"], "animal_id": sec["animal_id"], "animal_key": sec["animal_key"],
                "section_index": sec["section_index"], "tile": sec.get("tile", ""),
                "bin": b, "region": region_map.get(b, "OTHER"),
            })
        for r in region_order:
            section_region_groups.append({
                "sample": sec["sample"], "animal_id": sec["animal_id"], "animal_key": sec["animal_key"],
                "section_index": sec["section_index"], "tile": sec.get("tile", ""),
                "region": r,
            })

    section_region_area = build_section_region_area_table(section_df, region_order, args)
    section_region_area.to_csv(outdir / "per_reconstructed_section_region_area.csv", index=False)

    per_section_bin = summarize_cell_rows(cells_df, ["sample", "animal_id", "animal_key", "section_index", "tile", "bin", "region"], section_bin_groups)
    per_section_region = summarize_cell_rows(cells_df, ["sample", "animal_id", "animal_key", "section_index", "tile", "region"], section_region_groups)
    if not per_section_region.empty and not section_region_area.empty:
        merge_cols = ["sample", "animal_id", "animal_key", "section_index", "tile", "region"]
        area_cols = merge_cols + [
            "region_area_px", "region_area_um2", "pomc_assay_area_px", "pomc_assay_area_um2",
            "npy_assay_area_px", "npy_assay_area_um2",
            "has_pomc_channel", "pomc_channel_complete", "has_npy_channel", "npy_channel_complete",
        ]
        per_section_region = per_section_region.merge(section_region_area[area_cols], on=merge_cols, how="left")
        add_ratio_columns(per_section_region)
        add_area_normalized_columns(per_section_region)
    if not per_section_bin.empty:
        per_section_bin = attach_metadata_by_animal(per_section_bin, samplesheet)
    if not per_section_region.empty:
        per_section_region = attach_metadata_by_animal(per_section_region, samplesheet)
    per_section_bin.to_csv(outdir / "per_reconstructed_section_bin_summary.csv", index=False)
    per_section_region.to_csv(outdir / "per_reconstructed_section_region_summary.csv", index=False)

    # Per-animal total: source repetitions are summed into the same animal_key, not counted as separate N.
    if not cells_df.empty:
        animal_df = summarize_cell_rows(cells_df, ["animal_id", "animal_key"])
        src = section_df.groupby(["animal_id", "animal_key"], dropna=False).agg(
            source_samples=("sample", lambda x: ";".join(sorted(set(map(str, x))))),
            n_source_folders=("sample", lambda x: len(set(map(str, x)))),
            n_reconstructed_sections=("section_index", "count"),
            n_repetition_sections=("is_repetition", lambda x: int(pd.Series(x).astype(bool).sum())),
            n_pomc_sections=("has_pomc_channel", lambda x: int(pd.Series(x).astype(bool).sum())),
            n_npy_sections=("has_npy_channel", lambda x: int(pd.Series(x).astype(bool).sum())),
        ).reset_index()
        animal_df = animal_df.merge(src, on=["animal_id", "animal_key"], how="left")
        if not section_region_area.empty:
            total_area = section_region_area.groupby(["animal_id", "animal_key"], dropna=False).agg(
                analyzed_area_px=("region_area_px", "sum"),
                analyzed_area_um2=("region_area_um2", "sum"),
                pomc_assay_area_px=("pomc_assay_area_px", "sum"),
                pomc_assay_area_um2=("pomc_assay_area_um2", "sum"),
                npy_assay_area_px=("npy_assay_area_px", "sum"),
                npy_assay_area_um2=("npy_assay_area_um2", "sum"),
            ).reset_index()
            animal_df = animal_df.merge(total_area, on=["animal_id", "animal_key"], how="left")
        add_ratio_columns(animal_df)
        add_area_normalized_columns(animal_df, overall_area_col="analyzed_area_um2", pomc_area_col="pomc_assay_area_um2", npy_area_col="npy_assay_area_um2")
    else:
        animal_df = pd.DataFrame()
    if not animal_df.empty:
        animal_df = attach_metadata_by_animal(animal_df, samplesheet)
    animal_df.to_csv(outdir / "per_animal_total_summary.csv", index=False)
    animal_df.to_csv(outdir / "cfos_npy_pomc_summary_reconstructed_total_per_animal.csv", index=False)
    animal_df.to_csv(outdir / "cfos_npy_pomc_summary.csv", index=False)

    # Per-animal by bin/region.
    if not cells_df.empty:
        animal_bin_groups = []
        animal_region_groups = []
        animal_keys = cells_df[["animal_id", "animal_key"]].drop_duplicates().to_dict("records")
        for a in animal_keys:
            for b in range(int(args.n_bins)):
                animal_bin_groups.append({"animal_id": a["animal_id"], "animal_key": a["animal_key"], "bin": b, "region": region_map.get(b, "OTHER")})
            for r in region_order:
                animal_region_groups.append({"animal_id": a["animal_id"], "animal_key": a["animal_key"], "region": r})
        animal_bin_df = summarize_cell_rows(cells_df, ["animal_id", "animal_key", "bin", "region"], animal_bin_groups)
        animal_region_df = summarize_cell_rows(cells_df, ["animal_id", "animal_key", "region"], animal_region_groups)
        if not animal_region_df.empty and not section_region_area.empty:
            region_area = section_region_area.groupby(["animal_id", "animal_key", "region"], dropna=False).agg(
                region_area_px=("region_area_px", "sum"),
                region_area_um2=("region_area_um2", "sum"),
                pomc_assay_area_px=("pomc_assay_area_px", "sum"),
                pomc_assay_area_um2=("pomc_assay_area_um2", "sum"),
                npy_assay_area_px=("npy_assay_area_px", "sum"),
                npy_assay_area_um2=("npy_assay_area_um2", "sum"),
                n_sections_analyzed=("section_index", "count"),
                n_pomc_sections=("has_pomc_channel", lambda x: int(pd.Series(x).astype(bool).sum())),
                n_complete_pomc_sections=("pomc_channel_complete", lambda x: int(pd.Series(x).astype(bool).sum())),
                n_npy_sections=("has_npy_channel", lambda x: int(pd.Series(x).astype(bool).sum())),
                n_complete_npy_sections=("npy_channel_complete", lambda x: int(pd.Series(x).astype(bool).sum())),
            ).reset_index()
            animal_region_df = animal_region_df.merge(region_area, on=["animal_id", "animal_key", "region"], how="left")
            add_ratio_columns(animal_region_df)
            add_area_normalized_columns(animal_region_df)
    else:
        animal_bin_df = pd.DataFrame()
        animal_region_df = pd.DataFrame()
    if not animal_bin_df.empty:
        animal_bin_df = attach_metadata_by_animal(animal_bin_df, samplesheet)
    if not animal_region_df.empty:
        animal_region_df = attach_metadata_by_animal(animal_region_df, samplesheet)
    animal_bin_df.to_csv(outdir / "per_animal_bin_summary.csv", index=False)
    animal_region_df.to_csv(outdir / "per_animal_region_summary.csv", index=False)

    return cells_df, section_df, animal_df, animal_region_df, animal_bin_df


# -- statistics and plots for reconstructed total-per-animal data ----------------

LEGACY_REGION_METRIC_SPECS = [
    ("cfos_over_dapi", "c-FOS/DAPI ratio", "c-FOS/DAPI ratios per animal"),
    ("cfos_cells", "c-FOS+ cells", "c-FOS+ cells per animal"),
    ("cfos_only", "c-FOS-only cells", "c-FOS-only cells per animal"),
    ("cfos_pomc_all", "c-FOS∧POMC cells", "Activated POMC cells per animal"),
    ("cfos_pomc_over_pomc", "c-FOS∧POMC / total POMC", "POMC activation ratio per animal"),
    ("cfos_pomc_over_cfos", "c-FOS∧POMC / c-FOS", "Fraction of c-FOS cells that are POMC"),
]

# Exact compact two-panel ME/ARC endpoint suite requested for the final analysis.
FINAL_ARC_ME_METRIC_SPECS = [
    ("cfos_over_dapi", "c-FOS/DAPI ratio (per animal)", "c-FOS/DAPI Ratios per Animal"),
    ("cfos_pomc_over_dapi", "c-FOS∧POMC/DAPI ratio (per animal)", "c-FOS∧POMC/DAPI Ratios per Animal"),
    ("cfos_pomc_over_pomc", "c-FOS∧POMC/total POMC ratio (per animal)", "c-FOS∧POMC/Total POMC Ratios per Animal"),
    ("cfos_pomc_all", "c-FOS∧POMC cells (raw count per animal)", "c-FOS∧POMC Raw Counts per Animal"),
    ("cfos_pomc_per_um2", "c-FOS∧POMC density", "c-FOS∧POMC Density per Animal"),
    ("pomc_cells", "total POMC cells (raw count per animal)", "Total POMC Raw Counts per Animal"),
    ("pomc_per_um2", "total POMC density", "Total POMC Density per Animal"),
]

REGION_METRIC_SPECS = FINAL_ARC_ME_METRIC_SPECS
POMC_ASSAY_METRICS = {
    "cfos_pomc_over_dapi", "cfos_pomc_over_pomc", "cfos_pomc_over_cfos",
    "cfos_pomc_all", "cfos_pomc_per_um2", "pomc_cells", "pomc_per_um2",
}


def is_control_condition(c: Any) -> bool:
    return clean_condition(c) == "Control" or norm_token(c) in {"control", "negativectr", "negativecontrol", "ctrn", "noab"}


def _clean_conditions_for_plot(df: Any, exclude_control: bool = True) -> Any:
    out = df.copy()
    if "cond" not in out.columns:
        out["cond"] = "Unknown"
    out["cond_clean"] = out["cond"].map(clean_condition)
    if exclude_control:
        out = out[~out["cond_clean"].map(is_control_condition)].copy()
    return out


def _values_for_metric(df: Any, metric: str) -> Any:
    if metric not in df.columns:
        return pd.Series([np.nan] * len(df), index=df.index, dtype=float)
    values = pd.to_numeric(df[metric], errors="coerce").replace([np.inf, -np.inf], np.nan)

    # Do not treat animals without a POMC-assayed ROI as biological zeroes.
    if metric in POMC_ASSAY_METRICS:
        if "pomc_assay_area_um2" in df.columns:
            valid = pd.to_numeric(df["pomc_assay_area_um2"], errors="coerce").fillna(0) > 0
            values = values.where(valid)
        elif "n_pomc_sections" in df.columns:
            valid = pd.to_numeric(df["n_pomc_sections"], errors="coerce").fillna(0) > 0
            values = values.where(valid)
        if metric == "cfos_pomc_over_pomc" and "pomc_cells" in df.columns:
            values = values.where(pd.to_numeric(df["pomc_cells"], errors="coerce").fillna(0) > 0)
    elif metric == "cfos_over_dapi":
        if "region_area_um2" in df.columns:
            values = values.where(pd.to_numeric(df["region_area_um2"], errors="coerce").fillna(0) > 0)
        if "dapi_nuclei" in df.columns:
            values = values.where(pd.to_numeric(df["dapi_nuclei"], errors="coerce").fillna(0) > 0)
    return values


def _group_arrays_for_metric(df: Any, metric: str, conds: List[str], filtered: bool, args: argparse.Namespace) -> Tuple[List[Any], Dict[str, Any], Dict[str, Any]]:
    groups = []
    out_masks: Dict[str, Any] = {}
    raw_vals: Dict[str, Any] = {}
    for c in conds:
        vals = _values_for_metric(df[df["cond_clean"] == c], metric).dropna().to_numpy(dtype=float)
        out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
        raw_vals[c] = vals
        out_masks[c] = out
        groups.append(vals[~out] if filtered else vals)
    return groups, raw_vals, out_masks


def write_reconstructed_outlier_and_stats(animal_df: Any, animal_region_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    if not bool(args.write_stat_tables):
        return
    if animal_region_df is None or animal_region_df.empty:
        return

    rows_out: List[Dict[str, Any]] = []
    stats_all: List[Dict[str, Any]] = []
    stats_iqr: List[Dict[str, Any]] = []

    df = animal_region_df.copy()
    df = df[df["region"].astype(str).isin(["ME", "ARC"])].copy()
    df = _clean_conditions_for_plot(df, bool(args.exclude_control_plots))
    if df.empty:
        return

    conds = condition_order(df["cond_clean"])
    specs = FINAL_ARC_ME_METRIC_SPECS if str(getattr(args, "plot_suite", "final")) == "final" else LEGACY_REGION_METRIC_SPECS + FINAL_ARC_ME_METRIC_SPECS
    seen: Set[str] = set()
    specs = [s for s in specs if not (s[0] in seen or seen.add(s[0]))]

    for metric, _ylabel, _title in specs:
        for region in ["ME", "ARC"]:
            rdf = df[df["region"].astype(str) == region].copy()
            if rdf.empty:
                continue
            for c in conds:
                vals = _values_for_metric(rdf[rdf["cond_clean"] == c], metric).dropna().to_numpy(dtype=float)
                out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method))
                for v, o in zip(vals, out):
                    rows_out.append({"metric": metric, "region": region, "cond": c, "value": float(v), "is_outlier_iqr": bool(o)})
            for filtered, target in [(False, stats_all), (True, stats_iqr)]:
                groups, _raw_vals, _out_masks = _group_arrays_for_metric(rdf, metric, conds, filtered, args)
                target.append({"metric": metric, "region": region, "test": "ANOVA", "group_a": "ALL", "group_b": "", "pvalue": safe_anova(groups), "filtered_iqr": filtered})
                for i, c in enumerate(conds):
                    vals_i = groups[i]
                    target.append({
                        "metric": metric, "region": region, "test": "summary", "group_a": c, "group_b": "",
                        "n": int(vals_i.size), "mean": float(np.nanmean(vals_i)) if vals_i.size else float("nan"),
                        "sd": float(np.nanstd(vals_i, ddof=1)) if vals_i.size >= 2 else float("nan"),
                        "pvalue": float("nan"), "filtered_iqr": filtered,
                    })
                for i in range(len(conds)):
                    for j in range(i + 1, len(conds)):
                        target.append({
                            "metric": metric, "region": region, "test": "Mann-Whitney U",
                            "group_a": conds[i], "group_b": conds[j],
                            "pvalue": safe_mwu(groups[i], groups[j]), "filtered_iqr": filtered,
                        })

    pd.DataFrame(rows_out).to_csv(outdir / "outlier_report_iqr_reconstructed.csv", index=False)
    pd.DataFrame(stats_all).to_csv(outdir / "stats_condition_tests_all_values.csv", index=False)
    pd.DataFrame(stats_iqr).to_csv(outdir / "stats_condition_tests_iqr_filtered.csv", index=False)
    logger.warning("Final ARC/ME statistical tables written for %d metrics.", len(specs))



# -- definitive animal-level endpoint/statistics tables -------------------------

def _bh_adjust(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg FDR adjustment; NaNs remain NaN."""
    arr = np.asarray(pvals, dtype=float)
    out = np.full(arr.shape, np.nan, dtype=float)
    finite = np.isfinite(arr)
    if not np.any(finite):
        return out.tolist()
    idx = np.where(finite)[0]
    p = arr[idx]
    order = np.argsort(p)
    ranked = p[order]
    m = float(len(ranked))
    q = ranked * m / (np.arange(len(ranked), dtype=float) + 1.0)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    out[idx[order]] = q
    return out.tolist()


def _add_bh_columns(df: Any) -> Any:
    if df is None or df.empty or "p_value" not in df.columns:
        return df
    out = df.copy()
    out["q_value_bh_all_definitive_tests"] = np.nan
    out["q_value_bh_within_family_region"] = np.nan
    out["q_value_bh_within_endpoint_region"] = np.nan
    mask = pd.to_numeric(out["p_value"], errors="coerce").notna() & ~out["test"].astype(str).str.contains("summary", case=False, na=False)
    if mask.any():
        out.loc[mask, "q_value_bh_all_definitive_tests"] = _bh_adjust(pd.to_numeric(out.loc[mask, "p_value"], errors="coerce").tolist())
        for _keys, sub_idx in out[mask].groupby(["stats_values", "endpoint_family", "region"], dropna=False).groups.items():
            sub_idx = list(sub_idx)
            out.loc[sub_idx, "q_value_bh_within_family_region"] = _bh_adjust(pd.to_numeric(out.loc[sub_idx, "p_value"], errors="coerce").tolist())
        for _keys, sub_idx in out[mask].groupby(["stats_values", "endpoint", "region"], dropna=False).groups.items():
            sub_idx = list(sub_idx)
            out.loc[sub_idx, "q_value_bh_within_endpoint_region"] = _bh_adjust(pd.to_numeric(out.loc[sub_idx, "p_value"], errors="coerce").tolist())
    return out


def definitive_endpoint_specs(args: argparse.Namespace) -> List[Dict[str, Any]]:
    scale = float(max(getattr(args, "definitive_density_scale_um2", 100000.0), 1e-12))
    scale_label = f"per {scale:g} µm²"
    return [
        {"endpoint":"primary_pomc_activation_fraction","source_col":"cfos_pomc_over_pomc","endpoint_label":"Activated POMC fraction: c-FOS∧POMC / total POMC","endpoint_family":"POMC_activation","statistical_role":"PRIMARY","numerator_col":"cfos_pomc_all","denominator_col":"pomc_cells","area_col":"pomc_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"pomc_cells","recommended_model":"Primary animal-level fraction; confirm with binomial/quasi-binomial count model using activated POMC successes and total POMC denominator.","biological_question":"Among POMC neurons, what fraction is activated?"},
        {"endpoint":"activated_pomc_per_dapi","source_col":"cfos_pomc_over_dapi","endpoint_label":"Activated POMC / DAPI nuclei","endpoint_family":"POMC_activation","statistical_role":"SECONDARY_DAPI_NORMALIZED","numerator_col":"cfos_pomc_all","denominator_col":"pomc_assay_dapi_nuclei","area_col":"pomc_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"pomc_assay_dapi_nuclei","recommended_model":"Animal-level fraction; DAPI-normalized activated POMC abundance.","biological_question":"How much of the local cellular population is activated POMC?"},
        {"endpoint":"activated_pomc_density","source_col":"cfos_pomc_per_um2","endpoint_label":f"Activated POMC density ({scale_label})","endpoint_family":"POMC_activation","statistical_role":"SECONDARY_AREA_NORMALIZED","numerator_col":"cfos_pomc_all","denominator_col":"","area_col":"pomc_assay_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"pomc_assay_area_um2","recommended_model":"Count model with log(area_um2) offset; animal-level density shown for inference.","biological_question":"Activated POMC cells per analyzed tissue area."},
        {"endpoint":"total_cfos_activation_fraction","source_col":"cfos_over_dapi","endpoint_label":"Total c-FOS / DAPI nuclei","endpoint_family":"global_cFOS_activation","statistical_role":"SECONDARY_TOTAL_ACTIVATION","numerator_col":"cfos_cells","denominator_col":"dapi_nuclei","area_col":"region_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"dapi_nuclei","recommended_model":"Animal-level fraction; confirm with binomial/quasi-binomial count model using c-FOS successes and DAPI denominator.","biological_question":"Overall cell activation in the region."},
        {"endpoint":"total_cfos_density","source_col":"cfos_per_um2","endpoint_label":f"Total c-FOS density ({scale_label})","endpoint_family":"global_cFOS_activation","statistical_role":"SECONDARY_AREA_NORMALIZED","numerator_col":"cfos_cells","denominator_col":"","area_col":"region_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"region_area_um2","recommended_model":"Count model with log(area_um2) offset; animal-level density shown for inference.","biological_question":"Total activated cells per analyzed tissue area."},
        {"endpoint":"pomc_fraction_of_dapi","source_col":"pomc_over_dapi","endpoint_label":"Total POMC / DAPI nuclei","endpoint_family":"POMC_population_QC","statistical_role":"QC_DENOMINATOR_POPULATION","numerator_col":"pomc_cells","denominator_col":"pomc_assay_dapi_nuclei","area_col":"pomc_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"pomc_assay_dapi_nuclei","recommended_model":"QC endpoint; verify comparable POMC sampling across conditions.","biological_question":"Is the POMC population/sampling comparable?"},
        {"endpoint":"pomc_density","source_col":"pomc_per_um2","endpoint_label":f"Total POMC density ({scale_label})","endpoint_family":"POMC_population_QC","statistical_role":"QC_DENOMINATOR_POPULATION","numerator_col":"pomc_cells","denominator_col":"","area_col":"pomc_assay_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"pomc_assay_area_um2","recommended_model":"QC endpoint; count model with log(POMC-assayed area) offset if needed.","biological_question":"POMC sampling density."},
        {"endpoint":"activated_npy_fraction","source_col":"cfos_npy_over_npy","endpoint_label":"Activated NPY fraction: c-FOS∧NPY / total NPY","endpoint_family":"NPY_activation","statistical_role":"OPTIONAL_SECONDARY_IF_NPY_PRESENT","numerator_col":"cfos_npy_all","denominator_col":"npy_cells","area_col":"npy_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"npy_cells","recommended_model":"Animal-level fraction; only for animals/sections with NPY channel.","biological_question":"Among NPY neurons, what fraction is activated?"},
        {"endpoint":"activated_npy_per_dapi","source_col":"cfos_npy_over_dapi","endpoint_label":"Activated NPY / DAPI nuclei","endpoint_family":"NPY_activation","statistical_role":"OPTIONAL_SECONDARY_IF_NPY_PRESENT","numerator_col":"cfos_npy_all","denominator_col":"npy_assay_dapi_nuclei","area_col":"npy_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"npy_assay_dapi_nuclei","recommended_model":"Animal-level fraction; only for animals/sections with NPY channel.","biological_question":"How much of the local cellular population is activated NPY?"},
        {"endpoint":"activated_npy_density","source_col":"cfos_npy_per_um2","endpoint_label":f"Activated NPY density ({scale_label})","endpoint_family":"NPY_activation","statistical_role":"OPTIONAL_SECONDARY_IF_NPY_PRESENT","numerator_col":"cfos_npy_all","denominator_col":"","area_col":"npy_assay_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"npy_assay_area_um2","recommended_model":"Count model with log(NPY-assayed area) offset if needed.","biological_question":"Activated NPY cells per analyzed tissue area."},
        {"endpoint":"npy_fraction_of_dapi","source_col":"npy_over_dapi","endpoint_label":"Total NPY / DAPI nuclei","endpoint_family":"NPY_population_QC","statistical_role":"QC_DENOMINATOR_POPULATION","numerator_col":"npy_cells","denominator_col":"npy_assay_dapi_nuclei","area_col":"npy_assay_area_um2","display_scale":1.0,"unit":"fraction","valid_denominator_col":"npy_assay_dapi_nuclei","recommended_model":"QC endpoint; verify comparable NPY sampling across conditions.","biological_question":"NPY population/sampling density by nuclei."},
        {"endpoint":"npy_density","source_col":"npy_per_um2","endpoint_label":f"Total NPY density ({scale_label})","endpoint_family":"NPY_population_QC","statistical_role":"QC_DENOMINATOR_POPULATION","numerator_col":"npy_cells","denominator_col":"","area_col":"npy_assay_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"npy_assay_area_um2","recommended_model":"QC endpoint; count model with log(NPY-assayed area) offset if needed.","biological_question":"NPY sampling density."},
        {"endpoint":"dapi_density","source_col":"dapi_per_um2","endpoint_label":f"DAPI density ({scale_label})","endpoint_family":"tissue_sampling_QC","statistical_role":"QC_TISSUE_CELLULARITY","numerator_col":"dapi_nuclei","denominator_col":"","area_col":"region_area_um2","display_scale":scale,"unit":f"cells/{scale:g}um2","valid_denominator_col":"region_area_um2","recommended_model":"QC endpoint; should not be the primary activation endpoint.","biological_question":"Cellularity / sampling density."},
    ]


def build_definitive_endpoint_values(animal_region_df: Any, args: argparse.Namespace) -> Any:
    if animal_region_df is None or animal_region_df.empty:
        return pd.DataFrame()
    df = animal_region_df.copy()
    df = df[df["region"].astype(str).isin(["ME", "ARC"])].copy()
    df = _clean_conditions_for_plot(df, bool(args.exclude_control_plots))
    if df.empty:
        return pd.DataFrame()
    meta_cols = ["animal_id","animal_key","region","cond","cond_clean","genotype","sex","cage","source_samples","n_source_folders","n_sections_analyzed","n_reconstructed_sections","n_pomc_sections","n_npy_sections","n_complete_pomc_sections","n_complete_npy_sections"]
    count_cols = ["dapi_nuclei","cfos_cells","cfos_only","pomc_cells","npy_cells","cfos_pomc_all","cfos_npy_all","cfos_npy_pomc","pomc_assay_dapi_nuclei","npy_assay_dapi_nuclei","region_area_um2","pomc_assay_area_um2","npy_assay_area_um2"]
    rows: List[Dict[str, Any]] = []
    for rec in df.to_dict("records"):
        for spec in definitive_endpoint_specs(args):
            source = str(spec["source_col"])
            raw_val = rec.get(source, np.nan)
            try:
                raw_val = float(raw_val)
            except Exception:
                raw_val = float("nan")
            value = raw_val * float(spec.get("display_scale", 1.0)) if np.isfinite(raw_val) else float("nan")
            num_col = str(spec.get("numerator_col", "")); den_col = str(spec.get("denominator_col", "")); area_col = str(spec.get("area_col", "")); valid_col = str(spec.get("valid_denominator_col", ""))
            valid_basis = rec.get(valid_col, np.nan) if valid_col else np.nan
            try:
                valid_basis_f = float(valid_basis)
            except Exception:
                valid_basis_f = float("nan")
            valid = bool(np.isfinite(value) and np.isfinite(valid_basis_f) and valid_basis_f > 0)
            row: Dict[str, Any] = {c: rec.get(c, np.nan) for c in meta_cols if c in rec}
            for c in count_cols:
                row[c] = rec.get(c, np.nan)
            row.update({
                "endpoint": spec["endpoint"], "endpoint_label": spec["endpoint_label"], "endpoint_family": spec["endpoint_family"],
                "statistical_role": spec["statistical_role"], "biological_question": spec["biological_question"], "recommended_model": spec["recommended_model"],
                "source_col": source, "numerator_col": num_col, "denominator_col": den_col, "area_col": area_col,
                "numerator": rec.get(num_col, np.nan) if num_col else np.nan,
                "denominator": rec.get(den_col, np.nan) if den_col else np.nan,
                "area_um2": rec.get(area_col, np.nan) if area_col else np.nan,
                "raw_value": raw_val, "display_value": value, "display_unit": spec["unit"],
                "density_scale_um2": float(spec.get("display_scale", 1.0)),
                "valid_for_inference": valid,
                "exclusion_reason": "" if valid else f"invalid_or_missing_{valid_col or source}",
            })
            rows.append(row)
    return pd.DataFrame(rows)


def _endpoint_summary_rows(vals: Any, cond: str) -> Dict[str, Any]:
    arr = np.asarray(vals, dtype=float); arr = arr[np.isfinite(arr)]
    return {"condition_a": cond, "condition_b": "", "n_a": int(arr.size), "mean_a": float(np.nanmean(arr)) if arr.size else float("nan"), "sd_a": float(np.nanstd(arr, ddof=1)) if arr.size >= 2 else float("nan"), "sem_a": float(np.nanstd(arr, ddof=1) / math.sqrt(arr.size)) if arr.size >= 2 else float("nan"), "median_a": float(np.nanmedian(arr)) if arr.size else float("nan"), "q1_a": float(np.nanpercentile(arr, 25)) if arr.size else float("nan"), "q3_a": float(np.nanpercentile(arr, 75)) if arr.size else float("nan")}


def build_definitive_stats_table(endpoint_values: Any, args: argparse.Namespace, *, filtered: bool) -> Any:
    rows: List[Dict[str, Any]] = []
    if endpoint_values is None or endpoint_values.empty:
        return pd.DataFrame()
    ev = endpoint_values[endpoint_values["valid_for_inference"].astype(bool)].copy()
    if ev.empty:
        return pd.DataFrame()
    for (endpoint, region), sdf in ev.groupby(["endpoint", "region"], dropna=False):
        first = sdf.iloc[0].to_dict(); conds = condition_order(sdf["cond_clean"])
        groups_use: Dict[str, Any] = {}; out_count: Dict[str, int] = {}
        for c in conds:
            vals = pd.to_numeric(sdf.loc[sdf["cond_clean"] == c, "display_value"], errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy(dtype=float)
            out = iqr_outlier_mask(vals, k=float(args.iqr_k), method=str(args.outlier_method)) if str(args.outlier_method) == "iqr" else np.zeros(vals.shape, dtype=bool)
            groups_use[c] = vals[~out] if filtered else vals; out_count[c] = int(np.count_nonzero(out)) if vals.size else 0
        group_list = [groups_use[c] for c in conds]
        common = {"stats_values": "iqr_filtered" if filtered else "all_values", "endpoint": endpoint, "endpoint_label": first.get("endpoint_label", endpoint), "endpoint_family": first.get("endpoint_family", ""), "statistical_role": first.get("statistical_role", ""), "region": region, "display_unit": first.get("display_unit", ""), "recommended_model": first.get("recommended_model", ""), "biological_question": first.get("biological_question", "")}
        rows.append({**common, "test":"One-way ANOVA global", "p_value_method":"ordinary equal-variance parametric F distribution", "condition_a":"ALL", "condition_b":"", "n_groups":int(sum(1 for g in group_list if len(g) >= 2)), "p_value":safe_anova(group_list), "note":"Primary global test: animal is the experimental unit; interpret cautiously with small N."})
        for c in conds:
            rows.append({**common, "test":"summary", **_endpoint_summary_rows(groups_use[c], c), "n_outliers_removed_a":out_count.get(c, 0) if filtered else 0, "p_value":float("nan"), "note":"Animal-level endpoint summary."})
        for i in range(len(conds)):
            for j in range(i + 1, len(conds)):
                a, b = conds[i], conds[j]; va, vb = groups_use[a], groups_use[b]
                rows.append({**common, "test":"Mann-Whitney U pairwise", "p_value_method":"exhaustive two-sided independent-label permutation of rank-sum/U", "condition_a":a, "condition_b":b, "n_a":int(len(va)), "n_b":int(len(vb)), "mean_a":float(np.nanmean(va)) if len(va) else float("nan"), "mean_b":float(np.nanmean(vb)) if len(vb) else float("nan"), "median_a":float(np.nanmedian(va)) if len(va) else float("nan"), "median_b":float(np.nanmedian(vb)) if len(vb) else float("nan"), "median_difference_b_minus_a":(float(np.nanmedian(vb))-float(np.nanmedian(va))) if len(va) and len(vb) else float("nan"), "mean_difference_b_minus_a":(float(np.nanmean(vb))-float(np.nanmean(va))) if len(va) and len(vb) else float("nan"), "cliffs_delta_b_vs_a":cliffs_delta(va, vb), "n_outliers_removed_a":out_count.get(a,0) if filtered else 0, "n_outliers_removed_b":out_count.get(b,0) if filtered else 0, "p_value":safe_mwu(va, vb), "note":"Pairwise animal-level nonparametric comparison with observed ties retained."})
    return _add_bh_columns(pd.DataFrame(rows))


def write_definitive_statistics_csvs(animal_region_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    if not bool(getattr(args, "write_definitive_stats", True)):
        return
    if animal_region_df is None or animal_region_df.empty:
        logger.warning("Definitive stats skipped: no animal-region table.")
        return
    endpoint_values = build_definitive_endpoint_values(animal_region_df, args)
    endpoint_values.to_csv(outdir / "definitive_endpoint_values_long.csv", index=False)
    if not endpoint_values.empty:
        id_cols = [c for c in ["animal_id", "animal_key", "region", "cond", "cond_clean", "genotype", "sex", "cage"] if c in endpoint_values.columns]
        try:
            wide = endpoint_values.pivot_table(index=id_cols, columns="endpoint", values="display_value", aggfunc="first").reset_index()
            wide.columns = [str(c) for c in wide.columns]
            denom_cols = [c for c in ["dapi_nuclei","cfos_cells","pomc_cells","npy_cells","cfos_pomc_all","cfos_npy_all","pomc_assay_dapi_nuclei","npy_assay_dapi_nuclei","region_area_um2","pomc_assay_area_um2","npy_assay_area_um2","n_sections_analyzed","n_pomc_sections","n_npy_sections"] if c in endpoint_values.columns]
            denoms = endpoint_values[id_cols + denom_cols].drop_duplicates(id_cols)
            wide = wide.merge(denoms, on=id_cols, how="left")
            wide.to_csv(outdir / "definitive_endpoint_values_wide.csv", index=False)
        except Exception:
            logger.exception("Could not write definitive_endpoint_values_wide.csv")
    stats_all = build_definitive_stats_table(endpoint_values, args, filtered=False)
    stats_iqr = build_definitive_stats_table(endpoint_values, args, filtered=True)
    stats_all.to_csv(outdir / "definitive_statistics_all_values.csv", index=False)
    stats_iqr.to_csv(outdir / "definitive_statistics_iqr_filtered.csv", index=False)
    stats_combined = pd.concat([stats_all, stats_iqr], ignore_index=True) if not stats_all.empty or not stats_iqr.empty else pd.DataFrame()
    stats_combined.to_csv(outdir / "definitive_statistics.csv", index=False)
    plan = pd.DataFrame(definitive_endpoint_specs(args))
    plan = plan[["endpoint","endpoint_label","endpoint_family","statistical_role","numerator_col","denominator_col","area_col","unit","biological_question","recommended_model"]]
    plan.to_csv(outdir / "definitive_statistics_analysis_plan.csv", index=False)
    logger.warning("Definitive animal-level endpoint/statistics CSVs written to: %s", outdir)

def plot_reconstructed_metric_by_condition(df: Any, out_pdf: Path, metric: str, ylabel: str, title: str, args: argparse.Namespace, *, region: Optional[str] = None) -> None:
    plot_df = _clean_conditions_for_plot(df, bool(args.exclude_control_plots))
    if region is not None:
        plot_df = plot_df[plot_df["region"] == region].copy()
    if plot_df.empty:
        return
    conds = condition_order(plot_df["cond_clean"])
    if not conds:
        return
    groups, raw_vals, out_masks = _group_arrays_for_metric(plot_df, metric, conds, bool(args.exclude_outliers_stats), args)
    groups_all, raw_vals_all, out_masks_all = _group_arrays_for_metric(plot_df, metric, conds, False, args)
    colors = [COND_COLORS.get(c, "#cccccc") for c in conds]
    means = [float(np.nanmean(g)) if g.size else float("nan") for g in groups_all]
    errs = [float(np.nanstd(g, ddof=1)) if g.size >= 2 else 0.0 for g in groups_all]
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    x = np.arange(len(conds))
    ax.bar(x, means, yerr=errs, color=colors, edgecolor="black", linewidth=1.8, capsize=5, alpha=0.9)
    rng = np.random.default_rng(7)
    for i, c in enumerate(conds):
        vals = raw_vals_all.get(c, np.array([]))
        outs = out_masks_all.get(c, np.zeros(vals.shape, dtype=bool))
        jitter = rng.normal(0, 0.045, size=vals.size)
        ax.scatter(np.full(vals.size, x[i]) + jitter, vals, s=58, color=COND_COLORS.get(c, "#cccccc"), edgecolor="black", zorder=3)
        if bool(args.mark_outliers) and str(args.outlier_method) == "iqr" and vals.size:
            ax.scatter((np.full(vals.size, x[i]) + jitter)[outs], vals[outs], marker="x", s=105, linewidth=2.2, color="black", zorder=4)
        ax.text(x[i], means[i] + (errs[i] if np.isfinite(errs[i]) else 0), f"n={len(groups[i])}", ha="center", va="bottom", fontsize=10, weight="bold")
    p_anova = safe_anova(groups)
    stat_lines = ["stats: IQR-filtered" if bool(args.exclude_outliers_stats) else "stats: all values", f"ANOVA p={_fmt_p(p_anova)}"]
    for i in range(len(conds)):
        for j in range(i + 1, len(conds)):
            stat_lines.append(f"{conds[i]} vs {conds[j]} p={_fmt_p(safe_mwu(groups[i], groups[j]))}")
    ax.text(0.98, 0.98, "\n".join(stat_lines), transform=ax.transAxes, ha="right", va="top", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(conds)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}{' — ' + region if region else ''}", weight="bold")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def plot_reconstructed_region_facets(animal_region_df: Any, out_pdf: Path, metric: str, ylabel: str, title: str, args: argparse.Namespace, region_order: List[str]) -> None:
    if animal_region_df is None or animal_region_df.empty:
        return
    df = _clean_conditions_for_plot(animal_region_df, bool(args.exclude_control_plots))
    if df.empty:
        return
    regions = [r for r in region_order if r in set(df["region"].astype(str))]
    for r in sorted(set(df["region"].astype(str))):
        if r not in regions:
            regions.append(r)
    if not regions:
        return
    ncols = min(2, len(regions))
    nrows = int(math.ceil(len(regions) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 5.0 * nrows), squeeze=False)
    axes = axes.reshape(-1)
    fig.suptitle(title, fontsize=18, weight="bold", y=1.02)
    for ax, region in zip(axes, regions):
        sub = df[df["region"].astype(str) == region]
        conds = condition_order(sub["cond_clean"])
        groups, raw_vals, out_masks = _group_arrays_for_metric(sub, metric, conds, bool(args.exclude_outliers_stats), args)
        groups_all, raw_vals_all, out_masks_all = _group_arrays_for_metric(sub, metric, conds, False, args)
        x = np.arange(len(conds))
        means = [float(np.nanmean(g)) if g.size else float("nan") for g in groups_all]
        errs = [float(np.nanstd(g, ddof=1)) if g.size >= 2 else 0.0 for g in groups_all]
        colors = [COND_COLORS.get(c, "#cccccc") for c in conds]
        ax.bar(x, means, yerr=errs, capsize=5, color=colors, edgecolor="black", linewidth=1.5, alpha=0.9)
        rng = np.random.default_rng(17)
        for i, c in enumerate(conds):
            vals = raw_vals_all.get(c, np.array([]))
            outs = out_masks_all.get(c, np.zeros(vals.shape, dtype=bool))
            jitter = rng.normal(0, 0.045, size=vals.size)
            ax.scatter(np.full(vals.size, x[i]) + jitter, vals, s=50, color=COND_COLORS.get(c, "#cccccc"), edgecolor="black", zorder=3)
            if bool(args.mark_outliers) and str(args.outlier_method) == "iqr" and vals.size:
                ax.scatter((np.full(vals.size, x[i]) + jitter)[outs], vals[outs], marker="x", s=100, linewidth=2.1, color="black", zorder=4)
            if means[i] == means[i]:
                ax.text(x[i], means[i] + (errs[i] if np.isfinite(errs[i]) else 0), f"n={len(groups[i])}", ha="center", va="bottom", fontsize=9, weight="bold")
        p_anova = safe_anova(groups)
        stat_lines = [f"ANOVA p={_fmt_p(p_anova)}"]
        for i in range(len(conds)):
            for j in range(i + 1, len(conds)):
                stat_lines.append(f"{conds[i]} vs {conds[j]} p={_fmt_p(safe_mwu(groups[i], groups[j]))}")
        ax.text(0.98, 0.98, "\n".join(stat_lines), transform=ax.transAxes, ha="right", va="top", fontsize=8.5)
        ax.set_title(region, weight="bold")
        ax.set_ylabel(ylabel)
        ax.set_xticks(x); ax.set_xticklabels(conds)
        ax.grid(axis="y", alpha=0.22)
    for ax in axes[len(regions):]:
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def _density_scale_text(scale_um2: float) -> str:
    scale = float(scale_um2)
    if math.isclose(scale, 1.0):
        return "µm²"
    if math.isclose(scale, round(scale)):
        return f"{int(round(scale)):,} µm²"
    return f"{scale:g} µm²"


def _final_metric_display(metric: str, ylabel: str, title: str, args: argparse.Namespace) -> Tuple[str, str, float]:
    if metric.endswith("_per_um2"):
        scale = max(float(args.plot_density_area_um2), 1e-12)
        area_txt = _density_scale_text(scale)
        short = ylabel.replace(" density", "").replace("Density", "")
        return f"{short} / {area_txt}", title, scale
    return ylabel, title, 1.0


def _final_plot_formats(args: argparse.Namespace) -> List[str]:
    allowed = {"pdf", "png", "svg"}
    fmts = []
    for token in str(args.plot_formats).split(","):
        t = token.strip().lower().lstrip(".")
        if t in allowed and t not in fmts:
            fmts.append(t)
    return fmts or ["pdf", "png"]


def build_final_arc_me_figure(
    animal_region_df: Any,
    metric: str,
    ylabel: str,
    title: str,
    args: argparse.Namespace,
) -> Optional[Any]:
    """Create the compact, exact two-panel ME/ARC plot used for final reporting."""
    if animal_region_df is None or animal_region_df.empty:
        return None
    df = _clean_conditions_for_plot(animal_region_df, bool(args.exclude_control_plots))
    df = df[df["region"].astype(str).isin(["ME", "ARC"])].copy()
    if df.empty:
        return None

    conds = condition_order(df["cond_clean"])
    if not conds:
        return None

    ylabel_display, title_display, display_scale = _final_metric_display(metric, ylabel, title, args)
    axis_fs = float(args.plot_axis_fontsize)
    line_w = float(args.plot_linewidth)
    bar_width = float(args.plot_bar_width)
    err_mode = str(args.plot_error)

    plt.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.size": axis_fs,
        "axes.titlesize": axis_fs + 1.5,
        "axes.labelsize": axis_fs,
        "xtick.labelsize": axis_fs,
        "ytick.labelsize": axis_fs,
        "axes.linewidth": line_w,
    })

    fig, axes = plt.subplots(1, 2, figsize=(float(args.plot_figure_width), float(args.plot_figure_height)), squeeze=False)
    axes = axes.reshape(-1)
    fig.suptitle(title_display, fontsize=float(args.plot_title_fontsize), fontweight="bold", y=0.985)

    from matplotlib.lines import Line2D
    legend_handles: List[Any] = []
    for c in conds:
        legend_handles.append(Line2D(
            [0], [0], marker="o", linestyle="none",
            markerfacecolor=COND_COLORS.get(c, "#cccccc"), markeredgecolor="black",
            markeredgewidth=1.25, markersize=7.5, label=c,
        ))
    if bool(args.mark_outliers) and str(args.outlier_method) == "iqr":
        legend_handles.append(Line2D([0], [0], marker="x", linestyle="none", color="black", markeredgewidth=2.0, markersize=8.5, label="Outlier (IQR)"))
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.915),
        ncol=max(1, len(legend_handles)),
        frameon=False,
        handletextpad=0.35,
        columnspacing=0.75,
        prop={"size": float(args.plot_legend_fontsize), "weight": "bold"},
    )

    rng = np.random.default_rng(20260619)
    for ax, region in zip(axes, ["ME", "ARC"]):
        sub = df[df["region"].astype(str) == region].copy()
        groups_stats, _raw_stats, _out_stats = _group_arrays_for_metric(sub, metric, conds, bool(args.exclude_outliers_stats), args)
        _groups_all, raw_vals, out_masks = _group_arrays_for_metric(sub, metric, conds, False, args)

        x = np.arange(len(conds), dtype=float)
        display_vals: Dict[str, Any] = {c: raw_vals.get(c, np.array([], dtype=float)) * display_scale for c in conds}
        means = [float(np.nanmean(display_vals[c])) if display_vals[c].size else float("nan") for c in conds]
        if err_mode == "sem":
            errs = [float(np.nanstd(display_vals[c], ddof=1) / math.sqrt(display_vals[c].size)) if display_vals[c].size >= 2 else 0.0 for c in conds]
        else:
            errs = [float(np.nanstd(display_vals[c], ddof=1)) if display_vals[c].size >= 2 else 0.0 for c in conds]

        colors = [COND_COLORS.get(c, "#cccccc") for c in conds]
        ax.bar(
            x, means, width=bar_width, yerr=errs, capsize=4.5,
            color=colors, edgecolor="black", linewidth=line_w, alpha=0.92,
            error_kw={"elinewidth": line_w, "capthick": line_w},
            zorder=2,
        )

        finite_all: List[float] = []
        for i, c in enumerate(conds):
            vals = display_vals[c]
            outs = out_masks.get(c, np.zeros(vals.shape, dtype=bool))
            if vals.size:
                finite_all.extend(vals[np.isfinite(vals)].tolist())
                jitter = rng.normal(0.0, 0.035, size=vals.size)
                xs = np.full(vals.size, x[i], dtype=float) + jitter
                ax.scatter(
                    xs, vals, s=28, facecolor=colors[i], edgecolor="black",
                    linewidth=1.05, zorder=4,
                )
                if bool(args.mark_outliers) and str(args.outlier_method) == "iqr" and outs.size == vals.size and outs.any():
                    ax.scatter(xs[outs], vals[outs], marker="x", s=60, linewidth=1.9, color="black", zorder=5)

        data_max = max(finite_all) if finite_all else 1.0
        bar_max = max([m + e for m, e in zip(means, errs) if np.isfinite(m)] or [0.0])
        ymax = max(data_max, bar_max, 1e-12)
        upper = ymax * 1.55 + (0.02 * ymax if ymax > 0 else 0.1)
        if upper <= 0:
            upper = 1.0
        ax.set_ylim(0, upper)

        for i, c in enumerate(conds):
            vals = display_vals[c]
            if not vals.size or not np.isfinite(means[i]):
                continue
            y_n = means[i] + errs[i] + 0.035 * upper
            ax.text(x[i], y_n, f"n={len(groups_stats[i])}", ha="center", va="bottom", fontsize=axis_fs - 0.5, fontweight="bold")

        p_anova = safe_anova(groups_stats)
        stat_lines = [f"ANOVA p={_fmt_p(p_anova)}"]
        for i in range(len(conds)):
            for j in range(i + 1, len(conds)):
                stat_lines.append(f"MWU: {conds[i].lower()} vs {conds[j].lower()} p={_fmt_p(safe_mwu(groups_stats[i], groups_stats[j]))}")
        ax.text(
            0.985, 0.985, "\n".join(stat_lines), transform=ax.transAxes,
            ha="right", va="top", fontsize=float(args.plot_stats_fontsize), fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.0},
            zorder=10,
        )

        ax.set_title(region, fontweight="bold", pad=5)
        ax.set_xticks(x)
        ax.set_xticklabels(conds, fontweight="bold")
        ax.set_ylabel(ylabel_display, fontweight="bold")
        ax.grid(False)
        ax.tick_params(axis="both", width=line_w, length=4)
        for spine in ax.spines.values():
            spine.set_linewidth(line_w)
            spine.set_color("black")
        if metric.endswith("_per_um2") and display_scale == 1.0:
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)

    fig.subplots_adjust(left=0.115, right=0.985, bottom=0.19, top=0.72, wspace=0.30)
    return fig


def make_final_arc_me_plots(animal_region_df: Any, plot_dir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    ensure_dir(plot_dir)
    if animal_region_df is None or animal_region_df.empty:
        logger.warning("No per-animal ARC/ME rows available for final plots.")
        return

    plot_values = animal_region_df[animal_region_df["region"].astype(str).isin(["ME", "ARC"])].copy()
    scale = max(float(args.plot_density_area_um2), 1e-12)
    plot_values["cfos_endpoint_valid"] = pd.to_numeric(plot_values.get("region_area_um2", 0), errors="coerce").fillna(0) > 0
    plot_values["pomc_endpoint_valid"] = pd.to_numeric(plot_values.get("pomc_assay_area_um2", 0), errors="coerce").fillna(0) > 0
    for metric, _ylabel, _title in FINAL_ARC_ME_METRIC_SPECS:
        plot_values[f"plot_{metric}"] = _values_for_metric(plot_values, metric)
    plot_values[f"plot_cfos_pomc_per_{scale:g}_um2"] = plot_values["plot_cfos_pomc_per_um2"] * scale
    plot_values[f"plot_pomc_per_{scale:g}_um2"] = plot_values["plot_pomc_per_um2"] * scale
    plot_values.to_csv(plot_dir / "final_arc_me_plot_values.csv", index=False)

    from matplotlib.backends.backend_pdf import PdfPages
    combined_pdf = plot_dir / "ARC_ME_FINAL_ALL_ENDPOINTS.pdf"
    formats = _final_plot_formats(args)
    with PdfPages(combined_pdf) as pdf:
        for metric, ylabel, title in FINAL_ARC_ME_METRIC_SPECS:
            try:
                fig = build_final_arc_me_figure(animal_region_df, metric, ylabel, title, args)
                if fig is None:
                    continue
                slug_map = {
                    "cfos_over_dapi": "cfos_dapi_ratio",
                    "cfos_pomc_over_dapi": "cfos_pomc_dapi_ratio",
                    "cfos_pomc_over_pomc": "cfos_pomc_total_pomc_ratio",
                    "cfos_pomc_all": "cfos_pomc_raw_counts",
                    "cfos_pomc_per_um2": "cfos_pomc_density",
                    "pomc_cells": "total_pomc_raw_counts",
                    "pomc_per_um2": "total_pomc_density",
                }
                slug = slug_map.get(metric, re.sub(r"[^a-z0-9]+", "_", metric.lower()).strip("_"))
                for fmt in formats:
                    path = plot_dir / f"ARC_ME_{slug}.{fmt}"
                    save_kwargs = {"bbox_inches": "tight"}
                    if fmt == "png":
                        save_kwargs["dpi"] = int(args.plot_dpi)
                    fig.savefig(path, **save_kwargs)
                pdf.savefig(fig, bbox_inches="tight")
                plt.close(fig)
            except Exception:
                logger.exception("Failed final ARC/ME plot for metric %s", metric)
    logger.warning("Final compact ARC/ME plots written to %s", plot_dir)
    logger.warning("Combined multipage PDF: %s", combined_pdf)


def make_reconstructed_plots(animal_df: Any, animal_region_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    plot_dir = outdir / "plots"
    ensure_dir(plot_dir)

    # Default: only the compact two-panel ME/ARC figures requested for the final analysis.
    make_final_arc_me_plots(animal_region_df, plot_dir, args, logger)

    # Optional compatibility suite. This is OFF by default so the plots directory stays focused.
    if str(getattr(args, "plot_suite", "final")) == "all":
        legacy_dir = plot_dir / "archive_all_plots"
        ensure_dir(legacy_dir)
        region_map = parse_region_bins(str(args.region_bins), int(args.n_bins))
        reg_order = ["ME", "ARC"]
        for metric, ylabel, title in LEGACY_REGION_METRIC_SPECS:
            try:
                if animal_df is not None and not animal_df.empty:
                    plot_reconstructed_metric_by_condition(animal_df, legacy_dir / f"total_{metric}_by_condition.pdf", metric, ylabel, title + " — TOTAL", args)
                if animal_region_df is not None and not animal_region_df.empty:
                    plot_reconstructed_region_facets(animal_region_df, legacy_dir / f"region_{metric}_by_condition.pdf", metric, ylabel, title, args, reg_order)
            except Exception:
                logger.exception("Failed legacy reconstructed plot for metric %s", metric)

    # Group summary CSVs, including final ratios/densities.
    try:
        if animal_df is not None and not animal_df.empty and "cond" in animal_df.columns:
            candidates = [m for m, _, _ in FINAL_ARC_ME_METRIC_SPECS] + [
                "dapi_nuclei", "cfos_cells", "npy_cells", "pomc_cells", "cfos_only",
                "cfos_npy", "cfos_pomc", "cfos_npy_pomc", "cfos_pomc_all",
                "analyzed_area_um2", "pomc_assay_area_um2", "dapi_per_um2", "cfos_per_um2",
                "pomc_per_um2", "cfos_pomc_per_um2",
            ]
            numeric_cols = [c for c in candidates if c in animal_df.columns]
            group_cols = [c for c in ["cond", "genotype", "sex"] if c in animal_df.columns]
            if group_cols and numeric_cols:
                animal_df.groupby(group_cols, dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_group_summary_cond_genotype_sex.csv")
            group_cols2 = [c for c in ["cond", "genotype"] if c in animal_df.columns]
            if group_cols2 and numeric_cols:
                animal_df.groupby(group_cols2, dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_group_summary_cond_genotype.csv")
            if "cage" in animal_df.columns:
                animal_df.groupby(["cage"], dropna=False)[numeric_cols].agg(["mean", "std", "count"]).to_csv(outdir / "per_cage_summary.csv")
    except Exception:
        logger.exception("Failed group summaries for reconstructed analysis")


def legacy_main() -> None:
    ap = build_argparser()
    args = ap.parse_args()

    # Keep the full run as the default behavior. --full-run is a readable no-op flag.
    if args.full_run:
        args.migrate_only = False
        args.segment_only = False
        args.quantify_only = False

    # Optional conda/bootstrap happens before scientific imports.
    bootstrap_conda_if_needed(args)
    import_science_stack()

    raw_root = Path(args.raw_root).resolve()
    sanitized_root = Path(args.sanitized_root).resolve()
    args.sanitized_root = str(sanitized_root)
    outdir = Path(args.outdir).resolve() if args.outdir else sanitized_root / "analysis_arc_me_definitive_stats_v19"
    ensure_dir(outdir)
    logger = setup_logger(outdir, args.verbose)

    logger.warning("=== Keyence 20x reconstructed-section / projected-ventricle-wall ARC/ME pipeline starting ===")
    logger.warning("Sanitized/manual-seg input: %s", sanitized_root)
    logger.warning("Raw input, only if migration is needed: %s", raw_root)
    logger.warning("Output dir: %s", outdir)
    logger.warning("Existing *_seg.npy files are preserved unless --overwrite-seg is used.")
    logger.warning("No left/right statistical analysis will be produced; left/right halves are reconstructed before quantification.")
    logger.warning("Repetition folders such as FR6-2(R) are merged into the same animal N: %s", bool(args.keep_repetitions_as_same_animal))
    logger.warning("Region mode: %s | ARC/ME-only: %s | outside action: %s | n_bins=%d", args.region_mode, bool(args.arc_me_only), args.outside_arc_me_action, int(args.n_bins))
    logger.warning("ARC/ME rule: tissue/nuclei enclosed by projected left/right third-ventricle walls = ME; every other nucleus = ARC.")
    logger.warning("Anatomical wall fitting does not use c-FOS, NPY, or POMC signal; marker signals are quantified only after region masks are fixed.")

    samplesheet_path = Path(args.samplesheet).resolve() if args.samplesheet else None
    samplesheet = load_samplesheet(samplesheet_path, outdir, logger)

    # 1) Migrate only if needed/requested.
    migration_rows: List[MigrationRow] = []
    if (not sanitized_root.exists()) or args.force_migrate:
        migration_rows = migrate_raw_to_sanitized(raw_root, sanitized_root, args, logger)
        mig_df = pd.DataFrame([asdict(r) for r in migration_rows])
        mig_df.to_csv(outdir / "migration_report.csv", index=False)
        logger.warning("Migration rows: %d", len(migration_rows))
    else:
        logger.warning("Using existing sanitized folder; migration skipped.")

    if args.migrate_only:
        logger.warning("--migrate-only requested; stopping after migration.")
        return

    planes = filter_planes_by_args(discover_planes(sanitized_root, logger), args)
    disc_df = plane_discovery_dataframe(planes)
    if not disc_df.empty:
        disc_df = attach_sample_metadata(disc_df, samplesheet)
    disc_df.to_csv(outdir / "discovered_planes.csv", index=False)
    logger.warning("Discovered side plane directories: %d", len(planes))

    # 2) Segment only missing masks, while preserving manual masks.
    if args.segment_missing and not args.quantify_only:
        seg_report = segment_missing_masks(planes, args, logger, outdir)
        logger.warning("Segmentation report rows: %d", seg_report.shape[0])
        # Refresh planes after segmentation stage.
        planes = filter_planes_by_args(discover_planes(sanitized_root, logger), args)
    else:
        logger.warning("Skipping Cellpose segmentation stage.")

    if args.segment_only:
        logger.warning("--segment-only requested; stopping after segmentation.")
        return

    # 3) Reconstruct left/right section pairs, then quantify bins/regions.
    pairs = build_section_pairs(planes, samplesheet, args)
    pair_df = pd.DataFrame([{
        "sample": p.sample,
        "sample_key": p.sample_key,
        "animal_id": p.animal_id,
        "animal_key": p.animal_key,
        "is_repetition": p.is_repetition,
        "repetition_label": p.repetition_label,
        "section_index": p.section_index,
        "tile": p.tile,
        "has_left": p.left is not None,
        "has_right": p.right is not None,
        "left_dir": str(p.left.plane_dir) if p.left else "",
        "right_dir": str(p.right.plane_dir) if p.right else "",
    } for p in pairs])
    if not pair_df.empty:
        pair_df = attach_metadata_by_animal(pair_df, samplesheet)
    pair_df.to_csv(outdir / "discovered_reconstructed_sections.csv", index=False)
    logger.warning("Reconstructed section pairs: %d", len(pairs))

    cells_df, section_df, animal_df, animal_region_df, animal_bin_df = run_reconstructed_region_quantification(pairs, args, logger, outdir, samplesheet)
    logger.warning("Per-cell rows: %d", cells_df.shape[0] if cells_df is not None else 0)
    logger.warning("Per-reconstructed-section rows: %d", section_df.shape[0] if section_df is not None else 0)
    logger.warning("Per-animal rows, unique biological N: %d", animal_df.shape[0] if animal_df is not None else 0)

    # 4) Statistics and IQR outlier report on unique animals, total and region level.
    write_reconstructed_outlier_and_stats(animal_df, animal_region_df, outdir, args, logger)

    # 5) Base fallback plots use ME/ARC; the strict Figure 5 HIL path reaggregates ARC/ME/VMN.
    if args.make_plots:
        make_reconstructed_plots(animal_df, animal_region_df, outdir, args, logger)
        logger.warning("Final ARC/ME plots written to: %s", outdir / "plots")

    logger.warning("=== DONE ===")
    logger.warning("Main outputs:")
    logger.warning("  %s", outdir / "per_animal_total_summary.csv")
    logger.warning("  %s", outdir / "per_animal_region_summary.csv")
    logger.warning("  %s", outdir / "per_animal_bin_summary.csv")
    logger.warning("  %s", outdir / "per_reconstructed_section_summary.csv")
    logger.warning("  %s", outdir / "per_cell_reconstructed_measurements.csv")
    logger.warning("  %s", outdir / "bin_region_map.csv")
    if args.save_qc:
        logger.warning("  %s", outdir / args.reconstructed_qc_dirname)
    if args.save_reconstructed_tifs:
        logger.warning("  %s", outdir / args.reconstructed_dirname)



# =============================================================================
# Keyence 10x DAPI/GFAP/Iba1 ARC/ME + microglial morphology pipeline
# =============================================================================

import json
import zipfile
from collections import defaultdict

warnings.filterwarnings("ignore", category=FutureWarning, message=r".*deprecated.*")

GFAP_IBA1_VERSION = "1.2-reviewed-dapi-iba1-animal-blocked"
GFAP_IBA1_DEFAULT_RAW = DEFAULT_RAW_ROOT
GFAP_IBA1_DEFAULT_WORK = DEFAULT_SANITIZED_ROOT
GFAP_IBA1_DEFAULT_UM_PER_PX = 1.51
GFAP_IBA1_EXPECTED_REVIEW_SECTIONS = 61
GFAP_IBA1_ENV = "paper_apotome_repro"
GFAP_IBA1_MARKERS = ("DAPI", "GFAP", "Iba1")
MICROGLIA_STATES = ("Ramified", "Rod-like", "Activated", "Amoeboid")
MICROGLIA_STATE_TO_CODE = {s: i for i, s in enumerate(MICROGLIA_STATES)}
MICROGLIA_STATE_COLORS = {
    "Ramified": "#2c7fb8",
    "Rod-like": "#31a354",
    "Activated": "#f39c12",
    "Amoeboid": "#d7301f",
    "Uncertain": "#777777",
}
RE_GI_SECTION_DIR = re.compile(r"^S(?P<num>\d{1,3})$", re.IGNORECASE)
RE_GI_XY_DIR = re.compile(r"^XY(?P<num>\d{1,3})$", re.IGNORECASE)


# Override the old setup package checks. These definitions are intentionally
# later in the file, so the existing safe conda bootstrap calls these versions.
def _env_has_required_packages(conda_exe: str, env_name: str, need_cellpose: bool, need_region_dl: bool = False) -> bool:
    modules = [
        "numpy", "pandas", "scipy", "skimage", "matplotlib", "tifffile",
        "PIL", "imagecodecs", "sklearn", "torch",
    ]
    if need_cellpose:
        modules.append("cellpose")
    code = "\n".join([f"import {m}" for m in modules])
    proc = _conda_run_python(conda_exe, env_name, code)
    return proc.returncode == 0


def _install_pipeline_packages(conda_exe: str, env_name: str, *, need_cellpose: bool, need_region_dl: bool = False) -> None:
    core = [
        "numpy", "pandas", "scipy", "scikit-image", "matplotlib", "tifffile",
        "imagecodecs", "pillow", "scikit-learn", "torch",
    ]
    pkgs = core + (["cellpose"] if need_cellpose else [])
    subprocess.check_call([conda_exe, "run", "-n", env_name, "python", "-m", "pip", "install", "--upgrade", "pip"])
    subprocess.check_call([conda_exe, "run", "-n", env_name, "python", "-m", "pip", "install", "--upgrade"] + pkgs)


def build_argparser_gfap_iba1() -> argparse.ArgumentParser:
    """Extend the validated v19 parser while preserving every ARC/ME flag."""
    ap = build_argparser()
    ap.add_argument("--force", action="store_true", help="Safely clean only the designated generated analysis output before running; protected inputs and HIL annotations are never removed.")
    ap.description = (
        "Keyence 10x CH1=DAPI, CH3=GFAP, CH4=Iba1 full-section ARC/ME pipeline. "
        "It preserves the validated v17 projected-ventricle-wall anatomy, quantifies "
        "regional GFAP area/intensity, segments DAPI-associated Iba1 microglia, and "
        "classifies each microglial cell as Ramified, Rod-like, Activated, or Amoeboid."
    )
    ap.set_defaults(
        raw_root=str(GFAP_IBA1_DEFAULT_RAW),
        sanitized_root=str(GFAP_IBA1_DEFAULT_WORK),
        env_name=GFAP_IBA1_ENV,
        um_per_px=GFAP_IBA1_DEFAULT_UM_PER_PX,
        region_mode="arc-me-walls",
        outside_arc_me_action="arc",
        arc_me_only=True,
        reconstructed_dirname="processed_sections",
        reconstructed_qc_dirname="qc_arc_me_gfap_iba1",
        qc_mask_dirname="qc_arc_me_gfap_iba1",
        segment_missing=False,
        plot_suite="final",
        region_tissue_close_um=18.0,
    )

    src = ap.add_argument_group("Keyence 10x DAPI/GFAP/Iba1 source")
    src.add_argument(
        "--source-zip", default=None,
        help="Optional .zip containing Iba1_GFAP_final_10x. Overrides --raw-input when supplied.",
    )
    src.add_argument(
        "--full-section-xy", action="store_true", default=True,
        help="Treat each XY## folder as one complete bilateral section (the correct layout for this 10x dataset).",
    )
    src.add_argument(
        "--no-full-section-xy", dest="full_section_xy", action="store_false",
        help="Reserved compatibility switch; not recommended for this dataset.",
    )

    manual_regions = ap.add_argument_group("Manual ARC/ME region masks")
    manual_regions.add_argument(
        "--human-in-the-loop", "--human_in_the_loop",
        dest="human_in_the_loop", action="store_true", default=True,
        help=("Compatibility flag; every run requires one HIL-accepted ARC/ME/VMN TIFF+JSON review pair for every "
              "discovered image before quantification. If --manual-region-dir is "
              "omitted, use <Paper>/analyses/Fig5/hil/accepted_annotations. "
              "This mode never falls back to automated "
              "ARC/ME inference."),
    )
    manual_regions.add_argument("--human-review-host", default="127.0.0.1")
    manual_regions.add_argument(
        "--human-review-port", type=int, default=0,
        help=("Loopback reviewer port. The default 0 asks the OS for a free port "
              "and avoids colliding with an already-running review session."),
    )
    manual_regions.add_argument(
        "--manual-region-dir", default=None,
        help=("Optional root containing reviewed ARC/ME/VMN TIFF+JSON pairs. Exact layout: "
              "<root>/<animal>/<section>/<animal>_<section>_ARC_ME_labels.{tif,json}. "
              "The TIFF uses 0=background, 1=ARC, 2=ME, 3=VMN; the JSON must record human "
              "acceptance, reviewer, timestamp, dimensions, codes, and TIFF SHA-256."),
    )
    manual_regions.add_argument(
        "--require-manual-regions", action="store_true", default=True,
        help=("Compatibility flag; every analysis fails if --manual-region-dir is omitted or any "
              "discovered section has a missing or invalid manual ARC/ME mask. "
              "Automated ARC/ME fallback is disabled."),
    )
    manual_regions.add_argument(
        "--make-figure", action="store_true", default=False,
        help=("Render Figure 5 from this run once the analysis finishes. The interactive "
              "review prompt offers the same action as make_figure_5."),
    )
    manual_regions.add_argument(
        "--figure-output-root", default=str(FIG5_ANALYSIS_ROOT / "figure"),
        help=("Parent of the fresh, absent directory --make-figure renders Figure 5 into. "
              "It defaults inside the tree so the render can be found afterwards."),
    )

    seg = ap.add_argument_group("DAPI and Iba1 segmentation")
    seg.add_argument(
        "--reviewed-dapi-root", default=str(DEFAULT_REVIEWED_DAPI_ROOT),
        help=("Versioned reviewed DAPI masks used for all quantification. The root must "
              "contain one PASS receipt and segmentation for every discovered section."),
    )
    seg.add_argument(
        "--reviewed-iba1-root", default=str(DEFAULT_REVIEWED_IBA1_ROOT),
        help=("Versioned DAPI-seeded Iba1 masks used for all quantification. The root must "
              "contain one PASS receipt and segmentation for every discovered section."),
    )
    seg.add_argument(
        "--allow-historical-cell-masks", action="store_true", default=False,
        help=("Explicit maintenance escape hatch: use the immutable historical raw-local "
              "DAPI/Iba1 sidecars instead of the reviewed versioned masks."),
    )
    seg.add_argument(
        "--dapi-segmentation-mode", choices=["auto", "cellpose", "watershed"], default="auto",
        help="auto tries the supplied DAPI Cellpose model then falls back to marker-controlled watershed.",
    )
    seg.add_argument(
        "--iba1-segmentation-mode", choices=["seeded-watershed", "existing", "cellpose"], default="seeded-watershed",
        help="Instance segmentation strategy. Seeded watershed uses Iba1-associated DAPI nuclei and preserves processes.",
    )
    seg.add_argument("--model-iba1", default=None, help="Optional trained Cellpose Iba1 model used with --iba1-segmentation-mode cellpose.")
    seg.add_argument("--diameter-iba1", type=float, default=None, help="Optional Cellpose Iba1 diameter in pixels.")
    seg.add_argument("--dapi-min-area-px", type=int, default=5)
    seg.add_argument("--dapi-max-area-px", type=int, default=350)
    seg.add_argument("--dapi-watershed-min-distance-px", type=int, default=2)
    seg.add_argument("--iba1-background-radius-um", type=float, default=22.0)
    seg.add_argument("--iba1-threshold-factor", type=float, default=0.65)
    seg.add_argument("--iba1-nucleus-expand-um", type=float, default=7.5)
    seg.add_argument("--iba1-nucleus-min-score-quantile", type=float, default=0.60)
    seg.add_argument("--iba1-max-cell-radius-um", type=float, default=55.0)
    seg.add_argument("--iba1-min-cell-area-um2", type=float, default=12.0)
    seg.add_argument("--iba1-max-cell-area-um2", type=float, default=3500.0)
    seg.add_argument("--iba1-min-foreground-fraction", type=float, default=0.02)
    seg.add_argument(
        "--microglia-qc-max-nucleus-distance-um", type=float, default=4.5,
        help="Maximum Iba1-mask to DAPI-nucleus distance for biological-cell acceptance.",
    )
    seg.add_argument(
        "--microglia-qc-edge-margin-px", type=int, default=1,
        help="Reject truncated Iba1 objects touching this many pixels at an image edge.",
    )
    seg.add_argument(
        "--microglia-qc-min-anatomy-fraction", type=float, default=0.80,
        help="Minimum object fraction inside the ARC/ME anatomical tissue support.",
    )
    seg.add_argument(
        "--microglia-qc-min-largest-component-fraction", type=float, default=0.85,
        help="Minimum fraction belonging to the largest connected component of one instance label.",
    )
    seg.add_argument(
        "--microglia-qc-max-major-axis-um", type=float, default=130.0,
        help="Reject oversized Iba1 structures longer than this calibrated major-axis length.",
    )
    seg.add_argument(
        "--microglia-qc-linear-min-major-axis-um", type=float, default=60.0,
        help="Minimum length at which the compound line/vessel rule is evaluated.",
    )
    seg.add_argument(
        "--microglia-qc-linear-min-aspect-ratio", type=float, default=6.0,
        help="Aspect-ratio component of the compound line/vessel rule.",
    )
    seg.add_argument(
        "--microglia-qc-max-dapi-nuclei-touched", type=int, default=4,
        help="Maximum nearby DAPI nuclei before a long/large object is considered a non-cellular structure.",
    )
    seg.add_argument(
        "--accept-empty-iba1-seg", action="store_true", default=False,
        help=("Treat an existing zero-object Iba1 mask as a valid biological zero. "
              "Default: mark it invalid; regenerate it only when --segment-missing is enabled, "
              "and exclude its area from microglial density denominators otherwise."),
    )

    gf = ap.add_argument_group("GFAP and Iba1 regional intensity")
    gf.add_argument("--gfap-background-radius-um", type=float, default=35.0)
    gf.add_argument("--gfap-threshold-factor", type=float, default=0.80)
    gf.add_argument("--gfap-robust-k", type=float, default=3.0)
    gf.add_argument("--gfap-min-object-area-um2", type=float, default=8.0)
    gf.add_argument("--intensity-percentile-low", type=float, default=0.5)
    gf.add_argument("--intensity-percentile-high", type=float, default=99.8)

    clf = ap.add_argument_group("Microglial morphology deep learning")
    clf.add_argument(
        "--microglia-classifier-mode", choices=["hybrid", "cnn", "morphometry", "supervised"], default="hybrid",
        help=("hybrid: morphometric clustering creates high-confidence pseudo-labels, then a CNN learns image morphology; "
              "supervised requires --microglia-labels-csv; morphometry skips CNN."),
    )
    clf.add_argument(
        "--microglia-labels-csv", default=None,
        help="Optional reviewed labels. Accepted columns: cell_uid,state or patch_path,state. Manual labels override pseudo-labels.",
    )
    clf.add_argument("--cnn-device", choices=["auto", "cpu", "cuda"], default="auto")
    clf.add_argument("--cnn-epochs", type=int, default=25)
    clf.add_argument("--cnn-batch-size", type=int, default=64)
    clf.add_argument("--cnn-learning-rate", type=float, default=1e-3)
    clf.add_argument("--cnn-weight-decay", type=float, default=1e-4)
    clf.add_argument("--cnn-patience", type=int, default=6)
    clf.add_argument(
        "--cnn-calibrate-batchnorm", action="store_true",
        help="Recalculate batch-normalization statistics from unaugmented training cells before each validation pass; useful for small reviewed datasets.",
    )
    clf.add_argument("--cnn-patch-size", type=int, default=96)
    clf.add_argument("--cnn-min-training-cells", type=int, default=120)
    clf.add_argument("--cnn-min-class-cells", type=int, default=12)
    clf.add_argument(
        "--cnn-min-training-animals", type=int, default=3,
        help=("Minimum number of animals contributing CNN training labels. Validation "
              "is always animal-blocked; the pipeline retains morphometric classes "
              "instead of falling back to a leaking random cell split."),
    )
    clf.add_argument("--cnn-pseudo-confidence", type=float, default=0.12)
    clf.add_argument("--cnn-prediction-confidence", type=float, default=0.45)
    clf.add_argument("--cnn-seed", type=int, default=20260620)
    clf.add_argument("--cnn-model-out", default=None, help="Optional model output path. Default: <output>/models/microglia_morphology_cnn.pt")
    clf.add_argument("--cnn-model-in", default=None, help="Optional previously trained model to use without retraining.")
    patch_group = clf.add_mutually_exclusive_group()
    patch_group.add_argument("--save-microglia-patches", dest="save_microglia_patches", action="store_true", default=True)
    patch_group.add_argument("--no-microglia-patches", dest="save_microglia_patches", action="store_false")
    clf.add_argument("--prepare-labels-only", action="store_true", help="Segment/extract cells and write annotation template, but do not train or run statistics.")
    clf.add_argument(
        "--microglia-review-port", type=int, default=0,
        help=("Port for the microglia morphology reviewer, the second Figure 5 HIL "
              "server. It opens beside the ARC/ME/VMN region reviewer. 0 picks a free port."),
    )
    clf.add_argument(
        "--review-microglia", action="store_true",
        help=("Open the microglia morphology reviewer on this run's analysis directory "
              "once the analysis has written the cell patches, and wait there. Use it to "
              "turn model proposals into the reviewed labels --microglia-labels-csv wants."),
    )

    stats = ap.add_argument_group("GFAP/Iba1 definitive endpoints")
    stats.add_argument("--microglia-density-scale-um2", type=float, default=100000.0)
    stats.add_argument(
        "--activated-states", default="Activated,Amoeboid",
        help="Comma-separated microglial classes considered activated in the pooled activated fraction.",
    )
    stats.add_argument(
        "--resting-states", default="Ramified",
        help="Comma-separated classes considered resting. Rod-like remains separately reported by default.",
    )
    stats.add_argument("--gfap-intensity-scale", type=float, default=1.0)
    return ap


@dataclass
class GISection:
    sample: str
    sample_dir: Path
    section_index: int
    section: str
    section_dir: Path
    images: Dict[str, Optional[Path]]
    segs: Dict[str, Optional[Path]]


def _gi_signal_from_rgb(arr: Any, marker: str) -> Any:
    """Extract scalar fluorescence from Keyence pseudo-coloured RGB TIFFs."""
    x = np.asarray(arr)
    while x.ndim > 3:
        x = x.max(axis=0)
    if x.ndim == 3 and x.shape[-1] not in (3, 4):
        x = x.max(axis=0)
    if x.ndim == 2:
        return x.astype(np.float32)
    if x.ndim != 3:
        return np.squeeze(x).astype(np.float32)
    rgb = x[..., :3].astype(np.float32)
    m = str(marker).lower()
    if m == "dapi":
        return rgb[..., 2]
    if m == "gfap":
        return rgb[..., 0]
    if m == "iba1":
        # CH4 is magenta in this Keyence export; retain signal common to R and B.
        return 0.5 * (rgb[..., 0] + rgb[..., 2])
    return np.max(rgb, axis=-1)


def _gi_read_signal(path: Path, marker: str) -> Any:
    return _gi_signal_from_rgb(read_image(path), marker)


def _gi_robust_normalize(x: Any, args: argparse.Namespace) -> Any:
    a = np.asarray(x, dtype=np.float32)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a, dtype=np.float32)
    lo = float(np.percentile(finite, float(args.intensity_percentile_low)))
    hi = float(np.percentile(finite, float(args.intensity_percentile_high)))
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((a - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def _gi_save_seg(path: Path, labels_img: Any, *, source: str, metadata: Optional[Dict[str, Any]] = None) -> None:
    ensure_dir(path.parent)
    payload: Dict[str, Any] = {
        "masks": np.asarray(labels_img, dtype=np.int32),
        "source": source,
        "version": GFAP_IBA1_VERSION,
    }
    if metadata:
        payload.update(metadata)
    np.save(path, payload, allow_pickle=True)


def _gi_find_source_root(path: Path) -> Path:
    """Return the folder whose direct children are animal folders."""
    p = path.resolve()
    if not p.exists():
        raise FileNotFoundError(p)
    dirs = [d for d in p.iterdir() if d.is_dir() and not d.name.startswith(".")]
    if any(any(RE_GI_XY_DIR.match(c.name) for c in d.iterdir() if c.is_dir()) for d in dirs):
        return p
    if len(dirs) == 1:
        d = dirs[0]
        children = [x for x in d.iterdir() if x.is_dir()]
        if any(any(RE_GI_XY_DIR.match(c.name) for c in x.iterdir() if c.is_dir()) for x in children):
            return d
    # Search one additional level for the named extracted folder.
    hits = [d for d in p.rglob("Iba1_GFAP_final_10x") if d.is_dir()]
    if hits:
        return hits[0]
    raise RuntimeError(f"Could not find animal/XY## structure under {p}")


def _gi_prepare_raw_source(args: argparse.Namespace, work_root: Path, logger: logging.Logger) -> Path:
    source = Path(args.source_zip).resolve() if args.source_zip else Path(args.raw_root).resolve()
    if source.suffix.lower() == ".zip":
        extract_dir = work_root.parent / f".{work_root.name}_source_extracted"
        sentinel = extract_dir / ".extract_complete"
        if bool(args.force_migrate) and extract_dir.exists():
            shutil.rmtree(extract_dir)
        if not sentinel.exists():
            logger.warning("Extracting source ZIP: %s", source)
            ensure_dir(extract_dir)
            with zipfile.ZipFile(source, "r") as zf:
                zf.extractall(extract_dir)
            sentinel.write_text(f"source={source}\n", encoding="utf-8")
        else:
            logger.info("Using existing extracted ZIP cache: %s", extract_dir)
        return _gi_find_source_root(extract_dir)
    return _gi_find_source_root(source)


def _gi_copy_if_needed(src: Path, dst: Path, args: argparse.Namespace, logger: logging.Logger) -> str:
    return safe_copy(src, dst, overwrite=bool(args.overwrite_sanitized_files), dry_run=bool(args.dry_run), logger=logger)


def migrate_gfap_iba1_to_sanitized(raw_root: Path, work_root: Path, args: argparse.Namespace, logger: logging.Logger) -> Any:
    """Map each full bilateral XY image to one S## section. No left/right pairing is performed."""
    rows: List[Dict[str, Any]] = []
    if work_root.exists() and not bool(args.force_migrate):
        existing = list(work_root.glob("*/S*/Image_*_DAPI.tif")) + list(work_root.glob("*/S*/Image_*_DAPI.tiff"))
        if existing:
            logger.warning("Using existing sanitized 10x folder; source migration skipped: %s", work_root)
            return pd.DataFrame(rows)
    ensure_dir(work_root)
    for sample_dir in sorted([p for p in raw_root.iterdir() if p.is_dir() and not p.name.startswith(".")]):
        xy_dirs = [p for p in sample_dir.iterdir() if p.is_dir() and RE_GI_XY_DIR.match(p.name)]
        for xy_dir in sorted(xy_dirs, key=lambda p: int(RE_GI_XY_DIR.match(p.name).group("num"))):
            idx = int(RE_GI_XY_DIR.match(xy_dir.name).group("num"))
            sec = f"S{idx:02d}"
            dst_dir = work_root / sample_dir.name / sec
            mapping = {"CH1": "DAPI", "CH3": "GFAP", "CH4": "Iba1"}
            for ch, marker in mapping.items():
                hits = sorted(xy_dir.glob(f"*_{ch}.tif")) + sorted(xy_dir.glob(f"*_{ch}.tiff"))
                if not hits:
                    continue
                src = hits[0]
                dst = dst_dir / f"Image_{sec}_{marker}{src.suffix.lower()}"
                action = _gi_copy_if_needed(src, dst, args, logger)
                rows.append({"sample": sample_dir.name, "source_section": xy_dir.name, "section": sec, "marker": marker, "source": str(src), "destination": str(dst), "action": action})
                # Preserve channel- or marker-named manual segmentation if present.
                seg_hits = sorted(xy_dir.glob(f"*_{ch}_seg.npy")) + sorted(xy_dir.glob(f"*_{marker}_seg.npy"))
                if seg_hits:
                    sdst = dst_dir / f"Image_{sec}_{marker}_seg.npy"
                    saction = _gi_copy_if_needed(seg_hits[0], sdst, args, logger)
                    rows.append({"sample": sample_dir.name, "source_section": xy_dir.name, "section": sec, "marker": marker + "_seg", "source": str(seg_hits[0]), "destination": str(sdst), "action": saction})
            ov = sorted(xy_dir.glob("*Overlay.tif")) + sorted(xy_dir.glob("*Overlay.tiff"))
            if ov:
                dst = dst_dir / f"Image_{sec}_Overlay{ov[0].suffix.lower()}"
                action = _gi_copy_if_needed(ov[0], dst, args, logger)
                rows.append({"sample": sample_dir.name, "source_section": xy_dir.name, "section": sec, "marker": "Overlay", "source": str(ov[0]), "destination": str(dst), "action": action})
    return pd.DataFrame(rows)


def _gi_find_file(folder: Path, marker: str, seg: bool = False) -> Optional[Path]:
    ml = marker.lower()
    patterns = [f"*_{ml}_seg.npy"] if seg else [f"*_{ml}.tif", f"*_{ml}.tiff"]
    hits: List[Path] = []
    for pat in patterns:
        hits.extend(folder.glob(pat))
        hits.extend(folder.glob(pat.replace(ml, marker)))
    return sorted(set(hits))[0] if hits else None


def discover_gfap_iba1_sections(work_root: Path, args: argparse.Namespace, logger: logging.Logger) -> List[GISection]:
    sections: List[GISection] = []
    if not work_root.exists():
        return sections
    include = {canonical_sample_id(x) for x in (args.include_samples or [])}
    exclude = {canonical_sample_id(x) for x in (args.exclude_samples or [])}
    for sample_dir in sorted([p for p in work_root.iterdir() if p.is_dir() and not p.name.startswith(".")]):
        sk = canonical_sample_id(sample_dir.name)
        if include and sk not in include:
            continue
        if sk in exclude:
            continue
        for sec_dir in sorted([p for p in sample_dir.iterdir() if p.is_dir() and RE_GI_SECTION_DIR.match(p.name)], key=lambda p: int(RE_GI_SECTION_DIR.match(p.name).group("num"))):
            idx = int(RE_GI_SECTION_DIR.match(sec_dir.name).group("num"))
            images = {m: _gi_find_file(sec_dir, m, False) for m in GFAP_IBA1_MARKERS}
            segs = {m: _gi_find_file(sec_dir, m, True) for m in ("DAPI", "Iba1")}
            sections.append(GISection(sample_dir.name, sample_dir, idx, sec_dir.name, sec_dir, images, segs))
    if int(args.limit_planes or 0) > 0:
        sections = sections[:int(args.limit_planes)]
    logger.warning("Discovered full bilateral 10x sections: %d", len(sections))
    return sections


def _gi_require_raw_local_masks(sections: Sequence[GISection], work_root: Path) -> None:
    """Require each DAPI/Iba1 image to use its canonical sibling raw mask."""
    root = work_root.resolve()
    problems: List[str] = []
    for section in sections:
        for marker in ("DAPI", "Iba1"):
            image = section.images.get(marker)
            if image is None:
                problems.append(f"{section.sample}/{section.section}/{marker}: missing raw image")
                continue
            expected = image.with_name(f"{image.stem}_seg.npy").resolve()
            mask = section.segs.get(marker)
            if mask is None:
                problems.append(f"{section.sample}/{section.section}/{marker}: missing {expected.name}")
                continue
            mask_path = Path(mask).expanduser()
            resolved = mask_path.resolve()
            try:
                resolved.relative_to(root)
            except ValueError:
                problems.append(f"{section.sample}/{section.section}/{marker}: mask is outside raw root")
                continue
            if resolved != expected:
                problems.append(f"{section.sample}/{section.section}/{marker}: noncanonical mask {mask_path.name}")
            elif mask_path.is_symlink():
                problems.append(f"{section.sample}/{section.section}/{marker}: mask is a symlink")
            elif not resolved.is_file() or resolved.stat().st_size <= 0:
                problems.append(f"{section.sample}/{section.section}/{marker}: mask is missing or empty")
    if problems:
        detail = "\n".join(f"  - {item}" for item in problems[:20])
        if len(problems) > 20:
            detail += f"\n  - ... and {len(problems) - 20} more"
        raise ManualRegionMaskError(
            f"Raw-local Cellpose coverage is incomplete ({len(problems)} errors):\n{detail}"
        )


def _gi_apply_reviewed_cell_masks(
    sections: Sequence[GISection],
    args: argparse.Namespace,
    outdir: Path,
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    """Validate and select the versioned DAPI/Iba1 masks without touching raw_data."""
    rows: List[Dict[str, Any]] = []
    problems: List[str] = []
    ensure_dir(outdir)

    if bool(getattr(args, "allow_historical_cell_masks", False)):
        for section in sections:
            for marker in ("DAPI", "Iba1"):
                image = section.images.get(marker)
                mask = image.with_name(f"{image.stem}_seg.npy") if image is not None else None
                section.segs[marker] = mask
                rows.append({
                    "sample": section.sample,
                    "section": section.section,
                    "marker": marker,
                    "status": "historical_raw_local_explicitly_allowed",
                    "mask_path": str(mask or ""),
                })
        pd.DataFrame(rows).to_csv(outdir / "reviewed_cell_mask_coverage.csv", index=False)
        logger.warning("Explicit override active: using historical raw-local DAPI/Iba1 sidecars.")
        return rows

    roots = {
        "DAPI": Path(args.reviewed_dapi_root).expanduser().resolve(),
        "Iba1": Path(args.reviewed_iba1_root).expanduser().resolve(),
    }
    for marker, root in roots.items():
        if not root.is_dir() or root.is_symlink():
            problems.append(f"{marker}: reviewed mask root is missing, not physical, or not a directory: {root}")

    if len(sections) != GFAP_IBA1_EXPECTED_REVIEW_SECTIONS:
        problems.append(
            f"dataset: expected {GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} sections, found {len(sections)}"
        )

    for section in sections:
        dapi_image = section.images.get("DAPI")
        iba1_image = section.images.get("Iba1")
        source_dapi_sha = sha256_file(dapi_image) if dapi_image and dapi_image.is_file() else None
        source_iba1_sha = sha256_file(iba1_image) if iba1_image and iba1_image.is_file() else None
        reviewed_dapi_mask = (
            roots["DAPI"] / section.sample / section.section
            / f"{section.sample}_{section.section}_DAPI_seg.npy"
        )
        for marker in ("DAPI", "Iba1"):
            root = roots[marker]
            prefix = f"{section.sample}_{section.section}_{marker}"
            section_root = root / section.sample / section.section
            mask_path = section_root / f"{prefix}_seg.npy"
            label_path = section_root / f"{prefix}_labels.tif"
            qc_path = section_root / f"{prefix}_QC.png"
            receipt_path = section_root / f"{prefix}_receipt.json"
            row: Dict[str, Any] = {
                "sample": section.sample,
                "section": section.section,
                "marker": marker,
                "status": "pending",
                "root": str(root),
                "mask_path": str(mask_path),
                "receipt_path": str(receipt_path),
            }
            section_errors: List[str] = []
            expected_outputs = {
                "segmentation_npy": mask_path,
                "label_tiff": label_path,
                "qc_overlay_png": qc_path,
            }
            for output_name, output_path in expected_outputs.items():
                if not output_path.is_file() or output_path.is_symlink():
                    section_errors.append(f"missing physical {output_name}: {output_path}")
            if not receipt_path.is_file() or receipt_path.is_symlink():
                section_errors.append(f"missing physical receipt: {receipt_path}")
                receipt: Dict[str, Any] = {}
            else:
                try:
                    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                except Exception as exc:
                    receipt = {}
                    section_errors.append(f"cannot parse receipt {receipt_path}: {exc}")
            if receipt:
                if receipt.get("status") != "PASS":
                    section_errors.append(f"receipt status is not PASS: {receipt.get('status')!r}")
                if receipt.get("raw_data_modified") is not False:
                    section_errors.append("receipt does not assert raw_data_modified=false")
                if receipt.get("source_dapi_sha256") != source_dapi_sha:
                    section_errors.append("source DAPI SHA-256 mismatch")
                if marker == "Iba1":
                    if receipt.get("source_iba1_sha256") != source_iba1_sha:
                        section_errors.append("source Iba1 SHA-256 mismatch")
                    if reviewed_dapi_mask.is_file():
                        reviewed_dapi_sha = sha256_file(reviewed_dapi_mask)
                        if receipt.get("reviewed_dapi_segmentation_sha256") != reviewed_dapi_sha:
                            section_errors.append("reviewed DAPI segmentation SHA-256 mismatch")
                    else:
                        section_errors.append(f"reviewed DAPI segmentation is missing: {reviewed_dapi_mask}")
                output_hashes = receipt.get("output_hashes")
                if not isinstance(output_hashes, dict):
                    section_errors.append("receipt output_hashes is missing")
                    output_hashes = {}
                for output_name, output_path in expected_outputs.items():
                    expected_hash = output_hashes.get(f"{output_name}_sha256")
                    if output_path.is_file():
                        actual_hash = sha256_file(output_path)
                        row[f"{output_name}_sha256"] = actual_hash
                        if expected_hash != actual_hash:
                            section_errors.append(f"{output_name} SHA-256 mismatch")
            if mask_path.is_file():
                try:
                    labels = ensure_label_image(load_seg_mask(mask_path))
                    image = section.images.get(marker)
                    if image is None or not image.is_file():
                        section_errors.append(f"missing source {marker} image")
                    else:
                        expected_shape = tuple(np.asarray(_gi_read_signal(image, marker)).shape)
                        if tuple(labels.shape) != expected_shape:
                            section_errors.append(
                                f"mask shape {tuple(labels.shape)} != source shape {expected_shape}"
                            )
                    if int(labels.max()) <= 0:
                        section_errors.append("mask contains no instances")
                    row["instances"] = int(labels.max())
                    row["mask_shape"] = "x".join(str(value) for value in labels.shape)
                except Exception as exc:
                    section_errors.append(f"cannot load segmentation: {exc}")
            if section_errors:
                row["reason"] = " | ".join(section_errors)
                problems.append(
                    f"{section.sample}/{section.section}/{marker}: {row['reason']}"
                )
            else:
                row["status"] = "accepted"
                row["reason"] = ""
                section.segs[marker] = mask_path.resolve()
            rows.append(row)

    pd.DataFrame(rows).to_csv(outdir / "reviewed_cell_mask_coverage.csv", index=False)
    accepted = sum(row.get("status") == "accepted" for row in rows)
    logger.warning(
        "Versioned reviewed DAPI/Iba1 coverage: accepted=%d/%d masks",
        accepted,
        2 * GFAP_IBA1_EXPECTED_REVIEW_SECTIONS,
    )
    if problems:
        detail = "\n".join(f"  - {item}" for item in problems[:20])
        if len(problems) > 20:
            detail += f"\n  - ... and {len(problems) - 20} more"
        raise ManualRegionMaskError(
            f"Reviewed DAPI/Iba1 coverage is incomplete ({len(problems)} errors):\n{detail}"
        )
    return rows


def _gi_section_table(sections: List[GISection], samplesheet: Any) -> Any:
    rows = []
    for s in sections:
        base, rep, lab = strip_repetition_suffix(s.sample)
        rows.append({
            "sample": s.sample, "animal_id": base, "animal_key": canonical_sample_id(base),
            "is_repetition": rep, "repetition_label": lab, "section_index": s.section_index,
            "section": s.section, "section_dir": str(s.section_dir),
            "dapi_image": str(s.images.get("DAPI") or ""), "gfap_image": str(s.images.get("GFAP") or ""),
            "iba1_image": str(s.images.get("Iba1") or ""), "dapi_seg": str(s.segs.get("DAPI") or ""),
            "iba1_seg": str(s.segs.get("Iba1") or ""),
        })
    df = pd.DataFrame(rows)
    return attach_metadata_by_animal(df, samplesheet) if not df.empty else df


def _gi_watershed_dapi(signal: Any, args: argparse.Namespace) -> Any:
    from scipy import ndimage as ndi
    from skimage.feature import peak_local_max
    from skimage.filters import gaussian, threshold_otsu
    from skimage.segmentation import watershed
    from skimage.morphology import remove_small_objects as _rso, remove_small_holes as _rsh
    x = _gi_robust_normalize(signal, args)
    sm = gaussian(x, sigma=0.8, preserve_range=True)
    vals = sm[np.isfinite(sm) & (sm > 0)]
    thr = float(threshold_otsu(vals)) if vals.size > 10 else 0.10
    mask = sm > max(thr, 0.035)
    mask = _rso(mask, min_size=max(2, int(args.dapi_min_area_px)))
    mask = _rsh(mask, area_threshold=max(2, int(args.dapi_min_area_px)))
    dist = ndi.distance_transform_edt(mask)
    coords = peak_local_max(
        dist, labels=mask, min_distance=max(1, int(args.dapi_watershed_min_distance_px)),
        threshold_abs=0.8, exclude_border=False,
    )
    markers = np.zeros(mask.shape, dtype=np.int32)
    if coords.size:
        markers[tuple(coords.T)] = np.arange(1, coords.shape[0] + 1, dtype=np.int32)
    markers = label(markers > 0, connectivity=1)
    if int(markers.max()) == 0:
        markers = label(mask, connectivity=1)
    lbl = watershed(-dist, markers=markers, mask=mask)
    return filter_rois_by_area(lbl, int(args.dapi_min_area_px), int(args.dapi_max_area_px or 0)).astype(np.int32)


def _gi_resolve_torch_device(requested: str, logger: Optional[logging.Logger] = None) -> str:
    req = str(requested or "auto").lower()
    if req == "cpu":
        return "cpu"
    try:
        import torch
        if req in {"auto", "cuda"} and torch.cuda.is_available():
            try:
                _ = torch.zeros(1, device="cuda")
                return "cuda"
            except Exception as e:
                if logger:
                    logger.warning("CUDA was reported but could not initialize; using CPU: %s", e)
        if req == "cuda" and logger:
            logger.warning("CUDA explicitly requested but unavailable; using CPU.")
    except Exception:
        pass
    return "cpu"


def _gi_cellpose_dapi(signal: Any, args: argparse.Namespace, logger: logging.Logger) -> Optional[Any]:
    try:
        from cellpose import models
        use_gpu = bool(args.gpu)
        model_path = Path(str(args.model_dapi)).expanduser() if args.model_dapi else None
        if model_path and model_path.exists():
            model = models.CellposeModel(gpu=use_gpu, pretrained_model=str(model_path))
        else:
            try:
                model = models.CellposeModel(gpu=use_gpu, model_type="nuclei")
            except Exception:
                model = models.CellposeModel(gpu=use_gpu)
        img = _gi_robust_normalize(signal, args)
        kwargs: Dict[str, Any] = {"channels": [0, 0]}
        if args.diameter_dapi is not None:
            kwargs["diameter"] = float(args.diameter_dapi)
        if args.flow_threshold is not None:
            kwargs["flow_threshold"] = float(args.flow_threshold)
        if args.cellprob_threshold is not None:
            kwargs["cellprob_threshold"] = float(args.cellprob_threshold)
        out = model.eval(img, **kwargs)
        masks = out[0] if isinstance(out, tuple) else out
        return filter_rois_by_area(ensure_label_image(masks), int(args.dapi_min_area_px), int(args.dapi_max_area_px or 0))
    except Exception as e:
        logger.warning("DAPI Cellpose failed; watershed fallback will be used: %s", e)
        return None


def _gi_background_correct(signal: Any, radius_px: int) -> Tuple[Any, Any]:
    from scipy import ndimage as ndi
    x = np.asarray(signal, dtype=np.float32)
    size = int(max(3, 2 * int(radius_px) + 1))
    # Cap the footprint to keep the 960x720 10x images fast and memory-safe.
    size = min(size, 151)
    bg = ndi.grey_opening(x, size=(size, size))
    return np.clip(x - bg, 0, None).astype(np.float32), bg.astype(np.float32)


def _gi_robust_threshold(values: Any, factor: float = 1.0, k: float = 3.0) -> float:
    from skimage.filters import threshold_otsu
    v = np.asarray(values, dtype=np.float32)
    v = v[np.isfinite(v)]
    if v.size < 10:
        return float("inf")
    nz = v[v > 0]
    otsu = float(threshold_otsu(nz)) if nz.size >= 10 and float(np.max(nz)) > float(np.min(nz)) else float(np.percentile(v, 90))
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med))) * 1.4826
    robust = med + float(k) * mad
    return max(float(factor) * otsu, robust)


def _gi_segment_iba1_seeded(dapi_lbl: Any, iba1_signal: Any, args: argparse.Namespace, logger: logging.Logger) -> Tuple[Any, Any, Any, Dict[str, Any]]:
    from scipy import ndimage as ndi
    from skimage.filters import gaussian, sobel, threshold_otsu
    from skimage.morphology import closing, disk as _disk, remove_small_holes as _rsh, remove_small_objects as _rso
    from skimage.segmentation import watershed, expand_labels as _expand

    x = _gi_robust_normalize(iba1_signal, args)
    bg_radius = um_to_px(float(args.iba1_background_radius_um), args, min_px=3)
    corr, _bg = _gi_background_correct(x, bg_radius)
    thr = _gi_robust_threshold(corr, factor=float(args.iba1_threshold_factor), k=2.5)
    fg = corr > thr
    fg = _rso(fg, min_size=max(2, int(round(float(args.iba1_min_cell_area_um2) / max(float(args.um_per_px) ** 2, 1e-9) / 4.0))))
    fg = closing(fg, _disk(1))
    fg = _rsh(fg, area_threshold=5)

    expand_px = um_to_px(float(args.iba1_nucleus_expand_um), args, min_px=2)
    exp = _expand(np.asarray(dapi_lbl, dtype=np.int32), distance=expand_px)
    nuc_arr = np.unique(np.asarray(dapi_lbl, dtype=np.int32))
    nuc_arr = nuc_arr[nuc_arr > 0].astype(int, copy=False)
    if nuc_arr.size:
        score_arr = np.asarray(
            ndi.mean(corr, labels=np.asarray(exp, dtype=np.int32), index=nuc_arr),
            dtype=float,
        )
        score_arr[~np.isfinite(score_arr)] = 0.0
    else:
        score_arr = np.asarray([], dtype=float)
    if score_arr.size:
        try:
            otsu_n = float(threshold_otsu(score_arr)) if float(np.max(score_arr)) > float(np.min(score_arr)) else float(np.percentile(score_arr, 75))
        except Exception:
            otsu_n = float(np.percentile(score_arr, 75))
        q = float(np.quantile(score_arr, min(max(float(args.iba1_nucleus_min_score_quantile), 0.0), 1.0)))
        nucleus_thr = max(otsu_n, q)
        selected = nuc_arr[score_arr >= nucleus_thr]
    else:
        nucleus_thr = float("nan")
        selected = np.array([], dtype=int)

    marker_lookup = np.zeros(int(np.max(dapi_lbl)) + 1, dtype=np.int32)
    if selected.size:
        marker_lookup[selected] = np.arange(1, selected.size + 1, dtype=np.int32)
    markers = marker_lookup[np.asarray(dapi_lbl, dtype=np.int32)]
    seed_mask = markers > 0
    max_radius_px = um_to_px(float(args.iba1_max_cell_radius_um), args, min_px=8)
    near_seed = ndi.distance_transform_edt(~seed_mask) <= max_radius_px if np.any(seed_mask) else np.zeros_like(seed_mask)
    domain = (fg & near_seed) | seed_mask
    grad = sobel(gaussian(corr, sigma=0.8, preserve_range=True))
    if int(markers.max()) > 0:
        lbl = watershed(grad, markers=markers, mask=domain, compactness=0.0).astype(np.int32)
    else:
        lbl = np.zeros_like(dapi_lbl, dtype=np.int32)

    min_px = int(round(float(args.iba1_min_cell_area_um2) / max(float(args.um_per_px) ** 2, 1e-9)))
    max_px = int(round(float(args.iba1_max_cell_area_um2) / max(float(args.um_per_px) ** 2, 1e-9)))
    keep = np.zeros(int(lbl.max()) + 1, dtype=bool)
    if keep.size:
        keep[0] = False
    kept = 0
    for rp in regionprops(lbl, intensity_image=corr):
        area = int(rp.area)
        fgfrac = float(np.mean(fg[rp.slice][rp.image])) if area > 0 else 0.0
        if area >= max(1, min_px) and (max_px <= 0 or area <= max_px) and fgfrac >= float(args.iba1_min_foreground_fraction):
            keep[int(rp.label)] = True
            kept += 1
    # Preserve watershed instance boundaries. Binarizing retained labels here can
    # merge adjacent cells into long vessel/tissue-like objects in dense Iba1 fields.
    from skimage.segmentation import relabel_sequential
    out_lbl = np.asarray(lbl, dtype=np.int32).copy()
    out_lbl[~keep[out_lbl]] = 0
    out_lbl, _, _ = relabel_sequential(out_lbl)
    out_lbl = out_lbl.astype(np.int32, copy=False)
    info = {
        "iba1_threshold": float(thr), "iba1_nucleus_score_threshold": float(nucleus_thr),
        "iba1_candidate_nuclei": int(len(selected)), "iba1_cells_kept": int(kept),
        "iba1_foreground_fraction": float(np.mean(fg)),
    }
    return out_lbl, fg.astype(bool), corr, info


def _gi_cellpose_iba1(signal: Any, args: argparse.Namespace, logger: logging.Logger) -> Optional[Any]:
    if not args.model_iba1:
        logger.warning("--iba1-segmentation-mode cellpose requires --model-iba1; using seeded watershed.")
        return None
    path = Path(args.model_iba1).expanduser()
    if not path.exists():
        logger.warning("Iba1 model not found: %s", path)
        return None
    try:
        from cellpose import models
        model = models.CellposeModel(gpu=bool(args.gpu), pretrained_model=str(path))
        kwargs: Dict[str, Any] = {"channels": [0, 0]}
        if args.diameter_iba1 is not None:
            kwargs["diameter"] = float(args.diameter_iba1)
        out = model.eval(_gi_robust_normalize(signal, args), **kwargs)
        masks = out[0] if isinstance(out, tuple) else out
        return ensure_label_image(masks)
    except Exception as e:
        logger.warning("Iba1 Cellpose failed; seeded watershed fallback will be used: %s", e)
        return None


def segment_gfap_iba1_missing(sections: List[GISection], args: argparse.Namespace, logger: logging.Logger, outdir: Path) -> Any:
    rows: List[Dict[str, Any]] = []
    for i, s in enumerate(sections, start=1):
        logger.info("[seg %d/%d] %s %s", i, len(sections), s.sample, s.section)
        dapi_img = s.images.get("DAPI")
        iba_img = s.images.get("Iba1")
        if dapi_img is None or iba_img is None:
            rows.append({"sample": s.sample, "section": s.section, "marker": "required", "status": "missing_image"})
            continue
        dapi_seg_path = s.section_dir / f"Image_{s.section}_DAPI_seg.npy"
        iba_seg_path = s.section_dir / f"Image_{s.section}_Iba1_seg.npy"

        dapi_lbl: Optional[Any] = None
        if dapi_seg_path.exists() and not bool(args.overwrite_seg):
            dapi_lbl = ensure_label_image(load_seg_mask(dapi_seg_path))
            rows.append({"sample": s.sample, "section": s.section, "marker": "DAPI", "status": "existing_segmentation_kept", "objects": int(dapi_lbl.max()), "path": str(dapi_seg_path)})
        elif bool(args.segment_missing) or bool(args.overwrite_seg):
            sig = _gi_read_signal(dapi_img, "DAPI")
            mode = str(args.dapi_segmentation_mode)
            if mode in {"auto", "cellpose"}:
                dapi_lbl = _gi_cellpose_dapi(sig, args, logger)
            if dapi_lbl is None:
                dapi_lbl = _gi_watershed_dapi(sig, args)
                source = "watershed"
            else:
                source = "cellpose"
            if not args.dry_run:
                _gi_save_seg(dapi_seg_path, dapi_lbl, source=source)
            rows.append({"sample": s.sample, "section": s.section, "marker": "DAPI", "status": "segmented", "method": source, "objects": int(dapi_lbl.max()), "path": str(dapi_seg_path)})
        else:
            rows.append({"sample": s.sample, "section": s.section, "marker": "DAPI", "status": "missing_not_segmented", "path": str(dapi_seg_path)})
            continue

        reuse_iba = False
        empty_existing = False
        if iba_seg_path.exists() and not bool(args.overwrite_seg):
            lbl0 = ensure_label_image(load_seg_mask(iba_seg_path))
            empty_existing = int(lbl0.max()) == 0 and not bool(args.accept_empty_iba1_seg)
            if not empty_existing:
                rows.append({"sample": s.sample, "section": s.section, "marker": "Iba1", "status": "existing_segmentation_kept", "objects": int(lbl0.max()), "path": str(iba_seg_path)})
                reuse_iba = True
            else:
                logger.warning("Existing Iba1 mask is empty and considered invalid: %s", iba_seg_path)
        if not reuse_iba and (bool(args.segment_missing) or bool(args.overwrite_seg)):
            iba_sig = _gi_read_signal(iba_img, "Iba1")
            iba_lbl = None
            method = "seeded-watershed"
            if str(args.iba1_segmentation_mode) == "cellpose":
                iba_lbl = _gi_cellpose_iba1(iba_sig, args, logger)
                method = "cellpose" if iba_lbl is not None else "seeded-watershed"
            if iba_lbl is None:
                iba_lbl, _fg, _corr, info = _gi_segment_iba1_seeded(dapi_lbl, iba_sig, args, logger)
            else:
                info = {}
            if not args.dry_run:
                _gi_save_seg(iba_seg_path, iba_lbl, source=method, metadata=info)
            row = {"sample": s.sample, "section": s.section, "marker": "Iba1", "status": "empty_regenerated" if empty_existing else "segmented", "method": method, "objects": int(iba_lbl.max()), "path": str(iba_seg_path)}
            row.update(info)
            rows.append(row)
        elif not reuse_iba:
            rows.append({"sample": s.sample, "section": s.section, "marker": "Iba1", "status": "existing_empty_invalid_not_regenerated" if empty_existing else "missing_not_segmented", "objects": 0, "path": str(iba_seg_path)})
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "segmentation_report.csv", index=False)
    return df


def _gi_load_section_labels(s: GISection) -> Tuple[Optional[Any], Optional[Any], List[str]]:
    errs: List[str] = []
    dpath = s.segs.get("DAPI")
    ipath = s.segs.get("Iba1")
    dapi = ensure_label_image(load_seg_mask(dpath)) if dpath and dpath.exists() else None
    iba = ensure_label_image(load_seg_mask(ipath)) if ipath and ipath.exists() else None
    if dapi is None:
        errs.append("missing_DAPI_seg")
    if iba is None:
        errs.append("missing_Iba1_seg")
    if dapi is not None and iba is not None and dapi.shape != iba.shape:
        errs.append("shape_mismatch")
    return dapi, iba, errs


def _gi_fractal_dimension(mask: Any) -> float:
    z = np.asarray(mask, dtype=bool)
    if np.count_nonzero(z) < 4:
        return float("nan")
    p = min(z.shape)
    sizes = 2 ** np.arange(1, int(np.floor(np.log2(max(2, p)))) + 1)
    counts: List[float] = []
    valid_sizes: List[float] = []
    for size in sizes:
        if size > p:
            continue
        pad_y = (-z.shape[0]) % size
        pad_x = (-z.shape[1]) % size
        zp = np.pad(z, ((0, pad_y), (0, pad_x)), mode="constant")
        blocks = zp.reshape(zp.shape[0] // size, size, zp.shape[1] // size, size)
        n = int(np.count_nonzero(blocks.any(axis=(1, 3))))
        if n > 0:
            counts.append(float(n)); valid_sizes.append(float(size))
    if len(counts) < 2:
        return float("nan")
    slope = np.polyfit(np.log(1.0 / np.asarray(valid_sizes)), np.log(np.asarray(counts)), 1)[0]
    return float(slope)


def _gi_skeleton_stats(mask: Any, center_local: Tuple[float, float]) -> Dict[str, float]:
    from scipy import ndimage as ndi
    from skimage.morphology import skeletonize
    sk = skeletonize(np.asarray(mask, dtype=bool))
    neigh = ndi.convolve(sk.astype(np.uint8), np.ones((3, 3), dtype=np.uint8), mode="constant", cval=0) - sk.astype(np.uint8)
    endpoints = int(np.count_nonzero(sk & (neigh == 1)))
    branches = int(np.count_nonzero(sk & (neigh >= 3)))
    length = int(np.count_nonzero(sk))
    cy, cx = center_local
    yy, xx = np.indices(sk.shape)
    rr = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2)
    rmax = int(max(2, np.nanmax(rr[sk]) if np.any(sk) else 2))
    crossings: List[int] = []
    for rad in range(3, rmax + 1, 3):
        ann = sk & (np.abs(rr - rad) <= 0.75)
        crossings.append(int(label(ann, connectivity=2).max()))
    return {
        "skeleton_length_px": float(length), "skeleton_endpoints": float(endpoints),
        "skeleton_branchpoints": float(branches),
        "sholl_max_crossings": float(max(crossings) if crossings else 0),
        "sholl_mean_crossings": float(np.mean(crossings) if crossings else 0.0),
    }


def _gi_extract_patch(signal_corr: Any, cell_mask: Any, centroid: Tuple[float, float], patch_size: int) -> Any:
    from scipy import ndimage as ndi
    from skimage.transform import resize
    from skimage.morphology import skeletonize
    h, w = cell_mask.shape
    cy, cx = centroid
    half = int(max(16, patch_size // 2))
    y0, y1 = int(round(cy)) - half, int(round(cy)) + half
    x0, x1 = int(round(cx)) - half, int(round(cx)) + half
    pad_t, pad_b = max(0, -y0), max(0, y1 - h)
    pad_l, pad_r = max(0, -x0), max(0, x1 - w)
    yy0, yy1 = max(0, y0), min(h, y1)
    xx0, xx1 = max(0, x0), min(w, x1)
    sig = np.pad(signal_corr[yy0:yy1, xx0:xx1], ((pad_t, pad_b), (pad_l, pad_r)), mode="constant")
    m = np.pad(cell_mask[yy0:yy1, xx0:xx1], ((pad_t, pad_b), (pad_l, pad_r)), mode="constant")
    if sig.shape != (2 * half, 2 * half):
        sig = resize(sig, (2 * half, 2 * half), order=1, preserve_range=True, anti_aliasing=True)
        m = resize(m.astype(float), (2 * half, 2 * half), order=0, preserve_range=True) > 0.5
    vals = sig[m]
    if vals.size:
        lo, hi = np.percentile(vals, [1, 99])
    else:
        lo, hi = np.percentile(sig, [1, 99]) if sig.size else (0, 1)
    sig_n = np.clip((sig - lo) / (hi - lo + 1e-6), 0, 1)
    sk = skeletonize(m)
    dist = ndi.distance_transform_edt(m)
    dist = dist / (float(np.max(dist)) + 1e-6)
    ch3 = np.clip(0.65 * sk.astype(float) + 0.35 * dist, 0, 1)
    patch = np.stack([sig_n, m.astype(float), ch3], axis=-1)
    patch = resize(patch, (patch_size, patch_size, 3), order=1, preserve_range=True, anti_aliasing=True)
    return np.clip(patch * 255.0, 0, 255).astype(np.uint8)


def _gi_save_patch(path: Path, patch: Any) -> None:
    from PIL import Image
    ensure_dir(path.parent)
    Image.fromarray(np.asarray(patch, dtype=np.uint8), mode="RGB").save(path)


def _gi_region_for_point(region_info: Dict[str, Any], region_mask: Any, y: float, x: float) -> str:
    # Use the same final region mask for cell/nucleus assignment and for the
    # area/intensity denominator. This avoids numerator/denominator anatomy
    # mismatches. The projected-wall result remains available as an audit vote.
    try:
        r = region_name_from_mask(region_mask, y, x)
        return r if r in {"ARC", "ME", "VMN"} else "OUTSIDE_ARC_ME"
    except Exception:
        r = region_name_from_projected_walls(region_info, y, x)
        return r if r in {"ARC", "ME"} else "OUTSIDE_ARC_ME"


def _gi_quantify_gfap_region(gfap_signal: Any, region_mask: Any, args: argparse.Namespace) -> Tuple[Any, Any, Dict[str, Dict[str, float]]]:
    from skimage.morphology import remove_small_objects as _rso
    # Normalized signal is used only to define the positive mask. Quantitative
    # intensity remains on the raw exported channel scale after background
    # subtraction, so per-section percentile scaling cannot erase biological
    # between-image differences.
    x = _gi_robust_normalize(gfap_signal, args)
    radius = um_to_px(float(args.gfap_background_radius_um), args, min_px=3)
    corr, bg = _gi_background_correct(x, radius)
    raw = np.asarray(gfap_signal, dtype=np.float32)
    corr_raw, _bg_raw = _gi_background_correct(raw, radius)
    tissue = region_mask > 0
    threshold = _gi_robust_threshold(corr[tissue], factor=float(args.gfap_threshold_factor), k=float(args.gfap_robust_k)) if np.any(tissue) else float("inf")
    pos = (corr > threshold) & tissue
    min_px = int(round(float(args.gfap_min_object_area_um2) / max(float(args.um_per_px) ** 2, 1e-9)))
    pos = _rso(pos, min_size=max(1, min_px))
    out: Dict[str, Dict[str, float]] = {}
    delineated = [name for name in ("ARC", "ME", "VMN") if np.any(region_mask == REGION_NAME_TO_CODE[name])]
    for region in delineated:
        code = REGION_NAME_TO_CODE[region]
        rm = region_mask == code
        n = int(np.count_nonzero(rm))
        gp = pos & rm
        pn = int(np.count_nonzero(gp))
        out[region] = {
            "region_area_px": float(n),
            "gfap_positive_area_px": float(pn),
            "gfap_area_fraction": float(pn / n) if n else float("nan"),
            "gfap_raw_mean_intensity": float(np.mean(raw[rm])) if n else float("nan"),
            "gfap_corrected_mean_intensity": float(np.mean(corr_raw[rm])) if n else float("nan"),
            "gfap_corrected_integrated_intensity": float(np.sum(corr_raw[rm])) if n else 0.0,
            "gfap_positive_mean_intensity": float(np.mean(corr_raw[gp])) if pn else float("nan"),
            "gfap_reactivity_index": float((pn / n) * np.mean(corr_raw[gp])) if n and pn else 0.0,
            "gfap_threshold": float(threshold),
            "gfap_intensity_scale": "raw_channel_background_corrected_no_percentile_rescaling",
        }
    return pos.astype(bool), corr, out


def _gi_assign_best_dapi(iba_mask: Any, dapi_lbl: Any) -> int:
    vals = np.asarray(dapi_lbl)[iba_mask]
    vals = vals[vals > 0]
    if vals.size == 0:
        return 0
    c = np.bincount(vals.astype(int))
    return int(np.argmax(c)) if c.size > 1 else 0


MICROGLIA_CELL_QC_VERSION = "biological_cell_qc_v1"
MICROGLIA_CELL_QC_REASONS = (
    "area_below_min", "area_above_max", "weak_iba1_support",
    "no_dapi_association", "image_edge_truncated", "outside_anatomical_tissue",
    "outside_final_arc_me_mask",
    "fragmented_label", "oversized_structure", "linear_or_vascular_structure",
    "multi_nucleus_structure", "invalid_geometry",
)


def _gi_evaluate_microglia_object_qc(
    rp: Any,
    mask_global: Any,
    dapi_lbl: Any,
    iba_fg: Any,
    region_mask: Any,
    args: argparse.Namespace,
) -> Dict[str, Any]:
    """Treatment-blind biological-cell QC applied to every raw Iba1 instance."""
    from scipy import ndimage as ndi

    mask = np.asarray(mask_global, dtype=bool)
    dapi = np.asarray(dapi_lbl, dtype=np.int32)
    anatomy = np.asarray(region_mask) > 0
    um_per_px = max(float(args.um_per_px), 1e-9)
    area_px = int(rp.area)
    area_um2 = float(area_px) * um_per_px * um_per_px
    major_px = float(rp.major_axis_length)
    minor_px = float(rp.minor_axis_length)
    major_um = major_px * um_per_px
    minor_um = minor_px * um_per_px
    aspect = float(major_px / max(minor_px, 1e-6))
    foreground_fraction = float(np.mean(np.asarray(iba_fg, dtype=bool)[mask])) if area_px else 0.0
    anatomy_fraction = float(np.mean(anatomy[mask])) if area_px else 0.0

    minr, minc, maxr, maxc = (int(v) for v in rp.bbox)
    height, width = mask.shape
    edge_margin = max(0, int(args.microglia_qc_edge_margin_px))
    touches_edge = bool(
        minr <= edge_margin or minc <= edge_margin
        or maxr >= height - edge_margin or maxc >= width - edge_margin
    )

    local_components, n_components = ndi.label(
        np.asarray(rp.image, dtype=bool), structure=np.ones((3, 3), dtype=np.uint8)
    )
    component_counts = np.bincount(local_components.ravel())[1:]
    largest_component_fraction = (
        float(component_counts.max() / max(area_px, 1)) if component_counts.size else 0.0
    )

    max_distance_px = max(
        0, int(math.ceil(float(args.microglia_qc_max_nucleus_distance_um) / um_per_px))
    )
    y0, y1 = max(0, minr - max_distance_px), min(height, maxr + max_distance_px)
    x0, x1 = max(0, minc - max_distance_px), min(width, maxc + max_distance_px)
    local_mask = mask[y0:y1, x0:x1]
    local_dapi = dapi[y0:y1, x0:x1]
    direct_labels = dapi[mask]
    direct_labels = direct_labels[direct_labels > 0]
    distance_to_object = ndi.distance_transform_edt(~local_mask)
    near = distance_to_object <= float(max_distance_px)
    near_labels = local_dapi[near]
    near_labels = near_labels[near_labels > 0]
    dapi_nuclei_nearby = int(len(np.unique(near_labels))) if near_labels.size else 0
    label_source = direct_labels if direct_labels.size else near_labels
    if label_source.size:
        counts = np.bincount(label_source.astype(int))
        dapi_lab = int(np.argmax(counts)) if counts.size > 1 else 0
    else:
        dapi_lab = 0
    dapi_direct_overlap_px = int(np.count_nonzero(direct_labels == dapi_lab)) if dapi_lab > 0 else 0
    if dapi_lab > 0 and np.any(local_dapi == dapi_lab):
        dapi_distance_um = float(np.min(distance_to_object[local_dapi == dapi_lab])) * um_per_px
    else:
        dapi_distance_um = float("inf")

    local_center = (
        float(rp.centroid[0] - rp.bbox[0]),
        float(rp.centroid[1] - rp.bbox[1]),
    )
    skeleton_stats = _gi_skeleton_stats(rp.image, local_center)
    branchpoints = int(skeleton_stats["skeleton_branchpoints"])
    endpoints = int(skeleton_stats["skeleton_endpoints"])
    linear_structure = bool(
        major_um >= float(args.microglia_qc_linear_min_major_axis_um)
        and aspect >= float(args.microglia_qc_linear_min_aspect_ratio)
        and branchpoints <= 1 and endpoints <= 3
    )
    large_or_long = bool(
        major_um >= float(args.microglia_qc_linear_min_major_axis_um)
        or area_um2 >= 0.25 * float(args.iba1_max_cell_area_um2)
    )
    multi_nucleus_structure = bool(
        int(args.microglia_qc_max_dapi_nuclei_touched) >= 0
        and dapi_nuclei_nearby > int(args.microglia_qc_max_dapi_nuclei_touched)
        and large_or_long
    )
    finite_geometry = bool(
        np.isfinite(area_um2) and np.isfinite(major_um)
        and np.isfinite(minor_um) and np.isfinite(aspect)
    )

    reasons: List[str] = []
    if area_um2 < float(args.iba1_min_cell_area_um2): reasons.append("area_below_min")
    if area_um2 > float(args.iba1_max_cell_area_um2): reasons.append("area_above_max")
    if foreground_fraction < float(args.iba1_min_foreground_fraction): reasons.append("weak_iba1_support")
    if dapi_lab <= 0 or dapi_distance_um > float(args.microglia_qc_max_nucleus_distance_um): reasons.append("no_dapi_association")
    if touches_edge: reasons.append("image_edge_truncated")
    if anatomy_fraction < float(args.microglia_qc_min_anatomy_fraction): reasons.append("outside_anatomical_tissue")
    if largest_component_fraction < float(args.microglia_qc_min_largest_component_fraction): reasons.append("fragmented_label")
    if major_um > float(args.microglia_qc_max_major_axis_um): reasons.append("oversized_structure")
    if linear_structure: reasons.append("linear_or_vascular_structure")
    if multi_nucleus_structure: reasons.append("multi_nucleus_structure")
    if not finite_geometry: reasons.append("invalid_geometry")

    out: Dict[str, Any] = {
        "cell_qc_version": MICROGLIA_CELL_QC_VERSION,
        "cell_qc_pass": bool(not reasons),
        "cell_qc_rejection_reasons": ";".join(reasons),
        "object_iba1_foreground_fraction": foreground_fraction,
        "dapi_nucleus_label": int(dapi_lab),
        "dapi_direct_overlap_px": dapi_direct_overlap_px,
        "dapi_nucleus_distance_um": dapi_distance_um,
        "dapi_nuclei_nearby": dapi_nuclei_nearby,
        "mask_inside_anatomy_fraction": anatomy_fraction,
        "touches_image_edge": touches_edge,
        "connected_components": int(n_components),
        "largest_component_fraction": largest_component_fraction,
        "major_axis_length_um": major_um,
        "minor_axis_length_um": minor_um,
        "axis_aspect_ratio": aspect,
        "bbox_height_px": int(maxr - minr),
        "bbox_width_px": int(maxc - minc),
    }
    out.update(skeleton_stats)
    for reason in MICROGLIA_CELL_QC_REASONS:
        out[f"qc_reject_{reason}"] = bool(reason in reasons)
    return out


def _gi_require_qc_accepted(cell_df: Any, stage: str) -> None:
    """Fail closed if an object-level biological stage receives rejected cells."""
    if cell_df is None or cell_df.empty:
        return
    if "cell_qc_pass" not in cell_df.columns:
        raise RuntimeError(f"{stage}: required cell_qc_pass column is missing")
    accepted = cell_df["cell_qc_pass"].fillna(False).astype(bool)
    if not bool(accepted.all()):
        raise RuntimeError(
            f"{stage}: {int((~accepted).sum())} rejected Iba1 objects reached a biological-cell stage"
        )



GI_MANUAL_REGIONS = ("ARC", "ME", "VMN")
GI_MANUAL_LABEL_CODES = {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}
GI_LEGACY_LABEL_CODES = {"background": 0, "ARC": 1, "ME": 2}
GI_REGION_STATES = {"drawn", "explicitly_absent", "not_delineated"}


def _gi_review_region_status(
    review: Dict[str, Any], mask: Any
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Validate schema-v3 partial-region provenance, retaining legacy receipts."""
    receipt_codes = review.get("label_codes")
    legacy = receipt_codes == GI_LEGACY_LABEL_CODES
    if receipt_codes not in (GI_MANUAL_LABEL_CODES, GI_LEGACY_LABEL_CODES):
        raise ValueError(
            f"label_codes must equal {GI_MANUAL_LABEL_CODES} "
            f"(or legacy {GI_LEGACY_LABEL_CODES})"
        )
    raw_reasons = review.get("region_absence_reasons", {})
    if not isinstance(raw_reasons, dict):
        raise ValueError("region_absence_reasons must be an object")
    reasons = {
        name: str(raw_reasons.get(name, "")).strip()
        for name in GI_MANUAL_REGIONS
    }
    if not reasons["ME"]:
        reasons["ME"] = str(review.get("me_absent_reason", "")).strip()

    if legacy:
        me_absent = review.get("me_absent") is True
        presence = review.get("region_presence")
        if isinstance(presence, dict) and presence.get("ME") is False:
            me_absent = True
        status = {
            "ARC": "drawn" if np.any(mask == 1) else "not_delineated",
            "ME": (
                "drawn" if np.any(mask == 2) else
                "explicitly_absent" if me_absent else "not_delineated"
            ),
            "VMN": "not_delineated",
        }
    else:
        raw_status = review.get("region_status")
        if not isinstance(raw_status, dict):
            raise ValueError("schema-v3 receipt lacks region_status")
        status = {name: str(raw_status.get(name, "")) for name in GI_MANUAL_REGIONS}
        invalid = {k: v for k, v in status.items() if v not in GI_REGION_STATES}
        if invalid:
            raise ValueError(f"invalid region_status entries: {invalid}")
        expected_lists = {
            "delineated_regions": [n for n in GI_MANUAL_REGIONS if status[n] == "drawn"],
            "omitted_regions": [n for n in GI_MANUAL_REGIONS if status[n] == "not_delineated"],
            "explicitly_absent_regions": [n for n in GI_MANUAL_REGIONS if status[n] == "explicitly_absent"],
        }
        for field, expected in expected_lists.items():
            if review.get(field) != expected:
                raise ValueError(f"{field} must equal {expected}")
        presence = review.get("region_presence")
        expected_presence = {
            n: True if status[n] == "drawn" else False if status[n] == "explicitly_absent" else None
            for n in GI_MANUAL_REGIONS
        }
        if not isinstance(presence, dict) or {
            n: presence.get(n) for n in GI_MANUAL_REGIONS
        } != expected_presence:
            raise ValueError(f"region_presence must equal {expected_presence}")

    for name in GI_MANUAL_REGIONS:
        area = int(np.count_nonzero(mask == GI_MANUAL_LABEL_CODES[name]))
        state = status[name]
        if (state == "drawn") != (area > 0):
            raise ValueError(
                f"{name} status {state!r} conflicts with mask area {area}px"
            )
        if state == "explicitly_absent" and not reasons[name]:
            raise ValueError(f"explicit {name} absence requires a reason")
        if state != "explicitly_absent" and reasons[name]:
            raise ValueError(
                f"{name} absence reason is present but status is {state}"
            )
    if not any(status[name] == "drawn" for name in GI_MANUAL_REGIONS):
        raise ValueError("mask must delineate at least one of ARC, ME, or VMN")
    return status, reasons


class ManualRegionMaskError(RuntimeError):
    """A required or supplied manual ARC/ME/VMN review pair is unusable."""


def _gi_lexical_absolute_path(path: Path, *, base: Optional[Path] = None) -> Path:
    """Normalize historical receipt paths without filesystem probes."""
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = Path(base) / candidate if base is not None else Path.cwd() / candidate
    return Path(os.path.abspath(os.fspath(candidate)))


def _gi_resolve_relocated_figure_path(
    value: Any,
    *,
    repository_root: Path = PAPER_ROOT,
) -> Path:
    """Resolve an immutable pre-renumbering Figure 5 receipt inside Figure 5.

    The accepted GFAP/Iba1 receipts intentionally retain their historical
    absolute ``Fig5`` paths so their bytes and recorded receipt hashes are not
    rewritten during figure renumbering. Only that exact repository-local
    prefix is translated, and downstream code still verifies the referenced
    image or proposal by SHA-256 before accepting it.
    """
    root = Path(repository_root).expanduser().resolve()
    candidate = _gi_lexical_absolute_path(Path(str(value)), base=root)
    parts = candidate.parts
    figure_indices = [index for index, part in enumerate(parts) if part == "Fig5"]
    if len(figure_indices) != 1:
        return candidate
    relative = Path(*parts[figure_indices[0] + 1 :])
    if not relative.parts:
        return candidate
    figure_root = (root / "Fig5").resolve()
    relocated = _gi_lexical_absolute_path(figure_root / relative)
    if not relocated.is_relative_to(figure_root):
        return candidate
    # Always bind a historical Fig5 receipt to this repository. Returning an
    # existing path from the original checkout would let an adjacent clean-room
    # replay read outside its reconstructed tree. A missing local file must
    # therefore remain missing and fail the downstream hash/path checks.
    return relocated


def _gi_review_provenance_path_matches(
    value: Any,
    expected_path: Path,
    *,
    repository_root: Path = PAPER_ROOT,
) -> bool:
    """Match an HIL source path independently of the process working directory.

    Current receipts contain absolute source paths.  Earlier HIL receipts may
    contain only a leaf filename or a path relative to the Paper repository.
    Resolve those representations against the expected source section and the
    repository root, never against an arbitrary current working directory.

    The caller separately verifies the exact DAPI SHA-256, dimensions, section
    identity, and mask hash.  Accepting an equivalent path representation does
    not weaken the image-to-receipt binding.
    """
    text = str(value or "").strip()
    if not text:
        return False
    expected = Path(expected_path).expanduser().resolve()
    candidate = Path(text).expanduser()
    candidates = [candidate] if candidate.is_absolute() else [
        expected.parent / candidate,
        Path(repository_root).expanduser().resolve() / candidate,
    ]
    for path in candidates:
        try:
            if _gi_resolve_relocated_figure_path(path, repository_root=repository_root) == expected:
                return True
        except (OSError, RuntimeError):
            continue
    return False


def _gi_manual_region_review_candidates(
    s: GISection,
    root: Path,
) -> List[Tuple[str, Path, Path]]:
    """Return current and legacy receipt paths for one native source section.

    The reviewer writes annotations under the exact source folder name.  This
    matters for repetition folders such as ``FR6-2(R)``: they are merged only
    later at animal-level statistics and must not collide with ``FR6-2`` during
    image review or per-section provenance.  A legacy base-animal path remains
    readable when the exact-source pair is completely absent.
    """
    source_sample = str(s.sample).strip()
    section = str(s.section).strip()
    for value, field in ((source_sample, "source sample"), (section, "section")):
        if (
            not value
            or value in {".", ".."}
            or "/" in value
            or "\\" in value
            or Path(value).name != value
        ):
            raise ManualRegionMaskError(
                f"unsafe {field} in native GFAP review identity: {value!r}"
            )

    base, _is_rep, _rep_lab = strip_repetition_suffix(source_sample)
    animals = [source_sample]
    if base != source_sample:
        animals.append(base)
    candidates: List[Tuple[str, Path, Path]] = []
    for animal in animals:
        pair_dir = Path(root) / animal / section
        stem = f"{animal}_{section}_ARC_ME_labels"
        candidates.append(
            (animal, pair_dir / f"{stem}.tif", pair_dir / f"{stem}.json")
        )
    return candidates


def _gi_load_manual_region_mask(
    s: GISection,
    dapi_labels: Any,
    args: argparse.Namespace,
    logger: logging.Logger,
) -> Tuple[Optional[Any], Optional[Dict[str, Any]]]:
    """Load one native HIL-accepted 0/1/2/3 ARC/ME/VMN mask.

    The canonical pair preserves the exact source-folder identity::

        <manual-region-dir>/<source-sample>/<section>/<source-sample>_<section>_ARC_ME_labels.tif
        <manual-region-dir>/<source-sample>/<section>/<source-sample>_<section>_ARC_ME_labels.json

    A base-animal directory is accepted only as a backward-compatible fallback
    when the exact-source pair is completely absent.  The JSON receipt must
    contain ``accepted: true``, a nonblank reviewer, an ISO-8601 timestamp,
    image dimensions, the current 0/1/2/3 label contract (or its legacy 0/1/2
    form), and hashes matching the paired TIFF and exact native DAPI. Every
    failure raises; unreviewed automated fallback is forbidden for this assay.
    """
    dapi_array = np.asarray(dapi_labels)
    dapi_shape = tuple(int(v) for v in dapi_array.shape)
    root_arg = getattr(args, "manual_region_dir", None)
    human_in_the_loop = True
    required = True
    if not root_arg:
        raise ManualRegionMaskError(
            "strict HIL ARC/ME/VMN review requires --manual-region-dir"
        )

    source_sample = str(s.sample).strip()
    base, _is_rep, _rep_lab = strip_repetition_suffix(source_sample)
    candidates = _gi_manual_region_review_candidates(s, Path(root_arg))
    _current_animal, current_mask, current_review = candidates[0]

    # A partial current-format pair is an interrupted/invalid save and must not
    # be hidden by a legacy fallback.  Otherwise prefer the exact source-folder
    # pair and only then consider the pre-fix base-animal layout.
    if current_mask.is_file() != current_review.is_file():
        missing = current_review if current_mask.is_file() else current_mask
        raise ManualRegionMaskError(
            f"Manual ARC/ME/VMN review rejected for {s.sample} {s.section}: "
            f"incomplete exact-source receipt pair; missing {missing}"
        )
    selected = next(
        (
            (animal, mask, review)
            for animal, mask, review in candidates
            if mask.is_file() and review.is_file()
        ),
        None,
    )
    if selected is None:
        expected_pairs = "; ".join(
            f"{mask} + {review}" for _animal, mask, review in candidates
        )
        raise ManualRegionMaskError(
            f"Manual ARC/ME/VMN review rejected for {s.sample} {s.section}: "
            f"missing required paired TIFF+JSON; expected {expected_pairs}"
        )
    receipt_animal_path, mask_path, review_path = selected

    def reject(reason: str) -> Tuple[Optional[Any], Optional[Dict[str, Any]]]:
        message = f"Manual ARC/ME/VMN review rejected for {s.sample} {s.section}: {reason}"
        raise ManualRegionMaskError(message)

    try:
        review = json.loads(review_path.read_text(encoding="utf-8"))
    except Exception as exc:
        return reject(f"cannot parse review JSON {review_path}: {exc}")
    if not isinstance(review, dict):
        return reject(f"review JSON must contain one object: {review_path}")
    if review.get("accepted") is not True:
        return reject(f"review JSON does not have accepted=true: {review_path}")
    reviewer = str(review.get("reviewer", "")).strip()
    if not reviewer:
        return reject(f"review JSON has no nonblank reviewer: {review_path}")
    accepted_at = str(review.get("accepted_at", "")).strip()
    try:
        dt.datetime.fromisoformat(accepted_at.replace("Z", "+00:00"))
    except Exception:
        return reject(f"review JSON accepted_at is not ISO-8601: {review_path}")
    try:
        receipt_shape = (int(review["image_height_px"]), int(review["image_width_px"]))
    except Exception:
        return reject(f"review JSON lacks integer image_height_px/image_width_px: {review_path}")
    if receipt_shape != tuple(dapi_shape):
        return reject(
            f"review dimensions {receipt_shape} do not match DAPI {tuple(dapi_shape)}: {review_path}"
        )

    source_dapi_path = s.images.get("DAPI") if s.images else None
    resolved_source_dapi = Path(source_dapi_path).resolve() if source_dapi_path else None
    actual_source_dapi_sha256: Optional[str] = None
    source_dapi_hash_error: Optional[str] = None
    if resolved_source_dapi is not None and resolved_source_dapi.is_file():
        try:
            actual_source_dapi_sha256 = sha256_file(resolved_source_dapi)
        except Exception as exc:
            source_dapi_hash_error = str(exc)
    elif resolved_source_dapi is not None:
        source_dapi_hash_error = f"source DAPI does not exist: {resolved_source_dapi}"

    receipt_source_dapi_sha256: Optional[str] = None
    source_dapi_sha256_verified = False
    if human_in_the_loop and "source_dapi_sha256" not in review:
        return reject(f"native HIL review receipt lacks source_dapi_sha256: {review_path}")
    if "source_dapi_sha256" in review:
        receipt_source_dapi_sha256 = str(review.get("source_dapi_sha256", "")).strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", receipt_source_dapi_sha256):
            return reject(f"review JSON source_dapi_sha256 is invalid: {review_path}")
        if actual_source_dapi_sha256 is None:
            return reject(
                f"cannot verify source_dapi_sha256 against {resolved_source_dapi}: "
                f"{source_dapi_hash_error or 'source DAPI unavailable'}"
            )
        if receipt_source_dapi_sha256 != actual_source_dapi_sha256:
            return reject(
                "source_dapi_sha256 mismatch "
                f"(receipt={receipt_source_dapi_sha256}, actual={actual_source_dapi_sha256}): "
                f"{resolved_source_dapi}"
            )
        source_dapi_sha256_verified = True

    accepted_source = str(review.get("accepted_source", "")).strip()
    allowed_accepted_sources = {"automated_mask_unchanged", "manual_replacement_polygons"}
    if accepted_source not in allowed_accepted_sources:
        return reject(
            "native human-in-the-loop receipt requires accepted_source in "
            f"{sorted(allowed_accepted_sources)}: {review_path}"
        )
    if str(review.get("dataset", "")).strip().lower() != "gfap":
        return reject(f"native GFAP HIL review receipt requires dataset=gfap: {review_path}")
    receipt_animal = str(review.get("animal", "")).strip()
    allowed_receipt_animals = {source_sample, base, receipt_animal_path}
    if (
        receipt_animal not in allowed_receipt_animals
        or str(review.get("section", "")).strip() != s.section
    ):
        return reject(
            "review identity does not match native GFAP source section "
            f"{source_sample}/{s.section}: {review_path}"
        )
    source_paths = review.get("source_paths") if isinstance(review.get("source_paths"), dict) else {}
    receipt_dapi_text = str(source_paths.get("dapi", "") or "").strip()
    if (
        not receipt_dapi_text
        or resolved_source_dapi is None
        or not _gi_review_provenance_path_matches(
            receipt_dapi_text, resolved_source_dapi
        )
    ):
        return reject(f"review source_paths.dapi does not match the native source DAPI: {review_path}")
    source_manifest_path_text = str(source_paths.get("manifest", "") or "").strip()
    source_manifest_sha256 = str(review.get("source_manifest_sha256", "") or "").strip()
    if source_manifest_path_text or source_manifest_sha256 or review.get("registration") is not None:
        return reject(f"native GFAP review cannot cite registration or a review manifest: {review_path}")

    proposal_mask_path_text = str(source_paths.get("automated_mask", "") or "").strip()
    proposal_mask_sha256 = str(
        review.get("automated_mask_sha256")
        or review.get("proposal_mask_sha256")
        or ""
    ).strip().lower()
    proposal_mask_verified = False
    if accepted_source == "automated_mask_unchanged":
        confirmation = (
            review.get("automated_mask_confirmation")
            or review.get("proposal_confirmation")
        )
        automated_provenance = (
            review.get("automated_mask_provenance")
            if isinstance(review.get("automated_mask_provenance"), dict)
            else {}
        )
        if confirmation != "reviewed_image_by_image":
            return reject(f"unchanged proposal lacks image-by-image confirmation: {review_path}")
        if (
            str(automated_provenance.get("dataset", "")).strip().lower() != "gfap"
            or str(automated_provenance.get("kind", "")).strip()
            != "native_gfap_automated_mask"
        ):
            return reject(f"unchanged proposal is not explicitly native GFAP: {review_path}")
        if not proposal_mask_path_text or not re.fullmatch(r"[0-9a-f]{64}", proposal_mask_sha256):
            return reject(f"unchanged proposal has incomplete proposal provenance: {review_path}")
        proposal_path = _gi_resolve_relocated_figure_path(
            proposal_mask_path_text, repository_root=PAPER_ROOT
        )
        if not proposal_path.is_file():
            return reject(f"accepted unchanged proposal is missing: {proposal_path}")
        if sha256_file(proposal_path) != proposal_mask_sha256:
            return reject(f"accepted unchanged proposal SHA-256 mismatch: {proposal_path}")
        proposal_mask_verified = True
    else:
        if proposal_mask_path_text:
            proposal_path = _gi_resolve_relocated_figure_path(
                proposal_mask_path_text, repository_root=PAPER_ROOT
            )
        else:
            proposal_path = None

    try:
        raw = np.asarray(read_image(mask_path))
    except Exception as exc:
        return reject(f"cannot read manual TIFF {mask_path}: {exc}")
    mask = np.squeeze(raw)
    if mask.ndim != 2:
        return reject(f"mask must be one 2-D plane; got shape {raw.shape}: {mask_path}")
    if tuple(mask.shape) != tuple(dapi_shape):
        return reject(
            f"mask shape {tuple(mask.shape)} does not match DAPI {tuple(dapi_shape)}: {mask_path}"
        )
    if not np.issubdtype(mask.dtype, np.number) or not np.all(np.isfinite(mask)):
        return reject(f"mask must contain finite numeric codes: {mask_path}")
    if not np.all(mask == np.rint(mask)):
        return reject(f"mask contains non-integer values: {mask_path}")
    codes = {int(v) for v in np.unique(mask)}
    invalid_codes = sorted(codes.difference(set(GI_MANUAL_LABEL_CODES.values())))
    if invalid_codes:
        return reject(
            f"mask contains invalid codes {invalid_codes}; expected "
            f"{sorted(GI_MANUAL_LABEL_CODES.values())}: {mask_path}"
        )
    mask = np.asarray(mask, dtype=np.uint8)
    try:
        region_status, absence_reasons = _gi_review_region_status(review, mask)
    except ValueError as exc:
        return reject(str(exc))
    region_area_px = {
        name: int(np.count_nonzero(mask == GI_MANUAL_LABEL_CODES[name]))
        for name in GI_MANUAL_REGIONS
    }
    arc_px, me_px, vmn_px = (region_area_px[n] for n in GI_MANUAL_REGIONS)
    me_absent = region_status["ME"] == "explicitly_absent"
    me_absent_reason = absence_reasons["ME"]
    if accepted_source == "automated_mask_unchanged":
        try:
            proposal_mask = np.squeeze(np.asarray(read_image(proposal_path)))
        except Exception as exc:
            return reject(f"cannot read accepted native automated mask {proposal_path}: {exc}")
        if (
            proposal_mask.shape != mask.shape
            or not np.issubdtype(proposal_mask.dtype, np.number)
            or not np.all(np.isfinite(proposal_mask))
            or not np.all(proposal_mask == np.rint(proposal_mask))
            or not np.array_equal(proposal_mask.astype(np.uint8, copy=False), mask)
        ):
            return reject(
                f"accepted TIFF is not pixel-identical to unchanged native automated mask: {proposal_path}"
            )

    actual_mask_sha256 = hashlib.sha256(mask_path.read_bytes()).hexdigest()
    receipt_mask_sha256 = str(review.get("mask_sha256", "")).strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", receipt_mask_sha256):
        return reject(f"review JSON lacks a valid mask_sha256: {review_path}")
    if receipt_mask_sha256 != actual_mask_sha256:
        return reject(
            f"mask_sha256 mismatch (receipt={receipt_mask_sha256}, actual={actual_mask_sha256}): {mask_path}"
        )

    # Descriptive, non-destructive anatomy QC. These metrics make gross drawing
    # problems auditable without invalidating a legitimate human decision that
    # ME is absent in a particular section.
    roi = mask > 0
    region_masks = {
        name: mask == REGION_NAME_TO_CODE[name]
        for name in GI_MANUAL_REGIONS
    }
    dapi_binary = dapi_array > 0
    tissue = smooth_tissue_mask_from_dapi(dapi_array, args)

    def overlap_fraction(numerator_mask: Any, denominator_mask: Any) -> Optional[float]:
        denominator_n = int(np.count_nonzero(denominator_mask))
        if denominator_n == 0:
            return None
        return float(np.count_nonzero(np.asarray(numerator_mask) & np.asarray(denominator_mask)) / denominator_n)

    roi_px = int(np.count_nonzero(roi))
    image_px = int(mask.size)
    dapi_positive_px = int(np.count_nonzero(dapi_binary))
    tissue_px = int(np.count_nonzero(tissue))
    mask_tissue_overlap = overlap_fraction(tissue, roi)
    tissue_overlap = {
        name: overlap_fraction(tissue, region_masks[name])
        for name in GI_MANUAL_REGIONS
    }
    dapi_capture_fraction = overlap_fraction(roi, dapi_binary)
    dapi_occupancy_fraction = overlap_fraction(dapi_binary, roi)
    components = {
        name: (
            int(label(region_masks[name], connectivity=1).max())
            if label is not None and region_area_px[name] else 0
        )
        for name in GI_MANUAL_REGIONS
    }
    manual_qc_flags: List[str] = []
    if mask_tissue_overlap is not None and mask_tissue_overlap < 0.05:
        manual_qc_flags.append("low_mask_overlap_with_dapi_tissue")
    if dapi_positive_px and dapi_capture_fraction == 0.0:
        manual_qc_flags.append("no_dapi_label_pixels_inside_mask")
    for name in GI_MANUAL_REGIONS:
        threshold = 20 if name == "ARC" else 10
        if components[name] > threshold:
            manual_qc_flags.append(f"high_{name.lower()}_component_count")

    resolved_mask = mask_path.resolve()
    resolved_review = review_path.resolve()
    info: Dict[str, Any] = {
        "region_mode": "manual",
        "region_ai_model": "human_reviewed_arc_me_vmn_tiff",
        "region_source": "manual_region_tiff",
        "manual_region_review_status": "accepted",
        "manual_region_path": str(resolved_mask),
        "manual_region_review_path": str(resolved_review),
        "manual_region_sha256": actual_mask_sha256,
        "manual_region_review_sha256": hashlib.sha256(resolved_review.read_bytes()).hexdigest(),
        "manual_region_reviewer": reviewer,
        "manual_region_accepted_at": accepted_at,
        "manual_region_review_session": str(review.get("session", "")).strip() or None,
        "manual_region_receipt_schema_version": review.get("schema_version"),
        "manual_region_accepted_source": accepted_source or None,
        "manual_region_proposal_mask_path": str(proposal_path) if proposal_path else None,
        "manual_region_proposal_mask_sha256": proposal_mask_sha256 or None,
        "manual_region_proposal_mask_sha256_verified": proposal_mask_verified,
        "manual_region_code_contract": "0=background;1=ARC;2=ME;3=VMN",
        "manual_region_delineated_regions": ";".join(
            name for name in GI_MANUAL_REGIONS if region_status[name] == "drawn"
        ),
        "manual_region_omitted_regions": ";".join(
            name for name in GI_MANUAL_REGIONS if region_status[name] == "not_delineated"
        ),
        "manual_region_explicitly_absent_regions": ";".join(
            name for name in GI_MANUAL_REGIONS if region_status[name] == "explicitly_absent"
        ),
        "manual_region_required": required,
        "manual_region_human_in_the_loop": human_in_the_loop,
        "manual_region_me_absent": me_absent,
        "manual_region_me_absent_reason": me_absent_reason or None,
        "manual_region_source_dapi_path": str(resolved_source_dapi) if resolved_source_dapi else None,
        "manual_region_source_dapi_sha256": actual_source_dapi_sha256,
        "manual_region_receipt_source_dapi_sha256": receipt_source_dapi_sha256,
        "manual_region_source_dapi_sha256_verified": source_dapi_sha256_verified,
        "manual_region_source_dapi_hash_error": source_dapi_hash_error,
        "manual_region_source_paths_json": (
            json.dumps(review.get("source_paths"), sort_keys=True)
            if review.get("source_paths") is not None else None
        ),
        "manual_region_proposal_provenance_json": (
            json.dumps(review.get("proposal_provenance"), sort_keys=True)
            if review.get("proposal_provenance") is not None else None
        ),
        "arc_area_px": arc_px,
        "me_area_px": me_px,
        "vmn_area_px": vmn_px,
        "me_floor_present": bool(me_px > 0 and not me_absent),
        "manual_region_roi_area_px": roi_px,
        "manual_region_roi_area_fraction": float(roi_px / image_px) if image_px else None,
        "manual_region_arc_area_fraction": float(arc_px / image_px) if image_px else None,
        "manual_region_me_area_fraction": float(me_px / image_px) if image_px else None,
        "manual_region_vmn_area_fraction": float(vmn_px / image_px) if image_px else None,
        "manual_region_dapi_positive_px": dapi_positive_px,
        "manual_region_dapi_tissue_px": tissue_px,
        "manual_region_mask_tissue_overlap_fraction": mask_tissue_overlap,
        "manual_region_dapi_capture_fraction": dapi_capture_fraction,
        "manual_region_dapi_occupancy_fraction": dapi_occupancy_fraction,
        "manual_region_qc_status": "warning" if manual_qc_flags else "ok",
        "manual_region_qc_flags": ";".join(manual_qc_flags),
    }
    for name in GI_MANUAL_REGIONS:
        lower = name.lower()
        info[f"manual_region_{lower}_status"] = region_status[name]
        info[f"manual_region_{lower}_absence_reason"] = absence_reasons[name] or None
        info[f"manual_region_{lower}_tissue_overlap_fraction"] = tissue_overlap[name]
        info[f"manual_region_{lower}_components"] = components[name]
    logger.warning(
        "Using HIL-accepted ARC/ME/VMN mask for %s %s: %s "
        "(reviewer=%s, accepted_at=%s, statuses=%s, areas_px=%s, QC=%s)",
        s.sample, s.section, resolved_mask, reviewer, accepted_at,
        region_status, region_area_px, info["manual_region_qc_status"],
    )
    return mask, info


def _gi_audit_manual_region_coverage(
    sections: Sequence[GISection],
    args: argparse.Namespace,
    outdir: Path,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Validate the complete native 61-section review set, never an inferred fallback."""
    rows: List[Dict[str, Any]] = []
    pending: List[Dict[str, Any]] = []
    identities = [(str(section.sample), str(section.section)) for section in sections]
    unique_identities = set(identities)
    if (
        len(sections) != GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
        or len(unique_identities) != GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
    ):
        discovery_row = {
            "sample": "__dataset__",
            "section": "__discovery__",
            "review_status": "pending",
            "reason": (
                f"expected {GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} unique GFAP/Iba1 section "
                f"identities, but discovered rows={len(sections)}, "
                f"unique={len(unique_identities)}"
            ),
        }
        rows.append(discovery_row)
        pending.append(discovery_row)
    for section in sections:
        row: Dict[str, Any] = {"sample": section.sample, "section": section.section}
        dapi_seg = section.segs.get("DAPI")
        if dapi_seg is None or not dapi_seg.is_file():
            row.update({"review_status": "pending", "reason": "missing DAPI segmentation required for receipt/QC audit"})
            rows.append(row)
            pending.append(row)
            continue
        try:
            dapi_labels = ensure_label_image(load_seg_mask(dapi_seg))
            _mask, info = _gi_load_manual_region_mask(section, dapi_labels, args, logger)
            row.update({
                "review_status": "accepted",
                "reason": "",
                "accepted_source": (info or {}).get("manual_region_accepted_source"),
                "manual_region_path": (info or {}).get("manual_region_path"),
                "manual_region_sha256": (info or {}).get("manual_region_sha256"),
                "manual_region_review_path": (info or {}).get("manual_region_review_path"),
                "manual_region_review_sha256": (info or {}).get("manual_region_review_sha256"),
                "reviewer": (info or {}).get("manual_region_reviewer"),
                "accepted_at": (info or {}).get("manual_region_accepted_at"),
                "me_absent": (info or {}).get("manual_region_me_absent"),
                "me_absent_reason": (info or {}).get("manual_region_me_absent_reason"),
                "source_dapi_path": (info or {}).get("manual_region_source_dapi_path"),
                "source_dapi_sha256": (info or {}).get("manual_region_source_dapi_sha256"),
                "source_dapi_sha256_verified": (info or {}).get("manual_region_source_dapi_sha256_verified"),
                "proposal_mask_path": (info or {}).get("manual_region_proposal_mask_path"),
                "proposal_mask_sha256": (info or {}).get("manual_region_proposal_mask_sha256"),
                "proposal_mask_sha256_verified": (info or {}).get("manual_region_proposal_mask_sha256_verified"),
                "mask_tissue_overlap_fraction": (info or {}).get("manual_region_mask_tissue_overlap_fraction"),
                "manual_qc_status": (info or {}).get("manual_region_qc_status"),
                "manual_qc_flags": (info or {}).get("manual_region_qc_flags"),
            })
        except ManualRegionMaskError as exc:
            row.update({"review_status": "pending", "reason": str(exc)})
            pending.append(row)
        rows.append(row)
    ensure_dir(outdir)
    pd.DataFrame(rows).to_csv(outdir / "human_review_coverage.csv", index=False)
    logger.warning(
        "Native HIL ARC/ME/VMN review coverage: accepted=%d/%d required | pending=%d",
        sum(row.get("review_status") == "accepted" for row in rows),
        GFAP_IBA1_EXPECTED_REVIEW_SECTIONS,
        len(pending),
    )
    return rows, pending


def _gi_render_figure_5(args: argparse.Namespace, outdir: Path, work_root: Path,
                        logger: logging.Logger) -> None:
    """make_figure_5: render Figure 5 from this run into a fresh, absent directory."""
    renderer = PAPER_ROOT / "scripts" / "Fig5" / "02_make_figure_5_gfap_iba1_microglia.py"
    if not renderer.is_file():
        logger.warning("Figure 5 renderer is missing, so make_figure_5 was skipped: %s", renderer)
        return
    root = Path(args.figure_output_root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = root / f"paper_figure6_render_{stamp}_{os.getpid()}"
    if target.exists():
        logger.warning("Figure 5 render target already exists, so make_figure_5 was skipped: %s", target)
        return
    command = [
        sys.executable, str(renderer),
        "--analysis-dir", str(outdir),
        "--input-root", str(Path(work_root).expanduser().resolve()),
        "--manual-region-dir", str(Path(args.manual_region_dir).expanduser().resolve()),
        "--output-dir", str(target),
        "--um-per-px", str(float(args.um_per_px)),
        "--dpi", "600", "--formats", "pdf,png",
    ]
    logger.warning("make_figure_5: %s", shlex.join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode:
        logger.warning(
            "make_figure_5 failed with status %d; the analysis outputs in %s are unchanged.",
            completed.returncode, outdir,
        )
        return
    logger.warning("Figure 5 written to %s", target)
    for suffix in (".pdf", ".png"):
        master = target / f"Figure_GFAP_Iba1_microglia_ARC_ME{suffix}"
        if master.is_file():
            logger.warning("Figure 5 master graphic: %s", master)


def _gi_human_review_command(
    args: argparse.Namespace,
    work_root: Path,
    review_analysis_dir: Path,
) -> List[str]:
    annotator = PAPER_ROOT / "scripts" / "shared" / "01_annotate_regions.py"
    return [
        sys.executable, str(annotator),
        "--dataset", "gfap",
        "--input-root", str(work_root),
        "--analysis-dir", str(review_analysis_dir),
        "--output-dir", str(Path(args.manual_region_dir).resolve()),
        "--host", str(args.human_review_host),
        "--port", str(int(args.human_review_port)),
        # Reopening the reviewer means segmenting sections again, including ones
        # already accepted, so replacement is always enabled. The reviewer still
        # rewrites each JSON/TIFF pair atomically and deletes nothing.
        "--force",
        "--figure-output-root", str(Path(args.figure_output_root).expanduser()),
    ]


def _gi_microglia_review_command(args: argparse.Namespace, analysis_dir: Path) -> List[str]:
    """Command for the microglia morphology reviewer, the second Figure 5 HIL server."""
    reviewer = PAPER_ROOT / "scripts" / "shared" / "05_annotate_microglia.py"
    return [
        sys.executable, str(reviewer),
        "--analysis-dir", str(Path(analysis_dir).expanduser().resolve()),
        "--host", str(args.human_review_host),
        "--port", str(int(getattr(args, "microglia_review_port", 0))),
    ]


def _gi_launch_microglia_reviewer(args: argparse.Namespace, analysis_dir: Path,
                                  logger: logging.Logger):
    """Start the morphology reviewer beside the region reviewer, non-blocking.

    The two reviewers are different jobs on the same data: one delineates ARC, ME
    and VMN on whole sections, the other classifies individual Iba1 cells into the
    four morphological states that gate 4 asks for. They are meant to be open at
    the same time, in two browser tabs, so this one runs in the background while
    the region reviewer holds the foreground.

    It needs cell patches, which only exist once an analysis run has extracted
    them, so when the template is absent this says so and starts nothing.
    """
    analysis_dir = Path(analysis_dir).expanduser().resolve()
    template = analysis_dir / "microglia_annotation_template.csv"
    if not template.is_file():
        logger.warning(
            "Microglia morphology reviewer not started: no %s yet. It appears once "
            "an analysis run has extracted the cell patches; --prepare-labels-only "
            "writes it without training or statistics.", template,
        )
        return None
    command = _gi_microglia_review_command(args, analysis_dir)
    logger.warning("Starting microglia morphology reviewer: %s", shlex.join(command))
    try:
        return subprocess.Popen(command)
    except OSError as exc:
        logger.warning("Could not start the microglia morphology reviewer: %s", exc)
        return None


def _gi_prompt_native_review_choice(
    *,
    coverage_complete: bool,
    accepted_n: int,
    pending_n: int,
) -> str:
    """Return the interactive user's choice: reuse, review, or make_figure.

    ``reuse``       keep the accepted segmentations and continue the analysis.
    ``review``      reopen the reviewer to segment sections again.
    ``make_figure`` keep the accepted segmentations and render Figure 5 after
                    the analysis, which is make_figure_5 for this dataset.
    """
    if not sys.stdin.isatty():
        return "reuse" if coverage_complete else "stop"

    reuse_choices = {"", "r", "reuse", "s", "skip"}
    review_choices = {"o", "open", "review", "segment", "again"}
    figure_choices = {"m", "make", "make_figure", "make_figure_5", "figure"}

    if coverage_complete:
        prompt = (
            f"All {accepted_n}/{GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} native ARC/ME "
            "segmentations are accepted. Choose [Enter/r] reuse them and skip review "
            "(default), [o] open review/segment again, or [m] make_figure_5 from these "
            "segmentations: "
        )
        while True:
            try:
                choice = input(prompt).strip().lower()
            except (EOFError, KeyboardInterrupt):
                print(file=sys.stderr)
                return "reuse"
            if choice in reuse_choices:
                return "reuse"
            if choice in review_choices:
                return "review"
            if choice in figure_choices:
                return "make_figure"
            print(
                "Please press Enter (or r) to reuse/skip, o to open review again, "
                "or m to make_figure_5.",
                file=sys.stderr,
            )

    prompt = (
        f"Native ARC/ME segmentation is incomplete: {accepted_n}/"
        f"{GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} accepted, {pending_n} pending. "
        "Choose [y] start/continue image-by-image segmentation, [m] open the reviewer "
        "and use its Make Figure 5 button to accept the machine default segmentations, "
        "or [Enter/n] stop: "
    )
    while True:
        try:
            choice = input(prompt).strip().lower()
        except (EOFError, KeyboardInterrupt):
            print(file=sys.stderr)
            return "stop"
        if choice in {"y", "yes", "start", "continue"}:
            return "review"
        if choice in figure_choices:
            return "review"
        if choice in {"", "n", "no"}:
            return "stop"
        print(
            "Please enter y to start/continue, m to open the reviewer's Make Figure 5 "
            "button, or press Enter (n) to stop.",
            file=sys.stderr,
        )


def _gi_require_complete_human_review(
    sections: Sequence[GISection],
    args: argparse.Namespace,
    work_root: Path,
    coverage_outdir: Path,
    review_analysis_dir: Path,
    logger: logging.Logger,
) -> List[Dict[str, Any]]:
    rows, pending = _gi_audit_manual_region_coverage(sections, args, coverage_outdir, logger)
    command = _gi_human_review_command(args, work_root, review_analysis_dir)
    command_text = shlex.join(command)
    accepted_n = sum(row.get("review_status") == "accepted" for row in rows)
    coverage_complete = (
        not pending
        and len(rows) == GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
        and accepted_n == GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
    )
    if not coverage_complete:
        print(
            f"Pending native HIL ARC/ME/VMN reviews: {len(pending)}; "
            f"accepted {accepted_n} / "
            f"{GFAP_IBA1_EXPECTED_REVIEW_SECTIONS}",
            file=sys.stderr,
        )
        for row in pending:
            print(f"  - {row['sample']}/{row['section']}: {row['reason']}", file=sys.stderr)
    else:
        print(
            f"Complete native ARC/ME coverage found: {accepted_n}/"
            f"{GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} accepted segmentations.",
            file=sys.stderr,
        )
    print(
        "Standalone native review command (stop the reviewer, then rerun the "
        f"Figure 5 analyzer): {command_text}",
        file=sys.stderr,
    )

    choice = _gi_prompt_native_review_choice(
        coverage_complete=coverage_complete,
        accepted_n=accepted_n,
        pending_n=len(pending),
    )
    launch = choice == "review"
    if choice == "make_figure":
        # make_figure_5: keep the current segmentations and render after analysis.
        args.make_figure = True
    if coverage_complete and not launch:
        print(
            "Reusing all previously accepted ARC/ME segmentations; reviewer skipped."
            + (" Figure 5 will be rendered when this analysis finishes."
               if choice == "make_figure" else ""),
            file=sys.stderr,
        )
    elif not coverage_complete and not launch:
        session_kind = "Non-interactive session" if not sys.stdin.isatty() else "Interactive choice"
        print(
            f"{session_kind} with incomplete ARC/ME coverage: reviewer not started; "
            "failing closed before quantification or cleanup.",
            file=sys.stderr,
        )
    if launch:
        action = "re-opening" if coverage_complete else "starting/continuing"
        logger.warning("%s native image-by-image ARC/ME/VMN reviewer: %s", action.capitalize(), command_text)
        # Both Figure 5 reviewers come up together. The morphology one runs in the
        # background because the region reviewer below blocks until it is stopped.
        microglia_proc = _gi_launch_microglia_reviewer(args, review_analysis_dir, logger)
        try:
            completed = subprocess.run(command, check=False)
            if completed.returncode:
                logger.warning("Native reviewer exited with status %d; revalidating receipts.", completed.returncode)
        except KeyboardInterrupt:
            logger.warning("Reviewer stopped; revalidating all accepted masks.")
        finally:
            if microglia_proc is not None and microglia_proc.poll() is None:
                logger.warning("Stopping the microglia morphology reviewer.")
                microglia_proc.terminate()
                try:
                    microglia_proc.wait(timeout=20)
                except subprocess.TimeoutExpired:  # pragma: no cover
                    microglia_proc.kill()
        rows, pending = _gi_audit_manual_region_coverage(sections, args, coverage_outdir, logger)
        accepted_n = sum(row.get("review_status") == "accepted" for row in rows)
        coverage_complete = (
            not pending
            and len(rows) == GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
            and accepted_n == GFAP_IBA1_EXPECTED_REVIEW_SECTIONS
        )
        if not coverage_complete:
            print(
                f"HIL review remains incomplete: {len(pending)} pending; "
                f"accepted {accepted_n} / "
                f"{GFAP_IBA1_EXPECTED_REVIEW_SECTIONS}.",
                file=sys.stderr,
            )
            for row in pending:
                print(f"  - {row['sample']}/{row['section']}: {row['reason']}", file=sys.stderr)
    if not coverage_complete:
        raise ManualRegionMaskError(
            "native human-in-the-loop ARC/ME coverage is incomplete or invalid; "
            "quantification and output cleanup were not started. "
            f"Complete all {GFAP_IBA1_EXPECTED_REVIEW_SECTIONS} images with the "
            f"reviewer, stop it, and rerun this analyzer: {command_text}"
        )
    return rows


def process_gfap_iba1_section(s: GISection, args: argparse.Namespace, outdir: Path, logger: logging.Logger) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], Dict[str, Any]]:
    dapi_lbl, iba_lbl, errors = _gi_load_section_labels(s)
    base, is_rep, rep_lab = strip_repetition_suffix(s.sample)
    source_sample = str(s.sample)
    base_row: Dict[str, Any] = {
        "sample": s.sample, "animal_id": base, "animal_key": canonical_sample_id(base),
        "is_repetition": is_rep, "repetition_label": rep_lab,
        "section_index": int(s.section_index), "section": s.section,
        "status": "ok" if not errors else ";".join(errors),
    }
    if errors or dapi_lbl is None or iba_lbl is None:
        return [], [], base_row
    if s.images.get("GFAP") is None or s.images.get("Iba1") is None or s.images.get("DAPI") is None:
        base_row["status"] = "missing_required_image"
        return [], [], base_row

    iba1_segmentation_valid = bool(int(iba_lbl.max()) > 0 or args.accept_empty_iba1_seg)
    base_row["iba1_segmentation_valid"] = iba1_segmentation_valid
    base_row["iba1_segmentation_empty"] = bool(int(iba_lbl.max()) == 0)
    if not iba1_segmentation_valid:
        logger.warning(
            "Iba1 mask has zero objects; GFAP remains valid but microglial area is "
            "excluded from density denominators: %s %s", s.sample, s.section
        )

    dapi_signal = _gi_read_signal(s.images["DAPI"], "DAPI")
    gfap_signal = _gi_read_signal(s.images["GFAP"], "GFAP")
    iba_signal = _gi_read_signal(s.images["Iba1"], "Iba1")
    if dapi_lbl.shape != gfap_signal.shape or dapi_lbl.shape != iba_signal.shape:
        base_row["status"] = "shape_mismatch_images"
        return [], [], base_row

    region_mask, manual_region_info = _gi_load_manual_region_mask(
        s, dapi_lbl, args, logger
    )
    if region_mask is None:
        raise ManualRegionMaskError(
            f"native HIL ARC/ME/VMN mask unexpectedly unavailable for {s.sample} {s.section}; "
            "automated inference is disabled"
        )
    region_info = manual_region_info or {}
    base_row.update({
        k: v for k, v in region_info.items()
        if isinstance(v, (str, int, float, bool)) or v is None
    })
    gfap_pos, gfap_corr, gfap_by_region = _gi_quantify_gfap_region(gfap_signal, region_mask, args)
    iba_norm = _gi_robust_normalize(iba_signal, args)
    iba_corr, _iba_bg = _gi_background_correct(iba_norm, um_to_px(float(args.iba1_background_radius_um), args, min_px=3))
    iba_corr_raw, _iba_bg_raw = _gi_background_correct(
        np.asarray(iba_signal, dtype=np.float32),
        um_to_px(float(args.iba1_background_radius_um), args, min_px=3),
    )
    iba_fg_thr = _gi_robust_threshold(iba_corr[region_mask > 0], factor=float(args.iba1_threshold_factor), k=2.5)
    iba_fg = (iba_corr > iba_fg_thr) & (region_mask > 0)

    # Preserve the exact source folder in all per-section products. Repetition
    # folders are merged only when animal-level values are aggregated; using
    # ``base`` here would overwrite FR6-2/S## with FR6-2(R)/S## and reduce the
    # required 61 independent section geometries.
    processed_dir = outdir / args.reconstructed_dirname / source_sample / s.section
    ensure_dir(processed_dir)
    if bool(args.save_reconstructed_tifs):
        tifffile.imwrite(str(processed_dir / f"{source_sample}_{s.section}_DAPI_labels.tif"), np.asarray(dapi_lbl, dtype=np.uint32))
        region_output_path = (processed_dir / f"{source_sample}_{s.section}_ARC_ME_labels.tif").resolve()
        native_proposal_text = str(region_info.get("manual_region_proposal_mask_path") or "").strip()
        preserve_native_proposal = bool(
            region_info.get("manual_region_accepted_source") == "automated_mask_unchanged"
            and region_info.get("manual_region_proposal_mask_sha256_verified") is True
            and native_proposal_text
            and Path(native_proposal_text).expanduser().resolve() == region_output_path
        )
        if preserve_native_proposal:
            logger.warning("Preserving verified native proposal provenance file without overwriting: %s", region_output_path)
        else:
            tifffile.imwrite(str(region_output_path), np.asarray(region_mask, dtype=np.uint8))
        tifffile.imwrite(str(processed_dir / f"{source_sample}_{s.section}_GFAP_positive.tif"), np.asarray(gfap_pos, dtype=np.uint8))

    # Save wall/floor geometry in a compact NPZ for exact QC regeneration.
    geom_arrays = {k: np.asarray(v) for k, v in region_info.items() if str(k).startswith("_") and isinstance(v, np.ndarray)}
    np.savez_compressed(processed_dir / "region_geometry_arrays.npz", **geom_arrays)
    scalar_info = {k: v for k, v in region_info.items() if not str(k).startswith("_") and isinstance(v, (str, int, float, bool, type(None), np.integer, np.floating))}
    (processed_dir / "region_geometry.json").write_text(json.dumps(scalar_info, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o)), encoding="utf-8")

    pixel_area = float(args.um_per_px) ** 2
    region_rows: List[Dict[str, Any]] = []
    drawn_regions = [
        name for name in GI_MANUAL_REGIONS
        if region_info.get(f"manual_region_{name.lower()}_status") == "drawn"
    ]
    for region in drawn_regions:
        rm = region_mask == REGION_NAME_TO_CODE[region]
        area_px = int(np.count_nonzero(rm))
        dapi_n = 0
        for rp in regionprops(dapi_lbl):
            y, x = rp.centroid
            if _gi_region_for_point(region_info, region_mask, float(y), float(x)) == region:
                dapi_n += 1
        rr: Dict[str, Any] = {
            **base_row, "region": region, "region_area_px": area_px,
            "region_area_um2": float(area_px) * pixel_area, "dapi_nuclei": int(dapi_n),
            "dapi_density_per_um2": float(dapi_n / (area_px * pixel_area)) if area_px else float("nan"),
            "iba1_segmentation_valid": iba1_segmentation_valid,
            "iba1_assay_area_um2": float(area_px) * pixel_area if iba1_segmentation_valid else 0.0,
            "iba1_positive_area_px": int(np.count_nonzero(iba_fg & rm)),
            "iba1_area_fraction": float(np.count_nonzero(iba_fg & rm) / area_px) if area_px else float("nan"),
            "iba1_corrected_integrated_intensity": float(np.sum(iba_corr_raw[rm])) if area_px else 0.0,
            "iba1_corrected_mean_intensity": float(np.mean(iba_corr_raw[rm])) if area_px else float("nan"),
            "iba1_intensity_scale": "raw_channel_background_corrected_no_percentile_rescaling",
        }
        rr.update(gfap_by_region[region])
        rr["region_area_um2"] = float(area_px) * pixel_area
        rr["gfap_positive_area_um2"] = float(rr["gfap_positive_area_px"]) * pixel_area
        rr["gfap_corrected_intensity_per_um2"] = float(rr["gfap_corrected_integrated_intensity"] / rr["region_area_um2"]) if rr["region_area_um2"] > 0 else float("nan")
        rr["iba1_corrected_intensity_per_um2"] = float(rr["iba1_corrected_integrated_intensity"] / rr["region_area_um2"]) if rr["region_area_um2"] > 0 else float("nan")
        region_rows.append(rr)

    patch_size = int(args.cnn_patch_size)
    patch_root = outdir / "microglia_patches" / source_sample / s.section
    cell_rows: List[Dict[str, Any]] = []
    accepted_iba_lbl = np.zeros_like(iba_lbl, dtype=np.int32)
    rejected_iba_lbl = np.zeros_like(iba_lbl, dtype=np.int32)
    iba_props = list(regionprops(iba_lbl, intensity_image=iba_corr))
    dapi_props_by_label = {int(p.label): p for p in regionprops(dapi_lbl)}
    logger.info("Extracting morphology/QC for %s %s: %d raw Iba1 objects", s.sample, s.section, len(iba_props))
    for rp in iba_props:
        lab0 = int(rp.label)
        mask_global = iba_lbl == lab0
        object_qc = _gi_evaluate_microglia_object_qc(
            rp, mask_global, dapi_lbl, iba_fg, region_mask, args
        )
        dapi_lab = int(object_qc["dapi_nucleus_label"])
        nuc_prop = dapi_props_by_label.get(dapi_lab)
        nuc_centroid = nuc_prop.centroid if nuc_prop is not None else rp.centroid
        y, x = float(nuc_centroid[0]), float(nuc_centroid[1])
        region = _gi_region_for_point(region_info, region_mask, y, x)
        region_dl_vote = region_name_from_mask(region_mask, y, x)
        try:
            region_wall_vote = region_name_from_projected_walls(region_info, y, x)
        except Exception:
            region_wall_vote = "UNKNOWN"
        if region_dl_vote not in set(GI_MANUAL_REGIONS):
            region_dl_vote = "OUTSIDE_ARC_ME"
        if region not in set(GI_MANUAL_REGIONS):
            existing_reasons = [
                value for value in str(object_qc["cell_qc_rejection_reasons"]).split(";") if value
            ]
            if "outside_final_arc_me_mask" not in existing_reasons:
                existing_reasons.append("outside_final_arc_me_mask")
            object_qc["cell_qc_pass"] = False
            object_qc["cell_qc_rejection_reasons"] = ";".join(existing_reasons)
            object_qc["qc_reject_outside_final_arc_me_mask"] = True
        per = float(rp.perimeter) if float(rp.perimeter) > 0 else 1.0
        circ = float(4.0 * math.pi * float(rp.area) / (per * per))
        from scipy import ndimage as ndi
        dist_local = ndi.distance_transform_edt(rp.image)
        dmax = float(np.max(dist_local)) if dist_local.size else 0.0
        soma = rp.image & (dist_local >= max(1.0, 0.50 * dmax)) if dmax > 0 else np.zeros_like(rp.image)
        soma_fraction = float(np.count_nonzero(soma) / float(rp.area)) if rp.area else 0.0
        sk_len = float(object_qc["skeleton_length_px"])
        cell_uid = f"{canonical_sample_id(source_sample)}__{s.section.lower()}__iba1_{lab0:05d}"
        patch_path = patch_root / f"{cell_uid}.png"
        patch_path_value = ""
        if bool(object_qc["cell_qc_pass"]):
            accepted_iba_lbl[mask_global] = lab0
        else:
            rejected_iba_lbl[mask_global] = lab0
        if bool(object_qc["cell_qc_pass"]) and bool(args.save_microglia_patches):
            patch = _gi_extract_patch(iba_corr, mask_global, rp.centroid, patch_size)
            _gi_save_patch(patch_path, patch)
            patch_path_value = str(patch_path)
        row: Dict[str, Any] = {
            **base_row, "cell_uid": cell_uid, "iba1_label": lab0, "dapi_nucleus_label": int(dapi_lab),
            "centroid_y": y, "centroid_x": x, "region": region,
            "region_assignment_method": (
                "manual_arc_me_vmn_mask_centroid_consistent_with_area_denominator"
                if region_info.get("region_source") == "manual_region_tiff"
                else "automated_arc_me_mask_centroid_consistent_with_area_denominator"
            ),
            "region_assignment_source": str(region_info.get("region_source", "unknown")),
            "region_mask_vote": region_dl_vote,
            "region_dl_vote": region_dl_vote,
            "region_projected_wall_vote": region_wall_vote,
            "region_dl_agrees_with_anatomy": bool(region_dl_vote == region_wall_vote),
            "patch_path": patch_path_value,
            "area_px": float(rp.area), "area_um2": float(rp.area) * pixel_area,
            "perimeter_px": per, "perimeter_um": per * float(args.um_per_px),
            "circularity": circ, "solidity": float(rp.solidity), "eccentricity": float(rp.eccentricity),
            "extent": float(rp.extent), "major_axis_length_px": float(rp.major_axis_length),
            "minor_axis_length_px": float(rp.minor_axis_length), "equivalent_diameter_px": float(rp.equivalent_diameter_area),
            "convex_area_px": float(rp.convex_area), "soma_fraction": soma_fraction,
            "process_fraction": float(1.0 - soma_fraction),
            "mean_process_width_px": float(rp.area / max(sk_len, 1.0)),
            "iba1_mean_intensity": float(rp.mean_intensity),
            "iba1_max_intensity": float(np.max(iba_corr[mask_global])) if np.any(mask_global) else 0.0,
            "iba1_integrated_intensity": float(np.sum(iba_corr[mask_global])),
            "fractal_dimension": _gi_fractal_dimension(rp.image),
        }
        row.update(object_qc)
        row["skeleton_length_um"] = row["skeleton_length_px"] * float(args.um_per_px)
        row["skeleton_length_norm"] = row["skeleton_length_px"] / math.sqrt(max(row["area_px"], 1.0))
        row["branchpoints_per_100_skeleton_px"] = 100.0 * row["skeleton_branchpoints"] / max(row["skeleton_length_px"], 1.0)
        row["endpoints_per_100_skeleton_px"] = 100.0 * row["skeleton_endpoints"] / max(row["skeleton_length_px"], 1.0)
        cell_rows.append(row)

    if bool(args.save_reconstructed_tifs):
        tifffile.imwrite(str(processed_dir / f"{source_sample}_{s.section}_Iba1_cells_raw.tif"), np.asarray(iba_lbl, dtype=np.uint32))
        tifffile.imwrite(str(processed_dir / f"{source_sample}_{s.section}_Iba1_cells.tif"), np.asarray(accepted_iba_lbl, dtype=np.uint32))
        tifffile.imwrite(str(processed_dir / f"{source_sample}_{s.section}_Iba1_cells_rejected.tif"), np.asarray(rejected_iba_lbl, dtype=np.uint32))

    accepted_n = int(sum(bool(r["cell_qc_pass"]) for r in cell_rows))
    rejected_n = int(len(cell_rows) - accepted_n)
    for rr in region_rows:
        region_candidates = [r for r in cell_rows if r["region"] == rr["region"]]
        region_accepted = int(sum(bool(r["cell_qc_pass"]) for r in region_candidates))
        rr["iba1_objects_raw"] = int(len(region_candidates))
        rr["microglia_cells_qc_accepted"] = region_accepted
        rr["iba1_objects_qc_rejected"] = int(len(region_candidates) - region_accepted)

    base_row.update({
        "dapi_nuclei": int(len(dapi_props_by_label)), "iba1_cells": accepted_n,
        "iba1_objects_raw": int(len(cell_rows)), "microglia_cells_qc_accepted": accepted_n,
        "iba1_objects_qc_rejected": rejected_n, "iba1_labels_max": int(iba_lbl.max()),
        "image_height_px": int(dapi_lbl.shape[0]), "image_width_px": int(dapi_lbl.shape[1]),
        "um_per_px": float(args.um_per_px), "pixel_area_um2": pixel_area,
    })
    for k, v in region_info.items():
        if isinstance(v, (str, int, float, bool)) or v is None:
            base_row[k] = v
    return cell_rows, region_rows, base_row


def _gi_pseudo_labels(cell_df: Any, args: argparse.Namespace, logger: logging.Logger) -> Tuple[Any, Any, Dict[str, Any]]:
    """Create reproducible morphology clusters and map them to the four requested biological states."""
    _gi_require_qc_accepted(cell_df, "morphometric pseudo-label KMeans")
    from sklearn.cluster import KMeans
    from sklearn.preprocessing import StandardScaler
    features = [
        "circularity", "solidity", "soma_fraction", "mean_process_width_px",
        "skeleton_length_norm", "branchpoints_per_100_skeleton_px",
        "endpoints_per_100_skeleton_px", "sholl_max_crossings", "fractal_dimension",
        "eccentricity", "area_um2",
    ]
    xdf = cell_df[features].apply(pd.to_numeric, errors="coerce").copy()
    for c in features:
        med = float(xdf[c].median()) if xdf[c].notna().any() else 0.0
        xdf[c] = xdf[c].fillna(med)
    scaler = StandardScaler()
    X = scaler.fit_transform(xdf.to_numpy(dtype=float))
    k = min(4, max(1, len(cell_df)))
    km = KMeans(n_clusters=k, n_init=30, random_state=int(args.cnn_seed))
    cluster = km.fit_predict(X)
    dist = km.transform(X)
    order_dist = np.sort(dist, axis=1)
    conf = 1.0 - order_dist[:, 0] / np.maximum(order_dist[:, 1] if order_dist.shape[1] > 1 else order_dist[:, 0] + 1.0, 1e-9)
    zdf = pd.DataFrame(X, columns=features)
    zdf["cluster"] = cluster
    centers = zdf.groupby("cluster").mean()
    activation = (
        1.35 * centers["circularity"] + 1.10 * centers["solidity"] + 1.00 * centers["soma_fraction"]
        + 0.70 * centers["mean_process_width_px"] - 1.20 * centers["skeleton_length_norm"]
        - 0.90 * centers["branchpoints_per_100_skeleton_px"] - 0.45 * centers["sholl_max_crossings"]
    )
    mapping: Dict[int, str] = {}
    remaining = list(centers.index.astype(int))
    if remaining:
        amo = int(activation.idxmax()); mapping[amo] = "Amoeboid"; remaining.remove(amo)
    if remaining:
        # Rod-like microglia are elongated and bipolar: a stretched soma with few
        # processes. That is close to the opposite of the branching-dominated
        # score that used to sit here to pick Rod-like, so the rule was
        # rewritten with the nomenclature rather than merely relabelled.
        # eccentricity is already one of the clustering features.
        rod_score = (
            centers.loc[remaining, "eccentricity"]
            - 0.60 * centers.loc[remaining, "branchpoints_per_100_skeleton_px"]
            - 0.40 * centers.loc[remaining, "sholl_max_crossings"]
        )
        rod = int(rod_score.idxmax()); mapping[rod] = "Rod-like"; remaining.remove(rod)
    if remaining:
        ram_score = centers.loc[remaining, "skeleton_length_norm"] + centers.loc[remaining, "endpoints_per_100_skeleton_px"] - 0.4 * centers.loc[remaining, "circularity"]
        ram = int(ram_score.idxmax()); mapping[ram] = "Ramified"; remaining.remove(ram)
    for c in remaining:
        mapping[int(c)] = "Activated"
    # Degenerate datasets with fewer than four clusters use the nearest semantic state mapping above.
    labels_out = np.asarray([mapping.get(int(c), "Activated") for c in cluster], dtype=object)
    info = {"features": features, "cluster_to_state": mapping, "cluster_centers_z": centers.to_dict(orient="index")}
    logger.warning("Morphometry cluster mapping: %s", mapping)
    return labels_out, conf.astype(float), info


def _gi_load_manual_labels(path: Optional[str], cell_df: Any, logger: logging.Logger) -> Dict[str, str]:
    if not path:
        return {}
    p = Path(path).resolve()
    if not p.exists():
        raise FileNotFoundError(p)
    lab = pd.read_csv(p)
    if "state" not in lab.columns:
        raise ValueError("--microglia-labels-csv requires a state column")
    allowed = {norm_token(s): s for s in MICROGLIA_STATES}
    valid_uids = set(cell_df["cell_uid"].astype(str).tolist())
    out: Dict[str, str] = {}
    ignored = 0
    if "cell_uid" in lab.columns:
        for _, r in lab.iterrows():
            key = str(r.get("cell_uid", "")).strip(); st = allowed.get(norm_token(r.get("state", "")))
            if key and st and key in valid_uids:
                out[key] = st
            elif key and st:
                ignored += 1
    elif "patch_path" in lab.columns:
        patch_to_uid = {str(Path(p).resolve()): uid for p, uid in zip(cell_df["patch_path"], cell_df["cell_uid"]) if str(p)}
        for _, r in lab.iterrows():
            pp = str(Path(str(r.get("patch_path", ""))).resolve()); st = allowed.get(norm_token(r.get("state", "")))
            if pp in patch_to_uid and st:
                out[patch_to_uid[pp]] = st
            elif st:
                ignored += 1
    else:
        raise ValueError("Manual label CSV requires cell_uid or patch_path")
    logger.warning("Loaded reviewed microglia labels: %d", len(out))
    if ignored:
        logger.warning("Ignored %d reviewed labels absent from the biological-cell QC accepted set.", ignored)
    return out


def _gi_train_and_predict_cnn(cell_df: Any, pseudo_labels: Any, pseudo_conf: Any, manual: Dict[str, str], args: argparse.Namespace, outdir: Path, logger: logging.Logger) -> Tuple[Any, Dict[str, Any]]:
    """Train a compact CNN on reviewed labels when present, otherwise on high-confidence morphology pseudo-labels."""
    _gi_require_qc_accepted(cell_df, "microglia CNN training/inference")
    import torch
    import torch.nn as nn
    from torch.utils.data import Dataset, DataLoader
    from PIL import Image
    from sklearn.model_selection import GroupShuffleSplit
    from sklearn.metrics import confusion_matrix, balanced_accuracy_score

    state_to_idx = MICROGLIA_STATE_TO_CODE
    idx_to_state = {v: k for k, v in state_to_idx.items()}
    device = _gi_resolve_torch_device(args.cnn_device, logger)
    rng = np.random.default_rng(int(args.cnn_seed))
    torch.manual_seed(int(args.cnn_seed))

    labels = np.asarray(pseudo_labels, dtype=object).copy()
    sources = np.asarray(["pseudo"] * len(cell_df), dtype=object)
    train_mask = np.asarray(pseudo_conf, dtype=float) >= float(args.cnn_pseudo_confidence)
    for i, uid in enumerate(cell_df["cell_uid"].astype(str).tolist()):
        if uid in manual:
            labels[i] = manual[uid]
            sources[i] = "manual"
            train_mask[i] = True
    if str(args.microglia_classifier_mode) == "supervised":
        train_mask = sources == "manual"
    qc_ok = cell_df["cell_qc_pass"].fillna(False).astype(bool).to_numpy()
    valid_path = cell_df["patch_path"].astype(str).map(lambda x: bool(x) and Path(x).exists()).to_numpy()
    valid_path &= qc_ok
    train_mask &= valid_path & qc_ok

    class_counts = {s: int(np.sum(train_mask & (labels == s))) for s in MICROGLIA_STATES}
    usable_classes = [s for s, n in class_counts.items() if n >= int(args.cnn_min_class_cells)]
    training_animals = sorted(
        set(cell_df.loc[train_mask, "animal_key"].astype(str).tolist())
    )
    enough = (
        int(np.sum(train_mask)) >= int(args.cnn_min_training_cells)
        and len(usable_classes) >= 2
        and len(training_animals) >= int(args.cnn_min_training_animals)
    )
    logger.warning("CNN training candidates: %d | class counts: %s", int(np.sum(train_mask)), class_counts)
    logger.warning("CNN training animals: %d | %s", len(training_animals), training_animals)
    if str(args.microglia_classifier_mode) == "supervised" and not enough:
        raise ValueError(
            "Reviewed supervised labels do not meet the minimum training design: "
            f"candidates={int(np.sum(train_mask))}, classes={class_counts}, "
            f"animals={len(training_animals)}. Add reviewed cells/animals or use "
            "--microglia-classifier-mode hybrid for explicitly exploratory pseudo-labels."
        )

    # Defaults to morphometric result when CNN cannot be trained safely.
    result = pd.DataFrame({
        "cell_uid": cell_df["cell_uid"].astype(str),
        "microglia_state": labels,
        "classifier_source": "morphometry",
        "classifier_confidence": np.asarray(pseudo_conf, dtype=float),
    })
    for s in MICROGLIA_STATES:
        result[f"prob_{norm_token(s)}"] = (labels == s).astype(float)

    model_path = Path(args.cnn_model_out).resolve() if args.cnn_model_out else outdir / "models" / "microglia_morphology_cnn.pt"
    ensure_dir(model_path.parent)

    class PatchDataset(Dataset):
        def __init__(self, indices: Any, training: bool, target_labels: Any):
            self.indices = np.asarray(indices, dtype=int)
            self.training = training
            self.target_labels = target_labels
        def __len__(self): return len(self.indices)
        def __getitem__(self, j):
            idx = int(self.indices[j])
            a = np.asarray(Image.open(cell_df.iloc[idx]["patch_path"]).convert("RGB"), dtype=np.float32) / 255.0
            if self.training:
                if rng.random() < 0.5: a = np.flip(a, axis=1).copy()
                if rng.random() < 0.5: a = np.flip(a, axis=0).copy()
                krot = int(rng.integers(0, 4)); a = np.rot90(a, krot).copy()
                a[..., 0] = np.clip(a[..., 0] * float(rng.uniform(0.85, 1.15)), 0, 1)
            t = torch.from_numpy(np.transpose(a, (2, 0, 1)).copy()).float()
            y = int(state_to_idx[str(self.target_labels[idx])])
            return t, y, idx

    class TinyMorphCNN(nn.Module):
        def __init__(self):
            super().__init__()
            self.features = nn.Sequential(
                nn.Conv2d(3, 16, 3, padding=1), nn.BatchNorm2d(16), nn.ReLU(inplace=True), nn.MaxPool2d(2),
                nn.Conv2d(16, 32, 3, padding=1), nn.BatchNorm2d(32), nn.ReLU(inplace=True), nn.MaxPool2d(2),
                nn.Conv2d(32, 64, 3, padding=1), nn.BatchNorm2d(64), nn.ReLU(inplace=True), nn.MaxPool2d(2),
                nn.Conv2d(64, 96, 3, padding=1), nn.BatchNorm2d(96), nn.ReLU(inplace=True),
                nn.AdaptiveAvgPool2d((1, 1)),
            )
            self.head = nn.Sequential(nn.Flatten(), nn.Dropout(0.25), nn.Linear(96, 4))
        def forward(self, x): return self.head(self.features(x))

    model = TinyMorphCNN().to(device)
    history: List[Dict[str, Any]] = []
    loaded = False
    best_epoch = None
    calibrate_batchnorm = bool(getattr(args, "cnn_calibrate_batchnorm", False))
    if args.cnn_model_in:
        ckpt = torch.load(str(Path(args.cnn_model_in).resolve()), map_location=device)
        model.load_state_dict(ckpt["model_state"] if isinstance(ckpt, dict) and "model_state" in ckpt else ckpt)
        loaded = True
        logger.warning("Loaded existing microglia CNN: %s", args.cnn_model_in)
    elif enough and str(args.microglia_classifier_mode) in {"hybrid", "cnn", "supervised"}:
        idx_all = np.flatnonzero(train_mask)
        groups = cell_df.iloc[idx_all]["animal_key"].astype(str).to_numpy()
        splitter = GroupShuffleSplit(
            n_splits=1, test_size=0.20, random_state=int(args.cnn_seed)
        )
        tr_rel, va_rel = next(splitter.split(idx_all, labels[idx_all], groups))
        tr_idx, va_idx = idx_all[tr_rel], idx_all[va_rel]
        train_ds = PatchDataset(tr_idx, True, labels); val_ds = PatchDataset(va_idx, False, labels)
        train_loader = DataLoader(train_ds, batch_size=int(args.cnn_batch_size), shuffle=True, num_workers=0)
        val_loader = DataLoader(val_ds, batch_size=int(args.cnn_batch_size), shuffle=False, num_workers=0)
        if calibrate_batchnorm:
            from torch.optim.swa_utils import update_bn
            # Only fit-animal cells enter normalization calibration. Use the
            # inference representation, with no flips or intensity augmentation.
            calibration_loader = DataLoader(
                PatchDataset(tr_idx, False, labels), batch_size=int(args.cnn_batch_size),
                shuffle=False, num_workers=0,
            )
        counts = np.bincount([state_to_idx[str(labels[i])] for i in tr_idx], minlength=4).astype(float)
        weights = np.where(counts > 0, np.sum(counts) / (4.0 * np.maximum(counts, 1)), 0.0)
        criterion = nn.CrossEntropyLoss(weight=torch.tensor(weights, dtype=torch.float32, device=device))
        optim = torch.optim.AdamW(model.parameters(), lr=float(args.cnn_learning_rate), weight_decay=float(args.cnn_weight_decay))
        best_loss = float("inf"); best_state = None; bad = 0
        for epoch in range(1, int(args.cnn_epochs) + 1):
            model.train(); train_loss = 0.0; train_n = 0
            for xb, yb, _ in train_loader:
                xb, yb = xb.to(device), yb.to(device)
                optim.zero_grad(set_to_none=True); logits = model(xb); loss = criterion(logits, yb); loss.backward(); optim.step()
                train_loss += float(loss.item()) * len(yb); train_n += len(yb)
            if calibrate_batchnorm:
                # Preserve optimization RNG: calibration must not change the
                # subsequent training shuffles, augmentation, or dropout draws.
                cuda_devices = [] if str(device) == "cpu" else [torch.device(device).index or torch.cuda.current_device()]
                with torch.random.fork_rng(devices=cuda_devices):
                    update_bn(calibration_loader, model, device=device)
            model.eval(); val_loss = 0.0; val_n = 0; yp: List[int] = []; yt: List[int] = []
            with torch.no_grad():
                for xb, yb, _ in val_loader:
                    xb, yb = xb.to(device), yb.to(device); logits = model(xb); loss = criterion(logits, yb)
                    val_loss += float(loss.item()) * len(yb); val_n += len(yb)
                    yp.extend(torch.argmax(logits, 1).cpu().tolist()); yt.extend(yb.cpu().tolist())
            trl = train_loss / max(train_n, 1); val = val_loss / max(val_n, 1)
            bal = float(balanced_accuracy_score(yt, yp)) if yt and len(set(yt)) > 1 else float("nan")
            history.append({"epoch": epoch, "train_loss": trl, "val_loss": val, "val_balanced_accuracy": bal})
            logger.info("CNN epoch %d/%d train_loss=%.5f val_loss=%.5f val_bal_acc=%s", epoch, int(args.cnn_epochs), trl, val, f"{bal:.3f}" if np.isfinite(bal) else "NA")
            if val < best_loss - 1e-5:
                best_loss = val; best_epoch = epoch; best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}; bad = 0
            else:
                bad += 1
                if bad >= int(args.cnn_patience):
                    logger.info("CNN early stopping at epoch %d", epoch); break
        if best_state is not None:
            model.load_state_dict(best_state)
        torch.save({
            "model_state": model.state_dict(), "states": list(MICROGLIA_STATES), "patch_size": int(args.cnn_patch_size),
            "training_source": "manual" if manual else "pseudo_morphometry", "class_counts": class_counts,
            "best_epoch": best_epoch,
            "batchnorm_calibration": "training_cells_before_validation" if calibrate_batchnorm else "none",
            "version": GFAP_IBA1_VERSION,
        }, model_path)
        loaded = True
        pd.DataFrame(history).to_csv(outdir / "microglia_cnn_training_history.csv", index=False)
        if yt:
            cm = confusion_matrix(yt, yp, labels=list(range(4)))
            pd.DataFrame(cm, index=MICROGLIA_STATES, columns=MICROGLIA_STATES).to_csv(outdir / "microglia_cnn_validation_confusion_matrix.csv")
        logger.warning("Saved microglia morphology CNN: %s", model_path)
    else:
        logger.warning(
            "CNN training skipped; morphometric classes retained. Need >=%d cells, "
            ">=2 classes with >=%d cells, and >=%d animals for animal-blocked validation.",
            int(args.cnn_min_training_cells), int(args.cnn_min_class_cells),
            int(args.cnn_min_training_animals),
        )

    if loaded and str(args.microglia_classifier_mode) != "morphometry":
        all_idx = np.flatnonzero(valid_path & qc_ok)
        infer_ds = PatchDataset(all_idx, False, labels)
        infer_loader = DataLoader(infer_ds, batch_size=int(args.cnn_batch_size), shuffle=False, num_workers=0)
        probs_all = np.zeros((len(cell_df), 4), dtype=float)
        model.eval()
        with torch.no_grad():
            for xb, _yb, idxb in infer_loader:
                prob = torch.softmax(model(xb.to(device)), dim=1).cpu().numpy()
                probs_all[np.asarray(idxb, dtype=int)] = prob
        pred = np.argmax(probs_all, axis=1)
        conf = np.max(probs_all, axis=1)
        result["microglia_state"] = [idx_to_state[int(i)] for i in pred]
        result["classifier_source"] = "cnn_manual" if manual else "cnn_pseudo_morphometry"
        result["classifier_confidence"] = conf
        result["classifier_uncertain"] = conf < float(args.cnn_prediction_confidence)
        for s, idx in state_to_idx.items():
            result[f"prob_{norm_token(s)}"] = probs_all[:, idx]
    else:
        result["classifier_uncertain"] = result["classifier_confidence"] < float(args.cnn_prediction_confidence)

    meta = {
        "device": device, "model_path": str(model_path) if loaded else "", "cnn_used": bool(loaded),
        "manual_labels_n": len(manual), "training_candidates_n": int(np.sum(train_mask)),
        "class_counts": class_counts, "training_animals": training_animals,
        "best_epoch": best_epoch,
        "batchnorm_calibration": "training_cells_before_validation" if calibrate_batchnorm else "none",
        "validation_split": "whole_animal_blocked",
        "metric_interpretation": (
            "agreement with reviewed labels when supplied; otherwise whole-animal "
            "held-out morphology pseudo-label agreement, not biological truth"
        ),
    }
    return result, meta


def classify_microglia(cell_df: Any, args: argparse.Namespace, outdir: Path, logger: logging.Logger) -> Tuple[Any, Dict[str, Any]]:
    if cell_df.empty:
        return cell_df, {}
    _gi_require_qc_accepted(cell_df, "microglia classification")
    pseudo, pconf, cluster_info = _gi_pseudo_labels(cell_df, args, logger)
    manual = _gi_load_manual_labels(args.microglia_labels_csv, cell_df, logger)
    if str(args.microglia_classifier_mode) == "supervised" and not manual:
        raise ValueError(
            "Supervised microglia classification requires reviewed non-empty labels in "
            "--microglia-labels-csv. The blank annotation template is not a trained model."
        )
    template = cell_df[["cell_uid", "sample", "animal_id", "section", "region", "patch_path"]].copy()
    template["state"] = ""
    template["pseudo_state"] = pseudo
    template["pseudo_confidence"] = pconf
    template.to_csv(outdir / "microglia_annotation_template.csv", index=False)
    (outdir / "microglia_pseudo_cluster_mapping.json").write_text(json.dumps(cluster_info, indent=2, default=str), encoding="utf-8")
    if bool(args.prepare_labels_only):
        out = cell_df.copy(); out["microglia_state"] = pseudo; out["classifier_confidence"] = pconf; out["classifier_source"] = "morphometry_annotation_preparation"; out["classifier_uncertain"] = pconf < float(args.cnn_prediction_confidence)
        return out, {"prepare_labels_only": True}
    if str(args.microglia_classifier_mode) == "morphometry":
        pred = pd.DataFrame({"cell_uid": cell_df["cell_uid"], "microglia_state": pseudo, "classifier_source": "morphometry", "classifier_confidence": pconf, "classifier_uncertain": pconf < float(args.cnn_prediction_confidence)})
        meta = {"cnn_used": False, "cluster_info": cluster_info}
    else:
        pred, meta = _gi_train_and_predict_cnn(cell_df, pseudo, pconf, manual, args, outdir, logger)
        meta["cluster_info"] = cluster_info
    out = cell_df.merge(pred, on="cell_uid", how="left")
    return out, meta


def _gi_aggregate_animal_regions(cell_df: Any, section_region_df: Any, samplesheet: Any, args: argparse.Namespace) -> Any:
    _gi_require_qc_accepted(cell_df, "animal/region endpoint aggregation")
    if section_region_df.empty:
        return section_region_df
    sum_cols = [
        "region_area_px", "region_area_um2", "dapi_nuclei", "gfap_positive_area_px", "gfap_positive_area_um2",
        "gfap_corrected_integrated_intensity", "iba1_positive_area_px", "iba1_corrected_integrated_intensity",
        "iba1_assay_area_um2", "iba1_objects_raw", "microglia_cells_qc_accepted",
        "iba1_objects_qc_rejected",
    ]
    base_cols = ["animal_id", "animal_key", "region"]
    sr = section_region_df.copy()
    for c in sum_cols:
        sr[c] = pd.to_numeric(sr.get(c, 0), errors="coerce").fillna(0)
    agg = sr.groupby(base_cols, dropna=False)[sum_cols].sum().reset_index()
    agg["n_sections"] = sr.groupby(base_cols, dropna=False).size().values
    if "iba1_segmentation_valid" in sr.columns:
        agg["n_iba1_valid_sections"] = (
            sr.groupby(base_cols, dropna=False)["iba1_segmentation_valid"]
            .sum().astype(int).values
        )
    # Counts by state from per-cell table.
    if cell_df is not None and not cell_df.empty:
        cc = cell_df.copy()
        cc = cc[cc["region"].isin(GI_MANUAL_REGIONS)]
        total = cc.groupby(base_cols, dropna=False).size().rename("microglia_cells").reset_index()
        agg = agg.merge(total, on=base_cols, how="left")
        for state in MICROGLIA_STATES:
            col = f"microglia_{norm_token(state)}_cells"
            x = cc[cc["microglia_state"] == state].groupby(base_cols, dropna=False).size().rename(col).reset_index()
            agg = agg.merge(x, on=base_cols, how="left")
        uncertain = cc[cc["classifier_uncertain"].fillna(False)].groupby(base_cols, dropna=False).size().rename("microglia_uncertain_cells").reset_index()
        agg = agg.merge(uncertain, on=base_cols, how="left")
    else:
        agg["microglia_cells"] = 0
    count_cols = [c for c in agg.columns if c.startswith("microglia_") and c.endswith("_cells")]
    for c in count_cols + ["microglia_cells"]:
        if c in agg.columns: agg[c] = pd.to_numeric(agg[c], errors="coerce").fillna(0).astype(int)

    area = pd.to_numeric(agg["region_area_um2"], errors="coerce")
    agg["dapi_density_per_um2"] = agg["dapi_nuclei"] / area.replace(0, np.nan)
    agg["gfap_area_fraction"] = agg["gfap_positive_area_px"] / agg["region_area_px"].replace(0, np.nan)
    agg["gfap_corrected_intensity_per_um2"] = agg["gfap_corrected_integrated_intensity"] / area.replace(0, np.nan)
    agg["gfap_corrected_mean_intensity"] = agg["gfap_corrected_integrated_intensity"] / agg["region_area_px"].replace(0, np.nan)
    agg["gfap_reactivity_index"] = agg["gfap_area_fraction"] * agg["gfap_corrected_mean_intensity"]
    agg["iba1_area_fraction"] = agg["iba1_positive_area_px"] / agg["region_area_px"].replace(0, np.nan)
    agg["iba1_corrected_intensity_per_um2"] = agg["iba1_corrected_integrated_intensity"] / area.replace(0, np.nan)
    iba1_assay_area = pd.to_numeric(agg["iba1_assay_area_um2"], errors="coerce")
    agg["microglia_density_per_um2"] = agg["microglia_cells"] / iba1_assay_area.replace(0, np.nan)
    scale = float(args.microglia_density_scale_um2)
    agg["microglia_density_per_100k_um2"] = agg["microglia_density_per_um2"] * scale
    for state in MICROGLIA_STATES:
        token = norm_token(state)
        col = f"microglia_{token}_cells"
        if col not in agg: agg[col] = 0
        agg[f"microglia_{token}_fraction"] = agg[col] / agg["microglia_cells"].replace(0, np.nan)
        agg[f"microglia_{token}_density_per_um2"] = agg[col] / iba1_assay_area.replace(0, np.nan)
        agg[f"microglia_{token}_density_per_100k_um2"] = agg[f"microglia_{token}_density_per_um2"] * scale
    activated = [s.strip() for s in str(args.activated_states).split(",") if s.strip()]
    resting = [s.strip() for s in str(args.resting_states).split(",") if s.strip()]
    agg["activated_microglia_cells"] = sum(agg.get(f"microglia_{norm_token(s)}_cells", 0) for s in activated)
    agg["resting_microglia_cells"] = sum(agg.get(f"microglia_{norm_token(s)}_cells", 0) for s in resting)
    agg["activated_microglia_fraction"] = agg["activated_microglia_cells"] / agg["microglia_cells"].replace(0, np.nan)
    agg["resting_microglia_fraction"] = agg["resting_microglia_cells"] / agg["microglia_cells"].replace(0, np.nan)
    agg["activated_microglia_density_per_um2"] = agg["activated_microglia_cells"] / iba1_assay_area.replace(0, np.nan)
    agg["activated_microglia_density_per_100k_um2"] = agg["activated_microglia_density_per_um2"] * scale
    agg = attach_metadata_by_animal(agg, samplesheet)
    return agg


def _gi_endpoint_specs() -> List[Tuple[str, str, str, str]]:
    return [
        ("activated_microglia_fraction", "Activated + Amoeboid / total Iba1 cells", "Reactive Microglia Fraction per Animal (Activated + Amoeboid)", "primary"),
        ("microglia_amoeboid_fraction", "Amoeboid / total Iba1 cells", "Amoeboid Microglia Fraction per Animal", "primary"),
        ("microglia_density_per_um2", "Iba1 cells density", "Total Microglia Density per Animal", "primary"),
        ("gfap_area_fraction", "GFAP-positive area fraction", "GFAP Area Fraction per Animal", "primary"),
        ("gfap_corrected_intensity_per_um2", "GFAP corrected intensity density", "GFAP Corrected Intensity per Animal", "primary"),
        ("microglia_ramified_fraction", "Ramified / total Iba1 cells", "Ramified Microglia Fraction per Animal", "secondary"),
        ("microglia_rodlike_fraction", "Rod-like / total Iba1 cells", "Rod-like Microglia Fraction per Animal", "secondary"),
        ("microglia_activated_fraction", "Activated / total Iba1 cells", "Activated Microglia Fraction per Animal", "secondary"),
        ("iba1_area_fraction", "Iba1-positive area fraction", "Iba1 Area Fraction per Animal", "secondary"),
        ("iba1_corrected_intensity_per_um2", "Iba1 corrected intensity density", "Iba1 Corrected Intensity per Animal", "secondary"),
        ("dapi_density_per_um2", "DAPI nuclei density", "DAPI Density per Animal", "qc"),
    ]


def write_gfap_iba1_stats(animal_region_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    if animal_region_df is None or animal_region_df.empty:
        return
    values_rows: List[Dict[str, Any]] = []
    for metric, ylabel, title, tier in _gi_endpoint_specs():
        for _, r in animal_region_df.iterrows():
            values_rows.append({
                "animal_id": r.get("animal_id"), "animal_key": r.get("animal_key"), "region": r.get("region"),
                "cond": r.get("cond"), "genotype": r.get("genotype"), "sex": r.get("sex"), "cage": r.get("cage"),
                "endpoint": metric, "endpoint_tier": tier, "value": pd.to_numeric(pd.Series([r.get(metric)]), errors="coerce").iloc[0],
                "ylabel": ylabel, "title": title,
            })
    valdf = pd.DataFrame(values_rows)
    valdf.to_csv(outdir / "gfap_iba1_endpoint_values_long.csv", index=False)
    plan = pd.DataFrame([{"endpoint": m, "tier": t, "interpretation": y, "plot_title": ti, "statistical_unit": "animal", "regions": "ME and ARC separately"} for m, y, ti, t in _gi_endpoint_specs()])
    plan.to_csv(outdir / "gfap_iba1_definitive_analysis_plan.csv", index=False)

    def make_table(filtered: bool) -> Any:
        rows: List[Dict[str, Any]] = []
        d = valdf.copy(); d["cond_clean"] = d["cond"].map(clean_condition); d = d[~d["cond_clean"].map(is_control_condition)]
        for (metric, tier, region), sub in d.groupby(["endpoint", "endpoint_tier", "region"], dropna=False):
            conds = condition_order(sub["cond_clean"])
            arrays: Dict[str, Any] = {}
            for c in conds:
                v = pd.to_numeric(sub[sub["cond_clean"] == c]["value"], errors="coerce").dropna().to_numpy(float)
                outmask = iqr_outlier_mask(v, k=float(args.iqr_k), method=str(args.outlier_method))
                arrays[c] = v[~outmask] if filtered else v
                rows.append({"endpoint": metric, "endpoint_tier": tier, "region": region, "test": "descriptive", "comparison": c, "n": len(arrays[c]), "mean": float(np.mean(arrays[c])) if len(arrays[c]) else np.nan, "sd": float(np.std(arrays[c], ddof=1)) if len(arrays[c]) > 1 else np.nan, "median": float(np.median(arrays[c])) if len(arrays[c]) else np.nan, "iqr_filtered": filtered})
            groups = [arrays[c] for c in conds]
            rows.append({"endpoint": metric, "endpoint_tier": tier, "region": region, "test": "one-way ANOVA", "p_value_method": "ordinary equal-variance parametric F distribution", "comparison": "global", "p_value": safe_anova(groups), "iqr_filtered": filtered})
            for i in range(len(conds)):
                for j in range(i + 1, len(conds)):
                    a, b = arrays[conds[i]], arrays[conds[j]]
                    rows.append({"endpoint": metric, "endpoint_tier": tier, "region": region, "test": "Mann-Whitney U", "p_value_method": "exhaustive two-sided independent-label permutation of rank-sum/U", "comparison": f"{conds[i]} vs {conds[j]}", "p_value": safe_mwu(a, b), "effect_cliffs_delta": cliffs_delta(a, b), "n_group1": len(a), "n_group2": len(b), "iqr_filtered": filtered})
        tab = pd.DataFrame(rows)
        if not tab.empty and "p_value" in tab:
            tab["q_value_BH"] = np.nan
            idx = tab["p_value"].notna()
            tab.loc[idx, "q_value_BH"] = _bh_adjust(tab.loc[idx, "p_value"].astype(float).tolist())
        return tab
    allv = make_table(False); filtv = make_table(True)
    allv.to_csv(outdir / "gfap_iba1_definitive_statistics_all_values.csv", index=False)
    filtv.to_csv(outdir / "gfap_iba1_definitive_statistics_iqr_filtered.csv", index=False)
    primary = filtv if bool(args.exclude_outliers_stats) else allv
    primary.to_csv(outdir / "gfap_iba1_definitive_statistics.csv", index=False)
    logger.warning("Definitive GFAP/Iba1 animal-level statistics CSVs written.")


def make_gfap_iba1_plots(animal_region_df: Any, outdir: Path, args: argparse.Namespace, logger: logging.Logger) -> None:
    if animal_region_df is None or animal_region_df.empty:
        return
    from matplotlib.backends.backend_pdf import PdfPages
    plot_dir = outdir / "plots"; ensure_dir(plot_dir)
    specs = [
        ("activated_microglia_fraction", "Activated + Amoeboid / total Iba1", "Reactive Microglia Fraction per Animal (Activated + Amoeboid)", "activated_microglia_fraction"),
        ("microglia_density_per_um2", "Iba1 cells density", "Total Microglia Density per Animal", "total_microglia_density"),
        ("microglia_ramified_fraction", "Ramified / total Iba1", "Ramified Microglia Fraction per Animal", "ramified_fraction"),
        ("microglia_rodlike_fraction", "Rod-like / total Iba1", "Rod-like Microglia Fraction per Animal", "rod_like_fraction"),
        ("microglia_activated_fraction", "Activated / total Iba1", "Activated Microglia Fraction per Animal", "activated_fraction"),
        ("microglia_amoeboid_fraction", "Amoeboid / total Iba1", "Amoeboid Microglia Fraction per Animal", "amoeboid_fraction"),
        ("microglia_ramified_density_per_um2", "Ramified density", "Ramified Microglia Density per Animal", "ramified_density"),
        ("microglia_rodlike_density_per_um2", "Rod-like density", "Rod-like Microglia Density per Animal", "rod_like_density"),
        ("microglia_activated_density_per_um2", "Activated density", "Activated Microglia Density per Animal", "activated_density"),
        ("microglia_amoeboid_density_per_um2", "Amoeboid density", "Amoeboid Microglia Density per Animal", "amoeboid_density"),
        ("gfap_area_fraction", "GFAP-positive area fraction", "GFAP Area Fraction per Animal", "gfap_area_fraction"),
        ("gfap_corrected_intensity_per_um2", "GFAP corrected intensity density", "GFAP Corrected Intensity per Animal", "gfap_corrected_intensity"),
        ("gfap_reactivity_index", "GFAP area × intensity index", "GFAP Reactivity Index per Animal", "gfap_reactivity_index"),
        ("iba1_area_fraction", "Iba1-positive area fraction", "Iba1 Area Fraction per Animal", "iba1_area_fraction"),
        ("iba1_corrected_intensity_per_um2", "Iba1 corrected intensity density", "Iba1 Corrected Intensity per Animal", "iba1_corrected_intensity"),
    ]
    # build_final_arc_me_figure automatically applies the exact v18/v19 style.
    formats = _final_plot_formats(args)
    with PdfPages(plot_dir / "ARC_ME_GFAP_IBA1_ALL_ENDPOINTS.pdf") as pdf:
        for metric, ylabel, title, slug in specs:
            fig = build_final_arc_me_figure(animal_region_df, metric, ylabel, title, args)
            if fig is None:
                continue
            for fmt in formats:
                kw = {"bbox_inches": "tight"}
                if fmt == "png": kw["dpi"] = int(args.plot_dpi)
                fig.savefig(plot_dir / f"ARC_ME_{slug}.{fmt}", **kw)
            pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
    animal_region_df.to_csv(plot_dir / "final_arc_me_gfap_iba1_plot_values.csv", index=False)
    logger.warning("GFAP/Iba1 final ME/ARC figures written to %s", plot_dir)


def _gi_load_region_info(processed_dir: Path) -> Tuple[Dict[str, Any], Any]:
    info: Dict[str, Any] = {}
    jp = processed_dir / "region_geometry.json"
    if jp.exists(): info.update(json.loads(jp.read_text(encoding="utf-8")))
    npz = processed_dir / "region_geometry_arrays.npz"
    if npz.exists():
        z = np.load(npz)
        for k in z.files: info[k] = z[k]
    rm_path = next(iter(processed_dir.glob("*_ARC_ME_labels.tif")), None)
    rm = read_image(rm_path) if rm_path else None
    return info, np.asarray(rm).astype(np.uint8) if rm is not None else None


def _gi_draw_scale_bar(ax: Any, H: int, W: int, args: argparse.Namespace) -> None:
    from matplotlib.patches import Rectangle
    from matplotlib import patheffects as pe
    um = float(args.scalebar_um); px = int(round(um / max(float(args.um_per_px), 1e-9)))
    px = min(px, int(0.30 * W)); y = H - 28; x = W - px - 30
    ax.add_patch(Rectangle((x, y), px, 6, facecolor="white", edgecolor="black", linewidth=1.5))
    ax.text(x + px / 2, y - 8, f"{int(um)} µm", color="white", fontsize=10, ha="center", va="bottom", fontweight="bold", path_effects=[pe.withStroke(linewidth=2, foreground="black")])


def render_gfap_iba1_qc(s: GISection, cell_sub: Any, section_region_sub: Any, args: argparse.Namespace, outdir: Path, logger: logging.Logger) -> None:
    from skimage.measure import find_contours
    from matplotlib.lines import Line2D
    from matplotlib import patheffects as pe
    base, _rep, _lab = strip_repetition_suffix(s.sample)
    source_sample = str(s.sample)
    proc = outdir / args.reconstructed_dirname / source_sample / s.section
    dpath = next(iter(proc.glob("*_DAPI_labels.tif")), None); ipath = next(iter(proc.glob("*_Iba1_cells.tif")), None); gpath = next(iter(proc.glob("*_GFAP_positive.tif")), None)
    raw_ipath = next(iter(proc.glob("*_Iba1_cells_raw.tif")), None); rejected_ipath = next(iter(proc.glob("*_Iba1_cells_rejected.tif")), None)
    if not dpath or not ipath:
        return
    dapi = np.asarray(read_image(dpath)).astype(np.int32); iba = np.asarray(read_image(ipath)).astype(np.int32); gf = np.asarray(read_image(gpath)).astype(bool) if gpath else np.zeros_like(dapi, bool)
    raw_iba = np.asarray(read_image(raw_ipath)).astype(np.int32) if raw_ipath else iba
    rejected_iba = np.asarray(read_image(rejected_ipath)).astype(np.int32) if rejected_ipath else np.zeros_like(iba)
    info, rm = _gi_load_region_info(proc)
    if rm is None: return
    H, W = dapi.shape
    rgb = np.zeros((H, W, 3), dtype=float); rgb[dapi > 0] = 0.78
    # GFAP positive area as a subtle magenta underlay.
    rgb[gf] = 0.60 * rgb[gf] + 0.40 * np.array([1.0, 0.0, 0.65])
    state_by_label = {int(r.iba1_label): str(r.microglia_state) for r in cell_sub.itertuples()}
    for lab0, state in state_by_label.items():
        col = matplotlib.colors.to_rgb(MICROGLIA_STATE_COLORS.get(state, "#777777"))
        rgb[iba == lab0] = col
    fig = plt.figure(figsize=(14.5, 5.6), dpi=150)
    gs = fig.add_gridspec(1, 2, width_ratios=[5.3, 1.7]); ax = fig.add_subplot(gs[0, 0]); panel = fig.add_subplot(gs[0, 1])
    ax.imshow(rgb, interpolation="nearest"); ax.set_axis_off(); ax.set_xlim(0, W); ax.set_ylim(H, 0)
    for c in find_contours(rejected_iba > 0, 0.5):
        ax.plot(c[:, 1], c[:, 0], color="#ff5a5f", linewidth=0.45, alpha=0.80)
    for code, color, lw in [(REGION_NAME_TO_CODE["ME"], "#f39c12", 1.4), (REGION_NAME_TO_CODE["ARC"], "#12239e", 1.0)]:
        for c in find_contours(rm == code, 0.5): ax.plot(c[:, 1], c[:, 0], color=color, linewidth=lw)
    wy = np.asarray(info.get("_wall_y", [])); lx = np.asarray(info.get("_wall_left_x", [])); rx = np.asarray(info.get("_wall_right_x", []))
    if wy.size and lx.size == wy.size: ax.plot(lx, wy, color="white", linewidth=2.0, path_effects=[pe.withStroke(linewidth=3.5, foreground="#12239e")])
    if wy.size and rx.size == wy.size: ax.plot(rx, wy, color="white", linewidth=2.0, path_effects=[pe.withStroke(linewidth=3.5, foreground="#12239e")])
    fx = np.asarray(info.get("_me_floor_curve_x", [])); fy = np.asarray(info.get("_me_floor_curve_y", []))
    if fx.size and fy.size == fx.size: ax.plot(fx, fy, color="#f39c12", linewidth=1.8, linestyle="--")
    _gi_draw_scale_bar(ax, H, W, args)
    panel.set_axis_off(); panel.set_xlim(0, 1); panel.set_ylim(0, 1)
    raw_n = int(len(np.unique(raw_iba)) - (1 if np.any(raw_iba == 0) else 0))
    accepted_n = int(len(np.unique(iba)) - (1 if np.any(iba == 0) else 0))
    rejected_n = int(len(np.unique(rejected_iba)) - (1 if np.any(rejected_iba == 0) else 0))
    lines = [f"Animal: {base}", f"Source folder: {s.sample}", f"Section: {s.section}", "", f"DAPI nuclei: {int(dapi.max())}", f"Iba1 raw objects: {raw_n}", f"Iba1 QC accepted: {accepted_n}", f"Iba1 QC rejected: {rejected_n}"]
    for region in ("ME", "ARC"):
        rr = section_region_sub[section_region_sub["region"] == region]
        n = int((cell_sub["region"] == region).sum())
        lines += ["", region, f"  microglia: {n}"]
        for st in MICROGLIA_STATES:
            lines.append(f"  {st}: {int(((cell_sub['region']==region)&(cell_sub['microglia_state']==st)).sum())}")
        if not rr.empty:
            lines.append(f"  GFAP area: {100*float(rr.iloc[0].get('gfap_area_fraction', np.nan)):.2f}%")
            lines.append(f"  GFAP corr. int./µm²: {float(rr.iloc[0].get('gfap_corrected_intensity_per_um2', np.nan)):.4g}")
    lines += ["", f"Classifier: {cell_sub['classifier_source'].mode().iloc[0] if not cell_sub.empty else 'NA'}", f"Uncertain cells: {int(cell_sub['classifier_uncertain'].fillna(False).sum()) if not cell_sub.empty else 0}", "", f"ME floor detected: {bool(info.get('me_floor_present', False))}", f"Oval ventricle: {bool(info.get('wall_oval_detected', False))}"]
    panel.text(0.02, 0.98, "\n".join(lines), ha="left", va="top", fontsize=9.0, fontweight="bold")
    handles = [Line2D([0],[0], marker="s", linestyle="none", markerfacecolor=MICROGLIA_STATE_COLORS[s0], markeredgecolor="black", label=s0) for s0 in MICROGLIA_STATES]
    handles.append(Line2D([0],[0], marker="s", linestyle="none", markerfacecolor="#ff0090", markeredgecolor="none", alpha=0.5, label="GFAP+ area"))
    panel.legend(handles=handles, loc="lower left", frameon=True, fontsize=9, prop={"weight":"bold", "size":9})
    fig.tight_layout()
    qdir = outdir / args.reconstructed_qc_dirname / source_sample; ensure_dir(qdir)
    stem = qdir / f"{source_sample}_{s.section}_GFAP_Iba1_ARC_ME_qc"
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_gfap_iba1_analysis(sections: List[GISection], args: argparse.Namespace, outdir: Path, samplesheet: Any, logger: logging.Logger) -> Tuple[Any, Any, Any]:
    all_cells: List[Dict[str, Any]] = []; all_regions: List[Dict[str, Any]] = []; section_rows: List[Dict[str, Any]] = []
    for i, s in enumerate(sections, start=1):
        logger.warning("[quantify %d/%d] %s %s", i, len(sections), s.sample, s.section)
        try:
            c, r, sec = process_gfap_iba1_section(s, args, outdir, logger)
            all_cells.extend(c); all_regions.extend(r); section_rows.append(sec)
        except ManualRegionMaskError:
            logger.exception("Required/reviewed manual ARC/ME mask failed: %s %s", s.sample, s.section)
            raise
        except Exception:
            logger.exception("GFAP/Iba1 section failed: %s %s", s.sample, s.section)
            base, rep, lab = strip_repetition_suffix(s.sample)
            section_rows.append({"sample": s.sample, "animal_id": base, "animal_key": canonical_sample_id(base), "is_repetition": rep, "repetition_label": lab, "section": s.section, "section_index": s.section_index, "status": "exception"})
    candidate_df = pd.DataFrame(all_cells); section_region_df = pd.DataFrame(all_regions); section_df = pd.DataFrame(section_rows)
    candidate_df.to_csv(outdir / "microglia_object_qc.csv", index=False)
    qc_summary_rows: List[Dict[str, Any]] = []
    if not candidate_df.empty:
        if "cell_qc_pass" not in candidate_df.columns:
            raise RuntimeError("Object QC failed closed: candidate table lacks cell_qc_pass")
        qc_ok = candidate_df["cell_qc_pass"].fillna(False).astype(bool)
        rejected_df = candidate_df.loc[~qc_ok].copy()
        cell_df = candidate_df.loc[qc_ok].copy()
        group_cols = ["sample", "animal_id", "section", "region"]
        for keys, sub in candidate_df.groupby(group_cols, dropna=False, sort=True):
            key_values = keys if isinstance(keys, tuple) else (keys,)
            prefix = dict(zip(group_cols, key_values))
            sub_ok = sub["cell_qc_pass"].fillna(False).astype(bool)
            raw_n = int(len(sub)); accepted_n = int(sub_ok.sum()); rejected_n = int(raw_n - accepted_n)
            qc_summary_rows.append({
                **prefix, "raw_candidates": raw_n, "qc_accepted": accepted_n,
                "qc_rejected": rejected_n, "rejection_reason": "ALL",
                "reason_count": rejected_n,
                "rejection_fraction": float(rejected_n / raw_n) if raw_n else float("nan"),
            })
            reason_counts: Dict[str, int] = {}
            for reason_text in sub.loc[~sub_ok, "cell_qc_rejection_reasons"].fillna("").astype(str):
                for reason in (x for x in reason_text.split(";") if x):
                    reason_counts[reason] = reason_counts.get(reason, 0) + 1
            for reason, count in sorted(reason_counts.items()):
                qc_summary_rows.append({
                    **prefix, "raw_candidates": raw_n, "qc_accepted": accepted_n,
                    "qc_rejected": rejected_n, "rejection_reason": reason,
                    "reason_count": int(count),
                    "rejection_fraction": float(count / raw_n) if raw_n else float("nan"),
                })
    else:
        rejected_df = candidate_df.copy()
        cell_df = candidate_df.copy()
    rejected_df.to_csv(outdir / "per_microglia_rejected_objects.csv", index=False)
    pd.DataFrame(qc_summary_rows).to_csv(outdir / "microglia_object_qc_summary.csv", index=False)
    if not cell_df.empty:
        cell_df, classifier_meta = classify_microglia(cell_df, args, outdir, logger)
    else:
        classifier_meta = {}
    classifier_meta["biological_cell_qc"] = {
        "version": MICROGLIA_CELL_QC_VERSION,
        "raw_candidates_n": int(len(candidate_df)),
        "accepted_n": int(len(cell_df)),
        "rejected_n": int(len(rejected_df)),
        "reason_counts": (
            rejected_df["cell_qc_rejection_reasons"].fillna("").astype(str)
            .str.split(";").explode().loc[lambda x: x != ""].value_counts().to_dict()
            if not rejected_df.empty else {}
        ),
    }
    (outdir / "microglia_classifier_metadata.json").write_text(json.dumps(classifier_meta, indent=2, default=str), encoding="utf-8")
    cell_df.to_csv(outdir / "per_microglia_cell_measurements.csv", index=False)
    section_df = attach_metadata_by_animal(section_df, samplesheet) if not section_df.empty else section_df
    section_df.to_csv(outdir / "per_section_summary.csv", index=False)
    section_region_df = attach_metadata_by_animal(section_region_df, samplesheet) if not section_region_df.empty else section_region_df
    section_region_df.to_csv(outdir / "per_section_region_summary.csv", index=False)
    if bool(args.prepare_labels_only):
        return cell_df, section_region_df, pd.DataFrame()
    animal_region = _gi_aggregate_animal_regions(cell_df, section_region_df, samplesheet, args)
    animal_region.to_csv(outdir / "per_animal_region_summary.csv", index=False)
    # Animal total across ARC+ME, useful for global QC but not used in region figures.
    total_rows: List[Dict[str, Any]] = []
    if not animal_region.empty:
        for (aid, akey), sub in animal_region.groupby(["animal_id", "animal_key"], dropna=False):
            total_rows.append({"animal_id": aid, "animal_key": akey, "region_area_um2": float(sub["region_area_um2"].sum()), "microglia_cells": int(sub["microglia_cells"].sum()), "activated_microglia_cells": int(sub["activated_microglia_cells"].sum()), "n_regions_present": int((sub["region_area_um2"]>0).sum())})
    animal_total = pd.DataFrame(total_rows)
    animal_total = attach_metadata_by_animal(animal_total, samplesheet) if not animal_total.empty else animal_total
    animal_total.to_csv(outdir / "per_animal_total_summary.csv", index=False)
    # Render state-aware QC only after classification.
    if bool(args.save_qc) and not cell_df.empty:
        for s in sections:
            base, _r, _l = strip_repetition_suffix(s.sample)
            csub = cell_df[(cell_df["sample"] == s.sample) & (cell_df["section_index"] == s.section_index)]
            rsub = section_region_df[(section_region_df["sample"] == s.sample) & (section_region_df["section_index"] == s.section_index)]
            try: render_gfap_iba1_qc(s, csub, rsub, args, outdir, logger)
            except Exception: logger.exception("QC failed for %s %s", s.sample, s.section)
    return cell_df, section_region_df, animal_region


def main_gfap_iba1() -> None:
    ap = build_argparser_gfap_iba1()
    args = ap.parse_args()
    default_manual_root = DEFAULT_MANUAL_REGION_ROOT

    # Native image-by-image review is the only allowed ARC/ME source for this assay.
    args.human_in_the_loop = True
    args.require_manual_regions = True
    if not args.manual_region_dir:
        args.manual_region_dir = str(default_manual_root)
    if not 0 <= int(args.human_review_port) <= 65535:
        ap.error("--human-review-port must be between 0 and 65535")
    manual_root = Path(args.manual_region_dir).expanduser().resolve()
    manual_root.mkdir(parents=True, exist_ok=True)
    args.manual_region_dir = str(manual_root)
    reviewed_dapi_root = Path(args.reviewed_dapi_root).expanduser().resolve()
    reviewed_iba1_root = Path(args.reviewed_iba1_root).expanduser().resolve()
    args.reviewed_dapi_root = str(reviewed_dapi_root)
    args.reviewed_iba1_root = str(reviewed_iba1_root)
    if args.full_run:
        args.migrate_only = False
        args.segment_only = False
        args.quantify_only = False
    if bool(args.force) and bool(args.dry_run):
        ap.error("--force cannot be combined with --dry-run")

    bootstrap_conda_if_needed(args)
    import_science_stack()

    raw_arg = Path(args.source_zip).resolve() if args.source_zip else Path(args.raw_root).resolve()
    work_root = Path(args.sanitized_root).resolve()
    outdir_arg = Path(args.outdir) if args.outdir else DEFAULT_ANALYSIS_OUTPUT
    outdir_candidate = outdir_arg.expanduser().resolve()
    samplesheet_path = Path(args.samplesheet).resolve() if args.samplesheet else None
    migration_report = None

    def log_configuration(logger: logging.Logger, output: Path) -> None:
        logger.warning("=== Keyence 10x DAPI/GFAP/Iba1 ARC/ME microglia pipeline v%s ===", GFAP_IBA1_VERSION)
        logger.warning("Raw source or ZIP: %s", raw_arg)
        logger.warning("Sanitized working folder: %s", work_root)
        logger.warning("Output: %s", output)
        logger.warning("Channel map: CH1=DAPI, CH3=GFAP, CH4=Iba1")
        logger.warning("Each XY## is treated as one complete bilateral section; no left/right reconstruction or statistics.")
        logger.warning(
            "Native HIL ARC/ME/VMN policy: exactly %d/%d accepted TIFF+JSON receipts are required; "
            "automated fallback and registered POMC masks are disabled. Root=%s",
            GFAP_IBA1_EXPECTED_REVIEW_SECTIONS,
            GFAP_IBA1_EXPECTED_REVIEW_SECTIONS,
            manual_root,
        )
        logger.warning(
            "Accepted sources: manual_replacement_polygons or verified native automated_mask_unchanged; "
            "ME may be absent only with an explicit anatomical reason."
        )
        logger.warning(
            "Versioned cell masks: DAPI=%s | Iba1=%s | historical_override=%s",
            reviewed_dapi_root,
            reviewed_iba1_root,
            bool(args.allow_historical_cell_masks),
        )
        logger.warning("Microglial CNN mode=%s. Without reviewed cell-state labels, CNN targets remain high-confidence morphometric pseudo-labels.", args.microglia_classifier_mode)
        logger.warning("Pixel size: %.4f µm/px. Override --um-per-px if Keyence calibration differs.", float(args.um_per_px))

    # All discovery, segmentation preparation and review happen in a disposable
    # preflight area. A failed/incomplete review therefore cannot erase the
    # canonical output even when --force was requested.
    with tempfile.TemporaryDirectory(prefix="gfap-native-hil-preflight-") as preflight_text:
        preflight_dir = Path(preflight_text)
        preflight_logger = setup_logger(preflight_dir, bool(args.verbose))
        log_configuration(preflight_logger, outdir_candidate)
        samplesheet = load_samplesheet(samplesheet_path, preflight_dir, preflight_logger)

        existing_dapi = (
            list(work_root.glob("*/S*/Image_*_DAPI.tif"))
            + list(work_root.glob("*/S*/Image_*_DAPI.tiff"))
            if work_root.exists() else []
        )
        migration_needed = (not work_root.exists()) or bool(args.force_migrate) or not existing_dapi
        if migration_needed:
            source_root = _gi_prepare_raw_source(args, work_root, preflight_logger)
            migration_report = migrate_gfap_iba1_to_sanitized(source_root, work_root, args, preflight_logger)
            migration_report.to_csv(preflight_dir / "migration_report.csv", index=False)
        else:
            preflight_logger.warning("Existing sanitized folder detected; migration skipped and existing masks are preserved.")

        sections = discover_gfap_iba1_sections(work_root, args, preflight_logger)
        _gi_section_table(sections, samplesheet).to_csv(preflight_dir / "discovered_sections.csv", index=False)
        if bool(args.segment_missing) and not bool(args.quantify_only):
            segment_gfap_iba1_missing(sections, args, preflight_logger, preflight_dir)
            sections = discover_gfap_iba1_sections(work_root, args, preflight_logger)
        else:
            preflight_logger.warning("Missing-mask segmentation skipped.")

        try:
            _gi_require_raw_local_masks(sections, work_root)
            _gi_apply_reviewed_cell_masks(sections, args, preflight_dir, preflight_logger)
            coverage_rows = _gi_require_complete_human_review(
                sections,
                args,
                work_root,
                preflight_dir,
                outdir_candidate,
                preflight_logger,
            )
        except ManualRegionMaskError as exc:
            raise SystemExit(f"ERROR: {exc}") from exc

    preserved_native_proposals: List[Path] = []
    for row in coverage_rows:
        if (
            row.get("review_status") == "accepted"
            and row.get("accepted_source") == "automated_mask_unchanged"
            and row.get("proposal_mask_sha256_verified") is True
            and str(row.get("proposal_mask_path") or "").strip()
        ):
            proposal = Path(str(row["proposal_mask_path"])).expanduser().resolve()
            if proposal == outdir_candidate or _path_is_within(proposal, outdir_candidate):
                preserved_native_proposals.append(proposal)

    if bool(args.force):
        try:
            outdir = _force_clean_generated_dir(
                outdir_candidate,
                forbidden_roots=(work_root, raw_arg, PROJECT_ROOT, PROJECT_ROOT / "scripts"),
                preserved_roots=(
                    default_manual_root,
                    manual_root,
                    reviewed_dapi_root,
                    reviewed_iba1_root,
                    *preserved_native_proposals,
                ),
            )
        except ValueError as exc:
            ap.error(str(exc))
    else:
        outdir = outdir_candidate
    ensure_dir(outdir)
    logger = setup_logger(outdir, bool(args.verbose))
    log_configuration(logger, outdir)
    if preserved_native_proposals:
        logger.warning(
            "Preserved %d verified native proposal provenance file(s) through output cleanup.",
            len(preserved_native_proposals),
        )

    samplesheet = load_samplesheet(samplesheet_path, outdir, logger)
    if migration_report is not None:
        migration_report.to_csv(outdir / "migration_report.csv", index=False)
    try:
        _gi_apply_reviewed_cell_masks(sections, args, outdir, logger)
    except ManualRegionMaskError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    _gi_section_table(sections, samplesheet).to_csv(outdir / "discovered_sections.csv", index=False)
    final_coverage_rows, final_pending = _gi_audit_manual_region_coverage(sections, args, outdir, logger)
    if final_pending or len(final_coverage_rows) != GFAP_IBA1_EXPECTED_REVIEW_SECTIONS:
        raise SystemExit(
            "ERROR: native HIL review provenance changed after preflight; "
            "quantification was not started"
        )
    if args.migrate_only or args.segment_only:
        logger.warning("Native ARC/ME/VMN review preflight complete; requested preparation-only mode finished.")
        return

    cell_df, section_region_df, animal_region_df = run_gfap_iba1_analysis(
        sections, args, outdir, samplesheet, logger
    )
    logger.warning("Per-cell Iba1 rows: %d", len(cell_df))
    logger.warning("Per-animal ARC/ME rows: %d", len(animal_region_df))
    if bool(args.prepare_labels_only):
        logger.warning("Annotation preparation complete. Review microglia_annotation_template.csv, fill state, then rerun with --microglia-labels-csv.")
        return
    write_gfap_iba1_stats(animal_region_df, outdir, args, logger)
    if bool(args.make_plots):
        make_gfap_iba1_plots(animal_region_df, outdir, args, logger)
    logger.warning("=== DONE ===")
    logger.warning("Main outputs:")
    for name in [
        "human_review_coverage.csv",
        "reviewed_cell_mask_coverage.csv",
        "per_microglia_cell_measurements.csv",
        "per_section_region_summary.csv",
        "per_animal_region_summary.csv",
        "gfap_iba1_definitive_statistics.csv",
        "microglia_annotation_template.csv",
        "plots/ARC_ME_GFAP_IBA1_ALL_ENDPOINTS.pdf",
    ]:
        logger.warning("  %s", outdir / name)
    # The morphology reviewer is only useful once the patches exist, which is now.
    # Say how to open it every time, and open it here when asked, because this is
    # the moment a reviewer would want it.
    microglia_command = shlex.join(_gi_microglia_review_command(args, outdir))
    logger.warning("Microglia morphology review (turns model proposals into "
                   "reviewed labels): %s", microglia_command)
    if bool(getattr(args, "review_microglia", False)):
        logger.warning("Opening the microglia morphology reviewer; stop it with Ctrl-C.")
        try:
            subprocess.run(_gi_microglia_review_command(args, outdir), check=False)
        except KeyboardInterrupt:
            logger.warning("Microglia morphology reviewer stopped.")
    if bool(getattr(args, "make_figure", False)):
        _gi_render_figure_5(args, outdir, work_root, logger)


if __name__ == "__main__":
    main_gfap_iba1()
