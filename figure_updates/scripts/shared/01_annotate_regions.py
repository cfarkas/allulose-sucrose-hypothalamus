#!/usr/bin/env python3
"""Local human-in-the-loop ARC/ME/VMN polygon annotation server.

Accepted annotations are paired JSON and uint8 TIFF files. Coordinates remain in
original DAPI pixels; label codes are background=0, ARC=1, ME=2, VMN=3.

Any nonempty subset of the three regions is a valid review. Regions without a
polygon are recorded as ``not_delineated`` unless the reviewer explicitly marks
them anatomically absent and supplies a reason. This distinction is persisted
in the receipt and is never inferred from code 0 during downstream analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import ipaddress
import json
import math
import os
import re
import shlex
import subprocess
import sys
import threading
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from functools import lru_cache
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
import tifffile
from PIL import Image
from skimage.draw import polygon as draw_polygon
from skimage.morphology import binary_dilation, disk
from skimage.segmentation import find_boundaries

HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file()))), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")
ROOT = PAPER_ROOT
FIG4_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig4"
FIG5_ANALYSIS_ROOT = PAPER_ROOT / "analyses" / "Fig5"



def _first_existing(*candidates):
    """The first candidate that exists, else the first candidate.

    Directories under analyses/ are generated working output and are absent in a
    freshly unpacked copy of the tree, so every default that names one has to
    fall back to the equivalent input this tree actually ships. Returning the
    first candidate when nothing exists keeps the error message pointed at the
    working copy the author would have expected.
    """
    for candidate in candidates:
        if candidate is not None and candidate.exists():
            return candidate
    return candidates[0]


def _latest(parent, pattern):
    """The newest-named directory matching ``pattern``, or None."""
    if not parent.is_dir():
        return None
    matches = sorted((path for path in parent.glob(pattern) if path.is_dir()),
                     key=lambda path: path.name)
    return matches[-1] if matches else None


_SHIPPED_POMC_REVIEW_INPUT = _latest(
    PAPER_ROOT / "Fig4" / "hil_review", "prepared_review_input_*"
)
DEFAULT_INPUT = _first_existing(
    FIG5_ANALYSIS_ROOT / "raw" / "input", PAPER_ROOT / "Fig5" / "raw_data",
)
DEFAULT_ANALYSIS = _first_existing(
    FIG5_ANALYSIS_ROOT / "results" / "final_hybrid_run",
    ((_latest(PAPER_ROOT / "Fig5", "provisional_machine_hil_*") or PAPER_ROOT / "Fig5")
     / "review_input_analysis"),
)
DEFAULT_OUTPUT = FIG5_ANALYSIS_ROOT / "hil" / "accepted_annotations"
DEFAULT_POMC_INPUT = _first_existing(
    FIG4_ANALYSIS_ROOT / "results" / "final_run" / "reconstructed_sections",
    (_SHIPPED_POMC_REVIEW_INPUT / "reconstructed_sections"
     if _SHIPPED_POMC_REVIEW_INPUT is not None else None),
)
# The review images and the registration receipts that bind them have to come
# from the same stage, so this is always the chosen input root's own parent.
DEFAULT_POMC_ANALYSIS = DEFAULT_POMC_INPUT.parent
DEFAULT_POMC_OUTPUT = FIG4_ANALYSIS_ROOT / "hil" / "accepted_annotations"
REGIONS = ("ARC", "ME", "VMN")
LABEL_CODES = {"background": 0, "ARC": 1, "ME": 2, "VMN": 3}
LEGACY_LABEL_CODES = {"background": 0, "ARC": 1, "ME": 2}
COLORS = {
    "ARC": (0, 173, 181, 235),
    "ME": (230, 126, 34, 235),
    "VMN": (190, 62, 174, 235),
}
REGION_STATES = {"drawn", "explicitly_absent", "not_delineated"}
STITCH_REGISTRATION_SCHEMA = "fig4_stitch_registration_v1"
STITCH_REGISTRATION_ALGORITHM = "dapi_integer_overlap_v1"
STITCH_REGISTRATION_STATUSES = {
    "verified_overlap", "manual_review_required", "single_half_native"
}
STITCH_REGISTRATION_METHODS = {
    "sift_ransac_raw_ncc", "direct_ncc", "anatomical_y_fallback",
    "single_half_native",
}
STITCH_MANUAL_CONFIRMATION = "reviewed_registered_composite_image_by_image"
STITCH_AUTO_CONFIRMATIONS = {
    "verified_overlap": "automated_overlap_qc_verified",
    "single_half_native": "single_half_native_no_registration_needed",
}
MAX_BODY = 20 * 1024 * 1024
FIGURE_LABELS = {"pomc": "Figure 4", "gfap": "Figure 5"}
FIGURE_OUTPUT_TOKEN = "{output_dir}"
FIGURE_ANALYSIS_TOKEN = "{analysis_dir}"
FIGURE_OUTPUT_PREFIX = {"pomc": "paper_figure4_render", "gfap": "paper_figure6_render"}
FIGURE_ANALYSIS_PREFIX = {"pomc": "paper_figure4_analysis", "gfap": "paper_figure6_analysis"}
MACHINE_DEFAULT_ATTESTATION = "I REVIEWED THE MACHINE DEFAULT SEGMENTATIONS"
# The reviewer's --analysis-dir holds review inputs and automated proposals. The
# renderers instead need a completed analysis run, which is a different directory
# produced later by the figure's analyzer, so the Make Figure button takes it
# separately and refuses to start until these artifacts are actually present.
FIGURE_ANALYSIS_REQUIREMENTS = {
    "pomc": (
        "per_reconstructed_section_summary.csv",
        "per_cell_reconstructed_measurements.csv",
        "per_reconstructed_section_region_summary.csv",
        "per_animal_region_summary.csv",
        "reconstructed_sections",
    ),
    "gfap": (
        "human_review_coverage.csv",
        "discovered_sections.csv",
        "per_section_region_summary.csv",
        "per_animal_region_summary.csv",
        "per_microglia_cell_measurements.csv",
        "microglia_classifier_metadata.json",
        "processed_sections",
    ),
}
FIGURE_ANALYZERS = {
    "pomc": "scripts/Fig4/02_analyze_pomc_cfos.py",
    "gfap": "scripts/Fig5/01_analyze_gfap_iba1_microglia.py",
}
# Where the button's own analysis run lands when nothing was named. These are the
# directories each renderer already reads by default, so the analysis the button
# produces is also the one the documented plain renderer command picks up, and it
# stays inside this tree instead of a possibly small /tmp.
FIGURE_DEFAULT_ANALYSIS_OUTPUT = {
    "pomc": ROOT / "analyses" / "Fig4" / "results" / "human_final_run",
    "gfap": ROOT / "analyses" / "Fig5" / "results" / "human_final_run",
}
# Where a Make Figure run renders when no --figure-output-root was named. It is
# inside the tree, beside the analysis the button produces, because a render that
# lands in a random /tmp directory under a name nobody recognises reads to the
# reviewer as though the button had done nothing at all.
FIGURE_DEFAULT_OUTPUT_ROOT = {
    "pomc": ROOT / "analyses" / "Fig4" / "figure",
    "gfap": ROOT / "analyses" / "Fig5" / "figure",
}
# Parent of the fresh analysis directory the button falls back to when each
# figure's own default analysis output is already taken. Analysis output belongs
# beside the other analysis runs, not among the renders.
FIGURE_FALLBACK_ANALYSIS_ROOT = {
    "pomc": ROOT / "analyses" / "Fig4" / "results",
    "gfap": ROOT / "analyses" / "Fig5" / "results",
}
# The basename each renderer gives its master PDF/PNG. Figure 4's renderer takes
# it as --figure-name and would otherwise write Figure_ARC_ME_multipanel, which
# is not the name the figure is promoted and cited under; Figure 5's is fixed.
FIGURE_MASTER_NAME = {"pomc": "Figure_4", "gfap": "Figure_GFAP_Iba1_microglia_ARC_ME"}
FIGURE_MASTER_SUFFIXES = (".pdf", ".png")
LOG_TAIL_CHARS = 6000


def source_sha256():
    """SHA-256 of this file on disk, or "" if it cannot be read.

    A reviewer server keeps running the code it was started with. When this file
    is later corrected, the open browser tab still drives the old process and its
    old commands, which is indistinguishable from the fix not working. Comparing
    this against the hash captured at import time lets the page say so.
    """
    try:
        return hashlib.sha256(HERE.read_bytes()).hexdigest()
    except OSError:
        return ""


SOURCE_SHA256_AT_START = source_sha256()


@dataclass(frozen=True)
class Section:
    dataset: str
    animal: str
    section: str
    dapi: Path
    auto: Path | None
    width: int
    height: int
    dapi_sha256: str
    dapi_size: int
    dapi_mtime_ns: int
    auto_sha256: str | None = None
    auto_kind: str = "none"
    registration: Path | None = None
    registration_sha256: str | None = None
    registration_schema_version: str | None = None
    registration_status: str = "not_applicable"
    registration_requires_manual_review: bool = False
    registration_method: str = "not_applicable"
    registration_qc: dict | None = None
    registration_valid: bool = True
    registration_error: str | None = None

    @property
    def key(self):
        return self.animal, self.section


def plane(array):
    a = np.squeeze(np.asarray(array))
    if a.ndim < 2:
        raise ValueError("image has fewer than two dimensions")
    while a.ndim > 2:
        if a.shape[-1] <= 4:
            a = a.max(-1)
        elif a.shape[0] <= 4:
            a = a.max(0)
        else:
            a = a.max(0)
    return a


def dimensions(path):
    with tifffile.TiffFile(path) as tf:
        shape = tuple(int(v) for v in tf.series[0].shape if int(v) != 1)
    if len(shape) < 2:
        raise ValueError(f"cannot infer dimensions: {shape}")
    if len(shape) == 2:
        return shape[1], shape[0]
    if shape[-1] <= 4:
        return shape[-2], shape[-3]
    return shape[-1], shape[-2]


def safe_component(value, field):
    text = str(value).strip()
    if not text or text in (".", "..") or "/" in text or "\\" in text or "\x00" in text:
        raise ValueError(f"unsafe or empty {field}: {value!r}")
    return text


def _registration_document(path, width, height, expected_dapi_sha256=None):
    """Read and strictly validate one immutable Fig4 stitch receipt."""
    path = Path(path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"stitch registration JSON not found: {path}")
    raw = path.read_bytes()
    try:
        record = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid stitch registration JSON: {path}") from exc
    if not isinstance(record, dict):
        raise ValueError(f"stitch registration must be a JSON object: {path}")
    if record.get("schema_version") != STITCH_REGISTRATION_SCHEMA:
        raise ValueError(
            f"unsupported stitch registration schema {record.get('schema_version')!r}: {path}"
        )
    if record.get("algorithm") != STITCH_REGISTRATION_ALGORITHM:
        raise ValueError(
            f"unsupported stitch registration algorithm {record.get('algorithm')!r}: {path}"
        )
    status = str(record.get("status", ""))
    method = str(record.get("method", ""))
    requires_review = record.get("requires_manual_review")
    if status not in STITCH_REGISTRATION_STATUSES:
        raise ValueError(f"invalid stitch registration status {status!r}: {path}")
    if method not in STITCH_REGISTRATION_METHODS:
        raise ValueError(f"invalid stitch registration method {method!r}: {path}")
    if not isinstance(requires_review, bool):
        raise ValueError(f"requires_manual_review must be Boolean: {path}")
    if requires_review != (status == "manual_review_required"):
        raise ValueError(
            f"registration status/requires_manual_review disagree: {path}"
        )
    transform = record.get("transform")
    qc = record.get("qc")
    if not isinstance(transform, dict) or not isinstance(qc, dict):
        raise ValueError(f"registration transform and qc must be objects: {path}")
    output_shape = record.get("output_shape_px")
    if output_shape != [int(height), int(width)]:
        raise ValueError(
            f"registration output_shape_px {output_shape!r} != DAPI shape "
            f"{[int(height), int(width)]}: {path}"
        )
    output = record.get("output")
    if not isinstance(output, dict) or not str(output.get("dapi_path", "")).strip():
        raise ValueError(f"registration output DAPI provenance is missing: {path}")
    if expected_dapi_sha256 is not None and output.get("dapi_sha256") != expected_dapi_sha256:
        raise ValueError(
            "registration output DAPI SHA-256 does not match the displayed "
            f"registered composite: {path}"
        )
    expected_method = {
        "verified_overlap": {"sift_ransac_raw_ncc", "direct_ncc"},
        # A high-confidence overlap may still be escalated for explicit human
        # confirmation when conservative instance-seam QC finds an unmatched
        # object at the hard seam. Preserve and review that exact composite.
        "manual_review_required": {
            "sift_ransac_raw_ncc", "direct_ncc", "anatomical_y_fallback",
        },
        "single_half_native": {"single_half_native"},
    }[status]
    if method not in expected_method:
        raise ValueError(
            f"registration method {method!r} is incompatible with status {status!r}: {path}"
        )
    return record, hashlib.sha256(raw).hexdigest()


def make_section(dataset, animal, section, dapi, auto=None, *, auto_kind="none",
                 registration=None):
    if dataset not in {"gfap", "pomc"}:
        raise ValueError(f"native reviewer does not support dataset {dataset!r}")
    animal = safe_component(animal, "animal")
    section = safe_component(section, "section")
    dapi = Path(dapi).expanduser().resolve()
    if not dapi.is_file():
        raise FileNotFoundError(f"DAPI image not found: {dapi}")
    auto_path = Path(auto).expanduser().resolve() if auto else None
    if auto_path is not None and not auto_path.is_file():
        auto_path = None
    width, height = dimensions(dapi)
    dapi_hash = sha256(dapi)
    stat = dapi.stat()
    auto_hash = sha256(auto_path) if auto_path is not None else None
    registration_path = None
    registration_hash = None
    registration_schema = None
    registration_status = "not_applicable"
    registration_requires_review = False
    registration_method = "not_applicable"
    registration_qc = None
    registration_valid = True
    registration_error = None
    if dataset == "pomc":
        registration_path = Path(registration).expanduser().resolve() if registration else None
        try:
            if registration_path is None:
                raise FileNotFoundError("stitch_registration.json was not supplied")
            registration_record, registration_hash = _registration_document(
                registration_path, width, height, dapi_hash
            )
            registration_schema = str(registration_record["schema_version"])
            registration_status = str(registration_record["status"])
            registration_requires_review = bool(
                registration_record["requires_manual_review"]
            )
            registration_method = str(registration_record["method"])
            registration_qc = dict(registration_record.get("qc", {}))
        except (OSError, ValueError, TypeError) as exc:
            registration_valid = False
            registration_status = "missing_or_invalid"
            registration_requires_review = True
            registration_method = "unavailable"
            registration_error = str(exc)
    return Section(str(dataset), animal, section, dapi, auto_path, width, height,
                   dapi_hash, int(stat.st_size), int(stat.st_mtime_ns), auto_hash,
                   str(auto_kind), registration_path, registration_hash,
                   registration_schema, registration_status,
                   registration_requires_review, registration_method,
                   registration_qc, registration_valid, registration_error)


def discover_gfap(raw_root, analysis):
    if not raw_root.is_dir():
        raise FileNotFoundError(f"raw input directory not found: {raw_root}")
    processed = analysis / "processed_sections"
    rows = []
    for animal_dir in sorted((p for p in raw_root.iterdir() if p.is_dir()), key=lambda p: p.name):
        for section_dir in sorted((p for p in animal_dir.iterdir() if p.is_dir()), key=lambda p: p.name):
            expected = section_dir / f"Image_{section_dir.name}_DAPI.tif"
            candidates = [expected] if expected.is_file() else sorted(section_dir.glob("*DAPI*.tif*"))
            if not candidates:
                continue
            dapi = candidates[0].resolve()
            auto = processed / animal_dir.name / section_dir.name / f"{animal_dir.name}_{section_dir.name}_ARC_ME_labels.tif"
            try:
                rows.append(make_section("gfap", animal_dir.name, section_dir.name, dapi,
                                         auto if auto.is_file() else None,
                                         auto_kind="native_gfap_automated_mask"))
            except Exception as exc:
                print(f"Skipping {dapi}: {exc}", file=sys.stderr)
                continue
    if not rows:
        raise RuntimeError(f"no DAPI animal/section pairs found under {raw_root}")
    return rows


def discover_pomc(root):
    if not root.is_dir():
        raise FileNotFoundError(f"POMC reconstructed-section directory not found: {root}")
    rows = []
    for dapi in sorted(root.rglob("reconstructed_DAPI.tif")):
        if len(dapi.parents) < 2:
            continue
        auto = dapi.with_name("reconstructed_ARC_ME_region_mask.tif")
        registration = dapi.with_name("stitch_registration.json")
        try:
            rows.append(make_section("pomc", dapi.parents[1].name, dapi.parent.name, dapi,
                                     auto if auto.is_file() else None,
                                     auto_kind="native_pomc_automated_mask",
                                     registration=registration))
        except Exception as exc:
            print(f"Skipping {dapi}: {exc}", file=sys.stderr)
    if not rows:
        raise RuntimeError(f"no reconstructed POMC DAPI sections found under {root}")
    return rows


def discover(dataset, raw_root, analysis):
    if dataset == "gfap":
        return discover_gfap(raw_root, analysis)
    if dataset == "pomc":
        return discover_pomc(raw_root)
    raise ValueError(f"unknown dataset: {dataset}")


@lru_cache(maxsize=6)
def dapi_png(path_string, mtime_ns):
    del mtime_ns
    a = plane(tifffile.imread(path_string)).astype(np.float32, copy=False)
    finite = a[np.isfinite(a)]
    sample = finite[finite > 0]
    sample = sample if sample.size else finite
    if not sample.size:
        u8 = np.zeros(a.shape, np.uint8)
    else:
        lo, hi = np.percentile(sample, (0.5, 99.8))
        if not np.isfinite(hi) or hi <= lo:
            lo, hi = float(sample.min()), float(sample.max())
        u8 = np.zeros(a.shape, np.uint8) if hi <= lo else np.round(255 * np.nan_to_num(np.clip((a-lo)/(hi-lo), 0, 1))).astype(np.uint8)
    rgb = np.stack(((u8 * .10).astype(np.uint8), (u8 * .58).astype(np.uint8), u8), axis=-1)
    out = io.BytesIO()
    Image.fromarray(rgb, "RGB").save(out, "PNG")
    return out.getvalue()


@lru_cache(maxsize=6)
def auto_png(path_string, mtime_ns, width, height):
    del mtime_ns
    labels = validated_mask(Path(path_string), width, height, require_nonempty=True)
    rgba = np.zeros((height, width, 4), np.uint8)
    for name in REGIONS:
        code = LABEL_CODES[name]
        boundary = binary_dilation(find_boundaries(labels == code, connectivity=2, mode="thick"), disk(2))
        rgba[boundary] = COLORS[name]
    out = io.BytesIO()
    Image.fromarray(rgba, "RGBA").save(out, "PNG")
    return out.getvalue()


def paths(output, section):
    directory = output / section.animal / section.section
    stem = f"{section.animal}_{section.section}_ARC_ME_labels"
    return directory / f"{stem}.json", directory / f"{stem}.tif"


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def validated_mask(path, width, height, require_nonempty=True):
    labels = plane(tifffile.imread(path))
    if labels.shape != (height, width):
        raise ValueError(f"mask shape {labels.shape} != DAPI shape {(height, width)}: {path}")
    if not np.issubdtype(labels.dtype, np.integer):
        if not np.all(np.isfinite(labels)) or not np.all(labels == np.round(labels)):
            raise ValueError(f"mask is not integer-valued: {path}")
    codes = set(int(value) for value in np.unique(labels))
    allowed = set(LABEL_CODES.values())
    if not codes.issubset(allowed):
        raise ValueError(
            f"mask has invalid codes {sorted(codes)}; expected only "
            f"{sorted(allowed)}: {path}"
        )
    if require_nonempty and not any(code in codes for code in allowed if code):
        raise ValueError(f"mask contains no delineated ARC, ME, or VMN region: {path}")
    return np.asarray(labels, dtype=np.uint8)


def validate_source_dapi(section):
    stat = section.dapi.stat()
    if int(stat.st_size) != section.dapi_size or int(stat.st_mtime_ns) != section.dapi_mtime_ns:
        raise ValueError(f"source DAPI changed after server discovery: {section.dapi}")
    actual = sha256(section.dapi)
    if actual != section.dapi_sha256:
        raise ValueError(f"source DAPI SHA-256 changed after server discovery: {section.dapi}")
    return actual


def validate_stitch_registration(section):
    """Return the current POMC stitch receipt and reject post-discovery changes."""
    if section.dataset != "pomc":
        return None
    if (not section.registration_valid or section.registration is None
            or not section.registration_sha256):
        raise ValueError(
            "POMC registered composite has no valid stitch_registration.json: "
            f"{section.registration_error or section.registration}"
        )
    record, current_hash = _registration_document(
        section.registration, section.width, section.height,
        section.dapi_sha256,
    )
    if current_hash != section.registration_sha256:
        raise ValueError(
            "stitch_registration.json changed after server discovery; restart the "
            f"reviewer after rebuilding/reviewing the composite: {section.registration}"
        )
    if (record.get("status") != section.registration_status
            or bool(record.get("requires_manual_review"))
            != section.registration_requires_manual_review
            or record.get("method") != section.registration_method):
        raise ValueError("stitch registration metadata changed after server discovery")
    return record


def stitch_confirmation_fields(section, payload):
    """Bind an accepted POMC ROI receipt to the exact registered composite."""
    if section.dataset != "pomc":
        return {}
    registration = validate_stitch_registration(section)
    status = str(registration["status"])
    requires_review = bool(registration["requires_manual_review"])
    if requires_review:
        if (payload.get("registered_composite_reviewed") is not True
                or payload.get("stitch_registration_confirmation")
                != STITCH_MANUAL_CONFIRMATION):
            raise ValueError(
                "explicit image-by-image confirmation of the registered composite "
                "is required before saving ARC/ME/VMN ROIs"
            )
        confirmation = STITCH_MANUAL_CONFIRMATION
        reviewed = True
    else:
        confirmation = STITCH_AUTO_CONFIRMATIONS[status]
        reviewed = False
    return {
        "stitch_registration_schema_version": registration["schema_version"],
        "stitch_registration_file": section.registration.name,
        "stitch_registration_sha256": section.registration_sha256,
        "stitch_registration_status": status,
        "stitch_registration_method": str(registration["method"]),
        "stitch_registration_requires_manual_review": requires_review,
        "stitch_registration_confirmation": confirmation,
        "registered_composite_reviewed": reviewed,
    }


def region_absence(payload):
    """Return explicit absence flags/reasons, retaining legacy ME fields."""
    raw_flags = payload.get("region_absent", {})
    raw_reasons = payload.get("region_absence_reasons", {})
    if raw_flags is None:
        raw_flags = {}
    if raw_reasons is None:
        raw_reasons = {}
    if not isinstance(raw_flags, dict) or not isinstance(raw_reasons, dict):
        raise ValueError("region_absent and region_absence_reasons must be objects")
    flags, reasons = {}, {}
    for name in REGIONS:
        legacy_flag = payload.get("me_absent", False) if name == "ME" else False
        legacy_reason = payload.get("me_absent_reason", "") if name == "ME" else ""
        raw_flag = raw_flags.get(name, legacy_flag)
        if not isinstance(raw_flag, bool):
            raise ValueError(f"{name} absent must be true or false")
        reason = str(raw_reasons.get(name, legacy_reason) or "").strip()
        if raw_flag and not reason:
            raise ValueError(f"{name} absent requires a written anatomical reason")
        if len(reason) > 1000:
            raise ValueError(f"{name} absent reason must be 1000 characters or fewer")
        if not raw_flag and reason:
            raise ValueError(
                f"{name} absent reason is only valid when {name} absent is checked"
            )
        flags[name], reasons[name] = raw_flag, reason
    return flags, reasons


def region_status_for(labels, absent, reasons):
    status = {}
    for name in REGIONS:
        has_pixels = bool(np.any(labels == LABEL_CODES[name]))
        if has_pixels and absent[name]:
            raise ValueError(f"{name} cannot be both delineated and explicitly absent")
        status[name] = (
            "drawn" if has_pixels
            else "explicitly_absent" if absent[name]
            else "not_delineated"
        )
        if status[name] == "explicitly_absent" and not reasons[name]:
            raise ValueError(f"{name} explicit absence lacks a reason")
    return status


def accepted(output, section):
    jp, tp = paths(output, section)
    if not jp.is_file() or not tp.is_file():
        return False
    try:
        record = json.loads(jp.read_text(encoding="utf-8"))
        accepted_source = record.get("accepted_source")
        receipt_codes = record.get("label_codes")
        if not (record.get("accepted") is True and str(record.get("reviewer", "")).strip()
                and record.get("accepted_at") and record.get("image_width_px") == section.width
                and record.get("image_height_px") == section.height
                and receipt_codes in (LABEL_CODES, LEGACY_LABEL_CODES)
                and record.get("mask_sha256") == sha256(tp)
                and record.get("source_dapi_sha256") == section.dapi_sha256
                and accepted_source in {"automated_mask_unchanged", "manual_replacement_polygons"}):
            return False
        if section.dataset == "pomc":
            registration = validate_stitch_registration(section)
            status = str(registration["status"])
            requires_review = bool(registration["requires_manual_review"])
            expected_confirmation = (
                STITCH_MANUAL_CONFIRMATION if requires_review
                else STITCH_AUTO_CONFIRMATIONS[status]
            )
            if not (
                record.get("schema_version") == 4
                and record.get("stitch_registration_schema_version")
                == STITCH_REGISTRATION_SCHEMA
                and record.get("stitch_registration_sha256")
                == section.registration_sha256
                and record.get("stitch_registration_status") == status
                and record.get("stitch_registration_method")
                == registration["method"]
                and record.get("stitch_registration_requires_manual_review")
                is requires_review
                and record.get("stitch_registration_confirmation")
                == expected_confirmation
                and record.get("registered_composite_reviewed")
                is requires_review
            ):
                return False
        if accepted_source == "automated_mask_unchanged":
            if (record.get("automated_mask_confirmation") != "reviewed_image_by_image"
                    or section.auto is None
                    or record.get("automated_mask_sha256") != section.auto_sha256
                    or sha256(section.auto) != section.auto_sha256):
                return False
        validate_source_dapi(section)
        labels = validated_mask(tp, section.width, section.height, require_nonempty=True)
        if receipt_codes == LEGACY_LABEL_CODES:
            # Schema-v2 compatibility. New receipts always use the three-region
            # contract below and never fabricate an omitted region.
            me_absent = bool(record.get("me_absent", False))
            reason = str(record.get("me_absent_reason", "")).strip()
            if not np.any(labels == 1):
                return False
            return bool(
                (me_absent and reason and not np.any(labels == 2))
                or (not me_absent and np.any(labels == 2))
            )
        raw_status = record.get("region_status")
        raw_reasons = record.get("region_absence_reasons")
        if not isinstance(raw_status, dict) or not isinstance(raw_reasons, dict):
            return False
        for name in REGIONS:
            status = raw_status.get(name)
            reason = str(raw_reasons.get(name, "") or "").strip()
            has_pixels = bool(np.any(labels == LABEL_CODES[name]))
            if status not in REGION_STATES:
                return False
            if status == "drawn" and not has_pixels:
                return False
            if status != "drawn" and has_pixels:
                return False
            if status == "explicitly_absent" and not reason:
                return False
            if status != "explicitly_absent" and reason:
                return False
        return True
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return False


def orient(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])


def on_segment(a, b, p):
    return min(a[0], b[0])-1e-9 <= p[0] <= max(a[0], b[0])+1e-9 and min(a[1], b[1])-1e-9 <= p[1] <= max(a[1], b[1])+1e-9 and abs(orient(a, b, p)) <= 1e-9


def intersects(a, b, c, d):
    o1, o2, o3, o4 = orient(a,b,c), orient(a,b,d), orient(c,d,a), orient(c,d,b)
    if ((o1 > 1e-9 > o2) or (o2 > 1e-9 > o1)) and ((o3 > 1e-9 > o4) or (o4 > 1e-9 > o3)):
        return True
    return ((abs(o1)<=1e-9 and on_segment(a,b,c)) or (abs(o2)<=1e-9 and on_segment(a,b,d))
            or (abs(o3)<=1e-9 and on_segment(c,d,a)) or (abs(o4)<=1e-9 and on_segment(c,d,b)))


def valid_polygon(raw, name, index, width, height):
    if isinstance(raw, dict):
        raw = raw.get("vertices_px")
    if not isinstance(raw, list) or len(raw) > 10000:
        raise ValueError(f"{name} polygon {index+1} is not a valid vertex list")
    points = []
    for j, vertex in enumerate(raw):
        if not isinstance(vertex, (list, tuple)) or len(vertex) != 2:
            raise ValueError(f"{name} polygon {index+1} vertex {j+1} must be [x,y]")
        try:
            x, y = float(vertex[0]), float(vertex[1])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{name} polygon {index+1} has a nonnumeric vertex") from exc
        if not math.isfinite(x) or not math.isfinite(y) or x < 0 or y < 0 or x > width-1 or y > height-1:
            raise ValueError(f"{name} polygon {index+1} has an out-of-bounds/nonfinite vertex")
        points.append((x, y))
    if len(points) >= 2 and points[0] == points[-1]:
        points.pop()
    if len(points) < 3 or len(set(points)) < 3:
        raise ValueError(f"{name} polygon {index+1} needs at least 3 distinct vertices")
    n = len(points)
    for i in range(n):
        a, b = points[i], points[(i+1)%n]
        if a == b:
            raise ValueError(f"{name} polygon {index+1} has a zero-length edge")
        for j in range(i+1, n):
            if j == i or j == (i+1)%n or (i == 0 and j == n-1):
                continue
            if intersects(a, b, points[j], points[(j+1)%n]):
                raise ValueError(f"{name} polygon {index+1} self-intersects")
    area2 = abs(sum(points[i][0]*points[(i+1)%n][1]-points[(i+1)%n][0]*points[i][1] for i in range(n)))
    if area2 < 8:
        raise ValueError(f"{name} polygon {index+1} is too small or degenerate")
    return [[round(x,3), round(y,3)] for x,y in points]


def validate(payload, section):
    raw = payload.get("polygons") if isinstance(payload, dict) else None
    if not isinstance(raw, dict):
        raise ValueError("polygons must contain ARC, ME, and/or VMN lists")
    absent, absence_reasons = region_absence(payload)
    polygons, masks = {}, {}
    for name in REGIONS:
        source = raw.get(name, [])
        if not isinstance(source, list):
            raise ValueError(f"{name} polygons must be a list")
        if source and absent[name]:
            raise ValueError(f"{name} polygons must be empty when {name} absent is documented")
        polygons[name] = [valid_polygon(p, name, i, section.width, section.height) for i,p in enumerate(source)]
        mask = np.zeros((section.height, section.width), bool)
        for poly in polygons[name]:
            v = np.asarray(poly)
            rr, cc = draw_polygon(v[:,1], v[:,0], shape=mask.shape)
            if len(rr) < 3:
                raise ValueError(f"a {name} polygon covers no valid region")
            mask[rr,cc] = True
        masks[name] = mask
    for index, first in enumerate(REGIONS):
        for second in REGIONS[index + 1:]:
            overlap = int(np.count_nonzero(masks[first] & masks[second]))
            if overlap:
                raise ValueError(f"{first} and {second} overlap by {overlap:,} pixels")
    labels = np.zeros((section.height, section.width), np.uint8)
    for name in REGIONS:
        labels[masks[name]] = LABEL_CODES[name]
    if not np.any(labels):
        raise ValueError("draw at least one closed ARC, ME, or VMN polygon")
    status = region_status_for(labels, absent, absence_reasons)
    return polygons, labels, status, absence_reasons


def atomic_pair(json_path, tif_path, record, labels):
    json_path.parent.mkdir(parents=True, exist_ok=True)
    token = uuid.uuid4().hex
    jt = json_path.parent / f".{json_path.name}.{token}.tmp"
    tt = tif_path.parent / f".{tif_path.name}.{token}.tmp"
    try:
        buffer = io.BytesIO()
        tifffile.imwrite(buffer, labels, photometric="minisblack", metadata=None)
        tif_bytes = buffer.getvalue()
        record["mask_sha256"] = hashlib.sha256(tif_bytes).hexdigest()
        json_bytes = (json.dumps(record, indent=2, sort_keys=True)+"\n").encode()
        for target, data in ((tt,tif_bytes),(jt,json_bytes)):
            with target.open("xb") as handle:
                handle.write(data); handle.flush(); os.fsync(handle.fileno())
        os.replace(tt, tif_path)       # JSON is the commit marker; strict consumers verify the hash.
        os.replace(jt, json_path)
        try:
            fd = os.open(json_path.parent, os.O_RDONLY)
            try: os.fsync(fd)
            finally: os.close(fd)
        except OSError:
            pass
    finally:
        jt.unlink(missing_ok=True); tt.unlink(missing_ok=True)


class App:
    def __init__(self, dataset, raw, analysis, output, allow_replace=True, *,
                 figure_command=None, figure_output_root=None, figure_button=True,
                 figure_analysis=None, figure_input_root=None,
                 figure_prepared_dir=None, figure_analyze=True):
        self.dataset = str(dataset)
        self.raw, self.analysis, self.output = raw.resolve(), analysis.resolve(), output.resolve()
        # Re-segmenting an already accepted section is the normal reason to reopen
        # this reviewer, so replacement is permitted by default. Every replacement
        # still rewrites the durable JSON/TIFF pair atomically and rebinds the
        # reviewer, timestamp and source hashes; nothing is ever deleted.
        self.allow_replace = bool(allow_replace)
        self.sections = discover(self.dataset, self.raw, self.analysis)
        self.index = {s.key:s for s in self.sections}
        if len(self.index) != len(self.sections):
            raise ValueError("discovery produced duplicate animal/section keys")
        self.lock = threading.Lock()
        self.figure_label = FIGURE_LABELS.get(self.dataset, "the figure")
        self.figure_button = bool(figure_button)
        self.figure_output_root = Path(figure_output_root or "/tmp").expanduser().resolve()
        # A completed analysis directory if one was named, otherwise None: the
        # button then runs the figure's analyzer to produce one. The review has
        # to finish before the analysis can run, so no reviewer launch can be
        # expected to already know a completed analysis directory.
        self.figure_analysis = (
            Path(figure_analysis).expanduser().resolve() if figure_analysis else None
        )
        # Where the reviewed reconstructed sections live, which the analyzer
        # validates every accepted receipt against. It is this reviewer's own
        # --analysis-dir unless the launcher says otherwise.
        self.figure_prepared = (
            Path(figure_prepared_dir).expanduser().resolve()
            if figure_prepared_dir else self.analysis
        )
        self.figure_input_root = (
            Path(figure_input_root).expanduser().resolve() if figure_input_root else None
        )
        self.figure_analyze = bool(figure_analyze)
        self.figure_command = list(figure_command) if figure_command else None
        self.build_lock = threading.Lock()
        self._auto_analysis_lock = threading.Lock()
        self._auto_analysis_dir = None
        self.build = {
            "running": False, "started_at": None, "finished_at": None,
            "ok": None, "stage": "idle", "error": "",
            "machine_defaults_accepted": 0, "machine_defaults_skipped": [],
            "output_dir": "", "command": "", "stdout_tail": "", "stderr_tail": "",
            "analysis_dir": "", "steps": [], "figure_files": [],
        }

    def section(self, query):
        key = (query.get("animal",[""])[0], query.get("section",[""])[0])
        if key not in self.index:
            raise KeyError(f"unknown animal/section: {key[0]}/{key[1]}")
        return self.index[key]

    def rows(self):
        rows = []
        for s in self.sections:
            jp, tp = paths(self.output, s)
            valid = accepted(self.output, s)
            row = {
                "dataset":s.dataset,"animal":s.animal,"section":s.section,"width":s.width,"height":s.height,
                "auto_mask_available":s.auto is not None,"auto_mask_kind":s.auto_kind,
                "auto_mask_sha256":s.auto_sha256,"accepted":valid,
                "receipt_files_present":jp.exists() or tp.exists(),
                "invalid_existing_pair":(jp.exists() or tp.exists()) and not valid,
                "replacement_allowed":self.allow_replace,
                "stitch_registration_valid":s.registration_valid,
                "stitch_registration_status":s.registration_status,
                "stitch_registration_method":s.registration_method,
                "stitch_registration_requires_manual_review":s.registration_requires_manual_review,
                "stitch_registration_sha256":s.registration_sha256,
                "stitch_registration_error":s.registration_error,
                "stitch_registration_qc":dict(s.registration_qc or {}),
            }
            if s.dataset == "pomc" and s.registration_valid:
                try:
                    current = validate_stitch_registration(s)
                    row["stitch_registration_transform"] = dict(
                        current.get("transform", {})
                    )
                except (OSError, ValueError, TypeError) as exc:
                    row["stitch_registration_valid"] = False
                    row["stitch_registration_error"] = str(exc)
            rows.append(row)
        return rows

    def config(self):
        plan = self.figure_build_plan()
        return {
            "dataset": self.dataset, "figure_label": self.figure_label,
            "figure_button_enabled": self.figure_button,
            "figure_command": shlex.join(plan["render_command"]),
            "figure_analysis_command": (
                shlex.join(plan["analysis_command"]) if plan["analysis_command"] else ""
            ),
            "machine_default_attestation": MACHINE_DEFAULT_ATTESTATION,
            "replacement_allowed": self.allow_replace,
            "output_dir": str(self.output),
            "figure_analysis_dir": str(plan["analysis_dir"]),
            "figure_analysis_reused": bool(plan["analysis_reused"]),
            "figure_prepared_dir": str(self.figure_prepared),
            "figure_output_root": str(self.figure_output_root),
            "sections": len(self.sections),
            **self.source_state(),
        }

    def figure_command_for(self, target, analysis_dir=None):
        """Command that renders the figure from ``analysis_dir`` into ``target``."""
        if self.figure_command:
            return [
                part.replace(FIGURE_OUTPUT_TOKEN, str(target))
                for part in self.figure_command
            ]
        analysis = str(analysis_dir if analysis_dir is not None else FIGURE_ANALYSIS_TOKEN)
        if self.dataset == "pomc":
            return [
                sys.executable, str(ROOT / "scripts" / "Fig4" / "03_make_figure_4_pomc_cfos.py"),
                "--analysis-dir", analysis,
                "--manual-region-dir", str(self.output), "--require-manual-regions",
                "--outdir", str(target), "--figure-name", FIGURE_MASTER_NAME["pomc"],
                "--scalebar-um", "100",
                "--inset-field-um", "50", "--inset-scalebar-um", "20",
                "--um-per-px", "0.755", "--dpi", "600",
            ]
        return [
            sys.executable, str(ROOT / "scripts" / "Fig5" / "02_make_figure_5_gfap_iba1_microglia.py"),
            "--analysis-dir", analysis, "--input-root", str(self.raw),
            "--manual-region-dir", str(self.output), "--output-dir", str(target),
            "--um-per-px", "1.51", "--dpi", "600", "--formats", "pdf,png",
        ]

    def analysis_command_for(self, target):
        """Command that turns the accepted receipts into an analysis directory.

        This is the step the review used to dead-end on: the renderer needs an
        analysis run over the accepted masks, and only the figure's analyzer can
        produce one. It reads the prepared review images, never writes to them,
        and writes into ``target``, which the caller guarantees is absent.
        """
        analyzer = ROOT / FIGURE_ANALYZERS[self.dataset]
        command = [
            sys.executable, str(analyzer),
            "--output", str(target),
            "--manual-region-dir", str(self.output), "--require-manual-regions",
            "--no-auto-conda",
        ]
        if self.figure_input_root is not None:
            command += ["--input", str(self.figure_input_root)]
        if self.dataset == "pomc":
            command += ["--prepared-review-dir", str(self.figure_prepared)]
        return command

    def figure_build_plan(self):
        """The one or two commands the Make Figure button will run, in order."""
        analysis_dir = self.figure_analysis
        reused = bool(analysis_dir is not None and not self.missing_figure_analysis_inputs())
        analysis_command = None
        if self.figure_command:
            # A caller-supplied command owns the whole build; nothing is inferred.
            analysis_dir = analysis_dir or self.figure_prepared
        elif not reused:
            if analysis_dir is None:
                analysis_dir = self.auto_analysis_dir()
            analysis_command = self.analysis_command_for(analysis_dir)
        render_dir = self.fresh_output_dir(
            FIGURE_OUTPUT_PREFIX.get(self.dataset, "paper_figure_render")
        )
        return {
            "analysis_dir": analysis_dir,
            "analysis_reused": reused,
            "analysis_command": analysis_command,
            "render_dir": render_dir,
            "render_command": self.figure_command_for(render_dir, analysis_dir),
        }

    def auto_analysis_dir(self):
        """One stable, absent analysis target for this reviewer session.

        It is chosen once so the directory the page names is the directory the
        build writes, and it is never an existing path: analysis output is
        generated, and this reviewer replaces nothing.
        """
        with self._auto_analysis_lock:
            if self._auto_analysis_dir is None:
                preferred = FIGURE_DEFAULT_ANALYSIS_OUTPUT.get(self.dataset)
                if preferred is not None and not preferred.exists():
                    self._auto_analysis_dir = preferred.resolve()
                else:
                    self._auto_analysis_dir = self.fresh_output_dir(
                        FIGURE_ANALYSIS_PREFIX.get(self.dataset, "paper_figure_analysis"),
                        root=FIGURE_FALLBACK_ANALYSIS_ROOT.get(self.dataset),
                    )
            return self._auto_analysis_dir

    def fresh_output_dir(self, prefix, root=None):
        """An absolute, absent target under ``root``, default the figure output root."""
        parent = Path(root) if root is not None else self.figure_output_root
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        while True:
            target = parent / f"{prefix}_{stamp}_{uuid.uuid4().hex[:8]}"
            if not target.exists():
                return target

    def fresh_figure_output_dir(self):
        """An absolute, absent render target; existing outputs are never reused."""
        return self.fresh_output_dir(
            FIGURE_OUTPUT_PREFIX.get(self.dataset, "paper_figure_render")
        )

    def machine_default_candidates(self):
        """Split pending sections into batch-acceptable and blocked."""
        ready, blocked = [], []
        for section in self.sections:
            if accepted(self.output, section):
                continue
            if section.auto is None or not section.auto_kind.startswith("native_"):
                blocked.append((section, "no native automated mask to accept"))
                continue
            if section.dataset == "pomc":
                try:
                    registration = validate_stitch_registration(section)
                except (OSError, ValueError, TypeError) as exc:
                    blocked.append((section, f"invalid stitch registration: {exc}"))
                    continue
                if bool(registration["requires_manual_review"]):
                    blocked.append((
                        section,
                        "registered composite still needs image-by-image seam confirmation",
                    ))
                    continue
            ready.append(section)
        return ready, blocked

    def _accept_native_automated(self, section, reviewer, session, extra):
        """Write one accepted receipt that keeps the native automated mask."""
        validate_source_dapi(section)
        current_auto_hash = sha256(section.auto)
        if current_auto_hash != section.auto_sha256:
            raise ValueError("native automated mask changed after server discovery")
        labels = validated_mask(section.auto, section.width, section.height, require_nonempty=True)
        absent = {name: False for name in REGIONS}
        reasons = {name: "" for name in REGIONS}
        region_status = region_status_for(labels, absent, reasons)
        now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")
        registration_extra = stitch_confirmation_fields(section, {})
        registration_extra.update(extra)
        record, jp, tp = self._record(
            section, reviewer, session, now, "automated_mask_unchanged", labels,
            region_status, reasons, {name: [] for name in REGIONS}, registration_extra,
        )
        return record, jp, tp, now, current_auto_hash

    def accept_machine_defaults(self, reviewer, session, attestation):
        """Accept the native automated mask of every remaining eligible section.

        The receipts are ordinary ``automated_mask_unchanged`` receipts, so they
        are exactly as strong as the per-section button. They additionally record
        that the confirmation was given once for a batch, and the literal wording
        the reviewer attested to, so an auditor can always tell the two apart.
        """
        if str(attestation).strip().upper() != MACHINE_DEFAULT_ATTESTATION:
            raise ValueError(
                "machine default acceptance requires the exact attestation "
                f"{MACHINE_DEFAULT_ATTESTATION!r}"
            )
        ready, blocked = self.machine_default_candidates()
        batch_id = uuid.uuid4().hex
        extra_base = {
            "automated_mask_confirmation": "reviewed_image_by_image",
            "machine_default_batch": True,
            "machine_default_batch_id": batch_id,
            "machine_default_batch_size": len(ready),
            "machine_default_attestation": MACHINE_DEFAULT_ATTESTATION,
            "review_granularity": "batch_machine_default",
        }
        written, failed = [], []
        for section in ready:
            try:
                _record, jp, _tp, _now, _hash = self._accept_native_automated(
                    section, reviewer, session, dict(extra_base),
                )
                written.append(f"{section.animal}/{section.section}")
            except Exception as exc:
                failed.append(f"{section.animal}/{section.section}: {exc}")
        skipped = [f"{s.animal}/{s.section}: {why}" for s, why in blocked] + failed
        return {
            "batch_id": batch_id, "accepted": written, "skipped": skipped,
        }

    def missing_figure_analysis_inputs(self, analysis_dir=None):
        """Required analysis artifacts the renderer cannot run without."""
        root = self.figure_analysis if analysis_dir is None else Path(analysis_dir)
        required = FIGURE_ANALYSIS_REQUIREMENTS.get(self.dataset, ())
        if root is None:
            return list(required)
        return [name for name in required if not (root / name).exists()]

    def require_figure_analysis(self):
        """Refuse only when this session cannot reach a completed analysis."""
        if self.figure_command:
            return
        missing = self.missing_figure_analysis_inputs()
        if not missing:
            return
        analyzer = FIGURE_ANALYZERS.get(self.dataset, "the figure analyzer")
        if self.figure_analysis is not None and self.figure_analysis.exists():
            raise ValueError(
                f"{self.figure_analysis} exists but is not a completed analysis "
                f"directory: it is missing {', '.join(missing)}. Nothing is written "
                "over a directory that already holds something, so either delete it "
                "yourself, name an absent --figure-analysis-dir, or name a completed "
                f"one produced by {analyzer}."
            )
        if not self.figure_analyze:
            raise ValueError(
                "this session was started with --no-figure-analyze, so the button "
                f"only renders from a completed analysis directory. Run {analyzer} "
                "and relaunch with --figure-analysis-dir pointing at its output."
            )
        if not (ROOT / FIGURE_ANALYZERS[self.dataset]).is_file():
            raise ValueError(
                f"the figure analyzer is missing: {ROOT / FIGURE_ANALYZERS[self.dataset]}"
            )

    def start_figure_build(self, payload):
        """Queue one background render of Figure 4/6 from the accepted receipts."""
        if not self.figure_button:
            raise ValueError("the figure build button is disabled for this session")
        self.require_figure_analysis()
        reviewer, session = self._identity(payload)
        machine_defaults = bool(payload.get("machine_defaults"))
        attestation = str(payload.get("machine_default_attestation", ""))
        if machine_defaults and str(attestation).strip().upper() != MACHINE_DEFAULT_ATTESTATION:
            raise ValueError(
                "machine default acceptance requires the exact attestation "
                f"{MACHINE_DEFAULT_ATTESTATION!r}"
            )
        plan = self.figure_build_plan()
        target = plan["render_dir"]
        command = plan["render_command"]
        steps = []
        if plan["analysis_command"]:
            steps.append({
                "name": "analyzing the accepted segmentations",
                "command": plan["analysis_command"],
                "output_dir": str(plan["analysis_dir"]),
            })
        steps.append({
            "name": f"rendering {self.figure_label}",
            "command": command,
            "output_dir": str(target),
        })
        with self.build_lock:
            if self.build["running"]:
                raise ValueError("a figure build is already running in this session")
            self.build = {
                "running": True,
                "started_at": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z"),
                "finished_at": None, "ok": None,
                "stage": ("accepting machine defaults" if machine_defaults
                          else steps[0]["name"]),
                "error": "", "machine_defaults_accepted": 0, "machine_defaults_skipped": [],
                "output_dir": str(target), "command": shlex.join(command),
                "stdout_tail": "", "stderr_tail": "", "figure_files": [],
                "analysis_dir": str(plan["analysis_dir"]),
                "steps": [
                    {"name": step["name"], "command": shlex.join(step["command"]),
                     "output_dir": step["output_dir"]}
                    for step in steps
                ],
            }
        thread = threading.Thread(
            target=self._run_figure_build,
            args=(reviewer, session, machine_defaults, attestation, steps),
            daemon=True, name="make-figure",
        )
        thread.start()
        return {"ok": True, "started": True, "figure_label": self.figure_label,
                "output_dir": str(target), "command": shlex.join(command),
                "analysis_dir": str(plan["analysis_dir"]),
                "steps": [step["name"] for step in steps]}

    def figure_master_files(self, target):
        """The master graphics a finished render actually wrote, by full path.

        Naming the directory alone is what made a successful build look like
        nothing: the reviewer has no way to know which of its files is the
        figure. These are the paths to open.
        """
        target = Path(target)
        stem = FIGURE_MASTER_NAME.get(self.dataset)
        found = []
        if stem:
            found = [
                str(target / f"{stem}{suffix}") for suffix in FIGURE_MASTER_SUFFIXES
                if (target / f"{stem}{suffix}").is_file()
            ]
        if not found and target.is_dir():
            for suffix in FIGURE_MASTER_SUFFIXES:
                found += sorted(str(path) for path in target.glob(f"*{suffix}"))
        return found

    def source_state(self):
        """Whether this process is still running the reviewer that is on disk.

        A server keeps the code it was started with for its whole life, so a
        tab left open across a fix silently drives the old commands. That is
        indistinguishable from the fix not working, so the page has to be told.
        """
        current = source_sha256()
        return {
            "source_path": str(HERE),
            "source_sha256": SOURCE_SHA256_AT_START,
            "source_sha256_on_disk": current,
            "source_stale": bool(
                current and SOURCE_SHA256_AT_START and current != SOURCE_SHA256_AT_START
            ),
        }

    def _release_failed_auto_analysis(self):
        """Let a retry pick a new target when this session's own analysis failed.

        A half-written analysis directory would otherwise be handed to the next
        attempt, which correctly refuses to write into an existing directory. A
        completed one is kept so a render-only retry does not repeat the
        analysis. Nothing is deleted either way.
        """
        if self.figure_analysis is not None:
            return
        with self._auto_analysis_lock:
            target = self._auto_analysis_dir
            if target is not None and self.missing_figure_analysis_inputs(target):
                self._auto_analysis_dir = None

    def _finish_build(self, **fields):
        with self.build_lock:
            self.build.update(fields)
            self.build["running"] = False
            self.build["finished_at"] = (
                datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")
            )

    def _run_figure_build(self, reviewer, session, machine_defaults, attestation, steps):
        """Run the queued steps in order: analyze the accepted masks, then render.

        The analysis step is what the review used to stop short of. It reads the
        prepared review images and the accepted receipts and writes a completed
        analysis directory; only then can the renderer draw the figure.
        """
        try:
            if machine_defaults:
                batch = self.accept_machine_defaults(reviewer, session, attestation)
                with self.build_lock:
                    self.build["machine_defaults_accepted"] = len(batch["accepted"])
                    self.build["machine_defaults_skipped"] = batch["skipped"]
                    self.build["stage"] = steps[0]["name"]
            pending = [
                f"{s.animal}/{s.section}" for s in self.sections
                if not accepted(self.output, s)
            ]
            if pending:
                self._finish_build(
                    ok=False, stage="blocked",
                    error=(
                        f"{len(pending)} of {len(self.sections)} sections are still not "
                        f"accepted, so {self.figure_label} was not rendered: "
                        + ", ".join(pending[:8])
                        + (" ..." if len(pending) > 8 else "")
                    ),
                )
                return
            for index, step in enumerate(steps):
                with self.build_lock:
                    self.build["stage"] = step["name"]
                    self.build["command"] = shlex.join(step["command"])
                completed = subprocess.run(
                    [str(part) for part in step["command"]], cwd=str(ROOT), text=True,
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                    check=False,
                )
                stdout_tail = (completed.stdout or "")[-LOG_TAIL_CHARS:]
                stderr_tail = (completed.stderr or "")[-LOG_TAIL_CHARS:]
                if completed.returncode:
                    self._release_failed_auto_analysis()
                    self._finish_build(
                        ok=False, stage="failed",
                        error=(
                            f"step {index + 1} of {len(steps)}, {step['name']}, exited "
                            f"with status {completed.returncode}; see the log below"
                        ),
                        stdout_tail=stdout_tail, stderr_tail=stderr_tail,
                    )
                    return
            self._finish_build(
                ok=True, stage="finished", error="",
                stdout_tail=stdout_tail, stderr_tail=stderr_tail,
                figure_files=self.figure_master_files(steps[-1]["output_dir"]),
            )
        except Exception as exc:
            self._release_failed_auto_analysis()
            self._finish_build(ok=False, stage="failed", error=str(exc))

    def build_status(self):
        with self.build_lock:
            state = dict(self.build)
        state["figure_label"] = self.figure_label
        state["accepted"] = sum(1 for s in self.sections if accepted(self.output, s))
        state["total"] = len(self.sections)
        state.update(self.source_state())
        return state

    def _identity(self, payload):
        reviewer, session = str(payload.get("reviewer","")).strip(), str(payload.get("session","")).strip()
        if not reviewer or not session:
            raise ValueError("reviewer name and session identifier are required")
        if len(reviewer)>200 or len(session)>200:
            raise ValueError("reviewer/session must be 200 characters or fewer")
        return reviewer, session

    def _record(self, section, reviewer, session, now, accepted_source, labels,
                region_status, absence_reasons, polygons, extra=None):
        jp, tp = paths(self.output, section)
        source_dapi_hash = validate_source_dapi(section)
        source_match = re.match(r"^(?P<sample>.+)_S\d+(?:_T.*)?$", section.section, flags=re.IGNORECASE)
        source_sample = source_match.group("sample") if section.dataset == "pomc" and source_match else None
        record = {
            "schema_version":4 if section.dataset == "pomc" else 3,
            "accepted":True,"reviewer":reviewer,"session":session,"accepted_at":now,
            "accepted_source":accepted_source,"dataset":section.dataset,
            "animal":section.animal,"section":section.section,"source_sample":source_sample,
            "image_width_px":section.width,
            "image_height_px":section.height,"coordinate_system":"original DAPI pixels; top-left origin; x right; y down",
            "label_codes":LABEL_CODES,
            "region_status":dict(region_status),
            "region_presence":{
                name: True if region_status[name] == "drawn" else
                False if region_status[name] == "explicitly_absent" else None
                for name in REGIONS
            },
            "region_absence_reasons":dict(absence_reasons),
            "delineated_regions":[name for name in REGIONS if region_status[name] == "drawn"],
            "omitted_regions":[name for name in REGIONS if region_status[name] == "not_delineated"],
            "explicitly_absent_regions":[name for name in REGIONS if region_status[name] == "explicitly_absent"],
            "me_absent":region_status["ME"] == "explicitly_absent",
            "me_absent_reason":absence_reasons["ME"],
            "polygons":{name:[{"closed":True,"vertices_px":p} for p in polygons.get(name,[])] for name in REGIONS},
            "source_paths":{
                "dapi":str(section.dapi),
                "automated_mask":str(section.auto) if section.auto else None,
                **({"stitch_registration":str(section.registration)}
                   if section.dataset == "pomc" else {}),
            },
            "source_dapi_sha256":source_dapi_hash,
            "automated_mask_sha256":section.auto_sha256,
            "automated_mask_provenance":{"dataset":section.dataset,"kind":section.auto_kind},
            "mask_file":tp.name,"mask_sha256":"",
        }
        record.update(dict(extra or {}))
        with self.lock:
            if (jp.exists() or tp.exists()) and not self.allow_replace:
                state = "a valid accepted pair" if accepted(self.output, section) else "pre-existing output files"
                raise FileExistsError(
                    f"{section.animal}/{section.section} already has {state}; "
                    "this reviewer was started with --lock-accepted, so restart it "
                    "without that flag to segment the section again"
                )
            atomic_pair(jp,tp,record,labels)
            if not accepted(self.output, section):
                raise OSError(
                    "paired annotation failed immediate disk reload/hash validation"
                )
        return record, jp, tp

    def save(self, payload):
        if not isinstance(payload, dict):
            raise ValueError("request must be a JSON object")
        section = self.index.get((str(payload.get("animal","")),str(payload.get("section",""))))
        if section is None:
            raise ValueError("unknown animal/section")
        reviewer, session = self._identity(payload)
        polygons, labels, region_status, absence_reasons = validate(payload, section)
        now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")
        registration_extra = stitch_confirmation_fields(section, payload)
        record, jp, tp = self._record(section,reviewer,session,now,"manual_replacement_polygons",
                                      labels,region_status,absence_reasons,polygons,
                                      registration_extra)
        return {
            "ok":True,"accepted":True,"pair_valid":True,"accepted_at":now,
            "mask_sha256":record["mask_sha256"],"region_status":region_status,
            "json_path":str(jp),"mask_path":str(tp),
        }

    def accept_automated_mask(self, payload):
        if not isinstance(payload, dict):
            raise ValueError("request must be a JSON object")
        section = self.index.get((str(payload.get("animal","")),str(payload.get("section",""))))
        if section is None:
            raise ValueError("unknown animal/section")
        reviewer, session = self._identity(payload)
        if (payload.get("confirm_automated_mask_unchanged") is not True
                or payload.get("automated_mask_confirmation") != "reviewed_image_by_image"):
            raise ValueError("explicit image-by-image native automated-mask confirmation is required")
        if section.auto is None:
            raise ValueError("this native section has no automated mask to accept")
        if section.dataset not in {"gfap", "pomc"} or not section.auto_kind.startswith("native_"):
            raise ValueError("only a native GFAP or POMC automated mask can be accepted")
        validate_source_dapi(section)
        current_auto_hash = sha256(section.auto)
        if current_auto_hash != section.auto_sha256:
            raise ValueError("native automated mask changed after server discovery")
        labels = validated_mask(section.auto, section.width, section.height, require_nonempty=True)
        absent, absence_reasons = region_absence(payload)
        region_status = region_status_for(labels, absent, absence_reasons)
        accepted_source = "automated_mask_unchanged"
        now = datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z")
        registration_extra = stitch_confirmation_fields(section, payload)
        registration_extra["automated_mask_confirmation"] = "reviewed_image_by_image"
        record, jp, tp = self._record(section,reviewer,session,now,accepted_source,labels,
                                      region_status,absence_reasons,{name:[] for name in REGIONS},
                                      registration_extra)
        return {"ok":True,"accepted":True,"pair_valid":True,"accepted_at":now,"accepted_source":accepted_source,
                "automated_mask_sha256":current_auto_hash,"mask_sha256":record["mask_sha256"],
                "region_status":region_status,"json_path":str(jp),"mask_path":str(tp)}


class Server(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True
    def __init__(self,address,app):
        super().__init__(address,Handler); self.app=app


class Handler(BaseHTTPRequestHandler):
    server_version="ARCMEVMNAnnotator/3.0"
    @property
    def app(self): return self.server.app
    def log_message(self,fmt,*args): print(f"[{self.log_date_time_string()}] {self.client_address[0]} {fmt%args}")
    def send_bytes(self,data,ctype,status=200):
        self.send_response(status); self.send_header("Content-Type",ctype); self.send_header("Content-Length",str(len(data)))
        self.send_header("Cache-Control","no-store"); self.send_header("X-Content-Type-Options","nosniff")
        self.send_header("Content-Security-Policy","default-src 'self'; img-src 'self' blob: data:; style-src 'unsafe-inline'; script-src 'unsafe-inline'; connect-src 'self'")
        self.end_headers(); self.wfile.write(data)
    def send_json(self,obj,status=200): self.send_bytes((json.dumps(obj,separators=(',',':'))+'\n').encode(),"application/json; charset=utf-8",status)
    def error_json(self,status,message): self.send_json({"ok":False,"error":message},status)
    def do_GET(self):
        parsed=urlparse(self.path); q=parse_qs(parsed.query)
        try:
            if parsed.path=="/": self.send_bytes(HTML.encode(),"text/html; charset=utf-8")
            elif parsed.path in ("/health", "/api/health"):
                rows=self.app.rows()
                self.send_json({"ok":True,"status":"healthy","service":"arc-me-annotator","sections":len(rows),"accepted":sum(bool(x["accepted"]) for x in rows)})
            elif parsed.path=="/api/config":
                self.send_json({"ok":True,**self.app.config()})
            elif parsed.path=="/api/figure-status":
                self.send_json({"ok":True,**self.app.build_status()})
            elif parsed.path=="/api/sections":
                rows=self.app.rows(); self.send_json({"sections":rows,"accepted":sum(bool(x["accepted"]) for x in rows)})
            elif parsed.path=="/api/image":
                s=self.app.section(q); st=s.dapi.stat(); self.send_bytes(dapi_png(str(s.dapi),st.st_mtime_ns),"image/png")
            elif parsed.path=="/api/auto-mask":
                s=self.app.section(q)
                if s.auto is None: self.error_json(404,"no automated ARC/ME mask available"); return
                st=s.auto.stat(); self.send_bytes(auto_png(str(s.auto),st.st_mtime_ns,s.width,s.height),"image/png")
            elif parsed.path=="/api/annotation":
                s=self.app.section(q); jp,tp=paths(self.app.output,s)
                if not jp.is_file(): self.send_json({"found":False,"animal":s.animal,"section":s.section,"polygons":{"ARC":[],"ME":[]}}); return
                record=json.loads(jp.read_text(encoding="utf-8")); record["found"]=True; record["paired_mask_present"]=tp.is_file(); record["pair_valid"]=accepted(self.app.output,s); self.send_json(record)
            elif parsed.path=="/favicon.ico": self.send_bytes(b"","image/x-icon",204)
            else: self.error_json(404,"endpoint not found")
        except KeyError as exc: self.error_json(404,str(exc))
        except Exception as exc: self.error_json(500,str(exc))
    def do_POST(self):
        endpoint=urlparse(self.path).path
        if endpoint not in ("/api/annotation","/api/accept-automated-mask","/api/make-figure"): self.error_json(404,"endpoint not found"); return
        try:
            if self.headers.get_content_type()!="application/json": self.error_json(415,"Content-Type must be application/json"); return
            length=int(self.headers.get("Content-Length","0"))
            if length<=0 or length>MAX_BODY: self.error_json(413,"invalid or oversized request body"); return
            payload=json.loads(self.rfile.read(length).decode())
            if endpoint=="/api/make-figure": self.send_json(self.app.start_figure_build(payload))
            elif endpoint=="/api/accept-automated-mask": self.send_json(self.app.accept_automated_mask(payload))
            else: self.send_json(self.app.save(payload))
        except (ValueError,TypeError,json.JSONDecodeError) as exc: self.error_json(422,str(exc))
        except Exception as exc: self.error_json(500,str(exc))


HTML=r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>ARC/ME delineation</title><style>
:root{--ink:#18272d;--muted:#64737a;--paper:#f5f6f3;--rule:#d1d8d6;--arc:#00adb5;--me:#e67e22;--ok:#287a52;--bad:#a63d32;--canvas:#101619}*{box-sizing:border-box}html,body{height:100%;margin:0}body{font:14px/1.35 "Helvetica Neue",sans-serif;color:var(--ink);background:var(--paper);overflow:hidden}button,input,select,textarea{font:inherit;color:inherit}button{border:1px solid #abb6b5;background:#fff;border-radius:2px;padding:7px 9px;cursor:pointer;font-weight:650}button:hover{background:#edf1ef}button:disabled{opacity:.45}.app{height:100%;display:grid;grid-template-rows:54px minmax(0,1fr)}header{display:flex;align-items:center;gap:18px;padding:0 18px;background:#fff;border-bottom:1px solid var(--rule)}h1{font:700 18px Georgia,serif;margin:0}.sub{font-size:12px;color:var(--muted)}#progress{margin-left:auto;font-weight:700}.main{min-height:0;display:grid;grid-template-columns:306px minmax(0,1fr)}aside{overflow:auto;padding:14px;background:#fff;border-right:1px solid var(--rule)}.group{border-top:1px solid var(--rule);padding:12px 0}.group:first-child{border:0;padding-top:0}.title{font-size:11px;text-transform:uppercase;letter-spacing:.1em;color:var(--muted);font-weight:800;margin-bottom:7px}.field{display:block;margin:7px 0 3px;font-size:12px;font-weight:700}input[type=text],select,textarea{width:100%;border:1px solid #abb6b5;border-radius:2px;padding:7px}textarea{min-height:58px;resize:vertical}.grid{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin-top:7px}.class{border-width:2px}.class.arc.active{border-color:var(--arc);background:#e5f7f7}.class.me.active{border-color:var(--me);background:#fff2e6}.counts{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin-top:8px;font-size:12px}.count{padding:6px;border-left:4px solid}.count.arc{border-color:var(--arc);background:#edf8f8}.count.me{border-color:var(--me);background:#fff5eb}.range{display:grid;grid-template-columns:72px 1fr 38px;gap:6px;align-items:center;margin:8px 0;font-size:12px}.range output{text-align:right}.check{display:flex;gap:7px;align-items:center}.accept{width:100%;padding:10px;background:var(--ok);border-color:var(--ok);color:white}.autoaccept{width:100%;padding:10px;margin-bottom:8px;background:#0e6573;border-color:#0e6573;color:white}.warning{font-size:11px;color:var(--bad);margin:5px 0 9px}.status{font-size:12px;min-height:34px;margin-top:8px;color:var(--muted);overflow-wrap:anywhere}.status.error,.dirty{color:var(--bad);font-weight:800}.status.ok,.accepted{color:var(--ok);font-weight:800}.workspace{min-width:0;min-height:0;display:grid;grid-template-rows:40px minmax(0,1fr);background:var(--canvas)}.toolbar{display:flex;align-items:center;gap:7px;padding:5px 10px;background:#edf0ee}.toolbar .name{margin-right:auto;font-weight:800}.toolbar button{padding:4px 9px}.wrap{position:relative;min-height:0;overflow:hidden}canvas{position:absolute;inset:0;width:100%;height:100%;touch-action:none;cursor:crosshair}.hint{position:absolute;left:12px;bottom:10px;padding:6px 9px;background:#ffffffe8;font-size:11px;pointer-events:none}.makefig{width:100%;padding:10px;background:#4a3b8f;border-color:#4a3b8f;color:white}.buildlog{font:11px/1.35 ui-monospace,Menlo,monospace;white-space:pre-wrap;max-height:170px;overflow:auto;background:#eef1ee;border:1px solid var(--rule);padding:6px;margin-top:7px;display:none}@media(max-width:850px){.main{grid-template-columns:255px minmax(0,1fr)}aside{padding:9px}.sub{display:none}}
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:7px;margin-top:7px}.class.vmn.active{border-color:#be3eae;background:#faeafa}.counts3{display:grid;grid-template-columns:repeat(3,1fr);gap:5px;margin-top:8px;font-size:12px}.count.vmn{border-color:#be3eae;background:#fbedf9}.absence{margin-top:9px;padding-top:5px;border-top:1px dashed var(--rule)}.absence textarea{min-height:38px;margin:3px 0 7px}.registration{font-size:11px;padding:8px;border-left:4px solid var(--ok);background:#edf7f1;overflow-wrap:anywhere}.registration.manual{border-color:#bc711d;background:#fff3df}.registration.invalid{border-color:var(--bad);background:#fbeceb;color:var(--bad);font-weight:750}.registration-metrics{margin-top:5px;color:var(--muted)}#registrationReviewBox{margin-top:8px;font-size:12px;font-weight:750}
</style></head><body><div class="app"><header><h1>ARC / ME / VMN delineation</h1><span class="sub">HIL review on native DAPI coordinates; any nonempty region subset is accepted</span><span id="progress">Loading...</span></header><div class="main"><aside>
<div class="group"><div class="title">Review identity</div><label class="field">Reviewer name</label><input id="reviewer" type="text" maxlength="200"><label class="field">Session identifier</label><input id="session" type="text" maxlength="200" placeholder="e.g. 2026-08-14 pass 1"></div>
<div class="group"><div class="title">Section</div><select id="sectionSelect"></select><div class="grid"><button id="previous">Previous</button><button id="next">Next</button></div></div>
<div id="registrationGroup" class="group" hidden><div class="title">DAPI stitch registration</div><div id="registrationStatus" class="registration">Loading registration...</div><label id="registrationReviewBox" class="check" hidden><input id="registrationReviewed" type="checkbox"> I reviewed this registered composite and confirm that its medial seam and ventricular lumen are anatomically correct.</label></div>
<div class="group"><div class="title">Draw/edit replacement polygons</div><div class="grid3"><button id="arcClass" class="class arc active">ARC</button><button id="meClass" class="class me">ME</button><button id="vmnClass" class="class vmn">VMN</button></div><div class="grid"><button id="undo">Undo point</button><button id="close">Close polygon</button><button id="delete">Delete last ARC</button><button id="clear">Clear ARC</button></div><div class="counts3"><div class="count arc">ARC: <b id="arcCount">0</b></div><div class="count me">ME: <b id="meCount">0</b></div><div class="count vmn">VMN: <b id="vmnCount">0</b></div></div><div class="absence"><div class="title">Optional explicit absence</div><label class="check"><input id="arcAbsent" type="checkbox"> ARC absent / not identifiable</label><textarea id="arcAbsentReason" maxlength="1000" placeholder="Required only when explicitly absent" disabled></textarea><label class="check"><input id="meAbsent" type="checkbox"> ME absent / not identifiable</label><textarea id="meAbsentReason" maxlength="1000" placeholder="Required only when explicitly absent" disabled></textarea><label class="check"><input id="vmnAbsent" type="checkbox"> VMN absent / not identifiable</label><textarea id="vmnAbsentReason" maxlength="1000" placeholder="Required only when explicitly absent" disabled></textarea></div></div>
<div class="group"><div class="title">Display</div><label class="check"><input id="showAuto" type="checkbox" checked> Automated outline</label><div class="range"><span>Brightness</span><input id="brightness" type="range" min="40" max="260" value="100"><output id="brightnessOut">100%</output></div><div class="range"><span>Outline</span><input id="opacity" type="range" min="0" max="100" value="70"><output id="opacityOut">70%</output></div></div>
<div class="group"><button id="acceptAuto" class="autoaccept" disabled>Accept native automated mask unchanged</button><div class="warning">Inspect the native GFAP or POMC image itself before acceptance. Cross-assay masks are unavailable.</div><button id="accept" class="accept">Save drawn subset + accept</button><div id="status" class="status">Draw at least one ARC, ME, or VMN polygon. Other regions may remain not delineated.</div></div><div class="group"><div class="title">Figure output</div><button id="makeFigure" class="makefig" disabled>Make figure</button><label class="check"><input id="machineDefaults" type="checkbox"> Accept the machine default segmentation of every remaining section first</label><div id="buildStatus" class="status">Loading figure build configuration...</div><pre id="buildLog" class="buildlog"></pre></div></aside>
<section class="workspace"><div class="toolbar"><span id="sectionName" class="name">No section</span><span id="recordState"></span><button id="zoomOut">-</button><span id="zoomReadout">100%</span><button id="zoomIn">+</button><button id="fit">Fit</button></div><div id="wrap" class="wrap"><canvas id="canvas"></canvas><div class="hint">Click: vertex | Enter/click first: close | Space/middle-drag: pan | Wheel: zoom</div></div></section></div></div>
<script>'use strict';const $=id=>document.getElementById(id),canvas=$('canvas'),wrap=$('wrap'),ctx=canvas.getContext('2d'),regions=['ARC','ME','VMN'],colors={ARC:'#00adb5',ME:'#e67e22',VMN:'#be3eae'};const st={sections:[],index:0,currentSection:null,image:null,auto:null,polygons:{ARC:[],ME:[],VMN:[]},current:[],active:'ARC',scale:1,tx:0,ty:0,brightness:100,opacity:.7,space:false,panning:false,px:0,py:0,dirty:false,accepted:false,loading:false};
const endpoint=(path,s)=>`${path}?${new URLSearchParams({animal:s.animal,section:s.section})}`,lower=n=>n.toLowerCase();function status(msg,kind=''){$('status').textContent=msg;$('status').className=`status ${kind}`}function dirty(){st.dirty=true;st.accepted=false;record()}function record(){let e=$('recordState');e.textContent=st.dirty?'Unsaved edits':st.accepted?'HIL accepted':'Pending review';e.className=st.dirty?'dirty':st.accepted?'accepted':''}function counts(){for(let n of regions)$(lower(n)+'Count').textContent=st.polygons[n].length;$('delete').textContent=`Delete last ${st.active}`;$('clear').textContent=`Clear ${st.active}`}function refresh(){let s=$('sectionSelect');s.innerHTML='';st.sections.forEach((x,i)=>{let o=document.createElement('option');o.value=i;o.textContent=`${x.accepted?'Accepted':'Pending'} | ${x.animal} / ${x.section}`;s.append(o)});s.value=st.index;$('progress').textContent=`Accepted ${st.sections.filter(x=>x.accepted).length} / ${st.sections.length}`}function active(name){let key=lower(name);if($(key+'Absent').checked){status(`Uncheck ${name} absent before drawing ${name}.`,'error');return}st.active=name;st.current=[];for(let n of regions)$(lower(n)+'Class').classList.toggle('active',name===n);counts();draw()}function absenceUI(name,mark=false){let key=lower(name),absent=$(key+'Absent').checked;$(key+'AbsentReason').disabled=!absent;$(key+'Class').disabled=absent;if(absent&&st.active===name){let fallback=regions.find(n=>!$(lower(n)+'Absent').checked);if(fallback)active(fallback)}if(mark)dirty()}function loadImage(url){return new Promise((ok,no)=>{let i=new Image;i.onload=()=>ok(i);i.onerror=()=>no(new Error(`Could not load ${url}`));i.src=url})}function polys(r,n){return ((r.polygons&&r.polygons[n])||[]).map(p=>(p.vertices_px||p).map(v=>[+v[0],+v[1]]))}
function displayNumber(v,d=2){let n=Number(v);return Number.isFinite(n)?n.toFixed(d):'n/a'}
function registrationUI(s,rec){let group=$('registrationGroup'),box=$('registrationReviewBox'),check=$('registrationReviewed'),panel=$('registrationStatus');if(s.dataset!=='pomc'){group.hidden=true;check.checked=false;return}group.hidden=false;let valid=!!s.stitch_registration_valid,manual=!!s.stitch_registration_requires_manual_review,q=s.stitch_registration_qc||{},t=s.stitch_registration_transform||{};box.hidden=!manual;check.checked=manual&&!!(rec.pair_valid&&rec.registered_composite_reviewed&&rec.stitch_registration_sha256===s.stitch_registration_sha256);panel.className=`registration ${valid?(manual?'manual':''):'invalid'}`;if(!valid){panel.textContent=`Invalid or missing registration: ${s.stitch_registration_error||'rebuild this section before ROI review'}`;return}let details=[`method ${s.stitch_registration_method}`,`overlap ${displayNumber(t.overlap_width_px,0)} px`,`right y ${displayNumber(t.right_shift_y_px,0)} px`,`raw r ${displayNumber(q.raw_overlap_ncc,3)}`,`inliers ${displayNumber(q.sift_inlier_count,0)}`,`p95 ${displayNumber(q.sift_residual_p95_px,2)} px`].join(' · ');panel.innerHTML=`<b>${s.stitch_registration_status.replaceAll('_',' ')}</b><div class="registration-metrics">${details}<br>JSON SHA-256 ${String(s.stitch_registration_sha256).slice(0,16)}...</div>`}
function registrationFields(){let s=st.sections[st.index];if(s.dataset!=='pomc')return{};if(!s.stitch_registration_valid){status('This POMC composite has no valid stitch registration receipt. Rebuild it first.','error');return null}if(s.stitch_registration_requires_manual_review&&!$('registrationReviewed').checked){status('Review the registered composite and explicitly confirm its medial seam before saving ROIs.','error');return null}return{registered_composite_reviewed:!!s.stitch_registration_requires_manual_review,stitch_registration_confirmation:s.stitch_registration_requires_manual_review?'reviewed_registered_composite_image_by_image':null}}
function setAcceptState(){let s=st.sections[st.index],blocked=s?.dataset==='pomc'&&(!s.stitch_registration_valid||(s.stitch_registration_requires_manual_review&&!$('registrationReviewed').checked));$('accept').disabled=st.loading||!st.image||blocked;$('acceptAuto').disabled=st.loading||!st.image||!st.auto||blocked}
async function load(index,force=false){if(st.loading)return;if(!force&&st.dirty&&!confirm('Discard unsaved edits?')){refresh();return}st.loading=true;st.image=null;st.auto=null;st.index=Math.max(0,Math.min(index,st.sections.length-1));let s=st.sections[st.index];st.currentSection=s;refresh();status('Loading...');st.current=[];try{let ap=fetch(endpoint('/api/annotation',s)).then(async r=>{let j=await r.json();if(!r.ok)throw Error(j.error);return j}),ip=loadImage(endpoint('/api/image',s)),mp=s.auto_mask_available?loadImage(endpoint('/api/auto-mask',s)).catch(()=>null):Promise.resolve(null),[rec,img,mask]=await Promise.all([ap,ip,mp]);st.image=img;st.auto=mask;st.polygons=Object.fromEntries(regions.map(n=>[n,polys(rec,n)]));st.accepted=!!(rec.found&&rec.accepted&&rec.pair_valid);st.dirty=false;for(let n of regions){let key=lower(n),legacyAbsent=n==='ME'&&!!rec.me_absent,explicit=rec.region_status?.[n]==='explicitly_absent'||legacyAbsent;$(key+'Absent').checked=explicit;$(key+'AbsentReason').value=rec.region_absence_reasons?.[n]||(n==='ME'?rec.me_absent_reason:'')||'';absenceUI(n,false)}if(rec.reviewer&&!$('reviewer').value)$('reviewer').value=rec.reviewer;if(rec.session&&!$('session').value)$('session').value=rec.session;$('sectionName').textContent=`${s.dataset} | ${s.animal} / ${s.section} | ${s.width} x ${s.height} px`;$('showAuto').disabled=!mask;$('previous').disabled=st.index===0;$('next').disabled=st.index===st.sections.length-1;registrationUI(s,rec);counts();record();fit();let pending='Pending review: draw any nonempty ARC/ME/VMN subset; omitted regions remain not delineated.';if(s.dataset==='pomc'&&!s.stitch_registration_valid)pending='Registration receipt is invalid or missing; rebuild this composite before ROI review.';else if(s.dataset==='pomc'&&s.stitch_registration_requires_manual_review)pending='Pending registration review: inspect the medial seam and ventricular lumen, then check the confirmation box.';let acceptedMessage=s.replacement_allowed?'Accepted annotation loaded; segment it again and save to replace it.':'Accepted annotation loaded; replacement is locked by --lock-accepted (restart without it to segment again).';status(st.accepted?acceptedMessage:pending)}catch(e){status(e.message,'error')}finally{st.loading=false;setAcceptState()}}
function resize(){let d=devicePixelRatio||1,r=wrap.getBoundingClientRect();canvas.width=Math.max(1,Math.round(r.width*d));canvas.height=Math.max(1,Math.round(r.height*d));draw()}function fit(){if(!st.image)return;let w=canvas.clientWidth,h=canvas.clientHeight;st.scale=Math.min(w/st.image.width,h/st.image.height)*.96;st.tx=(w-st.image.width*st.scale)/2;st.ty=(h-st.image.height*st.scale)/2;draw()}function zoom(f,cx=canvas.clientWidth/2,cy=canvas.clientHeight/2){let old=st.scale;st.scale=Math.max(.02,Math.min(30,old*f));let r=st.scale/old;st.tx=cx-(cx-st.tx)*r;st.ty=cy-(cy-st.ty)*r;draw()}function poly(points,color,fill,open=false){if(!points.length)return;ctx.beginPath();ctx.moveTo(...points[0]);for(let i=1;i<points.length;i++)ctx.lineTo(...points[i]);if(!open)ctx.closePath();ctx.strokeStyle=color;ctx.lineWidth=4.5/st.scale;ctx.lineJoin='round';ctx.stroke();if(fill&&!open){ctx.fillStyle=color+'2b';ctx.fill()}for(let p of points){ctx.beginPath();ctx.arc(p[0],p[1],4/st.scale,0,7);ctx.fillStyle='#fff';ctx.fill();ctx.lineWidth=2.2/st.scale;ctx.strokeStyle=color;ctx.stroke()}}function draw(){let d=devicePixelRatio||1;ctx.setTransform(d,0,0,d,0,0);ctx.fillStyle='#101619';ctx.fillRect(0,0,canvas.clientWidth,canvas.clientHeight);if(!st.image)return;ctx.save();ctx.translate(st.tx,st.ty);ctx.scale(st.scale,st.scale);ctx.filter=`brightness(${st.brightness}%)`;ctx.drawImage(st.image,0,0);ctx.filter='none';if(st.auto&&$('showAuto').checked){ctx.globalAlpha=st.opacity;ctx.drawImage(st.auto,0,0);ctx.globalAlpha=1}for(let n of regions)for(let p of st.polygons[n])poly(p,colors[n],true);poly(st.current,colors[st.active],false,true);ctx.restore();$('zoomReadout').textContent=`${Math.round(st.scale*100)}%`}function point(e){let r=canvas.getBoundingClientRect();return[(e.clientX-r.left-st.tx)/st.scale,(e.clientY-r.top-st.ty)/st.scale]}function inside(p){return st.image&&p[0]>=0&&p[1]>=0&&p[0]<st.image.width&&p[1]<st.image.height}function close(){if(st.current.length<3){status('A polygon needs at least three vertices.','error');return}st.polygons[st.active].push(st.current.map(p=>[+p[0].toFixed(3),+p[1].toFixed(3)]));st.current=[];dirty();counts();status(`${st.active} polygon closed.`);draw()}
canvas.oncontextmenu=e=>e.preventDefault();canvas.onpointerdown=e=>{if(!st.image)return;if(e.button===1||e.button===2||st.space){st.panning=true;st.px=e.clientX;st.py=e.clientY;canvas.setPointerCapture(e.pointerId);canvas.style.cursor='grabbing';return}if(e.button)return;let p=point(e);if(!inside(p))return;if(st.current.length>=3&&Math.hypot((p[0]-st.current[0][0])*st.scale,(p[1]-st.current[0][1])*st.scale)<=12){close();return}st.current.push(p);dirty();draw()};canvas.onpointermove=e=>{if(st.panning){st.tx+=e.clientX-st.px;st.ty+=e.clientY-st.py;st.px=e.clientX;st.py=e.clientY;draw()}};function stop(e){st.panning=false;canvas.style.cursor='crosshair';try{canvas.releasePointerCapture(e.pointerId)}catch(_){}}canvas.onpointerup=stop;canvas.onpointercancel=stop;canvas.addEventListener('wheel',e=>{e.preventDefault();let r=canvas.getBoundingClientRect();zoom(e.deltaY<0?1.15:1/1.15,e.clientX-r.left,e.clientY-r.top)},{passive:false});window.onkeydown=e=>{if(['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName))return;if(e.code==='Space'){st.space=true;e.preventDefault()}if(e.key==='Enter'){close();e.preventDefault()}if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='z'){$('undo').click();e.preventDefault()}};window.onkeyup=e=>{if(e.code==='Space')st.space=false};for(let n of regions)$(lower(n)+'Class').onclick=()=>active(n);$('undo').onclick=()=>{if(st.current.length){st.current.pop();dirty();draw()}else status('No current vertex to undo.','error')};$('close').onclick=close;$('delete').onclick=()=>{let p=st.polygons[st.active];if(!p.length){status(`No closed ${st.active} polygon to delete.`,'error');return}p.pop();dirty();counts();draw()};$('clear').onclick=()=>{let n=st.active;if((st.polygons[n].length||st.current.length)&&confirm(`Clear all ${n} polygons and current vertices?`)){st.polygons[n]=[];st.current=[];dirty();counts();draw()}};$('brightness').oninput=e=>{st.brightness=+e.target.value;$('brightnessOut').textContent=`${st.brightness}%`;draw()};$('opacity').oninput=e=>{st.opacity=+e.target.value/100;$('opacityOut').textContent=`${e.target.value}%`;draw()};$('showAuto').onchange=draw;$('registrationReviewed').onchange=()=>{dirty();setAcceptState()};for(let n of regions){let key=lower(n);$(key+'Absent').onchange=()=>absenceUI(n,true);$(key+'AbsentReason').oninput=dirty}$('zoomIn').onclick=()=>zoom(1.25);$('zoomOut').onclick=()=>zoom(.8);$('fit').onclick=fit;$('previous').onclick=()=>load(st.index-1);$('next').onclick=()=>load(st.index+1);$('sectionSelect').onchange=e=>load(+e.target.value);for(let id of ['reviewer','session']){$(id).value=localStorage.getItem(`arcmevmn_${id}`)||'';$(id).oninput=e=>localStorage.setItem(`arcmevmn_${id}`,e.target.value)}
function identity(){let reviewer=$('reviewer').value.trim(),session=$('session').value.trim();if(!reviewer||!session){status('Reviewer name and session identifier are required.','error');return null}return{reviewer,session}}function absenceFields(){let region_absent={},region_absence_reasons={};for(let n of regions){let key=lower(n);region_absent[n]=$(key+'Absent').checked;region_absence_reasons[n]=$(key+'AbsentReason').value.trim()}return{region_absent,region_absence_reasons}}function validateAbsence(fields,checkPolygons=true){if(checkPolygons&&!regions.some(n=>st.polygons[n].length)){status('Draw at least one closed ARC, ME, or VMN polygon.','error');return false}for(let n of regions){if(fields.region_absent[n]&&checkPolygons&&st.polygons[n].length){status(`Clear ${n} polygons or uncheck ${n} absent.`,'error');return false}if(fields.region_absent[n]&&!fields.region_absence_reasons[n]){status(`Write a reason for explicit ${n} absence.`,'error');return false}}return true}async function confirmSaved(s,saved){let response=await fetch(endpoint('/api/annotation',s)),rec=await response.json();if(!response.ok||!rec.pair_valid||rec.mask_sha256!==saved.mask_sha256)throw Error(rec.error||'Disk reload/hash validation failed after save');return rec}
$('accept').onclick=async()=>{if(st.current.length){status('Close or undo the current polygon before acceptance.','error');return}let id=identity();if(!id)return;let absence=absenceFields();if(!validateAbsence(absence,true))return;let s=st.sections[st.index],registration=registrationFields();if(registration===null)return;st.loading=true;setAcceptState();status('Validating, atomically writing, and reloading paired files...');try{let response=await fetch('/api/annotation',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({animal:s.animal,section:s.section,...id,...absence,...registration,polygons:st.polygons})}),r=await response.json();if(!response.ok||!r.pair_valid)throw Error(r.error||'Acceptance failed durable-pair validation');await confirmSaved(s,r);st.dirty=false;st.accepted=true;s.accepted=true;refresh();record();let drawn=regions.filter(n=>r.region_status[n]==='drawn').join(', ');status(`Saved and reloaded ${drawn} subset ${r.accepted_at}. SHA-256 ${r.mask_sha256.slice(0,16)}...`,'ok')}catch(e){status(e.message,'error')}finally{st.loading=false;setAcceptState()}};
$('acceptAuto').onclick=async()=>{let id=identity();if(!id)return;let absence=absenceFields(),s=st.sections[st.index],registration=registrationFields();if(registration===null)return;if(!validateAbsence(absence,false))return;if(!st.auto||!s.auto_mask_available){status('No native automated mask is available for unchanged acceptance.','error');return}if(!$('showAuto').checked){status('Display the native automated outline before accepting it.','error');return}if(!confirm(`Confirm that you inspected ${s.animal} / ${s.section} image by image and accept the native automated mask unchanged.`))return;st.loading=true;setAcceptState();status('Revalidating, atomically writing, and reloading native DAPI/mask hashes...');try{let response=await fetch('/api/accept-automated-mask',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({animal:s.animal,section:s.section,...id,...absence,...registration,confirm_automated_mask_unchanged:true,automated_mask_confirmation:'reviewed_image_by_image'})}),r=await response.json();if(!response.ok||!r.pair_valid)throw Error(r.error||'Automated-mask acceptance failed');await confirmSaved(s,r);st.dirty=false;st.accepted=true;s.accepted=true;refresh();record();status(`Saved and reloaded unchanged ${r.accepted_source} ${r.accepted_at}. SHA-256 ${r.mask_sha256.slice(0,16)}...`,'ok')}catch(e){status(e.message,'error')}finally{st.loading=false;setAcceptState()}};
let figureConfig={figure_button_enabled:false,figure_label:'figure',machine_default_attestation:'',output_dir:'',figure_output_root:'',figure_analysis_dir:'',figure_analysis_command:'',source_stale:false,source_path:''};
function buildStatus(msg,kind=''){$('buildStatus').textContent=msg;$('buildStatus').className=`status ${kind}`}
function buildLog(text){let e=$('buildLog');e.textContent=text||'';e.style.display=text?'block':'none'}
function makeFigureEnabled(on){$('makeFigure').disabled=!(figureConfig.figure_button_enabled&&on)}
function staleNote(c){return c&&c.source_stale?`STALE SERVER: this process is running an older copy of ${c.source_path||'the reviewer'} than the one now on disk, so the button would run the superseded commands. Stop it and relaunch. `:''}
function buildPlanText(c){if(!c.figure_button_enabled)return 'The figure build button is disabled for this session.';let render=`renders ${c.figure_label} into a new directory under ${c.figure_output_root}`;return staleNote(c)+(c.figure_analysis_command?`Analyzes the accepted receipts in ${c.output_dir} into ${c.figure_analysis_dir}, then ${render}.`:`Reuses the completed analysis ${c.figure_analysis_dir} and ${render}.`)}
async function loadConfig(){try{let response=await fetch('/api/config'),c=await response.json();if(!response.ok)throw Error(c.error);figureConfig=c;$('makeFigure').textContent=`Make ${c.figure_label}`;makeFigureEnabled(true);buildStatus(buildPlanText(c),c.source_stale?'error':'')}catch(e){buildStatus(e.message,'error')}}
async function pollBuild(){let b;try{let response=await fetch('/api/figure-status');b=await response.json();if(!response.ok)throw Error(b.error)}catch(e){buildStatus(e.message,'error');return false}let head=`${staleNote(b)}${b.accepted} / ${b.total} accepted`;if(b.running){makeFigureEnabled(false);buildStatus(`${head} | ${b.stage}: ${b.command}`);return true}makeFigureEnabled(true);if(b.ok===true){let files=b.figure_files||[];buildStatus(`${head} | ${b.figure_label} written: ${files.length?files.join('   '):b.output_dir}${b.machine_defaults_accepted?`, after accepting ${b.machine_defaults_accepted} machine default sections`:''}.`,'ok');buildLog([files.length?`${b.figure_label} master files:\n`+files.join('\n')+`\n\nPanels, legends, source data and provenance for this run are in ${b.output_dir}\n`:'',b.stdout_tail,b.stderr_tail].filter(Boolean).join('\n'))}else if(b.ok===false){buildStatus(`${head} | ${b.error}`,'error');buildLog([(b.machine_defaults_skipped||[]).join('\n'),b.stderr_tail,b.stdout_tail].filter(Boolean).join('\n'))}else buildStatus(`${head} | ${figureConfig.figure_label} has not been built in this session.`);return false}
async function watchBuild(){while(await pollBuild()){await new Promise(r=>setTimeout(r,2000))}try{let response=await fetch('/api/sections'),r=await response.json();if(response.ok){st.sections=r.sections;refresh()}}catch(e){}}
$('makeFigure').onclick=async()=>{let id=identity();if(!id){buildStatus('Fill in the reviewer name and session identifier at the top of this sidebar first. Nothing was started.','error');return}let useMachine=$('machineDefaults').checked,attestation='';if(useMachine){attestation=(prompt(`Accepting the machine default segmentation of every remaining section writes one accepted receipt per section under your name, without opening each image.\n\nType exactly:\n${figureConfig.machine_default_attestation}`)||'').trim();if(attestation.toUpperCase()!==String(figureConfig.machine_default_attestation).toUpperCase()){buildStatus('Machine default acceptance was not confirmed; nothing was written.','error');return}}else if(!confirm(figureConfig.figure_analysis_command?`Analyze the currently accepted segmentations and then render ${figureConfig.figure_label}?\n\nThe analysis reads the prepared review images and writes ${figureConfig.figure_analysis_dir}. It takes a while; this page will follow it, and the render lands under ${figureConfig.figure_output_root}.`:`Render ${figureConfig.figure_label} now from the currently accepted segmentations, into a new directory under ${figureConfig.figure_output_root}?`)){buildStatus('Cancelled at the confirmation dialog; nothing was started.');return}makeFigureEnabled(false);buildStatus('Starting...');buildLog('');try{let response=await fetch('/api/make-figure',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({...id,machine_defaults:useMachine,machine_default_attestation:attestation})}),r=await response.json();if(!response.ok||!r.ok)throw Error(r.error||'Could not start the figure build')}catch(e){buildStatus(e.message,'error');makeFigureEnabled(true);return}watchBuild()};
window.onbeforeunload=e=>{if(st.dirty){e.preventDefault();e.returnValue=''}};new ResizeObserver(resize).observe(wrap);(async()=>{try{let response=await fetch('/api/sections'),r=await response.json();if(!response.ok)throw Error(r.error);st.sections=r.sections;if(!st.sections.length)throw Error('No sections found');refresh();await load(0,true);await loadConfig();await watchBuild()}catch(e){status(e.message,'error')}})();</script></body></html>'''


def parse_args(argv=None):
    p=argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Native-only browser UI for image-by-image ARC/ME/VMN review",
        epilog=("Datasets: gfap reviews native full sections; pomc reviews exact native reconstructed "
                "source-sample sections. Cross-assay manifests and registered masks are unsupported."))
    p.add_argument("--dataset",choices=("gfap","pomc"),default="gfap")
    p.add_argument("--input-root",type=Path,default=None,
                   help="GFAP sanitized root or POMC reconstructed_sections root")
    p.add_argument("--analysis-dir",type=Path,default=None,
                   help="GFAP analysis directory containing native automated masks; ignored for pomc")
    p.add_argument("--output-dir",type=Path,default=None,
                   help="Accepted native-coordinate JSON/TIFF root")
    p.add_argument("--host",default="127.0.0.1")
    p.add_argument("--port",type=int,default=8080)
    p.add_argument("--validate-only",action="store_true",help="Validate discovery plus every persisted JSON/TIFF/hash receipt without starting a server")
    p.add_argument("--force",action="store_true",
                   help="Accepted for compatibility. Replacement through Save + accept is now always allowed; no accepted review is ever deleted")
    p.add_argument("--lock-accepted","--no-force",dest="lock_accepted",action="store_true",
                   help="Restore the historical lock that refused to replace an already accepted review")
    p.add_argument("--figure-command",default=None,
                   help=("Shell-quoted command the Make Figure 4/6 button runs. "
                         "A literal {output_dir} in it is replaced by a fresh render directory. "
                         "The default renders the canonical Figure 4/6 for this dataset"))
    p.add_argument("--figure-analysis-dir",type=Path,default=None,
                   help=("Analysis directory the Make Figure 4/6 button renders from. This is "
                         "the analyzer's output, not this reviewer's --analysis-dir. A completed "
                         "one is reused; an absent one is where the button's own analysis run "
                         "writes. Omit it and the button picks a fresh directory under "
                         "--figure-output-root"))
    p.add_argument("--figure-input-root",type=Path,default=None,
                   help=("Raw/sanitized input the button's analysis run reads. Defaults to the "
                         "analyzer's own default, which is the raw data shipped with this tree"))
    p.add_argument("--figure-prepared-dir",type=Path,default=None,
                   help=("Prepared review input the button's analysis run validates the accepted "
                         "receipts against. Defaults to this reviewer's --analysis-dir, which is "
                         "the directory the review images were served from"))
    p.add_argument("--no-figure-analyze",action="store_true",
                   help=("Refuse to run the analyzer from the Make Figure 4/6 button. The button "
                         "then only renders, and only from a completed --figure-analysis-dir"))
    p.add_argument("--figure-output-root",type=Path,default=None,
                   help=("Parent of the fresh, absent directory each Make Figure run renders "
                         "into. Defaults to this figure's own analyses/FigN/figure directory "
                         "inside the tree, where the render can be found afterwards"))
    p.add_argument("--no-figure-button",action="store_true",help="Hide the Make Figure 4/6 button")
    args=p.parse_args(argv)
    if args.dataset=="gfap":
        args.input_root=args.input_root or DEFAULT_INPUT
        args.analysis_dir=args.analysis_dir or DEFAULT_ANALYSIS
        args.output_dir=args.output_dir or DEFAULT_OUTPUT
    else:
        args.input_root=args.input_root or DEFAULT_POMC_INPUT
        args.analysis_dir=args.analysis_dir or DEFAULT_POMC_ANALYSIS
        args.output_dir=args.output_dir or DEFAULT_POMC_OUTPUT
    args.figure_output_root=(
        args.figure_output_root or FIGURE_DEFAULT_OUTPUT_ROOT[args.dataset]
    )
    return args


def loopback(host):
    if host.lower()=="localhost": return True
    try: return ipaddress.ip_address(host).is_loopback
    except ValueError: return False


def main(argv=None):
    args=parse_args(argv)
    if not 0<=args.port<=65535: raise SystemExit("--port must be 0..65535")
    if not loopback(args.host):
        print("WARNING: non-loopback binding exposes microscopy and annotation writes without authentication or TLS.",file=sys.stderr)
    figure_command=shlex.split(args.figure_command) if args.figure_command else None
    app=App(args.dataset,args.input_root,args.analysis_dir,args.output_dir,
            allow_replace=not args.lock_accepted,figure_command=figure_command,
            figure_output_root=args.figure_output_root,figure_button=not args.no_figure_button,
            figure_analysis=args.figure_analysis_dir,
            figure_input_root=args.figure_input_root,
            figure_prepared_dir=args.figure_prepared_dir,
            figure_analyze=not args.no_figure_analyze)
    if args.validate_only:
        rows=app.rows()
        invalid = [
            f"{row['animal']}/{row['section']}"
            for row in rows if row["invalid_existing_pair"]
        ]
        accepted_n = sum(bool(row["accepted"]) for row in rows)
        print(json.dumps({"ok":not invalid,"dataset":args.dataset,"sections":len(rows),
                          "native_automated_mask_available":sum(bool(row["auto_mask_available"]) for row in rows),
                          "accepted":accepted_n,"pending":len(rows)-accepted_n,
                          "invalid_existing_pairs":len(invalid),
                          "invalid_existing_identities":invalid,
                          "output_dir":str(app.output)},indent=2,sort_keys=True))
        return 0
    server=Server((args.host,args.port),app)
    host,port=server.server_address[:2]; browser="127.0.0.1" if host in ("0.0.0.0","::") else host
    print(f"ARC/ME/VMN annotation server: http://{browser}:{port}"); print(f"Dataset: {args.dataset}"); print(f"Sections: {len(app.sections)}"); print(f"Accepted outputs: {app.output}")
    print(f"Re-segmenting an accepted section: {'locked by --lock-accepted' if app.allow_replace is False else 'allowed; Save + accept replaces the pair'}")
    if app.figure_button:
        plan=app.figure_build_plan()
        if plan["analysis_reused"]:
            print(f"Make {app.figure_label} reuses the completed analysis {plan['analysis_dir']}")
        elif plan["analysis_command"]:
            print(f"Make {app.figure_label} will analyze the accepted segmentations into "
                  f"{plan['analysis_dir']}, validating them against {app.figure_prepared}, then render")
            print(f"  {shlex.join(plan['analysis_command'])}")
        print(f"Make {app.figure_label} renders into a fresh directory under {app.figure_output_root}")
        print(f"  {shlex.join(plan['render_command'])}")
        stem = FIGURE_MASTER_NAME.get(args.dataset)
        if stem:
            print(f"  master graphics: {plan['render_dir']}/{stem}.pdf and .png")
    print("Press Ctrl-C to stop.")
    try: server.serve_forever(.25)
    except KeyboardInterrupt: print("\nStopping server.")
    finally: server.server_close()
    return 0


if __name__=="__main__": raise SystemExit(main())
