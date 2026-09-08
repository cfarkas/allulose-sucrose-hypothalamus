#!/usr/bin/env python3
"""
Human-guided July 2026 c-FOS and ACTH/CLIP analysis (ARC/ME endpoints).

This standalone pipeline consumes ``zstack_tiffs_flat/channel_map.csv`` and
``zstack_tiffs_flat/sample_manifest.csv`` produced by ``export_czi_zstacks.py``.
It deliberately keeps anatomical annotation DAPI-only: c-FOS and ACTH_CLIP are
never displayed to, or used by, the ARC/ME/VMN annotation step.

Stages
------
``--prepare``
    Validate manifests, refuse unresolved/registration-required acquisitions
    unless explicitly permitted, make one 2-D maximum-intensity projection per
    DAPI/cFOS/ACTH-CLIP channel, and write ``analysis_samplesheet.csv``.
``--segment``
    Run user-supplied Cellpose models on missing DAPI/cFOS/ACTH-CLIP MIPs. Existing
    ``*_seg.npy`` files are preserved unless ``--overwrite-seg`` is supplied.
``--annotate``
    Launch a resumable Flask polygon interface for HIL ARC, ME, and VMN
    annotation over DAPI segmentation or raw DAPI. Statistical endpoints remain
    ARC/ME, matching Figure 4.
``--status``
    Print preparation, segmentation, and annotation status.
``--quantify``
    Associate marker objects to unique DAPI nuclei, assign HIL regions,
    write nucleus/image/animal tables, check count invariants, render QC, run
    animal-level statistics, and make compact ME/ARC plots.
``--self-test``
    Run a small in-memory test of assignments and invariants.

Typical use
-----------
python july_cfos_acth_clip_human_in_loop.py --prepare

python july_cfos_acth_clip_human_in_loop.py --segment \
  --model-dapi /path/to/dapi_model \
  --model-cfos /path/to/cfos_model \
  --model-acth-clip /path/to/acth_clip_model

python july_cfos_acth_clip_human_in_loop.py --annotate --port 8050
python july_cfos_acth_clip_human_in_loop.py --quantify --qc-format both

Scientific assumptions
----------------------
* The fluor-to-marker mapping in channel_map.csv has been checked against the
  staining notebook. The exporter labels it provisional; this script records
  that status rather than silently claiming confirmation.
* Each acquisition is analyzed as a 2-D MIP for compatibility with the Nancy
  workflow. Z plane count/step/span are retained as QC covariates.
* A DAPI nucleus is the cell identity. A c-FOS ROI is accepted only by direct
  overlap with one best DAPI nucleus (default ROI overlap fraction >= 0.40).
* Each ACTH_CLIP ROI is assigned to at most one best expanded-DAPI territory.
* Animal is the inferential unit. Multiple acquisitions for an animal are
  summed before ratios are recomputed.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

from scipy.stats import f_oneway, kruskal, mannwhitneyu
from skimage.draw import polygon as draw_polygon
from skimage.measure import label, regionprops
from skimage.segmentation import expand_labels, find_boundaries
import tifffile

try:
    from PIL import Image
except Exception:  # only required for the annotation web server
    Image = None


SCRIPT_VERSION = "2.1.0"
ACTH_CLIP_DISPLAY = "ACTH/CLIP"
ACTH_CLIP_DAPI_AREA_QUANTILE = 0.05
EXPERIMENT_DESCRIPTION = "16 h fasting; first-ever acute solution exposure"
REGION_NAMES = ("ARC", "ME", "VMN")
REGION_CODES = {"OUTSIDE": 0, "ARC": 1, "ME": 2, "VMN": 3}
CODE_TO_REGION = {v: k for k, v in REGION_CODES.items()}
REGION_COLORS = {"ARC": "#00adb5", "ME": "#e67e22", "VMN": "#be3eae"}
REGION_RGB = {"ARC": (0, 173, 181), "ME": (230, 126, 34), "VMN": (190, 62, 174)}
COND_COLORS = {
    "Water": "#b9e3f2",
    "Sucrose": "#e31a1c",
    "Allulose": "#2ecc71",
    "Control": "#999999",
}
COND_ALIASES = {
    "water": "Water", "agua": "Water", "h2o": "Water",
    "sucrose": "Sucrose", "sacarosa": "Sucrose",
    "allulose": "Allulose", "alulosa": "Allulose", "psicose": "Allulose",
    "control": "Control",
}
BIOLOGICAL_CONDITIONS = ["Water", "Sucrose", "Allulose"]

COUNT_COLUMNS = [
    "n_dapi", "n_cfos_only", "n_acth_clip_only", "n_double", "n_neither",
    "n_total_cfos", "n_total_acth_clip",
]
ENDPOINTS: Dict[str, Tuple[str, str]] = {
    "double_over_total_acth_clip": ("c-FOS+ ACTH/CLIP+ / total ACTH/CLIP", "Activated ACTH/CLIP fraction"),
    "total_cfos_over_dapi": ("total c-FOS+ / DAPI", "Total c-FOS activation"),
    "total_acth_clip_over_dapi": ("total ACTH/CLIP+ / DAPI", "ACTH/CLIP fraction of DAPI"),
    "double_over_dapi": ("c-FOS+ ACTH/CLIP+ / DAPI", "Double-positive fraction of DAPI"),
    "cfos_only_over_dapi": ("c-FOS-only / DAPI", "c-FOS-only fraction"),
    "acth_clip_only_over_dapi": ("ACTH/CLIP-only / DAPI", "ACTH/CLIP-only fraction"),
    "double_over_total_cfos": ("c-FOS+ ACTH/CLIP+ / total c-FOS", "ACTH/CLIP fraction of c-FOS"),
    "double_per_100k_um2": ("c-FOS+ ACTH/CLIP+ / 100,000 µm²", "Double-positive density"),
    "total_cfos_per_100k_um2": ("total c-FOS+ / 100,000 µm²", "Total c-FOS density"),
    "total_acth_clip_per_100k_um2": ("total ACTH/CLIP+ / 100,000 µm²", "Total ACTH/CLIP density"),
    "n_double": ("c-FOS+ ACTH/CLIP+ nuclei", "Double-positive raw counts"),
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_name(value: Any) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value).strip())
    out = re.sub(r"_+", "_", out).strip("_")
    return out or "unnamed"


def token(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


def bool_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y", "t"}


def canonical_condition(value: Any) -> str:
    raw = "" if value is None else str(value).strip()
    return COND_ALIASES.get(token(raw), raw or "REVIEW")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def atomic_tiff_write(path: Path, arr: np.ndarray, **kwargs: Any) -> None:
    ensure_dir(path.parent)
    tmp = path.with_name(path.name + ".tmp.tif")
    tifffile.imwrite(str(tmp), np.asarray(arr), **kwargs)
    os.replace(tmp, path)


def read_stack(path: Path) -> np.ndarray:
    arr = np.asarray(tifffile.imread(str(path)))
    # OME exports are Z,Y,X. Remove singleton axes, but retain RGB only when
    # it is genuinely a final color axis (the channel-separated exporter is gray).
    arr = np.squeeze(arr)
    return arr


def maximum_projection(arr: np.ndarray) -> np.ndarray:
    x = np.asarray(arr)
    if x.ndim == 2:
        return x
    if x.ndim == 3 and x.shape[-1] in (3, 4):
        return x[..., :3].max(axis=-1)
    while x.ndim > 2:
        x = np.max(x, axis=0)
    if x.ndim != 2:
        raise ValueError(f"Could not reduce image to 2-D; final shape={x.shape}")
    return x


def normalize_u8(arr: np.ndarray, low: float = 0.2, high: float = 99.8) -> np.ndarray:
    x = np.asarray(arr, dtype=np.float32)
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return np.zeros(x.shape, dtype=np.uint8)
    lo, hi = np.percentile(finite, [low, high])
    if not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    return np.clip((x - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)


def load_segmentation(path: Path) -> np.ndarray:
    obj = np.load(str(path), allow_pickle=True)
    if isinstance(obj, np.ndarray) and obj.dtype == object and obj.shape == ():
        obj = obj.item()
    if isinstance(obj, Mapping) and "masks" in obj:
        obj = obj["masks"]
    elif isinstance(obj, (tuple, list)) and obj:
        obj = obj[0]
    arr = np.asarray(obj)
    arr = np.squeeze(arr)
    while arr.ndim > 2:
        arr = arr[0]
    if arr.ndim != 2:
        raise ValueError(f"Segmentation {path} is not 2-D: {arr.shape}")
    if np.issubdtype(arr.dtype, np.integer) and int(np.max(arr)) > 1:
        return arr.astype(np.int32, copy=False)
    return label(arr > 0, connectivity=1).astype(np.int32)


def save_segmentation(path: Path, masks: np.ndarray, **metadata: Any) -> None:
    ensure_dir(path.parent)
    obj = {"masks": np.asarray(masks, dtype=np.int32), "pipeline": SCRIPT_VERSION, **metadata}
    tmp = path.with_name(path.name + ".tmp.npy")
    np.save(str(tmp), obj, allow_pickle=True)
    os.replace(tmp, path)


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def source_segmentation_path(source_tiff: Path) -> Path:
    name = source_tiff.name
    if not name.lower().endswith(".tif"):
        raise ValueError(f"Expected .tif source for Cellpose sidecar lookup: {source_tiff}")
    return source_tiff.with_name(name[:-4] + "_seg.npy")


def segmentation_provenance(path: Path) -> Dict[str, Any]:
    raw = np.load(str(path), allow_pickle=True)
    obj = raw.item() if isinstance(raw, np.ndarray) and raw.dtype == object and raw.shape == () else {}
    masks = load_segmentation(path)
    return {
        "source_segmentation": str(path),
        "source_segmentation_sha256": sha256_file(path),
        "source_segmentation_shape": f"{masks.shape[0]}x{masks.shape[1]}",
        "source_segmentation_objects": int(masks.max()),
        "source_segmentation_model_path": str(obj.get("model_path", "")) if isinstance(obj, Mapping) else "",
        "source_segmentation_flow_threshold": obj.get("flow_threshold", "") if isinstance(obj, Mapping) else "",
        "source_segmentation_cellprob_threshold": obj.get("cellprob_threshold", "") if isinstance(obj, Mapping) else "",
        "source_segmentation_manual_changes": len(obj.get("manual_changes", [])) if isinstance(obj, Mapping) else 0,
    }


def stage_source_segmentation(source_tiff: Path, target: Path, expected_shape: Tuple[int, int]) -> Dict[str, Any]:
    source = source_segmentation_path(source_tiff)
    if not source.is_file():
        raise FileNotFoundError(f"Completed Cellpose sidecar missing: {source}")
    masks = load_segmentation(source)
    if tuple(map(int, masks.shape)) != tuple(map(int, expected_shape)):
        raise ValueError(
            f"Cellpose/image shape mismatch for {source}: {masks.shape} vs {expected_shape}"
        )
    source_hash = sha256_file(source)
    if target.exists():
        target_hash = sha256_file(target)
        if target_hash != source_hash:
            raise FileExistsError(
                f"Prepared segmentation collision: {target} differs from source {source}"
            )
    else:
        ensure_dir(target.parent)
        incoming = target.with_name(target.name + ".incoming")
        if incoming.exists():
            raise FileExistsError(f"Stale incoming segmentation blocks safe staging: {incoming}")
        shutil.copy2(source, incoming)
        if sha256_file(incoming) != source_hash:
            raise RuntimeError(f"Segmentation copy verification failed: {source} -> {incoming}")
        os.replace(incoming, target)
    provenance = segmentation_provenance(source)
    provenance["prepared_segmentation"] = str(target)
    provenance["prepared_segmentation_sha256"] = sha256_file(target)
    return provenance


def paths(args: argparse.Namespace) -> Dict[str, Path]:
    export_dir = Path(args.export_dir).expanduser().resolve()
    workdir = Path(args.workdir).expanduser().resolve()
    return {
        "export": export_dir,
        "channel_map": Path(args.channel_map).expanduser().resolve() if args.channel_map else export_dir / "channel_map.csv",
        "manifest": Path(args.sample_manifest).expanduser().resolve() if args.sample_manifest else export_dir / "sample_manifest.csv",
        "work": workdir,
        "prepared": workdir / "prepared",
        "sheet": workdir / "analysis_samplesheet.csv",
        "prepare_audit": workdir / "prepare_audit.csv",
        "annotations": workdir / "annotations",
        "region_masks": workdir / "region_masks",
        "results": workdir / "results",
    }


def portable_work_path(workdir: Path, path: Path) -> str:
    resolved = path.expanduser().resolve()
    try:
        return resolved.relative_to(workdir.resolve()).as_posix()
    except ValueError as exc:
        raise ValueError(f"Generated path is outside the work tree: {path}") from exc


def resolve_work_path(workdir: Path, value: Any) -> Path:
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = workdir / path
    return path.resolve()


def parse_requested_acquisitions(value: str) -> Set[str]:
    return {x.strip() for x in str(value or "").split(",") if x.strip()}


def _find_channel_rows(channel_map: pd.DataFrame, acquisition_id: str, marker: str) -> pd.DataFrame:
    sub = channel_map[
        (channel_map["acquisition_id"].astype(str) == str(acquisition_id))
        & (channel_map["marker"].astype(str).map(token) == token(marker))
    ].copy()
    if "include_channel" in sub.columns:
        sub = sub[sub["include_channel"].map(bool_value)].copy()
    return sub


def _json_scalar(value: Any, default: float = float("nan")) -> float:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return default
    try:
        decoded = json.loads(str(value))
        if isinstance(decoded, list):
            vals = [float(x) for x in decoded if x is not None]
            return float(vals[0]) if vals else default
        return float(decoded)
    except Exception:
        try:
            return float(value)
        except Exception:
            return default


def center_crop(arr: np.ndarray, shape: Tuple[int, int]) -> np.ndarray:
    h, w = arr.shape; th, tw = shape
    if th > h or tw > w:
        raise ValueError(f"Cannot center crop {arr.shape} to {shape}")
    y0=(h-th)//2; x0=(w-tw)//2
    return np.asarray(arr[y0:y0+th, x0:x0+tw])


def dapi_cv_vector(path: Path, size: int = 96) -> np.ndarray:
    """Small contrast-normalized DAPI thumbnail used only for duplicate screening."""
    arr=maximum_projection(read_stack(path)).astype(np.float32)
    ys=np.linspace(0,arr.shape[0]-1,size).round().astype(int); xs=np.linspace(0,arr.shape[1]-1,size).round().astype(int)
    thumb=arr[np.ix_(ys,xs)]; lo,hi=np.percentile(thumb,[1,99]); thumb=np.clip((thumb-lo)/max(hi-lo,1e-6),0,1)
    vec=thumb.ravel(); vec-=float(vec.mean()); norm=float(np.linalg.norm(vec))
    return vec/norm if norm>0 else vec


def registration_overlay(mips: Mapping[str, Path], output: Path) -> None:
    if Image is None:
        return
    d=normalize_u8(read_stack(mips["DAPI"])); c=normalize_u8(read_stack(mips["cFOS"])); p=normalize_u8(read_stack(mips["ACTH_CLIP"])); rgb=np.stack([p,c,d],axis=-1)
    ensure_dir(output.parent); Image.fromarray(rgb).save(output)


def prepare(args: argparse.Namespace) -> int:
    p = paths(args)
    ensure_dir(p["work"])
    ensure_dir(p["prepared"])
    if not p["channel_map"].exists() or not p["manifest"].exists():
        raise FileNotFoundError(f"Missing manifest(s): {p['channel_map']} or {p['manifest']}")
    channel_map = pd.read_csv(p["channel_map"])
    manifest = pd.read_csv(p["manifest"])
    required_map = {"acquisition_id", "marker", "output_tiff"}
    required_manifest = {"acquisition_id", "sample_id", "animal_id", "condition"}
    if required_map - set(channel_map.columns):
        raise ValueError(f"channel_map.csv missing: {sorted(required_map - set(channel_map.columns))}")
    if required_manifest - set(manifest.columns):
        raise ValueError(f"sample_manifest.csv missing: {sorted(required_manifest - set(manifest.columns))}")

    requested = parse_requested_acquisitions(args.include_acquisitions)
    audit: List[Dict[str, Any]] = []
    sheet_rows: List[Dict[str, Any]] = []
    for _, mrow in manifest.iterrows():
        acq = str(mrow["acquisition_id"])
        if requested and acq not in requested:
            continue
        # Reused animal IDs are allowed: acquisition_id is the analysis unit.
        # Therefore an explicit filename condition wins even when the exporter
        # marked an animal-level cross-acquisition conflict for review.
        cond_file = canonical_condition(mrow.get("condition_filename", "REVIEW"))
        cond_manifest = canonical_condition(mrow.get("condition", "REVIEW"))
        if cond_file in {"Water", "Sucrose", "Allulose", "Control"}:
            cond = cond_file
            condition_resolution = "EXPLICIT_FILENAME_ACQUISITION_UNIT"
        elif str(mrow.get("condition_source", "")).startswith("historical") and cond_manifest in {"Water", "Sucrose", "Allulose", "Control"}:
            cond = cond_manifest
            condition_resolution = "HISTORICAL_EXACT_MATCH_REVIEW"
        else:
            cond = cond_manifest
            condition_resolution = "UNRESOLVED"
        include_main = bool_value(mrow.get("include_main_analysis", False))
        condition_status = str(mrow.get("condition_status", ""))
        unresolved = cond not in {"Water", "Sucrose", "Allulose", "Control"}
        registration_required = bool_value(mrow.get("registration_required", False))
        reasons: List[str] = []
        if not include_main and not args.include_nonmain_cohorts:
            reasons.append("NON_MAIN_COHORT_REFUSED")
        if unresolved and not args.include_unresolved:
            reasons.append("UNRESOLVED_METADATA_REFUSED")
        if registration_required and not args.include_registration_required:
            reasons.append("REGISTRATION_REQUIRED_REFUSED")

        selected: Dict[str, pd.Series] = {}
        for marker in ["DAPI", "cFOS", "ACTH_CLIP"]:
            rows = _find_channel_rows(channel_map, acq, marker)
            if len(rows) != 1:
                reasons.append(f"{marker}_CHANNEL_ROWS_{len(rows)}")
            elif not str(rows.iloc[0].get("output_tiff", "")).strip():
                reasons.append(f"{marker}_OUTPUT_MISSING")
            else:
                selected[marker] = rows.iloc[0]

        channel_statuses = sorted({
            str(row.get("channel_map_status", "")).strip()
            for row in selected.values()
        })
        channel_map_unconfirmed = bool(channel_statuses) and any(
            not status.upper().startswith("CONFIRMED")
            for status in channel_statuses
        )
        if channel_map_unconfirmed and not args.include_unresolved:
            reasons.append(
                "CHANNEL_MAP_UNCONFIRMED:" + ",".join(channel_statuses)
            )

        if reasons:
            audit.append({
                "acquisition_id": acq, "sample_id": mrow.get("sample_id", ""),
                "animal_id": mrow.get("animal_id", ""), "condition": cond,
                "decision": "REFUSED", "reason": ";".join(reasons),
                "unresolved": unresolved, "registration_required": registration_required,
            })
            continue

        acq_dir = ensure_dir(p["prepared"] / safe_name(acq))
        mip_paths: Dict[str, Path] = {}
        source_paths: Dict[str, Path] = {}
        source_shapes: Dict[str, Tuple[int, int]] = {}
        for marker, crow in selected.items():
            source = Path(str(crow["output_tiff"]))
            if not source.is_absolute():
                source = p["export"] / source
            if not source.exists():
                reasons.append(f"{marker}_SOURCE_NOT_FOUND:{source}")
                continue
            source_paths[marker] = source
            out = acq_dir / f"{safe_name(acq)}_{marker}_mip.tif"
            if out.exists() and not args.overwrite_prepare:
                mip = maximum_projection(read_stack(out))
            else:
                mip = maximum_projection(read_stack(source))
                atomic_tiff_write(out, mip, photometric="minisblack", metadata={
                    "axes": "YX", "source_zstack": str(source), "marker": marker,
                    "pipeline": SCRIPT_VERSION,
                })
            mip_paths[marker] = out
            source_shapes[marker] = tuple(map(int, mip.shape))

        if reasons or set(mip_paths) != {"DAPI", "cFOS", "ACTH_CLIP"}:
            audit.append({
                "acquisition_id": acq, "sample_id": mrow.get("sample_id", ""),
                "animal_id": mrow.get("animal_id", ""), "condition": cond,
                "decision": "REFUSED", "reason": ";".join(reasons or ["MIP_PREPARATION_INCOMPLETE"]),
                "unresolved": unresolved, "registration_required": registration_required,
            })
            continue
        registration_qc = ""
        registration_status = str(mrow.get("registration_status", ""))
        if len(set(source_shapes.values())) != 1:
            if registration_required and args.include_registration_required and args.align_center_crop_registration:
                common=(min(x[0] for x in source_shapes.values()),min(x[1] for x in source_shapes.values()))
                for marker,out in mip_paths.items(): atomic_tiff_write(out,center_crop(maximum_projection(read_stack(out)),common),photometric="minisblack",metadata={"axes":"YX","registration":"center_crop_same_stage_qc_pending"})
                source_shapes={marker:common for marker in source_shapes}; registration_status="CENTER_CROP_SAME_STAGE_QC_PENDING"; registration_qc=str(acq_dir/f"{safe_name(acq)}_registration_overlay.png"); registration_overlay(mip_paths,Path(registration_qc))
            else:
                audit.append({"acquisition_id":acq,"sample_id":mrow.get("sample_id",""),"animal_id":mrow.get("animal_id",""),"condition":cond,"decision":"REFUSED","reason":f"MIP_SHAPE_MISMATCH_REGISTRATION_NOT_APPROVED:{source_shapes}","unresolved":unresolved,"registration_required":registration_required})
                continue

        seg_paths = {
            marker: acq_dir / f"{safe_name(acq)}_{marker}_seg.npy"
            for marker in ("DAPI", "cFOS", "ACTH_CLIP")
        }
        seg_provenance: Dict[str, Dict[str, Any]] = {}
        for marker in ("DAPI", "cFOS", "ACTH_CLIP"):
            try:
                seg_provenance[marker] = stage_source_segmentation(
                    source_paths[marker], seg_paths[marker], source_shapes[marker]
                )
            except Exception as exc:
                reasons.append(f"{marker}_SEGMENTATION_IMPORT_FAILED:{exc}")
        if reasons:
            audit.append({
                "acquisition_id": acq,
                "sample_id": mrow.get("sample_id", ""),
                "animal_id": mrow.get("animal_id", ""),
                "condition": cond,
                "decision": "REFUSED",
                "reason": ";".join(reasons),
                "unresolved": unresolved,
                "registration_required": registration_required,
            })
            continue

        px_x = float(pd.to_numeric(pd.Series([mrow.get("pixel_size_x_um", np.nan)]), errors="coerce").iloc[0])
        px_y = float(pd.to_numeric(pd.Series([mrow.get("pixel_size_y_um", np.nan)]), errors="coerce").iloc[0])
        z_planes = _json_scalar(mrow.get("z_planes", np.nan))
        z_step = _json_scalar(mrow.get("z_steps_um", np.nan))
        z_span = (z_planes - 1) * z_step if np.isfinite(z_planes) and np.isfinite(z_step) and z_planes > 0 else np.nan
        unsafe = unresolved or registration_required or channel_map_unconfirmed
        sheet_rows.append({
            "analysis_unit_id": acq,
            "acquisition_id": acq,
            "sample_id": str(mrow.get("sample_id", acq)),
            "animal_id": str(mrow.get("animal_id", mrow.get("sample_id", acq))),
            "genotype": str(mrow.get("genotype", "")),
            "fasting_hours": mrow.get("fasting_hours", ""),
            "prior_sugar_exposure": str(mrow.get("prior_sugar_exposure", "")),
            "exposure_history": "first-ever solution exposure",
            "experiment_description": EXPERIMENT_DESCRIPTION,
            "antibody_display_name": ACTH_CLIP_DISPLAY,
            "condition": cond,
            "condition_status": condition_status,
            "condition_resolution": condition_resolution,
            "include_main_analysis_manifest": include_main,
            "unsafe_explicit_inclusion": unsafe,
            "registration_required": registration_required,
            "registration_status": registration_status,
            "registration_qc_overlay": registration_qc,
            "channel_map_status": ";".join(channel_statuses),
            "pixel_size_x_um": px_x, "pixel_size_y_um": px_y,
            "z_planes": z_planes, "z_step_um": z_step, "z_span_um": z_span,
            "dapi_mip": portable_work_path(p["work"], mip_paths["DAPI"]),
            "cfos_mip": portable_work_path(p["work"], mip_paths["cFOS"]),
            "acth_clip_mip": portable_work_path(p["work"], mip_paths["ACTH_CLIP"]),
            "dapi_seg": portable_work_path(p["work"], seg_paths["DAPI"]),
            "cfos_seg": portable_work_path(p["work"], seg_paths["cFOS"]),
            "acth_clip_seg": portable_work_path(p["work"], seg_paths["ACTH_CLIP"]),
            "dapi_seg_source": seg_provenance["DAPI"]["source_segmentation"],
            "dapi_seg_sha256": seg_provenance["DAPI"]["source_segmentation_sha256"],
            "dapi_seg_model_path": seg_provenance["DAPI"]["source_segmentation_model_path"],
            "cfos_seg_source": seg_provenance["cFOS"]["source_segmentation"],
            "cfos_seg_sha256": seg_provenance["cFOS"]["source_segmentation_sha256"],
            "cfos_seg_model_path": seg_provenance["cFOS"]["source_segmentation_model_path"],
            "acth_clip_seg_source": seg_provenance["ACTH_CLIP"]["source_segmentation"],
            "acth_clip_seg_sha256": seg_provenance["ACTH_CLIP"]["source_segmentation_sha256"],
            "acth_clip_seg_model_path": seg_provenance["ACTH_CLIP"]["source_segmentation_model_path"],
            "acth_clip_seg_manual_changes": seg_provenance["ACTH_CLIP"]["source_segmentation_manual_changes"],
            "source_relpaths": str(mrow.get("source_relpaths", "")),
            "prepared_at_utc": utc_now(),
        })
        audit.append({
            "acquisition_id": acq, "sample_id": mrow.get("sample_id", ""),
            "animal_id": mrow.get("animal_id", ""), "condition": cond,
            "decision": "INCLUDED_EXPLICIT_UNSAFE" if unsafe else "INCLUDED",
            "reason": "EXPLICIT_OVERRIDE" if unsafe else "",
            "unresolved": unresolved, "registration_required": registration_required,
        })

    sheet = pd.DataFrame(sheet_rows)
    pd.DataFrame(audit).to_csv(p["prepare_audit"], index=False)
    if sheet.empty:
        raise RuntimeError(
            "No safe acquisitions were prepared. See prepare_audit.csv. "
            "Confirm channel mappings first; explicit override flags are for reviewed QC exceptions only."
        )
    # Gate exact/near-exact DAPI duplicates at the acquisition level. This
    # screen does not equate reused animal IDs; only image similarity matters.
    vectors: Dict[str,np.ndarray]={}; kept=[]; duplicate_rows=[]
    for rec in sheet.to_dict("records"):
        vec=dapi_cv_vector(resolve_work_path(p["work"], rec["dapi_mip"])); best_id=""; best_corr=-1.0
        for other_id,other in vectors.items():
            corr=float(np.dot(vec,other));
            if corr>best_corr: best_corr=corr; best_id=other_id
        if best_corr>=float(args.duplicate_correlation_threshold):
            rec["duplicate_check_status"]="REVIEW_POTENTIAL_DAPI_DUPLICATE"; rec["duplicate_match_analysis_unit_id"]=best_id; rec["duplicate_correlation"]=best_corr
            if not args.include_duplicate_candidates:
                duplicate_rows.append({"acquisition_id":rec["acquisition_id"],"sample_id":rec["sample_id"],"animal_id":rec["animal_id"],"condition":rec["condition"],"decision":"REFUSED","reason":f"DAPI_CV_DUPLICATE_CANDIDATE:{best_id}:{best_corr:.6f}","unresolved":False,"registration_required":rec["registration_required"]})
                continue
        else:
            rec["duplicate_check_status"]="DISTINCT_DAPI_CV" if vectors else "FIRST_DAPI_REFERENCE"; rec["duplicate_match_analysis_unit_id"]=best_id; rec["duplicate_correlation"]=best_corr if vectors else np.nan
        vectors[rec["analysis_unit_id"]]=vec; kept.append(rec)
    audit.extend(duplicate_rows); sheet=pd.DataFrame(kept)
    pd.DataFrame(audit).to_csv(p["prepare_audit"], index=False)
    if sheet.empty: raise RuntimeError("All prepared acquisitions failed the DAPI CV duplicate gate")
    sheet.sort_values(["condition", "analysis_unit_id"]).to_csv(p["sheet"], index=False)
    print(f"Prepared {len(sheet)} acquisition(s): {p['sheet']}")
    print(f"Audit: {p['prepare_audit']}")
    return 0


def read_sheet(args: argparse.Namespace) -> pd.DataFrame:
    p = paths(args)
    if not p["sheet"].exists():
        raise FileNotFoundError(f"Run --prepare first; missing {p['sheet']}")
    df = pd.read_csv(p["sheet"])
    required = {"acquisition_id", "animal_id", "condition", "dapi_mip", "cfos_mip", "acth_clip_mip", "dapi_seg", "cfos_seg", "acth_clip_seg"}
    if required - set(df.columns):
        raise ValueError(f"analysis_samplesheet.csv missing: {sorted(required - set(df.columns))}")
    return df


def _cellpose_eval(model: Any, image: np.ndarray, args: argparse.Namespace, diameter: Optional[float]) -> np.ndarray:
    kwargs: Dict[str, Any] = {}
    if diameter is not None:
        kwargs["diameter"] = diameter
    if args.flow_threshold is not None:
        kwargs["flow_threshold"] = args.flow_threshold
    if args.cellprob_threshold is not None:
        kwargs["cellprob_threshold"] = args.cellprob_threshold
    if args.min_size_seg is not None:
        kwargs["min_size"] = args.min_size_seg
    try:
        result = model.eval(image, channels=[0, 0], **kwargs)
    except TypeError:
        result = model.eval(image, **kwargs)
    masks = result[0] if isinstance(result, tuple) else result
    return np.asarray(masks, dtype=np.int32)


def segment(args: argparse.Namespace) -> int:
    sheet = read_sheet(args)
    workdir = paths(args)["work"]
    models_given = {"DAPI": args.model_dapi, "cFOS": args.model_cfos, "ACTH_CLIP": args.model_acth_clip}
    missing_models = [m for m, value in models_given.items() if not value]
    if missing_models:
        raise SystemExit(f"--segment requires model paths for: {', '.join(missing_models)}")
    try:
        from cellpose import models as cp_models
    except Exception as exc:
        raise RuntimeError("Cellpose is unavailable. Run in a Cellpose environment.") from exc
    cp: Dict[str, Any] = {}
    for marker, model_path in models_given.items():
        if not Path(str(model_path)).expanduser().exists():
            raise FileNotFoundError(f"{marker} model not found: {model_path}")
        try:
            cp[marker] = cp_models.CellposeModel(gpu=bool(args.gpu), pretrained_model=str(model_path))
        except TypeError:
            cp[marker] = cp_models.CellposeModel(pretrained_model=str(model_path), gpu=bool(args.gpu))
    diameters = {"DAPI": args.diameter_dapi, "cFOS": args.diameter_cfos, "ACTH_CLIP": args.diameter_acth_clip}
    report: List[Dict[str, Any]] = []
    for _, row in sheet.iterrows():
        acq = str(row["acquisition_id"])
        for marker, mip_col, seg_col in [
            ("DAPI", "dapi_mip", "dapi_seg"), ("cFOS", "cfos_mip", "cfos_seg"), ("ACTH_CLIP", "acth_clip_mip", "acth_clip_seg")
        ]:
            mip = resolve_work_path(workdir, row[mip_col])
            seg_path = resolve_work_path(workdir, row[seg_col])
            if seg_path.exists() and not args.overwrite_seg:
                report.append({"acquisition_id": acq, "marker": marker, "status": "EXISTING_MASK_PRESERVED", "segmentation": str(seg_path)})
                continue
            image = maximum_projection(read_stack(mip))
            masks = _cellpose_eval(cp[marker], image, args, diameters[marker])
            if masks.shape != image.shape:
                raise RuntimeError(f"Cellpose shape mismatch for {acq} {marker}: {masks.shape} vs {image.shape}")
            save_segmentation(seg_path, masks, marker=marker, source_mip=str(mip), model=str(models_given[marker]), segmented_at_utc=utc_now())
            report.append({"acquisition_id": acq, "marker": marker, "status": "SEGMENTED", "n_objects": int(masks.max()), "segmentation": str(seg_path)})
            print(f"[SEGMENTED] {acq} {marker}: {int(masks.max())} objects")
    pd.DataFrame(report).to_csv(paths(args)["work"] / "segmentation_report.csv", index=False)
    return 0


def annotation_json_path(args: argparse.Namespace, acquisition_id: str) -> Path:
    return paths(args)["annotations"] / f"{safe_name(acquisition_id)}.json"


def region_mask_path(args: argparse.Namespace, acquisition_id: str) -> Path:
    return paths(args)["region_masks"] / f"{safe_name(acquisition_id)}_regions.tif"


def read_annotation(args: argparse.Namespace, acquisition_id: str) -> Optional[Dict[str, Any]]:
    path = annotation_json_path(args, acquisition_id)
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def rasterize_polygons(polygons: Sequence[Mapping[str, Any]], shape: Tuple[int, int]) -> np.ndarray:
    mask = np.zeros(shape, dtype=np.uint8)
    for item in polygons:
        region = str(item.get("region", "")).upper()
        if region not in REGION_NAMES:
            continue
        points = item.get("points", [])
        if len(points) < 3:
            continue
        xs = np.asarray([float(p[0]) for p in points], dtype=float)
        ys = np.asarray([float(p[1]) for p in points], dtype=float)
        rr, cc = draw_polygon(ys, xs, shape=shape)
        mask[rr, cc] = REGION_CODES[region]
    return mask


def save_annotation(args: argparse.Namespace, row: Mapping[str, Any], decision: str, polygons: Sequence[Mapping[str, Any]], notes: str, shape: Tuple[int, int]) -> None:
    acq = str(row["acquisition_id"])
    ensure_dir(paths(args)["annotations"]); ensure_dir(paths(args)["region_masks"])
    clean_polys: List[Dict[str, Any]] = []
    if decision == "include":
        for item in polygons:
            region = str(item.get("region", "")).upper()
            pts = [[float(x), float(y)] for x, y in item.get("points", [])]
            if region in REGION_NAMES and len(pts) >= 3:
                clean_polys.append({"region": region, "points": pts})
    payload = {
        "schema_version": 2, "pipeline_version": SCRIPT_VERSION,
        "acquisition_id": acq, "sample_id": str(row.get("sample_id", acq)),
        "animal_id": str(row.get("animal_id", "")), "condition": str(row.get("condition", "")),
        "decision": decision, "notes": notes, "image_height": int(shape[0]), "image_width": int(shape[1]),
        "polygons": clean_polys, "arc_present": any(x["region"] == "ARC" for x in clean_polys),
        "me_present": any(x["region"] == "ME" for x in clean_polys),
        "vmn_present": any(x["region"] == "VMN" for x in clean_polys),
        "region_codes": REGION_CODES,
        "anatomy_channels_used": ["DAPI"], "saved_at_utc": utc_now(),
    }
    annotation_json_path(args, acq).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    mask = rasterize_polygons(clean_polys, shape) if decision == "include" else np.zeros(shape, dtype=np.uint8)
    atomic_tiff_write(region_mask_path(args, acq), mask, photometric="minisblack", compression="zlib", metadata={"axes": "YX"})


def annotation_states(args: argparse.Namespace, sheet: pd.DataFrame) -> List[str]:
    states: List[str] = []
    for _, row in sheet.iterrows():
        ann = read_annotation(args, str(row["acquisition_id"]))
        if ann is None:
            states.append("pending")
        elif ann.get("decision") == "exclude":
            states.append("excluded")
        elif ann.get("decision") == "include":
            states.append("annotated")
        else:
            states.append("pending")
    return states

# -----------------------------------------------------------------------------
# DAPI-only HIL ARC/ME annotation
# -----------------------------------------------------------------------------

ANNOTATION_HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Figure S5 ARC/ME/VMN review</title><style>
:root{--ink:#18272d;--paper:#f5f6f3;--rule:#c8d0ce;--arc:#00adb5;--me:#e67e22;--vmn:#be3eae;--ok:#238636;--bad:#b3261e;--build:#5d42a8}
*{box-sizing:border-box}html,body{height:100%;margin:0}body{font:14px/1.35 Arial,sans-serif;background:#222;color:#fff;overflow:hidden}button,select,textarea{font:inherit}button,select{padding:7px;font-weight:700}button{cursor:pointer}button:disabled{opacity:.48;cursor:not-allowed}
#top{height:58px;background:#111;display:flex;align-items:center;gap:8px;padding:5px 10px;border-bottom:2px solid #333}#progress{font-weight:800;white-space:nowrap}#meta{flex:1;min-width:100px;font-weight:700;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.top-build{background:var(--build);border:2px solid #b9abef;color:#fff;padding:9px 13px;white-space:nowrap}.top-build:disabled,.build-side:disabled{opacity:.78}
#main{display:grid;grid-template-columns:minmax(0,1fr) 336px;height:calc(100vh - 58px)}#wrap{background:#000;overflow:hidden}svg{width:100%;height:100%;touch-action:none;cursor:crosshair}#side{background:var(--paper);color:var(--ink);padding:12px;overflow:auto}.section{padding-bottom:11px;margin-bottom:11px;border-bottom:1px solid var(--rule)}.small{font-size:12px;color:#47565c}.regions{display:grid;grid-template-columns:repeat(3,1fr);gap:6px}.region{padding:10px 4px;border:4px solid transparent}.region.active{border-color:#111;box-shadow:0 0 0 2px #fff inset}#arc{background:var(--arc);color:#fff}#me{background:var(--me);color:#111}#vmn{background:var(--vmn);color:#fff}.action{display:block;width:100%;margin:5px 0}.save{background:var(--ok);color:#fff}.exclude{background:var(--bad);color:#fff}.build-side{background:var(--build);color:#fff;border:2px solid #493379;padding:11px}textarea{width:100%;min-height:55px;box-sizing:border-box}.key{font-weight:800;color:#111}#buildStatus{overflow-wrap:anywhere}
</style></head><body>
<div id="top"><button onclick="loadItem(index-1)">◀</button><button onclick="loadItem(index+1)">▶</button><span id="progress">Loading…</span><select id="mode" onchange="reloadImage()"><option value="seg">DAPI segmentation</option><option value="raw">Raw DAPI</option></select><button onclick="resetView()">Reset view</button><div id="meta"></div><button id="buildBtnTop" class="top-build" onclick="buildFigure()" disabled>Make Supplementary Figure S5 · 600 dpi</button><button id="stopBtn" onclick="stopServer()" disabled>STOP SERVER</button></div>
<div id="main"><div id="wrap"><svg id="viewer"><g id="viewport"><image id="bg" x="0" y="0"></image><g id="saved"></g><polyline id="working" fill="none" stroke="#fff" stroke-width="6" stroke-linejoin="round" vector-effect="non-scaling-stroke"></polyline><g id="workingDots"></g></g></svg></div><div id="side">
<div class="section"><b>Status:</b> <span id="state"></span><p class="small">Anatomy is DAPI-only. c-FOS and ACTH/CLIP are intentionally unavailable during delineation.</p></div>
<div class="section"><b>Draw/edit region</b><div class="regions"><button id="arc" class="region active" onclick="setRegion('ARC')">ARC</button><button id="me" class="region" onclick="setRegion('ME')">ME</button><button id="vmn" class="region" onclick="setRegion('VMN')">VMN</button></div><p class="small"><span class="key">Auto-close:</span> after ≥3 vertices, click within 16 screen pixels of the first white dot. Enter or double-click also closes. Lines remain bold while zooming.</p></div>
<div class="section"><button class="action" onclick="finishPoly()">Close polygon</button><button class="action" onclick="undo()">Undo point / last selected polygon</button><button class="action" onclick="clearRegion()">Clear selected region</button><button class="action" onclick="clearAll()">Clear all</button><button class="action" onclick="copyPrevious()">Copy previous</button><div id="polyList" class="small"></div></div>
<div class="section"><textarea id="notes" placeholder="Review notes"></textarea></div><button class="action save" onclick="saveDecision('include',false)">Save included</button><button class="action save" onclick="saveDecision('include',true)">Save & next</button><button class="action exclude" onclick="saveDecision('exclude',true)">Exclude & next</button>
<div class="section"><button id="buildBtnSide" class="action build-side" onclick="buildFigure()" disabled>Make Supplementary Figure S5 · 600 dpi</button><div id="buildStatus" class="small">Loading review state…</div></div><p class="small">At least one ARC, ME, or VMN polygon is required for an included field. Supplementary Figure S5 statistics remain ARC/ME. Shift, middle, or right drag pans; wheel zooms.</p></div></div>
<script>
'use strict';
const regionNames=['ARC','ME','VMN'],colors={ARC:'#00adb5',ME:'#e67e22',VMN:'#be3eae'};
let items=[],index=0,meta=null,polygons=[],points=[],region='ARC',zoom=1,panX=0,panY=0,drag=false,lx=0,ly=0,canBuild=false,building=false;
const svg=document.getElementById('viewer'),vp=document.getElementById('viewport');
function setBuildButtons(enabled){for(const id of ['buildBtnTop','buildBtnSide'])document.getElementById(id).disabled=!(enabled&&!building)}
function updateBuildState(counts){canBuild=counts.pending===0;setBuildButtons(canBuild);document.getElementById('stopBtn').disabled=!canBuild;document.getElementById('buildStatus').textContent=canBuild?'All decisions complete. Ready to quantify and render at 600 dpi.':`${counts.pending} review decision${counts.pending===1?'':'s'} remain. The build control is visible and will enable when complete.`}
async function init(){let response=await fetch('/api/list',{cache:'no-store'});items=await response.json();if(!response.ok)throw Error(items.error||'Could not load review list');let p=items.findIndex(x=>x.status==='pending');await loadItem(p>=0?p:0)}
function setRegion(r){region=r;for(const n of regionNames)document.getElementById(n.toLowerCase()).classList.toggle('active',r===n);points=[];render()}
function transform(){vp.setAttribute('transform',`translate(${panX} ${panY}) scale(${zoom})`);render()}
function resetView(){zoom=1;panX=0;panY=0;transform()}
function svgPoint(e){let p=svg.createSVGPoint();p.x=e.clientX;p.y=e.clientY;return p.matrixTransform(vp.getScreenCTM().inverse())}
function screenPoint(p){let q=svg.createSVGPoint();q.x=p[0];q.y=p[1];return q.matrixTransform(vp.getScreenCTM())}
function nearFirst(e){if(points.length<3)return false;let q=screenPoint(points[0]);return Math.hypot(e.clientX-q.x,e.clientY-q.y)<=16}
function inside(p){return meta&&p.x>=0&&p.y>=0&&p.x<meta.display_width&&p.y<meta.display_height}
function finishPoly(){if(!points.length)return true;if(points.length<3){document.getElementById('state').textContent='A polygon needs at least three vertices.';return false}polygons.push({region,points:points.map(p=>[+p[0].toFixed(3),+p[1].toFixed(3)])});points=[];render();return true}
function undo(){if(points.length)points.pop();else{let i=-1;for(let j=polygons.length-1;j>=0;j--)if(polygons[j].region===region){i=j;break}if(i>=0)polygons.splice(i,1)}render()}
function clearRegion(){polygons=polygons.filter(p=>p.region!==region);points=[];render()}
function clearAll(){if(confirm('Clear all ARC, ME, and VMN polygons?')){polygons=[];points=[];render()}}
function dot(parent,p,color){let c=document.createElementNS('http://www.w3.org/2000/svg','circle');c.setAttribute('cx',p[0]);c.setAttribute('cy',p[1]);let m=vp.getScreenCTM(),screenScale=m?Math.hypot(m.a,m.b):zoom;c.setAttribute('r',String(5/Math.max(screenScale,.001)));c.setAttribute('fill','#fff');c.setAttribute('stroke',color);c.setAttribute('stroke-width','2.5');c.setAttribute('vector-effect','non-scaling-stroke');parent.appendChild(c)}
function render(){let g=document.getElementById('saved');g.innerHTML='';for(const p of polygons){let q=document.createElementNS('http://www.w3.org/2000/svg','polygon');q.setAttribute('points',p.points.map(x=>x.join(',')).join(' '));q.setAttribute('fill',colors[p.region]);q.setAttribute('fill-opacity','.22');q.setAttribute('stroke',colors[p.region]);q.setAttribute('stroke-width','6');q.setAttribute('stroke-linejoin','round');q.setAttribute('vector-effect','non-scaling-stroke');g.appendChild(q);for(const v of p.points)dot(g,v,colors[p.region])}let w=document.getElementById('working');w.setAttribute('points',points.map(x=>x.join(',')).join(' '));w.setAttribute('stroke',colors[region]);let d=document.getElementById('workingDots');d.innerHTML='';for(const p of points)dot(d,p,colors[region]);document.getElementById('polyList').innerHTML=regionNames.map(n=>`${n}: <b>${polygons.filter(p=>p.region===n).length}</b>`).join(' · ')}
svg.addEventListener('contextmenu',e=>e.preventDefault());
svg.addEventListener('wheel',e=>{e.preventDefault();let p=svgPoint(e),f=e.deltaY<0?1.15:1/1.15,o=zoom;zoom=Math.min(15,Math.max(.25,zoom*f));panX+=p.x*(o-zoom);panY+=p.y*(o-zoom);transform()},{passive:false});
svg.addEventListener('pointerdown',e=>{if(e.button===1||e.button===2||e.shiftKey){drag=true;lx=e.clientX;ly=e.clientY;svg.setPointerCapture(e.pointerId);return}if(e.button!==0)return;if(nearFirst(e)){finishPoly();return}let p=svgPoint(e);if(!inside(p))return;if(points.length){let q=screenPoint(points[points.length-1]);if(Math.hypot(e.clientX-q.x,e.clientY-q.y)<2)return}points.push([p.x,p.y]);render()});
svg.addEventListener('pointermove',e=>{if(drag){let r=svg.getBoundingClientRect();panX+=(e.clientX-lx)*(meta.display_width/r.width)/zoom;panY+=(e.clientY-ly)*(meta.display_height/r.height)/zoom;lx=e.clientX;ly=e.clientY;transform()}});
svg.addEventListener('pointerup',e=>{drag=false;try{svg.releasePointerCapture(e.pointerId)}catch(_){}});svg.addEventListener('pointercancel',()=>drag=false);svg.addEventListener('dblclick',e=>{e.preventDefault();finishPoly()});
document.addEventListener('keydown',e=>{if(['TEXTAREA','INPUT','SELECT'].includes(document.activeElement.tagName))return;if(e.key==='Enter'){finishPoly();e.preventDefault()}if(e.key==='Escape'){points=[];render()}if(e.key==='Backspace'){e.preventDefault();undo()}});
async function loadItem(i){if(!items.length)return;index=Math.max(0,Math.min(items.length-1,i));let response=await fetch(`/api/item/${index}`,{cache:'no-store'});meta=await response.json();if(!response.ok)throw Error(meta.error||'Could not load item');polygons=(meta.polygons||[]).filter(p=>regionNames.includes(p.region));points=[];document.getElementById('notes').value=meta.notes||'';document.getElementById('state').textContent=meta.status;document.getElementById('progress').textContent=`${index+1}/${items.length} | annotated ${meta.counts.annotated}, excluded ${meta.counts.excluded}, pending ${meta.counts.pending}`;document.getElementById('meta').textContent=`${meta.acquisition_id} | ${meta.condition}`;updateBuildState(meta.counts);svg.setAttribute('viewBox',`0 0 ${meta.display_width} ${meta.display_height}`);let bg=document.getElementById('bg');bg.setAttribute('width',meta.display_width);bg.setAttribute('height',meta.display_height);reloadImage();resetView();render()}
function reloadImage(){document.getElementById('bg').setAttribute('href',`/api/image/${index}?mode=${document.getElementById('mode').value}&t=${Date.now()}`)}
async function saveDecision(decision,next){if(!finishPoly())return;if(decision==='include'&&!polygons.some(p=>regionNames.includes(p.region))){alert('Draw at least one ARC, ME, or VMN polygon.');return}let body={decision,polygons:decision==='include'?polygons:[],notes:document.getElementById('notes').value,display_width:meta.display_width,display_height:meta.display_height};let r=await fetch(`/api/save/${index}`,{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)});let t=await r.text();if(!r.ok){alert(t);return}items=await(await fetch('/api/list',{cache:'no-store'})).json();await loadItem(next?Math.min(index+1,items.length-1):index)}
async function copyPrevious(){let r=await fetch(`/api/previous/${index}`,{cache:'no-store'});if(!r.ok){alert('No prior included annotation');return}polygons=((await r.json()).polygons||[]).filter(p=>regionNames.includes(p.region));points=[];render()}
async function buildFigure(){if(!canBuild){document.getElementById('buildStatus').textContent='Complete all review decisions before building.';return}building=true;setBuildButtons(false);let s=document.getElementById('buildStatus');s.textContent='Building statistics and bilingual 600-dpi panels…';let r=await fetch('/api/build',{method:'POST'}),t=await r.text();if(!r.ok){s.textContent=t;building=false;setBuildButtons(canBuild);return}let d=JSON.parse(t);s.textContent='Built: '+d.output_dir;building=false;setBuildButtons(canBuild)}
async function stopServer(){let r=await fetch('/api/stop',{method:'POST'});if(r.ok)document.body.innerHTML='<h1 style="padding:40px">Stopped. Supplementary Figure S5 review is complete.</h1>';else alert(await r.text())}
init().catch(e=>{document.getElementById('state').textContent=e.message;document.getElementById('buildStatus').textContent=e.message});
</script></body></html>'''


def _display_image(row: Mapping[str, Any], mode: str, max_dim: int, workdir: Path) -> Tuple[bytes, Tuple[int, int], Tuple[int, int]]:
    if Image is None:
        raise RuntimeError("Pillow is required for --annotate")
    dapi = maximum_projection(read_stack(resolve_work_path(workdir, row["dapi_mip"])))
    shape = tuple(map(int, dapi.shape))
    if mode == "raw":
        u8 = normalize_u8(dapi)
        rgb = np.repeat(u8[..., None], 3, axis=-1)
    else:
        seg = load_segmentation(resolve_work_path(workdir, row["dapi_seg"]))
        if seg.shape != shape:
            raise ValueError(f"DAPI mask/image mismatch: {seg.shape} vs {shape}")
        rgb = np.zeros((*shape, 3), dtype=np.uint8)
        rgb[seg > 0] = 190
        rgb[find_boundaries(seg, mode="inner")] = 255
    scale = min(1.0, float(max_dim) / max(shape))
    display = (max(1, int(round(shape[0] * scale))), max(1, int(round(shape[1] * scale))))
    image = Image.fromarray(rgb)
    if display != shape:
        image = image.resize((display[1], display[0]), Image.Resampling.LANCZOS)
    buf = io.BytesIO()
    image.save(buf, format="PNG", optimize=True)
    return buf.getvalue(), shape, display


def annotate(args: argparse.Namespace) -> int:
    sheet = read_sheet(args)
    p = paths(args)
    for _, row in sheet.iterrows():
        if not resolve_work_path(p["work"], row["dapi_seg"]).exists():
            raise FileNotFoundError(f"DAPI segmentation missing for {row['acquisition_id']}; run --segment")
    try:
        from flask import Flask, Response, jsonify, request
    except Exception as exc:
        raise RuntimeError("Flask is required for --annotate") from exc
    ensure_dir(p["annotations"])
    ensure_dir(p["region_masks"])
    records = sheet.to_dict("records")
    app = Flask(__name__)
    cache: Dict[Tuple[int, str], Tuple[bytes, Tuple[int, int], Tuple[int, int]]] = {}

    def rendered(i: int, mode: str) -> Tuple[bytes, Tuple[int, int], Tuple[int, int]]:
        if (i, mode) not in cache:
            cache[(i, mode)] = _display_image(records[i], mode, args.web_max_dim, p["work"])
        return cache[(i, mode)]

    @app.get("/")
    def root_route():
        response = Response(ANNOTATION_HTML, mimetype="text/html")
        response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
        return response

    @app.get("/api/list")
    def list_route():
        states = annotation_states(args, sheet)
        return jsonify([{"index": i, "acquisition_id": str(r["acquisition_id"]), "status": states[i]} for i, r in enumerate(records)])

    @app.get("/api/item/<int:i>")
    def item_route(i: int):
        if i < 0 or i >= len(records):
            return jsonify({"error": "index"}), 404
        _, native, display = rendered(i, "seg")
        row = records[i]
        acq = str(row["acquisition_id"])
        ann = read_annotation(args, acq) or {}
        sx, sy = display[1] / native[1], display[0] / native[0]
        polygons = [{"region": q.get("region"), "points": [[float(x)*sx, float(y)*sy] for x, y in q.get("points", [])]} for q in ann.get("polygons", [])]
        states = annotation_states(args, sheet)
        counts = {s: states.count(s) for s in ["annotated", "excluded", "pending"]}
        return jsonify({"index": i, "acquisition_id": acq, "condition": str(row.get("condition", "")), "status": states[i], "notes": ann.get("notes", ""), "polygons": polygons, "display_width": display[1], "display_height": display[0], "counts": counts})

    @app.get("/api/image/<int:i>")
    def image_route(i: int):
        if i < 0 or i >= len(records):
            return Response(status=404)
        mode = request.args.get("mode", "seg")
        mode = mode if mode in {"seg", "raw"} else "seg"
        return Response(rendered(i, mode)[0], mimetype="image/png")

    @app.post("/api/save/<int:i>")
    def save_route(i: int):
        incoming = request.get_json(force=True)
        decision = str(incoming.get("decision", "")).lower()
        if decision not in {"include", "exclude"}:
            return jsonify({"error": "decision"}), 400
        _, native, display = rendered(i, "seg")
        sx = native[1] / max(1.0, float(incoming.get("display_width", display[1])))
        sy = native[0] / max(1.0, float(incoming.get("display_height", display[0])))
        polygons = [{"region": q.get("region"), "points": [[float(x)*sx, float(y)*sy] for x, y in q.get("points", [])]} for q in incoming.get("polygons", [])]
        if decision == "include" and not any(str(q.get("region", "")).upper() in REGION_NAMES for q in polygons):
            return jsonify({"error": "At least one ARC, ME, or VMN polygon is required"}), 400
        save_annotation(args, records[i], decision, polygons, str(incoming.get("notes", "")), native)
        return jsonify({"ok": True})

    @app.get("/api/previous/<int:i>")
    def previous_route(i: int):
        _, native, display = rendered(i, "seg")
        for j in range(i-1, -1, -1):
            ann = read_annotation(args, str(records[j]["acquisition_id"]))
            if ann and ann.get("decision") == "include" and ann.get("polygons"):
                sh, sw = float(ann.get("image_height", native[0])), float(ann.get("image_width", native[1]))
                out = [{"region": q.get("region"), "points": [[float(x)/sw*display[1], float(y)/sh*display[0]] for x, y in q.get("points", [])]} for q in ann["polygons"]]
                return jsonify({"polygons": out})
        return jsonify({"error": "none"}), 404

    @app.post("/api/build")
    def build_route():
        states = annotation_states(args, sheet)
        pending = states.count("pending")
        if pending:
            return jsonify({"error": f"{pending} review decisions remain"}), 409
        p = paths(args)
        result_sentinel = p["results"] / "per_analysis_unit_region_summary.csv"
        if result_sentinel.exists():
            return jsonify({
                "error": (
                    f"Regional results already exist and were preserved: {result_sentinel}. "
                    "For a figure-only retry, complete the DAPI third-ventricle review "
                    "at http://127.0.0.1:34250/ "
                    "and run FigS5/04_run_figure_s5.sh render with a fresh FIGS5_OUTPUT_DIR. "
                    "A new quantification requires a fresh FIGS5_WORKDIR."
                )
            }), 409
        paper_root = Path(__file__).resolve().parent.parent
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        output_dir = paper_root / "analyses" / "FigS5" / "figure" / f"browser_candidate_{stamp}"
        quantify_command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--quantify",
            "--export-dir", str(p["export"]),
            "--channel-map", str(p["channel_map"]),
            "--sample-manifest", str(p["manifest"]),
            "--workdir", str(p["work"]),
            "--qc-format", "both",
            "--qc-dpi", "600",
            "--plot-dpi", "600",
            "--acth-clip-dapi-area-quantile", str(args.acth_clip_dapi_area_quantile),
        ]
        render_command = [
            sys.executable,
            str(paper_root / "FigS5" / "03_make_figure_s5_acth_clip_cfos.py"),
            "--workdir", str(p["work"]),
            "--output-dir", str(output_dir),
            "--dpi", "600",
        ]
        logs: List[str] = []
        for command in (quantify_command, render_command):
            completed = subprocess.run(
                command,
                cwd=str(paper_root),
                text=True,
                capture_output=True,
                timeout=1800,
            )
            logs.append(completed.stdout)
            logs.append(completed.stderr)
            if completed.returncode:
                return jsonify({
                    "error": "Supplementary Figure S5 build failed",
                    "command": command,
                    "returncode": completed.returncode,
                    "log": "\n".join(logs)[-12000:],
                }), 500
        return jsonify({
            "ok": True,
            "output_dir": str(output_dir),
            "figure_png": str(output_dir / "Figure_S5.png"),
            "figure_pdf": str(output_dir / "Figure_S5.pdf"),
            "log": "\n".join(logs)[-12000:],
        })

    @app.post("/api/stop")
    def stop_route():
        states = annotation_states(args, sheet)
        if states.count("pending"):
            return jsonify({"error": f"{states.count('pending')} pending"}), 409
        (paths(args)["work"] / "ANNOTATION_COMPLETE.txt").write_text(f"completed_at_utc={utc_now()}\n", encoding="utf-8")
        def shutdown() -> None:
            time.sleep(0.35)
            os.kill(os.getpid(), signal.SIGINT)
        threading.Thread(target=shutdown, daemon=True).start()
        return jsonify({"ok": True})

    print(f"DAPI-only annotation server: http://{args.host}:{args.port}")
    try:
        app.run(host=args.host, port=args.port, threaded=True, use_reloader=False)
    except KeyboardInterrupt:
        print("Annotation server stopped")
    return 0


def status(args: argparse.Namespace) -> int:
    p = paths(args)
    payload: Dict[str, Any] = {"workdir": str(p["work"]), "prepared": p["sheet"].exists()}
    if p["sheet"].exists():
        sheet = read_sheet(args)
        payload["analysis_units"] = int(len(sheet))
        payload["segmentations"] = {m: int(sum(resolve_work_path(p["work"], r[c]).exists() for _, r in sheet.iterrows())) for m, c in [("DAPI", "dapi_seg"), ("cFOS", "cfos_seg"), ("ACTH_CLIP", "acth_clip_seg")]}
        states = annotation_states(args, sheet)
        payload["annotations"] = {s: states.count(s) for s in ["annotated", "excluded", "pending"]}
    print(json.dumps(payload, indent=2))
    return 0


# -----------------------------------------------------------------------------
# Nucleus-centric quantification
# -----------------------------------------------------------------------------

def assign_cfos_rois(dapi: np.ndarray, cfos: np.ndarray, threshold: float, min_area: int) -> Tuple[Set[int], List[Dict[str, Any]]]:
    """Strictly assign each c-FOS ROI to its one best directly-overlapped nucleus."""
    positives: Set[int] = set()
    rows: List[Dict[str, Any]] = []
    max_dapi = int(dapi.max())
    for prop in regionprops(cfos):
        area = int(prop.area)
        best = overlap_px = 0
        fraction = 0.0
        accepted = False
        reason = "too_small" if area < min_area else "no_dapi_overlap"
        if area >= min_area and max_dapi > 0:
            vals = dapi[prop.coords[:, 0], prop.coords[:, 1]]
            counts = np.bincount(vals.astype(np.int64), minlength=max_dapi + 1)
            counts[0] = 0
            best = int(np.argmax(counts)) if counts.size else 0
            overlap_px = int(counts[best]) if best else 0
            fraction = overlap_px / float(max(1, area))
            accepted = bool(best > 0 and fraction >= threshold)
            reason = "accepted_strict_direct_overlap" if accepted else "overlap_below_threshold"
        if accepted:
            positives.add(best)
        rows.append({
            "marker": "cFOS", "roi_label": int(prop.label), "roi_area_px": area,
            "assigned_nucleus": best if accepted else 0, "best_overlap_px": overlap_px,
            "roi_overlap_fraction": fraction, "accepted": accepted, "reason": reason,
        })
    return positives, rows


def assign_acth_clip_rois(
    dapi: np.ndarray,
    acth_clip: np.ndarray,
    expansion_px: int,
    min_pixels: int,
    dapi_area_quantile: float = ACTH_CLIP_DAPI_AREA_QUANTILE,
) -> Tuple[Set[int], List[Dict[str, Any]]]:
    """Filter likely fragments by DAPI size, then assign each ACTH/CLIP ROI once."""
    if not 0.0 <= float(dapi_area_quantile) <= 0.5:
        raise ValueError("ACTH/CLIP DAPI-area quantile must be between 0 and 0.5")
    dapi_areas = np.asarray([int(prop.area) for prop in regionprops(dapi)], dtype=float)
    if dapi_areas.size == 0:
        raise ValueError("Cannot define ACTH/CLIP size filter without DAPI nuclei")
    dapi_area_threshold = float(
        np.quantile(dapi_areas, float(dapi_area_quantile), method="higher")
    )
    dapi_area_median = float(np.median(dapi_areas))

    positives: Set[int] = set()
    rows: List[Dict[str, Any]] = []
    territories = expand_labels(dapi.astype(np.int32), distance=max(0, int(expansion_px)))
    max_dapi = int(dapi.max())
    for prop in regionprops(acth_clip):
        area = int(prop.area)
        size_filter_pass = bool(area >= dapi_area_threshold)
        best = 0
        overlap_px = 0
        if size_filter_pass:
            vals = territories[prop.coords[:, 0], prop.coords[:, 1]]
            counts = np.bincount(vals.astype(np.int64), minlength=max_dapi + 1)
            if counts.size:
                counts[0] = 0
            best = int(np.argmax(counts)) if counts.size else 0
            overlap_px = int(counts[best]) if best else 0
        accepted = bool(size_filter_pass and best > 0 and overlap_px >= int(min_pixels))
        if accepted:
            positives.add(best)
        if not size_filter_pass:
            reason = "below_acquisition_dapi_area_p05_likely_noncell"
        elif accepted:
            reason = "accepted_size_and_best_expanded_dapi"
        else:
            reason = "insufficient_expanded_dapi_pixels"
        rows.append({
            "marker": "ACTH_CLIP",
            "marker_display": ACTH_CLIP_DISPLAY,
            "roi_label": int(prop.label),
            "roi_area_px": area,
            "dapi_area_reference_n": int(dapi_areas.size),
            "dapi_area_quantile": float(dapi_area_quantile),
            "dapi_area_threshold_px": dapi_area_threshold,
            "dapi_area_median_px": dapi_area_median,
            "size_filter_pass": size_filter_pass,
            "likely_noncell_by_size": not size_filter_pass,
            "assigned_nucleus": best if accepted else 0,
            "best_overlap_px": overlap_px,
            "roi_overlap_fraction": overlap_px / float(max(1, area)),
            "accepted": accepted,
            "reason": reason,
            "dapi_expansion_px": int(expansion_px),
        })
    return positives, rows


def assign_nucleus_region(prop: Any, region_mask: np.ndarray, minimum_fraction: float) -> Tuple[str, float, bool]:
    vals = region_mask[prop.coords[:, 0], prop.coords[:, 1]].astype(np.int32)
    counts = np.bincount(vals, minlength=max(REGION_CODES.values()) + 1)
    counts[0] = 0
    best = int(np.argmax(counts)) if counts.sum() else 0
    fraction = float(counts[best]) / float(max(1, prop.area)) if best else 0.0
    boundary = bool(np.unique(vals[vals > 0]).size > 1 or (best > 0 and fraction < minimum_fraction))
    if best and fraction >= minimum_fraction:
        return CODE_TO_REGION[best], fraction, boundary
    y, x = prop.centroid
    code = int(region_mask[int(np.clip(round(y), 0, region_mask.shape[0]-1)), int(np.clip(round(x), 0, region_mask.shape[1]-1))])
    return CODE_TO_REGION.get(code, "OUTSIDE"), fraction, boundary


def phenotype(cfos: bool, acth_clip: bool) -> str:
    if cfos and acth_clip:
        return "cFOS_ACTH_CLIP"
    if cfos:
        return "cFOS_only"
    if acth_clip:
        return "ACTH_CLIP_only"
    return "neither"


def safe_ratio(num: Any, den: Any) -> float:
    try:
        n, d = float(num), float(den)
    except Exception:
        return float("nan")
    return n / d if np.isfinite(n) and np.isfinite(d) and d > 0 else float("nan")


def add_endpoints(row: Dict[str, Any]) -> Dict[str, Any]:
    if not bool(row.get("region_present", True)):
        for endpoint in ENDPOINTS:
            if endpoint not in COUNT_COLUMNS:
                row[endpoint] = np.nan
        return row
    row["double_over_total_acth_clip"] = safe_ratio(row["n_double"], row["n_total_acth_clip"])
    row["total_cfos_over_dapi"] = safe_ratio(row["n_total_cfos"], row["n_dapi"])
    row["total_acth_clip_over_dapi"] = safe_ratio(row["n_total_acth_clip"], row["n_dapi"])
    row["double_over_dapi"] = safe_ratio(row["n_double"], row["n_dapi"])
    row["cfos_only_over_dapi"] = safe_ratio(row["n_cfos_only"], row["n_dapi"])
    row["acth_clip_only_over_dapi"] = safe_ratio(row["n_acth_clip_only"], row["n_dapi"])
    row["double_over_total_cfos"] = safe_ratio(row["n_double"], row["n_total_cfos"])
    row["double_per_100k_um2"] = safe_ratio(row["n_double"] * 100000.0, row["area_um2"])
    row["total_cfos_per_100k_um2"] = safe_ratio(row["n_total_cfos"] * 100000.0, row["area_um2"])
    row["total_acth_clip_per_100k_um2"] = safe_ratio(row["n_total_acth_clip"] * 100000.0, row["area_um2"])
    return row


def count_summary(nuclei: pd.DataFrame, region: str, present: bool, area_px: int, px_x: float, px_y: float) -> Dict[str, Any]:
    sub = nuclei[nuclei["region"] == region] if not nuclei.empty else nuclei
    states = sub["phenotype"].value_counts().to_dict() if not sub.empty else {}
    row: Dict[str, Any] = {
        "region": region, "region_present": bool(present), "area_px": int(area_px),
        "area_um2": float(area_px) * float(px_x) * float(px_y),
        "n_dapi": int(len(sub)), "n_cfos_only": int(states.get("cFOS_only", 0)),
        "n_acth_clip_only": int(states.get("ACTH_CLIP_only", 0)), "n_double": int(states.get("cFOS_ACTH_CLIP", 0)),
        "n_neither": int(states.get("neither", 0)),
    }
    row["n_total_cfos"] = row["n_cfos_only"] + row["n_double"]
    row["n_total_acth_clip"] = row["n_acth_clip_only"] + row["n_double"]
    row["invariant_partition_ok"] = row["n_cfos_only"] + row["n_acth_clip_only"] + row["n_double"] + row["n_neither"] == row["n_dapi"]
    row["invariant_total_cfos_ok"] = row["n_total_cfos"] == row["n_cfos_only"] + row["n_double"]
    row["invariant_total_acth_clip_ok"] = row["n_total_acth_clip"] == row["n_acth_clip_only"] + row["n_double"]
    row["invariant_double_bounds_ok"] = row["n_double"] <= min(row["n_total_cfos"], row["n_total_acth_clip"])
    row["invariants_ok"] = all(row[x] for x in ["invariant_partition_ok", "invariant_total_cfos_ok", "invariant_total_acth_clip_ok", "invariant_double_bounds_ok"])
    return add_endpoints(row)


def render_quant_qc(args: argparse.Namespace, row: Mapping[str, Any], dapi: np.ndarray, region_mask: np.ndarray, nuclei: pd.DataFrame, outdir: Path) -> None:
    shape = dapi.shape
    rgb = np.zeros((*shape, 3), dtype=np.float32)
    rgb[dapi > 0] = 0.65
    rgb[find_boundaries(dapi, mode="inner")] = 0.92
    colors = {"neither": (0.72, 0.72, 0.72), "cFOS_only": (0.92, 0.05, 0.75), "ACTH_CLIP_only": (0.95, 0.12, 0.05), "cFOS_ACTH_CLIP": (0.0, 0.95, 0.95)}
    for state, color in colors.items():
        labs = nuclei.loc[nuclei["phenotype"] == state, "nucleus_label"].to_numpy(dtype=int) if not nuclei.empty else np.array([], dtype=int)
        if labs.size:
            rgb[np.isin(dapi, labs)] = color
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.imshow(rgb)
    for region in REGION_NAMES:
        mask = region_mask == REGION_CODES[region]
        if mask.any():
            ax.contour(mask, [0.5], colors=[REGION_COLORS[region]], linewidths=2.2)
    px = float(row.get("pixel_size_x_um", np.nan))
    if np.isfinite(px) and px > 0:
        bar_px = int(round(args.scalebar_um / px))
        if 0 < bar_px < shape[1] * 0.5:
            x0, y0 = shape[1] - bar_px - 35, shape[0] - 45
            ax.add_patch(Rectangle((x0, y0), bar_px, 8, facecolor="white", edgecolor="black"))
            ax.text(x0 + bar_px/2, y0-5, f"{args.scalebar_um:g} µm", color="white", ha="center", va="bottom", fontweight="bold")
    ax.set_title(f"{row['analysis_unit_id']} | {row['condition']} | DAPI-only HIL anatomy", fontweight="bold")
    ax.axis("off")
    ax.legend(handles=[Patch(facecolor=colors[x], label=x) for x in colors] + [Patch(facecolor="none", edgecolor=REGION_COLORS[x], linewidth=2, label=x) for x in REGION_NAMES], loc="upper right", fontsize=8)
    ensure_dir(outdir)
    stem = safe_name(row["analysis_unit_id"]) + "_quantification_qc"
    if args.qc_format in {"png", "both"}:
        fig.savefig(outdir / f"{stem}.png", dpi=args.qc_dpi, bbox_inches="tight")
    if args.qc_format in {"pdf", "both"}:
        fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def cliffs_delta(a: Sequence[float], b: Sequence[float]) -> float:
    x = np.asarray(a, dtype=float); y = np.asarray(b, dtype=float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if not x.size or not y.size:
        return float("nan")
    return float((sum(v > w for v in x for w in y) - sum(v < w for v in x for w in y)) / (x.size * y.size))


def safe_anova(groups: Sequence[np.ndarray]) -> float:
    usable = [np.asarray(g, dtype=float)[np.isfinite(g)] for g in groups]
    usable = [g for g in usable if g.size >= 2]
    try:
        return float(f_oneway(*usable).pvalue) if len(usable) >= 2 else np.nan
    except Exception:
        return np.nan


def safe_kruskal(groups: Sequence[np.ndarray]) -> float:
    usable = [np.asarray(g, dtype=float)[np.isfinite(g)] for g in groups]
    usable = [g for g in usable if g.size >= 2]
    try:
        return float(kruskal(*usable).pvalue) if len(usable) >= 2 else np.nan
    except Exception:
        return np.nan


def safe_mwu(a: Sequence[float], b: Sequence[float]) -> float:
    x = np.asarray(a, dtype=float); y = np.asarray(b, dtype=float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    try:
        return float(mannwhitneyu(x, y, alternative="two-sided").pvalue) if x.size >= 2 and y.size >= 2 else np.nan
    except Exception:
        return np.nan


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float); out = np.full(p.shape, np.nan); valid = np.where(np.isfinite(p))[0]
    if not valid.size:
        return out
    order = valid[np.argsort(p[valid])]; ranked = p[order] * len(order) / np.arange(1, len(order)+1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.clip(ranked, 0, 1)
    return out


def iqr_mask(values: Sequence[float], k: float = 1.5) -> np.ndarray:
    x = np.asarray(values, dtype=float); out = np.zeros(x.shape, dtype=bool); finite = np.isfinite(x)
    if finite.sum() < 4:
        return out
    q1, q3 = np.percentile(x[finite], [25, 75]); span = q3-q1
    if not np.isfinite(span) or span <= 0:
        return out
    out[finite] = (x[finite] < q1-k*span) | (x[finite] > q3+k*span)
    return out


def build_statistics(unit_df: pd.DataFrame, endpoints: Sequence[str]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    conds = [c for c in BIOLOGICAL_CONDITIONS if c in set(unit_df["condition"])]
    for endpoint in endpoints:
        for region in ["ME", "ARC"]:
            sub = unit_df[(unit_df["region"] == region) & unit_df["region_present"].map(bool_value)]
            groups = {c: pd.to_numeric(sub.loc[sub["condition"] == c, endpoint], errors="coerce").dropna().to_numpy(float) for c in conds}
            rows.append({"endpoint": endpoint, "region": region, "test": "Kruskal-Wallis", "group1": "ALL", "group2": "", "n1": sum(len(x) for x in groups.values()), "n2": np.nan, "p_value": safe_kruskal(list(groups.values())), "cliffs_delta": np.nan, "analysis_unit": "acquisition"})
            rows.append({"endpoint": endpoint, "region": region, "test": "ANOVA", "group1": "ALL", "group2": "", "n1": sum(len(x) for x in groups.values()), "n2": np.nan, "p_value": safe_anova(list(groups.values())), "cliffs_delta": np.nan, "analysis_unit": "acquisition"})
            for a, b in combinations(conds, 2):
                rows.append({"endpoint": endpoint, "region": region, "test": "Mann-Whitney U", "group1": a, "group2": b, "n1": len(groups[a]), "n2": len(groups[b]), "p_value": safe_mwu(groups[a], groups[b]), "cliffs_delta": cliffs_delta(groups[a], groups[b]), "analysis_unit": "acquisition", "insufficient_n": len(groups[a]) < 3 or len(groups[b]) < 3})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value_bh"] = bh_adjust(out["p_value"])
    return out


def plot_endpoint(unit_df: pd.DataFrame, endpoint: str, outdir: Path, dpi: int, iqr_k: float) -> None:
    ylabel, title = ENDPOINTS[endpoint]
    conds = [c for c in BIOLOGICAL_CONDITIONS if c in set(unit_df["condition"])]
    if not conds:
        return
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 5.4), squeeze=False)
    axes = axes.ravel(); rng = np.random.default_rng(17)
    all_vals = pd.to_numeric(unit_df[endpoint], errors="coerce").dropna().to_numpy(float)
    global_max = float(np.max(all_vals)) if all_vals.size else 1.0
    if not np.isfinite(global_max) or global_max <= 0:
        global_max = 1.0
    for ax, region in zip(axes, ["ME", "ARC"]):
        sub = unit_df[(unit_df["region"] == region) & unit_df["region_present"].map(bool_value)]
        groups = [pd.to_numeric(sub.loc[sub["condition"] == c, endpoint], errors="coerce").dropna().to_numpy(float) for c in conds]
        means = [float(np.mean(x)) if x.size else np.nan for x in groups]
        sds = [float(np.std(x, ddof=1)) if x.size >= 2 else 0.0 for x in groups]
        x = np.arange(len(conds))
        ax.bar(x, means, yerr=sds, width=0.52, color=[COND_COLORS[c] for c in conds], edgecolor="black", linewidth=2, capsize=5, zorder=2)
        data_max = 0.0
        for i, (cond, values) in enumerate(zip(conds, groups)):
            outliers = iqr_mask(values, iqr_k); jitter = rng.uniform(-0.07, 0.07, values.size)
            for j, value in enumerate(values):
                ax.scatter(i+jitter[j], value, marker="x" if outliers[j] else "o", s=55, color="black" if outliers[j] else COND_COLORS[cond], edgecolor="black", linewidth=1.5, zorder=4)
            if values.size:
                data_max = max(data_max, float(np.max(values)))
                ax.text(i, max(float(np.max(values)), means[i]+sds[i]) + 0.04*global_max, f"n={values.size}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        stats_lines = [f"KW p={safe_kruskal(groups):.3g}" if np.isfinite(safe_kruskal(groups)) else "KW p=NA", f"ANOVA p={safe_anova(groups):.3g}" if np.isfinite(safe_anova(groups)) else "ANOVA p=NA"]
        for i, j in combinations(range(len(conds)), 2):
            pv = safe_mwu(groups[i], groups[j]); stats_lines.append(f"MWU {conds[i]} vs {conds[j]} p={pv:.3g}" if np.isfinite(pv) else f"MWU {conds[i]} vs {conds[j]} p=NA")
        ax.set_ylim(0, max(global_max*1.65, data_max*1.35, 0.1))
        ax.text(0.98, 0.98, "\n".join(stats_lines), transform=ax.transAxes, ha="right", va="top", fontsize=8.2, fontweight="bold", bbox={"facecolor":"white","edgecolor":"none","alpha":0.72})
        ax.set_title(region, fontsize=16, fontweight="bold"); ax.set_xticks(x); ax.set_xticklabels(conds, fontweight="bold"); ax.set_ylabel(ylabel, fontweight="bold")
        for spine in ax.spines.values(): spine.set_linewidth(2)
    handles = [Line2D([0],[0],marker="o",linestyle="none",markerfacecolor=COND_COLORS[c],markeredgecolor="black",label=c) for c in conds]
    handles.append(Line2D([0],[0],marker="x",linestyle="none",color="black",label="IQR outlier"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5,0.92), ncol=len(handles), frameon=False)
    fig.suptitle(f"{title} per acquisition", fontsize=20, fontweight="bold", y=1.01)
    fig.tight_layout(rect=[0,0,1,0.84]); ensure_dir(outdir); stem=safe_name(endpoint)
    fig.savefig(outdir/f"ME_ARC_{stem}.pdf", bbox_inches="tight"); fig.savefig(outdir/f"ME_ARC_{stem}.png", dpi=dpi, bbox_inches="tight"); plt.close(fig)


def audit_segmentations(args: argparse.Namespace) -> int:
    """Audit imported masks and the acquisition-specific ACTH/CLIP size rule."""
    sheet = read_sheet(args)
    rows: List[Dict[str, Any]] = []
    for _, srow in sheet.iterrows():
        acq = str(srow["acquisition_id"])
        dapi = load_segmentation(resolve_work_path(paths(args)["work"], srow["dapi_seg"]))
        acth_clip = load_segmentation(resolve_work_path(paths(args)["work"], srow["acth_clip_seg"]))
        if dapi.shape != acth_clip.shape:
            raise ValueError(f"Segmentation shape mismatch for {acq}: {dapi.shape} vs {acth_clip.shape}")
        px_x = float(srow.get("pixel_size_x_um", np.nan))
        expansion_px = int(round(args.acth_clip_expand_um / px_x)) if np.isfinite(px_x) and px_x > 0 else int(round(args.acth_clip_expand_um))
        _positives, assignments = assign_acth_clip_rois(
            dapi,
            acth_clip,
            expansion_px,
            args.acth_clip_min_pixels,
            args.acth_clip_dapi_area_quantile,
        )
        rows.append({
            "analysis_unit_id": acq,
            "acquisition_id": acq,
            "sample_id": srow.get("sample_id", acq),
            "animal_id": srow.get("animal_id", ""),
            "condition": canonical_condition(srow.get("condition", "")),
            "genotype": srow.get("genotype", ""),
            "dapi_reference_n": assignments[0]["dapi_area_reference_n"] if assignments else int(dapi.max()),
            "dapi_area_quantile": args.acth_clip_dapi_area_quantile,
            "dapi_area_threshold_px": assignments[0]["dapi_area_threshold_px"] if assignments else np.nan,
            "dapi_area_median_px": assignments[0]["dapi_area_median_px"] if assignments else np.nan,
            "n_acth_clip_rois_raw": len(assignments),
            "n_rejected_likely_noncell_by_size": sum(bool(item["likely_noncell_by_size"]) for item in assignments),
            "n_passing_size_rule": sum(bool(item["size_filter_pass"]) for item in assignments),
            "n_accepted_after_size_and_dapi_assignment": sum(bool(item["accepted"]) for item in assignments),
            "source_segmentation": srow.get("acth_clip_seg_source", ""),
            "source_segmentation_sha256": srow.get("acth_clip_seg_sha256", ""),
            "source_model_path": srow.get("acth_clip_seg_model_path", ""),
            "raw_masks_modified": False,
            "rule": "reject ROI area below acquisition-specific DAPI nuclear-area fifth percentile",
        })
    audit = pd.DataFrame(rows)
    audit["source_model_consistent_across_cohort"] = audit["source_model_path"].nunique() == 1
    output = paths(args)["work"] / "acth_clip_segmentation_size_audit.csv"
    audit.to_csv(output, index=False)
    print(audit[[
        "sample_id", "condition", "dapi_area_threshold_px", "n_acth_clip_rois_raw",
        "n_rejected_likely_noncell_by_size", "n_passing_size_rule", "source_model_path",
    ]].to_string(index=False))
    print(
        f"ACTH/CLIP size audit: rejected "
        f"{int(audit['n_rejected_likely_noncell_by_size'].sum())}/"
        f"{int(audit['n_acth_clip_rois_raw'].sum())} likely non-cell objects; {output}"
    )
    return 0


def quantify(args: argparse.Namespace) -> int:
    sheet = read_sheet(args)
    p = paths(args); result_dir = ensure_dir(p["results"]); qc_dir = ensure_dir(result_dir/"qc_overlays"); plot_dir = ensure_dir(result_dir/"plots")
    nucleus_rows: List[Dict[str, Any]] = []; image_rows: List[Dict[str, Any]] = []; roi_rows: List[Dict[str, Any]] = []; annotation_rows: List[Dict[str, Any]] = []
    acth_clip_size_audit_rows: List[Dict[str, Any]] = []
    for _, srow in sheet.iterrows():
        acq = str(srow["acquisition_id"]); ann = read_annotation(args, acq)
        ann_status = "pending" if ann is None else str(ann.get("decision", "pending"))
        annotation_rows.append({"analysis_unit_id":acq,"acquisition_id":acq,"animal_id":srow.get("animal_id",""),"condition":srow.get("condition",""),"status":ann_status,"annotation_json":portable_work_path(p["work"], annotation_json_path(args,acq)),"region_mask":portable_work_path(p["work"], region_mask_path(args,acq))})
        if ann_status != "include":
            continue
        required = [resolve_work_path(p["work"], srow[c]) for c in ["dapi_seg","cfos_seg","acth_clip_seg"]]
        if not all(x.exists() for x in required):
            raise FileNotFoundError(f"Missing segmentation for {acq}: {required}")
        dapi, cfos, acth_clip = [load_segmentation(x) for x in required]
        if len({dapi.shape,cfos.shape,acth_clip.shape}) != 1:
            raise ValueError(f"Segmentation shape mismatch for {acq}: {dapi.shape}, {cfos.shape}, {acth_clip.shape}")
        mask_path=region_mask_path(args,acq)
        region_mask=np.asarray(tifffile.imread(str(mask_path)),dtype=np.uint8) if mask_path.exists() else rasterize_polygons(ann.get("polygons",[]),dapi.shape)
        if region_mask.shape != dapi.shape:
            raise ValueError(f"Region-mask shape mismatch for {acq}: {region_mask.shape} vs {dapi.shape}")
        cfos_pos, cfos_assign = assign_cfos_rois(dapi,cfos,args.cfos_dapi_overlap,args.cfos_min_area_px)
        px_x=float(srow.get("pixel_size_x_um",np.nan)); px_y=float(srow.get("pixel_size_y_um",px_x)); expand_px=int(round(args.acth_clip_expand_um/px_x)) if np.isfinite(px_x) and px_x>0 else int(round(args.acth_clip_expand_um))
        acth_clip_pos, acth_clip_assign = assign_acth_clip_rois(
            dapi,
            acth_clip,
            expand_px,
            args.acth_clip_min_pixels,
            args.acth_clip_dapi_area_quantile,
        )
        for item in cfos_assign+acth_clip_assign:
            roi_rows.append({"analysis_unit_id":acq,"acquisition_id":acq,**item})
        acth_rows = [item for item in acth_clip_assign if item["marker"] == "ACTH_CLIP"]
        acth_clip_size_audit_rows.append({
            "analysis_unit_id": acq,
            "acquisition_id": acq,
            "animal_id": srow.get("animal_id", ""),
            "condition": canonical_condition(srow.get("condition", "")),
            "marker": ACTH_CLIP_DISPLAY,
            "rule": "reject ROI area below acquisition DAPI area 5th percentile",
            "dapi_area_quantile": args.acth_clip_dapi_area_quantile,
            "dapi_area_threshold_px": acth_rows[0]["dapi_area_threshold_px"] if acth_rows else np.nan,
            "dapi_area_median_px": acth_rows[0]["dapi_area_median_px"] if acth_rows else np.nan,
            "dapi_reference_n": acth_rows[0]["dapi_area_reference_n"] if acth_rows else int(dapi.max()),
            "n_acth_clip_rois_raw": len(acth_rows),
            "n_rejected_likely_noncell_by_size": sum(bool(item["likely_noncell_by_size"]) for item in acth_rows),
            "n_passing_size_rule": sum(bool(item["size_filter_pass"]) for item in acth_rows),
            "n_accepted_after_size_and_dapi_assignment": sum(bool(item["accepted"]) for item in acth_rows),
            "source_segmentation": srow.get("acth_clip_seg_source", ""),
            "source_segmentation_sha256": srow.get("acth_clip_seg_sha256", ""),
            "source_model_path": srow.get("acth_clip_seg_model_path", ""),
        })
        local_nuclei: List[Dict[str,Any]]=[]
        for prop in regionprops(dapi):
            lab=int(prop.label); region,fraction,boundary=assign_nucleus_region(prop,region_mask,args.region_nucleus_overlap); is_cfos=lab in cfos_pos; is_acth_clip=lab in acth_clip_pos
            rec={"analysis_unit_id":acq,"acquisition_id":acq,"sample_id":srow.get("sample_id",acq),"animal_id":srow.get("animal_id",""),"condition":canonical_condition(srow.get("condition","")),"nucleus_label":lab,"centroid_y":float(prop.centroid[0]),"centroid_x":float(prop.centroid[1]),"nucleus_area_px":int(prop.area),"region":region,"region_overlap_fraction":fraction,"boundary_nucleus":boundary,"is_cfos":is_cfos,"is_acth_clip":is_acth_clip,"phenotype":phenotype(is_cfos,is_acth_clip)}
            local_nuclei.append(rec); nucleus_rows.append(rec)
        local_df=pd.DataFrame(local_nuclei)
        metadata={"analysis_unit_id":acq,"acquisition_id":acq,"sample_id":srow.get("sample_id",acq),"animal_id":srow.get("animal_id",""),"condition":canonical_condition(srow.get("condition","")),"fasting_hours":srow.get("fasting_hours",16),"prior_sugar_exposure":srow.get("prior_sugar_exposure","No"),"exposure_history":srow.get("exposure_history","first-ever solution exposure"),"antibody_display_name":ACTH_CLIP_DISPLAY,"pixel_size_x_um":px_x,"pixel_size_y_um":px_y,"z_planes":srow.get("z_planes",np.nan),"z_step_um":srow.get("z_step_um",np.nan),"z_span_um":srow.get("z_span_um",np.nan),"duplicate_check_status":srow.get("duplicate_check_status","NOT_RUN"),"registration_status":srow.get("registration_status","")}
        for region in ["ARC","ME"]:
            present=bool(ann.get("arc_present" if region=="ARC" else "me_present",False)); summary=count_summary(local_df,region,present,int(np.count_nonzero(region_mask==REGION_CODES[region])),px_x,px_y); image_rows.append({**metadata,**summary,"n_boundary_nuclei":int(local_df.loc[local_df["region"]==region,"boundary_nucleus"].sum()) if not local_df.empty else 0,"n_accepted_cfos_rois":sum(bool(x["accepted"]) for x in cfos_assign),"n_accepted_acth_clip_rois":sum(bool(x["accepted"]) for x in acth_clip_assign)})
        render_quant_qc(args,metadata,dapi,region_mask,local_df,qc_dir)
    pd.DataFrame(annotation_rows).to_csv(result_dir/"annotation_status.csv",index=False)
    nuclei_df=pd.DataFrame(nucleus_rows); images_df=pd.DataFrame(image_rows); rois_df=pd.DataFrame(roi_rows)
    nuclei_df.to_csv(result_dir/"per_nucleus_assignments.csv",index=False); images_df.to_csv(result_dir/"per_image_region_summary.csv",index=False); rois_df.to_csv(result_dir/"marker_roi_assignments.csv",index=False)
    pd.DataFrame(acth_clip_size_audit_rows).to_csv(
        result_dir/"acth_clip_dapi_size_filter_audit.csv", index=False
    )
    if images_df.empty:
        raise RuntimeError("No included annotated acquisitions were available for quantification")
    # acquisition_id is intentionally the analysis unit. Reused animal IDs are metadata only.
    unit_df=images_df.copy(); unit_df.to_csv(result_dir/"per_analysis_unit_region_summary.csv",index=False)
    # Compatibility filename: rows remain acquisition-level and are never collapsed by animal_id.
    unit_df.assign(summary_unit="acquisition_not_collapsed_by_animal_id").to_csv(result_dir/"per_animal_region_summary.csv",index=False)
    invariants=unit_df[["analysis_unit_id","acquisition_id","animal_id","condition","region","region_present","invariant_partition_ok","invariant_total_cfos_ok","invariant_total_acth_clip_ok","invariant_double_bounds_ok","invariants_ok"]].copy(); invariants.to_csv(result_dir/"count_invariants.csv",index=False)
    if not invariants["invariants_ok"].all():
        raise RuntimeError(f"Count invariants failed; see {result_dir/'count_invariants.csv'}")
    endpoints=[x.strip() for x in args.endpoints.split(",") if x.strip()]
    unknown=[x for x in endpoints if x not in ENDPOINTS]
    if unknown: raise ValueError(f"Unknown endpoint(s): {unknown}; choices={sorted(ENDPOINTS)}")
    long_rows=[]
    for _,r in unit_df.iterrows():
        for endpoint in endpoints:
            long_rows.append({"analysis_unit_id":r["analysis_unit_id"],"acquisition_id":r["acquisition_id"],"animal_id":r["animal_id"],"condition":r["condition"],"region":r["region"],"region_present":r["region_present"],"endpoint":endpoint,"endpoint_label":ENDPOINTS[endpoint][0],"value":r.get(endpoint,np.nan),"analysis_unit":"acquisition"})
    pd.DataFrame(long_rows).to_csv(result_dir/"endpoint_values_long.csv",index=False)
    stats=build_statistics(unit_df,endpoints); stats.to_csv(result_dir/"statistics_all_values.csv",index=False)
    for endpoint in endpoints: plot_endpoint(unit_df,endpoint,plot_dir,args.plot_dpi,args.iqr_k)
    (result_dir/"analysis_readme.txt").write_text("Statistical unit: acquisition_id. animal_id is metadata only and reused IDs are not collapsed.\nExperiment: 16 h fasting followed by first-ever acute solution exposure.\nAnatomical annotation used DAPI only. c-FOS and ACTH/CLIP never influenced ARC/ME/VMN placement; inferential endpoints remain ARC/ME.\nACTH/CLIP objects below the acquisition-specific DAPI nuclear-area 5th percentile are classified as likely non-cell fragments and excluded; raw Cellpose masks are unchanged.\n",encoding="utf-8")
    print(f"Quantified {unit_df['analysis_unit_id'].nunique()} acquisition(s). Results: {result_dir}")
    return 0


def self_test() -> int:
    dapi=np.zeros((32,32),dtype=np.int32); dapi[3:10,3:10]=1; dapi[3:10,18:25]=2; dapi[18:25,3:10]=3; dapi[18:25,18:25]=4
    cfos=np.zeros_like(dapi); cfos[4:9,4:9]=1; cfos[19:24,19:24]=2
    acth_clip=np.zeros_like(dapi); acth_clip[2:11,17:26]=1; acth_clip[23:27,19:24]=2
    cset,crows=assign_cfos_rois(dapi,cfos,0.40,5); pset,prows=assign_acth_clip_rois(dapi,acth_clip,3,2)
    assert cset=={1,4},cset; assert pset=={2},pset
    assert prows[0]["size_filter_pass"] and not prows[1]["size_filter_pass"],prows
    rows=[]
    for lab in range(1,5): rows.append({"region":"ARC","phenotype":phenotype(lab in cset,lab in pset)})
    summary=count_summary(pd.DataFrame(rows),"ARC",True,100,1.0,1.0)
    test_regions=rasterize_polygons([{"region":"VMN","points":[[2,2],[10,2],[6,10]]}],(16,16))
    assert int(test_regions.max())==REGION_CODES["VMN"]
    assert summary["invariants_ok"] and summary["n_dapi"]==4
    print(json.dumps({"self_test":"PASS","cfos_positive":sorted(cset),"acth_clip_positive":sorted(pset),"summary":{k:summary[k] for k in COUNT_COLUMNS}},indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser=argparse.ArgumentParser(description="July 2026 DAPI-only HIL ARC/ME/VMN annotation and ARC/ME c-FOS and ACTH/CLIP acquisition-level quantification.",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare",action="store_true",help="Build DAPI/cFOS/ACTH-CLIP MIPs and analysis_samplesheet.csv")
    mode.add_argument("--segment",action="store_true",help="Run Cellpose on missing prepared MIPs")
    mode.add_argument("--annotate",action="store_true",help="Launch resumable DAPI-only ARC/ME/VMN polygon UI")
    mode.add_argument("--status",action="store_true",help="Print stage status as JSON")
    mode.add_argument("--audit-segmentations",action="store_true",help="Audit imported masks and ACTH/CLIP DAPI-size filtering without HIL regions")
    mode.add_argument("--quantify",action="store_true",help="Quantify annotated acquisitions and write QC/statistics/plots")
    mode.add_argument("--self-test",action="store_true",help="Run in-memory assignment/invariant test")
    parser.add_argument("--export-dir",default="zstack_tiffs_flat")
    parser.add_argument("--channel-map",default="",help="Default: <export-dir>/channel_map.csv")
    parser.add_argument("--sample-manifest",default="",help="Default: <export-dir>/sample_manifest.csv")
    parser.add_argument("--workdir",default="july_cfos_acth_clip_work")
    parser.add_argument("--include-acquisitions",default="",help="Optional comma-separated acquisition IDs")
    parser.add_argument("--include-unresolved",action="store_true",help="Explicitly permit unresolved metadata or unconfirmed channel-map rows for QC-only preparation")
    parser.add_argument("--include-nonmain-cohorts",action="store_true",help="Explicitly include acquisitions whose manifest marks include_main_analysis=false")
    parser.add_argument("--include-registration-required",action="store_true",help="Explicitly permit rows flagged registration_required; shape checks still apply")
    parser.add_argument("--include-duplicate-candidates",action="store_true",help="Explicitly permit DAPI CV duplicate candidates")
    parser.add_argument("--duplicate-correlation-threshold",type=float,default=0.995,help="DAPI CV duplicate gate threshold; acquisitions with duplicate correlation >= this value are flagged")
    parser.add_argument("--align-center-crop-registration",action="store_true",help="For explicitly allowed same-stage split-channel rows, center-crop MIPs to common shape and mark QC pending")
    parser.add_argument("--overwrite-prepare",action="store_true")
    parser.add_argument("--model-dapi",default=""); parser.add_argument("--model-cfos",default=""); parser.add_argument("--model-acth-clip",default="")
    parser.add_argument("--gpu",action="store_true"); parser.add_argument("--overwrite-seg",action="store_true")
    parser.add_argument("--diameter-dapi",type=float,default=None); parser.add_argument("--diameter-cfos",type=float,default=None); parser.add_argument("--diameter-acth-clip",type=float,default=None)
    parser.add_argument("--flow-threshold",type=float,default=None); parser.add_argument("--cellprob-threshold",type=float,default=None); parser.add_argument("--min-size-seg",type=int,default=None)
    parser.add_argument("--host",default="127.0.0.1"); parser.add_argument("--port",type=int,default=8050); parser.add_argument("--web-max-dim",type=int,default=2200)
    parser.add_argument("--cfos-dapi-overlap",type=float,default=0.40); parser.add_argument("--cfos-min-area-px",type=int,default=8)
    parser.add_argument("--acth-clip-expand-um",type=float,default=7.5); parser.add_argument("--acth-clip-min-pixels",type=int,default=5)
    parser.add_argument("--acth-clip-dapi-area-quantile",type=float,default=ACTH_CLIP_DAPI_AREA_QUANTILE,help="Reject ACTH/CLIP ROIs smaller than this acquisition-specific DAPI-area quantile")
    parser.add_argument("--region-nucleus-overlap",type=float,default=0.25)
    parser.add_argument("--qc-format",choices=["png","pdf","both"],default="both"); parser.add_argument("--qc-dpi",type=int,default=180); parser.add_argument("--scalebar-um",type=float,default=200.0)
    parser.add_argument("--endpoints",default="double_over_total_acth_clip,total_cfos_over_dapi,total_acth_clip_over_dapi,double_over_dapi,double_per_100k_um2,n_double")
    parser.add_argument("--iqr-k",type=float,default=1.5); parser.add_argument("--plot-dpi",type=int,default=600)
    return parser


def main() -> int:
    args=build_parser().parse_args()
    if args.self_test: return self_test()
    if args.prepare: return prepare(args)
    if args.segment: return segment(args)
    if args.annotate: return annotate(args)
    if args.status: return status(args)
    if args.audit_segmentations: return audit_segmentations(args)
    if args.quantify: return quantify(args)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
