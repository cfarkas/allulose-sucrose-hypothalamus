#!/usr/bin/env python3
"""Regenerate Figure 5 Iba1 instances from the reviewed DAPI masks.

This maintenance utility is raw-safe. It reads immutable microscopy TIFFs and
versioned DAPI masks, applies the same seeded-watershed algorithm used by the
Figure 5 analyzer, and writes relocatable NPY/TIFF/QC/JSON artifacts under a
separate HIL root. The normal reproduction workflow consumes these reviewed
outputs; it never rewrites raw_data sidecars.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import inspect
import json
import logging
import os
import platform
import sys
import time
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import tifffile
from PIL import Image
from skimage.segmentation import find_boundaries


HERE = Path(__file__).resolve()
PAPER_ROOT = next(
    (
        path
        for path in HERE.parents
        if path.name == "Paper"
        or ((path / "scripts" / "setup").is_dir() and (path / "README.txt").is_file())
    ),
    None,
)
if PAPER_ROOT is None:
    raise RuntimeError(f"Could not locate the Paper repository above {HERE}")

DEFAULT_RAW_ROOT = PAPER_ROOT / "Fig5" / "raw_data"
DEFAULT_DAPI_ROOT = PAPER_ROOT / "Fig5" / "hil_review" / "dapi_final_20260830_v1"
DEFAULT_OUTPUT_ROOT = PAPER_ROOT / "Fig5" / "hil_review" / "iba1_final_20260830_v1"
DEFAULT_ANALYZER = PAPER_ROOT / "Fig5" / "01_analyze_gfap_iba1_microglia.py"
VERSION = "fig5_dapi_seeded_iba1_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Regenerate all Figure 5 Iba1 masks from reviewed DAPI labels.",
    )
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--dapi-root", type=Path, default=DEFAULT_DAPI_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--analyzer", type=Path, default=DEFAULT_ANALYZER)
    parser.add_argument("--expected-sections", type=int, default=61)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_npy(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.npy")
    with temporary.open("wb") as handle:
        np.save(handle, payload, allow_pickle=True)
    os.replace(temporary, path)


def atomic_tiff(path: Path, labels: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.tif")
    tifffile.imwrite(temporary, labels, compression="deflate", metadata=None)
    os.replace(temporary, path)


def atomic_png(path: Path, image: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.png")
    Image.fromarray(image).save(temporary, format="PNG", optimize=True)
    os.replace(temporary, path)


def load_labels(path: Path) -> np.ndarray:
    payload = np.load(path, allow_pickle=True)
    if isinstance(payload, np.ndarray) and payload.dtype == object and payload.shape == ():
        payload = payload.item()
    if isinstance(payload, dict):
        payload = payload.get("masks")
    labels = np.asarray(payload)
    if labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer):
        raise ValueError(f"Invalid reviewed DAPI labels at {path}: {labels.shape}, {labels.dtype}")
    return labels.astype(np.int32, copy=False)


def load_analyzer(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("fig5_analyzer_for_iba1", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import analyzer: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.import_science_stack()
    return module


def algorithm_args() -> SimpleNamespace:
    return SimpleNamespace(
        intensity_percentile_low=0.5,
        intensity_percentile_high=99.8,
        um_per_px=1.51,
        iba1_background_radius_um=22.0,
        iba1_threshold_factor=0.65,
        iba1_min_cell_area_um2=12.0,
        iba1_nucleus_expand_um=7.5,
        iba1_nucleus_min_score_quantile=0.60,
        iba1_max_cell_radius_um=55.0,
        iba1_max_cell_area_um2=3500.0,
        iba1_min_foreground_fraction=0.02,
    )


def parameter_dict(args: SimpleNamespace) -> dict[str, Any]:
    return {key: value for key, value in vars(args).items()}


def qc_overlay(signal: np.ndarray, dapi: np.ndarray, iba1: np.ndarray) -> np.ndarray:
    finite = signal[np.isfinite(signal)]
    lo, hi = np.percentile(finite, [0.5, 99.8]) if finite.size else (0.0, 1.0)
    if hi <= lo:
        hi = lo + 1.0
    gray = np.clip((signal.astype(np.float32) - lo) / (hi - lo), 0.0, 1.0)
    gray = np.round(gray * 255).astype(np.uint8)
    rgb = np.repeat(gray[..., None], 3, axis=2)
    rgb[find_boundaries(dapi, mode="outer")] = np.array([32, 255, 255], dtype=np.uint8)
    rgb[find_boundaries(iba1, mode="outer")] = np.array([255, 32, 32], dtype=np.uint8)
    return rgb


def mask_stats(labels: np.ndarray) -> dict[str, Any]:
    counts = np.bincount(labels.ravel())[1:]
    counts = counts[counts > 0]
    return {
        "instances": int(counts.size),
        "positive_pixels": int(np.count_nonzero(labels)),
        "area_px_q10": float(np.percentile(counts, 10)) if counts.size else None,
        "area_px_median": float(np.median(counts)) if counts.size else None,
        "area_px_q90": float(np.percentile(counts, 90)) if counts.size else None,
        "area_px_max": int(counts.max()) if counts.size else None,
    }


def receipt_matches(
    path: Path,
    source_dapi_sha: str,
    source_iba1_sha: str,
    dapi_mask_sha: str,
    algorithm_sha: str,
    parameters: dict[str, Any],
    expected_outputs: dict[str, Path],
) -> bool:
    try:
        receipt = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    if not (
        receipt.get("status") == "PASS"
        and receipt.get("raw_data_modified") is False
        and receipt.get("source_dapi_sha256") == source_dapi_sha
        and receipt.get("source_iba1_sha256") == source_iba1_sha
        and receipt.get("reviewed_dapi_segmentation_sha256") == dapi_mask_sha
        and receipt.get("algorithm_sha256") == algorithm_sha
        and receipt.get("parameters") == parameters
    ):
        return False
    hashes = receipt.get("output_hashes")
    if not isinstance(hashes, dict):
        return False
    return all(
        output.is_file() and hashes.get(f"{name}_sha256") == sha256(output)
        for name, output in expected_outputs.items()
    )


def write_summary(output_root: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = sorted({key for row in rows for key in row})
    destination = output_root / "iba1_regeneration_summary.csv"
    temporary = destination.with_suffix(".tmp.csv")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, destination)


def main() -> int:
    cli = parse_args()
    raw_root = cli.raw_root.expanduser().resolve()
    dapi_root = cli.dapi_root.expanduser().resolve()
    output_root = cli.output_root.expanduser().resolve()
    analyzer_path = cli.analyzer.expanduser().resolve()
    if not raw_root.is_dir() or not dapi_root.is_dir() or not analyzer_path.is_file():
        raise FileNotFoundError(
            f"Required input missing: raw={raw_root}, DAPI={dapi_root}, analyzer={analyzer_path}"
        )
    if output_root == raw_root or raw_root in output_root.parents or output_root in raw_root.parents:
        raise RuntimeError("--output-root must be separate from Fig5/raw_data")

    images = sorted(raw_root.glob("*/S*/Image_*_Iba1.tif"))
    images += sorted(raw_root.glob("*/S*/Image_*_Iba1.tiff"))
    images = sorted(set(images))
    if cli.limit > 0:
        images = images[: cli.limit]
    elif len(images) != cli.expected_sections:
        raise RuntimeError(f"Expected {cli.expected_sections} Iba1 sections, found {len(images)}")

    analyzer = load_analyzer(analyzer_path)
    algo_args = algorithm_args()
    parameters = parameter_dict(algo_args)
    algorithm_source = inspect.getsource(analyzer._gi_segment_iba1_seeded).encode("utf-8")
    algorithm_sha = hashlib.sha256(algorithm_source).hexdigest()
    logger = logging.getLogger("fig5_iba1_regeneration")
    logger.handlers.clear()
    logger.addHandler(logging.StreamHandler(sys.stdout))
    logger.setLevel(logging.WARNING)

    output_root.mkdir(parents=True, exist_ok=True)
    started_at = dt.datetime.now(dt.timezone.utc)
    metadata: dict[str, Any] = {
        "version": VERSION,
        "started_at": started_at.isoformat(),
        "raw_root": str(raw_root),
        "dapi_root": str(dapi_root),
        "output_root": str(output_root),
        "analyzer_path": str(analyzer_path),
        "algorithm_sha256": algorithm_sha,
        "parameters": parameters,
        "python": sys.version,
        "platform": platform.platform(),
        "sections": len(images),
        "raw_data_modified": False,
    }
    atomic_json(output_root / "RUN_METADATA.json", metadata)

    rows: list[dict[str, Any]] = []
    for index, iba1_path in enumerate(images, start=1):
        section_started = time.monotonic()
        section_dir = iba1_path.parent
        animal = section_dir.parent.name
        section = section_dir.name
        dapi_image = iba1_path.with_name(iba1_path.name.replace("_Iba1.", "_DAPI."))
        dapi_mask = (
            dapi_root / animal / section / f"{animal}_{section}_DAPI_seg.npy"
        )
        if not dapi_image.is_file() or not dapi_mask.is_file():
            raise FileNotFoundError(
                f"Missing DAPI source or reviewed mask for {animal} {section}: "
                f"{dapi_image}, {dapi_mask}"
            )

        prefix = f"{animal}_{section}_Iba1"
        section_out = output_root / animal / section
        mask_path = section_out / f"{prefix}_seg.npy"
        label_path = section_out / f"{prefix}_labels.tif"
        qc_path = section_out / f"{prefix}_QC.png"
        receipt_path = section_out / f"{prefix}_receipt.json"
        expected_outputs = {
            "segmentation_npy": mask_path,
            "label_tiff": label_path,
            "qc_overlay_png": qc_path,
        }
        source_dapi_sha = sha256(dapi_image)
        source_iba1_sha = sha256(iba1_path)
        dapi_mask_sha = sha256(dapi_mask)

        if not cli.force and receipt_matches(
            receipt_path,
            source_dapi_sha,
            source_iba1_sha,
            dapi_mask_sha,
            algorithm_sha,
            parameters,
            expected_outputs,
        ):
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            row = dict(receipt["summary_row"])
            row["action"] = "verified_resume_skip"
            rows.append(row)
            write_summary(output_root, rows)
            print(f"[{index:02d}/{len(images):02d}] {animal} {section}: verified, skipped", flush=True)
            continue

        dapi_labels = load_labels(dapi_mask)
        raw_iba1 = tifffile.imread(iba1_path)
        iba1_signal = np.asarray(analyzer._gi_signal_from_rgb(raw_iba1, "Iba1"), dtype=np.float32)
        if iba1_signal.shape != dapi_labels.shape:
            raise RuntimeError(
                f"Iba1/DAPI shape mismatch for {animal} {section}: "
                f"{iba1_signal.shape} != {dapi_labels.shape}"
            )
        labels, foreground, _corrected, info = analyzer._gi_segment_iba1_seeded(
            dapi_labels,
            iba1_signal,
            algo_args,
            logger,
        )
        labels = np.asarray(labels, dtype=np.int32)
        stats = mask_stats(labels)
        status = "PASS" if labels.shape == dapi_labels.shape and stats["instances"] > 0 else "REVIEW"

        payload = {
            "masks": labels,
            "source": "dapi_seeded_iba1_watershed",
            "version": VERSION,
            "source_dapi_path": str(dapi_image),
            "source_dapi_sha256": source_dapi_sha,
            "source_iba1_path": str(iba1_path),
            "source_iba1_sha256": source_iba1_sha,
            "reviewed_dapi_segmentation_path": str(dapi_mask),
            "reviewed_dapi_segmentation_sha256": dapi_mask_sha,
            "algorithm_sha256": algorithm_sha,
            "parameters": parameters,
            "segmentation_info": info,
            "qc_status": status,
        }
        atomic_npy(mask_path, payload)
        label_dtype = np.uint16 if int(labels.max()) <= np.iinfo(np.uint16).max else np.uint32
        atomic_tiff(label_path, labels.astype(label_dtype))
        atomic_png(qc_path, qc_overlay(iba1_signal, dapi_labels, labels))
        output_hashes = {
            f"{name}_sha256": sha256(output)
            for name, output in expected_outputs.items()
        }
        elapsed = time.monotonic() - section_started
        summary_row = {
            "animal": animal,
            "section": section,
            "status": status,
            "action": "segmented",
            "instances": stats["instances"],
            "positive_pixels": stats["positive_pixels"],
            "area_px_median": stats["area_px_median"],
            "candidate_nuclei": info.get("iba1_candidate_nuclei"),
            "foreground_fraction": info.get("iba1_foreground_fraction"),
            "elapsed_seconds": round(elapsed, 3),
            "source_iba1_sha256": source_iba1_sha,
            "reviewed_dapi_segmentation_sha256": dapi_mask_sha,
            "label_tiff_sha256": output_hashes["label_tiff_sha256"],
        }
        receipt = {
            "status": status,
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "raw_data_modified": False,
            "source_dapi_path": str(dapi_image),
            "source_dapi_sha256": source_dapi_sha,
            "source_iba1_path": str(iba1_path),
            "source_iba1_sha256": source_iba1_sha,
            "reviewed_dapi_segmentation_path": str(dapi_mask),
            "reviewed_dapi_segmentation_sha256": dapi_mask_sha,
            "algorithm_sha256": algorithm_sha,
            "parameters": parameters,
            "mask_statistics": stats,
            "segmentation_info": info,
            "foreground_fraction": float(np.mean(foreground)),
            "outputs": {name: str(output) for name, output in expected_outputs.items()},
            "output_hashes": output_hashes,
            "summary_row": summary_row,
        }
        atomic_json(receipt_path, receipt)
        rows.append(summary_row)
        write_summary(output_root, rows)
        remaining = len(images) - index
        mean_elapsed = float(np.mean([
            float(row["elapsed_seconds"]) for row in rows if row.get("action") == "segmented"
        ]))
        eta_minutes = remaining * mean_elapsed / 60.0
        print(
            f"[{index:02d}/{len(images):02d}] {animal} {section}: {status}; "
            f"n={stats['instances']}; {elapsed:.1f}s; ETA={eta_minutes:.1f}m",
            flush=True,
        )

    metadata.update({
        "finished_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pass_sections": sum(row["status"] == "PASS" for row in rows),
        "review_sections": sum(row["status"] == "REVIEW" for row in rows),
        "summary_csv": str(output_root / "iba1_regeneration_summary.csv"),
    })
    atomic_json(output_root / "RUN_METADATA.json", metadata)
    print(json.dumps(
        {key: metadata[key] for key in ("sections", "pass_sections", "review_sections", "summary_csv")},
        indent=2,
    ))
    return 0 if metadata["review_sections"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
