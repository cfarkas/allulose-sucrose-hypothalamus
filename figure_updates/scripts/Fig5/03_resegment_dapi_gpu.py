#!/usr/bin/env python3
"""Versioned, raw-safe GPU DAPI re-segmentation for all Figure 5 sections.

The source TIFFs and their historical ``*_seg.npy`` siblings are immutable.
This utility reads them, applies the shipped ``dapi_channel`` Cellpose model,
and writes masks, label TIFFs, QC overlays, and receipts below a separate HIL
review root.  Runs are resumable: a section is skipped only when its receipt
matches the current source image, model, and parameters and all outputs exist.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import tifffile
from PIL import Image
from scipy import ndimage as ndi
from skimage.measure import regionprops_table
from skimage.segmentation import find_boundaries, relabel_sequential


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
DEFAULT_OUTPUT_ROOT = PAPER_ROOT / "Fig5" / "hil_review" / "dapi_final_20260830_v1"
DEFAULT_REGION_ROOT = PAPER_ROOT / "Fig5" / "hil_review" / "human_final_20260823_v1"
DEFAULT_MODEL = (
    PAPER_ROOT
    / "Fig2"
    / "raw_data"
    / "Apotome"
    / "24_06_2025"
    / "analysis"
    / "models"
    / "dapi_channel"
)
VERSION = "fig5_dapi_channel_gpu_v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Re-segment Figure 5 DAPI images with dapi_channel on GPU without altering raw_data.",
    )
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--region-root", type=Path, default=DEFAULT_REGION_ROOT)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--diameter", type=float, default=6.5)
    parser.add_argument("--flow-threshold", type=float, default=0.4)
    parser.add_argument("--cellprob-threshold", type=float, default=0.0)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-size", type=int, default=350)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--gpu-device", type=int, default=0)
    parser.add_argument("--expected-sections", type=int, default=61)
    parser.add_argument("--limit", type=int, default=0, help="Debug: process only the first N sections.")
    parser.add_argument("--force", action="store_true", help="Recompute sections with matching receipts.")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def atomic_npy(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npy")
    with tmp.open("wb") as handle:
        np.save(handle, payload, allow_pickle=True)
    os.replace(tmp, path)


def atomic_tiff(path: Path, labels: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.tif")
    tifffile.imwrite(tmp, labels, compression="deflate", metadata=None)
    os.replace(tmp, path)


def atomic_png(path: Path, image: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.png")
    Image.fromarray(image).save(tmp, format="PNG", optimize=True)
    os.replace(tmp, path)


def source_signal(image: np.ndarray) -> np.ndarray:
    arr = np.asarray(image)
    while arr.ndim > 3:
        arr = arr.max(axis=0)
    if arr.ndim == 2:
        return arr.astype(np.uint8, copy=False)
    if arr.ndim != 3 or arr.shape[-1] < 3:
        raise ValueError(f"Unexpected DAPI image shape: {arr.shape}")
    return arr[..., 2].astype(np.uint8, copy=False)


def cyan_input(signal: np.ndarray) -> np.ndarray:
    # dapi_channel was trained on Keyence cyan exports: R=0 and G=B=DAPI.
    cyan = np.zeros((*signal.shape, 3), dtype=np.uint8)
    cyan[..., 1] = signal
    cyan[..., 2] = signal
    return cyan


def filter_and_relabel(labels: np.ndarray, min_size: int, max_size: int) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.int32)
    counts = np.bincount(labels.ravel())
    keep = np.ones(counts.size, dtype=bool)
    keep[0] = False
    keep &= counts >= int(min_size)
    if int(max_size) > 0:
        keep &= counts <= int(max_size)
    filtered = labels.copy()
    filtered[~keep[filtered]] = 0
    filtered, _, _ = relabel_sequential(filtered)
    return np.asarray(filtered, dtype=np.int32)


def mask_stats(labels: np.ndarray, signal: np.ndarray) -> dict[str, Any]:
    counts = np.bincount(labels.ravel())[1:]
    counts = counts[counts > 0]
    inside = signal[labels > 0].astype(np.float32)
    outside = signal[labels == 0].astype(np.float32)
    return {
        "instances": int(counts.size),
        "positive_pixels": int(np.count_nonzero(labels)),
        "area_px_q10": float(np.percentile(counts, 10)) if counts.size else None,
        "area_px_median": float(np.median(counts)) if counts.size else None,
        "area_px_q90": float(np.percentile(counts, 90)) if counts.size else None,
        "area_px_max": int(counts.max()) if counts.size else None,
        "raw_intensity_inside_median": float(np.median(inside)) if inside.size else None,
        "raw_intensity_outside_median": float(np.median(outside)) if outside.size else None,
        "raw_intensity_inside_mean": float(np.mean(inside)) if inside.size else None,
        "raw_intensity_outside_mean": float(np.mean(outside)) if outside.size else None,
    }


def region_qc(labels: np.ndarray, signal: np.ndarray, region_path: Path) -> dict[str, Any]:
    result: dict[str, Any] = {
        "region_mask_path": str(region_path),
        "region_mask_sha256": sha256(region_path) if region_path.is_file() else None,
        "largest_raw_bright_label_empty_area_px": None,
        "largest_dropout_region": None,
        "largest_dropout_bbox_y0_x0_y1_x1": None,
        "left_vmn_instances": None,
        "right_vmn_instances": None,
        "left_vmn_positive_fraction": None,
        "right_vmn_positive_fraction": None,
    }
    if not region_path.is_file():
        return result
    regions = tifffile.imread(region_path)
    if regions.shape != labels.shape:
        raise ValueError(f"Region/DAPI shape mismatch: {regions.shape} != {labels.shape}")

    density = ndi.gaussian_filter((labels > 0).astype(np.float32), sigma=18.0)
    largest: dict[str, Any] | None = None
    for code, name in ((1, "ARC"), (2, "ME"), (3, "VMN")):
        components, n_components = ndi.label(regions == code)
        for component_id in range(1, n_components + 1):
            component = components == component_id
            if int(component.sum()) < 3000:
                continue
            interior = ndi.binary_erosion(component, iterations=8)
            holes, n_holes = ndi.label(interior & (density < 0.02))
            component_mean = float(np.mean(signal[component]))
            for hole_id in range(1, n_holes + 1):
                hole = holes == hole_id
                area = int(hole.sum())
                if area < 2500:
                    continue
                hole_mean = float(np.mean(signal[hole]))
                raw_ratio = hole_mean / component_mean if component_mean > 0 else 0.0
                # Only retain holes with raw signal comparable to their anatomical component.
                if raw_ratio < 0.6:
                    continue
                ys, xs = np.nonzero(hole)
                candidate = {
                    "area": area,
                    "region": name,
                    "bbox": f"{ys.min()}:{xs.min()}:{ys.max()+1}:{xs.max()+1}",
                    "raw_ratio": raw_ratio,
                }
                if largest is None or candidate["area"] > largest["area"]:
                    largest = candidate
    if largest is not None:
        result["largest_raw_bright_label_empty_area_px"] = largest["area"]
        result["largest_dropout_region"] = largest["region"]
        result["largest_dropout_bbox_y0_x0_y1_x1"] = largest["bbox"]
        result["largest_dropout_raw_intensity_ratio"] = largest["raw_ratio"]

    vmn, n_vmn = ndi.label(regions == 3)
    components = []
    for component_id in range(1, n_vmn + 1):
        yy, xx = np.nonzero(vmn == component_id)
        if yy.size >= 3000:
            components.append((component_id, float(xx.mean())))
    if len(components) >= 2:
        left_id = min(components, key=lambda item: item[1])[0]
        right_id = max(components, key=lambda item: item[1])[0]
        left = vmn == left_id
        right = vmn == right_id
        props = regionprops_table(labels, properties=("centroid",))
        cy = np.rint(props["centroid-0"]).astype(int).clip(0, labels.shape[0] - 1)
        cx = np.rint(props["centroid-1"]).astype(int).clip(0, labels.shape[1] - 1)
        result.update(
            {
                "left_vmn_instances": int(np.sum(left[cy, cx])),
                "right_vmn_instances": int(np.sum(right[cy, cx])),
                "left_vmn_positive_fraction": float(np.mean((labels > 0)[left])),
                "right_vmn_positive_fraction": float(np.mean((labels > 0)[right])),
            }
        )
    return result


def qc_overlay(signal: np.ndarray, labels: np.ndarray, regions: np.ndarray | None) -> np.ndarray:
    lo, hi = np.percentile(signal, [1.0, 99.5])
    if hi <= lo:
        hi = lo + 1.0
    gray = np.clip((signal.astype(np.float32) - lo) / (hi - lo), 0, 1)
    gray = np.round(gray * 255).astype(np.uint8)
    rgb = np.repeat(gray[..., None], 3, axis=2)
    borders = find_boundaries(labels, mode="outer")
    rgb[borders] = np.array([255, 32, 32], dtype=np.uint8)
    if regions is not None:
        region_border = find_boundaries(regions, mode="outer")
        rgb[region_border] = np.array([32, 255, 255], dtype=np.uint8)
    return rgb


def receipt_matches(receipt_path: Path, source_sha: str, model_sha: str, params: dict[str, Any]) -> bool:
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return (
        receipt.get("source_dapi_sha256") == source_sha
        and receipt.get("model_sha256") == model_sha
        and receipt.get("parameters") == params
        and receipt.get("status") in {"PASS", "REVIEW"}
        and all(Path(path).is_file() for path in receipt.get("outputs", {}).values())
    )


def write_summary(output_root: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fields = sorted({key for row in rows for key in row})
    path = output_root / "dapi_resegmentation_summary.csv"
    tmp = path.with_suffix(".tmp.csv")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()
    region_root = args.region_root.expanduser().resolve()
    model_path = args.model.expanduser().resolve()
    if not raw_root.is_dir():
        raise FileNotFoundError(raw_root)
    if not model_path.is_file():
        raise FileNotFoundError(model_path)
    if output_root == raw_root or raw_root in output_root.parents or output_root in raw_root.parents:
        raise RuntimeError("--output-root must be separate from Fig5/raw_data")

    images = sorted(raw_root.glob("*/S*/Image_*_DAPI.tif"))
    images += sorted(raw_root.glob("*/S*/Image_*_DAPI.tiff"))
    images = sorted(set(images))
    if args.limit > 0:
        images = images[: args.limit]
    elif len(images) != args.expected_sections:
        raise RuntimeError(f"Expected {args.expected_sections} DAPI sections, found {len(images)}")

    import torch
    from cellpose import models

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for this batch; activate the working cellpose_env")
    torch.cuda.set_device(args.gpu_device)
    output_root.mkdir(parents=True, exist_ok=True)
    model_sha = sha256(model_path)
    params: dict[str, Any] = {
        "diameter_px": float(args.diameter),
        "flow_threshold": float(args.flow_threshold),
        "cellprob_threshold": float(args.cellprob_threshold),
        "min_size_px": int(args.min_size),
        "max_size_px": int(args.max_size),
        "normalization_percentiles": [1.0, 99.0],
        "input_conversion": "Keyence blue DAPI plane -> cyan RGB (R=0,G=DAPI,B=DAPI)",
        "resample": True,
    }
    run_started = dt.datetime.now(dt.timezone.utc)
    run_metadata = {
        "version": VERSION,
        "started_at": run_started.isoformat(),
        "raw_root": str(raw_root),
        "output_root": str(output_root),
        "region_root": str(region_root),
        "model_path": str(model_path),
        "model_sha256": model_sha,
        "parameters": params,
        "python": sys.version,
        "platform": platform.platform(),
        "cellpose": importlib.metadata.version("cellpose"),
        "torch": torch.__version__,
        "cuda": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(args.gpu_device),
        "sections": len(images),
        "raw_data_modified": False,
    }
    atomic_json(output_root / "RUN_METADATA.json", run_metadata)

    print(f"Loading dapi_channel on {run_metadata['gpu']} (CUDA {run_metadata['cuda']})", flush=True)
    model = models.CellposeModel(gpu=True, pretrained_model=str(model_path))
    summary_rows: list[dict[str, Any]] = []
    for index, image_path in enumerate(images, start=1):
        section_started = time.monotonic()
        section_dir = image_path.parent
        animal = section_dir.parent.name
        section = section_dir.name
        prefix = f"{animal}_{section}_DAPI"
        section_out = output_root / animal / section
        labels_npy = section_out / f"{prefix}_seg.npy"
        labels_tif = section_out / f"{prefix}_labels.tif"
        overlay_png = section_out / f"{prefix}_QC.png"
        receipt_path = section_out / f"{prefix}_receipt.json"
        source_sha = sha256(image_path)
        historical_seg = image_path.with_name(image_path.stem + "_seg.npy")
        historical_sha = sha256(historical_seg) if historical_seg.is_file() else None

        if not args.force and receipt_matches(receipt_path, source_sha, model_sha, params):
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            row = dict(receipt["summary_row"])
            row["action"] = "verified_resume_skip"
            summary_rows.append(row)
            write_summary(output_root, summary_rows)
            print(f"[{index:02d}/{len(images):02d}] {animal} {section}: verified, skipped", flush=True)
            continue

        raw = tifffile.imread(image_path)
        signal = source_signal(raw)
        masks, _flows, _styles = model.eval(
            cyan_input(signal),
            batch_size=int(args.batch_size),
            resample=True,
            normalize=True,
            diameter=float(args.diameter),
            flow_threshold=float(args.flow_threshold),
            cellprob_threshold=float(args.cellprob_threshold),
            min_size=int(args.min_size),
        )
        labels = filter_and_relabel(masks, int(args.min_size), int(args.max_size))
        if labels.shape != signal.shape or int(labels.max()) == 0:
            raise RuntimeError(f"Invalid segmentation for {animal} {section}: shape={labels.shape}, max={labels.max()}")

        region_path = region_root / animal / section / f"{animal}_{section}_ARC_ME_labels.tif"
        regions = tifffile.imread(region_path) if region_path.is_file() else None
        stats = mask_stats(labels, signal)
        anatomical_qc = region_qc(labels, signal, region_path)
        dropout_area = anatomical_qc.get("largest_raw_bright_label_empty_area_px") or 0
        status = "PASS" if int(dropout_area) < 5000 else "REVIEW"

        payload = {
            "masks": labels.astype(np.int32),
            "source": "cellpose_dapi_channel",
            "version": VERSION,
            "source_dapi_path": str(image_path),
            "source_dapi_sha256": source_sha,
            "historical_segmentation_path": str(historical_seg),
            "historical_segmentation_sha256": historical_sha,
            "model_path": str(model_path),
            "model_sha256": model_sha,
            "parameters": params,
            "qc_status": status,
        }
        atomic_npy(labels_npy, payload)
        out_dtype = np.uint16 if int(labels.max()) <= np.iinfo(np.uint16).max else np.uint32
        atomic_tiff(labels_tif, labels.astype(out_dtype))
        atomic_png(overlay_png, qc_overlay(signal, labels, regions))
        elapsed = time.monotonic() - section_started
        outputs = {
            "segmentation_npy": str(labels_npy),
            "label_tiff": str(labels_tif),
            "qc_overlay_png": str(overlay_png),
        }
        output_sha = {key + "_sha256": sha256(Path(path)) for key, path in outputs.items()}
        summary_row = {
            "animal": animal,
            "section": section,
            "status": status,
            "action": "segmented",
            "instances": stats["instances"],
            "positive_pixels": stats["positive_pixels"],
            "area_px_median": stats["area_px_median"],
            "largest_dropout_area_px": dropout_area,
            "left_vmn_instances": anatomical_qc.get("left_vmn_instances"),
            "right_vmn_instances": anatomical_qc.get("right_vmn_instances"),
            "elapsed_seconds": round(elapsed, 3),
            "source_dapi_sha256": source_sha,
            "label_tiff_sha256": output_sha["label_tiff_sha256"],
        }
        receipt = {
            "status": status,
            "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "raw_data_modified": False,
            "source_dapi_path": str(image_path),
            "source_dapi_sha256": source_sha,
            "historical_segmentation_path": str(historical_seg),
            "historical_segmentation_sha256": historical_sha,
            "model_path": str(model_path),
            "model_sha256": model_sha,
            "parameters": params,
            "mask_statistics": stats,
            "anatomical_qc": anatomical_qc,
            "outputs": outputs,
            "output_hashes": output_sha,
            "summary_row": summary_row,
        }
        atomic_json(receipt_path, receipt)
        summary_rows.append(summary_row)
        write_summary(output_root, summary_rows)
        remaining = len(images) - index
        mean_elapsed = float(np.mean([float(row["elapsed_seconds"]) for row in summary_rows if row.get("action") == "segmented"]))
        eta_min = remaining * mean_elapsed / 60.0
        print(
            f"[{index:02d}/{len(images):02d}] {animal} {section}: {status}; "
            f"n={stats['instances']}; dropout={dropout_area}px; {elapsed:.1f}s; ETA={eta_min:.1f}m",
            flush=True,
        )

    run_metadata.update(
        {
            "finished_at": dt.datetime.now(dt.timezone.utc).isoformat(),
            "pass_sections": sum(row["status"] == "PASS" for row in summary_rows),
            "review_sections": sum(row["status"] == "REVIEW" for row in summary_rows),
            "summary_csv": str(output_root / "dapi_resegmentation_summary.csv"),
        }
    )
    atomic_json(output_root / "RUN_METADATA.json", run_metadata)
    print(json.dumps({key: run_metadata[key] for key in ("sections", "pass_sections", "review_sections", "summary_csv")}, indent=2))
    return 0 if run_metadata["review_sections"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
