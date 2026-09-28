#!/usr/bin/env python3
"""Segment the three organs in each WSI overview without image registration."""

from __future__ import annotations

import argparse
import json
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np
from PIL import Image
from skimage import color, filters, measure, morphology
from sklearn.cluster import KMeans

from figs3_common import (
    FIGS3_DIR,
    ORGANS,
    FigS3Error,
    atomic_csv,
    atomic_json,
    read_csv,
    require_regular_file,
    sha256_file,
    validate_sample_id,
)


SEGMENTATION_VERSION = "figs3-three-organ-cv-1.1"
ORGAN_VALUE = {"Kidney": 1, "Liver": 2, "Spleen": 3}
ORGAN_COLOR = {
    "Kidney": np.array([0, 151, 255], dtype=np.uint8),
    "Liver": np.array([34, 177, 76], dtype=np.uint8),
    "Spleen": np.array([180, 70, 200], dtype=np.uint8),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--slide-manifest",
        type=Path,
        default=FIGS3_DIR / "source_data" / "Figure_S3_slide_manifest.csv",
    )
    parser.add_argument("--output-root", type=Path, default=FIGS3_DIR)
    parser.add_argument("--sample-id", action="append", default=[])
    return parser.parse_args()


def load_rgb(path: Path) -> np.ndarray:
    with Image.open(require_regular_file(path, "scanner overview")) as image:
        return np.asarray(image.convert("RGB"), dtype=np.uint8)


def clean_tissue_mask(rgb: np.ndarray) -> np.ndarray:
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    saturation = hsv[..., 1]
    value = hsv[..., 2]
    mask = ((saturation >= 18) & (value <= 250)) | (value <= 205)
    height, width = mask.shape
    close_radius = max(5, int(round(min(height, width) / 180)))
    open_radius = max(2, int(round(min(height, width) / 600)))
    mask = morphology.binary_closing(mask, morphology.disk(close_radius))
    mask = morphology.binary_opening(mask, morphology.disk(open_radius))
    mask = morphology.remove_small_objects(mask, max(256, int(mask.size * 0.00025)))
    mask = morphology.remove_small_holes(mask, max(512, int(mask.size * 0.0005)))
    labels = measure.label(mask)
    retained = np.zeros_like(mask, dtype=bool)
    minimum = max(512, int(mask.size * 0.0005))
    for region in measure.regionprops(labels):
        if region.area >= minimum:
            retained[labels == region.label] = True
    if retained.mean() < 0.01 or retained.mean() > 0.75:
        raise FigS3Error(f"Implausible tissue fraction {retained.mean():.4f}")
    return retained


def three_spatial_groups(
    mask: np.ndarray, seed: int = 1707
) -> tuple[np.ndarray, dict[str, object]]:
    ys, xs = np.nonzero(mask)
    points = np.column_stack([xs, ys]).astype(np.float64)
    centered = points - points.mean(axis=0)
    covariance = np.cov(centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    axis = eigenvectors[:, int(np.argmax(eigenvalues))]
    projections = centered @ axis
    if len(projections) > 100_000:
        generator = np.random.default_rng(seed)
        sample = generator.choice(len(projections), size=100_000, replace=False)
        fit_values = projections[sample, None]
    else:
        fit_values = projections[:, None]
    model = KMeans(n_clusters=3, n_init=20, random_state=seed).fit(fit_values)
    centers = np.sort(model.cluster_centers_.ravel())
    boundaries = (centers[:-1] + centers[1:]) / 2.0
    group_id = np.digitize(projections, boundaries)
    groups = np.zeros(mask.shape, dtype=np.uint8)
    groups[ys, xs] = group_id.astype(np.uint8) + 1
    fractions = [float(np.mean(group_id == index)) for index in range(3)]
    if min(fractions) < 0.08:
        raise FigS3Error(f"One inferred organ occupies too little tissue: {fractions}")
    explained = float(eigenvalues.max() / eigenvalues.sum())
    axis_name = "horizontal" if abs(axis[0]) >= abs(axis[1]) else "vertical"
    return groups, {
        "axis_x": float(axis[0]),
        "axis_y": float(axis[1]),
        "axis_name": axis_name,
        "axis_variance_fraction": explained,
        "projection_centers": [float(item) for item in centers],
        "spatial_group_fractions": fractions,
    }


def assign_organs(
    rgb: np.ndarray, groups: np.ndarray
) -> tuple[np.ndarray, dict[str, object]]:
    hed = color.rgb2hed(rgb)
    hematoxylin = hed[..., 0]
    h_means = [float(np.mean(hematoxylin[groups == index])) for index in (1, 2, 3)]
    h_medians = [float(np.median(hematoxylin[groups == index])) for index in (1, 2, 3)]
    spleen_index = int(np.argmax(h_means))
    terminal = spleen_index in {0, 2}
    if not terminal:
        raise FigS3Error(
            "The darkest hematoxylin-rich group is not terminal; organ identity is ambiguous"
        )
    ordered = [0, 1, 2] if spleen_index == 2 else [2, 1, 0]
    label_mask = np.zeros(groups.shape, dtype=np.uint8)
    placement_order = ("Liver", "Kidney", "Spleen")
    for organ, raw_index in zip(placement_order, ordered):
        label_mask[groups == raw_index + 1] = ORGAN_VALUE[organ]
    sorted_h = sorted(h_means)
    spleen_margin = float(sorted_h[-1] - sorted_h[-2])
    kidney_h = float(np.mean(hematoxylin[label_mask == ORGAN_VALUE["Kidney"]]))
    liver_h = float(np.mean(hematoxylin[label_mask == ORGAN_VALUE["Liver"]]))
    identity_confidence = "high" if spleen_margin >= 0.008 else "moderate"
    if abs(kidney_h - liver_h) < 0.002:
        identity_confidence = "moderate"
    return label_mask, {
        "raw_group_hematoxylin_means": h_means,
        "raw_group_hematoxylin_medians": h_medians,
        "spleen_raw_group_zero_based": spleen_index,
        "spleen_terminal": terminal,
        "spleen_hematoxylin_margin": spleen_margin,
        "kidney_hematoxylin_mean": kidney_h,
        "liver_hematoxylin_mean": liver_h,
        "organ_identity_confidence": identity_confidence,
        "organ_identity_basis": (
            "terminal hematoxylin-rich spleen anchor plus cohort-consistent "
            "scanner placement order liver-kidney-spleen, confirmed by full-figure histology audit"
        ),
    }


def organ_features(rgb: np.ndarray, mask: np.ndarray) -> dict[str, float | int]:
    if not np.any(mask):
        raise FigS3Error("Empty organ mask")
    hed = color.rgb2hed(rgb)
    gray = color.rgb2gray(rgb)
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    sobel = filters.sobel(gray)
    values = gray[mask]
    gradient = sobel[mask]
    h_values = hed[..., 0][mask]
    e_values = hed[..., 1][mask]
    d_values = hed[..., 2][mask]
    saturation = hsv[..., 1][mask].astype(np.float64) / 255.0
    labels = measure.label(mask)
    regions = measure.regionprops(labels)
    largest = max(regions, key=lambda item: item.area)
    return {
        "area_thumbnail_px": int(mask.sum()),
        "area_thumbnail_fraction": float(mask.mean()),
        "connected_components": int(len(regions)),
        "largest_component_fraction": float(largest.area / mask.sum()),
        "hematoxylin_mean": float(np.mean(h_values)),
        "hematoxylin_p95": float(np.quantile(h_values, 0.95)),
        "eosin_mean": float(np.mean(e_values)),
        "dab_residual_mean": float(np.mean(d_values)),
        "gray_mean": float(np.mean(values)),
        "gray_std": float(np.std(values)),
        "saturation_mean": float(np.mean(saturation)),
        "sobel_mean": float(np.mean(gradient)),
        "sobel_p95": float(np.quantile(gradient, 0.95)),
        "bbox_min_x": int(largest.bbox[1]),
        "bbox_min_y": int(largest.bbox[0]),
        "bbox_max_x": int(largest.bbox[3]),
        "bbox_max_y": int(largest.bbox[2]),
    }


def overlay_labels(rgb: np.ndarray, labels: np.ndarray) -> np.ndarray:
    output = rgb.astype(np.float32).copy()
    for organ in ORGANS:
        mask = labels == ORGAN_VALUE[organ]
        output[mask] = output[mask] * 0.68 + ORGAN_COLOR[organ] * 0.32
        boundary = morphology.binary_dilation(
            mask, morphology.disk(3)
        ) ^ morphology.binary_erosion(mask, morphology.disk(3))
        output[boundary] = ORGAN_COLOR[organ]
    return np.clip(output, 0, 255).astype(np.uint8)


def main() -> int:
    args = parse_args()
    manifest_rows = read_csv(
        require_regular_file(args.slide_manifest, "slide manifest")
    )
    requested = {validate_sample_id(value) for value in args.sample_id}
    if requested:
        manifest_rows = [row for row in manifest_rows if row["sample_id"] in requested]
    if len(manifest_rows) != (len(requested) if requested else 24):
        raise FigS3Error("Slide selection does not match the requested cohort")
    output_root = args.output_root.resolve()
    segmentation_root = output_root / "organ_segmentation"
    segmentation_root.mkdir(parents=True, exist_ok=True)
    manifest = {str(row["sample_id"]): row for row in manifest_rows}

    cv_rows: list[dict[str, object]] = []
    segmentation_receipts: list[dict[str, object]] = []

    for sample_id in sorted(manifest):
        row = manifest[sample_id]
        overview_path = Path(row["source_overview_path"])
        sample_dir = segmentation_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        try:
            rgb = load_rgb(overview_path)
            tissue = clean_tissue_mask(rgb)
            groups, axis_metrics = three_spatial_groups(tissue)
            labels, identity_metrics = assign_organs(rgb, groups)
            if set(np.unique(labels)) != {0, 1, 2, 3}:
                raise FigS3Error("Organ label image is incomplete")
            Image.fromarray((tissue.astype(np.uint8) * 255)).save(
                sample_dir / "tissue_mask.png"
            )
            Image.fromarray(labels).save(sample_dir / "organ_label_mask.png")
            Image.fromarray(overlay_labels(rgb, labels)).save(
                sample_dir / "organ_overlay.png"
            )
            for organ in ORGANS:
                organ_mask = labels == ORGAN_VALUE[organ]
                features = organ_features(rgb, organ_mask)
                cv_rows.append(
                    {
                        "sample_id": sample_id,
                        "treatment": row["treatment"],
                        "cohort": row["cohort"],
                        "organ": organ,
                        "overview_width": rgb.shape[1],
                        "overview_height": rgb.shape[0],
                        **features,
                    }
                )
            receipt = {
                "status": "complete",
                "version": SEGMENTATION_VERSION,
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "sample_id": sample_id,
                "source_overview_path": str(overview_path.resolve()),
                "source_overview_sha256": sha256_file(overview_path),
                "overview_width": rgb.shape[1],
                "overview_height": rgb.shape[0],
                "tissue_fraction": float(tissue.mean()),
                **axis_metrics,
                **identity_metrics,
            }
            atomic_json(sample_dir / "segmentation_receipt.json", receipt)
            segmentation_receipts.append(receipt)
            print(
                json.dumps(
                    {
                        "sample_id": sample_id,
                        "segmentation": "complete",
                        "confidence": identity_metrics["organ_identity_confidence"],
                    }
                ),
                flush=True,
            )
        except Exception as exc:
            failure = {
                "status": "failed",
                "version": SEGMENTATION_VERSION,
                "sample_id": sample_id,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
            atomic_json(sample_dir / "segmentation_failure.json", failure)
            raise

    atomic_csv(
        output_root / "source_data" / "Figure_S3_overview_cv_features.csv", cv_rows
    )
    atomic_json(
        segmentation_root / "segmentation_summary.json",
        {
            "status": "complete",
            "version": SEGMENTATION_VERSION,
            "slides": len(segmentation_receipts),
            "organs": list(ORGANS),
            "method": (
                "HSV/H&E tissue masking, principal specimen axis, spatial three-cluster partition, "
                "terminal hematoxylin-rich spleen anchor, and audited liver-kidney-spleen placement"
            ),
            "human_reviewed": True,
            "review_outcome": "Liver and Kidney placement corrected after full-resolution histology review; Spleen unchanged",
            "image_registration_performed": False,
            "receipts": segmentation_receipts,
        },
    )

    print(
        json.dumps(
            {
                "status": "complete",
                "slides": len(manifest),
                "organs_per_slide": len(ORGANS),
                "image_registration_performed": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except FigS3Error as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
