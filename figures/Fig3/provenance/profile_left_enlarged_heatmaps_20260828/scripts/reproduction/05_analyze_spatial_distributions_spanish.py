#!/usr/bin/env python3
"""Exploratory animal-level spatial-distribution analysis for Figure 3.

This companion analysis uses the existing Cellpose-derived DAPI, c-FOS, and
NPY instance masks.  It does not train a condition classifier.  For each
animal, DAPI centroids define six covariance-normalized anillo radials that
contain approximately equal numbers of DAPI nuclei.  c-FOS/DAPI and
c-FOS+NPY/DAPI occurrence rates are calculated in every shell, transformed by
arcsin(sqrt(rate)), and compared across Water, Sucrose, and Allulose using an
exhaustive animal-label PERMANOVA (all 1,680 allocations of 3/3/3 animals).

Cells and shells are feature-extraction units only.  The biological animal is
the sole inferential unit.  This is explicitly exploratory because n=3 animals
per condition, one field per animal, and acquisition/model heterogeneity limit
generalization.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import json
import math
import os
import platform
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy
from scipy.ndimage import center_of_mass
import skimage


HERE = Path(__file__).resolve()
PAPER_ROOT = next((path for path in HERE.parents if path.name == "Paper"), None)
if PAPER_ROOT is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate Paper above {HERE}")

BASE_ANALYSIS_SCRIPT = HERE.with_name("01_analyze_cfos_npy.py")
RAW_ROOT = PAPER_ROOT / "Fig3" / "raw" / "legacy_experiment_2025_07_28"
DEFAULT_OUTPUT = PAPER_ROOT / "analyses" / "Fig3" / "results" / "spatial_distribution_run"

CONDITIONS = ("Water", "Sucrose", "Allulose")
CONDITION_COLORS = {"Water": "#b9e3f2", "Sucrose": "#e31a1c", "Allulose": "#2ecc71"}
DISPLAY_CONDITION_ES = {"Water": "Agua", "Sucrose": "Sacarosa", "Allulose": "Alulosa"}
MEAN_LINE_COLORS = {"Water": "#2f91b8", "Sucrose": "#c51b1f", "Allulose": "#159447"}
N_RADIAL_SHELLS = 6
EXPECTED_ANIMALS_PER_CONDITION = 3
EXPECTED_LABELINGS = 1680
ANALYSIS_VERSION = "fig3_spatial_distribution_v1.3_2026-08-28"
# Display-only F/G source layout. Keep the inferential radial profile on the
# left and use the available right-side space for a larger descriptive heatmap.
SPATIAL_PANEL_FIGSIZE = (6.40, 3.12)
SPATIAL_PANEL_GRID_WIDTH_RATIOS = (2.55, 1.45)
SPATIAL_PANEL_GRID_WSPACE = 0.26
SPATIAL_PANEL_MARGINS = {"left": 0.112, "right": 0.985, "bottom": 0.205, "top": 0.79}
SPATIAL_PROFILE_MIN_WIDTH_INCHES = 3.00
SPATIAL_HEATMAP_MIN_WIDTH_INCHES = 1.50

OUTPUT_NAMES = (
    "animal_spatial_counts.csv",
    "animal_radial_occurrence_features.csv",
    "exact_spatial_permanova.csv",
    "exact_spatial_permutation_receipts.csv",
    "exact_spatial_dispersion.csv",
    "exact_spatial_dispersion_permutation_receipts.csv",
    "segmentation_model_setting_audit.csv",
    "Figure3_Spatial_cFOS_occurrence.png",
    "Figure3_Spatial_cFOS_occurrence.pdf",
    "Figure3_Spatial_cFOS_NPY_occurrence.png",
    "Figure3_Spatial_cFOS_NPY_occurrence.pdf",
    "spatial_distribution_provenance.json",
    "output_manifest.csv",
)


@dataclass(frozen=True)
class EndpointSpec:
    endpoint: str
    label: str
    plot_title: str
    y_label: str
    output_stem: str


ENDPOINTS = (
    EndpointSpec(
        "cfos_occurrence",
        "Ocurrencia espacial de núcleos c-FOS positivos / núcleos DAPI",
        "Ocurrencia espacial de c-FOS",
        "Núcleos c-FOS⁺ / núcleos DAPI (%)",
        "Figure3_Spatial_cFOS_occurrence",
    ),
    EndpointSpec(
        "double_occurrence",
        "Ocurrencia espacial de núcleos dobles c-FOS y NPY positivos / núcleos DAPI",
        "Ocurrencia espacial c-FOS⁺NPY⁺",
        "Núcleos c-FOS⁺NPY⁺ / núcleos DAPI (%)",
        "Figure3_Spatial_cFOS_NPY_occurrence",
    ),
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace a prior output bundle only when it contains this script's known files.",
    )
    return parser


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_to_paper(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(PAPER_ROOT.resolve()))
    except ValueError:
        return str(resolved)


def load_base_analysis() -> Any:
    if not BASE_ANALYSIS_SCRIPT.is_file():
        raise FileNotFoundError(f"Missing Figure 3 base analysis: {BASE_ANALYSIS_SCRIPT}")
    module_name = "_fig3_base_analysis_for_spatial_distribution"
    spec = importlib.util.spec_from_file_location(module_name, BASE_ANALYSIS_SCRIPT)
    if spec is None or spec.loader is None:  # pragma: no cover
        raise RuntimeError(f"Could not load module specification: {BASE_ANALYSIS_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def validate_output_target(path: Path) -> Path:
    requested = path.expanduser()
    if requested.exists() and requested.is_symlink():
        raise RuntimeError(f"Refusing output through symbolic link: {requested}")
    resolved = requested.resolve()
    allowed = (PAPER_ROOT / "analyses" / "Fig3" / "results").resolve()
    if resolved == allowed or allowed not in resolved.parents:
        raise RuntimeError(f"Output must be a child of {allowed}: {resolved}")
    return resolved


def validate_existing_output(path: Path, force: bool) -> None:
    if not path.exists():
        return
    if not path.is_dir() or path.is_symlink():
        raise RuntimeError(f"Output exists but is not a normal directory: {path}")
    unknown = sorted(item.name for item in path.iterdir() if item.name not in OUTPUT_NAMES)
    if unknown:
        raise RuntimeError("Refusing output bundle with unknown files: " + ", ".join(unknown))
    if any(path.iterdir()) and not force:
        raise RuntimeError(f"Output already exists; rerun with --force: {path}")


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temp, index=False, float_format="%.12g", lineterminator="\n")
    os.replace(temp, path)


def atomic_json(payload: Mapping[str, Any], path: Path) -> None:
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, path)


def replace_bundle(stage: Path, target: Path) -> None:
    backup: Path | None = None
    if target.exists():
        backup = target.with_name(target.name + ".backup")
        if backup.exists():
            raise RuntimeError(f"Refusing stale backup directory: {backup}")
        os.replace(target, backup)
    try:
        os.replace(stage, target)
    except Exception:
        if backup is not None and backup.exists() and not target.exists():
            os.replace(backup, target)
        raise
    if backup is not None and backup.exists():
        shutil.rmtree(backup)


def label_centroids(labels: np.ndarray) -> np.ndarray:
    maximum = int(labels.max())
    if maximum <= 0:
        raise ValueError("DAPI label field is empty")
    positive = np.unique(labels[labels > 0])
    expected = np.arange(1, maximum + 1, dtype=positive.dtype)
    if not np.array_equal(positive, expected):
        raise ValueError("DAPI labels must be contiguous for label-to-centroid indexing")
    label_ids = np.arange(1, maximum + 1, dtype=int)
    points = np.asarray(
        center_of_mass(np.ones(labels.shape, dtype=np.uint8), labels, label_ids),
        dtype=np.float64,
    )
    if points.shape != (maximum, 2) or not np.isfinite(points).all():
        raise RuntimeError("Could not derive one finite centroid for every DAPI label")
    return points


def covariance_normalized_radial_shells(
    dapi_points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    centered = dapi_points - np.mean(dapi_points, axis=0, keepdims=True)
    covariance = np.cov(centered, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    if not np.isfinite(eigenvalues).all() or float(np.max(eigenvalues)) <= 0:
        raise ValueError("Degenerate DAPI centroid covariance")
    eigenvalues = np.maximum(eigenvalues, float(np.max(eigenvalues)) * 1e-8)
    whitened = centered @ eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues))
    radius = np.sqrt(np.sum(whitened * whitened, axis=1))
    finite_edges = np.quantile(radius, np.linspace(0.0, 1.0, N_RADIAL_SHELLS + 1))
    if np.any(np.diff(finite_edges) <= 0):
        raise ValueError("Repeated DAPI radial quantiles prevent six-shell construction")
    assignment_edges = finite_edges.copy()
    assignment_edges[0] = -np.inf
    assignment_edges[-1] = np.inf
    shell_index = np.searchsorted(assignment_edges[1:-1], radius, side="right")
    if int(shell_index.min()) != 0 or int(shell_index.max()) != N_RADIAL_SHELLS - 1:
        raise RuntimeError("Radial shell indexing did not cover all six shells")
    return shell_index.astype(np.int16), radius, finite_edges, eigenvalues


def positive_indicator(positive_ids: set[int], nucleus_count: int) -> np.ndarray:
    indicator = np.zeros(nucleus_count, dtype=bool)
    if positive_ids:
        indices = np.fromiter((identifier - 1 for identifier in sorted(positive_ids)), dtype=int)
        if int(indices.min()) < 0 or int(indices.max()) >= nucleus_count:
            raise ValueError("Positive nucleus identifier is outside DAPI label range")
        indicator[indices] = True
    return indicator


def channel_audit_row(
    module: Any,
    sample: Any,
    channel: str,
    segmentation_relative: str,
    segmentation_path: Path,
    metadata: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "animal_id": sample.animal_id,
        "sample_id": sample.sample_id,
        "condition": sample.condition,
        "acquisition_batch": sample.acquisition_batch,
        "channel": channel,
        "segmentation_path": relative_to_paper(segmentation_path),
        "segmentation_relative_mapping": segmentation_relative,
        "segmentation_sha256": sha256_file(segmentation_path),
        "segmentation_bytes": int(segmentation_path.stat().st_size),
        "mask_shape_y": int(metadata["mask_shape_y"]),
        "mask_shape_x": int(metadata["mask_shape_x"]),
        "mask_dtype": metadata["mask_dtype"],
        "mask_objects": int(metadata["object_count"]),
        "labels_contiguous": bool(metadata["labels_contiguous"]),
        "model_path_recorded": str(metadata["model_path"]),
        "model_basename_recorded": Path(str(metadata["model_path"])).name,
        "flow_threshold_recorded": metadata["flow_threshold"],
        "cellprob_threshold_recorded": metadata["cellprob_threshold"],
        "diameter_recorded": metadata["diameter"],
        "manual_label_count": int(metadata["manual_label_count"]),
        "manual_flag_count": int(metadata["manual_flag_count"]),
        "manual_change_record_count": int(metadata["manual_change_record_count"]),
        "cellpose_source_filename_recorded": str(metadata["cellpose_source_filename"]),
        "classification_rule": (
            "c-FOS ROI area >=5 px; best DAPI overlap fraction >=0.40"
            if channel == "c-FOS"
            else (
                "NPY union overlap >=1 px and >=0.05 of DAPI nucleus"
                if channel == "NPY"
                else "all positive contiguous DAPI labels"
            )
        ),
        "used_for_spatial_analysis": True,
        "base_analysis_version_source": relative_to_paper(BASE_ANALYSIS_SCRIPT),
        "base_cfos_min_roi_area_px": int(module.CFOS_MIN_ROI_AREA_PX),
        "base_cfos_dapi_overlap_fraction": float(module.CFOS_DAPI_OVERLAP_FRACTION),
        "base_npy_min_overlap_px": int(module.NPY_MIN_OVERLAP_PX),
        "base_npy_dapi_overlap_fraction": float(module.NPY_DAPI_OVERLAP_FRACTION),
    }


def extract_animal_features(
    module: Any,
    raw_root: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    count_rows: list[dict[str, Any]] = []
    feature_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []

    if tuple(module.CONDITIONS) != CONDITIONS:
        raise RuntimeError(f"Base analysis condition order changed: {module.CONDITIONS}")
    samples = tuple(module.SAMPLES)
    condition_counts = {
        condition: sum(sample.condition == condition for sample in samples)
        for condition in CONDITIONS
    }
    expected_counts = {condition: EXPECTED_ANIMALS_PER_CONDITION for condition in CONDITIONS}
    if condition_counts != expected_counts:
        raise RuntimeError(f"Expected 3 animals per condition; observed {condition_counts}")

    for animal_order, sample in enumerate(samples, start=1):
        channel_specs = (
            ("DAPI", sample.dapi_seg),
            ("c-FOS", sample.cfos_seg),
            ("NPY", sample.npy_seg),
        )
        labels: dict[str, np.ndarray] = {}
        metadata_by_channel: dict[str, Mapping[str, Any]] = {}
        for channel, relative_path in channel_specs:
            segmentation_path = raw_root / relative_path
            mask, metadata = module.load_cellpose_mask(segmentation_path)
            labels[channel] = mask
            metadata_by_channel[channel] = metadata
            audit_rows.append(
                channel_audit_row(
                    module,
                    sample,
                    channel,
                    relative_path,
                    segmentation_path,
                    metadata,
                )
            )

        shapes = {channel: tuple(mask.shape) for channel, mask in labels.items()}
        if len(set(shapes.values())) != 1:
            raise ValueError(f"Within-animal segmentation shape mismatch for {sample.animal_id}: {shapes}")

        dapi = labels["DAPI"]
        cfos = labels["c-FOS"]
        npy = labels["NPY"]
        cfos_ids, accepted_cfos_rois = module.cfos_positive_nuclei(dapi, cfos)
        npy_ids = module.marker_positive_nuclei(
            dapi,
            npy,
            min_overlap_px=module.NPY_MIN_OVERLAP_PX,
            min_fraction=module.NPY_DAPI_OVERLAP_FRACTION,
        )
        double_ids = cfos_ids & npy_ids

        dapi_points = label_centroids(dapi)
        shell_index, _radius, shell_edges, covariance_eigenvalues = (
            covariance_normalized_radial_shells(dapi_points)
        )
        indicators = {
            "cfos_occurrence": positive_indicator(cfos_ids, dapi_points.shape[0]),
            "double_occurrence": positive_indicator(double_ids, dapi_points.shape[0]),
        }

        shell_denominators: list[int] = []
        for shell_zero in range(N_RADIAL_SHELLS):
            shell_selected = shell_index == shell_zero
            denominator = int(np.count_nonzero(shell_selected))
            if denominator <= 0:
                raise RuntimeError(f"Empty anillo radial {shell_zero + 1}: {sample.animal_id}")
            shell_denominators.append(denominator)
            for endpoint in ENDPOINTS:
                positive_count = int(np.count_nonzero(indicators[endpoint.endpoint] & shell_selected))
                raw_rate = float(positive_count / denominator)
                transformed = float(np.arcsin(np.sqrt(np.clip(raw_rate, 0.0, 1.0))))
                feature_rows.append({
                    "animal_order": animal_order,
                    "animal_id": sample.animal_id,
                    "sample_id": sample.sample_id,
                    "condition": sample.condition,
                    "condition_order": CONDITIONS.index(sample.condition) + 1,
                    "experimental_unit": "biological animal",
                    "field_count": 1,
                    "endpoint": endpoint.endpoint,
                    "endpoint_label": endpoint.label,
                    "radial_shell": shell_zero + 1,
                    "shell_dapi_quantile_start": shell_zero / N_RADIAL_SHELLS,
                    "shell_dapi_quantile_end": (shell_zero + 1) / N_RADIAL_SHELLS,
                    "shell_dapi_nuclei": denominator,
                    "shell_positive_nuclei": positive_count,
                    "occurrence_rate": raw_rate,
                    "occurrence_percent": 100.0 * raw_rate,
                    "arcsin_sqrt_occurrence": transformed,
                    "shell_inner_whitened_radius": float(shell_edges[shell_zero]),
                    "shell_outer_whitened_radius": float(shell_edges[shell_zero + 1]),
                    "inference_role": "one component of one animal-level six-shell vector",
                })

        if max(shell_denominators) - min(shell_denominators) > 1:
            raise RuntimeError(
                f"DAPI quantile shells are not equal-density for {sample.animal_id}: {shell_denominators}"
            )
        count_rows.append({
            "animal_order": animal_order,
            "animal_id": sample.animal_id,
            "sample_id": sample.sample_id,
            "condition": sample.condition,
            "condition_order": CONDITIONS.index(sample.condition) + 1,
            "experimental_unit": "biological animal",
            "field_count": 1,
            "acquisition_batch": sample.acquisition_batch,
            "image_height_px": int(dapi.shape[0]),
            "image_width_px": int(dapi.shape[1]),
            "dapi_nuclei": int(dapi_points.shape[0]),
            "accepted_cfos_roi_objects": int(accepted_cfos_rois),
            "cfos_positive_nuclei": int(len(cfos_ids)),
            "npy_positive_nuclei": int(len(npy_ids)),
            "cfos_npy_double_positive_nuclei": int(len(double_ids)),
            "cfos_over_dapi": float(len(cfos_ids) / dapi_points.shape[0]),
            "cfos_npy_over_dapi": float(len(double_ids) / dapi_points.shape[0]),
            "cfos_npy_over_npy": float(len(double_ids) / len(npy_ids)),
            "radial_shell_count": N_RADIAL_SHELLS,
            "radial_shell_min_dapi_nuclei": min(shell_denominators),
            "radial_shell_max_dapi_nuclei": max(shell_denominators),
            "dapi_centroid_covariance_eigenvalue_small": float(covariance_eigenvalues[0]),
            "dapi_centroid_covariance_eigenvalue_large": float(covariance_eigenvalues[1]),
            "animal_id_status": sample.animal_id_status,
        })

    counts = pd.DataFrame(count_rows).sort_values("animal_order").reset_index(drop=True)
    features = pd.DataFrame(feature_rows).sort_values(
        ["animal_order", "endpoint", "radial_shell"]
    ).reset_index(drop=True)
    audit = pd.DataFrame(audit_rows).sort_values(["animal_id", "channel"]).reset_index(drop=True)

    model_modes = (
        audit.groupby("channel", sort=False)["model_basename_recorded"]
        .agg(lambda values: values.value_counts().index[0])
        .to_dict()
    )
    audit["channel_model_mode"] = audit["channel"].map(model_modes)
    audit["matches_channel_model_mode"] = (
        audit["model_basename_recorded"] == audit["channel_model_mode"]
    )
    audit["model_setting_caveat"] = ""
    dapi_mismatch = audit["channel"].eq("DAPI") & ~audit["matches_channel_model_mode"]
    audit.loc[dapi_mismatch, "model_setting_caveat"] = (
        "DAPI segmentation used dapi_Keyence4 instead of modal dapi_channel; DAPI labels "
        "define every shell and denominator, so both spatial endpoints remain exploratory"
    )
    npy_mismatch = audit["channel"].eq("NPY") & ~audit["matches_channel_model_mode"]
    audit.loc[npy_mismatch, "model_setting_caveat"] = (
        "NPY segmentation used dapi_channel instead of modal cfos_channel3; retain as recorded "
        "mask and interpret c-FOS+NPY spatial result as exploratory"
    )
    return counts, features, audit


def feature_matrix(
    features: pd.DataFrame,
    animals: pd.DataFrame,
    endpoint: str,
) -> np.ndarray:
    selected = features.loc[features["endpoint"].eq(endpoint)]
    pivot = selected.pivot(index="animal_id", columns="radial_shell", values="arcsin_sqrt_occurrence")
    pivot = pivot.reindex(animals["animal_id"].tolist())
    expected_columns = list(range(1, N_RADIAL_SHELLS + 1))
    if list(pivot.columns) != expected_columns or pivot.isna().any().any():
        raise RuntimeError(f"Incomplete six-shell feature matrix for {endpoint}")
    matrix = pivot.to_numpy(dtype=np.float64)
    if matrix.shape != (len(animals), N_RADIAL_SHELLS) or not np.isfinite(matrix).all():
        raise RuntimeError(f"Invalid feature matrix for {endpoint}: {matrix.shape}")
    return matrix


def pseudo_f_and_r_squared(matrix: np.ndarray, labels: np.ndarray) -> tuple[float, float]:
    grand_mean = np.mean(matrix, axis=0)
    total_sum_squares = float(np.sum((matrix - grand_mean) ** 2))
    between_sum_squares = 0.0
    for group_index in range(len(CONDITIONS)):
        group = matrix[labels == group_index]
        if group.shape[0] == 0:
            raise ValueError("A PERMANOVA allocation contains an empty condition")
        between_sum_squares += float(
            group.shape[0] * np.sum((np.mean(group, axis=0) - grand_mean) ** 2)
        )
    within_sum_squares = total_sum_squares - between_sum_squares
    df_between = len(CONDITIONS) - 1
    df_within = matrix.shape[0] - len(CONDITIONS)
    if within_sum_squares <= 0:
        pseudo_f = math.inf
    else:
        pseudo_f = (between_sum_squares / df_between) / (within_sum_squares / df_within)
    r_squared = between_sum_squares / total_sum_squares if total_sum_squares > 0 else math.nan
    return float(pseudo_f), float(r_squared)


def enumerate_condition_allocations(animal_count: int) -> Iterable[np.ndarray]:
    if animal_count != 9:
        raise ValueError(f"The exhaustive 3/3/3 enumeration requires 9 animals, got {animal_count}")
    indices = tuple(range(animal_count))
    for water_indices in itertools.combinations(indices, EXPECTED_ANIMALS_PER_CONDITION):
        remaining = tuple(index for index in indices if index not in water_indices)
        for sucrose_indices in itertools.combinations(
            remaining, EXPECTED_ANIMALS_PER_CONDITION
        ):
            labels = np.full(animal_count, 2, dtype=np.int8)
            labels[list(water_indices)] = 0
            labels[list(sucrose_indices)] = 1
            yield labels


def allocation_text(labels: np.ndarray, animal_ids: Sequence[str], group_index: int) -> str:
    return ";".join(
        animal_id for animal_id, label in zip(animal_ids, labels) if int(label) == group_index
    )


def benjamini_hochberg(raw_p_values: Sequence[float]) -> np.ndarray:
    values = np.asarray(raw_p_values, dtype=np.float64)
    if values.ndim != 1 or values.size == 0 or np.any((values < 0) | (values > 1)):
        raise ValueError("BH adjustment requires a non-empty one-dimensional p-value vector")
    order = np.argsort(values, kind="mergesort")
    ranked = values[order]
    adjusted_ranked = ranked * values.size / np.arange(1, values.size + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.minimum(adjusted_ranked, 1.0)
    return adjusted


def exact_permanova(
    counts: pd.DataFrame,
    features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    animal_ids = counts["animal_id"].astype(str).tolist()
    condition_to_index = {condition: index for index, condition in enumerate(CONDITIONS)}
    observed_labels = counts["condition"].map(condition_to_index).to_numpy(dtype=np.int8)
    if tuple(np.bincount(observed_labels, minlength=3)) != (3, 3, 3):
        raise RuntimeError("Observed Figure 3 labels are not a 3/3/3 allocation")
    allocations = list(enumerate_condition_allocations(len(counts)))
    if len(allocations) != EXPECTED_LABELINGS:
        raise RuntimeError(f"Expected 1,680 allocations; enumerated {len(allocations)}")

    summary_rows: list[dict[str, Any]] = []
    receipt_rows: list[dict[str, Any]] = []
    for endpoint in ENDPOINTS:
        matrix = feature_matrix(features, counts, endpoint.endpoint)
        observed_f, observed_r_squared = pseudo_f_and_r_squared(matrix, observed_labels)
        endpoint_receipts: list[dict[str, Any]] = []
        for labeling_index, labels in enumerate(allocations, start=1):
            permuted_f, permuted_r_squared = pseudo_f_and_r_squared(matrix, labels)
            endpoint_receipts.append({
                "endpoint": endpoint.endpoint,
                "endpoint_label": endpoint.label,
                "labeling_index": labeling_index,
                "water_animals": allocation_text(labels, animal_ids, 0),
                "sucrose_animals": allocation_text(labels, animal_ids, 1),
                "allulose_animals": allocation_text(labels, animal_ids, 2),
                "pseudo_f": permuted_f,
                "r_squared": permuted_r_squared,
                "is_observed_allocation": bool(np.array_equal(labels, observed_labels)),
                "is_as_or_more_extreme": bool(permuted_f >= observed_f - 1e-12),
                "permutation_unit": "biological animal",
            })
        extreme = sum(bool(row["is_as_or_more_extreme"]) for row in endpoint_receipts)
        observed_receipts = sum(bool(row["is_observed_allocation"]) for row in endpoint_receipts)
        if observed_receipts != 1:
            raise RuntimeError(f"Observed allocation appears {observed_receipts} times for {endpoint.endpoint}")
        p_value = float(extreme / len(endpoint_receipts))
        summary_rows.append({
            "endpoint": endpoint.endpoint,
            "endpoint_label": endpoint.label,
            "test": "exact animal-label PERMANOVA",
            "comparison": "Water|Sucrose|Allulose",
            "n_animals": len(counts),
            "n_water": 3,
            "n_sucrose": 3,
            "n_allulose": 3,
            "feature_count": N_RADIAL_SHELLS,
            "feature_transform": "arcsin(sqrt(shell positive nuclei / shell DAPI nuclei))",
            "distance": "Euclidean on six-shell animal vectors",
            "pseudo_f": observed_f,
            "r_squared": observed_r_squared,
            "p_value_exact": p_value,
            "extreme_labelings": extreme,
            "enumerated_labelings": len(endpoint_receipts),
            "permutation_exhaustive": True,
            "experimental_unit": "biological animal",
            "cells_as_independent_replicates": False,
            "analysis_tier": "exploratory",
        })
        receipt_rows.extend(endpoint_receipts)

    summary = pd.DataFrame(summary_rows)
    summary["bh_q_value_two_endpoint_family"] = benjamini_hochberg(
        summary["p_value_exact"].to_numpy(dtype=float)
    )
    summary["multiple_testing_family"] = (
        "two declared global spatial occurrence profiles: c-FOS/DAPI and c-FOS+NPY/DAPI"
    )
    receipts = pd.DataFrame(receipt_rows)
    return summary, receipts



def dispersion_f_and_distances(
    matrix: np.ndarray,
    labels: np.ndarray,
) -> tuple[float, np.ndarray]:
    """One-way F on Euclidean distances to allocation-specific centroids."""
    distances = np.empty(matrix.shape[0], dtype=np.float64)
    for group_index in range(len(CONDITIONS)):
        selected = labels == group_index
        group = matrix[selected]
        if group.shape[0] != EXPECTED_ANIMALS_PER_CONDITION:
            raise ValueError("A dispersion allocation does not contain three animals per group")
        centroid = np.mean(group, axis=0)
        distances[selected] = np.linalg.norm(group - centroid, axis=1)
    grand_mean = float(np.mean(distances))
    between_sum_squares = 0.0
    within_sum_squares = 0.0
    for group_index in range(len(CONDITIONS)):
        group_distances = distances[labels == group_index]
        group_mean = float(np.mean(group_distances))
        between_sum_squares += float(
            group_distances.size * (group_mean - grand_mean) ** 2
        )
        within_sum_squares += float(np.sum((group_distances - group_mean) ** 2))
    if within_sum_squares <= 0:
        statistic = math.inf
    else:
        statistic = (
            (between_sum_squares / (len(CONDITIONS) - 1))
            / (within_sum_squares / (matrix.shape[0] - len(CONDITIONS)))
        )
    return float(statistic), distances


def exact_multivariate_dispersion(
    counts: pd.DataFrame,
    features: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Exact allocation diagnostic for heterogeneity of multivariate dispersion.

    For every one of the 1,680 labeled 3/3/3 allocations, group centroids are
    recomputed in the six-shell arcsin-square-root feature space. Euclidean
    distances from animals to their allocated centroid are then compared with
    a one-way ANOVA F statistic. The exact p value includes the observed
    allocation and is extreme/1,680.
    """
    animal_ids = counts["animal_id"].astype(str).tolist()
    condition_to_index = {condition: index for index, condition in enumerate(CONDITIONS)}
    observed_labels = counts["condition"].map(condition_to_index).to_numpy(dtype=np.int8)
    allocations = list(enumerate_condition_allocations(len(counts)))
    if len(allocations) != EXPECTED_LABELINGS:
        raise RuntimeError(f"Expected 1,680 dispersion allocations; got {len(allocations)}")

    summary_rows: list[dict[str, Any]] = []
    receipt_rows: list[dict[str, Any]] = []
    for endpoint in ENDPOINTS:
        matrix = feature_matrix(features, counts, endpoint.endpoint)
        observed_f, observed_distances = dispersion_f_and_distances(matrix, observed_labels)
        endpoint_receipts: list[dict[str, Any]] = []
        for labeling_index, labels in enumerate(allocations, start=1):
            statistic, distances = dispersion_f_and_distances(matrix, labels)
            group_means = [
                float(np.mean(distances[labels == group_index]))
                for group_index in range(len(CONDITIONS))
            ]
            endpoint_receipts.append({
                "endpoint": endpoint.endpoint,
                "endpoint_label": endpoint.label,
                "labeling_index": labeling_index,
                "water_animals": allocation_text(labels, animal_ids, 0),
                "sucrose_animals": allocation_text(labels, animal_ids, 1),
                "allulose_animals": allocation_text(labels, animal_ids, 2),
                "dispersion_f": statistic,
                "water_mean_distance": group_means[0],
                "sucrose_mean_distance": group_means[1],
                "allulose_mean_distance": group_means[2],
                "is_observed_allocation": bool(np.array_equal(labels, observed_labels)),
                "is_as_or_more_extreme": bool(statistic >= observed_f - 1e-12),
                "permutation_unit": "biological animal",
            })
        extreme = sum(bool(row["is_as_or_more_extreme"]) for row in endpoint_receipts)
        if sum(bool(row["is_observed_allocation"]) for row in endpoint_receipts) != 1:
            raise RuntimeError(
                f"Observed dispersion allocation is not unique for {endpoint.endpoint}"
            )
        observed_groups = [
            observed_distances[observed_labels == group_index]
            for group_index in range(len(CONDITIONS))
        ]
        summary_rows.append({
            "endpoint": endpoint.endpoint,
            "endpoint_label": endpoint.label,
            "test": "exact animal-label multivariate dispersion diagnostic",
            "comparison": "Water|Sucrose|Allulose",
            "n_animals": len(counts),
            "n_water": 3,
            "n_sucrose": 3,
            "n_allulose": 3,
            "feature_count": N_RADIAL_SHELLS,
            "feature_transform": "arcsin(sqrt(shell positive nuclei / shell DAPI nuclei))",
            "distance": "Euclidean in six-shell animal feature space",
            "centroid_definition": "allocation-specific condition centroid recomputed for every labeling",
            "statistic": "one-way ANOVA F on Euclidean distances to allocated group centroids",
            "dispersion_f": observed_f,
            "p_value_exact": float(extreme / len(endpoint_receipts)),
            "extreme_labelings": extreme,
            "enumerated_labelings": len(endpoint_receipts),
            "permutation_exhaustive": True,
            "experimental_unit": "biological animal",
            "cells_as_independent_replicates": False,
            "analysis_tier": "exploratory assumption diagnostic",
            "water_mean_distance": float(np.mean(observed_groups[0])),
            "water_sd_distance": float(np.std(observed_groups[0], ddof=1)),
            "sucrose_mean_distance": float(np.mean(observed_groups[1])),
            "sucrose_sd_distance": float(np.std(observed_groups[1], ddof=1)),
            "allulose_mean_distance": float(np.mean(observed_groups[2])),
            "allulose_sd_distance": float(np.std(observed_groups[2], ddof=1)),
            "small_sample_bias_adjustment": (
                "not applied; every allocated group has n=3, so the common sqrt(3/2) "
                "distance factor leaves F and p unchanged"
            ),
        })
        receipt_rows.extend(endpoint_receipts)

    summary = pd.DataFrame(summary_rows)
    summary["bh_q_value_two_endpoint_diagnostic_family"] = benjamini_hochberg(
        summary["p_value_exact"].to_numpy(dtype=float)
    )
    summary["multiple_testing_family"] = (
        "two declared multivariate-dispersion diagnostics corresponding to the two spatial endpoints"
    )
    return summary, pd.DataFrame(receipt_rows)

def p_text(value: float) -> str:
    if value < 0.001:
        return f"{value:.2e}"
    if value < 0.1:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    return f"{value:.3f}".rstrip("0").rstrip(".")


def configure_plot_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.0,
        "axes.titlesize": 10.0,
        "axes.labelsize": 8.6,
        "xtick.labelsize": 7.4,
        "ytick.labelsize": 7.4,
        "legend.fontsize": 7.2,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })


def validate_panel_layout(
    fig: plt.Figure,
    heat_axis: plt.Axes,
    profile_axis: plt.Axes,
) -> None:
    """Fail unless the profile is left of an enlarged condition-means heatmap."""
    fig.canvas.draw()
    heat_position = heat_axis.get_position()
    profile_position = profile_axis.get_position()
    profile_width_inches = profile_position.width * fig.get_figwidth()
    heat_width_inches = heat_position.width * fig.get_figwidth()
    if profile_position.x1 >= heat_position.x0:
        raise RuntimeError("Radial profile is not strictly left of the condition-means heatmap")
    if profile_width_inches < SPATIAL_PROFILE_MIN_WIDTH_INCHES:
        raise RuntimeError(
            f"Radial-profile axis is only {profile_width_inches:.3f} in wide; "
            f"at least {SPATIAL_PROFILE_MIN_WIDTH_INCHES:.3f} in is required"
        )
    if heat_width_inches < SPATIAL_HEATMAP_MIN_WIDTH_INCHES:
        raise RuntimeError(
            f"Condition-means heatmap is only {heat_width_inches:.3f} in wide; "
            f"at least {SPATIAL_HEATMAP_MIN_WIDTH_INCHES:.3f} in is required"
        )


def save_panel(
    features: pd.DataFrame,
    summary: pd.DataFrame,
    endpoint: EndpointSpec,
    output_dir: Path,
    dpi: int,
) -> None:
    configure_plot_style()
    selected = features.loc[features["endpoint"].eq(endpoint.endpoint)].copy()
    statistic = summary.loc[summary["endpoint"].eq(endpoint.endpoint)].iloc[0]

    fig = plt.figure(figsize=SPATIAL_PANEL_FIGSIZE, facecolor="white")
    grid = fig.add_gridspec(
        1, 2,
        width_ratios=SPATIAL_PANEL_GRID_WIDTH_RATIOS,
        wspace=SPATIAL_PANEL_GRID_WSPACE,
    )
    axis = fig.add_subplot(grid[0, 0])
    heat_axis = fig.add_subplot(grid[0, 1])
    fig.subplots_adjust(**SPATIAL_PANEL_MARGINS)
    x_values = np.arange(1, N_RADIAL_SHELLS + 1, dtype=float)

    for condition in CONDITIONS:
        condition_rows = selected.loc[selected["condition"].eq(condition)]
        animal_profiles: list[np.ndarray] = []
        for animal_id in condition_rows["animal_id"].drop_duplicates().tolist():
            profile = (
                condition_rows.loc[condition_rows["animal_id"].eq(animal_id)]
                .sort_values("radial_shell")["occurrence_percent"]
                .to_numpy(dtype=float)
            )
            if profile.size != N_RADIAL_SHELLS:
                raise RuntimeError(f"Incomplete plot profile for {animal_id}/{endpoint.endpoint}")
            animal_profiles.append(profile)
            axis.plot(
                x_values,
                profile,
                color=MEAN_LINE_COLORS[condition],
                alpha=0.20,
                linewidth=0.8,
                marker="o",
                markersize=2.5,
                markeredgewidth=0.0,
                zorder=2,
            )
        matrix = np.stack(animal_profiles)
        mean_profile = np.mean(matrix, axis=0)
        axis.plot(
            x_values,
            mean_profile,
            color=MEAN_LINE_COLORS[condition],
            linewidth=2.25,
            marker="o",
            markersize=4.5,
            markerfacecolor=CONDITION_COLORS[condition],
            markeredgecolor="#172126",
            markeredgewidth=0.55,
            zorder=4,
        )

    ymax = float(selected["occurrence_percent"].max())
    axis.set_ylim(0.0, max(1.0, ymax * 1.20))
    axis.set_xlim(0.72, N_RADIAL_SHELLS + 0.28)
    axis.set_xticks(x_values)
    axis.set_xticklabels(["1\ninterior", "2", "3", "4", "5", "6\nexterior"])
    axis.set_xlabel("Anillo radial de igual densidad DAPI", labelpad=3.0)
    axis.set_ylabel(endpoint.y_label, labelpad=4.0)
    axis.set_title(endpoint.plot_title, loc="left", fontweight="bold", pad=8.0)
    axis.text(
        0.98,
        0.955,
        (
            f"PERMANOVA exacta p = {p_text(float(statistic['p_value_exact']))}"
            f"  ·  BH q = {p_text(float(statistic['bh_q_value_two_endpoint_family']))}"
        ),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=7.3,
        fontweight="bold",
        color="#26353b",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.86, "pad": 1.7},
        zorder=8,
    )
    axis.grid(axis="y", color="#d9dee1", linewidth=0.6, alpha=0.9)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(width=0.75, length=3.0, color="#30393d")

    legend_handles = [
        Line2D(
            [0], [0],
            color=MEAN_LINE_COLORS[condition],
            marker="o",
            markerfacecolor=CONDITION_COLORS[condition],
            markeredgecolor="#172126",
            markeredgewidth=0.45,
            linewidth=2.0,
            markersize=4.0,
            label=DISPLAY_CONDITION_ES[condition],
        )
        for condition in CONDITIONS
    ]
    axis.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.19),
        ncol=3,
        frameon=False,
        columnspacing=1.3,
        handlelength=1.6,
    )
    condition_means = np.vstack([
        selected.loc[selected["condition"].eq(condition)]
        .groupby("radial_shell", sort=True)["occurrence_percent"]
        .mean()
        .reindex(range(1, N_RADIAL_SHELLS + 1))
        .to_numpy(dtype=float)
        for condition in CONDITIONS
    ])
    if condition_means.shape != (len(CONDITIONS), N_RADIAL_SHELLS):
        raise RuntimeError(f"Invalid descriptive heatmap for {endpoint.endpoint}")
    heat_max = float(np.max(condition_means))
    heat_image = heat_axis.imshow(
        condition_means,
        cmap="inferno",
        vmin=0.0,
        vmax=max(heat_max, 1e-12),
        aspect="auto",
        interpolation="nearest",
    )
    heat_axis.set_title("Medias por condición", fontsize=8.7, fontweight="bold", pad=7.0)
    heat_axis.set_xticks(np.arange(N_RADIAL_SHELLS))
    heat_axis.set_xticklabels(range(1, N_RADIAL_SHELLS + 1), fontsize=6.8)
    heat_axis.set_xlabel("anillo radial", fontsize=7.2, labelpad=2.0)
    heat_axis.set_yticks(np.arange(len(CONDITIONS)))
    heat_axis.set_yticklabels(("Ag", "Sa", "Al"), fontsize=7.4, fontweight="bold")
    for tick, condition in zip(heat_axis.get_yticklabels(), CONDITIONS):
        tick.set_color(MEAN_LINE_COLORS[condition])
    for row_index in range(len(CONDITIONS)):
        for column_index in range(N_RADIAL_SHELLS):
            value = float(condition_means[row_index, column_index])
            red, green, blue, _alpha = heat_image.cmap(heat_image.norm(value))
            luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
            text_color = "#172126" if luminance >= 0.56 else "white"
            heat_axis.text(
                column_index,
                row_index,
                f"{value:.1f}",
                ha="center",
                va="center",
                fontsize=5.8,
                fontweight="bold",
                color=text_color,
            )
    heat_axis.set_xticks(np.arange(-0.5, N_RADIAL_SHELLS, 1), minor=True)
    heat_axis.set_yticks(np.arange(-0.5, len(CONDITIONS), 1), minor=True)
    heat_axis.grid(which="minor", color="white", linewidth=0.55, alpha=0.72)
    heat_axis.tick_params(which="minor", bottom=False, left=False)
    heat_axis.tick_params(which="major", length=0)
    for spine in heat_axis.spines.values():
        spine.set_color("#45565d")
        spine.set_linewidth(0.65)
    colorbar = fig.colorbar(heat_image, ax=heat_axis, fraction=0.085, pad=0.055)
    colorbar.ax.set_title("%", fontsize=7.0, pad=2.0)
    colorbar.ax.tick_params(labelsize=6.2, width=0.6, length=2.0)
    colorbar.outline.set_linewidth(0.55)
    validate_panel_layout(fig, heat_axis, axis)

    png_path = output_dir / f"{endpoint.output_stem}.png"
    pdf_path = output_dir / f"{endpoint.output_stem}.pdf"
    fig.savefig(
        png_path,
        dpi=dpi,
        facecolor="white",
        bbox_inches=None,
        metadata={"Software": ANALYSIS_VERSION, "Title": endpoint.plot_title},
    )
    fig.savefig(
        pdf_path,
        dpi=dpi,
        facecolor="white",
        bbox_inches=None,
        metadata={
            "Title": endpoint.plot_title,
            "Author": "Figure 3 reproducible analysis",
            "Creator": ANALYSIS_VERSION,
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)


def validate_results(
    counts: pd.DataFrame,
    features: pd.DataFrame,
    summary: pd.DataFrame,
    receipts: pd.DataFrame,
    dispersion_summary: pd.DataFrame,
    dispersion_receipts: pd.DataFrame,
    audit: pd.DataFrame,
    output_dir: Path,
) -> None:
    if len(counts) != 9 or counts["animal_id"].nunique() != 9:
        raise RuntimeError("Spatial analysis must contain exactly nine unique animals")
    if counts.groupby("condition")["animal_id"].nunique().to_dict() != {
        condition: 3 for condition in CONDITIONS
    }:
        raise RuntimeError("Spatial analysis must contain three animals per condition")
    expected_feature_rows = 9 * len(ENDPOINTS) * N_RADIAL_SHELLS
    if len(features) != expected_feature_rows:
        raise RuntimeError(f"Expected {expected_feature_rows} feature rows; got {len(features)}")
    if not np.allclose(
        np.sin(features["arcsin_sqrt_occurrence"].to_numpy(dtype=float)) ** 2,
        features["occurrence_rate"].to_numpy(dtype=float),
        atol=1e-12,
        rtol=1e-10,
    ):
        raise RuntimeError("Stored arcsin-square-root features do not invert to occurrence rates")
    if len(summary) != len(ENDPOINTS):
        raise RuntimeError("Expected exactly two global PERMANOVA rows")
    if not (summary["enumerated_labelings"] == EXPECTED_LABELINGS).all():
        raise RuntimeError("A global test lacks all 1,680 allocation receipts")
    if len(receipts) != EXPECTED_LABELINGS * len(ENDPOINTS):
        raise RuntimeError("Full permutation receipt table has the wrong row count")
    if not (receipts.groupby("endpoint")["is_observed_allocation"].sum() == 1).all():
        raise RuntimeError("Every endpoint must have exactly one observed allocation receipt")
    expected_endpoints = {endpoint.endpoint for endpoint in ENDPOINTS}
    if (
        len(dispersion_summary) != len(ENDPOINTS)
        or set(dispersion_summary["endpoint"].astype(str)) != expected_endpoints
    ):
        raise RuntimeError("Expected exactly two declared dispersion-diagnostic rows")
    if not (dispersion_summary["enumerated_labelings"] == EXPECTED_LABELINGS).all():
        raise RuntimeError("A dispersion diagnostic lacks all 1,680 allocation receipts")
    if len(dispersion_receipts) != EXPECTED_LABELINGS * len(ENDPOINTS):
        raise RuntimeError("Full dispersion allocation receipt table has the wrong row count")
    if not (
        dispersion_receipts.groupby("endpoint")["is_observed_allocation"].sum() == 1
    ).all():
        raise RuntimeError("Every dispersion endpoint must have one observed receipt")
    observed_extreme = (
        dispersion_receipts.groupby("endpoint")["is_as_or_more_extreme"].sum()
    )
    summary_by_endpoint = dispersion_summary.set_index("endpoint")
    for endpoint in expected_endpoints:
        row = summary_by_endpoint.loc[endpoint]
        extreme = int(observed_extreme.loc[endpoint])
        if extreme != int(row["extreme_labelings"]):
            raise RuntimeError(f"Dispersion extreme-count mismatch for {endpoint}")
        if not math.isclose(
            float(row["p_value_exact"]), extreme / EXPECTED_LABELINGS,
            rel_tol=0.0, abs_tol=1e-15,
        ):
            raise RuntimeError(f"Dispersion exact-p mismatch for {endpoint}")
    expected_dispersion_q = benjamini_hochberg(
        dispersion_summary["p_value_exact"].to_numpy(dtype=float)
    )
    if not np.allclose(
        dispersion_summary["bh_q_value_two_endpoint_diagnostic_family"].to_numpy(dtype=float),
        expected_dispersion_q, atol=1e-15, rtol=0.0,
    ):
        raise RuntimeError("Dispersion BH q values do not reproduce")
    if len(audit) != 27 or audit["segmentation_sha256"].str.len().ne(64).any():
        raise RuntimeError("Expected 27 hashed segmentation audit rows")
    cfos_models = audit.loc[audit["channel"].eq("c-FOS"), "model_basename_recorded"].unique()
    if list(cfos_models) != ["cfos_channel3"]:
        raise RuntimeError(f"Unexpected c-FOS Cellpose model inventory: {cfos_models}")
    npy_nonmodal = audit.loc[
        audit["channel"].eq("NPY") & ~audit["matches_channel_model_mode"], "animal_id"
    ].tolist()
    if sorted(npy_nonmodal) != ["E7_FR7-5", "E8_FR6-1"]:
        raise RuntimeError(f"Unexpected NPY model heterogeneity inventory: {npy_nonmodal}")
    dapi_nonmodal = audit.loc[
        audit["channel"].eq("DAPI") & ~audit["matches_channel_model_mode"], "animal_id"
    ].tolist()
    if dapi_nonmodal != []:
        raise RuntimeError(f"Unexpected DAPI model heterogeneity inventory: {dapi_nonmodal}")
    for endpoint in ENDPOINTS:
        for suffix in (".png", ".pdf"):
            path = output_dir / f"{endpoint.output_stem}{suffix}"
            if not path.is_file() or path.stat().st_size <= 0:
                raise RuntimeError(f"Missing or empty panel asset: {path}")


def build_manifest(output_dir: Path, names: Sequence[str], recorded_output_dir: Path) -> pd.DataFrame:
    rows = []
    for name in names:
        path = output_dir / name
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Manifest source missing or empty: {path}")
        suffix = path.suffix.lower()
        role = (
            "panel_asset"
            if suffix in {".png", ".pdf"}
            else ("analysis_table" if suffix == ".csv" else "provenance")
        )
        rows.append({
            "role": role,
            "path": relative_to_paper(recorded_output_dir / name),
            "bytes": int(path.stat().st_size),
            "sha256": sha256_file(path),
        })
    return pd.DataFrame(rows)


def validate_promoted_manifest(output_dir: Path) -> None:
    manifest_path = output_dir / "output_manifest.csv"
    manifest = pd.read_csv(manifest_path)
    required = {"role", "path", "bytes", "sha256"}
    if not required.issubset(manifest.columns):
        raise RuntimeError(f"Output manifest lacks required columns: {required - set(manifest.columns)}")
    expected_paths = {
        (output_dir / name).resolve()
        for name in OUTPUT_NAMES
        if name != "output_manifest.csv"
    }
    observed_paths: set[Path] = set()
    for row in manifest.itertuples(index=False):
        recorded = Path(str(row.path))
        path = recorded.resolve() if recorded.is_absolute() else (PAPER_ROOT / recorded).resolve()
        observed_paths.add(path)
        if path.parent != output_dir.resolve():
            raise RuntimeError(f"Manifest path escapes the promoted bundle: {path}")
        if not path.is_file() or path.stat().st_size <= 0:
            raise RuntimeError(f"Manifest path is missing or empty after promotion: {path}")
        if int(row.bytes) != int(path.stat().st_size):
            raise RuntimeError(f"Manifest byte count mismatch after promotion: {path}")
        if str(row.sha256) != sha256_file(path):
            raise RuntimeError(f"Manifest SHA-256 mismatch after promotion: {path}")
    if observed_paths != expected_paths:
        raise RuntimeError(
            "Promoted manifest inventory mismatch: "
            f"missing={sorted(map(str, expected_paths - observed_paths))}, "
            f"unexpected={sorted(map(str, observed_paths - expected_paths))}"
        )


def main() -> None:
    args = build_parser().parse_args()
    if args.dpi < 150 or args.dpi > 1200:
        raise ValueError("--dpi must be between 150 and 1200")
    output_target = validate_output_target(args.output_dir)
    validate_existing_output(output_target, args.force)
    raw_root = args.raw_root.expanduser().resolve()
    if not raw_root.is_dir():
        raise FileNotFoundError(f"Missing raw root: {raw_root}")

    output_target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".fig3-spatial-stage-", dir=output_target.parent))
    try:
        module = load_base_analysis()
        counts, features, audit = extract_animal_features(module, raw_root)
        summary, receipts = exact_permanova(counts, features)
        dispersion_summary, dispersion_receipts = exact_multivariate_dispersion(
            counts, features
        )

        atomic_csv(counts, stage / "animal_spatial_counts.csv")
        atomic_csv(features, stage / "animal_radial_occurrence_features.csv")
        atomic_csv(summary, stage / "exact_spatial_permanova.csv")
        atomic_csv(receipts, stage / "exact_spatial_permutation_receipts.csv")
        atomic_csv(dispersion_summary, stage / "exact_spatial_dispersion.csv")
        atomic_csv(
            dispersion_receipts,
            stage / "exact_spatial_dispersion_permutation_receipts.csv",
        )
        atomic_csv(audit, stage / "segmentation_model_setting_audit.csv")
        for endpoint in ENDPOINTS:
            save_panel(features, summary, endpoint, stage, args.dpi)

        validate_results(
            counts, features, summary, receipts,
            dispersion_summary, dispersion_receipts, audit, stage,
        )
        npy_nonmodal = audit.loc[
            audit["channel"].eq("NPY") & ~audit["matches_channel_model_mode"]
        ]
        dapi_nonmodal = audit.loc[
            audit["channel"].eq("DAPI") & ~audit["matches_channel_model_mode"]
        ]
        provenance = {
            "analysis_version": ANALYSIS_VERSION,
            "analysis_tier": "exploratory",
            "analysis_description": (
                "Cellpose-derived animal-level spatial occurrence profiles; no newly trained "
                "condition classifier"
            ),
            "statistical_unit": "biological animal",
            "cells_as_independent_replicates": False,
            "animal_count": int(len(counts)),
            "condition_counts": counts.groupby("condition")["animal_id"].nunique().to_dict(),
            "endpoints": [endpoint.endpoint for endpoint in ENDPOINTS],
            "method": {
                "coordinate_source": "centroids of contiguous Cellpose DAPI instance labels",
                "coordinate_normalization": (
                    "center DAPI centroid cloud and whiten by its two covariance eigenvalues; "
                    "use orientation-invariant whitened radial distance"
                ),
                "radial_shells": N_RADIAL_SHELLS,
                "shell_definition": (
                    "within-animal quantiles of whitened DAPI radial distance; approximately "
                    "equal DAPI nuclei per shell"
                ),
                "occurrence_denominator": "DAPI nuclei in the same animal and anillo radial",
                "feature_transform": "arcsin(sqrt(positive nuclei / DAPI nuclei))",
                "distance": "Euclidean between six-shell animal vectors",
                "test": "PERMANOVA pseudo-F with exhaustive animal-label allocation",
                "enumerated_labelings_per_endpoint": EXPECTED_LABELINGS,
                "multiple_testing": (
                    "Benjamini-Hochberg across the two declared global spatial occurrence profiles"
                ),
            },
            "global_results": summary.to_dict(orient="records"),
            "multivariate_dispersion_diagnostic": {
                "summary_table": "exact_spatial_dispersion.csv",
                "allocation_receipts": (
                    "exact_spatial_dispersion_permutation_receipts.csv"
                ),
                "method": (
                    "For every 3/3/3 allocation, recompute each allocated group centroid "
                    "in six-shell feature space, calculate animal Euclidean distances to "
                    "that centroid, and apply a one-way ANOVA F to the distances"
                ),
                "exact_p": "allocations with dispersion F >= observed F / 1,680",
                "multiple_testing": (
                    "Benjamini-Hochberg across the two endpoint-matched diagnostics"
                ),
                "results": dispersion_summary.to_dict(orient="records"),
                "interpretation": (
                    "The double-occurrence diagnostic detects unequal multivariate "
                    "dispersion; its PERMANOVA result cannot be interpreted as a pure "
                    "condition-centroid difference"
                ),
            },
            "base_analysis_script": relative_to_paper(BASE_ANALYSIS_SCRIPT),
            "base_analysis_script_sha256": sha256_file(BASE_ANALYSIS_SCRIPT),
            "spatial_analysis_script": relative_to_paper(HERE),
            "spatial_analysis_script_sha256": sha256_file(HERE),
            "raw_root": relative_to_paper(raw_root),
            "segmentation_inputs": {
                "count": int(len(audit)),
                "audit_table": "segmentation_model_setting_audit.csv",
                "dapi_modal_model": str(
                    audit.loc[audit["channel"].eq("DAPI"), "channel_model_mode"].iloc[0]
                ),
                "dapi_nonmodal_animals": dapi_nonmodal["animal_id"].tolist(),
                "dapi_nonmodal_models": dapi_nonmodal[
                    ["animal_id", "condition", "model_basename_recorded"]
                ].to_dict(orient="records"),
                "all_cfos_model": bool(
                    audit.loc[audit["channel"].eq("c-FOS"), "model_basename_recorded"]
                    .eq("cfos_channel3")
                    .all()
                ),
                "npy_modal_model": str(
                    audit.loc[audit["channel"].eq("NPY"), "channel_model_mode"].iloc[0]
                ),
                "npy_nonmodal_animals": npy_nonmodal["animal_id"].tolist(),
                "npy_nonmodal_models": npy_nonmodal[
                    ["animal_id", "condition", "model_basename_recorded"]
                ].to_dict(orient="records"),
            },
            "caveats": [
                "Only three animals and one analyzed field per condition are available.",
                (
                    "Water acquisition batch differs from the Sucrose/Allulose legacy Apotome "
                    "batch; Water3 is a 2026 Leica SP8 acquisition. Mask-derived ratios reduce "
                    "but do not eliminate acquisition confounding."
                ),
                (
                    "The mandatory human-accepted WATER_NPY3 DAPI mask records the modal "
                    "dapi_channel model; its accepted labels define all shells and denominators."
                ),
                (
                    "NPY masks for E7_FR7-5 and E8_FR6-1 record dapi_channel rather than the "
                    "modal cfos_channel3 Cellpose model."
                ),
                (
                    "The double-occurrence endpoint has unequal multivariate dispersion by "
                    "the exact endpoint-matched diagnostic; PERMANOVA therefore does not "
                    "distinguish condition-centroid separation from dispersion heterogeneity."
                ),
                (
                    "The analysis was developed after data inspection and is exploratory, not "
                    "preregistered confirmatory inference."
                ),
            ],
            "registered_microscopy_insets": {
                "included": False,
                "reason": (
                    "Only one representative field per condition has an audited display "
                    "registration; using those three insets beside all-nine-animal inference "
                    "would imply unsupported section-level registration."
                ),
            },
            "descriptive_panel_insets": {
                "included": True,
                "type": "condition-by-six-shell mean occurrence heatmap",
                "source": "same animal_radial_occurrence_features.csv values shown in line profiles",
                "rows": ["Water", "Sucrose", "Allulose"],
                "columns": [1, 2, 3, 4, 5, 6],
                "cell_value": "arithmetic mean animal occurrence_percent within condition and shell",
                "color_scale": "inferno; zero to endpoint-specific maximum condition mean",
                "inferential_role": "descriptive only; no additional statistical test or replicate",
                "position": "right of the radial-profile plot",
                "panel_figure_inches": list(SPATIAL_PANEL_FIGSIZE),
                "grid_width_ratios_radial_profile_to_condition_means": list(SPATIAL_PANEL_GRID_WIDTH_RATIOS),
                "grid_wspace": SPATIAL_PANEL_GRID_WSPACE,
                "radial_profile_minimum_width_inches": SPATIAL_PROFILE_MIN_WIDTH_INCHES,
                "condition_means_heatmap_minimum_width_inches": SPATIAL_HEATMAP_MIN_WIDTH_INCHES,
            },

            "condition_colors": CONDITION_COLORS,
            "mean_line_colors": MEAN_LINE_COLORS,
            "environment": {
                "python": platform.python_version(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "scipy": scipy.__version__,
                "scikit_image": skimage.__version__,
                "matplotlib": matplotlib.__version__,
                "platform": platform.platform(),
            },
            "determinism": {
                "randomness_used": False,
                "allocation_order": "lexicographic combinations in fixed base sample order",
                "pdf_creation_and_modification_dates": None,
                "dpi": int(args.dpi),
            },
        }
        atomic_json(provenance, stage / "spatial_distribution_provenance.json")
        manifest_source_names = [
            name for name in OUTPUT_NAMES
            if name not in {"output_manifest.csv"}
        ]
        manifest = build_manifest(stage, manifest_source_names, output_target)
        atomic_csv(manifest, stage / "output_manifest.csv")
        missing = sorted(set(OUTPUT_NAMES) - {path.name for path in stage.iterdir()})
        unexpected = sorted({path.name for path in stage.iterdir()} - set(OUTPUT_NAMES))
        if missing or unexpected:
            raise RuntimeError(f"Bundle inventory mismatch; missing={missing}, unexpected={unexpected}")
        replace_bundle(stage, output_target)
        validate_promoted_manifest(output_target)
    except Exception:
        if stage.exists():
            shutil.rmtree(stage)
        raise

    printed = summary[[
        "endpoint", "pseudo_f", "r_squared", "p_value_exact", "bh_q_value_two_endpoint_family",
        "extreme_labelings", "enumerated_labelings",
    ]]
    print(printed.to_string(index=False))
    printed_dispersion = dispersion_summary[[
        "endpoint", "dispersion_f", "p_value_exact",
        "bh_q_value_two_endpoint_diagnostic_family",
        "extreme_labelings", "enumerated_labelings",
    ]]
    print(printed_dispersion.to_string(index=False))
    print(f"Wrote {output_target}")


if __name__ == "__main__":
    main()
