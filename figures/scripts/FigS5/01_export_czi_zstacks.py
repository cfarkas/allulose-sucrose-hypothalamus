#!/usr/bin/env python3
"""
Export every source CZI acquisition below a directory as channel-separated,
flat, OME-TIFF Z stacks.

This exporter is deliberately conservative:

* Existing "Maximum intensity projection" derivatives are excluded.
* CZI channels are named by fluor/display metadata (DAPI, AF488, AF546,
  AF633), not by an unverified antibody assignment.
* All TIFFs are written into one output directory with collision-safe names.
* Zeiss mosaic acquisitions are stitched at native resolution, one Z plane
  at a time.
* Source and output SHA-256 hashes, physical scale, dimensions, inferred
  conditions, and provisional marker mappings are recorded in CSV manifests.
* Writes are atomic and resumable. Existing valid outputs are kept unless
  --overwrite is supplied.

The generated ``channel_map.csv`` intentionally marks biological marker
assignments as PROVISIONAL. Confirm the staining panel in that CSV before
running the c-FOS/POMC analysis.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import tifffile
from aicspylibczi import CziFile


SCRIPT_VERSION = "1.1.0-figure5-16h"
PROJECTION_RE = re.compile(r"maximum[\s_-]*intensity[\s_-]*projection", re.I)
CONDITION_ALIASES = {
    "agua": "Water",
    "water": "Water",
    "h2o": "Water",
    "sacarosa": "Sucrose",
    "sucrose": "Sucrose",
    "alulosa": "Allulose",
    "allulose": "Allulose",
    "psicose": "Allulose",
}
ALLOWED_CONDITIONS = {"Water", "Sucrose", "Allulose", "Control", "REVIEW"}

# Exact historical condition supplied by the user in the prior Apotome
# samplesheet. It is retained as an auditable inference, not silently treated
# as filename-confirmed metadata.
HISTORICAL_CONDITION = {
    "FR3-3": "Water",
}


@dataclass
class SourceRecord:
    source: Path
    source_relpath: str
    date: str
    raw_label: str
    sample_id: str
    animal_id: str
    genotype: str
    fasting_hours: int
    prior_sugar_exposure: str
    acquisition_id: str
    condition_raw: str
    condition_filename: str
    condition_canonical: str
    condition_source: str
    condition_status: str
    include_main_analysis: bool
    exclusion_reason: str
    assay_hint: str
    dims: str
    c_size: int
    z_size: int
    m_size: int
    y_size: int
    x_size: int
    pixel_type: str
    pixel_size_x_um: float
    pixel_size_y_um: float
    z_step_um: float
    fluor_labels: List[str]
    channel_metadata: List[Dict[str, str]]
    is_mosaic: bool
    source_size_bytes: int
    source_mtime_ns: int
    source_sha256: str = ""


def safe_token(value: Any, *, keep_hyphen: bool = True) -> str:
    pattern = r"[^A-Za-z0-9_-]+" if keep_hyphen else r"[^A-Za-z0-9_]+"
    out = re.sub(pattern, "_", str(value).strip())
    out = re.sub(r"_+", "_", out).strip("_")
    return out or "unnamed"


def canonical_key(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def canonical_date(folder_name: str) -> str:
    match = re.fullmatch(r"(\d{4})[._-](\d{2})[._-](\d{2})", folder_name)
    if not match:
        return safe_token(folder_name)
    return "-".join(match.groups())


def clean_raw_label(stem: str) -> str:
    text = PROJECTION_RE.sub("", stem)
    text = re.sub(r"\s*[-_ ]*0*001\s*$", "", text, flags=re.I)
    text = re.sub(r"\s*25x\s*$", "", text, flags=re.I)
    text = re.sub(r"\s*[-_ ]*0*001\s*$", "", text, flags=re.I)
    text = re.sub(r"\s+", " ", text).strip(" -_")
    return text


def infer_condition(raw_label: str, animal_id: str) -> Tuple[str, str, str, str]:
    tokens = re.findall(r"[A-Za-z0-9]+", raw_label.lower())
    explicit: List[Tuple[str, str]] = []
    for token in tokens:
        if token in CONDITION_ALIASES:
            explicit.append((token, CONDITION_ALIASES[token]))
    unique = sorted({x[1] for x in explicit})
    if len(unique) == 1:
        return explicit[0][0], unique[0], "filename", "CONFIRMED_FILENAME"
    if len(unique) > 1:
        return ",".join(x[0] for x in explicit), "REVIEW", "filename", "REVIEW_MULTIPLE_CONDITIONS"
    if animal_id in HISTORICAL_CONDITION:
        return "", HISTORICAL_CONDITION[animal_id], "historical_user_samplesheet", "INFERRED_HISTORICAL_REVIEW"
    return "", "REVIEW", "none", "REVIEW_MISSING_CONDITION"


def parse_identity(path: Path, root: Path) -> Dict[str, str]:
    label = clean_raw_label(path.stem)
    # These final words describe a separately acquired color channel, not the
    # biological sample. Removing them groups the three FR4-3 Water CZIs.
    grouped_label = re.sub(r"\s+(AZUL|VERDE|ROJO)\s*$", "", label, flags=re.I).strip()
    sanitized = re.fullmatch(
        r"(?P<sample>FR\d+(?:-\d+)?|AGUA-M|NPY-M)_"
        r"(?P<genotype>WT|NPY-Tg)_"
        r"(?P<condition>Water|Sucrose|Allulose)_16h_NoPriorSugar",
        grouped_label,
        flags=re.I,
    )
    if sanitized:
        sample_id = sanitized.group("sample").upper()
        genotype = "WT" if sanitized.group("genotype").lower() == "wt" else "NPY-Tg"
        condition = sanitized.group("condition").title()
        date = canonical_date(path.parent.name)
        return {
            "date": date,
            "raw_label": label,
            "grouped_label": grouped_label,
            "sample_id": sample_id,
            "animal_id": sample_id,
            "genotype": genotype,
            "fasting_hours": "16",
            "prior_sugar_exposure": "No",
            "condition_raw": condition.lower(),
            "condition_canonical": condition,
            "condition_source": "sanitized_filename_and_source_readme",
            "condition_status": "CONFIRMED_FILENAME",
            "acquisition_id": safe_token(f"{date}__{sample_id}__{condition}"),
            "assay_hint": "cFOS_POMC",
        }

    animal_match = re.search(r"\bFR\s*-?\s*(\d+)(?:\s*-\s*(\d+))?\b", grouped_label, flags=re.I)
    if animal_match:
        animal_id = f"FR{animal_match.group(1)}"
        if animal_match.group(2):
            animal_id += f"-{animal_match.group(2)}"
    else:
        animal_id = safe_token(grouped_label)

    prefix_match = re.match(r"^(E\d+)\b", grouped_label, flags=re.I)
    if prefix_match and animal_match:
        sample_id = f"{prefix_match.group(1).upper()}_{animal_id}"
    elif animal_match:
        sample_id = animal_id
    else:
        sample_id = safe_token(grouped_label)

    date = canonical_date(path.parent.name)
    condition_raw, condition, source, status = infer_condition(grouped_label, animal_id)
    acquisition_id = safe_token(f"{date}__{sample_id}__{condition}")
    assay_hint = "NPY" if re.search(r"\bNPY\b", grouped_label, re.I) else (
        "POMC" if re.search(r"\bPO[MN]C\b", grouped_label, re.I) else "cFOS_POMC"
    )
    return {
        "date": date,
        "raw_label": label,
        "grouped_label": grouped_label,
        "sample_id": sample_id,
        "animal_id": animal_id,
        "genotype": "",
        "fasting_hours": "",
        "prior_sugar_exposure": "",
        "condition_raw": condition_raw,
        "condition_canonical": condition,
        "condition_source": source,
        "condition_status": status,
        "acquisition_id": acquisition_id,
        "assay_hint": assay_hint,
    }


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def first_descendant_text(element: Any, name: str) -> str:
    for child in element.iter():
        if local_name(child.tag) == name and child.text:
            return str(child.text).strip()
    return ""


def channel_metadata(czi: CziFile, c_size: int) -> List[Dict[str, str]]:
    by_id: Dict[str, Dict[str, str]] = {}
    order: List[str] = []
    for element in czi.meta.iter():
        if local_name(element.tag) != "Channel":
            continue
        channel_id = str(element.get("Id", ""))
        if not channel_id:
            continue
        if channel_id not in by_id:
            by_id[channel_id] = {
                "id": channel_id,
                "name": str(element.get("Name", "")),
                "fluor": "",
                "color": "",
                "excitation_nm": "",
                "emission_nm": "",
            }
            order.append(channel_id)
        row = by_id[channel_id]
        for field, xml_name in [
            ("fluor", "Fluor"),
            ("color", "Color"),
            ("excitation_nm", "ExcitationWavelength"),
            ("emission_nm", "EmissionWavelength"),
        ]:
            value = first_descendant_text(element, xml_name)
            if value:
                row[field] = value
    rows = [by_id[key] for key in order]
    # Metadata normally repeats each channel once for acquisition and once for
    # display. The ID merge above preserves acquisition order.
    if len(rows) < c_size:
        for index in range(len(rows), c_size):
            rows.append({
                "id": f"index:{index}", "name": f"C{index}", "fluor": "",
                "color": "", "excitation_nm": "", "emission_nm": "",
            })
    return rows[:c_size]


def fluor_label(meta: Dict[str, str], index: int) -> str:
    fluor = meta.get("fluor", "").lower()
    name = meta.get("name", "").lower()
    color = meta.get("color", "").upper()
    try:
        excitation = float(meta.get("excitation_nm", "") or "nan")
    except Exception:
        excitation = float("nan")
    if "dapi" in fluor or color == "#0000FF" or (math.isfinite(excitation) and excitation < 450):
        return "DAPI"
    if "488" in fluor or color == "#00FF00" or (math.isfinite(excitation) and 450 <= excitation < 530):
        return "AF488"
    if "546" in fluor or color == "#FF0000" or (math.isfinite(excitation) and 530 <= excitation < 620):
        return "AF546"
    if "633" in fluor or color == "#FFFFFF" or (math.isfinite(excitation) and excitation >= 620):
        return "AF633"
    return f"C{index:02d}"


def scaling_um(czi: CziFile) -> Dict[str, float]:
    out = {"X": float("nan"), "Y": float("nan"), "Z": float("nan")}
    for element in czi.meta.iter():
        if local_name(element.tag) != "Distance":
            continue
        axis = str(element.get("Id", "")).upper()
        if axis not in out:
            continue
        raw = first_descendant_text(element, "Value")
        try:
            out[axis] = float(raw) * 1_000_000.0
        except Exception:
            pass
    return out


def dim_size(dims_shape: Dict[str, Tuple[int, int]], dim: str, default: int = 1) -> int:
    value = dims_shape.get(dim)
    return int(value[1]) if value else int(default)


def sha256_file(path: Path, block_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def provisional_marker(fluor: str, available_fluors: Sequence[str]) -> Tuple[str, str]:
    available = set(available_fluors)
    if fluor == "DAPI":
        return "DAPI", "metadata_DAPI"
    if fluor == "AF546":
        return "POMC", "panel_consistency_inference"
    if fluor == "AF633":
        return "cFOS", "panel_consistency_inference"
    if fluor == "AF488":
        if "AF633" in available:
            return "NPY", "four_channel_panel_inference"
        return "cFOS", "three_channel_panel_inference"
    return "REVIEW", "unresolved_fluor"


def read_source_record(path: Path, root: Path) -> SourceRecord:
    czi = CziFile(path)
    dims_shape = czi.get_dims_shape()[0]
    c_size = dim_size(dims_shape, "C")
    z_size = dim_size(dims_shape, "Z")
    m_size = dim_size(dims_shape, "M")
    y_size = dim_size(dims_shape, "Y")
    x_size = dim_size(dims_shape, "X")
    meta = channel_metadata(czi, c_size)
    fluors = [fluor_label(row, i) for i, row in enumerate(meta)]
    scale = scaling_um(czi)
    ident = parse_identity(path, root)
    stat = path.stat()
    if ident["genotype"]:
        include = (
            ident["condition_status"] == "CONFIRMED_FILENAME"
            and ident["genotype"] == "WT"
        )
    else:
        include = (
            ident["condition_status"] == "CONFIRMED_FILENAME"
            and bool(re.match(r"^FR\d", ident["animal_id"]))
        )
    reason = ""
    if not include:
        if ident["genotype"] and ident["genotype"] != "WT":
            reason = "NON_WT_COHORT_REQUIRES_SEPARATE_ANALYSIS"
        elif not re.match(r"^FR\d", ident["animal_id"]):
            reason = "REVIEW_NO_UNAMBIGUOUS_ANIMAL_ID"
        elif ident["condition_status"] != "CONFIRMED_FILENAME":
            reason = ident["condition_status"]
    return SourceRecord(
        source=path.resolve(),
        source_relpath=str(path.resolve().relative_to(root.resolve())),
        date=ident["date"],
        raw_label=ident["raw_label"],
        sample_id=ident["sample_id"],
        animal_id=ident["animal_id"],
        genotype=ident["genotype"],
        fasting_hours=int(ident["fasting_hours"] or 0),
        prior_sugar_exposure=ident["prior_sugar_exposure"],
        acquisition_id=ident["acquisition_id"],
        condition_raw=ident["condition_raw"],
        condition_filename=ident["condition_canonical"],
        condition_canonical=ident["condition_canonical"],
        condition_source=ident["condition_source"],
        condition_status=ident["condition_status"],
        include_main_analysis=include,
        exclusion_reason=reason,
        assay_hint=ident["assay_hint"],
        dims=str(czi.dims),
        c_size=c_size,
        z_size=z_size,
        m_size=m_size,
        y_size=y_size,
        x_size=x_size,
        pixel_type=str(czi.pixel_type),
        pixel_size_x_um=scale["X"],
        pixel_size_y_um=scale["Y"],
        z_step_um=scale["Z"],
        fluor_labels=fluors,
        channel_metadata=meta,
        is_mosaic=m_size > 1 or "M" in str(czi.dims),
        source_size_bytes=int(stat.st_size),
        source_mtime_ns=int(stat.st_mtime_ns),
    )


def mark_condition_conflicts(records: List[SourceRecord]) -> None:
    by_animal: Dict[str, set[str]] = {}
    for record in records:
        if re.match(r"^FR\d", record.animal_id) and record.condition_filename != "REVIEW":
            by_animal.setdefault(record.animal_id, set()).add(record.condition_filename)
    for record in records:
        values = by_animal.get(record.animal_id, set())
        if len(values) > 1:
            record.condition_canonical = "REVIEW"
            record.condition_status = "REVIEW_CONDITION_CONFLICT"
            record.include_main_analysis = False
            record.exclusion_reason = "CONDITION_CONFLICT:" + "|".join(sorted(values))
            record.acquisition_id = safe_token(
                f"{record.date}__{record.sample_id}__{record.condition_filename}"
            )


def discover_sources(root: Path, folders: Optional[Sequence[str]] = None) -> List[Path]:
    normalized_root = root.resolve()
    if not folders:
        search_roots = [normalized_root]
    else:
        search_roots: List[Path] = []
        seen: set[Path] = set()
        for raw_folder in folders:
            candidate_tokens = [token.strip() for token in str(raw_folder).split(",") if token.strip()]
            for token in candidate_tokens:
                folder = Path(token)
                if not folder.is_absolute():
                    folder = (normalized_root / folder).resolve()
                else:
                    folder = folder.resolve()
                if not folder.exists():
                    raise SystemExit(f"Requested folder does not exist: {token}")
                if not folder.is_dir():
                    raise SystemExit(f"Requested folder is not a directory: {token}")
                if folder != normalized_root and normalized_root not in folder.parents:
                    raise SystemExit(f"Requested folder is outside --input root: {token}")
                if folder not in seen:
                    search_roots.append(folder)
                    seen.add(folder)

    discovered: set[Path] = set()
    for base in search_roots:
        for path in base.rglob("*.czi"):
            if path.is_file() and not PROJECTION_RE.search(path.name):
                discovered.add(path)
    return sorted(discovered)


def read_channel_stack(record: SourceRecord, channel: int, cores: int) -> np.ndarray:
    czi = CziFile(record.source)
    if record.is_mosaic:
        first = np.squeeze(czi.read_mosaic(C=channel, Z=0, scale_factor=1.0))
        if first.ndim != 2:
            raise RuntimeError(f"Unexpected mosaic plane shape: {first.shape}")
        stack = np.empty((record.z_size, *first.shape), dtype=first.dtype)
        stack[0] = first
        for z_index in range(1, record.z_size):
            plane = np.squeeze(czi.read_mosaic(C=channel, Z=z_index, scale_factor=1.0))
            if plane.shape != first.shape:
                raise RuntimeError(
                    f"Mosaic plane shape changed at Z={z_index}: {plane.shape} vs {first.shape}"
                )
            stack[z_index] = plane
        return stack

    array, _dims = czi.read_image(C=channel, cores=max(1, int(cores)))
    stack = np.squeeze(np.asarray(array))
    if stack.ndim == 2:
        stack = stack[np.newaxis, ...]
    if stack.ndim != 3:
        raise RuntimeError(f"Expected ZYX after selecting C={channel}; got {stack.shape}")
    if int(stack.shape[0]) != int(record.z_size):
        raise RuntimeError(
            f"Unexpected Z size after selecting C={channel}: {stack.shape[0]} vs {record.z_size}"
        )
    return stack


def write_ome_stack(
    path: Path,
    stack: np.ndarray,
    record: SourceRecord,
    fluor: str,
    compression_level: int,
) -> None:
    metadata: Dict[str, Any] = {
        "axes": "ZYX",
        "Name": f"{record.acquisition_id} {fluor} Z stack",
        "PhysicalSizeX": float(record.pixel_size_x_um),
        "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": float(record.pixel_size_y_um),
        "PhysicalSizeYUnit": "µm",
        "Description": json.dumps({
            "exporter": f"export_czi_zstacks.py {SCRIPT_VERSION}",
            "source_relpath": record.source_relpath,
            "acquisition_id": record.acquisition_id,
            "fluor": fluor,
            "condition": record.condition_canonical,
            "condition_status": record.condition_status,
            "is_mosaic": record.is_mosaic,
        }, sort_keys=True),
    }
    if math.isfinite(record.z_step_um) and record.z_step_um > 0:
        metadata["PhysicalSizeZ"] = float(record.z_step_um)
        metadata["PhysicalSizeZUnit"] = "µm"
    tifffile.imwrite(
        str(path),
        np.asarray(stack),
        ome=True,
        bigtiff=bool(stack.nbytes >= 3_500_000_000),
        photometric="minisblack",
        compression="zlib",
        compressionargs={"level": int(compression_level)},
        predictor=True,
        metadata=metadata,
    )


def validate_output(path: Path, expected_shape: Tuple[int, int, int], expected_dtype: np.dtype) -> None:
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        if tuple(series.shape) != tuple(expected_shape):
            raise RuntimeError(f"Output shape {series.shape} != expected {expected_shape}")
        if np.dtype(series.dtype) != np.dtype(expected_dtype):
            raise RuntimeError(f"Output dtype {series.dtype} != expected {expected_dtype}")
        if not tif.is_ome:
            raise RuntimeError("Output is not OME-TIFF")


def csv_value(value: Any) -> Any:
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, (list, dict)):
        return json.dumps(value, sort_keys=True)
    return value


def write_csv(path: Path, rows: List[Dict[str, Any]], preferred: Optional[Sequence[str]] = None) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    keys: List[str] = []
    for key in preferred or []:
        if any(key in row for row in rows) and key not in keys:
            keys.append(key)
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    tmp = path.with_suffix(path.suffix + ".part")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: csv_value(row.get(key, "")) for key in keys})
    os.replace(tmp, path)


def group_sample_rows(records: Sequence[SourceRecord]) -> List[Dict[str, Any]]:
    grouped: Dict[str, List[SourceRecord]] = {}
    for record in records:
        grouped.setdefault(record.acquisition_id, []).append(record)
    rows: List[Dict[str, Any]] = []
    for acquisition_id, items in sorted(grouped.items()):
        first = items[0]
        fluors: List[str] = []
        sources: List[str] = []
        for item in items:
            sources.append(item.source_relpath)
            for fluor in item.fluor_labels:
                if fluor not in fluors:
                    fluors.append(fluor)
        marker_map = {
            fluor: provisional_marker(fluor, fluors)[0]
            for fluor in fluors
        }
        rows.append({
            "acquisition_id": acquisition_id,
            "sample_id": first.sample_id,
            "animal_id": first.animal_id,
            "genotype": first.genotype,
            "fasting_hours": first.fasting_hours,
            "prior_sugar_exposure": first.prior_sugar_exposure,
            "date": first.date,
            "condition_raw": first.condition_raw,
            "condition_filename": first.condition_filename,
            "condition": first.condition_canonical,
            "condition_source": first.condition_source,
            "condition_status": first.condition_status,
            "include_main_analysis": first.include_main_analysis,
            "exclusion_reason": first.exclusion_reason,
            "assay_hint": first.assay_hint,
            "source_czi_count": len(items),
            "source_relpaths": sources,
            "fluors": fluors,
            "provisional_marker_map": marker_map,
            "channel_map_status": "PROVISIONAL_NEEDS_CONFIRMATION",
            "pixel_size_x_um": first.pixel_size_x_um,
            "pixel_size_y_um": first.pixel_size_y_um,
            "z_steps_um": sorted({round(x.z_step_um, 9) for x in items if math.isfinite(x.z_step_um)}),
            "z_planes": sorted({x.z_size for x in items}),
            "registration_required": len(items) > 1,
            "registration_status": "REVIEW" if len(items) > 1 else "NOT_REQUIRED",
            "notes": (
                "Split-channel sources require registration/crop approval before colocalization."
                if len(items) > 1 else ""
            ),
        })
    return rows


def build_channel_map_rows(records: Sequence[SourceRecord], export_rows: Sequence[Dict[str, Any]]) -> List[Dict[str, Any]]:
    output_by_key = {
        (str(row.get("source_relpath")), int(row.get("channel_index", -1))): row
        for row in export_rows
    }
    grouped_fluors: Dict[str, List[str]] = {}
    for record in records:
        grouped_fluors.setdefault(record.acquisition_id, [])
        for fluor in record.fluor_labels:
            if fluor not in grouped_fluors[record.acquisition_id]:
                grouped_fluors[record.acquisition_id].append(fluor)
    rows: List[Dict[str, Any]] = []
    for record in records:
        available = grouped_fluors[record.acquisition_id]
        for index, fluor in enumerate(record.fluor_labels):
            marker, rule = provisional_marker(fluor, available)
            export = output_by_key.get((record.source_relpath, index), {})
            rows.append({
                "acquisition_id": record.acquisition_id,
                "sample_id": record.sample_id,
                "animal_id": record.animal_id,
                "genotype": record.genotype,
                "fasting_hours": record.fasting_hours,
                "prior_sugar_exposure": record.prior_sugar_exposure,
                "condition": record.condition_canonical,
                "source_relpath": record.source_relpath,
                "channel_index": index,
                "fluor": fluor,
                "marker": marker,
                "mapping_rule": rule,
                "channel_map_status": "PROVISIONAL_NEEDS_CONFIRMATION",
                "output_tiff": export.get("output_tiff", ""),
                "include_channel": marker in {"DAPI", "cFOS", "POMC"},
                "review_notes": "Confirm antibody/fluor assignment from staining notebook before analysis.",
            })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export source CZI files to a flat folder of channel-separated OME-TIFF Z stacks.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=Path("."), help="Root containing dated CZI subfolders.")
    parser.add_argument(
        "--folders",
        nargs="+",
        default=None,
        help="Optional list of folder names/paths under --input to process (e.g. \"2026.07.10\" \"2026.07.09\").",
    )
    parser.add_argument(
        "--outdir",
        "--output",
        type=Path,
        default=Path("zstack_tiffs_flat"),
        help="Single flat output folder.",
    )
    parser.add_argument("--cores", type=int, default=4, help="libCZI read threads for non-mosaic images.")
    parser.add_argument("--compression-level", type=int, default=4, choices=range(0, 10))
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-source-hash", action="store_true", help="Skip source SHA-256 hashing.")
    parser.add_argument("--no-output-hash", action="store_true", help="Skip output SHA-256 hashing.")
    parser.add_argument("--inventory-only", action="store_true", help="Write source/sample manifests without TIFF export.")
    args = parser.parse_args()

    root = args.input.expanduser().resolve()
    outdir = args.outdir.expanduser().resolve()
    if not root.is_dir():
        raise SystemExit(f"Input directory does not exist: {root}")
    outdir.mkdir(parents=True, exist_ok=True)
    folder_inputs: List[str] = []
    if args.folders:
        for value in args.folders:
            folder_inputs.extend(token.strip() for token in str(value).split(",") if token.strip())
    sources = discover_sources(root, folder_inputs)
    if not sources:
        raise SystemExit(f"No source CZI files found below {root}")

    print(f"[INFO] Reading metadata for {len(sources)} source CZI file(s)", flush=True)
    records = [read_source_record(path, root) for path in sources]
    mark_condition_conflicts(records)
    if not args.no_source_hash:
        for number, record in enumerate(records, start=1):
            print(f"[HASH source {number}/{len(records)}] {record.source_relpath}", flush=True)
            record.source_sha256 = sha256_file(record.source)

    export_rows: List[Dict[str, Any]] = []
    total_channels = sum(record.c_size for record in records)
    done = 0
    for record in records:
        for channel_index, fluor in enumerate(record.fluor_labels):
            done += 1
            filename = safe_token(
                f"{record.acquisition_id}__{fluor}__zstack", keep_hyphen=True
            ) + ".ome.tif"
            output = outdir / filename
            base_row: Dict[str, Any] = {
                "exporter_version": SCRIPT_VERSION,
                "source_relpath": record.source_relpath,
                "source_path": str(record.source),
                "source_size_bytes": record.source_size_bytes,
                "source_mtime_ns": record.source_mtime_ns,
                "source_sha256": record.source_sha256,
                "date": record.date,
                "raw_label": record.raw_label,
                "acquisition_id": record.acquisition_id,
                "sample_id": record.sample_id,
                "animal_id": record.animal_id,
                "condition_raw": record.condition_raw,
                "condition_filename": record.condition_filename,
                "condition": record.condition_canonical,
                "condition_source": record.condition_source,
                "condition_status": record.condition_status,
                "include_main_analysis": record.include_main_analysis,
                "exclusion_reason": record.exclusion_reason,
                "assay_hint": record.assay_hint,
                "channel_index": channel_index,
                "fluor": fluor,
                "channel_metadata": record.channel_metadata[channel_index],
                "dims": record.dims,
                "is_mosaic": record.is_mosaic,
                "c_size": record.c_size,
                "z_size": record.z_size,
                "m_size": record.m_size,
                "source_y_size": record.y_size,
                "source_x_size": record.x_size,
                "pixel_type": record.pixel_type,
                "pixel_size_x_um": record.pixel_size_x_um,
                "pixel_size_y_um": record.pixel_size_y_um,
                "z_step_um": record.z_step_um,
                "z_span_um": max(0, record.z_size - 1) * record.z_step_um,
                "output_tiff": filename,
                "output_path": str(output),
                "output_sha256": "",
                "output_size_bytes": 0,
                "output_shape_zyx": "",
                "status": "inventory_only" if args.inventory_only else "pending",
                "error": "",
            }
            if args.inventory_only:
                export_rows.append(base_row)
                continue
            print(
                f"[EXPORT {done}/{total_channels}] {record.source_relpath} "
                f"C={channel_index} {fluor} -> {filename}",
                flush=True,
            )
            try:
                if output.exists() and not args.overwrite:
                    with tifffile.TiffFile(output) as tif:
                        shape = tuple(int(x) for x in tif.series[0].shape)
                        dtype = np.dtype(tif.series[0].dtype)
                        is_ome = bool(tif.is_ome)
                    if len(shape) == 3 and shape[0] == record.z_size and dtype == np.dtype("uint16") and is_ome:
                        base_row["status"] = "existing_valid_kept"
                        base_row["output_shape_zyx"] = "x".join(map(str, shape))
                        base_row["output_size_bytes"] = output.stat().st_size
                        if not args.no_output_hash:
                            base_row["output_sha256"] = sha256_file(output)
                        export_rows.append(base_row)
                        continue
                    raise RuntimeError(
                        f"Existing output is invalid; use --overwrite after review: {output}"
                    )
                stack = read_channel_stack(record, channel_index, args.cores)
                temp = output.with_suffix(output.suffix + ".part")
                if temp.exists():
                    temp.unlink()
                write_ome_stack(temp, stack, record, fluor, args.compression_level)
                validate_output(temp, tuple(int(x) for x in stack.shape), stack.dtype)
                os.replace(temp, output)
                base_row["status"] = "exported"
                base_row["output_shape_zyx"] = "x".join(map(str, stack.shape))
                base_row["output_size_bytes"] = output.stat().st_size
                if not args.no_output_hash:
                    base_row["output_sha256"] = sha256_file(output)
                del stack
            except Exception as exc:
                base_row["status"] = "failed"
                base_row["error"] = f"{type(exc).__name__}: {exc}"
                print(f"[ERROR] {base_row['error']}", file=sys.stderr, flush=True)
            export_rows.append(base_row)
            write_csv(outdir / "zstack_export_manifest.partial.csv", export_rows)

    source_rows: List[Dict[str, Any]] = []
    for record in records:
        row = asdict(record)
        row["source"] = str(record.source)
        source_rows.append(row)
    sample_rows = group_sample_rows(records)
    channel_map_rows = build_channel_map_rows(records, export_rows)
    exclusions = [
        row for row in sample_rows
        if str(row.get("include_main_analysis")).upper() not in {"TRUE", "1"}
        or str(row.get("condition")) == "REVIEW"
        or str(row.get("registration_status")) == "REVIEW"
    ]
    write_csv(outdir / "source_czi_inventory.csv", source_rows)
    write_csv(outdir / "zstack_export_manifest.csv", export_rows)
    write_csv(outdir / "sample_manifest.csv", sample_rows)
    write_csv(outdir / "channel_map.csv", channel_map_rows)
    write_csv(outdir / "exclusions_and_conflicts.csv", exclusions)
    partial = outdir / "zstack_export_manifest.partial.csv"
    if partial.exists():
        partial.unlink()

    failed = sum(row.get("status") == "failed" for row in export_rows)
    exported = sum(row.get("status") in {"exported", "existing_valid_kept"} for row in export_rows)
    print(f"[DONE] {exported}/{len(export_rows)} channel stacks available; failed={failed}", flush=True)
    print(f"[DONE] Output: {outdir}", flush=True)
    print(f"[REVIEW REQUIRED] Confirm {outdir / 'channel_map.csv'} before biological analysis.", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
