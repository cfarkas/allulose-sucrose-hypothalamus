#!/usr/bin/env python3
"""Add the 2026 Zeiss LSM 780 Water/NPY controls to the Figure 3 cell tables.

The nuclear objects, the NPY association rule (>=1 overlapping pixel and >=5%
of the nucleus) and the c-FOS assignment rule (ROI >=5 px whose best-matching
DAPI nucleus covers >=40% of the ROI) are the ones the original nine animals
use.  Nothing about those nine animals is recomputed: their recovered rows are
copied verbatim and their hashes are recorded.  The command writes a new,
self-contained source tree; it never edits the 2026-09-19 reanalysis package.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

PAPER = Path(__file__).resolve().parents[2]
RECOVERED = PAPER / "thesis_and_manuscript/Scientific_Reports_Segura_et_al_2026/spatial_reanalysis_20260919"
DERIVED = PAPER / "Fig3/raw/new_water_control_LSM780_24_SEPT_2026/derived"

CONDITION = "Water"
PIXEL_UM = 0.664212669961432
# Ventricular axis fitted to each acquisition's DAPI support gap, as fractions
# of image width/height. SIDE_A/B are image sides, not registered anatomy.
ADDED = {
    "WATER_NPY4": dict(label="N132-5", side_partition=[[0.47, 0.00], [0.51, 1.00]]),
    "WATER_NPY5": dict(label="N132-3", side_partition=[[0.39, 0.00], [0.45, 1.00]]),
}
TISSUE_SUPPORT_UM = 15.0
NPY_MIN_OVERLAP_PX = 1
NPY_MIN_FRACTION = 0.05
CFOS_MIN_AREA_PX = 5
CFOS_MIN_OVERLAP_FRACTION = 0.40
LEGACY_ZEISS_PIXEL_UM = 0.664212669961432


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def recovered_module():
    """Reuse the delivered recovery code so the rules cannot drift."""
    spec = importlib.util.spec_from_file_location("fig3_extract_cells", RECOVERED / "src/extract_cells.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["fig3_extract_cells"] = module
    sys.path.insert(0, str(RECOVERED / "src"))
    spec.loader.exec_module(module)
    return module


def read_mask(path: Path) -> np.ndarray:
    record = np.load(path, allow_pickle=True).item()
    masks = np.asarray(record["masks"])
    if masks.ndim != 2 or not np.issubdtype(masks.dtype, np.integer) or masks.min() < 0:
        raise ValueError(f"Unusable Cellpose mask: {path}")
    return masks


def read_display(path: Path) -> np.ndarray:
    import tifffile

    array = tifffile.imread(path)
    if array.ndim == 3 and array.shape[-1] in (3, 4):
        array = array[..., :3].max(axis=-1)
    if array.ndim != 2:
        raise ValueError(f"Expected a 2D image: {path} {array.shape}")
    return array


def marker_positive_nuclei(dapi: np.ndarray, marker: np.ndarray) -> set[int]:
    """The original Figure 3 NPY rule, evaluated on the accepted masks."""
    positive = np.asarray(marker > 0, dtype=bool)
    top = int(dapi.max())
    if top == 0:
        return set()
    labels = np.arange(1, top + 1)
    area = np.bincount(dapi.ravel(), minlength=top + 1)[1:]
    overlap = np.bincount(dapi[positive].ravel(), minlength=top + 1)[1:]
    keep = (overlap >= NPY_MIN_OVERLAP_PX) & (overlap / np.maximum(area, 1) >= NPY_MIN_FRACTION) & (area > 0)
    return {int(label) for label in labels[keep]}


def build_animal(module, animal: str, output: Path) -> tuple[pd.DataFrame, dict, dict, list, list]:
    spec = ADDED[animal]
    sample_id = animal + "_S01"
    acquisition_id = animal + "__S01"
    dapi = read_mask(DERIVED / f"{sample_id}_dapi_seg.npy")
    npy = read_mask(DERIVED / f"{sample_id}_GFP_seg.npy")
    fos = read_mask(DERIVED / f"{sample_id}_fos_seg.npy")
    if not (dapi.shape == npy.shape == fos.shape):
        raise ValueError(f"Channel masks disagree in shape for {animal}")
    sx = sy = PIXEL_UM
    height, width = dapi.shape

    legacy = marker_positive_nuclei(dapi, npy)
    cfos, _, cfos_rows = module.classify_cfos(dapi, fos, CFOS_MIN_AREA_PX, CFOS_MIN_OVERLAP_FRACTION)
    nuc = module.props(dapi).rename(columns={"label": "nucleus_label", "centroid-0": "y_px", "centroid-1": "x_px", "area": "area_px"})
    roi = np.ones(dapi.shape, np.uint8)
    (ax, ay), (bx, by) = spec["side_partition"]
    anchor = np.array([ax * width, ay * height])
    direction = np.array([bx * width, by * height]) - anchor
    yy, xx = np.indices(dapi.shape)
    side = (((xx - anchor[0]) * direction[1] - (yy - anchor[1]) * direction[0]) >= 0).astype(np.uint8) + 1
    tissue = ndi.distance_transform_edt(dapi == 0, sampling=(sy, sx)) <= TISSUE_SUPPORT_UM
    ventricle = np.zeros(dapi.shape, bool)
    landmark_source = "DAPI support gaps; no validated ventricular boundary"

    assigns, ledger, details = module.marker_assignments(dapi, npy, np.zeros_like(npy), roi, sx, sy, False)
    labels = nuc.nucleus_label.astype(int)
    nuc["cohort"] = "NPY"
    nuc["marker"] = "NPY"
    nuc["animal_id"] = animal
    nuc["condition"] = CONDITION
    nuc["cage"] = "unknown"
    nuc["genotype"] = "NPY-GFP"
    nuc["sex"] = "unknown"
    nuc["acquisition_id"] = acquisition_id
    nuc["section"] = "S01"
    nuc["region"] = "FIELD"
    nuc["region_code"] = 1
    nuc["cell_id"] = "NPY|" + acquisition_id + "|N" + labels.astype(str)
    nuc["x_um"] = nuc.x_px * sx
    nuc["y_um"] = nuc.y_px * sy
    nuc["pixel_x_um"] = sx
    nuc["pixel_y_um"] = sy
    iy = np.clip(nuc.y_px.round().astype(int), 0, height - 1)
    ix = np.clip(nuc.x_px.round().astype(int), 0, width - 1)
    nuc["hemifield"] = np.where(side[iy, ix] == 1, "SIDE_A", "SIDE_B")
    nuc["hemisphere_evidence"] = "DAPI-only computational partition; anatomical side unregistered"
    nuc["in_tissue"] = tissue[iy, ix]
    nuc["is_cfos"] = labels.isin(cfos)
    nuc["accepted_legacy"] = labels.isin(legacy)
    nuc["marker_channel_present"] = True
    nuc["npy_channel_informative"] = True
    nuc["pomc_channel_present"] = False
    nuc["is_npy_legacy"] = labels.isin(legacy)
    nuc["is_pomc_legacy"] = pd.NA
    nuc["source_dapi_mask"] = str(DERIVED / f"{sample_id}_dapi_seg.npy")
    nuc["source_marker_mask"] = str(DERIVED / f"{sample_id}_GFP_seg.npy")
    nuc["source_cfos_mask"] = str(DERIVED / f"{sample_id}_fos_seg.npy")
    nuc["include_acquisition"] = True
    for variant in module.VARIANTS:
        nuc[variant] = labels.isin(assigns[variant])
    nuc["marker_label"] = labels.map({key: value["marker_label"] for key, value in details.items()})
    for key in ("association_overlap_fraction", "association_distance_um"):
        nuc[key] = labels.map({label: value[key] for label, value in details.items()})
    for threshold in (0.3, 0.5):
        nuc[f"cfos_overlap_{threshold}"] = labels.isin(module.classify_cfos(dapi, fos, CFOS_MIN_AREA_PX, threshold)[0])
    nuc["cfos_physical_size"] = labels.isin(
        module.classify_cfos(dapi, fos, CFOS_MIN_AREA_PX * LEGACY_ZEISS_PIXEL_UM**2 / (sx * sy), CFOS_MIN_OVERLAP_FRACTION)[0])
    quality = []
    for channel, suffix in (("dapi", "dapi"), ("npy", "GFP"), ("cfos", "fos")):
        image = read_display(DERIVED / f"{sample_id}_{suffix}.tif")
        nuc[channel + "_nuclear_mean"] = ndi.mean(image, dapi, labels)
        mask = {"dapi": dapi, "npy": npy, "cfos": fos}[channel]
        quality.append(dict(cohort="NPY", animal_id=animal, acquisition_id=acquisition_id, channel=channel,
                            **module.intensity_qc(image, mask)))
    nuc["ventricular_distance_um"] = np.nan
    nuc["arc_me_boundary_distance_um"] = np.nan
    nuc["opposite_region_distance_um"] = np.nan
    nuc["primary_marker_status"] = np.where(~nuc.marker_channel_present, "missing_channel",
                                            np.where(nuc.compact_unique, "eligible_mask_association", "excluded_by_technical_rule"))

    geometry = output / "cache" / f"{acquisition_id}_geometry.npz"
    geometry.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(geometry, roi=roi, tissue=tissue, side=side, ventricle=ventricle)

    manifest_row = dict(
        cohort="NPY", animal_id=animal, condition=CONDITION, acquisition_id=acquisition_id, section="S01",
        cage="unknown", genotype="NPY-GFP", sex="unknown",
        dapi_path=str(DERIVED / f"{sample_id}_dapi.tif"),
        marker_path=str(DERIVED / f"{sample_id}_GFP.tif"),
        cfos_path=str(DERIVED / f"{sample_id}_fos.tif"),
        dapi_mask=str(DERIVED / f"{sample_id}_dapi_seg.npy"),
        marker_mask=str(DERIVED / f"{sample_id}_GFP_seg.npy"),
        cfos_mask=str(DERIVED / f"{sample_id}_fos_seg.npy"),
        geometry_path=str(geometry), pixel_x_um=sx, pixel_y_um=sy, include_acquisition=True,
        acquisition_reason="single accepted field",
        region_source="sampled field; no transferable complete ARC annotation",
        landmark_source=landmark_source,
        source_record=str(DERIVED / f"{sample_id}_dapi_seg.npy"),
    )
    receipt = dict(
        animal_id=animal, native_animal_label=spec["label"], acquisition_id=acquisition_id, condition=CONDITION,
        pixel_um=PIXEL_UM, image_shape_yx=[height, width],
        field_size_um=[round(height * sy, 3), round(width * sx, 3)],
        side_partition_fractions=spec["side_partition"],
        counts=dict(dapi_nuclei=int(len(nuc)), npy_positive_nuclei=int(nuc.accepted_legacy.sum()),
                    cfos_positive_nuclei=int(nuc.is_cfos.sum()),
                    double_positive_nuclei=int((nuc.accepted_legacy & nuc.is_cfos).sum()),
                    primary_marker_nuclei=int(nuc.compact_unique.sum()),
                    npy_objects=int(npy.max()), cfos_objects=int(fos.max()),
                    in_tissue_nuclei=int(nuc.in_tissue.sum()),
                    side_a_nuclei=int((nuc.hemifield == "SIDE_A").sum()),
                    side_b_nuclei=int((nuc.hemifield == "SIDE_B").sum()),
                    tissue_support_fraction=float(tissue.mean())),
        image_quality=quality,
        source_masks={name: sha256(DERIVED / f"{sample_id}_{name}_seg.npy") for name in ("dapi", "GFP", "fos")},
        geometry_sha256=sha256(geometry),
    )
    for row in cfos_rows:
        row.update(cohort="NPY", animal_id=animal, acquisition_id=acquisition_id)
    for row in ledger:
        row.update(cohort="NPY", animal_id=animal, acquisition_id=acquisition_id,
                   source_mask=str(DERIVED / f"{sample_id}_GFP_seg.npy"))
    return nuc, manifest_row, receipt, cfos_rows, ledger


def build(output: Path, animals: list[str], overwrite: bool = False) -> dict:
    module = recovered_module()
    output = output.resolve()
    if output.exists() and not overwrite:
        raise RuntimeError(f"Fresh source tree already exists: {output}")
    if output.exists():
        shutil.rmtree(output)
    for folder in ("data", "audit", "cache", "src"):
        (output / folder).mkdir(parents=True)

    frames, manifest_rows, receipts, cfos_rows, ledger_rows = [], [], {}, [], []
    for animal in animals:
        nuc, row, receipt, rows, ledger = build_animal(module, animal, output)
        frames.append(nuc)
        manifest_rows.append(row)
        receipts[animal] = receipt
        cfos_rows.extend(rows)
        ledger_rows.extend(ledger)
        print(f"{animal}: " + json.dumps(receipt["counts"]), flush=True)

    original_cells = pd.read_csv(RECOVERED / "data/NPY_cells.csv.gz", low_memory=False)
    overlap = set(original_cells.animal_id) & set(animals)
    if overlap:
        raise ValueError(f"The recovered table already contains {sorted(overlap)}")
    added = pd.concat(frames, ignore_index=True)
    missing = [column for column in original_cells.columns if column not in added.columns]
    extra = [column for column in added.columns if column not in original_cells.columns]
    if missing or extra:
        raise ValueError(f"Column mismatch; missing={missing} extra={extra}")
    combined = pd.concat([original_cells, added[original_cells.columns]], ignore_index=True)
    if combined.cell_id.duplicated().any():
        raise ValueError("Nonunique cell identifiers")
    combined.to_csv(output / "data/NPY_cells.csv.gz", index=False, compression="gzip")

    manifest = pd.read_csv(RECOVERED / "audit/NPY_acquisition_manifest.csv")
    manifest = pd.concat([manifest, pd.DataFrame(manifest_rows)[manifest.columns]], ignore_index=True)
    manifest.to_csv(output / "audit/NPY_acquisition_manifest.csv", index=False)

    for source in sorted((RECOVERED / "cache").glob("*_geometry.npz")):
        if source.name.startswith(tuple(f"{animal}__" for animal in manifest.animal_id if animal not in animals)):
            target = output / "cache" / source.name
            if not target.exists():
                shutil.copy2(source, target)
    for name in ("spatial.py", "extract_cells.py"):
        shutil.copy2(RECOVERED / "src" / name, output / "src" / name)
    pd.DataFrame(cfos_rows).to_csv(output / "audit/new_water_cfos_assignment_audit.csv", index=False)
    pd.DataFrame(ledger_rows).to_csv(output / "audit/new_water_marker_object_audit.csv", index=False)

    receipt = dict(
        schema="fig3-new-water-cells-v1", added_animals=animals, per_animal=receipts,
        rules=dict(npy_min_overlap_px=NPY_MIN_OVERLAP_PX, npy_min_fraction=NPY_MIN_FRACTION,
                   cfos_min_roi_area_px=CFOS_MIN_AREA_PX, cfos_dapi_overlap_fraction=CFOS_MIN_OVERLAP_FRACTION,
                   tissue_support_um=TISSUE_SUPPORT_UM),
        recovered_inputs={
            "NPY_cells.csv.gz": sha256(RECOVERED / "data/NPY_cells.csv.gz"),
            "NPY_acquisition_manifest.csv": sha256(RECOVERED / "audit/NPY_acquisition_manifest.csv"),
        },
        original_animal_rows_unchanged=True,
        outputs={
            "NPY_cells.csv.gz": sha256(output / "data/NPY_cells.csv.gz"),
            "NPY_acquisition_manifest.csv": sha256(output / "audit/NPY_acquisition_manifest.csv"),
        },
        script_sha256=sha256(Path(__file__)),
    )
    (output / "new_water_cells_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True, default=float) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, required=True, help="New source tree for the extended cohort")
    parser.add_argument("--animal", action="append", choices=sorted(ADDED), help="Repeatable; default is every 2026 Water animal")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    build(args.output, args.animal or sorted(ADDED), args.force)


if __name__ == "__main__":
    main()
