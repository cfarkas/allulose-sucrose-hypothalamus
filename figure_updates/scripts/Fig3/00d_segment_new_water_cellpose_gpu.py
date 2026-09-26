#!/usr/bin/env python3
"""GPU Cellpose segmentation of the 2026 Zeiss LSM 780 Water/NPY controls.

The models, flow threshold, cell-probability threshold and normalization are
the ones recorded for the original nine Figure 3 animals: ``dapi_channel`` for
DAPI and ``cfos_channel3`` for the NPY/GFP and c-FOS channels.  Inputs are the
legacy-hue RGB16 exports written by ``00c_extract_new_water_channels.py`` at
the legacy Zeiss pixel pitch, so the models see objects at the same physical
scale as in the sugar cohort.

The script never edits the microscopy source.  It saves Cellpose
GUI-compatible ``*_seg.npy`` dictionaries and a hash-bound receipt.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import cellpose
import numpy as np
import tifffile
import torch
from cellpose import io, models

PAPER_ROOT = next(
    path
    for path in Path(__file__).resolve().parents
    if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir() and (path / "README.txt").is_file()))
)
DERIVED = PAPER_ROOT / "Fig3/raw/new_water_control_LSM780_24_SEPT_2026/derived"
ANIMALS = {"WATER_NPY4": "N132-5", "WATER_NPY5": "N132-3"}
ANALYSIS_PIXEL_UM = 0.664212669961432
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0
CHANNEL_MODELS = (("dapi_channel", ("dapi",)), ("cfos_channel3", ("GFP", "fos")))
MODEL_SHA256 = {
    "dapi_channel": "ae49db6b54dd729502d9f7fa03e1ddd1e6e6cfc78e97ed6b9ffdfe9db8799e24",
    "cfos_channel3": "2fb8312633be96f4a5d1bf48fd9a629b9fe3d0d2ca7a38de34b778300dc0549c",
}
NORMALIZE_PARAMS = {
    "lowhigh": None, "percentile": [1.0, 99.0], "normalize": True, "norm3D": True,
    "sharpen_radius": 0.0, "smooth_radius": 0.0, "tile_norm_blocksize": 0.0,
    "tile_norm_smooth3D": 0.0, "invert": False,
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pixel_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def require_gpu() -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; this workflow is GPU-only")
    free_bytes, total_bytes = torch.cuda.mem_get_info()
    if free_bytes < 12 * 1024**3:
        raise RuntimeError(f"Less than 12 GiB GPU memory is free: {free_bytes / 1024**3:.2f} GiB")
    print(f"gpu={torch.cuda.get_device_name(0)} free_gib={free_bytes / 1024**3:.2f} total_gib={total_bytes / 1024**3:.2f}", flush=True)


def read_rgb_input(path: Path) -> np.ndarray:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Missing/non-regular input: {path}")
    image = np.asarray(tifffile.imread(path))
    if image.ndim != 3 or image.shape[-1] != 3 or image.dtype != np.uint16:
        raise RuntimeError(f"Unexpected input {path}: {image.shape}/{image.dtype}")
    return image


def save_gui_compatible(image, masks, flows, source_path: Path, output_path: Path, model_path: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="new_water_seg_save_", dir=output_path.parent) as tmp:
        temporary_image = Path(tmp) / source_path.name
        shutil.copy2(source_path, temporary_image)
        io.masks_flows_to_seg(image, masks, flows, str(temporary_image), channels=None)
        generated = temporary_image.with_name(temporary_image.stem + "_seg.npy")
        data = np.load(generated, allow_pickle=True).item()
        data.update({
            "colors": np.zeros((int(masks.max()), 3), dtype=np.uint8),
            "ismanual": np.zeros(int(masks.max()), dtype=bool),
            "manual_changes": [], "model_path": str(model_path),
            "flow_threshold": FLOW_THRESHOLD, "cellprob_threshold": CELLPROB_THRESHOLD,
            "normalize_params": NORMALIZE_PARAMS, "restore": None, "ratio": 1.0,
            "diameter": None, "filename": str(source_path),
        })
        staged = output_path.with_suffix(output_path.suffix + ".partial")
        with staged.open("wb") as handle:
            np.save(handle, data, allow_pickle=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(staged, output_path)


def segment_one(directory: Path, sample_id: str, suffix: str, model, model_path: Path) -> dict[str, object]:
    source = directory / f"{sample_id}_{suffix}.tif"
    image = read_rgb_input(source)
    result = model.eval(image, channels=None, channel_axis=-1, diameter=None,
                        flow_threshold=FLOW_THRESHOLD, cellprob_threshold=CELLPROB_THRESHOLD,
                        normalize=NORMALIZE_PARAMS)
    masks, flows = np.asarray(result[0]), result[1]
    if masks.shape != image.shape[:2] or masks.dtype.kind not in "ui":
        raise RuntimeError(f"Unexpected {suffix} masks: {masks.shape}/{masks.dtype}")
    output = directory / f"{sample_id}_{suffix}_seg.npy"
    save_gui_compatible(image, masks, flows, source, output, model_path)
    loaded = np.load(output, allow_pickle=True).item()["masks"]
    if not np.array_equal(masks, loaded):
        raise RuntimeError(f"Saved {suffix} mask array differs from inference result")
    return {
        "suffix": suffix, "source": str(source), "source_sha256": sha256_file(source),
        "source_pixel_sha256": pixel_sha256(image.max(axis=-1)),
        "model": str(model_path), "model_sha256": sha256_file(model_path),
        "output": str(output), "output_bytes": output.stat().st_size, "output_sha256": sha256_file(output),
        "shape": list(masks.shape), "dtype": str(masks.dtype), "objects": int(masks.max()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--animal", action="append", choices=sorted(ANIMALS), help="Repeatable; default is every 2026 Water animal")
    parser.add_argument("--derived-dir", type=Path, default=DERIVED)
    parser.add_argument("--model-root", type=Path, default=Path("/home/server/.cellpose/models"))
    parser.add_argument("--force", action="store_true", help="Replace existing masks for the selected animals")
    args = parser.parse_args()
    require_gpu()
    for model_name, _ in CHANNEL_MODELS:
        path = args.model_root / model_name
        if path.is_symlink() or not path.is_file() or sha256_file(path) != MODEL_SHA256[model_name]:
            raise RuntimeError(f"Model missing or hash mismatch: {path}")
    animals = args.animal or sorted(ANIMALS)
    directory = args.derived_dir
    directory.mkdir(parents=True, exist_ok=True)
    for animal in animals:
        for _, suffixes in CHANNEL_MODELS:
            for suffix in suffixes:
                existing = directory / f"{animal}_S01_{suffix}_seg.npy"
                if existing.exists() and not args.force:
                    raise RuntimeError(f"Fresh output already exists: {existing}")
    records: list[dict[str, object]] = []
    for model_name, suffixes in CHANNEL_MODELS:
        model_path = args.model_root / model_name
        model = models.CellposeModel(gpu=True, pretrained_model=str(model_path))
        for animal in animals:
            for suffix in suffixes:
                print(f"{animal}: segmenting={suffix} model={model_name}", flush=True)
                record = segment_one(directory, animal + "_S01", suffix, model, model_path)
                record["animal_id"] = animal
                records.append(record)
                print(f"  objects={record['objects']}", flush=True)
        del model
        gc.collect()
        torch.cuda.empty_cache()
    receipt = {
        "schema_version": "figure3-new-water-cellpose-gpu-v1",
        "animals": {animal: ANIMALS[animal] for animal in animals},
        "cellpose_version": getattr(cellpose, "version", getattr(cellpose, "__version__", "unknown")),
        "torch_version": torch.__version__, "cuda_device": torch.cuda.get_device_name(0),
        "analysis_pixel_um": ANALYSIS_PIXEL_UM,
        "parameters": {"diameter": None, "flow_threshold": FLOW_THRESHOLD,
                       "cellprob_threshold": CELLPROB_THRESHOLD, "channels": None,
                       "channel_axis": -1, "normalization": NORMALIZE_PARAMS},
        "model_assignment_matches_original_nine": {"DAPI": "dapi_channel", "NPY": "cfos_channel3", "c-FOS": "cfos_channel3"},
        "records": records,
    }
    (directory / "NEW_WATER_CELLPOSE_GPU_RECEIPT.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    for record in records:
        print(f"{record['animal_id']} {record['suffix']}: {record['objects']} objects", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise
