#!/usr/bin/env python3
"""GPU-only, fresh-output Cellpose reconstruction for the three Water3 masks.

This script never edits the microscopy source or an existing output.  It saves
Cellpose GUI-compatible ``*_seg.npy`` dictionaries and a hash-bound receipt.
Frozen pre-incident file hashes and object counts are used as comparison gates;
a regenerated mismatch is retained but explicitly marked pending HIL review.
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
from cellpose import io, models, utils


EXPECTED_SHAPE = (1419, 1898)
FLOW_THRESHOLD = 0.4
CELLPROB_THRESHOLD = 0.0
MODEL_SHA256 = {
    "dapi_Keyence4": "d973df68cf36f226f4df97c106d82fbeedfc1ea5a0ddd427a5a8a00f0657d004",
    "cfos_channel3": "2fb8312633be96f4a5d1bf48fd9a629b9fe3d0d2ca7a38de34b778300dc0549c",
}
INPUT_PIXEL_SHA256 = {
    "dapi": "c9bb43c80f07a019d24e557922d6a78fb50d3f627b54fc217aa6494cd72d239c",
    "GFP": "28fd400469901436b0a5be093d5fb4afe49c56f4905951d51e9e8b4ce23cf10e",
    "fos": "33bc1f0fb9f625d04cef15bb711da44555a554a2c5e4df5a55d4558c08f10aef",
}
FROZEN_SEG = {
    "dapi": {
        "model": "dapi_Keyence4",
        "objects": 3128,
        "bytes": 53_884_707,
        "sha256": "58be1372ada83b0f799c27a32037e95e8ee28b4eeaaba92825b545fc91abad87",
    },
    "GFP": {
        "model": "cfos_channel3",
        "objects": 188,
        "bytes": 53_947_212,
        "sha256": "f662cf3ce31f564c4f760d4617a614107916f53d4ae9fc620114480baec50c94",
    },
    "fos": {
        "model": "cfos_channel3",
        "objects": 53,
        "bytes": 53_872_584,
        "sha256": "91761029120ce769c213fe9e305f5e12cebf43d5a4d7ea4af0ef7a608919c061",
    },
}
NORMALIZE_PARAMS = {
    "lowhigh": None,
    "percentile": [1.0, 99.0],
    "normalize": True,
    "norm3D": True,
    "sharpen_radius": 0.0,
    "smooth_radius": 0.0,
    "tile_norm_blocksize": 0.0,
    "tile_norm_smooth3D": 0.0,
    "invert": False,
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
        raise RuntimeError(
            f"Less than 12 GiB GPU memory is free: {free_bytes / 1024**3:.2f} GiB"
        )
    print(
        f"gpu={torch.cuda.get_device_name(0)} "
        f"free_gib={free_bytes / 1024**3:.2f} total_gib={total_bytes / 1024**3:.2f}",
        flush=True,
    )


def read_rgb_input(path: Path, suffix: str) -> np.ndarray:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Missing/non-regular input: {path}")
    image = np.asarray(tifffile.imread(path))
    if image.shape != (*EXPECTED_SHAPE, 3) or image.dtype != np.uint16:
        raise RuntimeError(f"Unexpected input {path}: {image.shape}/{image.dtype}")
    signal = image.max(axis=-1)
    observed = pixel_sha256(signal)
    if observed != INPUT_PIXEL_SHA256[suffix]:
        raise RuntimeError(
            f"Frozen {suffix} pixel hash mismatch: {observed} != {INPUT_PIXEL_SHA256[suffix]}"
        )
    return image


def save_gui_compatible(
    image: np.ndarray,
    masks: np.ndarray,
    flows: list[np.ndarray],
    source_path: Path,
    output_path: Path,
    model_path: Path,
) -> None:
    with tempfile.TemporaryDirectory(prefix="water3_seg_save_", dir=output_path.parent) as tmp:
        temporary_image = Path(tmp) / source_path.name
        shutil.copy2(source_path, temporary_image)
        io.masks_flows_to_seg(
            image,
            masks,
            flows,
            str(temporary_image),
            channels=None,
        )
        generated = temporary_image.with_name(temporary_image.stem + "_seg.npy")
        data = np.load(generated, allow_pickle=True).item()
        data.update(
            {
                "colors": np.zeros((int(masks.max()), 3), dtype=np.uint8),
                "ismanual": np.zeros(int(masks.max()), dtype=bool),
                "manual_changes": [],
                "model_path": str(model_path),
                "flow_threshold": FLOW_THRESHOLD,
                "cellprob_threshold": CELLPROB_THRESHOLD,
                "normalize_params": NORMALIZE_PARAMS,
                "restore": None,
                "ratio": 1.0,
                "diameter": None,
                "filename": str(source_path),
            }
        )
        staged = output_path.with_suffix(output_path.suffix + ".partial")
        with staged.open("wb") as handle:
            np.save(handle, data, allow_pickle=True)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(staged, output_path)


def segment_one(
    input_dir: Path,
    output_dir: Path,
    suffix: str,
    model: models.CellposeModel,
    model_path: Path,
) -> dict[str, object]:
    source = input_dir / f"WATER_NPY3_S01_{suffix}.tif"
    image = read_rgb_input(source, suffix)
    result = model.eval(
        image,
        channels=None,
        channel_axis=-1,
        diameter=None,
        flow_threshold=FLOW_THRESHOLD,
        cellprob_threshold=CELLPROB_THRESHOLD,
        normalize=NORMALIZE_PARAMS,
    )
    masks, flows = np.asarray(result[0]), result[1]
    if masks.shape != EXPECTED_SHAPE or masks.dtype.kind not in "ui":
        raise RuntimeError(f"Unexpected {suffix} masks: {masks.shape}/{masks.dtype}")
    output = output_dir / f"WATER_NPY3_S01_{suffix}_seg.npy"
    save_gui_compatible(image, masks, flows, source, output, model_path)
    loaded = np.load(output, allow_pickle=True).item()["masks"]
    if not np.array_equal(masks, loaded):
        raise RuntimeError(f"Saved {suffix} mask array differs from inference result")
    observed_hash = sha256_file(output)
    object_count = int(masks.max())
    expected = FROZEN_SEG[suffix]
    return {
        "suffix": suffix,
        "source": str(source),
        "source_sha256": sha256_file(source),
        "model": str(model_path),
        "model_sha256": sha256_file(model_path),
        "output": str(output),
        "output_bytes": output.stat().st_size,
        "output_sha256": observed_hash,
        "shape": list(masks.shape),
        "dtype": str(masks.dtype),
        "objects": object_count,
        "frozen_objects": expected["objects"],
        "frozen_bytes": expected["bytes"],
        "frozen_sha256": expected["sha256"],
        "object_count_match": object_count == expected["objects"],
        "byte_exact_match": observed_hash == expected["sha256"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--model-root", type=Path, default=Path("/home/server/.cellpose/models")
    )
    args = parser.parse_args()
    if args.output_dir.exists() or args.output_dir.is_symlink():
        raise RuntimeError(f"Fresh output already exists: {args.output_dir}")
    if args.input_dir.resolve() == args.output_dir.resolve():
        raise RuntimeError("Input and output directories overlap")
    require_gpu()
    for name, expected in MODEL_SHA256.items():
        path = args.model_root / name
        if path.is_symlink() or not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Model missing or hash mismatch: {path}")
    args.output_dir.mkdir(parents=False, mode=0o755)
    records: list[dict[str, object]] = []
    for model_name, suffixes in (
        ("dapi_Keyence4", ("dapi",)),
        ("cfos_channel3", ("GFP", "fos")),
    ):
        model_path = args.model_root / model_name
        model = models.CellposeModel(gpu=True, pretrained_model=str(model_path))
        for suffix in suffixes:
            print(f"segmenting={suffix} model={model_name}", flush=True)
            records.append(
                segment_one(args.input_dir, args.output_dir, suffix, model, model_path)
            )
        del model
        gc.collect()
        torch.cuda.empty_cache()
    exact = all(bool(item["byte_exact_match"]) for item in records)
    counts = all(bool(item["object_count_match"]) for item in records)
    receipt = {
        "schema_version": "figure3-water3-cellpose-gpu-reconstruction-v1",
        "cellpose_version": getattr(cellpose, "version", getattr(cellpose, "__version__", "unknown")),
        "torch_version": torch.__version__,
        "cuda_device": torch.cuda.get_device_name(0),
        "parameters": {
            "diameter": None,
            "flow_threshold": FLOW_THRESHOLD,
            "cellprob_threshold": CELLPROB_THRESHOLD,
            "channels": None,
            "channel_axis": -1,
            "normalization": NORMALIZE_PARAMS,
        },
        "records": records,
        "all_frozen_object_counts_match": counts,
        "all_frozen_files_byte_exact": exact,
        "human_review_required": not exact,
        "publication_final": exact,
    }
    receipt_path = args.output_dir / "WATER3_CELLPOSE_GPU_RECEIPT.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True), flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise
