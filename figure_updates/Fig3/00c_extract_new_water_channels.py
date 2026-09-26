#!/usr/bin/env python3
"""Extract the 2026 Zeiss LSM 780 Water/NPY controls from their CZI files.

The 24 September 2026 acquisitions use the same LSM 780 and the same
LD LCI Plan-Apochromat 25x/0.8 objective as the original Sucrose and Allulose
animals, but sample the field at twice their spatial frequency
(0.332106 um/pixel instead of 0.664213 um/pixel).  This command writes the
native maximum-Z projections and an exact 2x2 block-mean reduction whose pixel
size equals the legacy Zeiss pitch bit-for-bit, so that the accepted Cellpose
models, the >=1 pixel NPY rule and the >=5 pixel c-FOS rule apply at the same
physical scale as they do for the sugar cohort.  An acquisition with an odd
pixel extent is trimmed by one row or column first; the trim is recorded.

The command is fresh-output-only unless --force is given.  It never edits a
source CZI and never performs segmentation or inference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

import aicsimageio
import numpy as np
import tifffile
from aicsimageio import AICSImage

PAPER_ROOT = next(
    path
    for path in Path(__file__).resolve().parents
    if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir() and (path / "README.txt").is_file()))
)
RAW_DIR_NAME = "new_water_control_LSM780_24_SEPT_2026"
DEFAULT_RAW = PAPER_ROOT / "Fig3/raw" / RAW_DIR_NAME
FALLBACK_RAW = PAPER_ROOT / "Fig3/raw_data" / RAW_DIR_NAME
DEFAULT_OUTPUT = DEFAULT_RAW / "derived"

# Animal identity, specimen label recorded on the slide, and its acquisition.
SPECS: dict[str, dict[str, str]] = {
    "WATER_NPY4": {"label": "N132-5", "czi": "N132-5 CTRL NPY GFP 555 647 25X - 001.czi"},
    "WATER_NPY5": {"label": "N132-3", "czi": "N132-3 CTRL NPY GFP 555 647 25X - 001.czi"},
}
EXPECTED_CHANNELS = 4
EXPECTED_Z = 9
NATIVE_PIXEL_UM = 0.332106334980716
LEGACY_ZEISS_PIXEL_UM = 0.664212669961432
BIN = 2

# Physical channel order read from the CZI detection windows, not from the
# generic Ch1/ChS1 display names.
CHANNELS = (
    ("dapi", 0, "DAPI", "410.00-487.57 nm detection window", (0, 0, 1), True),
    ("GFP", 1, "NPY/GFP", "490.30-560.32 nm detection window", (0, 1, 0), True),
    ("fos", 2, "c-FOS", "638.01-747.00 nm detection window", (1, 0, 1), True),
    ("NeuN", 3, "NeuN", "569.07-621.58 nm detection window", (1, 0, 0), False),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pixel_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def require_environment() -> None:
    if aicsimageio.__version__ != "4.14.0":
        raise RuntimeError(f"AICSImageIO {aicsimageio.__version__} is not the recorded 4.14.0")


def block_mean_uint16(array: np.ndarray, factor: int) -> np.ndarray:
    """Exact 2x2 block mean with round-half-away-from-zero, kept in uint16.

    Averaging four uint16 samples in float64 is exact, so the only
    approximation is rounding the mean back onto the integer grid.
    """
    height, width = array.shape
    if height % factor or width % factor:
        raise RuntimeError(f"Shape {array.shape} is not divisible by {factor}")
    blocks = array.reshape(height // factor, factor, width // factor, factor).astype(np.float64)
    return np.floor(blocks.mean(axis=(1, 3)) + 0.5).astype(np.uint16)


def describe(array: np.ndarray) -> dict[str, Any]:
    return {
        "shape_yx": list(array.shape),
        "dtype": str(array.dtype),
        "min": int(array.min()),
        "max": int(array.max()),
        "mean": float(array.mean()),
        "saturated_65535_pixels": int(np.count_nonzero(array == 65535)),
        "pixel_sha256": pixel_sha256(array),
    }


def write_gray(path: Path, array: np.ndarray, name: str, pixel_um: float) -> dict[str, Any]:
    tifffile.imwrite(
        path, array, ome=True, photometric="minisblack", compression="zlib",
        metadata={"axes": "YX", "Name": name, "PhysicalSizeX": pixel_um, "PhysicalSizeXUnit": "µm",
                  "PhysicalSizeY": pixel_um, "PhysicalSizeYUnit": "µm"},
    )
    if not np.array_equal(array, tifffile.imread(path)):
        raise RuntimeError(f"Grayscale TIFF round-trip failed: {path}")
    return {"file": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size, **describe(array)}


def write_rgb(path: Path, gray: np.ndarray, color: tuple[int, int, int], name: str, pixel_um: float) -> dict[str, Any]:
    rgb = np.zeros(gray.shape + (3,), dtype=np.uint16)
    for index, enabled in enumerate(color):
        if enabled:
            rgb[..., index] = gray
    tifffile.imwrite(
        path, rgb, ome=True, photometric="rgb", compression="zlib",
        metadata={"axes": "YXS", "Name": name, "PhysicalSizeX": pixel_um, "PhysicalSizeXUnit": "µm",
                  "PhysicalSizeY": pixel_um, "PhysicalSizeYUnit": "µm"},
    )
    if not np.array_equal(rgb, tifffile.imread(path)):
        raise RuntimeError(f"RGB TIFF round-trip failed: {path}")
    if not np.array_equal(rgb.max(axis=-1), gray):
        raise RuntimeError(f"RGB encoding is not lossless: {path}")
    return {"file": path.name, "sha256": sha256_file(path), "bytes": path.stat().st_size, **describe(gray)}


def extract_one(animal: str, czi: Path, output: Path, force: bool) -> dict[str, Any]:
    spec = SPECS[animal]
    sample_id = animal + "_S01"
    existing = output / f"{sample_id}_extraction_provenance.json"
    if existing.exists() and not force:
        raise RuntimeError(f"Fresh output already exists: {existing}")
    image = AICSImage(str(czi))
    if len(image.shape) != 5 or image.shape[1] != EXPECTED_CHANNELS or image.shape[2] != EXPECTED_Z:
        raise RuntimeError(f"Unexpected source shape {image.shape} for {animal}")
    px = image.physical_pixel_sizes
    for value in (px.X, px.Y):
        if not np.isclose(value, NATIVE_PIXEL_UM, rtol=0, atol=1e-12):
            raise RuntimeError(f"Unexpected native pixel size {px} for {animal}")
    if not np.isclose(px.X * BIN, LEGACY_ZEISS_PIXEL_UM, rtol=0, atol=1e-12):
        raise RuntimeError("The binned pitch does not reproduce the legacy Zeiss pixel size")
    native_shape = tuple(int(value) for value in image.shape[3:])
    trimmed_shape = tuple(value - value % BIN for value in native_shape)
    trim = [native - kept for native, kept in zip(native_shape, trimmed_shape)]

    records: list[dict[str, Any]] = []
    output.mkdir(parents=True, exist_ok=True)
    (output / "native_uint16").mkdir(exist_ok=True)
    for suffix, index, marker, basis, color, cellpose_input in CHANNELS:
        native = np.asarray(image.get_image_dask_data("ZYX", T=0, C=index).max(axis=0).compute())
        if native.shape != native_shape or native.dtype != np.uint16:
            raise RuntimeError(f"Unexpected {marker} projection: {native.shape}/{native.dtype}")
        binned = block_mean_uint16(native[: trimmed_shape[0], : trimmed_shape[1]], BIN)
        native_record = write_gray(output / "native_uint16" / f"{sample_id}_{suffix}_native16.tif", native,
                                   f"{sample_id} {marker} native max-Z", NATIVE_PIXEL_UM)
        name = f"{sample_id}_{suffix}.tif" if cellpose_input else f"{sample_id}_{suffix}_display_only.tif"
        analysis_record = write_rgb(output / name, binned, color,
                                    f"{sample_id} {suffix} legacy-color RGB16 2x2 block-mean max-Z",
                                    LEGACY_ZEISS_PIXEL_UM)
        records.append({
            "marker": marker, "source_channel_index": index, "mapping_basis": basis,
            "projection": f"maximum across all {EXPECTED_Z} Z planes",
            "legacy_color_rgb": list(color), "used_as_cellpose_input": cellpose_input,
            "native_uint16_master": {"pixel_um": NATIVE_PIXEL_UM, **native_record},
            "analysis_rgb16": {"pixel_um": LEGACY_ZEISS_PIXEL_UM,
                               "reduction": "2x2 block mean, round half away from zero", **analysis_record},
        })
        print(f"{animal}: extracted={marker} native={native.shape} analysis={binned.shape}", flush=True)

    provenance = {
        "schema": "figure3-new-water-lsm780-extraction-v1",
        "animal_id": animal, "native_animal_label": spec["label"], "sample_id": sample_id,
        "condition": "Water", "acquisition_date": "2026-09-24",
        "source_czi_path": str(czi), "source_czi_sha256": sha256_file(czi), "source_czi_bytes": czi.stat().st_size,
        "scene": image.scenes[0], "source_shape_TCZYX": [int(value) for value in image.shape],
        "instrument": "Zeiss LSM 780, AxioObserver",
        "objective": "LD LCI Plan-Apochromat 25x/0.8 Imm Korr DIC M27",
        "same_instrument_and_objective_as_sugar_cohort": True,
        "native_physical_size_x_um": float(px.X), "native_physical_size_y_um": float(px.Y),
        "z_step_um": float(px.Z), "analysis_physical_size_um": LEGACY_ZEISS_PIXEL_UM,
        "native_shape_yx": list(native_shape), "analysis_input_shape_yx": list(trimmed_shape),
        "trimmed_rows_columns_before_binning": trim,
        "trim_rationale": ("A 2x2 block mean needs an even extent; at most one row and one column, "
                           "0.332106 um each, are dropped from the high edge before reduction."),
        "analysis_scale_rationale": (
            "The native acquisition samples twice as finely as the legacy Zeiss sugar images. "
            "An exact 2x2 block mean restores the identical 0.664212669961432 um pitch, so the "
            "accepted Cellpose models and the pixel-count association thresholds act at the same "
            "physical scale in every animal of the comparison."),
        "channel_order_validation": "CZI detection windows and laser lines; generic Ch/ChS display names were not used",
        "full_physical_channel_order": {"C0": "DAPI", "C1": "NPY/GFP", "C2": "c-FOS", "C3": "NeuN"},
        "legacy_color_encoding": {"DAPI": "RGB=(0,0,I)", "NPY/GFP": "RGB=(0,I,0)", "c-FOS": "RGB=(I,0,I)", "NeuN": "RGB=(I,0,0)"},
        "processing": ["channel separation", f"maximum projection across {EXPECTED_Z} Z planes",
                       "even-extent trim when required", "exact 2x2 block-mean reduction to the legacy Zeiss pitch"],
        "not_performed": ["intensity rescaling", "8-bit conversion", "LUT baking", "denoising",
                          "deconvolution", "segmentation", "registration"],
        "outputs": records,
        "software": {"aicsimageio": aicsimageio.__version__, "tifffile": tifffile.__version__, "numpy": np.__version__},
    }
    (output / f"{sample_id}_extraction_provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--animal", action="append", choices=sorted(SPECS), help="Repeatable; default is every 2026 Water animal")
    parser.add_argument("--raw-dir", type=Path, default=None, help="Directory holding the CZI files")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--force", action="store_true", help="Rewrite existing outputs for the selected animals")
    args = parser.parse_args()

    require_environment()
    raw = args.raw_dir or (DEFAULT_RAW if DEFAULT_RAW.is_dir() else FALLBACK_RAW)
    animals = args.animal or sorted(SPECS)
    summary = []
    for animal in animals:
        czi = raw / SPECS[animal]["czi"]
        if not czi.is_file():
            raise RuntimeError(f"Source acquisition is missing: {czi}")
        provenance = extract_one(animal, czi, args.output_dir, args.force)
        summary.append({"animal_id": animal, "native_animal_label": provenance["native_animal_label"],
                        "analysis_input_shape_yx": provenance["analysis_input_shape_yx"],
                        "trimmed_rows_columns_before_binning": provenance["trimmed_rows_columns_before_binning"]})
    print(json.dumps({"status": "PASS", "output_dir": str(args.output_dir), "animals": summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
