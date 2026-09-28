#!/usr/bin/env python3
"""Recreate the audited Water3 channel TIFFs from the Leica LIF.

This is a source reconstruction of the pre-incident extraction command.  It
reproduces the three native uint16 maximum projections and their lossless
legacy-hue RGB16 Cellpose inputs pixel-for-pixel when run in the recorded
``microscopy_extract`` environment (AICSImageIO 4.14.0, tifffile 2023.2.28).

The command is fresh-output-only.  It never edits the source LIF, never
replaces an output directory, and never performs segmentation or inference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Sequence

import aicsimageio
import numpy as np
import tifffile
from aicsimageio import AICSImage


PAPER_ROOT = next(path for path in Path(__file__).resolve().parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
DEFAULT_LIF = (
    PAPER_ROOT
    / "Fig3/raw_data/new_water_control_sp8_2026_07_20/"
    "cFos-NeuN-Npy-Ctrl-Agua-40X.lif"
)
DEFAULT_LIFEXT = DEFAULT_LIF.with_suffix(".lifext")

SOURCE_LIF_SHA256 = "db4cbdb72c8e8af6a5ab264d1526083d4bb70c4ece54e9c9ee59a8c90d7671ab"
SOURCE_LIFEXT_SHA256 = "1359dd14c59ad5d558312db0a6ff71bbf1ecaf2660e126fabba039903eb22fc7"
SCENE = "TileScan 1 Merged"
SERIES_INDEX = 12
SOURCE_SHAPE = (1, 4, 21, 1419, 1898)
OUTPUT_SHAPE = (1419, 1898)
SAMPLE_ID = "WATER_NPY3_S01"
HISTORICAL_FIELD_ID = "WATER_NPY3-0001"
ANIMAL_ID = "WATER_NPY3"

CHANNELS = (
    ("dapi", 0, "DAPI", "410.96–468.17 nm detector band", (0, 0, 1)),
    ("GFP", 1, "NPY/GFP", "496.01–543.80 nm detector band", (0, 1, 0)),
    ("fos", 3, "c-FOS", "644.29–742.10 nm detector band", (1, 0, 1)),
)

EXPECTED_NATIVE_SHA256 = {
    "dapi": "5fa403c79840778b86386815f992cdcc17953abaf857fe819534a52c7f0fc012",
    "GFP": "f32d2b74a0ce59ce8952d3b676351819fd1ffc7a10711587761c0bc50f37ce5e",
    "fos": "7a61c31dfc7c6ce5d69162767ec48e23373290c9b57a592caf3c37671a32dd2c",
}
EXPECTED_RGB_SHA256 = {
    "dapi": "4fb94b76b67425cff4b95c1195c096b366a316c65ddef46d56196f9f4f941880",
    "GFP": "c8f5f75d558f2011e9b44354d2886a42ab1a8e8938e47fa17f2632c5b8ce7dc8",
    "fos": "d0b399c3eb46541a9eed3cd6aa54b0122523ecb495b6cfdca941b4b7c7cee9b7",
}
EXPECTED_PROVENANCE_SHA256 = "22871c210374bcbf529dc85adcb0d2d490aa88d564438c4e691df6325c4d43b6"
EXPECTED_RGB_BYTES = {"dapi": 6_044_540, "GFP": 6_565_299, "fos": 8_413_489}
EXPECTED_PIXEL_SHA256 = {
    "dapi": "c9bb43c80f07a019d24e557922d6a78fb50d3f627b54fc217aa6494cd72d239c",
    "GFP": "28fd400469901436b0a5be093d5fb4afe49c56f4905951d51e9e8b4ce23cf10e",
    "fos": "33bc1f0fb9f625d04cef15bb711da44555a554a2c5e4df5a55d4558c08f10aef",
}
OME_UUID = {
    "dapi_native": "00000000-0000-5000-8000-000000000001",
    "GFP_native": "00000000-0000-5000-8000-000000000002",
    "fos_native": "00000000-0000-5000-8000-000000000003",
    "dapi_rgb": "00000000-0000-5000-8000-000000000011",
    "GFP_rgb": "00000000-0000-5000-8000-000000000012",
    "fos_rgb": "00000000-0000-5000-8000-000000000013",
}

# Both now resolve inside the tree: the original SP8 acquisition was copied
# into the raw bundle, so this step runs from any unpacked copy.
HISTORICAL_WATER_DIR = PAPER_ROOT / "Fig3/raw/legacy_experiment_2025_07_28/water"
ORIGINAL_LIF_PATH = HISTORICAL_WATER_DIR / "cFos-NeuN-Npy-Ctrl-Agua-40X.lif"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_regular_hash(path: Path, expected: str, role: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size <= 0:
        raise RuntimeError(f"Missing/non-regular {role}: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise RuntimeError(f"{role} SHA256 mismatch: expected {expected}, observed {observed}")


def require_environment() -> None:
    if aicsimageio.__version__ != "4.14.0":
        raise RuntimeError(
            f"AICSImageIO {aicsimageio.__version__} is not the recorded 4.14.0"
        )
    if tifffile.__version__ != "2023.2.28":
        raise RuntimeError(
            f"tifffile {tifffile.__version__} is not the recorded 2023.2.28"
        )


def write_native(
    stage: Path,
    image: AICSImage,
    suffix: str,
    channel_index: int,
    marker: str,
    mapping_basis: str,
) -> tuple[Path, np.ndarray, dict[str, Any]]:
    array = np.asarray(
        image.get_image_dask_data("ZYX", T=0, C=channel_index).max(axis=0).compute()
    )
    if array.shape != OUTPUT_SHAPE or array.dtype != np.uint16:
        raise RuntimeError(f"Unexpected {marker} projection: {array.shape}/{array.dtype}")
    temporary_name = stage / f"{HISTORICAL_FIELD_ID}_{suffix}.tif"
    px = image.physical_pixel_sizes
    tifffile.imwrite(
        temporary_name,
        array,
        ome=True,
        photometric="minisblack",
        compression="zlib",
        metadata={
            "UUID": OME_UUID[f"{suffix}_native"],
            "axes": "YX",
            "Name": f"{HISTORICAL_FIELD_ID} {marker} max-Z",
            "PhysicalSizeX": float(px.X),
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": float(px.Y),
            "PhysicalSizeYUnit": "µm",
        },
    )
    if not np.array_equal(array, tifffile.imread(temporary_name)):
        raise RuntimeError(f"Native TIFF round-trip failed: {temporary_name}")
    observed = sha256_file(temporary_name)
    pixel_sha256 = hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()
    if pixel_sha256 != EXPECTED_PIXEL_SHA256[suffix]:
        raise RuntimeError(f"Native {marker} pixels differ from frozen receipt: {pixel_sha256}")
    record = {
        "file": f"{HISTORICAL_FIELD_ID}_{suffix}.tif",
        "marker": marker,
        "source_channel_index": channel_index,
        "mapping_basis": mapping_basis,
        "projection": "maximum across all 21 Z planes",
        "shape_yx": list(array.shape),
        "dtype": str(array.dtype),
        "min": int(array.min()),
        "max": int(array.max()),
        "mean": float(array.mean()),
        "saturated_65535_pixels": int(np.count_nonzero(array == 65535)),
        "sha256": observed,
        "pixel_sha256": pixel_sha256,
    }
    return temporary_name, array, record


def build_initial_provenance(
    lif: Path,
    lifext: Path,
    image: AICSImage,
    records: list[dict[str, Any]],
) -> dict[str, Any]:
    px = image.physical_pixel_sizes
    return {
        "sample_id": HISTORICAL_FIELD_ID,
        "condition": "Water",
        "biological_replicate": "yes (author-confirmed)",
        "animal_id_status": "provisional; replace WATER_NPY3 with true animal ID before inference",
        "source_lif_original_path": str(ORIGINAL_LIF_PATH),
        "source_lif_staged_path": str(lif),
        "source_lif_sha256": sha256_file(lif),
        "source_lifext_staged_path": str(lifext),
        "source_lifext_sha256": sha256_file(lifext),
        "selected_scene_name": SCENE,
        "bioformats_series_index": SERIES_INDEX,
        "excluded_scene": "TileScan 1 Merged_Lng (processed intensities; not used)",
        "source_shape_TCZYX": list(SOURCE_SHAPE),
        "physical_size_x_um": float(px.X),
        "physical_size_y_um": float(px.Y),
        "z_step_um": float(px.Z),
        "channel_order_validation": (
            "detector emission windows, Leica LUTs, and visual morphology; "
            "OME excitation fields are mis-associated and were not used"
        ),
        "full_physical_channel_order": {
            "C0": "DAPI",
            "C1": "NPY/GFP",
            "C2": "NeuN (failed)",
            "C3": "c-FOS",
        },
        "processing": [
            "ordinary stitched merged scene",
            "channel separation",
            "maximum projection across 21 Z planes",
        ],
        "not_performed": [
            "intensity rescaling",
            "8-bit conversion",
            "LUT baking",
            "denoising",
            "LIGHTNING/deconvolution",
            "segmentation",
            "registration",
        ],
        "failed_channel": "NeuN / C2 retained only in source LIF; no Cellpose TIFF exported",
        "outputs": records,
        "software": {
            "aicsimageio": aicsimageio.__version__,
            "tifffile": tifffile.__version__,
        },
    }


def write_rgb(
    stage: Path,
    suffix: str,
    marker: str,
    color: tuple[int, int, int],
    gray: np.ndarray,
) -> dict[str, Any]:
    rgb = np.zeros(gray.shape + (3,), dtype=np.uint16)
    for channel, enabled in enumerate(color):
        if enabled:
            rgb[..., channel] = gray
    path = stage / f"{SAMPLE_ID}_{suffix}.tif"
    tifffile.imwrite(
        path,
        rgb,
        ome=True,
        photometric="rgb",
        compression="zlib",
        metadata={
            "UUID": OME_UUID[f"{suffix}_rgb"],
            "axes": "YXS",
            "Name": f"{SAMPLE_ID} {suffix} legacy-color RGB16 max-Z",
            "PhysicalSizeX": 0.5687380073800737,
            "PhysicalSizeXUnit": "µm",
            "PhysicalSizeY": 0.5687377997179125,
            "PhysicalSizeYUnit": "µm",
        },
    )
    if not np.array_equal(rgb, tifffile.imread(path)):
        raise RuntimeError(f"RGB TIFF round-trip failed: {path}")
    observed = sha256_file(path)
    pixel_sha256 = hashlib.sha256(np.ascontiguousarray(gray).tobytes()).hexdigest()
    if pixel_sha256 != EXPECTED_PIXEL_SHA256[suffix]:
        raise RuntimeError(f"RGB {marker} source pixels differ from frozen receipt: {pixel_sha256}")
    return {
        "file": path.name,
        "marker": marker,
        "role": "Cellpose input with legacy Water hue encoding",
        "source_native_master": f"native_uint16/{SAMPLE_ID}_{suffix}_native16.tif",
        "legacy_rgb_signal_channels": list(color),
        "color": {"DAPI": "blue", "NPY/GFP": "green", "c-FOS": "magenta"}[marker],
        "shape_yxs": list(rgb.shape),
        "dtype": str(rgb.dtype),
        "intensity_rescaled": False,
        "sha256": observed,
        "source_pixel_sha256": pixel_sha256,
    }


def finalize_current_provenance(
    initial: dict[str, Any],
    native_records: list[dict[str, Any]],
    rgb_records: list[dict[str, Any]],
) -> dict[str, Any]:
    provenance = dict(initial)
    provenance["sample_id"] = SAMPLE_ID
    provenance["animal_id"] = ANIMAL_ID
    provenance["animal_id_status"] = (
        "provisional; replace WATER_NPY3 with true animal ID before inference"
    )
    masters: list[dict[str, Any]] = []
    for record in native_records:
        updated = dict(record)
        suffix = {"DAPI": "dapi", "NPY/GFP": "GFP", "c-FOS": "fos"}[record["marker"]]
        updated["file"] = f"native_uint16/{SAMPLE_ID}_{suffix}_native16.tif"
        updated["role"] = "lossless native grayscale max-projection master"
        masters.append(updated)
    provenance["native_uint16_masters"] = masters
    provenance["outputs"] = rgb_records
    provenance["cellpose_input_format"] = (
        "single-page RGB uint16; signal copied losslessly into legacy Water color samples"
    )
    provenance["legacy_color_encoding"] = {
        "DAPI": "RGB=(0,0,I)",
        "NPY/GFP": "RGB=(0,I,0)",
        "c-FOS": "RGB=(I,0,I)",
    }
    provenance["processing"] = [
        value for value in provenance["processing"] if value != "channel separation"
    ] + [
        "channel separation to native uint16 grayscale masters",
        "lossless legacy-hue RGB16 encoding for Cellpose inputs",
    ]
    provenance["not_performed"] = [
        value
        for value in provenance["not_performed"]
        if value not in {"LUT baking", "8-bit conversion", "contrast/gamma adjustment"}
    ]
    provenance["not_performed"].extend(["8-bit conversion", "contrast/gamma adjustment"])
    return provenance


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_readme(path: Path) -> None:
    path.write_text(
        """WATER_NPY3_S01 — Cellpose input

This is one additional Water biological replicate acquired on Leica SP8 40X.
The main TIFFs use the same hue convention as WATER_NPY1/2:

  WATER_NPY3_S01_dapi.tif  — blue
  WATER_NPY3_S01_fos.tif   — magenta
  WATER_NPY3_S01_GFP.tif   — green (NPY/GFP)

They are lossless RGB uint16 images: native intensity was copied into the legacy
color sample(s) without contrast stretching, gamma adjustment, clipping, or
8-bit conversion. The original grayscale uint16 projections are preserved under
native_uint16/ and should not be separately segmented.

Verified channel order: C0 DAPI, C1 NPY/GFP, C2 NeuN, C3 c-FOS. NeuN failed and
has no quantitative TIFF; its original pixels remain in the LIF.
""",
        encoding="utf-8",
    )


def reconstruct(lif: Path, lifext: Path, output_dir: Path) -> None:
    require_environment()
    lif = lif.expanduser().resolve()
    lifext = lifext.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists() or output_dir.is_symlink():
        raise RuntimeError(f"Output must be absent: {output_dir}")
    if not output_dir.parent.is_dir() or output_dir.parent.is_symlink():
        raise RuntimeError(f"Output parent must be an existing regular directory: {output_dir.parent}")
    require_regular_hash(lif, SOURCE_LIF_SHA256, "Water3 LIF")
    require_regular_hash(lifext, SOURCE_LIFEXT_SHA256, "Water3 LIFEXT")

    stage = Path(tempfile.mkdtemp(prefix=f".{output_dir.name}.stage-", dir=output_dir.parent))
    try:
        image = AICSImage(lif)
        image.set_scene(SCENE)
        if tuple(image.shape) != SOURCE_SHAPE:
            raise RuntimeError(f"Unexpected Water3 source shape: {image.shape}")
        native_records: list[dict[str, Any]] = []
        arrays: dict[str, np.ndarray] = {}
        native_paths: dict[str, Path] = {}
        for suffix, channel, marker, basis, _color in CHANNELS:
            path, array, record = write_native(stage, image, suffix, channel, marker, basis)
            native_records.append(record)
            arrays[suffix] = array
            native_paths[suffix] = path
        initial = build_initial_provenance(lif, lifext, image, native_records)

        native_dir = stage / "native_uint16"
        native_dir.mkdir()
        for suffix, path in native_paths.items():
            os.replace(path, native_dir / f"{SAMPLE_ID}_{suffix}_native16.tif")

        rgb_records = [
            write_rgb(stage, suffix, marker, color, arrays[suffix])
            for suffix, _channel, marker, _basis, color in CHANNELS
        ]
        provenance = finalize_current_provenance(initial, native_records, rgb_records)
        provenance_path = stage / f"{SAMPLE_ID}_extraction_provenance.json"
        write_json(provenance_path, provenance)
        observed_provenance = sha256_file(provenance_path)
        write_readme(stage / f"{SAMPLE_ID}_CELLPOSE_README.txt")
        receipt = {
            "schema_version": "figure3-water3-extraction-reconstruction-v1",
            "source_lif_path": str(lif),
            "source_lif_sha256": SOURCE_LIF_SHA256,
            "source_lifext_path": str(lifext),
            "source_lifext_sha256": SOURCE_LIFEXT_SHA256,
            "current_provenance_sha256": observed_provenance,
            "preincident_container_sha256_receipts": {
                "provenance": EXPECTED_PROVENANCE_SHA256,
                "native_tiff": EXPECTED_NATIVE_SHA256,
                "rgb_tiff": EXPECTED_RGB_SHA256,
            },
            "pixel_sha256": EXPECTED_PIXEL_SHA256,
            "output_is_fresh": True,
            "segmentation_performed": False,
            "generator_path": str(Path(__file__).resolve()),
            "generator_sha256": sha256_file(Path(__file__).resolve()),
        }
        write_json(stage / "WATER3_EXTRACTION_RECONSTRUCTION_RECEIPT.json", receipt)
        os.replace(stage, output_dir)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-lif", type=Path, default=DEFAULT_LIF)
    parser.add_argument("--source-lifext", type=Path, default=DEFAULT_LIFEXT)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reconstruct(args.source_lif, args.source_lifext, args.output_dir)
    print(f"output_dir\t{args.output_dir.expanduser().resolve()}")
    print("frozen_pixel_hashes\tPASS")
    print("segmentation_performed\tFalse")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
