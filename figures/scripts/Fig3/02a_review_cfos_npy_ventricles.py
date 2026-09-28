#!/usr/bin/env python3
"""Serve a source-bound HIL editor for Figure 3 cartoon ventricular lumens.

The reviewer draws one closed ventricular-lumen polygon in the exact registered
866 x 1374 common frame for Water and Allulose. Acceptance writes a binary PNG
mask plus a JSON receipt bound to the raw DAPI TIFF, registered DAPI labels,
raw-DAPI lumen proposal, reviewer identity, timestamp and polygon vertices.
Existing accepted pairs are copied into a timestamped history directory before
replacement; nothing is deleted.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import importlib.util
import json
import mimetypes
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
from typing import Any, Mapping, Sequence
from urllib.parse import parse_qs, urlparse
from uuid import uuid4

import cv2
import numpy as np
from PIL import Image
from pypdf import PdfReader
from scipy.ndimage import binary_fill_holes, label


DEFAULT_PORT = 34247
HIL_RELATIVE = Path("Fig3/hil_review/ventricle_cartoon_20260827_v1")
SCHEMA = "fig3_cartoon_ventricle_hil_v1"
CONDITIONS = ("Water", "Allulose")
FIGURE_STEM = "Figure_3_cFos_NPY"
FIGURE_DPI = 600
FIGURE_LOG_TAIL_CHARS = 24_000
HIL_FILENAMES = (
    "water_ventricle_mask.png",
    "water_ventricle_receipt.json",
    "allulose_ventricle_mask.png",
    "allulose_ventricle_receipt.json",
)
MACHINE_FIGURE_PYTHON = Path(
    "/home/server/anaconda3/envs/paper_apotome_repro/bin/python"
)


def select_figure_python() -> Path:
    configured = os.environ.get("FIG3_PYTHON") or os.environ.get("PAPER_PYTHON")
    if configured:
        return Path(configured).expanduser().resolve()
    if MACHINE_FIGURE_PYTHON.is_file():
        return MACHINE_FIGURE_PYTHON
    return Path(sys.executable).resolve()


FIGURE_PYTHON = select_figure_python()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover_paper_root(explicit: Path | None) -> Path:
    if explicit is not None:
        root = explicit.expanduser().resolve()
    else:
        root = next(
            parent
            for parent in Path(__file__).resolve().parents
            if parent.name == "Paper"
            or (
                (parent / "README.txt").is_file()
                and (parent / "scripts" / "setup").is_dir()
            )
        )
    if not root.is_dir() or root.is_symlink():
        raise RuntimeError(f"Paper root is missing/non-normal: {root}")
    return root


def load_renderer(root: Path):
    path = root / "Fig3/02_make_cfos_npy_cartoons.py"
    if not path.is_file() or path.is_symlink():
        raise RuntimeError(
            f"Canonical Figure 3 cartoon renderer is missing/non-normal: {path}"
        )
    name = "fig3_cartoon_renderer_for_ventricle_hil"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import Figure 3 cartoon renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def encode_png(array: np.ndarray) -> bytes:
    ok, encoded = cv2.imencode(".png", np.ascontiguousarray(array))
    if not ok:
        raise RuntimeError("OpenCV failed to encode PNG")
    return encoded.tobytes()


def immutable_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        if path.is_file() and not path.is_symlink() and path.read_bytes() == payload:
            return
        raise RuntimeError(f"Refusing to replace changed review input: {path}")
    with path.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def validate_or_create_review_manifest(path: Path, manifest: Mapping[str, Any]) -> None:
    """Preserve the true input-generation renderer while validating content identity."""
    payload = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if not path.exists() and not path.is_symlink():
        immutable_write(path, payload)
        return
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"Review input manifest is missing/non-normal: {path}")
    try:
        existing = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"Cannot parse immutable review input manifest: {path}"
        ) from exc
    volatile = {"renderer_path", "renderer_sha256"}
    stable_existing = {
        key: value for key, value in existing.items() if key not in volatile
    }
    stable_current = {
        key: value for key, value in manifest.items() if key not in volatile
    }
    if stable_existing != stable_current:
        raise RuntimeError(
            "Review input content/source identity changed; use a new versioned HIL directory"
        )
    recorded_hash = str(existing.get("renderer_sha256", ""))
    if len(recorded_hash) != 64 or any(
        character not in "0123456789abcdef" for character in recorded_hash
    ):
        raise RuntimeError(
            "Immutable review input manifest has an invalid generator hash"
        )


def atomic_replace(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def make_background(dapi: np.ndarray, tissue: np.ndarray) -> bytes:
    height, width = dapi.shape
    image = np.full((height, width, 3), (248, 246, 245), dtype=np.uint8)
    image[tissue] = (242, 237, 232)
    nuclei = dapi > 0
    nucleus_bgr = np.asarray((239, 185, 158), dtype=np.float32)
    image[nuclei] = np.rint(
        image[nuclei].astype(np.float32) * 0.58 + nucleus_bgr * 0.42
    ).astype(np.uint8)
    outline = np.asarray(binary_fill_holes(tissue), dtype=np.uint8)
    contours, _ = cv2.findContours(outline, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    cv2.drawContours(image, contours, -1, (82, 64, 53), 3, lineType=cv2.LINE_AA)
    return encode_png(image)


def make_source_overlay(source: np.ndarray) -> bytes:
    overlay = np.zeros((*source.shape, 4), dtype=np.uint8)
    overlay[source] = (210, 40, 220, 48)
    contours, _ = cv2.findContours(
        np.asarray(source, dtype=np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE
    )
    cv2.drawContours(
        overlay, contours, -1, (210, 40, 220, 210), 3, lineType=cv2.LINE_AA
    )
    return encode_png(overlay)


class ReviewState:
    def __init__(
        self,
        root: Path,
        output_dir: Path,
        figure_build_root: Path | None = None,
    ) -> None:
        self.root = root
        self.output_dir = output_dir.expanduser().resolve()
        if self.output_dir.is_symlink():
            raise RuntimeError(
                f"HIL output directory cannot be a symlink: {self.output_dir}"
            )
        self.output_dir.mkdir(parents=True, exist_ok=True)
        requested_build_root = (
            figure_build_root.expanduser()
            if figure_build_root is not None
            else self.output_dir / "figure_builds"
        )
        if requested_build_root.is_symlink():
            raise RuntimeError(
                f"Figure build root cannot be a symlink: {requested_build_root}"
            )
        requested_build_root.mkdir(parents=True, exist_ok=True)
        self.figure_build_root = requested_build_root.resolve()
        if not self.figure_build_root.is_dir() or self.figure_build_root.is_symlink():
            raise RuntimeError(
                f"Figure build root is missing/non-normal: {self.figure_build_root}"
            )
        if not FIGURE_PYTHON.is_file() or not os.access(FIGURE_PYTHON, os.X_OK):
            raise RuntimeError(f"Figure Python is not executable: {FIGURE_PYTHON}")
        self._hil_lock = threading.RLock()
        self._build_lock = threading.Lock()
        self._build: dict[str, Any] = {
            "schema": "fig3_hil_figure_build_status_v1",
            "status": "idle",
            "running": False,
            "ok": None,
            "step": "Accept both HIL polygons to enable rendering.",
            "completed_steps": 0,
            "total_steps": 9,
            "build_dir": None,
            "started_at_utc": None,
            "finished_at_utc": None,
            "error": None,
            "log": "",
            "outputs": [],
        }
        self.renderer = load_renderer(root)
        self.samples: dict[str, dict[str, Any]] = {}
        specs = (
            self.renderer.WATER_CARTOON_SPEC,
            self.renderer.ALLULOSE_CARTOON_SPEC,
        )
        review_inputs = self.output_dir / "review_inputs"
        for spec in specs:
            data = self.renderer._matched_load_inputs(spec, root)
            dapi = self.renderer._matched_display_transform(data.dapi_labels, spec)
            tissue = self.renderer._tissue_envelope(dapi)
            source, source_receipt = self.renderer._matched_lumen_source_result(data)
            slug = spec.condition.lower()
            background_payload = make_background(dapi, tissue)
            overlay_payload = make_source_overlay(source)
            background_path = review_inputs / f"{slug}_registered_dapi_background.png"
            overlay_path = review_inputs / f"{slug}_raw_dapi_lumen_proposal.png"
            immutable_write(background_path, background_payload)
            immutable_write(overlay_path, overlay_payload)
            self.samples[spec.condition] = {
                "data": data,
                "dapi": dapi,
                "source": source,
                "source_receipt": source_receipt,
                "background_path": background_path,
                "overlay_path": overlay_path,
                "background_sha256": sha256_bytes(background_payload),
                "overlay_sha256": sha256_bytes(overlay_payload),
            }
        manifest = {
            "schema": "fig3_cartoon_ventricle_hil_review_inputs_v1",
            "renderer_path": str(
                (root / "Fig3/02_make_cfos_npy_cartoons.py").resolve()
            ),
            "renderer_sha256": sha256_file(root / "Fig3/02_make_cfos_npy_cartoons.py"),
            "coordinate_frame_yx": [866, 1374],
            "conditions": {
                condition: {
                    "sample_id": item["data"].spec.sample_id,
                    "animal_id": item["data"].spec.animal_id,
                    "background": item["background_path"].name,
                    "background_sha256": item["background_sha256"],
                    "source_overlay": item["overlay_path"].name,
                    "source_overlay_sha256": item["overlay_sha256"],
                    "raw_dapi_lumen_source_binary_sha256": self.renderer.binary_sha256(
                        item["source"]
                    ),
                }
                for condition, item in self.samples.items()
            },
        }
        validate_or_create_review_manifest(
            review_inputs / "REVIEW_INPUT_MANIFEST.json", manifest
        )

    def paths(self, condition: str) -> tuple[Path, Path]:
        slug = condition.lower()
        return (
            self.output_dir / f"{slug}_ventricle_mask.png",
            self.output_dir / f"{slug}_ventricle_receipt.json",
        )

    def annotation(self, condition: str) -> dict[str, Any]:
        with self._hil_lock:
            return self._annotation_unlocked(condition)

    def _annotation_unlocked(self, condition: str) -> dict[str, Any]:
        mask_path, receipt_path = self.paths(condition)
        if not mask_path.exists() and not receipt_path.exists():
            return {"found": False, "condition": condition}
        if not mask_path.is_file() or not receipt_path.is_file():
            raise RuntimeError(f"Incomplete accepted HIL pair for {condition}")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        data = self.samples[condition]["data"]
        source = self.samples[condition]["source"]
        self.renderer._load_hil_lumen_mask(data, source, self.output_dir)
        return {"found": True, **receipt}

    def config(self) -> dict[str, Any]:
        rows = []
        for condition in CONDITIONS:
            item = self.samples[condition]
            annotation = self.annotation(condition)
            rows.append(
                {
                    "condition": condition,
                    "sample_id": item["data"].spec.sample_id,
                    "animal_id": item["data"].spec.animal_id,
                    "width": int(item["dapi"].shape[1]),
                    "height": int(item["dapi"].shape[0]),
                    "accepted": bool(annotation.get("found")),
                }
            )
        accepted = sum(int(row["accepted"]) for row in rows)
        build = self.figure_status()
        return {
            "schema": "fig3_cartoon_ventricle_hil_config_v1",
            "sections": rows,
            "accepted": accepted,
            "total": len(rows),
            "output_dir": str(self.output_dir),
            "figure_build_root": str(self.figure_build_root),
            "figure_button_enabled": accepted == len(rows) and not build["running"],
            "figure_dpi": FIGURE_DPI,
        }

    def _archive_existing(self, condition: str) -> Path | None:
        mask_path, receipt_path = self.paths(condition)
        existing = [
            path
            for path in (mask_path, receipt_path)
            if path.exists() or path.is_symlink()
        ]
        if not existing:
            return None
        if any(not path.is_file() or path.is_symlink() for path in existing):
            raise RuntimeError(
                f"Cannot archive non-normal existing HIL pair for {condition}"
            )
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        archive = (
            self.output_dir
            / "history"
            / condition.lower()
            / f"{stamp}_{uuid4().hex[:8]}"
        )
        archive.mkdir(parents=True, exist_ok=False)
        records = []
        for path in existing:
            target = archive / path.name
            shutil.copy2(path, target)
            if sha256_file(target) != sha256_file(path):
                raise RuntimeError(f"Archived HIL copy mismatch: {path}")
            records.append({"filename": path.name, "sha256": sha256_file(path)})
        immutable_write(
            archive / "ARCHIVE_RECEIPT.json",
            (
                json.dumps(
                    {
                        "schema": "fig3_cartoon_ventricle_hil_archive_v1",
                        "archived_at_utc": datetime.now(timezone.utc).isoformat(),
                        "condition": condition,
                        "files": records,
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            ).encode("utf-8"),
        )
        return archive

    def accept(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        with self._hil_lock:
            return self._accept_unlocked(payload)

    def _accept_unlocked(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        condition = str(payload.get("condition", ""))
        if condition not in self.samples:
            raise ValueError("condition must be Water or Allulose")
        reviewer = str(payload.get("reviewer", "")).strip()
        session = str(payload.get("session", "")).strip()
        notes = str(payload.get("notes", "")).strip()
        if not reviewer or not session:
            raise ValueError("reviewer name and session identifier are required")
        if max(len(reviewer), len(session), len(notes)) > 1000:
            raise ValueError("review fields are too long")
        vertices = np.asarray(payload.get("polygon_vertices_xy", []), dtype=np.float64)
        if (
            vertices.ndim != 2
            or vertices.shape[1:] != (2,)
            or not 3 <= vertices.shape[0] <= 2000
        ):
            raise ValueError("one closed polygon with 3-2000 vertices is required")
        if not np.all(np.isfinite(vertices)):
            raise ValueError("polygon vertices must be finite")
        data = self.samples[condition]["data"]
        height, width = data.spec.output_shape_yx
        vertices = np.round(vertices, 3)
        raster_vertices = np.rint(vertices).astype(np.int32)
        if (
            np.any(raster_vertices[:, 0] < 0)
            or np.any(raster_vertices[:, 0] >= width)
            or np.any(raster_vertices[:, 1] < 0)
            or np.any(raster_vertices[:, 1] >= height)
        ):
            raise ValueError("polygon vertices must stay inside the common frame")
        mask_u8 = np.zeros((height, width), dtype=np.uint8)
        cv2.fillPoly(mask_u8, [raster_vertices.reshape(-1, 1, 2)], 255)
        mask = mask_u8 == 255
        if int(mask.sum()) < 50:
            raise ValueError("accepted ventricle polygon is too small")
        if not np.array_equal(mask, binary_fill_holes(mask)):
            raise ValueError("accepted ventricle polygon contains an internal hole")
        component_count = int(label(mask, structure=np.ones((3, 3), dtype=np.uint8))[1])
        path_counts = self.renderer._contour_path_counts(mask)
        if component_count != 1 or path_counts != {"external": 1, "tree": 1}:
            raise ValueError(
                "accepted ventricle must rasterize as exactly one filled component"
            )
        mask_payload = encode_png(mask_u8)
        mask_sha = sha256_bytes(mask_payload)
        ys, xs = np.where(mask)
        mask_path, receipt_path = self.paths(condition)
        item = self.samples[condition]
        accepted_at = datetime.now(timezone.utc).isoformat()
        receipt = {
            "schema": SCHEMA,
            "status": "accepted",
            "condition": condition,
            "sample_id": data.spec.sample_id,
            "animal_id": data.spec.animal_id,
            "coordinate_frame_yx": [height, width],
            "reviewer": reviewer,
            "session": session,
            "accepted_at_utc": accepted_at,
            "notes": notes,
            "polygon_vertices_xy": vertices.tolist(),
            "raster_vertices_sha256_i32_le": hashlib.sha256(
                np.ascontiguousarray(raster_vertices, dtype="<i4").tobytes()
            ).hexdigest(),
            "mask_filename": mask_path.name,
            "mask_file_sha256": mask_sha,
            "mask_geometry": {
                "area_px": int(mask.sum()),
                "bbox_yxyx_half_open": [
                    int(ys.min()),
                    int(xs.min()),
                    int(ys.max()) + 1,
                    int(xs.max()) + 1,
                ],
                "centroid_yx": [float(ys.mean()), float(xs.mean())],
                "binary_sha256": self.renderer.binary_sha256(mask),
                "connected_components_8": component_count,
                "contour_path_counts": path_counts,
            },
            "source_binding": {
                "dapi_tiff_sha256": data.sources["dapi_tiff"].sha256,
                "registered_dapi_binary_sha256": self.renderer.binary_sha256(
                    item["dapi"] > 0
                ),
                "raw_dapi_lumen_source_binary_sha256": self.renderer.binary_sha256(
                    item["source"]
                ),
                "review_background_sha256": item["background_sha256"],
                "raw_dapi_lumen_overlay_sha256": item["overlay_sha256"],
            },
            "semantics": {
                "role": "display-only ventricular lumen for Figure 3 cartoons B/C",
                "analysis_mask": False,
                "changes_cell_calls_or_counts": False,
                "render": "opaque white interior with one black outer boundary",
            },
            "generator": {
                "path": str(Path(__file__).resolve()),
                "sha256": sha256_file(Path(__file__).resolve()),
            },
        }
        receipt_payload = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode(
            "utf-8"
        )
        archive = self._archive_existing(condition)
        atomic_replace(mask_path, mask_payload)
        atomic_replace(receipt_path, receipt_payload)
        self.renderer._load_hil_lumen_mask(data, item["source"], self.output_dir)
        return {
            "status": "accepted",
            "condition": condition,
            "mask_sha256": mask_sha,
            "receipt_sha256": sha256_bytes(receipt_payload),
            "archived_previous_to": str(archive) if archive else None,
        }

    def figure_status(self) -> dict[str, Any]:
        with self._build_lock:
            return json.loads(json.dumps(self._build))

    def _set_build(self, **updates: Any) -> None:
        with self._build_lock:
            self._build.update(updates)

    def _append_build_log(self, text: str) -> None:
        if not text:
            return
        with self._build_lock:
            build_dir_value = self._build.get("build_dir")
            combined = (str(self._build.get("log", "")) + text)[-FIGURE_LOG_TAIL_CHARS:]
            self._build["log"] = combined
        if build_dir_value:
            log_path = Path(str(build_dir_value)) / "BUILD_LOG.txt"
            with log_path.open("a", encoding="utf-8") as handle:
                handle.write(text)
                handle.flush()
                os.fsync(handle.fileno())

    def _new_build_dir(self) -> Path:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        for _ in range(20):
            target = self.figure_build_root / f"figure3_hil_{stamp}_{uuid4().hex[:8]}"
            try:
                target.mkdir(mode=0o700)
                return target
            except FileExistsError:
                continue
        raise RuntimeError(
            "Could not reserve a fresh versioned Figure 3 build directory"
        )

    def _snapshot_hil(
        self,
        build_dir: Path,
        reviewer: str,
        session: str,
    ) -> tuple[Path, dict[str, Any]]:
        snapshot = build_dir / "accepted_ventricle_hil"
        snapshot.mkdir(mode=0o700)
        records = []
        accepted_receipts: dict[str, Any] = {}
        with self._hil_lock:
            for condition in CONDITIONS:
                annotation = self._annotation_unlocked(condition)
                if not annotation.get("found"):
                    raise RuntimeError(
                        "Make Figure requires accepted Water and Allulose HIL polygons"
                    )
                accepted_receipts[condition] = {
                    "reviewer": annotation.get("reviewer"),
                    "session": annotation.get("session"),
                    "accepted_at_utc": annotation.get("accepted_at_utc"),
                    "mask_file_sha256": annotation.get("mask_file_sha256"),
                }
            for filename in HIL_FILENAMES:
                source = self.output_dir / filename
                target = snapshot / filename
                if source.is_symlink() or not source.is_file():
                    raise RuntimeError(
                        f"Accepted HIL source is missing/non-normal: {source}"
                    )
                with source.open("rb") as src, target.open("xb") as dst:
                    shutil.copyfileobj(src, dst, 4 * 1024 * 1024)
                    dst.flush()
                    os.fsync(dst.fileno())
                source_hash = sha256_file(source)
                target_hash = sha256_file(target)
                if source_hash != target_hash:
                    raise RuntimeError(f"Accepted HIL snapshot mismatch: {filename}")
                records.append(
                    {
                        "filename": filename,
                        "bytes": target.stat().st_size,
                        "sha256": target_hash,
                    }
                )
            for condition in CONDITIONS:
                item = self.samples[condition]
                self.renderer._load_hil_lumen_mask(
                    item["data"], item["source"], snapshot
                )
        receipt = {
            "schema": "fig3_hil_build_input_snapshot_v1",
            "status": "PASS",
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "button_reviewer": reviewer,
            "button_session": session,
            "source_hil_dir": str(self.output_dir),
            "snapshot_dir": str(snapshot),
            "accepted_receipts": accepted_receipts,
            "files": records,
        }
        immutable_write(
            build_dir / "HIL_INPUT_SNAPSHOT_RECEIPT.json",
            (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode("utf-8"),
        )
        return snapshot, receipt

    def _figure_commands(
        self, build_dir: Path, snapshot: Path
    ) -> list[tuple[str, list[str]]]:
        bundle = self.root / "Fig3/reproduce_final_20260823"
        scripts = bundle / "scripts"
        raw_root = self.root / "Fig3/raw/legacy_experiment_2025_07_28"
        for required in (
            scripts / "01_analyze_cfos_npy.py",
            scripts / "02_make_cfos_npy_cartoons.py",
            scripts / "05_analyze_spatial_distributions.py",
            scripts / "03_make_figure_3_cfos_npy.py",
            scripts / "02_make_cfos_npy_cartoons_spanish.py",
            scripts / "05_analyze_spatial_distributions_spanish.py",
            scripts / "03_make_figure_3_cfos_npy_spanish.py",
            bundle / "assemble_final.py",
        ):
            if required.is_symlink() or not required.is_file():
                raise RuntimeError(
                    f"Frozen Figure 3 build input is missing/non-normal: {required}"
                )
        if raw_root.is_symlink() or not raw_root.is_dir():
            raise RuntimeError(f"Figure 3 raw root is missing/non-normal: {raw_root}")

        analysis = build_dir / "analysis"
        spatial_en = build_dir / "spatial_english"
        publication_en = build_dir / "english_publication"
        panels_es = build_dir / "panels_spanish"
        spatial_es = build_dir / "spatial_spanish"
        publication_es = build_dir / "spanish_publication"
        final = build_dir / "final"
        work = build_dir / "normalization_work"
        py = str(FIGURE_PYTHON)
        utility_py = str(Path(sys.executable).resolve())
        return [
            (
                "Analyze animal-level Figure 3 data",
                [
                    py,
                    str(scripts / "01_analyze_cfos_npy.py"),
                    "--raw-root",
                    str(raw_root),
                    "--output-dir",
                    str(analysis),
                ],
            ),
            (
                "Render English HIL cartoons at 600 dpi",
                [
                    py,
                    str(scripts / "02_make_cfos_npy_cartoons.py"),
                    "--root",
                    str(self.root),
                    "--output-dir",
                    str(analysis / "panels"),
                    "--dpi",
                    str(FIGURE_DPI),
                    "--ventricle-hil-dir",
                    str(snapshot),
                ],
            ),
            (
                "Render English spatial panels at 600 dpi",
                [
                    py,
                    str(scripts / "05_analyze_spatial_distributions.py"),
                    "--raw-root",
                    str(raw_root),
                    "--output-dir",
                    str(spatial_en),
                    "--dpi",
                    str(FIGURE_DPI),
                ],
            ),
            (
                "Compose English Figure 3 and A-G panels at 600 dpi",
                [
                    py,
                    str(scripts / "03_make_figure_3_cfos_npy.py"),
                    "--analysis-dir",
                    str(analysis),
                    "--output-dir",
                    str(publication_en),
                    "--raw-root",
                    str(raw_root),
                    "--spatial-dir",
                    str(spatial_en),
                    "--per-animal",
                    str(analysis / "per_animal_cfos_npy.csv"),
                    "--statistics",
                    str(analysis / "anova_mwu_statistics.csv"),
                    "--panel-b",
                    str(
                        analysis
                        / "panels/Figure3_Panel_B_Water_marker_positive_nuclei.png"
                    ),
                    "--panel-c",
                    str(
                        analysis
                        / "panels/Figure3_Panel_C_Allulose_marker_positive_nuclei.png"
                    ),
                    "--panel-f",
                    str(spatial_en / "Figure3_Spatial_cFOS_occurrence.png"),
                    "--panel-g",
                    str(spatial_en / "Figure3_Spatial_cFOS_NPY_occurrence.png"),
                    "--figure-name",
                    FIGURE_STEM,
                    "--dpi",
                    str(FIGURE_DPI),
                    "--panel-a-dpi",
                    str(FIGURE_DPI),
                    "--seed",
                    "31",
                ],
            ),
            (
                "Render Spanish HIL cartoons at 600 dpi",
                [
                    py,
                    str(scripts / "02_make_cfos_npy_cartoons_spanish.py"),
                    "--root",
                    str(self.root),
                    "--output-dir",
                    str(panels_es),
                    "--dpi",
                    str(FIGURE_DPI),
                    "--ventricle-hil-dir",
                    str(snapshot),
                ],
            ),
            (
                "Render Spanish spatial panels at 600 dpi",
                [
                    py,
                    str(scripts / "05_analyze_spatial_distributions_spanish.py"),
                    "--raw-root",
                    str(raw_root),
                    "--output-dir",
                    str(spatial_es),
                    "--dpi",
                    str(FIGURE_DPI),
                ],
            ),
            (
                "Compose Spanish Figure 3 and A-G panels at 600 dpi",
                [
                    py,
                    str(scripts / "03_make_figure_3_cfos_npy_spanish.py"),
                    "--analysis-dir",
                    str(analysis),
                    "--output-dir",
                    str(publication_es),
                    "--raw-root",
                    str(raw_root),
                    "--spatial-dir",
                    str(spatial_es),
                    "--per-animal",
                    str(analysis / "per_animal_cfos_npy.csv"),
                    "--statistics",
                    str(analysis / "anova_mwu_statistics.csv"),
                    "--panel-b",
                    str(panels_es / "Figure3_Panel_B_Water_marker_positive_nuclei.png"),
                    "--panel-c",
                    str(
                        panels_es
                        / "Figure3_Panel_C_Allulose_marker_positive_nuclei.png"
                    ),
                    "--panel-f",
                    str(spatial_es / "Figure3_Spatial_cFOS_occurrence.png"),
                    "--panel-g",
                    str(spatial_es / "Figure3_Spatial_cFOS_NPY_occurrence.png"),
                    "--figure-name",
                    f"{FIGURE_STEM}_spanish",
                    "--dpi",
                    str(FIGURE_DPI),
                    "--panel-a-dpi",
                    str(FIGURE_DPI),
                    "--seed",
                    "31",
                ],
            ),
            (
                "Assemble English master and bilingual A-G deliverables",
                [
                    utility_py,
                    str(bundle / "assemble_final.py"),
                    "--english-publication",
                    str(publication_en),
                    "--spanish-publication",
                    str(publication_es),
                    "--output",
                    str(final),
                    "--work-dir",
                    str(work),
                ],
            ),
        ]

    def _validate_figure_candidate(self, build_dir: Path) -> dict[str, Any]:
        final = build_dir / "final"
        if final.is_symlink() or not final.is_dir():
            raise RuntimeError("Fresh Figure 3 final candidate is missing/non-normal")
        expected = {f"{FIGURE_STEM}.pdf", f"{FIGURE_STEM}.png"}
        for letter in "ABCDEFG":
            for language in ("", "_spanish"):
                for suffix in ("pdf", "png"):
                    expected.add(
                        f"panels/{FIGURE_STEM}_Panel_{letter}{language}.{suffix}"
                    )
        actual = {
            path.relative_to(final).as_posix()
            for path in final.rglob("*")
            if path.is_file() and path.suffix.lower() in {".pdf", ".png"}
        }
        if actual != expected:
            raise RuntimeError(
                f"Figure graphic inventory mismatch: missing={sorted(expected-actual)} "
                f"extra={sorted(actual-expected)}"
            )

        Image.MAX_IMAGE_PIXELS = None
        dpi_records: dict[str, list[float]] = {}
        for png in sorted(final.rglob("*.png")):
            if png.is_symlink() or not png.is_file():
                raise RuntimeError(f"Non-normal Figure PNG: {png}")
            with Image.open(png) as image:
                dpi = image.info.get("dpi", (0.0, 0.0))
                if not isinstance(dpi, (tuple, list)) or len(dpi) < 2:
                    raise RuntimeError(
                        f"Figure PNG has no two-axis DPI metadata: {png}"
                    )
                axes = [float(dpi[0]), float(dpi[1])]
                if min(axes) < 599.0:
                    raise RuntimeError(
                        f"Figure PNG is below 600-dpi tolerance: {png}: {axes}"
                    )
                dpi_records[png.relative_to(final).as_posix()] = axes

        pair_geometry: dict[str, Any] = {}
        for letter in "ABCDEFG":
            en_png = final / "panels" / f"{FIGURE_STEM}_Panel_{letter}.png"
            es_png = final / "panels" / f"{FIGURE_STEM}_Panel_{letter}_spanish.png"
            with Image.open(en_png) as en_image, Image.open(es_png) as es_image:
                if en_image.size != es_image.size:
                    raise RuntimeError(
                        f"Bilingual PNG geometry mismatch: Panel {letter}"
                    )
                pixels = list(en_image.size)
            en_pdf = final / "panels" / f"{FIGURE_STEM}_Panel_{letter}.pdf"
            es_pdf = final / "panels" / f"{FIGURE_STEM}_Panel_{letter}_spanish.pdf"
            en_page = PdfReader(str(en_pdf)).pages
            es_page = PdfReader(str(es_pdf)).pages
            if len(en_page) != 1 or len(es_page) != 1:
                raise RuntimeError(
                    f"Bilingual panel PDF is not single-page: Panel {letter}"
                )
            en_points = [
                float(en_page[0].mediabox.width),
                float(en_page[0].mediabox.height),
            ]
            es_points = [
                float(es_page[0].mediabox.width),
                float(es_page[0].mediabox.height),
            ]
            if (
                max(abs(en_points[0] - es_points[0]), abs(en_points[1] - es_points[1]))
                > 1e-3
            ):
                raise RuntimeError(f"Bilingual PDF geometry mismatch: Panel {letter}")
            pair_geometry[letter] = {"pixels": pixels, "pdf_points": en_points}

        raster_rows = 0
        minimum_pdf_raster_ppi: float | None = None
        for pdf in sorted(final.rglob("*.pdf")):
            if pdf.is_symlink() or not pdf.is_file():
                raise RuntimeError(f"Non-normal Figure PDF: {pdf}")
            if len(PdfReader(str(pdf)).pages) != 1:
                raise RuntimeError(f"Figure PDF is not single-page: {pdf}")
            fonts = subprocess.run(
                ["pdffonts", str(pdf)],
                check=True,
                text=True,
                stdin=subprocess.DEVNULL,
                capture_output=True,
            ).stdout
            if "Type 3" in fonts:
                raise RuntimeError(f"Type 3 font found in Figure PDF: {pdf}")
            for row in fonts.splitlines()[2:]:
                if row.strip() and " yes " not in f" {row} ":
                    raise RuntimeError(f"Unembedded font in Figure PDF: {pdf}: {row}")
            images = subprocess.run(
                ["pdfimages", "-list", str(pdf)],
                check=True,
                text=True,
                stdin=subprocess.DEVNULL,
                capture_output=True,
            ).stdout
            for row in images.splitlines()[2:]:
                fields = row.split()
                if len(fields) < 14 or not fields[0].isdigit():
                    continue
                try:
                    x_ppi, y_ppi = float(fields[12]), float(fields[13])
                except ValueError as exc:
                    raise RuntimeError(
                        f"Cannot parse embedded raster PPI: {pdf}: {row}"
                    ) from exc
                minimum_pdf_raster_ppi = (
                    min(x_ppi, y_ppi)
                    if minimum_pdf_raster_ppi is None
                    else min(minimum_pdf_raster_ppi, x_ppi, y_ppi)
                )
                raster_rows += 1
                if min(x_ppi, y_ppi) < 599.0:
                    raise RuntimeError(
                        f"Embedded PDF raster is below 600-dpi tolerance: {pdf}: "
                        f"{x_ppi} x {y_ppi} ppi"
                    )

        spanish_png = build_dir / "spanish_publication" / f"{FIGURE_STEM}_spanish.png"
        spanish_pdf = build_dir / "spanish_publication" / f"{FIGURE_STEM}_spanish.pdf"
        with Image.open(spanish_png) as image:
            dpi = image.info.get("dpi", (0.0, 0.0))
            if min(float(dpi[0]), float(dpi[1])) < 599.0:
                raise RuntimeError("Spanish Figure 3 master is below 600-dpi tolerance")
        if len(PdfReader(str(spanish_pdf)).pages) != 1:
            raise RuntimeError("Spanish Figure 3 master PDF is not single-page")

        return {
            "status": "PASS",
            "graphics": len(actual),
            "expected_graphics": 30,
            "png_dpi_floor_tolerance": 599.0,
            "png_dpi_records": dpi_records,
            "bilingual_geometry_pairs": "7/7",
            "pair_geometry": pair_geometry,
            "pdf_single_page_embedded_fonts_no_type3": True,
            "embedded_pdf_raster_rows": raster_rows,
            "minimum_pdf_raster_ppi": minimum_pdf_raster_ppi,
            "spanish_master_600_dpi": True,
        }

    def start_figure_build(self, payload: Mapping[str, Any]) -> dict[str, Any]:
        reviewer = str(payload.get("reviewer", "")).strip()
        session = str(payload.get("session", "")).strip()
        if not reviewer or not session:
            raise ValueError("reviewer name and session identifier are required")
        if max(len(reviewer), len(session)) > 1000:
            raise ValueError("review fields are too long")
        with self._build_lock:
            if self._build["running"]:
                raise RuntimeError("A Figure 3 build is already running")
            self._build = {
                "schema": "fig3_hil_figure_build_status_v1",
                "status": "snapshotting",
                "running": True,
                "ok": None,
                "step": "Snapshotting both accepted HIL pairs…",
                "completed_steps": 0,
                "total_steps": 9,
                "build_dir": None,
                "started_at_utc": datetime.now(timezone.utc).isoformat(),
                "finished_at_utc": None,
                "error": None,
                "log": "",
                "outputs": [],
            }
        try:
            build_dir = self._new_build_dir()
            self._set_build(build_dir=str(build_dir))
            snapshot, snapshot_receipt = self._snapshot_hil(
                build_dir, reviewer, session
            )
            commands = self._figure_commands(build_dir, snapshot)
        except Exception as exc:
            self._set_build(
                status="failed",
                running=False,
                ok=False,
                step="Build could not start.",
                error=str(exc),
                finished_at_utc=datetime.now(timezone.utc).isoformat(),
            )
            raise
        worker = threading.Thread(
            target=self._run_figure_build,
            args=(build_dir, snapshot_receipt, commands, reviewer, session),
            name=f"fig3-build-{build_dir.name}",
            daemon=True,
        )
        worker.start()
        return self.figure_status()

    def _run_figure_build(
        self,
        build_dir: Path,
        snapshot_receipt: Mapping[str, Any],
        commands: Sequence[tuple[str, list[str]]],
        reviewer: str,
        session: str,
    ) -> None:
        command_records: list[dict[str, Any]] = []
        environment = os.environ.copy()
        for name in (
            "PYTHONPATH",
            "PYTHONHOME",
            "BASH_ENV",
            "ENV",
            "CDPATH",
            "GLOBIGNORE",
        ):
            environment.pop(name, None)
        environment.update(
            {
                "PYTHONDONTWRITEBYTECODE": "1",
                "MPLBACKEND": "Agg",
            }
        )
        try:
            self._set_build(status="running")
            for index, (label_text, command) in enumerate(commands, start=1):
                self._set_build(
                    step=f"{index}/9 · {label_text}", completed_steps=index - 1
                )
                self._append_build_log(f"\n[{index}/9] {label_text}\n")
                started = datetime.now(timezone.utc)
                result = subprocess.run(
                    command,
                    cwd=self.root,
                    env=environment,
                    check=False,
                    text=True,
                    stdin=subprocess.DEVNULL,
                    capture_output=True,
                )
                finished = datetime.now(timezone.utc)
                combined = result.stdout + result.stderr
                self._append_build_log(combined)
                command_records.append(
                    {
                        "step": index,
                        "label": label_text,
                        "command": command,
                        "started_at_utc": started.isoformat(),
                        "finished_at_utc": finished.isoformat(),
                        "returncode": result.returncode,
                        "stdout_sha256": sha256_bytes(result.stdout.encode("utf-8")),
                        "stderr_sha256": sha256_bytes(result.stderr.encode("utf-8")),
                    }
                )
                if result.returncode != 0:
                    raise RuntimeError(
                        f"Step {index} failed ({label_text}); see BUILD_LOG.txt"
                    )
                self._set_build(completed_steps=index)

            self._set_build(step="9/9 · Validate 600-dpi bilingual deliverables")
            self._append_build_log("\n[9/9] Validate 600-dpi bilingual deliverables\n")
            validation = self._validate_figure_candidate(build_dir)
            final = build_dir / "final"
            spanish = build_dir / "spanish_publication"
            output_paths = [
                ("English master PDF", final / f"{FIGURE_STEM}.pdf"),
                ("English master PNG", final / f"{FIGURE_STEM}.png"),
                ("Spanish master PDF", spanish / f"{FIGURE_STEM}_spanish.pdf"),
                ("Spanish master PNG", spanish / f"{FIGURE_STEM}_spanish.png"),
            ]
            outputs = []
            for label_text, path in output_paths:
                relative = path.relative_to(self.figure_build_root).as_posix()
                outputs.append(
                    {
                        "label": label_text,
                        "relative": relative,
                        "path": str(path),
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                )
            receipt = {
                "schema": "fig3_hil_button_figure_build_v1",
                "status": "PASS",
                "created_at_utc": datetime.now(timezone.utc).isoformat(),
                "reviewer": reviewer,
                "session": session,
                "build_dir": str(build_dir),
                "figure_dpi": FIGURE_DPI,
                "review_candidate_only": True,
                "promoted": False,
                "certified_closed_world_raw_project_validation": False,
                "hil_input_snapshot": snapshot_receipt,
                "commands": command_records,
                "validation": validation,
                "outputs": outputs,
            }
            immutable_write(
                build_dir / "FIGURE_BUILD_RECEIPT.json",
                (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode("utf-8"),
            )
            receipt_relative = (
                (build_dir / "FIGURE_BUILD_RECEIPT.json")
                .relative_to(self.figure_build_root)
                .as_posix()
            )
            outputs.append(
                {
                    "label": "Build validation receipt",
                    "relative": receipt_relative,
                    "path": str(build_dir / "FIGURE_BUILD_RECEIPT.json"),
                    "bytes": (build_dir / "FIGURE_BUILD_RECEIPT.json").stat().st_size,
                    "sha256": sha256_file(build_dir / "FIGURE_BUILD_RECEIPT.json"),
                }
            )
            self._append_build_log(
                "[PASS] 30 final graphics; 7/7 bilingual pairs; 600 dpi.\n"
            )
            self._set_build(
                status="passed",
                running=False,
                ok=True,
                step="Figure 3 candidate passed all build checks.",
                completed_steps=9,
                outputs=outputs,
                finished_at_utc=datetime.now(timezone.utc).isoformat(),
            )
        except Exception as exc:
            self._append_build_log(f"\n[FAIL] {exc}\n")
            self._set_build(
                status="failed",
                running=False,
                ok=False,
                step="Figure 3 build failed; accepted HIL files were not changed.",
                error=str(exc),
                finished_at_utc=datetime.now(timezone.utc).isoformat(),
            )

    def figure_file(self, relative: str) -> Path:
        candidate = Path(relative)
        if not relative or candidate.is_absolute() or ".." in candidate.parts:
            raise ValueError("invalid figure file path")
        resolved = (self.figure_build_root / candidate).resolve()
        try:
            resolved.relative_to(self.figure_build_root)
        except ValueError as exc:
            raise ValueError("figure file path escapes the build root") from exc
        if resolved.is_symlink() or not resolved.is_file():
            raise ValueError("figure file is missing/non-normal")
        return resolved


class Server(ThreadingHTTPServer):
    state: ReviewState
    daemon_threads = True
    allow_reuse_address = True


class Handler(BaseHTTPRequestHandler):
    server: Server

    def log_message(self, fmt: str, *args: object) -> None:
        sys.stderr.write("[fig3-hil] " + fmt % args + "\n")

    def send_payload(self, status: int, payload: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(payload)

    def send_json(self, status: int, payload: Mapping[str, Any]) -> None:
        self.send_payload(
            status,
            json.dumps(payload, separators=(",", ":")).encode("utf-8"),
            "application/json; charset=utf-8",
        )

    def send_file(self, path: Path, content_type: str) -> None:
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(path.stat().st_size))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        with path.open("rb") as handle:
            shutil.copyfileobj(handle, self.wfile, 1024 * 1024)

    def condition(self, query: Mapping[str, list[str]]) -> str:
        condition = str(query.get("condition", [""])[0])
        if condition not in CONDITIONS:
            raise ValueError("condition must be Water or Allulose")
        return condition

    def do_GET(self) -> None:
        try:
            parsed = urlparse(self.path)
            query = parse_qs(parsed.query)
            if parsed.path == "/":
                self.send_payload(200, HTML.encode("utf-8"), "text/html; charset=utf-8")
            elif parsed.path == "/api/config":
                self.send_json(200, self.server.state.config())
            elif parsed.path == "/api/annotation":
                self.send_json(200, self.server.state.annotation(self.condition(query)))
            elif parsed.path == "/api/figure-status":
                self.send_json(200, self.server.state.figure_status())
            elif parsed.path == "/api/figure-file":
                relative = str(query.get("relative", [""])[0])
                path = self.server.state.figure_file(relative)
                content_type = (
                    mimetypes.guess_type(path.name)[0] or "application/octet-stream"
                )
                self.send_file(path, content_type)
            elif parsed.path in {"/api/image", "/api/source-overlay"}:
                item = self.server.state.samples[self.condition(query)]
                path = (
                    item["background_path"]
                    if parsed.path == "/api/image"
                    else item["overlay_path"]
                )
                self.send_payload(200, path.read_bytes(), "image/png")
            else:
                self.send_json(404, {"error": "not found"})
        except Exception as exc:
            self.send_json(400, {"error": str(exc)})

    def do_POST(self) -> None:
        try:
            route = urlparse(self.path).path
            if route not in {"/api/accept", "/api/make-figure"}:
                self.send_json(404, {"error": "not found"})
                return
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 2_000_000:
                raise ValueError("invalid JSON body size")
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
            if not isinstance(payload, dict):
                raise ValueError("JSON body must be an object")
            if route == "/api/accept":
                result = self.server.state.accept(payload)
            else:
                result = self.server.state.start_figure_build(payload)
            self.send_json(200, result)
        except Exception as exc:
            self.send_json(400, {"error": str(exc)})


HTML = r"""<!doctype html>
<html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Figure 3 ventricle HIL</title>
<style>
:root{--ink:#18222d;--muted:#65727c;--paper:#f4f5f2;--rule:#cad1cf;--ok:#21734b;--bad:#a13b35;--accent:#6b3fb3}*{box-sizing:border-box}html,body{height:100%;margin:0}body{font:14px/1.35 Arial,sans-serif;color:var(--ink);background:var(--paper);overflow:hidden}button,input,select,textarea{font:inherit}button{border:1px solid #aab4b2;background:#fff;padding:7px 9px;border-radius:3px;cursor:pointer;font-weight:700}button:hover{background:#edf0ee}button:disabled{cursor:not-allowed;opacity:.48}.app{height:100%;display:grid;grid-template-rows:56px minmax(0,1fr)}header{display:flex;align-items:center;gap:16px;background:#fff;border-bottom:1px solid var(--rule);padding:0 16px}h1{font:700 18px Georgia,serif;margin:0}.sub{color:var(--muted);font-size:12px}#progress{margin-left:auto;font-weight:800}.main{min-height:0;display:grid;grid-template-columns:330px minmax(0,1fr)}aside{overflow:auto;background:#fff;border-right:1px solid var(--rule);padding:13px}.group{border-top:1px solid var(--rule);padding:11px 0}.group:first-child{border:0;padding-top:0}.title{font-size:11px;text-transform:uppercase;letter-spacing:.09em;color:var(--muted);font-weight:800;margin-bottom:6px}.field{display:block;font-size:12px;font-weight:700;margin:7px 0 3px}input,select,textarea{width:100%;padding:7px;border:1px solid #aab4b2;border-radius:3px}textarea{height:65px;resize:vertical}.grid{display:grid;grid-template-columns:1fr 1fr;gap:7px;margin-top:7px}.accept{width:100%;padding:11px;color:#fff;background:var(--ok);border-color:var(--ok)}.makefig{width:100%;padding:12px;color:#fff;background:var(--accent);border-color:var(--accent)}.status{min-height:42px;margin-top:8px;color:var(--muted);overflow-wrap:anywhere}.status.ok{color:var(--ok);font-weight:800}.status.error{color:var(--bad);font-weight:800}.buildlog{display:none;max-height:135px;overflow:auto;white-space:pre-wrap;background:#151b20;color:#dbe5e2;padding:7px;border-radius:3px;font:10px/1.35 monospace;margin:7px 0 0}.figurelinks a{display:block;margin-top:5px;color:var(--accent);font-weight:800;overflow-wrap:anywhere}.workspace{min-width:0;min-height:0;display:grid;grid-template-rows:42px minmax(0,1fr);background:#101619}.toolbar{display:flex;align-items:center;gap:7px;background:#edf0ee;padding:5px 10px}.name{margin-right:auto;font-weight:800}.wrap{position:relative;min-height:0;overflow:hidden}canvas{position:absolute;inset:0;width:100%;height:100%;touch-action:none;cursor:crosshair}.hint{position:absolute;left:10px;bottom:9px;background:#ffffffe8;padding:6px 8px;font-size:11px;pointer-events:none}.accepted{color:var(--ok);font-weight:800}.dirty{color:var(--bad);font-weight:800}@media(max-width:850px){.main{grid-template-columns:275px minmax(0,1fr)}.sub{display:none}}
</style></head><body><div class="app"><header><h1>Figure 3 ventricle HIL</h1><span class="sub">Draw one exact closed lumen polygon for Water and Allulose</span><span id="progress">0 / 2 accepted</span></header><div class="main"><aside>
<div class="group"><div class="title">Review identity</div><label class="field">Reviewer</label><input id="reviewer" maxlength="200"><label class="field">Session</label><input id="session" maxlength="200" placeholder="e.g. final ventricle pass"><label class="field">Notes</label><textarea id="notes"></textarea></div>
<div class="group"><div class="title">Section</div><select id="condition"></select><label class="field"><input id="showSource" type="checkbox" checked style="width:auto"> Show raw-DAPI proposal</label><div class="sub">The magenta proposal is guidance only. The accepted polygon fully controls the white lumen and black wall in all three cartoon views.</div></div>
<div class="group"><div class="title">Polygon</div><div class="grid"><button id="undo">Undo point</button><button id="reopen">Reopen</button><button id="close">Close polygon</button><button id="clear">Clear canvas</button></div><div id="vertexCount" class="sub" style="margin-top:8px">0 vertices</div></div>
<div class="group"><button id="accept" class="accept">Accept this HIL polygon</button><div id="status" class="status">Loading…</div></div>
<div class="group"><div class="title">Final render</div><button id="makeFigure" class="makefig" disabled>Make Figure 3 · 600 dpi</button><div id="buildStatus" class="status">Accept both HIL polygons to enable rendering.</div><div id="figureLinks" class="figurelinks"></div><pre id="buildLog" class="buildlog"></pre></div>
</aside><section class="workspace"><div class="toolbar"><span id="sectionName" class="name">No section</span><span id="recordState"></span><button id="zoomOut">−</button><span id="zoomReadout">100%</span><button id="zoomIn">+</button><button id="fit">Fit</button></div><div id="wrap" class="wrap"><canvas id="canvas"></canvas><div class="hint">Click: vertex | Drag vertex: refine | Enter/double-click: close | Space/middle-drag: pan | Wheel: zoom</div></div></section></div></div>
<script>'use strict';const $=id=>document.getElementById(id),canvas=$('canvas'),wrap=$('wrap'),ctx=canvas.getContext('2d');const st={sections:[],condition:'Water',image:null,source:null,vertices:[],closed:false,scale:1,tx:0,ty:0,space:false,panning:false,drag:-1,px:0,py:0,dirty:false,loading:false,canBuild:false,buildRunning:false,poll:null};
function endpoint(path){return `${path}?condition=${encodeURIComponent(st.condition)}`}function status(msg,kind=''){$('status').textContent=msg;$('status').className=`status ${kind}`}function loadImage(url){return new Promise((ok,bad)=>{let im=new Image();im.onload=()=>ok(im);im.onerror=()=>bad(Error(`Cannot load ${url}`));im.src=url+`&v=${Date.now()}`})}function refresh(){let n=st.vertices.length;$('vertexCount').textContent=`${n} vertices · ${st.closed?'closed':'open'}`;$('recordState').textContent=st.dirty?'UNSAVED':st.closed?'ACCEPTED/LOADED':'PENDING';$('recordState').className=st.dirty?'dirty':st.closed?'accepted':'';draw()}
function setMakeButton(){let button=$('makeFigure');button.disabled=!st.canBuild||st.buildRunning;button.textContent=st.buildRunning?'Making Figure 3…':'Make Figure 3 · 600 dpi'}
function showBuild(build){st.buildRunning=!!build.running;setMakeButton();let box=$('buildStatus'),kind=build.ok===true?'ok':build.ok===false?'error':'';box.className=`status ${kind}`;let progress=build.running?` (${build.completed_steps}/${build.total_steps})`:'';let message=build.status==='idle'&&st.canBuild?'Ready: both HIL polygons are accepted.':(build.step||build.status||'No build started.');box.textContent=message+progress+(build.error?` — ${build.error}`:'');let log=$('buildLog');log.textContent=build.log||'';log.style.display=build.log?'block':'none';if(build.log)log.scrollTop=log.scrollHeight;let links=$('figureLinks');links.replaceChildren();for(let item of (build.outputs||[])){let a=document.createElement('a');a.textContent=`Open ${item.label}`;a.href=`/api/figure-file?relative=${encodeURIComponent(item.relative)}`;a.target='_blank';a.rel='noopener';links.appendChild(a)}}
async function pollBuild(){try{let r=await fetch('/api/figure-status'),j=await r.json();if(!r.ok)throw Error(j.error);showBuild(j)}catch(e){let box=$('buildStatus');box.textContent=`Cannot read build status: ${e.message}`;box.className='status error'}finally{clearTimeout(st.poll);st.poll=setTimeout(pollBuild,2500)}}
function resize(){let d=devicePixelRatio||1,r=wrap.getBoundingClientRect();canvas.width=Math.max(1,Math.round(r.width*d));canvas.height=Math.max(1,Math.round(r.height*d));draw()}function fit(){if(!st.image)return;let w=canvas.clientWidth,h=canvas.clientHeight;st.scale=Math.min(w/st.image.width,h/st.image.height)*.96;st.tx=(w-st.image.width*st.scale)/2;st.ty=(h-st.image.height*st.scale)/2;draw()}function zoom(f,cx=canvas.clientWidth/2,cy=canvas.clientHeight/2){let old=st.scale;st.scale=Math.max(.03,Math.min(30,old*f));let r=st.scale/old;st.tx=cx-(cx-st.tx)*r;st.ty=cy-(cy-st.ty)*r;draw()}
function draw(){let d=devicePixelRatio||1;ctx.setTransform(d,0,0,d,0,0);ctx.fillStyle='#101619';ctx.fillRect(0,0,canvas.clientWidth,canvas.clientHeight);if(!st.image)return;ctx.save();ctx.translate(st.tx,st.ty);ctx.scale(st.scale,st.scale);ctx.drawImage(st.image,0,0);if(st.source&&$('showSource').checked){ctx.globalAlpha=.75;ctx.drawImage(st.source,0,0);ctx.globalAlpha=1}if(st.vertices.length){ctx.beginPath();ctx.moveTo(...st.vertices[0]);for(let i=1;i<st.vertices.length;i++)ctx.lineTo(...st.vertices[i]);if(st.closed)ctx.closePath();if(st.closed){ctx.fillStyle='rgba(255,255,255,.66)';ctx.fill()}ctx.strokeStyle='#6b3fb3';ctx.lineWidth=4.5/st.scale;ctx.lineJoin='round';ctx.stroke();for(let i=0;i<st.vertices.length;i++){let p=st.vertices[i];ctx.beginPath();ctx.arc(p[0],p[1],5/st.scale,0,7);ctx.fillStyle=i===st.drag?'#ffda45':'#fff';ctx.fill();ctx.lineWidth=2/st.scale;ctx.strokeStyle='#6b3fb3';ctx.stroke()}}ctx.restore();$('zoomReadout').textContent=`${Math.round(st.scale*100)}%`}
function point(e){let r=canvas.getBoundingClientRect();return[(e.clientX-r.left-st.tx)/st.scale,(e.clientY-r.top-st.ty)/st.scale]}function inside(p){return st.image&&p[0]>=0&&p[1]>=0&&p[0]<st.image.width&&p[1]<st.image.height}function mark(){st.dirty=true;refresh()}function closePoly(){if(st.vertices.length<3){status('A polygon needs at least three vertices.','error');return}st.closed=true;mark();status('Polygon closed. Drag vertices to refine, then accept.')}
async function load(condition,force=false){if(st.loading)return;if(!force&&st.dirty&&!confirm('Discard unsaved canvas edits?')){$('condition').value=st.condition;return}st.loading=true;st.condition=condition;st.image=null;st.source=null;st.vertices=[];st.closed=false;st.dirty=false;status('Loading registered field…');try{let [rec,img,src]=await Promise.all([fetch(endpoint('/api/annotation')).then(async r=>{let j=await r.json();if(!r.ok)throw Error(j.error);return j}),loadImage(endpoint('/api/image')),loadImage(endpoint('/api/source-overlay'))]);st.image=img;st.source=src;if(rec.found){st.vertices=rec.polygon_vertices_xy.map(p=>[+p[0],+p[1]]);st.closed=true;if(rec.reviewer&&!$('reviewer').value)$('reviewer').value=rec.reviewer;if(rec.session&&!$('session').value)$('session').value=rec.session;$('notes').value=rec.notes||'';status('Accepted polygon loaded. Editing and re-accepting archives the previous pair.','ok')}else status('Draw the complete ventricular lumen as one closed polygon.');let s=st.sections.find(x=>x.condition===condition);$('sectionName').textContent=`${s.condition} · ${s.animal_id} · ${s.width} × ${s.height}`;refresh();fit()}catch(e){status(e.message,'error')}finally{st.loading=false}}
canvas.oncontextmenu=e=>e.preventDefault();canvas.onpointerdown=e=>{if(!st.image)return;if(e.button===1||e.button===2||st.space){st.panning=true;st.px=e.clientX;st.py=e.clientY;canvas.setPointerCapture(e.pointerId);canvas.style.cursor='grabbing';return}if(e.button)return;let p=point(e);if(!inside(p))return;let nearest=-1,best=12/st.scale;for(let i=0;i<st.vertices.length;i++){let q=st.vertices[i],d=Math.hypot(q[0]-p[0],q[1]-p[1]);if(d<best){best=d;nearest=i}}if(nearest>=0){st.drag=nearest;canvas.setPointerCapture(e.pointerId);draw();return}if(st.closed){status('Use Reopen before adding new points, or drag an existing vertex.','error');return}st.vertices.push(p.map(v=>+v.toFixed(3)));mark()};canvas.onpointermove=e=>{if(st.panning){st.tx+=e.clientX-st.px;st.ty+=e.clientY-st.py;st.px=e.clientX;st.py=e.clientY;draw();return}if(st.drag>=0){let p=point(e);if(inside(p)){st.vertices[st.drag]=p.map(v=>+v.toFixed(3));mark()}}};function stop(e){st.panning=false;st.drag=-1;canvas.style.cursor='crosshair';try{canvas.releasePointerCapture(e.pointerId)}catch(_){ }draw()}canvas.onpointerup=stop;canvas.onpointercancel=stop;canvas.ondblclick=e=>{e.preventDefault();closePoly()};canvas.addEventListener('wheel',e=>{e.preventDefault();let r=canvas.getBoundingClientRect();zoom(e.deltaY<0?1.15:1/1.15,e.clientX-r.left,e.clientY-r.top)},{passive:false});window.onkeydown=e=>{if(['INPUT','SELECT','TEXTAREA'].includes(document.activeElement.tagName))return;if(e.code==='Space'){st.space=true;e.preventDefault()}if(e.key==='Enter'){closePoly();e.preventDefault()}if((e.ctrlKey||e.metaKey)&&e.key.toLowerCase()==='z'){$('undo').click();e.preventDefault()}};window.onkeyup=e=>{if(e.code==='Space')st.space=false};
$('undo').onclick=()=>{if(!st.vertices.length){status('No vertex to undo.','error');return}st.closed=false;st.vertices.pop();mark()};$('reopen').onclick=()=>{if(!st.closed){status('Polygon is already open.');return}st.closed=false;mark();status('Polygon reopened; add or remove vertices, then close it again.')};$('close').onclick=closePoly;$('clear').onclick=()=>{if(st.vertices.length&&confirm('Clear the browser canvas? Accepted disk files are not deleted.')){st.vertices=[];st.closed=false;mark();status('Canvas cleared; accepted disk files remain until a new polygon is accepted.')}};$('showSource').onchange=draw;$('zoomIn').onclick=()=>zoom(1.25);$('zoomOut').onclick=()=>zoom(.8);$('fit').onclick=fit;$('condition').onchange=e=>load(e.target.value);for(let id of ['reviewer','session']){$(id).value=localStorage.getItem(`fig3vent_${id}`)||'';$(id).oninput=e=>localStorage.setItem(`fig3vent_${id}`,e.target.value)}
$('accept').onclick=async()=>{let reviewer=$('reviewer').value.trim(),session=$('session').value.trim();if(!reviewer||!session){status('Reviewer and session are required.','error');return}if(!st.closed||st.vertices.length<3){status('Close one polygon before accepting.','error');return}if(!confirm(`Accept ${st.condition} ventricle with ${st.vertices.length} vertices?`))return;status('Saving and hash-validating…');try{let r=await fetch('/api/accept',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({condition:st.condition,reviewer,session,notes:$('notes').value,polygon_vertices_xy:st.vertices})}),j=await r.json();if(!r.ok)throw Error(j.error);st.dirty=false;status(`Accepted and validated. Mask SHA ${j.mask_sha256.slice(0,16)}…`,'ok');await bootstrap(st.condition)}catch(e){status(e.message,'error')}};
$('makeFigure').onclick=async()=>{let reviewer=$('reviewer').value.trim(),session=$('session').value.trim();if(!reviewer||!session){let box=$('buildStatus');box.textContent='Reviewer and session are required before rendering.';box.className='status error';return}if(st.dirty&&!confirm('The canvas has unsaved edits. Make the figure from the two accepted disk polygons anyway?'))return;if(!confirm('Create a new versioned 600-dpi English/Spanish Figure 3 candidate? The promoted figure will not be overwritten.'))return;st.buildRunning=true;setMakeButton();let box=$('buildStatus');box.textContent='Starting a source-bound Figure 3 build…';box.className='status';try{let r=await fetch('/api/make-figure',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({reviewer,session})}),j=await r.json();if(!r.ok)throw Error(j.error);showBuild(j)}catch(e){st.buildRunning=false;setMakeButton();box.textContent=e.message;box.className='status error'}};
async function bootstrap(keep='Water'){try{let r=await fetch('/api/config'),j=await r.json();if(!r.ok)throw Error(j.error);st.sections=j.sections;st.canBuild=j.accepted===j.total;$('progress').textContent=`${j.accepted} / ${j.total} accepted`;$('condition').innerHTML=j.sections.map(s=>`<option value="${s.condition}">${s.condition}${s.accepted?' ✓':''}</option>`).join('');$('condition').value=keep;setMakeButton();await load(keep,true)}catch(e){status(e.message,'error')}}window.addEventListener('resize',resize);window.addEventListener('beforeunload',()=>clearTimeout(st.poll));resize();bootstrap().then(pollBuild);
</script></body></html>"""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--figure-build-root",
        type=Path,
        default=None,
        help="Fresh versioned candidates are created here (default: OUTPUT_DIR/figure_builds).",
    )
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT)
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = discover_paper_root(args.root)
    output_dir = (
        args.output_dir.expanduser().resolve()
        if args.output_dir is not None
        else (root / HIL_RELATIVE).resolve()
    )
    state = ReviewState(root, output_dir, args.figure_build_root)
    print(f"review_output\t{output_dir}", flush=True)
    print(f"figure_build_root\t{state.figure_build_root}", flush=True)
    if args.prepare_only:
        print(json.dumps(state.config(), indent=2, sort_keys=True), flush=True)
        return 0
    if not 1 <= args.port <= 65535:
        raise ValueError("port must be in 1..65535")
    server = Server((args.host, args.port), Handler)
    server.state = state
    print(f"review_url\thttp://{args.host}:{args.port}/", flush=True)
    print(
        "No accepted HIL file is deleted; replacements are archived under history/.",
        flush=True,
    )
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
