#!/usr/bin/env python3
"""Shared, dependency-light helpers for the Figure S3 WSI workflow."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping


FIGS3_DIR = Path(__file__).resolve().parent
PAPER_DIR = FIGS3_DIR.parent
DEFAULT_SOURCE_ROOT = Path(
    "/media/server/STORAGE/Motic_AnatomiaPatologica_2025/"
    "HE_Liver_Spleen_Kidney_Nancy_2026"
)
SOURCE_MPP_DEFAULT = 0.261780
ORGANS = ("Kidney", "Liver", "Spleen")
TREATMENTS = ("Water", "Sucrose", "Allulose")
SAMPLE_RE = re.compile(r"^m26-(\d{3})$")


class FigS3Error(RuntimeError):
    """Raised when an invariant of the Figure S3 workflow is violated."""


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def atomic_json(path: Path, payload: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def atomic_csv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Iterable[str] | None = None,
) -> None:
    rows = list(rows)
    if fieldnames is None:
        if not rows:
            raise FigS3Error(f"Cannot infer CSV columns for empty table: {path}")
        fieldnames = list(rows[0])
    fieldnames = list(fieldnames)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def validate_sample_id(sample_id: str) -> str:
    value = str(sample_id).strip().lower().replace("_", "-")
    if not SAMPLE_RE.fullmatch(value):
        raise FigS3Error(f"Unsafe or unexpected sample ID: {sample_id!r}")
    return value


def sample_source_folder(sample_id: str) -> str:
    value = validate_sample_id(sample_id)
    return value.upper()


def canonical_json_sha256(payload: Any) -> str:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return sha256_text(text)


def require_regular_file(path: Path, description: str = "file") -> Path:
    candidate = Path(path).expanduser().absolute()
    if candidate.is_symlink() or not candidate.is_file():
        raise FigS3Error(f"Missing or unsafe {description}: {candidate}")
    return candidate.resolve()


def require_directory(path: Path, description: str = "directory") -> Path:
    candidate = Path(path).expanduser().absolute()
    if candidate.is_symlink() or not candidate.is_dir():
        raise FigS3Error(f"Missing or unsafe {description}: {candidate}")
    return candidate.resolve()


def file_identity(path: Path) -> dict[str, Any]:
    candidate = require_regular_file(path)
    stat = candidate.stat()
    return {
        "path": str(candidate),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def parse_float(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise FigS3Error(f"Invalid {name}: {value!r}") from exc
    if not (result > 0):
        raise FigS3Error(f"{name} must be greater than zero")
    return result
