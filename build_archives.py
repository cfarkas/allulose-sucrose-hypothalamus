#!/usr/bin/env python3
"""Build deterministic, manifest-bound ZIP shards for the public Paper tree.

The canonical public exclusions are those in
``scripts/utilities/01_validate_bundle.py``.  Files are assigned to exactly one
of four logical Zenodo records: software/core, main-figure source data,
supplementary-figure source data, or the full-resolution Figure S3 whole-slide
images.  A record may contain several ZIP shards; a shard is never modelled as
a separate record.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.util
import json
import os
import re
import shutil
import stat
import sys
import uuid
import zipfile
from dataclasses import dataclass, replace
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, Sequence


SCHEMA = "apotome-paper-archive-manifest-v1"
RECORDS = (
    "software_core",
    "main_figure_data",
    "supplementary_figure_data",
    "figs3_wsi",
)
RECORD_PREFIX = {
    "software_core": "paper-software-core",
    "main_figure_data": "paper-main-figure-data",
    "supplementary_figure_data": "paper-supplementary-figure-data",
    "figs3_wsi": "paper-figs3-wsi",
}
LOCAL_ONLY = {
    "analyses",
    "legacy",
    "quarantine",
    "recovery_archive",
    "recovery_support",
    "docs",
    ".claude",
    ".mypy_cache",
    "__pycache__",
}
LOCAL_ONLY_ROOT_PREFIXES = ("apotome_rebuild_",)
LOCAL_ONLY_FILES = {"README2.txt"}
LOCAL_ONLY_FILE_SUFFIXES = (".orig", ".rej")
PUBLIC_EXCLUDED_RELATIVE_PREFIXES = (
    "Fig4/figure_bundle_v3/raw_data/allen_p56_snapshot",
    "FigS3/registration",
    "literature_webscrap_30_08_2026",
)
REDISTRIBUTION_EXCLUDED_BASENAMES = {"full_text_open_access.pdf"}
MAIN_FIGURE_ROOTS = {
    "Fig1",
    "Fig2",
    "Fig3",
    "Fig4",
    "Fig5",
}
SUPPLEMENTARY_FIGURE_ROOTS = {
    "FigS1",
    "FigS2",
    "FigS3",
    "FigS4",
    "FigS5",
}
DATA_COMPONENTS = {"raw", "raw_data", "channel_tiffs"}
FIGS3_DATA_PREFIXES = (
    "FigS3/exports",
    "FigS3/histoplus",
    "FigS3/organ_segmentation",
    "FigS3/raw",
    "FigS3/raw_data",
)
FIGS3_WSI_RE = re.compile(
    r"^FigS3/exports/m26-(?:00[1-9]|01[0-9]|02[0-4])/1_L0_rgb\.tif$"
)
FIGURE_ROOT_RE = re.compile(r"^Fig(?:S)?[0-9]+$")
DEFAULT_SHARD_PAYLOAD_BYTES = 43_000_000_000
DEFAULT_MAX_ARCHIVE_BYTES = 45_000_000_000
ZENODO_INDIVIDUAL_FILE_LIMIT = 50_000_000_000
ZENODO_MAX_FILES_PER_RECORD = 100
ZENODO_MAX_RECORD_BYTES = 200_000_000_000
ZENODO_DEFAULT_RECORD_BYTES = 50_000_000_000
ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE = 150_000_000_000
ZIP_DATETIME = (1980, 1, 1, 0, 0, 0)
CHUNK = 8 * 1024 * 1024
SHA_RE = re.compile(r"^[0-9a-f]{64}$")


class ArchiveError(RuntimeError):
    """Fail-closed archive contract violation."""


@dataclass(frozen=True)
class DirectoryEntry:
    path: str
    mode: int
    empty: bool = False


@dataclass(frozen=True)
class FileEntry:
    path: str
    source: Path
    mode: int
    size: int
    sha256: str
    record: str
    shard: str = ""


@dataclass(frozen=True)
class Inventory:
    root: Path
    root_mode: int
    directories: tuple[DirectoryEntry, ...]
    files: tuple[FileEntry, ...]


@dataclass(frozen=True)
class Shard:
    record: str
    index: int
    name: str
    files: tuple[FileEntry, ...]
    payload_bytes: int


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_text(mode: int) -> str:
    return f"{mode & 0o7777:04o}"


def safe_relative(value: str) -> str:
    if not value or "\\" in value or "\x00" in value:
        raise ArchiveError(f"Unsafe relative path: {value!r}")
    pure = PurePosixPath(value)
    if pure.is_absolute() or any(part in {"", ".", ".."} for part in pure.parts):
        raise ArchiveError(f"Unsafe relative path: {value!r}")
    if pure.as_posix() != value:
        raise ArchiveError(f"Non-normalized relative path: {value!r}")
    return value


def is_public_relative_path(relative: str, *, is_directory: bool) -> bool:
    """Return whether one normalized Paper-relative path is public."""
    safe_relative(relative)
    parts = PurePosixPath(relative).parts
    directory_parts = parts if is_directory else parts[:-1]
    if any(part in LOCAL_ONLY for part in directory_parts):
        return False
    if (
        is_directory
        and len(parts) == 1
        and any(parts[0].startswith(prefix) for prefix in LOCAL_ONLY_ROOT_PREFIXES)
    ):
        return False
    if not is_directory:
        if len(parts) == 1 and parts[0] in LOCAL_ONLY_FILES:
            return False
        if parts[-1].endswith(LOCAL_ONLY_FILE_SUFFIXES):
            return False
        if parts[-1] in REDISTRIBUTION_EXCLUDED_BASENAMES:
            return False
    return not any(
        relative == prefix or relative.startswith(prefix + "/")
        for prefix in PUBLIC_EXCLUDED_RELATIVE_PREFIXES
    )


def public_walk(root: Path):
    """Yield the same paths selected by canonical bundle_entries()."""
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        base = Path(directory)
        dirnames[:] = [
            name
            for name in dirnames
            if is_public_relative_path(
                (base / name).relative_to(root).as_posix(), is_directory=True
            )
        ]
        for name in dirnames:
            yield base / name
        for name in filenames:
            relative = (base / name).relative_to(root).as_posix()
            if not is_public_relative_path(relative, is_directory=False):
                continue
            yield base / name


def assert_canonical_exclusions(contract_script: Path) -> str:
    """Fail if the copied selection constants drift from the canonical script."""
    if not contract_script.is_file() or contract_script.is_symlink():
        raise ArchiveError(
            f"Canonical bundle validator is missing or linked: {contract_script}"
        )
    spec = importlib.util.spec_from_file_location(
        f"paper_bundle_contract_{uuid.uuid4().hex}",
        contract_script,
    )
    if spec is None or spec.loader is None:
        raise ArchiveError(
            f"Cannot import canonical bundle validator: {contract_script}"
        )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    observed = (
        set(module.LOCAL_ONLY),
        tuple(module.LOCAL_ONLY_ROOT_PREFIXES),
        set(module.LOCAL_ONLY_FILES),
        tuple(module.LOCAL_ONLY_FILE_SUFFIXES),
        tuple(module.PUBLIC_EXCLUDED_RELATIVE_PREFIXES),
        set(module.REDISTRIBUTION_EXCLUDED_BASENAMES),
    )
    expected = (
        LOCAL_ONLY,
        LOCAL_ONLY_ROOT_PREFIXES,
        LOCAL_ONLY_FILES,
        LOCAL_ONLY_FILE_SUFFIXES,
        PUBLIC_EXCLUDED_RELATIVE_PREFIXES,
        REDISTRIBUTION_EXCLUDED_BASENAMES,
    )
    if observed != expected:
        raise ArchiveError("Release exclusions differ from 01_validate_bundle.py")
    return sha256_file(contract_script)


def assign_record(relative: str, size: int) -> str:
    safe_relative(relative)
    parts = PurePosixPath(relative).parts
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        raise ArchiveError(f"Invalid file size for record assignment: {size!r}")
    first = parts[0]
    candidates: set[str] = set()
    if FIGS3_WSI_RE.fullmatch(relative):
        candidates.add("figs3_wsi")
    elif any(
        relative == prefix or relative.startswith(prefix + "/")
        for prefix in FIGS3_DATA_PREFIXES
    ):
        candidates.add("supplementary_figure_data")
    if first in MAIN_FIGURE_ROOTS and any(
        part in DATA_COMPONENTS for part in parts[1:-1]
    ):
        candidates.add("main_figure_data")
    if first in SUPPLEMENTARY_FIGURE_ROOTS and any(
        part in DATA_COMPONENTS for part in parts[1:-1]
    ):
        candidates.add("supplementary_figure_data")
    if (
        FIGURE_ROOT_RE.fullmatch(first)
        and first not in MAIN_FIGURE_ROOTS
        and first not in SUPPLEMENTARY_FIGURE_ROOTS
        and any(part in DATA_COMPONENTS for part in parts[1:-1])
    ):
        raise ArchiveError(f"Unknown figure data root: {first}")
    if not candidates:
        candidates.add("software_core")
    if len(candidates) != 1:
        raise ArchiveError(
            f"Duplicate logical-record assignment for {relative}: {sorted(candidates)}"
        )
    return next(iter(candidates))


def inventory_tree(
    root: Path, *, contract_script: Path | None = None
) -> tuple[Inventory, str | None]:
    lexical_root = root.expanduser().absolute()
    if lexical_root.is_symlink():
        raise ArchiveError(f"Paper root may not be a symlink: {lexical_root}")
    paper = lexical_root.resolve()
    if not paper.is_dir():
        raise ArchiveError(f"Paper root is missing or not a directory: {paper}")
    contract_hash = (
        assert_canonical_exclusions(contract_script) if contract_script else None
    )
    directories: list[DirectoryEntry] = []
    files: list[FileEntry] = []
    seen: set[str] = set()
    for path in public_walk(paper):
        relative = safe_relative(path.relative_to(paper).as_posix())
        if relative in seen:
            raise ArchiveError(f"Duplicate public-tree entry: {relative}")
        seen.add(relative)
        observed = path.lstat()
        if stat.S_ISLNK(observed.st_mode):
            raise ArchiveError(f"Symlink is forbidden in the public tree: {relative}")
        mode = stat.S_IMODE(observed.st_mode)
        if stat.S_ISDIR(observed.st_mode):
            directories.append(DirectoryEntry(relative, mode))
        elif stat.S_ISREG(observed.st_mode):
            files.append(
                FileEntry(
                    path=relative,
                    source=path,
                    mode=mode,
                    size=int(observed.st_size),
                    sha256=sha256_file(path),
                    record=assign_record(relative, int(observed.st_size)),
                )
            )
        else:
            raise ArchiveError(
                f"Non-regular public-tree entry is forbidden: {relative}"
            )
    file_parents = {str(PurePosixPath(item.path).parent) for item in files}
    directory_parents = {str(PurePosixPath(item.path).parent) for item in directories}
    directories = [
        replace(
            item,
            empty=item.path not in file_parents and item.path not in directory_parents,
        )
        for item in directories
    ]
    files.sort(key=lambda item: item.path)
    directories.sort(key=lambda item: item.path)
    validate_inventory(files)
    if len({item.path for item in directories}) != len(directories):
        raise ArchiveError("Duplicate directory inventory")
    return (
        Inventory(
            root=paper,
            root_mode=stat.S_IMODE(paper.stat().st_mode),
            directories=tuple(directories),
            files=tuple(files),
        ),
        contract_hash,
    )


def validate_inventory(files: Sequence[FileEntry]) -> None:
    paths: set[str] = set()
    for item in files:
        safe_relative(item.path)
        if item.path in paths:
            raise ArchiveError(f"File is assigned more than once: {item.path}")
        paths.add(item.path)
        if item.record not in RECORDS or not SHA_RE.fullmatch(item.sha256):
            raise ArchiveError(f"Invalid record or SHA-256 for {item.path}")
        if item.size < 0:
            raise ArchiveError(f"Negative file size for {item.path}")
        if item.record != assign_record(item.path, item.size):
            raise ArchiveError(f"Logical-record classification differs: {item.path}")


def plan_shards(
    files: Sequence[FileEntry],
    *,
    payload_limit: int = DEFAULT_SHARD_PAYLOAD_BYTES,
) -> tuple[Shard, ...]:
    validate_inventory(files)
    if payload_limit <= 0 or payload_limit > DEFAULT_SHARD_PAYLOAD_BYTES:
        raise ArchiveError(f"Invalid uncompressed shard payload limit: {payload_limit}")
    plans: list[Shard] = []
    assigned: set[str] = set()
    for record in RECORDS:
        current: list[FileEntry] = []
        payload = 0
        index = 1
        for item in sorted(
            (value for value in files if value.record == record),
            key=lambda value: value.path,
        ):
            if item.size > payload_limit:
                raise ArchiveError(
                    f"One file exceeds the shard payload target: {item.path}"
                )
            if current and payload + item.size > payload_limit:
                name = f"{RECORD_PREFIX[record]}.part{index:03d}.zip"
                plans.append(Shard(record, index, name, tuple(current), payload))
                index += 1
                current, payload = [], 0
            if item.path in assigned:
                raise ArchiveError(f"Duplicate shard assignment: {item.path}")
            assigned.add(item.path)
            current.append(item)
            payload += item.size
        if current:
            name = f"{RECORD_PREFIX[record]}.part{index:03d}.zip"
            plans.append(Shard(record, index, name, tuple(current), payload))
    if assigned != {item.path for item in files}:
        raise ArchiveError("Not every file was assigned to exactly one shard")
    for record in RECORDS:
        if (
            sum(plan.record == record for plan in plans)
            + (1 if record == "software_core" else 0)
        ) > ZENODO_MAX_FILES_PER_RECORD:
            raise ArchiveError(f"{record} would exceed Zenodo's 100-file record limit")
    return tuple(plans)


def zip_info(item: FileEntry, compression_level: int) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo("Paper/" + item.path, ZIP_DATETIME)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.create_version = 45
    info.extract_version = 45
    info.external_attr = ((stat.S_IFREG | item.mode) & 0xFFFF) << 16
    info.file_size = item.size
    info._compresslevel = compression_level  # type: ignore[attr-defined]
    return info


def write_shard(
    plan: Shard,
    target: Path,
    *,
    compression_level: int,
    max_archive_bytes: int,
) -> dict[str, Any]:
    if target.exists() or target.is_symlink():
        raise ArchiveError(f"Refusing to replace archive: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        target,
        mode="x",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=compression_level,
        allowZip64=True,
        strict_timestamps=False,
    ) as archive:
        for item in plan.files:
            with item.source.open("rb") as source, archive.open(
                zip_info(item, compression_level),
                mode="w",
                force_zip64=True,
            ) as destination:
                for block in iter(lambda: source.read(CHUNK), b""):
                    destination.write(block)
    archive_size = target.stat().st_size
    if archive_size > max_archive_bytes or archive_size > ZENODO_INDIVIDUAL_FILE_LIMIT:
        raise ArchiveError(
            f"Compressed shard exceeds its hard ceiling: {target.name} ({archive_size})"
        )
    verify_shard(target, plan)
    return {
        "name": target.name,
        "relative_archive_path": f"{plan.record}/{target.name}",
        "member_count": len(plan.files),
        "uncompressed_bytes": plan.payload_bytes,
        "compressed_bytes": archive_size,
        "sha256": sha256_file(target),
    }


def verify_shard(path: Path, plan: Shard) -> None:
    expected = ["Paper/" + item.path for item in plan.files]
    with zipfile.ZipFile(path, mode="r") as archive:
        infos = archive.infolist()
        if [item.filename for item in infos] != expected or len(infos) != len(
            set(expected)
        ):
            raise ArchiveError(f"ZIP member inventory/order differs: {path}")
        for info, item in zip(infos, plan.files):
            if (
                info.date_time != ZIP_DATETIME
                or info.compress_type != zipfile.ZIP_DEFLATED
            ):
                raise ArchiveError(f"ZIP determinism metadata differs: {info.filename}")
            if info.file_size != item.size:
                raise ArchiveError(f"ZIP member size differs: {info.filename}")
            if ((info.external_attr >> 16) & 0o7777) != item.mode:
                raise ArchiveError(f"ZIP member POSIX mode differs: {info.filename}")
            digest = hashlib.sha256()
            count = 0
            with archive.open(info) as handle:
                for block in iter(lambda: handle.read(CHUNK), b""):
                    digest.update(block)
                    count += len(block)
            if count != item.size or digest.hexdigest() != item.sha256:
                raise ArchiveError(f"ZIP member content differs: {info.filename}")


def canonical_json(document: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def tree_digest(
    root_mode: str,
    directories: Sequence[Mapping[str, Any]],
    files: Sequence[Mapping[str, Any]],
) -> str:
    payload = {
        "root_mode": root_mode,
        "directories": list(directories),
        "files": list(files),
    }
    return hashlib.sha256(canonical_json(payload)).hexdigest()


def make_manifest(
    inventory: Inventory,
    plans: Sequence[Shard],
    archive_rows: Mapping[str, Mapping[str, Any]],
    *,
    contract_hash: str | None,
    payload_limit: int,
    max_archive_bytes: int,
    compression_level: int,
) -> dict[str, Any]:
    assignment = {item.path: plan.name for plan in plans for item in plan.files}
    directories = [
        {"path": item.path, "mode": mode_text(item.mode), "empty": item.empty}
        for item in inventory.directories
    ]
    files = [
        {
            "path": item.path,
            "mode": mode_text(item.mode),
            "size_bytes": item.size,
            "sha256": item.sha256,
            "record": item.record,
            "shard": assignment[item.path],
            "member": "Paper/" + item.path,
        }
        for item in inventory.files
    ]
    record_rows: dict[str, Any] = {}
    compressed_totals: dict[str, int] = {}
    for record in RECORDS:
        shards = [archive_rows[plan.name] for plan in plans if plan.record == record]
        compressed = sum(int(row["compressed_bytes"]) for row in shards)
        compressed_totals[record] = compressed
        record_rows[record] = {
            "one_logical_zenodo_record": True,
            "shard_count": len(shards),
            "zenodo_uploaded_shard_files": len(shards),
            "tree_file_count": sum(item.record == record for item in inventory.files),
            "uncompressed_bytes": sum(
                item.size for item in inventory.files if item.record == record
            ),
            "compressed_archive_bytes": compressed,
            "shards": shards,
        }
    if compressed_totals["software_core"] >= ZENODO_DEFAULT_RECORD_BYTES:
        raise ArchiveError(
            "software/core record must remain below 50,000,000,000 bytes"
        )
    if any(total > ZENODO_MAX_RECORD_BYTES for total in compressed_totals.values()):
        raise ArchiveError("A logical record exceeds the 200 GB quota-raised ceiling")
    additional_quota = sum(
        max(0, total - ZENODO_DEFAULT_RECORD_BYTES)
        for total in compressed_totals.values()
    )
    if additional_quota > ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE:
        raise ArchiveError(
            "The records exceed the account-level +150 GB quota allowance"
        )
    document: dict[str, Any] = {
        "schema": SCHEMA,
        "paper_root_name": "Paper",
        "archive_member_prefix": "Paper/",
        "compression": {"method": "ZIP_DEFLATED", "level": compression_level},
        "archive_byte_identity_scope": "deterministic within one Python/zlib runtime",
        "reconstruction_authority": "tree paths, sizes, SHA-256 values, POSIX modes, and directories",
        "canonical_exclusion_contract": {
            "script": "scripts/utilities/01_validate_bundle.py",
            "script_sha256": contract_hash,
            "local_only_directories": sorted(LOCAL_ONLY),
            "local_only_root_prefixes": list(LOCAL_ONLY_ROOT_PREFIXES),
            "local_only_root_files": sorted(LOCAL_ONLY_FILES),
            "local_only_file_suffixes": list(LOCAL_ONLY_FILE_SUFFIXES),
            "public_excluded_relative_prefixes": list(
                PUBLIC_EXCLUDED_RELATIVE_PREFIXES
            ),
            "redistribution_excluded_basenames": sorted(
                REDISTRIBUTION_EXCLUDED_BASENAMES
            ),
        },
        "limits": {
            "uncompressed_payload_per_shard": payload_limit,
            "compressed_bytes_per_shard": max_archive_bytes,
            "zenodo_individual_file_bytes": ZENODO_INDIVIDUAL_FILE_LIMIT,
            "zenodo_files_per_record": ZENODO_MAX_FILES_PER_RECORD,
            "zenodo_quota_raised_record_bytes": ZENODO_MAX_RECORD_BYTES,
            "zenodo_account_additional_allowance_bytes": ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE,
            "calculated_additional_quota_bytes": additional_quota,
        },
        "tree": {
            "root_mode": mode_text(inventory.root_mode),
            "directory_count": len(directories),
            "file_count": len(files),
            "payload_bytes": sum(item.size for item in inventory.files),
            "directories": directories,
            "files": files,
        },
        "records": record_rows,
    }
    document["tree"]["sha256"] = tree_digest(
        document["tree"]["root_mode"], directories, files
    )
    return document


def validate_manifest(document: Mapping[str, Any]) -> None:
    records = document.get("records", {})
    if document.get("schema") != SCHEMA or set(records) != set(RECORDS):
        raise ArchiveError("Manifest schema or logical-record inventory differs")
    tree = document.get("tree")
    limits = document.get("limits")
    if not isinstance(tree, dict) or not isinstance(limits, dict):
        raise ArchiveError("Manifest tree/limits are invalid")
    if not re.fullmatch(r"[0-7]{4}", str(tree.get("root_mode", ""))):
        raise ArchiveError("Manifest root mode is invalid")
    directories = tree.get("directories")
    files = tree.get("files")
    if not isinstance(directories, list) or not isinstance(files, list):
        raise ArchiveError("Manifest path inventories are invalid")
    payload_limit = limits.get("uncompressed_payload_per_shard")
    archive_limit = limits.get("compressed_bytes_per_shard")
    if (
        not isinstance(payload_limit, int)
        or not 0 < payload_limit <= DEFAULT_SHARD_PAYLOAD_BYTES
        or not isinstance(archive_limit, int)
        or not 0 < archive_limit <= DEFAULT_MAX_ARCHIVE_BYTES
    ):
        raise ArchiveError("Manifest shard ceilings are invalid")

    shard_record: dict[str, str] = {}
    shard_rows: dict[str, Mapping[str, Any]] = {}
    for record in RECORDS:
        row = records[record]
        shards = row.get("shards") if isinstance(row, dict) else None
        if not isinstance(shards, list) or len(shards) > ZENODO_MAX_FILES_PER_RECORD:
            raise ArchiveError(f"Invalid shard inventory for {record}")
        if (
            row.get("one_logical_zenodo_record") is not True
            or row.get("shard_count") != len(shards)
            or row.get("zenodo_uploaded_shard_files") != len(shards)
        ):
            raise ArchiveError(f"Logical-record/shard count differs for {record}")
        compressed_total = 0
        for shard in shards:
            if not isinstance(shard, dict):
                raise ArchiveError(f"Invalid shard row for {record}")
            name = safe_relative(shard["name"])
            if "/" in name or name in shard_record:
                raise ArchiveError(f"Duplicate/unsafe shard name: {name}")
            if shard["relative_archive_path"] != f"{record}/{name}":
                raise ArchiveError(f"Shard relative path differs: {name}")
            if not SHA_RE.fullmatch(shard["sha256"]):
                raise ArchiveError(f"Shard SHA-256 differs: {name}")
            if (
                not isinstance(shard.get("member_count"), int)
                or shard["member_count"] < 0
                or not isinstance(shard.get("uncompressed_bytes"), int)
                or shard["uncompressed_bytes"] < 0
                or not isinstance(shard.get("compressed_bytes"), int)
                or shard["compressed_bytes"] < 0
            ):
                raise ArchiveError(f"Invalid shard counts/sizes: {name}")
            if shard["uncompressed_bytes"] > payload_limit:
                raise ArchiveError(f"Shard payload exceeds plan: {name}")
            if shard["compressed_bytes"] > archive_limit:
                raise ArchiveError(f"Shard compressed size exceeds ceiling: {name}")
            compressed_total += shard["compressed_bytes"]
            shard_record[name] = record
            shard_rows[name] = shard
        if compressed_total != row.get("compressed_archive_bytes"):
            raise ArchiveError(f"Record compressed total differs: {record}")

    paths: set[str] = set()
    assignments: set[str] = set()
    per_shard_count = {name: 0 for name in shard_record}
    per_shard_payload = {name: 0 for name in shard_record}
    per_record_count = {record: 0 for record in RECORDS}
    per_record_payload = {record: 0 for record in RECORDS}
    for row in directories:
        path = safe_relative(row["path"])
        if (
            path in paths
            or not re.fullmatch(r"[0-7]{4}", str(row.get("mode", "")))
            or not isinstance(row.get("empty"), bool)
        ):
            raise ArchiveError(f"Duplicate/invalid directory row: {path}")
        if not is_public_relative_path(path, is_directory=True):
            raise ArchiveError(f"Manifest contains excluded directory: {path}")
        paths.add(path)
    for row in files:
        path = safe_relative(row["path"])
        record = row.get("record")
        shard = row.get("shard")
        if path in paths or path in assignments or record not in RECORDS:
            raise ArchiveError(f"Duplicate/invalid file row: {path}")
        if not is_public_relative_path(path, is_directory=False):
            raise ArchiveError(f"Manifest contains excluded file: {path}")
        if assign_record(path, row.get("size_bytes")) != record:
            raise ArchiveError(f"Logical-record classification differs: {path}")
        if shard_record.get(shard) != record:
            raise ArchiveError(f"Cross-record or missing shard assignment: {path}")
        if (
            not SHA_RE.fullmatch(str(row.get("sha256", "")))
            or not re.fullmatch(r"[0-7]{4}", str(row.get("mode", "")))
            or not isinstance(row.get("size_bytes"), int)
            or row["size_bytes"] < 0
            or row.get("member") != "Paper/" + path
        ):
            raise ArchiveError(f"Invalid content/mode/member row: {path}")
        assignments.add(path)
        per_shard_count[shard] += 1
        per_shard_payload[shard] += row["size_bytes"]
        per_record_count[record] += 1
        per_record_payload[record] += row["size_bytes"]
    for name, shard in shard_rows.items():
        if (
            shard.get("member_count") != per_shard_count[name]
            or shard.get("uncompressed_bytes") != per_shard_payload[name]
        ):
            raise ArchiveError(f"Shard manifest totals differ: {name}")
    for record in RECORDS:
        row = records[record]
        if (
            row.get("tree_file_count") != per_record_count[record]
            or row.get("uncompressed_bytes") != per_record_payload[record]
        ):
            raise ArchiveError(f"Logical-record tree totals differ: {record}")
    if tree.get("directory_count") != len(directories) or tree.get("file_count") != len(
        files
    ):
        raise ArchiveError("Manifest counts differ")
    if tree.get("payload_bytes") != sum(row["size_bytes"] for row in files):
        raise ArchiveError("Manifest payload total differs")
    if tree.get("sha256") != tree_digest(tree["root_mode"], directories, files):
        raise ArchiveError("Manifest tree digest differs")


def build_archives(
    paper_root: Path,
    output: Path,
    *,
    contract_script: Path | None = None,
    payload_limit: int = DEFAULT_SHARD_PAYLOAD_BYTES,
    max_archive_bytes: int = DEFAULT_MAX_ARCHIVE_BYTES,
    compression_level: int = 6,
    jobs: int = 1,
) -> Path:
    if not (0 < payload_limit <= DEFAULT_SHARD_PAYLOAD_BYTES):
        raise ArchiveError("Shard planning payload may not exceed 43,000,000,000 bytes")
    if not (0 < max_archive_bytes <= DEFAULT_MAX_ARCHIVE_BYTES):
        raise ArchiveError(
            "Compressed shard ceiling may not exceed 45,000,000,000 bytes"
        )
    if not 1 <= compression_level <= 9:
        raise ArchiveError("ZIP deflate compression level must be 1..9")
    if not isinstance(jobs, int) or isinstance(jobs, bool) or jobs < 1:
        raise ArchiveError("Parallel archive jobs must be a positive integer")
    lexical_paper = paper_root.expanduser().absolute()
    if lexical_paper.is_symlink():
        raise ArchiveError(f"Paper root may not be a symlink: {lexical_paper}")
    paper = lexical_paper.resolve()
    target = output.expanduser().absolute()
    resolved_target = target.resolve(strict=False)
    if (
        resolved_target == paper
        or resolved_target.is_relative_to(paper)
        or paper.is_relative_to(resolved_target)
    ):
        raise ArchiveError("Archive output must be outside and must not contain Paper")
    if target.exists() or target.is_symlink():
        raise ArchiveError(f"Archive output must be absent: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = target.with_name(f".{target.name}.building.{os.getpid()}")
    if stage.exists() or stage.is_symlink():
        raise ArchiveError(f"Temporary output already exists: {stage}")
    inventory, contract_hash = inventory_tree(paper, contract_script=contract_script)
    plans = plan_shards(inventory.files, payload_limit=payload_limit)
    archive_rows: dict[str, Mapping[str, Any]] = {}
    try:
        stage.mkdir(mode=0o700)
        for record in RECORDS:
            (stage / record).mkdir(mode=0o755)
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
            futures = {
                executor.submit(
                    write_shard,
                    plan,
                    stage / plan.record / plan.name,
                    compression_level=compression_level,
                    max_archive_bytes=max_archive_bytes,
                ): plan
                for plan in plans
            }
            for future in concurrent.futures.as_completed(futures):
                plan = futures[future]
                archive_rows[plan.name] = future.result()
        manifest = make_manifest(
            inventory,
            plans,
            archive_rows,
            contract_hash=contract_hash,
            payload_limit=payload_limit,
            max_archive_bytes=max_archive_bytes,
            compression_level=compression_level,
        )
        validate_manifest(manifest)
        manifest_path = stage / "manifest.json"
        manifest_path.write_bytes(canonical_json(manifest))
        os.chmod(manifest_path, 0o644)
        os.chmod(stage, 0o755)
        os.replace(stage, target)
    except Exception:
        if stage.exists() and not stage.is_symlink():
            shutil.rmtree(stage)
        raise
    return target


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--contract-script",
        type=Path,
        help="Defaults to <Paper>/scripts/utilities/01_validate_bundle.py.",
    )
    parser.add_argument(
        "--compression-level", type=int, default=6, choices=range(1, 10)
    )
    parser.add_argument("--jobs", type=int, default=1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paper = args.paper_root.expanduser().absolute()
    contract = (
        args.contract_script.expanduser().resolve()
        if args.contract_script
        else paper / "scripts/utilities/01_validate_bundle.py"
    )
    destination = build_archives(
        paper,
        args.output,
        contract_script=contract,
        compression_level=args.compression_level,
        jobs=args.jobs,
    )
    manifest = json.loads((destination / "manifest.json").read_text(encoding="utf-8"))
    print(f"[PASS] release archive set: {destination}")
    for record in RECORDS:
        row = manifest["records"][record]
        print(
            f"[RECORD] {record}: {row['shard_count']} shards, "
            f"{row['uncompressed_bytes']} raw bytes, "
            f"{row['compressed_archive_bytes']} ZIP bytes"
        )
    print(
        "[QUOTA] additional bytes: "
        f"{manifest['limits']['calculated_additional_quota_bytes']} / "
        f"{ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ArchiveError as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        raise SystemExit(2)
