#!/usr/bin/env python3
"""Production-only, fail-closed executor for the four-record Zenodo release.

This module is deliberately separate from :mod:`release_control`, whose remote
actions remain disabled.  Nothing is sent to Zenodo unless the corresponding
command receives both its explicit execution flag and ``--production``.
Publishing additionally requires a literal irreversible-action confirmation,
the current state-file SHA-256, verified downloads, and a successful all-figure
reproduction launched from the remotely reconstructed tree.
"""

from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
import queue
import re
import stat
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Callable, Iterable, Mapping, MutableMapping, Sequence


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(HERE))

import build_archives as archives  # noqa: E402
import reconstruct_paper as reconstruction  # noqa: E402
import release_control as control  # noqa: E402


STATE_SCHEMA = "apotome-zenodo-production-state-v1"
REQUEST_SCHEMA = control.REQUEST_SCHEMA
PRODUCTION_API = "https://zenodo.org/api"
PRODUCTION_HOST = "zenodo.org"
TOKEN_ENV = "ZENODO_TOKEN"
PUBLISH_CONFIRMATION = "PUBLISH_ALL_FOUR_PUBLIC_RECORDS"
RECORDS = tuple(archives.RECORDS)
EXPECTED_LICENSES = {
    "software_core": "other-open",
    "main_figure_data": "cc-by-4.0",
    "supplementary_figure_data": "cc-by-4.0",
    "figs3_wsi": "cc-by-4.0",
}
EXPECTED_RIGHTS_CONFIRMATIONS = (
    "all_authors_approved_public_release",
    "software_manuscript_and_figure_redistribution_rights_confirmed",
    "main_figure_data_redistribution_rights_confirmed",
    "supplementary_figure_data_redistribution_rights_confirmed",
    "histology_wsi_redistribution_rights_confirmed",
    "third_party_materials_and_model_outputs_reviewed",
    "no_sensitive_or_restricted_data",
)
MUTATING_COMMAND_FLAGS = {
    "drafts": "--execute-drafts",
    "upload": "--execute-uploads",
    "publish": "--execute-publish",
}
SHA256_RE = re.compile(r"[0-9a-f]{64}")
MD5_RE = re.compile(r"[0-9a-f]{32}")
DOI_RE = re.compile(r"10\.[0-9]{4,9}/[-._;()/:A-Za-z0-9]+")
RELEASE_ID_RE = re.compile(r"apotome-[0-9a-f]{24}")
RETRYABLE_HTTP = {408, 409, 425, 429, 500, 502, 503, 504}
CHUNK = 8 * 1024 * 1024
MAX_JSON_BYTES = 32 * 1024 * 1024
DEFAULT_DOWNLOAD_WORKERS = 4
MAX_DOWNLOAD_WORKERS = 8
SEGMENTED_TRANSPORT_SCHEMA_V1 = "allulose-zenodo-segmented-transport-v1"
SEGMENTED_TRANSPORT_SCHEMA = "allulose-zenodo-segmented-transport-v2"
SEGMENTED_TRANSPORT_SCHEMAS = {
    SEGMENTED_TRANSPORT_SCHEMA_V1,
    SEGMENTED_TRANSPORT_SCHEMA,
}
TRANSPORT_MANIFEST_FILENAME = "transport_manifest.json"
TRANSPORT_REASSEMBLY_INSTRUCTION = (
    "Concatenate each archive's segments in ascending number order; "
    "the resulting ZIP must match the recorded size and SHA-256."
)


class ZenodoReleaseError(RuntimeError):
    """The production release executor failed closed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def canonical_json(document: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def hashes_file(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5(usedforsecurity=False)
    count = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            count += len(block)
            sha256.update(block)
            md5.update(block)
    return sha256.hexdigest(), md5.hexdigest(), count


def _unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ZenodoReleaseError(f"Duplicate JSON key: {key!r}")
        result[key] = value
    return result


def read_json(path: Path, label: str) -> dict[str, Any]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise ZenodoReleaseError(f"Missing or linked {label}: {lexical}")
    if lexical.stat().st_size > MAX_JSON_BYTES:
        raise ZenodoReleaseError(f"Unexpectedly large {label}: {lexical}")
    try:
        value = json.loads(
            lexical.read_text(encoding="utf-8"), object_pairs_hook=_unique_object
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ZenodoReleaseError(f"Invalid {label}: {lexical}: {exc}") from exc
    if not isinstance(value, dict):
        raise ZenodoReleaseError(f"{label} must contain one JSON object")
    return value


def state_sha256(state: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_json(state))


def _assert_no_secret_keys(value: Any, path: str = "state") -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            if re.search(r"(?i)(token|authorization|password|secret)", str(key)):
                raise ZenodoReleaseError(
                    f"Credential-like key is forbidden in persisted state: {path}.{key}"
                )
            _assert_no_secret_keys(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _assert_no_secret_keys(child, f"{path}[{index}]")


def write_state(path: Path, state: Mapping[str, Any]) -> None:
    _assert_no_secret_keys(state)
    target = path.expanduser().absolute()
    if target.is_symlink():
        raise ZenodoReleaseError(f"State path may not be a symlink: {target}")
    if not target.parent.is_dir() or target.parent.is_symlink():
        raise ZenodoReleaseError(
            f"State parent must be a regular directory: {target.parent}"
        )
    payload = canonical_json(state)
    temporary = target.with_name(f".{target.name}.tmp.{os.getpid()}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(temporary, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, target)
        os.chmod(target, 0o600)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def load_state(path: Path) -> dict[str, Any]:
    target = path.expanduser().absolute()
    state = read_json(target, "Zenodo release state")
    mode = stat.S_IMODE(target.stat().st_mode)
    if mode & 0o077:
        raise ZenodoReleaseError(
            f"State file must not be group/world accessible (expected mode 0600): {target}"
        )
    validate_state(state)
    return state


def _plain_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _creator_names(value: Any) -> list[str] | None:
    if not isinstance(value, list) or any(
        not isinstance(row, Mapping) for row in value
    ):
        return None
    names = [row.get("name") for row in value]
    if any(not control.valid_zenodo_name(name) for name in names):
        return None
    return names


def _segmented_binding(state: Mapping[str, Any]) -> Mapping[str, Any] | None:
    binding = state.get("segmented_transport")
    if binding is None:
        return None
    expected_keys = {"schema", "manifest_path", "size_bytes", "sha256", "md5"}
    if not isinstance(binding, Mapping) or set(binding) != expected_keys:
        raise ZenodoReleaseError("Segmented-transport state binding is invalid")
    path = binding.get("manifest_path")
    if (
        binding.get("schema") not in SEGMENTED_TRANSPORT_SCHEMAS
        or not isinstance(path, str)
        or not Path(path).is_absolute()
        or not _plain_int(binding.get("size_bytes"))
        or binding["size_bytes"] <= 0
        or not SHA256_RE.fullmatch(str(binding.get("sha256", "")))
        or not MD5_RE.fullmatch(str(binding.get("md5", "")))
    ):
        raise ZenodoReleaseError("Segmented-transport state binding is invalid")
    return binding


def _canonical_files_by_name(
    state: Mapping[str, Any], record: str
) -> dict[str, Mapping[str, Any]]:
    return {
        item["name"]: item
        for item in state["records"][record]["files"]
        if item["name"] != "manifest.json"
    }


def load_segmented_transport(
    state: Mapping[str, Any],
) -> dict[str, Any] | None:
    """Load and strictly validate the optional segmented-transport manifest."""
    binding = _segmented_binding(state)
    if binding is None:
        return None
    path = Path(binding["manifest_path"])
    document = read_json(path, "segmented-transport manifest")
    sha256, md5, size = hashes_file(path)
    if (
        size != binding["size_bytes"]
        or sha256 != binding["sha256"]
        or md5 != binding["md5"]
    ):
        raise ZenodoReleaseError("Segmented-transport manifest file binding differs")
    expected_root_keys = {
        "schema",
        "release_id",
        "canonical_manifest",
        "records",
        "reassembly",
    }
    schema = document.get("schema")
    records = document.get("records")
    if (
        set(document) != expected_root_keys
        or schema not in SEGMENTED_TRANSPORT_SCHEMAS
        or schema != binding["schema"]
        or document.get("release_id") != state.get("release_id")
        or document.get("canonical_manifest") != state.get("manifest")
        or document.get("reassembly") != TRANSPORT_REASSEMBLY_INSTRUCTION
        or not isinstance(records, Mapping)
        or set(records) != set(RECORDS)
    ):
        raise ZenodoReleaseError("Segmented-transport manifest binding differs")

    physical_names: dict[str, set[str]] = {record: set() for record in RECORDS}
    physical_bytes: dict[str, int] = {record: 0 for record in RECORDS}
    for record in RECORDS:
        record_document = records[record]
        if (
            not isinstance(record_document, Mapping)
            or set(record_document) != {"archives"}
            or not isinstance(record_document["archives"], list)
        ):
            raise ZenodoReleaseError(f"Segmented-transport record is invalid: {record}")
        manifest_rows = [
            item
            for item in state["records"][record]["files"]
            if item["name"] == "manifest.json"
        ]
        if record == "software_core":
            if (
                len(manifest_rows) != 1
                or manifest_rows[0]["relative_local_path"] != "manifest.json"
                or manifest_rows[0]["size_bytes"] != state["manifest"]["size_bytes"]
                or manifest_rows[0]["sha256"] != state["manifest"]["sha256"]
            ):
                raise ZenodoReleaseError(
                    "Canonical manifest inventory differs in segmented transport"
                )
        elif manifest_rows:
            raise ZenodoReleaseError(
                f"Canonical manifest belongs only to software_core, not {record}"
            )
        canonical = _canonical_files_by_name(state, record)
        observed_archives = record_document["archives"]
        if [
            row.get("name") for row in observed_archives if isinstance(row, Mapping)
        ] != list(canonical) or len(observed_archives) != len(canonical):
            raise ZenodoReleaseError(
                f"Segmented-transport archive inventory differs: {record}"
            )
        segment_names: set[str] = set()
        for archive in observed_archives:
            expected_archive_keys = {
                "name",
                "relative_local_path",
                "size_bytes",
                "sha256",
                "md5",
                "segment_bytes",
                "segments",
            }
            if (
                not isinstance(archive, Mapping)
                or set(archive) != expected_archive_keys
            ):
                raise ZenodoReleaseError(
                    f"Segmented-transport archive entry is invalid: {record}"
                )
            expected = canonical.get(archive.get("name"))
            segments = archive.get("segments")
            segment_bytes = archive.get("segment_bytes")
            if (
                expected is None
                or archive.get("relative_local_path")
                != expected.get("relative_local_path")
                or archive.get("size_bytes") != expected.get("size_bytes")
                or archive.get("sha256") != expected.get("sha256")
                or not MD5_RE.fullmatch(str(archive.get("md5", "")))
                or (
                    expected.get("md5") is not None
                    and archive.get("md5") != expected.get("md5")
                )
                or not _plain_int(segment_bytes)
                or segment_bytes <= 0
                or not isinstance(segments, list)
                or not segments
                or len(segments) > 999
            ):
                raise ZenodoReleaseError(
                    f"Segmented-transport archive binding differs: "
                    f"{record}/{archive.get('name')}"
                )
            offset = 0
            count = len(segments)
            for number, segment in enumerate(segments, start=1):
                expected_segment_keys = {
                    "name",
                    "number",
                    "offset_bytes",
                    "size_bytes",
                    "sha256",
                    "md5",
                }
                observed_segment_keys = (
                    set(segment) if isinstance(segment, Mapping) else set()
                )
                storage_record = (
                    segment.get("storage_record", record)
                    if isinstance(segment, Mapping)
                    else None
                )
                expected_name = f"{archive['name']}.chunk{number:03d}-of-{count:03d}"
                if (
                    not isinstance(segment, Mapping)
                    or (
                        observed_segment_keys != expected_segment_keys
                        and not (
                            schema == SEGMENTED_TRANSPORT_SCHEMA
                            and observed_segment_keys
                            == expected_segment_keys | {"storage_record"}
                        )
                    )
                    or storage_record not in RECORDS
                    or (
                        "storage_record" in observed_segment_keys
                        and storage_record == record
                    )
                    or segment.get("name") != expected_name
                    or PurePosixPath(expected_name).name != expected_name
                    or expected_name in segment_names
                    or not _plain_int(segment.get("number"))
                    or segment.get("number") != number
                    or not _plain_int(segment.get("offset_bytes"))
                    or segment.get("offset_bytes") != offset
                    or not _plain_int(segment.get("size_bytes"))
                    or segment["size_bytes"] <= 0
                    or segment["size_bytes"] > segment_bytes
                    or (number < count and segment["size_bytes"] != segment_bytes)
                    or not SHA256_RE.fullmatch(str(segment.get("sha256", "")))
                    or not MD5_RE.fullmatch(str(segment.get("md5", "")))
                ):
                    raise ZenodoReleaseError(
                        f"Segmented-transport segment differs: "
                        f"{record}/{archive['name']} #{number}"
                    )
                segment_names.add(expected_name)
                if expected_name in physical_names[storage_record]:
                    raise ZenodoReleaseError(
                        f"Segmented-transport physical filename collides: "
                        f"{storage_record}/{expected_name}"
                    )
                physical_names[storage_record].add(expected_name)
                physical_bytes[storage_record] += segment["size_bytes"]
                offset += segment["size_bytes"]
            if offset != archive["size_bytes"]:
                raise ZenodoReleaseError(
                    f"Segmented-transport segment sizes differ: "
                    f"{record}/{archive['name']}"
                )
    if physical_names["software_core"] & {
        "manifest.json",
        TRANSPORT_MANIFEST_FILENAME,
    }:
        raise ZenodoReleaseError(
            "Segmented transport collides with an ancillary filename"
        )
    physical_names["software_core"].update(
        {"manifest.json", TRANSPORT_MANIFEST_FILENAME}
    )
    physical_bytes["software_core"] += (
        state["manifest"]["size_bytes"] + binding["size_bytes"]
    )
    for record in RECORDS:
        if len(physical_names[record]) > archives.ZENODO_MAX_FILES_PER_RECORD:
            raise ZenodoReleaseError(
                f"Zenodo file-count limit exceeded by segmented transport: {record}"
            )
        allocated = state["records"][record]["quota"]["allocated_record_bytes"]
        physical_limit = (
            min(archives.ZENODO_MAX_RECORD_BYTES, allocated)
            if _plain_int(allocated)
            else archives.ZENODO_MAX_RECORD_BYTES
        )
        if physical_bytes[record] > physical_limit:
            raise ZenodoReleaseError(
                f"Zenodo record-byte limit exceeded by segmented transport: {record}"
            )
    return document


def validate_segmented_transport(state: Mapping[str, Any]) -> None:
    """Validate an optional segmented transport binding without changing state."""
    load_segmented_transport(state)


def bind_segmented_transport(
    state: MutableMapping[str, Any], manifest_path: Path
) -> dict[str, Any]:
    """Bind a sealed transport manifest and invalidate stale remote gates."""
    existing = _segmented_binding(state)
    target = manifest_path.expanduser().absolute()
    if target.is_symlink() or not target.is_file():
        raise ZenodoReleaseError(
            f"Segmented-transport manifest is missing or linked: {target}"
        )
    sha256, md5, size = hashes_file(target)
    document = read_json(target, "segmented-transport manifest")
    schema = document.get("schema")
    if schema not in SEGMENTED_TRANSPORT_SCHEMAS:
        raise ZenodoReleaseError("Segmented-transport manifest schema differs")
    binding = {
        "schema": schema,
        "manifest_path": str(target),
        "size_bytes": size,
        "sha256": sha256,
        "md5": md5,
    }
    if existing is not None:
        if dict(existing) != binding:
            raise ZenodoReleaseError(
                "Segmented transport is already bound to a different manifest"
            )
        validate_state(state)
        return binding
    validate_state(state)
    if any(state["records"][record].get("published") for record in RECORDS):
        raise ZenodoReleaseError("Cannot bind segmented transport after publication")
    state["segmented_transport"] = binding
    load_segmented_transport(state)
    for record in RECORDS:
        row = state["records"][record]
        row["remote_verified"] = False
        for item in row["files"]:
            item["uploaded"] = False
            item["remote_verified"] = False
    state["remote_download"] = {
        "status": "NOT_STARTED",
        "root": None,
        "completed_files": [],
        "completed_segments": [],
        "reconstructed_root": None,
    }
    state["reproduction"] = {
        "status": "NOT_STARTED",
        "log_path": None,
        "log_sha256": None,
        "command": None,
    }
    _touch(state)
    validate_state(state)
    return binding


def expected_remote_items(
    state: Mapping[str, Any], record: str
) -> list[dict[str, Any]]:
    """Return the exact Zenodo inventory for canonical or segmented transport."""
    if record not in RECORDS:
        raise ZenodoReleaseError(f"Unknown logical record: {record}")
    transport = load_segmented_transport(state)
    if transport is None:
        return [dict(item) for item in state["records"][record]["files"]]
    result: list[dict[str, Any]] = []
    for logical_record in RECORDS:
        for archive in transport["records"][logical_record]["archives"]:
            for segment in archive["segments"]:
                storage_record = segment.get("storage_record", logical_record)
                if storage_record != record:
                    continue
                result.append(
                    {
                        "name": segment["name"],
                        "size_bytes": segment["size_bytes"],
                        "sha256": segment["sha256"],
                        "md5": segment["md5"],
                        "transport_kind": "segment",
                        "logical_record": logical_record,
                        "storage_record": storage_record,
                        "canonical_name": archive["name"],
                        "canonical_relative_local_path": archive["relative_local_path"],
                        "canonical_size_bytes": archive["size_bytes"],
                        "canonical_sha256": archive["sha256"],
                        "canonical_md5": archive["md5"],
                        "number": segment["number"],
                        "offset_bytes": segment["offset_bytes"],
                    }
                )
    if record == "software_core":
        manifest_item = next(
            item
            for item in state["records"][record]["files"]
            if item["name"] == "manifest.json"
        )
        path = _local_path(state, manifest_item)
        sha256, md5, size = hashes_file(path)
        if (
            size != manifest_item["size_bytes"]
            or sha256 != manifest_item["sha256"]
            or (manifest_item.get("md5") is not None and manifest_item["md5"] != md5)
        ):
            raise ZenodoReleaseError("Canonical manifest local binding differs")
        result.append(
            {
                **dict(manifest_item),
                "md5": md5,
                "transport_kind": "canonical_manifest",
            }
        )
        binding = _segmented_binding(state)
        if binding is None:
            raise ZenodoReleaseError("Segmented-transport state binding disappeared")
        result.append(
            {
                "name": TRANSPORT_MANIFEST_FILENAME,
                "relative_local_path": TRANSPORT_MANIFEST_FILENAME,
                "size_bytes": binding["size_bytes"],
                "sha256": binding["sha256"],
                "md5": binding["md5"],
                "transport_kind": "transport_manifest",
            }
        )
    return result


def _validated_quota_state(record: str, row: Mapping[str, Any]) -> tuple[int, bool]:
    """Return quota requirements derived from the immutable file inventory.

    Pre-draft approvals are deliberately not quota authority. The persisted
    quota fields must remain an exact, redundant representation of the file
    inventory and one of the two valid post-initialization phases below.
    """
    files = row.get("files")
    if (
        not isinstance(files, list)
        or not files
        or any(
            not isinstance(item, Mapping)
            or not _plain_int(item.get("size_bytes"))
            or item["size_bytes"] < 0
            for item in files
        )
    ):
        raise ZenodoReleaseError(
            f"Cannot derive Zenodo quota from the file inventory for {record}"
        )
    required = sum(item["size_bytes"] for item in files)
    if required > archives.ZENODO_MAX_RECORD_BYTES:
        raise ZenodoReleaseError(
            f"Zenodo record {record} requires {required:,} bytes, exceeding the "
            f"{archives.ZENODO_MAX_RECORD_BYTES:,}-byte record limit"
        )
    confirmation_required = required > archives.ZENODO_DEFAULT_RECORD_BYTES
    quota = row.get("quota")
    expected_keys = {
        "required_record_bytes",
        "confirmation_required",
        "confirmed_for_deposition_id",
        "allocated_record_bytes",
    }
    if not isinstance(quota, Mapping) or set(quota) != expected_keys:
        raise ZenodoReleaseError(f"Zenodo quota state is invalid for {record}")
    if (
        not _plain_int(quota.get("required_record_bytes"))
        or quota["required_record_bytes"] != required
    ):
        raise ZenodoReleaseError(
            f"Zenodo quota required bytes differ from the file inventory for {record}"
        )
    if quota.get("confirmation_required") is not confirmation_required:
        raise ZenodoReleaseError(
            f"Zenodo quota threshold state differs from the file inventory for {record}"
        )

    confirmed_for = quota.get("confirmed_for_deposition_id")
    allocated = quota.get("allocated_record_bytes")
    if not confirmation_required:
        if (
            confirmed_for is not None
            or allocated != archives.ZENODO_DEFAULT_RECORD_BYTES
        ):
            raise ZenodoReleaseError(
                f"Small Zenodo record quota state is inconsistent for {record}"
            )
        return required, False

    if confirmed_for is None and allocated is None:
        return required, True
    deposition_id = row.get("deposition_id")
    if (
        not _plain_int(confirmed_for)
        or confirmed_for <= 0
        or not _plain_int(deposition_id)
        or deposition_id <= 0
        or confirmed_for != deposition_id
        or not _plain_int(allocated)
        or allocated < required
        or allocated > archives.ZENODO_MAX_RECORD_BYTES
    ):
        raise ZenodoReleaseError(
            f"Large Zenodo record quota is not bound to its current draft for {record}"
        )
    return required, True


def validate_state(state: Mapping[str, Any]) -> None:
    _assert_no_secret_keys(state)
    if state.get("schema") != STATE_SCHEMA:
        raise ZenodoReleaseError("Zenodo state schema differs")
    if (
        state.get("api_base") != PRODUCTION_API
        or state.get("environment") != "production"
    ):
        raise ZenodoReleaseError("Only the production Zenodo endpoint is supported")
    if not RELEASE_ID_RE.fullmatch(str(state.get("release_id", ""))):
        raise ZenodoReleaseError("Release execution identifier is invalid")
    if not SHA256_RE.fullmatch(str(state.get("request_sha256", ""))):
        raise ZenodoReleaseError("Request SHA-256 binding is invalid")
    gates = state.get("reviewed_gates")
    rights = gates.get("rights_confirmations") if isinstance(gates, Mapping) else None
    if (
        not isinstance(gates, Mapping)
        or gates.get("request_operation") != "zenodo-production-draft"
        or gates.get("public_visibility_acknowledged") is not True
        or gates.get("local_clean_room_reproduction_passed") is not True
        or not isinstance(rights, Mapping)
        or set(rights) != set(EXPECTED_RIGHTS_CONFIRMATIONS)
        or any(rights[key] is not True for key in EXPECTED_RIGHTS_CONFIRMATIONS)
    ):
        raise ZenodoReleaseError(
            "Persisted rights/public/local-clean-room gates are incomplete"
        )
    github_url = state.get("github_url")
    if not isinstance(github_url, str) or not re.fullmatch(
        r"https://github\.com/[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})/" r"[A-Za-z0-9._-]+",
        github_url,
    ):
        raise ZenodoReleaseError("State GitHub URL is invalid")
    if state.get("related_identifier_policy") != {
        "mutual_reserved_dois": "isSupplementTo",
        "github_repository": "isSupplementTo",
        "article_doi": "isSupplementTo",
    }:
        raise ZenodoReleaseError("Related-identifier policy differs")
    manifest = state.get("manifest")
    if (
        not isinstance(manifest, Mapping)
        or not SHA256_RE.fullmatch(str(manifest.get("sha256", "")))
        or not SHA256_RE.fullmatch(str(manifest.get("tree_sha256", "")))
        or not _plain_int(manifest.get("size_bytes"))
    ):
        raise ZenodoReleaseError("Manifest state binding is invalid")
    records = state.get("records")
    if not isinstance(records, Mapping) or set(records) != set(RECORDS):
        raise ZenodoReleaseError("State must contain exactly the four release records")
    relative_paths: set[str] = set()
    for name in RECORDS:
        row = records[name]
        if (
            not isinstance(row, Mapping)
            or row.get("license") != EXPECTED_LICENSES[name]
        ):
            raise ZenodoReleaseError(f"License binding differs for {name}")
        if _creator_names(row.get("creators")) != list(control.EXPECTED_ZENODO_NAMES):
            raise ZenodoReleaseError(
                f"Creator names/order differ from the confirmed list for {name}"
            )
        files = row.get("files")
        if not isinstance(files, list) or not files:
            raise ZenodoReleaseError(f"Expected file inventory is absent for {name}")
        names: set[str] = set()
        for file_row in files:
            if not isinstance(file_row, Mapping):
                raise ZenodoReleaseError(f"Invalid expected-file entry for {name}")
            filename = file_row.get("name")
            if (
                not isinstance(filename, str)
                or filename in names
                or PurePosixPath(filename).name != filename
                or not _plain_int(file_row.get("size_bytes"))
                or file_row["size_bytes"] < 0
                or not SHA256_RE.fullmatch(str(file_row.get("sha256", "")))
            ):
                raise ZenodoReleaseError(
                    f"Invalid expected file for {name}: {filename!r}"
                )
            md5 = file_row.get("md5")
            if md5 is not None and not MD5_RE.fullmatch(str(md5)):
                raise ZenodoReleaseError(f"Invalid local MD5 state for {filename}")
            relative = file_row.get("relative_local_path")
            try:
                archives.safe_relative(relative)
            except (TypeError, archives.ArchiveError) as exc:
                raise ZenodoReleaseError(
                    f"Invalid local archive path for {filename}: {relative!r}"
                ) from exc
            if relative in relative_paths:
                raise ZenodoReleaseError(
                    f"Duplicate local archive path across records: {relative}"
                )
            relative_paths.add(relative)
            names.add(filename)
        _validated_quota_state(name, row)
    validate_segmented_transport(state)


def require_mutation_gate(command: str, *, execute: bool, production: bool) -> None:
    flag = MUTATING_COMMAND_FLAGS[command]
    if not execute or not production:
        raise ZenodoReleaseError(
            f"Remote mutation refused: {command} requires both {flag} and --production"
        )


def _record_resource_type(upload_type: str) -> str:
    if upload_type not in {"software", "dataset"}:
        raise ZenodoReleaseError(f"Unsupported Zenodo upload type: {upload_type!r}")
    return upload_type


def _normal_doi(value: Any, *, optional: bool = False) -> str | None:
    if value is None and optional:
        return None
    if not isinstance(value, str):
        raise ZenodoReleaseError(f"Invalid DOI: {value!r}")
    candidate = value.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if candidate.lower().startswith(prefix):
            candidate = candidate[len(prefix) :]
            break
    if not DOI_RE.fullmatch(candidate):
        raise ZenodoReleaseError(f"Invalid DOI: {value!r}")
    return candidate


def _remote_checksum(value: Any) -> str:
    checksum = str(value or "").lower()
    if checksum.startswith("md5:"):
        checksum = checksum[4:]
    if not MD5_RE.fullmatch(checksum):
        raise ZenodoReleaseError(f"Zenodo returned an invalid MD5 checksum: {value!r}")
    return checksum


def normal_remote_files(deposition: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    files = deposition.get("files", [])
    if not isinstance(files, list):
        raise ZenodoReleaseError("Zenodo deposition files are not a list")
    result: dict[str, dict[str, Any]] = {}
    for raw in files:
        if not isinstance(raw, Mapping):
            raise ZenodoReleaseError("Zenodo returned an invalid file resource")
        name = raw.get("filename", raw.get("name", raw.get("key")))
        size = raw.get("filesize", raw.get("size"))
        if not isinstance(name, str) or PurePosixPath(name).name != name:
            raise ZenodoReleaseError(f"Zenodo returned an unsafe filename: {name!r}")
        try:
            size_int = int(size)
        except (TypeError, ValueError) as exc:
            raise ZenodoReleaseError(
                f"Zenodo returned an invalid size for {name}"
            ) from exc
        if size_int < 0 or name in result:
            raise ZenodoReleaseError(
                f"Zenodo returned a duplicate/invalid file: {name}"
            )
        links = raw.get("links", {})
        result[name] = {
            "name": name,
            "size_bytes": size_int,
            "md5": _remote_checksum(raw.get("checksum")),
            "links": dict(links) if isinstance(links, Mapping) else {},
        }
    return result


def _validate_request(request: Mapping[str, Any], plan: Mapping[str, Any]) -> None:
    if (
        request.get("schema") != REQUEST_SCHEMA
        or request.get("operation") != "zenodo-production-draft"
        or request.get("status") != "READY_FOR_SEPARATELY_REVIEWED_EXECUTOR"
        or request.get("local_only") is not True
        or request.get("remote_execution_enabled") is not False
        or request.get("credentials_embedded") is not False
    ):
        raise ZenodoReleaseError(
            "Release request is not an approved production-draft request"
        )
    visibility = request.get("visibility", {})
    if visibility != control.REQUIRED_VISIBILITY:
        raise ZenodoReleaseError("Release request visibility is not public/open")
    approvals = request.get("approvals", {})
    if (
        not isinstance(approvals, Mapping)
        or approvals.get("public_visibility_acknowledged") is not True
        or approvals.get("clean_room_reproduction_passed") is not True
    ):
        raise ZenodoReleaseError(
            "Public visibility and local clean-room gates are not approved"
        )
    rights = approvals.get("rights_confirmations")
    if not isinstance(rights, Mapping) or any(
        rights.get(key) is not True
        for key in plan["required_inputs"]["rights_confirmations"]
    ):
        raise ZenodoReleaseError("Every rights confirmation must be explicitly true")
    records = request.get("zenodo_records")
    if not isinstance(records, Mapping) or set(records) != set(RECORDS):
        raise ZenodoReleaseError(
            "Release request must contain exactly four Zenodo records"
        )
    try:
        expected_authors = [row["zenodo_name"] for row in plan["manuscript"]["authors"]]
    except (KeyError, TypeError) as exc:
        raise ZenodoReleaseError(
            "Reviewed plan is missing structured Zenodo creator names"
        ) from exc
    if expected_authors != list(control.EXPECTED_ZENODO_NAMES) or any(
        not control.valid_zenodo_name(name) for name in expected_authors
    ):
        raise ZenodoReleaseError(
            "Reviewed-plan Zenodo creator names/order differ from the confirmed list"
        )
    for name in RECORDS:
        row = records[name]
        metadata = row.get("metadata", {}) if isinstance(row, Mapping) else {}
        creators = metadata.get("creators", []) if isinstance(metadata, Mapping) else []
        if (
            metadata.get("title") != plan["records"][name]["title"]
            or metadata.get("description") != plan["records"][name]["description"]
            or metadata.get("resource_type") != plan["records"][name]["resource_type"]
            or metadata.get("license") != EXPECTED_LICENSES[name]
            or row.get("access") != control.REQUIRED_VISIBILITY["zenodo"]
            or _creator_names(creators) != expected_authors
        ):
            raise ZenodoReleaseError(
                f"Request metadata differs from the reviewed plan: {name}"
            )
    quotas = approvals.get("quota_confirmations")
    expected_quota_keys = set(plan["required_inputs"]["quota_confirmations"])
    if (
        not isinstance(quotas, Mapping)
        or set(quotas) != expected_quota_keys
        or any(not isinstance(quotas[key], bool) for key in expected_quota_keys)
    ):
        raise ZenodoReleaseError(
            "Pre-draft quota status must cover exactly all three data records"
        )


def _expected_files(
    request: Mapping[str, Any], archive_root: Path, summary: Mapping[str, Any]
) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    request_records = request["zenodo_records"]
    for record in RECORDS:
        rows: list[dict[str, Any]] = []
        request_files = list(request_records[record]["archive_shards"]) + list(
            request_records[record]["ancillary_files"]
        )
        expected_summary = summary["records"][record]["shards"]
        if request_records[record]["archive_shards"] != expected_summary:
            raise ZenodoReleaseError(f"Request shard inventory differs: {record}")
        for item in request_files:
            name = item.get("name")
            if name == "manifest.json":
                relative = "manifest.json"
            else:
                relative = f"{record}/{name}"
            path = archive_root.joinpath(*PurePosixPath(relative).parts)
            if path.is_symlink() or not path.is_file():
                raise ZenodoReleaseError(f"Missing or linked upload source: {path}")
            if path.stat().st_size != item.get("size_bytes"):
                raise ZenodoReleaseError(f"Upload source size differs: {path}")
            if sha256_file(path) != item.get("sha256"):
                raise ZenodoReleaseError(f"Upload source SHA-256 differs: {path}")
            rows.append(
                {
                    "name": name,
                    "relative_local_path": relative,
                    "size_bytes": item["size_bytes"],
                    "sha256": item["sha256"],
                    "md5": None,
                    "uploaded": False,
                    "remote_verified": False,
                }
            )
        if len(rows) > archives.ZENODO_MAX_FILES_PER_RECORD:
            raise ZenodoReleaseError(f"Zenodo file-count limit exceeded: {record}")
        total = sum(item["size_bytes"] for item in rows)
        if total > archives.ZENODO_MAX_RECORD_BYTES:
            raise ZenodoReleaseError(f"Zenodo 200 GB record limit exceeded: {record}")
        result[record] = rows
    return result


def _validated_adopted_drafts(
    adopted: Mapping[str, int] | None,
) -> dict[str, int]:
    if adopted is None:
        return {}
    if not isinstance(adopted, Mapping):
        raise ZenodoReleaseError("Adopted drafts must be a record-to-deposition map")
    result: dict[str, int] = {}
    for record, deposition_id in adopted.items():
        if record not in RECORDS:
            raise ZenodoReleaseError(f"Unknown adopted-draft record: {record}")
        if (
            not _plain_int(deposition_id)
            or deposition_id <= 0
            or deposition_id in result.values()
        ):
            raise ZenodoReleaseError(
                "Adopted draft deposition IDs must be positive and unique"
            )
        result[record] = deposition_id
    return result


def parse_adopted_drafts(values: Sequence[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    for value in values:
        match = re.fullmatch(r"([a-z0-9_]+)=([1-9][0-9]*)", value)
        if match is None:
            raise ZenodoReleaseError(
                "--adopt-draft must be exactly RECORD=POSITIVE_DEPOSITION_ID"
            )
        record, raw_id = match.groups()
        if record in result:
            raise ZenodoReleaseError(f"Adopted draft repeated for record: {record}")
        result[record] = int(raw_id)
    return _validated_adopted_drafts(result)


def initialize_state(
    request_path: Path,
    archive_root: Path,
    *,
    plan_path: Path = HERE / "release_plan.json",
    audit_path: Path = PAPER / control.EXPECTED_AUTHOR_AUDIT,
    contract_script: Path = PAPER / "scripts/utilities/01_validate_bundle.py",
    adopt_drafts: Mapping[str, int] | None = None,
) -> dict[str, Any]:
    request = read_json(request_path, "release request")
    plan, _ = control.load_validated_plan(plan_path, audit_path)
    _validate_request(request, plan)
    lexical_root = archive_root.expanduser().absolute()
    if lexical_root.is_symlink() or not lexical_root.is_dir():
        raise ZenodoReleaseError(f"Archive root is missing or linked: {lexical_root}")
    root = lexical_root.resolve()
    summary = control.verify_archive_set(root, contract_script=contract_script)
    binding = request.get("archive_manifest", {})
    if (
        binding.get("sha256") != summary["manifest_sha256"]
        or binding.get("tree_sha256") != summary["tree_sha256"]
        or binding.get("tree_files") != summary["tree_files"]
    ):
        raise ZenodoReleaseError(
            "Release request is bound to a different archive manifest"
        )
    expected_files = _expected_files(request, root, summary)
    # A quota cannot be allocated or truthfully confirmed until Zenodo has
    # created the draft and assigned its deposition ID. Derive fresh,
    # deposition-bound quota state below; never treat a pre-draft approval flag
    # as upload authority.
    request_payload = canonical_json(request)
    release_id = "apotome-" + sha256_bytes(request_payload)[:24]
    github = request["github"]
    owner = github.get("owner")
    if not isinstance(owner, str) or not re.fullmatch(
        r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})", owner
    ):
        raise ZenodoReleaseError(
            "GitHub owner must be fixed in the reviewed request before initialization"
        )
    github_url = f"https://github.com/{github['owner']}/{github['repository_name']}"
    adopted = _validated_adopted_drafts(adopt_drafts)
    records: dict[str, Any] = {}
    for name in RECORDS:
        source = request["zenodo_records"][name]
        input_metadata = source["metadata"]
        records[name] = {
            "title": input_metadata["title"],
            "description": input_metadata["description"],
            "upload_type": _record_resource_type(input_metadata["resource_type"]),
            "publication_date": input_metadata["publication_date"],
            "creators": [{"name": row["name"]} for row in input_metadata["creators"]],
            "license": EXPECTED_LICENSES[name],
            "article_doi": _normal_doi(source.get("article_doi"), optional=True),
            "files": expected_files[name],
            "deposition_id": adopted.get(name),
            "reserved_doi": None,
            "concept_record_id": None,
            "bucket_url": None,
            "metadata_synced": False,
            "remote_verified": False,
            "published": False,
            "public_url": None,
            "quota": {
                "required_record_bytes": sum(
                    item["size_bytes"] for item in expected_files[name]
                ),
                "confirmation_required": sum(
                    item["size_bytes"] for item in expected_files[name]
                )
                > archives.ZENODO_DEFAULT_RECORD_BYTES,
                "confirmed_for_deposition_id": None,
                "allocated_record_bytes": (
                    archives.ZENODO_DEFAULT_RECORD_BYTES
                    if sum(item["size_bytes"] for item in expected_files[name])
                    <= archives.ZENODO_DEFAULT_RECORD_BYTES
                    else None
                ),
            },
        }
    state = {
        "schema": STATE_SCHEMA,
        "environment": "production",
        "api_base": PRODUCTION_API,
        "release_id": release_id,
        "request_sha256": sha256_bytes(request_payload),
        "reviewed_gates": {
            "request_operation": request["operation"],
            "rights_confirmations": dict(request["approvals"]["rights_confirmations"]),
            "public_visibility_acknowledged": request["approvals"][
                "public_visibility_acknowledged"
            ],
            "local_clean_room_reproduction_passed": request["approvals"][
                "clean_room_reproduction_passed"
            ],
        },
        "archive_root": str(root),
        "manifest": {
            "sha256": summary["manifest_sha256"],
            "size_bytes": summary["manifest_size_bytes"],
            "tree_sha256": summary["tree_sha256"],
            "tree_files": summary["tree_files"],
        },
        "github_url": github_url,
        "related_identifier_policy": {
            "mutual_reserved_dois": "isSupplementTo",
            "github_repository": "isSupplementTo",
            "article_doi": "isSupplementTo",
        },
        "records": records,
        "remote_download": {
            "status": "NOT_STARTED",
            "root": None,
            "completed_files": [],
            "reconstructed_root": None,
        },
        "reproduction": {
            "status": "NOT_STARTED",
            "log_path": None,
            "log_sha256": None,
            "command": None,
        },
        "created_at": utc_now(),
        "updated_at": utc_now(),
    }
    validate_state(state)
    return state


def _safe_api_url(url: str, *, allow_query: bool = True) -> str:
    try:
        parsed = urllib.parse.urlsplit(url)
        port = parsed.port
    except (TypeError, ValueError) as exc:
        raise ZenodoReleaseError("Zenodo returned a malformed URL") from exc
    if (
        parsed.scheme != "https"
        or parsed.hostname != PRODUCTION_HOST
        or port not in (None, 443)
        or parsed.username is not None
        or parsed.password is not None
        or parsed.fragment
        or not parsed.path.startswith("/api/")
        or (not allow_query and parsed.query)
    ):
        raise ZenodoReleaseError(
            "Refusing a non-production or credential-bearing API URL"
        )
    query_keys = {key.lower() for key, _ in urllib.parse.parse_qsl(parsed.query)}
    if query_keys & reconstruction.SECRET_QUERY_KEYS:
        raise ZenodoReleaseError(
            "Refusing an API URL containing credential query parameters"
        )
    return url


def _safe_public_url(url: str) -> str:
    try:
        parsed = urllib.parse.urlsplit(url)
        port = parsed.port
    except (TypeError, ValueError) as exc:
        raise ZenodoReleaseError("Zenodo returned a malformed public URL") from exc
    if (
        parsed.scheme != "https"
        or parsed.hostname != PRODUCTION_HOST
        or port not in (None, 443)
        or parsed.username is not None
        or parsed.password is not None
        or parsed.fragment
    ):
        raise ZenodoReleaseError("Refusing an invalid public Zenodo URL")
    query_keys = {key.lower() for key, _ in urllib.parse.parse_qsl(parsed.query)}
    if query_keys & reconstruction.SECRET_QUERY_KEYS:
        raise ZenodoReleaseError(
            "Refusing a public Zenodo URL containing credential query parameters"
        )
    return url


class _ProductionRedirect(urllib.request.HTTPRedirectHandler):
    """Prevent urllib from forwarding a bearer token to any other origin."""

    def redirect_request(
        self,
        req: urllib.request.Request,
        fp: BinaryIO,
        code: int,
        msg: str,
        headers: Mapping[str, str],
        newurl: str,
    ) -> urllib.request.Request | None:
        _safe_api_url(newurl)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


@dataclass(frozen=True)
class HTTPResult:
    status: int
    payload: Any
    headers: Mapping[str, str]
    url: str


class ZenodoClient:
    """Small production-only Zenodo client that never serializes/logs its token."""

    def __init__(
        self,
        token: str,
        *,
        timeout_seconds: float = 300.0,
        retries: int = 4,
        opener: Any | None = None,
        sleeper: Callable[[float], None] = time.sleep,
    ) -> None:
        if not isinstance(token, str) or not token.strip() or "\n" in token:
            raise ZenodoReleaseError(
                "Production Zenodo credential is absent or malformed"
            )
        if timeout_seconds <= 0 or retries < 0:
            raise ZenodoReleaseError("HTTP timeout/retry configuration is invalid")
        self.__token = token.strip()
        self.timeout_seconds = timeout_seconds
        self.retries = retries
        self._custom_opener = opener is not None
        self.opener = opener or urllib.request.build_opener(_ProductionRedirect())
        self.sleeper = sleeper

    @property
    def _authorization(self) -> str:
        return "Bearer " + self.__token

    def _redact(self, text: str) -> str:
        return text.replace(self.__token, "[REDACTED]")[:2000]

    def new_download_worker(self) -> "ZenodoClient":
        """Return an isolated real-network client for one downloader thread."""
        if self._custom_opener:
            # Injected openers are test/application-owned; callers can instead
            # supply an explicit download_client_factory when isolation matters.
            return self
        return ZenodoClient(
            self.__token,
            timeout_seconds=self.timeout_seconds,
            retries=self.retries,
            sleeper=self.sleeper,
        )

    def _request(
        self,
        method: str,
        url: str,
        *,
        data: Any = None,
        headers: Mapping[str, str] | None = None,
        retryable: bool,
        expected: Iterable[int],
    ) -> HTTPResult:
        _safe_api_url(url)
        allowed = set(expected)
        request_headers = {
            "Authorization": self._authorization,
            "User-Agent": "Apotome-Zenodo-Production-Executor/1",
            "Accept": "application/json",
        }
        request_headers.update(dict(headers or {}))
        for attempt in range(self.retries + 1):
            request = urllib.request.Request(
                url, data=data, headers=request_headers, method=method
            )
            try:
                response_context = self.opener.open(
                    request, timeout=self.timeout_seconds
                )
                with response_context as response:
                    status = int(getattr(response, "status", response.getcode()))
                    final_url = _safe_api_url(response.geturl())
                    body = response.read(MAX_JSON_BYTES + 1)
                    if len(body) > MAX_JSON_BYTES:
                        raise ZenodoReleaseError(
                            "Zenodo JSON response exceeds safety limit"
                        )
                    payload = json.loads(body.decode("utf-8")) if body else None
                    if status not in allowed:
                        raise ZenodoReleaseError(
                            f"Unexpected Zenodo HTTP status {status}"
                        )
                    return HTTPResult(
                        status, payload, dict(response.headers), final_url
                    )
            except urllib.error.HTTPError as exc:
                status = int(exc.code)
                body = exc.read(16384).decode("utf-8", errors="replace")
                if retryable and status in RETRYABLE_HTTP and attempt < self.retries:
                    self.sleeper(self._retry_delay(attempt, exc.headers))
                    continue
                message = self._redact(body)
                raise ZenodoReleaseError(
                    f"Zenodo HTTP {status}: {message or exc.reason}"
                ) from None
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                if retryable and attempt < self.retries:
                    self.sleeper(min(2**attempt, 30))
                    continue
                raise ZenodoReleaseError(
                    "Zenodo network request failed: " + self._redact(str(exc))
                ) from None
            except (UnicodeError, json.JSONDecodeError) as exc:
                raise ZenodoReleaseError("Zenodo returned invalid JSON") from exc
        raise AssertionError("retry loop exhausted")

    @staticmethod
    def _retry_delay(attempt: int, headers: Mapping[str, str]) -> float:
        raw = headers.get("Retry-After") if headers is not None else None
        if raw is not None:
            try:
                return min(max(float(raw), 0.0), 300.0)
            except ValueError:
                pass
        return float(min(2**attempt, 30))

    def get_deposition(self, deposition_id: int | str) -> Mapping[str, Any]:
        result = self._request(
            "GET",
            f"{PRODUCTION_API}/deposit/depositions/{int(deposition_id)}",
            retryable=True,
            expected={200},
        )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo deposition response is not an object")
        return result.payload

    def get_public_record(self, record_id: int | str) -> Mapping[str, Any]:
        result = self._request(
            "GET",
            f"{PRODUCTION_API}/records/{int(record_id)}",
            retryable=True,
            expected={200},
        )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo public-record response is invalid")
        return result.payload

    def find_drafts(self, release_id: str) -> list[Mapping[str, Any]]:
        query = urllib.parse.urlencode(
            {"status": "draft", "size": 100, "q": f'keywords:"{release_id}"'}
        )
        result = self._request(
            "GET",
            f"{PRODUCTION_API}/deposit/depositions?{query}",
            retryable=True,
            expected={200},
        )
        if not isinstance(result.payload, list) or any(
            not isinstance(row, Mapping) for row in result.payload
        ):
            raise ZenodoReleaseError("Zenodo draft search response is invalid")
        return list(result.payload)

    def create_deposition(self, metadata: Mapping[str, Any]) -> Mapping[str, Any]:
        payload = canonical_json({"metadata": metadata})
        result = self._request(
            "POST",
            f"{PRODUCTION_API}/deposit/depositions",
            data=payload,
            headers={"Content-Type": "application/json"},
            retryable=False,
            expected={201},
        )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo draft creation response is invalid")
        return result.payload

    def update_deposition(
        self, deposition_id: int | str, metadata: Mapping[str, Any]
    ) -> Mapping[str, Any]:
        payload = canonical_json({"metadata": metadata})
        result = self._request(
            "PUT",
            f"{PRODUCTION_API}/deposit/depositions/{int(deposition_id)}",
            data=payload,
            headers={"Content-Type": "application/json"},
            retryable=True,
            expected={200},
        )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo metadata response is invalid")
        return result.payload

    def upload_file(
        self, bucket_url: str, filename: str, path: Path, size_bytes: int
    ) -> Mapping[str, Any]:
        _safe_api_url(bucket_url, allow_query=False)
        if PurePosixPath(filename).name != filename:
            raise ZenodoReleaseError(f"Unsafe upload filename: {filename!r}")
        url = bucket_url.rstrip("/") + "/" + urllib.parse.quote(filename, safe="")
        with path.open("rb") as handle:
            result = self._request(
                "PUT",
                url,
                data=handle,
                headers={
                    "Content-Type": "application/octet-stream",
                    "Content-Length": str(size_bytes),
                },
                retryable=False,
                expected={200, 201},
            )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo upload response is invalid")
        return result.payload

    def publish(self, deposition_id: int | str) -> Mapping[str, Any]:
        result = self._request(
            "POST",
            f"{PRODUCTION_API}/deposit/depositions/{int(deposition_id)}/actions/publish",
            data=b"",
            retryable=False,
            expected={202},
        )
        if not isinstance(result.payload, Mapping):
            raise ZenodoReleaseError("Zenodo publish response is invalid")
        return result.payload

    def download_file(
        self,
        url: str,
        target: Path,
        *,
        expected_size: int,
        expected_sha256: str,
    ) -> None:
        """Download to a .part file and resume only on a valid HTTP byte range."""
        _safe_api_url(url)
        partial = target.with_name(target.name + ".part")
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if (
                target.is_symlink()
                or target.stat().st_size != expected_size
                or sha256_file(target) != expected_sha256
            ):
                raise ZenodoReleaseError(f"Existing downloaded file differs: {target}")
            return
        if partial.is_symlink() or (partial.exists() and not partial.is_file()):
            raise ZenodoReleaseError(f"Unsafe partial download path: {partial}")
        start = partial.stat().st_size if partial.exists() else 0
        if start > expected_size:
            raise ZenodoReleaseError(
                f"Partial download exceeds expected size: {partial}"
            )
        headers = {
            "Authorization": self._authorization,
            "User-Agent": "Apotome-Zenodo-Production-Executor/1",
            "Accept-Encoding": "identity",
        }
        if start:
            headers["Range"] = f"bytes={start}-"
        request = urllib.request.Request(url, headers=headers, method="GET")
        try:
            response_context = self.opener.open(request, timeout=self.timeout_seconds)
            with response_context as response:
                status = int(getattr(response, "status", response.getcode()))
                _safe_api_url(response.geturl())
                if start and status == 206:
                    content_range = response.headers.get("Content-Range", "")
                    if not content_range.startswith(f"bytes {start}-"):
                        raise ZenodoReleaseError(
                            "Zenodo returned an invalid resume range"
                        )
                    mode = "ab"
                elif status == 200:
                    start = 0
                    mode = "wb"
                else:
                    raise ZenodoReleaseError(
                        f"Unexpected Zenodo download HTTP status {status}"
                    )
                count = start
                with partial.open(mode) as output:
                    while True:
                        block = response.read(CHUNK)
                        if not block:
                            break
                        count += len(block)
                        if count > expected_size:
                            raise ZenodoReleaseError(
                                "Downloaded file exceeds expected size"
                            )
                        output.write(block)
                    output.flush()
                    os.fsync(output.fileno())
        except ZenodoReleaseError:
            raise
        except urllib.error.HTTPError as exc:
            body = exc.read(16384).decode("utf-8", errors="replace")
            raise ZenodoReleaseError(
                f"Zenodo download HTTP {exc.code}: {self._redact(body or str(exc.reason))}"
            ) from None
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise ZenodoReleaseError(
                "Zenodo download failed; rerun to resume: " + self._redact(str(exc))
            ) from None
        if (
            partial.stat().st_size != expected_size
            or sha256_file(partial) != expected_sha256
        ):
            raise ZenodoReleaseError(
                f"Downloaded file size/SHA-256 differs; partial retained for diagnosis: {target.name}"
            )
        os.chmod(partial, 0o644)
        os.replace(partial, target)


def load_token(
    token_file: Path | None, environ: Mapping[str, str] | None = None
) -> str:
    environment = os.environ if environ is None else environ
    from_environment = environment.get(TOKEN_ENV, "").strip()
    if token_file is not None and from_environment:
        raise ZenodoReleaseError(
            f"Choose one credential source: --token-file or {TOKEN_ENV}, not both"
        )
    if token_file is None:
        if not from_environment:
            raise ZenodoReleaseError(
                f"Production credential is absent; use a mode-0600 --token-file or {TOKEN_ENV}"
            )
        return from_environment
    path = token_file.expanduser().absolute()
    if path.is_symlink() or not path.is_file():
        raise ZenodoReleaseError(f"Credential file is missing or linked: {path}")
    if stat.S_IMODE(path.stat().st_mode) & 0o077:
        raise ZenodoReleaseError("Credential file must not be group/world accessible")
    if path.stat().st_size > 16384:
        raise ZenodoReleaseError("Credential file is unexpectedly large")
    token = path.read_text(encoding="utf-8").strip()
    if not token or "\n" in token:
        raise ZenodoReleaseError("Credential file must contain exactly one token")
    return token


def _bucket_from(deposition: Mapping[str, Any]) -> str:
    links = deposition.get("links")
    value = links.get("bucket") if isinstance(links, Mapping) else None
    if not isinstance(value, str):
        raise ZenodoReleaseError("Zenodo draft has no large-file bucket URL")
    return _safe_api_url(value, allow_query=False)


def _deposition_id(deposition: Mapping[str, Any]) -> int:
    value = deposition.get("id")
    if not _plain_int(value) or value <= 0:
        raise ZenodoReleaseError("Zenodo returned an invalid deposition ID")
    return value


def _reserved_doi(deposition: Mapping[str, Any]) -> str:
    metadata = deposition.get("metadata", {})
    reserve = (
        metadata.get("prereserve_doi", {}) if isinstance(metadata, Mapping) else {}
    )
    value = reserve.get("doi") if isinstance(reserve, Mapping) else None
    if value is None:
        value = deposition.get("doi")
    return str(_normal_doi(value))


Persist = Callable[[Mapping[str, Any]], None]


def _touch(state: MutableMapping[str, Any]) -> None:
    state["updated_at"] = utc_now()


def _create_metadata(state: Mapping[str, Any], record: str) -> dict[str, Any]:
    row = state["records"][record]
    return {
        "title": row["title"],
        "upload_type": row["upload_type"],
        "description": row["description"],
        "creators": list(row["creators"]),
        "publication_date": row["publication_date"],
        "access_right": "open",
        "license": row["license"],
        "prereserve_doi": True,
        "keywords": [
            "allulose",
            "sucrose",
            "hypothalamus",
            "reproducible research",
            state["release_id"],
        ],
        "notes": (
            "Four-record reproducibility release. Executor marker: "
            + state["release_id"]
        ),
    }


def effective_metadata(state: Mapping[str, Any], record: str) -> dict[str, Any]:
    """Render the exact reviewed metadata after every reserved DOI is known."""
    if record not in RECORDS:
        raise ZenodoReleaseError(f"Unknown logical record: {record}")
    metadata = _create_metadata(state, record)
    identifiers: list[dict[str, str]] = []
    resource_types = {
        "software_core": "software",
        "main_figure_data": "dataset",
        "supplementary_figure_data": "dataset",
        "figs3_wsi": "dataset",
    }
    for other in RECORDS:
        if other == record:
            continue
        doi = state["records"][other].get("reserved_doi")
        if not isinstance(doi, str):
            raise ZenodoReleaseError(
                "All four DOIs must be reserved before metadata sync"
            )
        identifiers.append(
            {
                "identifier": str(_normal_doi(doi)),
                "relation": state["related_identifier_policy"]["mutual_reserved_dois"],
                "resource_type": resource_types[other],
            }
        )
    identifiers.append(
        {
            "identifier": state["github_url"],
            "relation": state["related_identifier_policy"]["github_repository"],
            "resource_type": "software",
        }
    )
    article = state["records"][record].get("article_doi")
    if article is not None:
        identifiers.append(
            {
                "identifier": str(_normal_doi(article)),
                "relation": state["related_identifier_policy"]["article_doi"],
                "resource_type": "publication-article",
            }
        )
    metadata["related_identifiers"] = identifiers
    return metadata


def _metadata_license(metadata: Mapping[str, Any]) -> Any:
    value = metadata.get("license")
    return value.get("id") if isinstance(value, Mapping) else value


def _description_matches(actual: Any, expected: Any) -> bool:
    """Accept only Zenodo's reversible ampersand serialization.

    Zenodo stores descriptions as HTML and returns a literal ampersand as
    ``&amp;``. Keep this compatibility rule deliberately narrower than a
    general HTML unescape so alternate entities, double escaping, markup, and
    semantic text changes continue to fail closed.
    """
    if not isinstance(actual, str) or not isinstance(expected, str):
        return False
    return actual == expected or (
        "&" in expected
        and html.unescape(expected) == expected
        and actual == expected.replace("&", "&amp;")
        and html.unescape(actual) == expected
    )


def _normalized_identifier(value: Any) -> str:
    if not isinstance(value, str):
        return ""
    lowered = value.lower()
    if lowered.startswith("https://doi.org/"):
        return value[len("https://doi.org/") :]
    return value


def verify_metadata(
    deposition: Mapping[str, Any],
    expected: Mapping[str, Any],
    *,
    require_related: bool,
) -> None:
    actual = deposition.get("metadata")
    if not isinstance(actual, Mapping):
        raise ZenodoReleaseError("Zenodo deposition metadata is absent")
    scalar = ("title", "upload_type", "publication_date", "access_right")
    if any(
        actual.get(key) != expected.get(key) for key in scalar
    ) or not _description_matches(
        actual.get("description"), expected.get("description")
    ):
        raise ZenodoReleaseError("Zenodo metadata differs from the reviewed release")
    if _metadata_license(actual) != expected.get("license"):
        raise ZenodoReleaseError("Zenodo license differs from the reviewed release")
    actual_creators = actual.get("creators", [])
    expected_creators = _creator_names(expected.get("creators"))
    if (
        expected_creators != list(control.EXPECTED_ZENODO_NAMES)
        or _creator_names(actual_creators) != expected_creators
    ):
        raise ZenodoReleaseError("Zenodo creator names/order differ")
    if state_marker := expected.get("keywords", [])[-1:]:
        if state_marker[0] not in actual.get("keywords", []):
            raise ZenodoReleaseError("Zenodo release marker is absent")
    if require_related:
        actual_related = actual.get("related_identifiers")
        if not isinstance(actual_related, list):
            raise ZenodoReleaseError("Zenodo related identifiers are absent")
        normalized_actual = {
            (
                _normalized_identifier(row.get("identifier")),
                row.get("relation"),
                row.get("resource_type"),
            )
            for row in actual_related
            if isinstance(row, Mapping)
        }
        normalized_expected = {
            (
                _normalized_identifier(row["identifier"]),
                row["relation"],
                row["resource_type"],
            )
            for row in expected["related_identifiers"]
        }
        if normalized_actual != normalized_expected:
            raise ZenodoReleaseError("Zenodo DOI/repository cross-links differ")


def _verify_explicit_adoption(
    deposition: Mapping[str, Any], state: Mapping[str, Any], record: str
) -> None:
    row = state["records"][record]
    if deposition.get("submitted") is not False:
        raise ZenodoReleaseError(
            f"Adopted deposition is not an unpublished draft: {record}"
        )
    if _deposition_id(deposition) != row["deposition_id"]:
        raise ZenodoReleaseError(f"Adopted deposition ID differs for {record}")
    if normal_remote_files(deposition):
        raise ZenodoReleaseError(f"Adopted draft must be empty: {record}")
    actual = deposition.get("metadata")
    expected = _create_metadata(state, record)
    if not isinstance(actual, Mapping):
        raise ZenodoReleaseError(f"Adopted draft metadata is absent: {record}")
    stable = ("title", "upload_type", "publication_date", "access_right")
    if (
        any(actual.get(key) != expected.get(key) for key in stable)
        or _metadata_license(actual) != expected["license"]
        or _creator_names(actual.get("creators"))
        != _creator_names(expected.get("creators"))
    ):
        raise ZenodoReleaseError(f"Adopted draft identity differs for {record}")
    expected_keywords = expected["keywords"][:-1]
    actual_keywords = actual.get("keywords")
    if (
        not isinstance(actual_keywords, list)
        or actual_keywords[:-1] != expected_keywords
        or len(actual_keywords) != len(expected["keywords"])
        or not RELEASE_ID_RE.fullmatch(str(actual_keywords[-1]))
    ):
        raise ZenodoReleaseError(
            f"Adopted draft release marker differs unsafely: {record}"
        )
    if not re.fullmatch(
        r"(?:Three|Four)-record reproducibility release\. Executor marker: "
        r"apotome-[0-9a-f]{24}",
        str(actual.get("notes", "")),
    ):
        raise ZenodoReleaseError(f"Adopted draft notes differ unsafely: {record}")
    _reserved_doi(deposition)
    _bucket_from(deposition)


def _adoptable_drafts(
    drafts: Sequence[Mapping[str, Any]], state: Mapping[str, Any], record: str
) -> list[Mapping[str, Any]]:
    title = state["records"][record]["title"]
    marker = state["release_id"]
    matches: list[Mapping[str, Any]] = []
    for draft in drafts:
        metadata = draft.get("metadata", {})
        if (
            isinstance(metadata, Mapping)
            and metadata.get("title") == title
            and marker in metadata.get("keywords", [])
            and draft.get("submitted") is False
        ):
            matches.append(draft)
    return matches


def create_and_link_drafts(
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    *,
    execute: bool,
    production: bool,
) -> None:
    require_mutation_gate("drafts", execute=execute, production=production)
    validate_state(state)
    drafts = client.find_drafts(state["release_id"])
    fully_linked_at_start = all(
        state["records"][record].get("reserved_doi") is not None for record in RECORDS
    )
    explicit_adoptions: dict[str, Mapping[str, Any]] = {}
    for record in RECORDS:
        row = state["records"][record]
        if row["deposition_id"] is not None and row.get("reserved_doi") is None:
            deposition = client.get_deposition(row["deposition_id"])
            _verify_explicit_adoption(deposition, state, record)
            explicit_adoptions[record] = deposition

    for record in RECORDS:
        row = state["records"][record]
        deposition: Mapping[str, Any]
        if row["deposition_id"] is not None:
            deposition = explicit_adoptions.get(record) or client.get_deposition(
                row["deposition_id"]
            )
            if record in explicit_adoptions:
                expected_base = _create_metadata(state, record)
                try:
                    verify_metadata(deposition, expected_base, require_related=False)
                except ZenodoReleaseError:
                    deposition = client.update_deposition(
                        row["deposition_id"], expected_base
                    )
                    verify_metadata(deposition, expected_base, require_related=False)
                row["reserved_doi"] = _reserved_doi(deposition)
                row["concept_record_id"] = deposition.get("conceptrecid")
                row["bucket_url"] = _bucket_from(deposition)
                _touch(state)
                persist(state)
        else:
            matches = _adoptable_drafts(drafts, state, record)
            if len(matches) > 1:
                raise ZenodoReleaseError(
                    f"Multiple matching Zenodo drafts exist for {record}; adopt manually"
                )
            if matches:
                deposition = matches[0]
            else:
                # Creation POST is deliberately never retried automatically.
                deposition = client.create_deposition(_create_metadata(state, record))
            row["deposition_id"] = _deposition_id(deposition)
            row["reserved_doi"] = _reserved_doi(deposition)
            row["concept_record_id"] = deposition.get("conceptrecid")
            row["bucket_url"] = _bucket_from(deposition)
            _touch(state)
            persist(state)
        if deposition.get("submitted") is not False:
            raise ZenodoReleaseError(f"Expected an unpublished draft for {record}")
        if _deposition_id(deposition) != row["deposition_id"]:
            raise ZenodoReleaseError(f"Zenodo deposition ID differs for {record}")
        row["reserved_doi"] = _reserved_doi(deposition)
        if fully_linked_at_start:
            expected_current = effective_metadata(state, record)
            require_related = True
        else:
            expected_current = _create_metadata(state, record)
            require_related = False
        try:
            verify_metadata(
                deposition, expected_current, require_related=require_related
            )
        except ZenodoReleaseError:
            deposition = client.update_deposition(
                row["deposition_id"], expected_current
            )
            verify_metadata(
                deposition, expected_current, require_related=require_related
            )
        _touch(state)
        persist(state)

    dois = [state["records"][name]["reserved_doi"] for name in RECORDS]
    if len(set(dois)) != len(RECORDS):
        raise ZenodoReleaseError("Zenodo reserved duplicate DOIs")
    for record in RECORDS:
        row = state["records"][record]
        expected = effective_metadata(state, record)
        deposition = client.get_deposition(row["deposition_id"])
        try:
            verify_metadata(deposition, expected, require_related=True)
        except ZenodoReleaseError:
            deposition = client.update_deposition(row["deposition_id"], expected)
            verify_metadata(deposition, expected, require_related=True)
        row["metadata_synced"] = True
        row["metadata_sha256"] = sha256_bytes(canonical_json(expected))
        _touch(state)
        persist(state)


def confirm_quota(
    state: MutableMapping[str, Any],
    record: str,
    *,
    deposition_id: int,
    allocated_record_bytes: int,
    confirmed: bool,
) -> None:
    if not confirmed:
        raise ZenodoReleaseError(
            "Quota confirmation requires --confirm-quota-allocation"
        )
    if record not in RECORDS:
        raise ZenodoReleaseError(f"Unknown logical record: {record}")
    validate_state(state)
    row = state["records"][record]
    if (
        not _plain_int(deposition_id)
        or deposition_id <= 0
        or row["deposition_id"] != deposition_id
    ):
        raise ZenodoReleaseError("Quota confirmation deposition ID differs from state")
    required, confirmation_required = _validated_quota_state(record, row)
    if not confirmation_required:
        raise ZenodoReleaseError(
            f"Zenodo quota confirmation is not required for small record {record}"
        )
    if (
        not _plain_int(allocated_record_bytes)
        or allocated_record_bytes < required
        or allocated_record_bytes > archives.ZENODO_MAX_RECORD_BYTES
    ):
        raise ZenodoReleaseError(
            f"Allocated quota must cover {required:,} bytes and not exceed 200,000,000,000"
        )
    row["quota"]["confirmed_for_deposition_id"] = deposition_id
    row["quota"]["allocated_record_bytes"] = allocated_record_bytes
    _validated_quota_state(record, row)
    _touch(state)


def _quota_ready(record: str, row: Mapping[str, Any]) -> None:
    required, confirmation_required = _validated_quota_state(record, row)
    quota = row["quota"]
    if confirmation_required and quota["confirmed_for_deposition_id"] is None:
        raise ZenodoReleaseError(
            f"Zenodo storage quota is not confirmed for {record} draft "
            f"{row['deposition_id']}: allocate at least {required:,} total bytes "
            "with Manage storage, then run confirm-quota with that deposition ID"
        )


def _local_path(state: Mapping[str, Any], file_row: Mapping[str, Any]) -> Path:
    root = Path(state["archive_root"]).resolve()
    target = root.joinpath(*PurePosixPath(file_row["relative_local_path"]).parts)
    if (
        target.is_symlink()
        or not target.is_file()
        or not target.resolve().is_relative_to(root)
    ):
        raise ZenodoReleaseError(f"Missing or unsafe upload source: {target}")
    return target


def _remote_matches(remote: Mapping[str, Any], expected: Mapping[str, Any]) -> bool:
    return (
        remote["size_bytes"] == expected["size_bytes"]
        and expected.get("md5") is not None
        and remote["md5"] == expected["md5"]
    )


def _verify_remote_inventory(
    state: Mapping[str, Any], record: str, deposition: Mapping[str, Any]
) -> None:
    row = state["records"][record]
    actual = normal_remote_files(deposition)
    expected = {item["name"]: item for item in expected_remote_items(state, record)}
    if set(actual) != set(expected):
        raise ZenodoReleaseError(
            "Zenodo file inventory differs: "
            f"missing={sorted(set(expected)-set(actual))} "
            f"extra={sorted(set(actual)-set(expected))}"
        )
    for name, item in expected.items():
        if not _remote_matches(actual[name], item):
            raise ZenodoReleaseError(f"Zenodo size/MD5 differs: {name}")
        item["uploaded"] = True
        item["remote_verified"] = True
    for item in row["files"]:
        item["uploaded"] = True
        item["remote_verified"] = True
    row["remote_verified"] = True


def _check_upload_response(
    response: Mapping[str, Any], expected: Mapping[str, Any]
) -> None:
    wrapped = normal_remote_files({"files": [response]})
    if set(wrapped) != {expected["name"]}:
        raise ZenodoReleaseError("Zenodo upload response filename differs")
    if not _remote_matches(wrapped[expected["name"]], expected):
        raise ZenodoReleaseError("Zenodo upload response size/MD5 differs")


def upload_all(
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    *,
    execute: bool,
    production: bool,
) -> None:
    require_mutation_gate("upload", execute=execute, production=production)
    validate_state(state)
    if load_segmented_transport(state) is not None:
        raise ZenodoReleaseError(
            "The canonical uploader cannot write byte-range segments; use the "
            "reviewed segmented uploader, then run verify-remote"
        )
    # Preflight every record before the first network request or file upload so
    # a later data-record quota failure cannot leave an avoidable partial upload.
    for record in RECORDS:
        row = state["records"][record]
        if not row.get("metadata_synced"):
            raise ZenodoReleaseError(f"Zenodo metadata is not synchronized: {record}")
        _quota_ready(record, row)
    for record in RECORDS:
        row = state["records"][record]
        deposition = client.get_deposition(row["deposition_id"])
        if deposition.get("submitted") is not False:
            raise ZenodoReleaseError(f"Cannot upload to a published record: {record}")
        verify_metadata(
            deposition, effective_metadata(state, record), require_related=True
        )
        bucket = _bucket_from(deposition)
        remote = normal_remote_files(deposition)
        expected_names = {item["name"] for item in row["files"]}
        extra = set(remote) - expected_names
        if extra:
            raise ZenodoReleaseError(
                f"Unexpected existing Zenodo files for {record}: {sorted(extra)}"
            )
        for item in row["files"]:
            path = _local_path(state, item)
            sha256, md5, size = hashes_file(path)
            if size != item["size_bytes"] or sha256 != item["sha256"]:
                raise ZenodoReleaseError(
                    f"Local upload source size/SHA-256 changed: {path}"
                )
            item["md5"] = md5
            current = remote.get(item["name"])
            if current is not None:
                if not _remote_matches(current, item):
                    raise ZenodoReleaseError(
                        f"Existing Zenodo file conflicts; it will not be overwritten: {item['name']}"
                    )
                item["uploaded"] = True
                item["remote_verified"] = True
                _touch(state)
                persist(state)
                continue
            try:
                response = client.upload_file(
                    bucket, item["name"], path, item["size_bytes"]
                )
                _check_upload_response(response, item)
            except ZenodoReleaseError as exc:
                # A connection may fail after Zenodo committed the PUT. Reconcile once;
                # never blindly retry a multi-gigabyte write or overwrite a remote name.
                reconciled = client.get_deposition(row["deposition_id"])
                observed = normal_remote_files(reconciled).get(item["name"])
                if observed is None or not _remote_matches(observed, item):
                    required = row["quota"]["required_record_bytes"]
                    raise ZenodoReleaseError(
                        f"Upload failed for {record}/{item['name']}. The record needs "
                        f"{required:,} bytes total; verify its Manage storage allocation. "
                        f"Underlying error: {exc}"
                    ) from None
            item["uploaded"] = True
            item["remote_verified"] = True
            _touch(state)
            persist(state)
            deposition = client.get_deposition(row["deposition_id"])
            remote = normal_remote_files(deposition)
        deposition = client.get_deposition(row["deposition_id"])
        _verify_remote_inventory(state, record, deposition)
        verify_metadata(
            deposition, effective_metadata(state, record), require_related=True
        )
        _touch(state)
        persist(state)


def verify_remote(
    state: MutableMapping[str, Any], client: Any, persist: Persist
) -> None:
    """Read-only server reconciliation using exact filenames, bytes, and MD5."""
    validate_state(state)
    for record in RECORDS:
        row = state["records"][record]
        if row.get("deposition_id") is None:
            raise ZenodoReleaseError(f"Zenodo draft is absent: {record}")
        if any(
            item.get("md5") is None for item in expected_remote_items(state, record)
        ):
            raise ZenodoReleaseError(
                f"Local MD5 inventory is incomplete; run the upload phase: {record}"
            )
        deposition = client.get_deposition(row["deposition_id"])
        verify_metadata(
            deposition, effective_metadata(state, record), require_related=True
        )
        _verify_remote_inventory(state, record, deposition)
        _touch(state)
        persist(state)


def _download_url(
    deposition: Mapping[str, Any], remote: Mapping[str, Any], filename: str
) -> str:
    links = remote.get("links", {})
    if isinstance(links, Mapping):
        for key in ("download", "self", "content"):
            value = links.get(key)
            if isinstance(value, str):
                return _safe_api_url(value)
    return (
        _bucket_from(deposition).rstrip("/")
        + "/"
        + urllib.parse.quote(filename, safe="")
    )


def _download_relative(file_row: Mapping[str, Any]) -> str:
    return str(file_row["relative_local_path"])


def _closed_download_inventory(root: Path, expected: set[str]) -> None:
    observed: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ZenodoReleaseError(
                f"Downloaded archive set contains a symlink: {path}"
            )
        if path.is_file() and path.name.endswith(".part"):
            raise ZenodoReleaseError(
                f"Partial download remains after verification: {path}"
            )
        if path.is_file():
            observed.add(path.relative_to(root).as_posix())
        elif not path.is_file() and not path.is_dir():
            raise ZenodoReleaseError(
                f"Downloaded archive set has an unsafe entry: {path}"
            )
    if observed != expected:
        raise ZenodoReleaseError(
            "Downloaded archive inventory differs: "
            f"missing={sorted(expected-observed)} extra={sorted(observed-expected)}"
        )


def _verified_download(path: Path, expected: Mapping[str, Any], label: str) -> None:
    if (
        path.is_symlink()
        or not path.is_file()
        or path.stat().st_size != expected["size_bytes"]
        or sha256_file(path) != expected["sha256"]
    ):
        raise ZenodoReleaseError(f"Downloaded size/SHA-256 differs: {label}")


def _sha256_file_range(path: Path, offset: int, size: int) -> str:
    digest = hashlib.sha256()
    remaining = size
    with path.open("rb") as handle:
        handle.seek(offset)
        while remaining:
            block = handle.read(min(CHUNK, remaining))
            if not block:
                raise ZenodoReleaseError(f"Unexpected EOF in partial archive: {path}")
            remaining -= len(block)
            digest.update(block)
    return digest.hexdigest()


def _download_direct_item(
    *,
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    target: Path,
    deposition: Mapping[str, Any],
    remote: Mapping[str, Mapping[str, Any]],
    item: Mapping[str, Any],
    completed: set[str],
) -> None:
    relative = str(item["relative_local_path"])
    try:
        safe = archives.safe_relative(relative)
    except (TypeError, archives.ArchiveError) as exc:
        raise ZenodoReleaseError(f"Unsafe roundtrip path: {relative!r}") from exc
    destination = target.joinpath(*PurePosixPath(safe).parts)
    if relative in completed:
        _verified_download(destination, item, relative)
        return
    client.download_file(
        _download_url(deposition, remote[item["name"]], item["name"]),
        destination,
        expected_size=item["size_bytes"],
        expected_sha256=item["sha256"],
    )
    _verified_download(destination, item, relative)
    completed.add(relative)
    state["remote_download"]["completed_files"] = sorted(completed)
    _touch(state)
    persist(state)


def _reconcile_partial_archive(
    part: Path,
    segments: Sequence[Mapping[str, Any]],
    completed_count: int,
) -> int:
    expected_size = sum(item["size_bytes"] for item in segments[:completed_count])
    if not part.exists():
        if expected_size:
            raise ZenodoReleaseError(
                f"State records missing partial archive bytes: {part}"
            )
        return completed_count
    if part.is_symlink() or not part.is_file():
        raise ZenodoReleaseError(f"Partial archive is unsafe: {part}")
    observed_size = part.stat().st_size
    if observed_size < expected_size:
        raise ZenodoReleaseError(f"Partial archive is shorter than state: {part}")
    offset = expected_size
    reconciled = completed_count
    while offset < observed_size and reconciled < len(segments):
        segment = segments[reconciled]
        end = offset + segment["size_bytes"]
        if end > observed_size:
            # Appending is fsynced before its completion enters state, so a
            # crash or ENOSPC can leave only the uncommitted trailing segment
            # incomplete.  Roll that tail back to the last verified boundary;
            # the retained segment cache (or a fresh download) can replay it.
            with part.open("r+b") as handle:
                handle.truncate(offset)
                handle.flush()
                os.fsync(handle.fileno())
            observed_size = offset
            break
        if _sha256_file_range(part, offset, segment["size_bytes"]) != segment["sha256"]:
            raise ZenodoReleaseError(f"Unrecorded partial archive bytes differ: {part}")
        offset = end
        reconciled += 1
    if offset != observed_size:
        raise ZenodoReleaseError(f"Partial archive size is not segment-aligned: {part}")
    return reconciled


def _append_verified_segment(part: Path, segment_path: Path, offset: int) -> None:
    part.parent.mkdir(parents=True, exist_ok=True)
    if part.exists() and part.stat().st_size != offset:
        raise ZenodoReleaseError(f"Partial archive offset changed: {part}")
    if not part.exists() and offset != 0:
        raise ZenodoReleaseError(f"Partial archive is absent at nonzero offset: {part}")
    with part.open("ab") as output, segment_path.open("rb") as source:
        for block in iter(lambda: source.read(CHUNK), b""):
            output.write(block)
        output.flush()
        os.fsync(output.fileno())


def _discard_completed_segment_cache(
    record_temp: Path, segments: Sequence[Mapping[str, Any]]
) -> None:
    if not record_temp.exists():
        return
    if record_temp.is_symlink() or not record_temp.is_dir():
        raise ZenodoReleaseError(f"Transport cache is unsafe: {record_temp}")
    for segment in segments:
        segment_path = record_temp / segment["name"]
        if segment_path.exists():
            _verified_download(segment_path, segment, segment["name"])
            segment_path.unlink()
        partial = segment_path.with_name(segment_path.name + ".part")
        if partial.exists():
            if partial.is_symlink() or not partial.is_file():
                raise ZenodoReleaseError(
                    f"Transport partial cache is unsafe: {partial}"
                )
            partial.unlink()


@dataclass(frozen=True)
class _SegmentedArchivePlan:
    record: str
    archive: Mapping[str, Any]
    destination: Path
    part: Path
    record_temp: Path
    segments: tuple[Mapping[str, Any], ...]
    segment_ids: tuple[str, ...]
    prefix: int
    offset: int
    expected_items: Mapping[str, Mapping[str, Mapping[str, Any]]]
    depositions: Mapping[str, Mapping[str, Any]]
    remotes: Mapping[str, Mapping[str, Mapping[str, Any]]]


def _run_segmented_archive_plan(
    plan: _SegmentedArchivePlan,
    client: Any,
    events: "queue.Queue[tuple[str, str]]",
    stop: threading.Event,
) -> None:
    """Download one archive; all shared-state commits remain with the caller."""
    relative = str(plan.archive["relative_local_path"])
    offset = plan.offset
    for index in range(plan.prefix, len(plan.segments)):
        if stop.is_set():
            return
        segment = plan.segments[index]
        storage_record = segment.get("storage_record", plan.record)
        expected = plan.expected_items[storage_record][segment["name"]]
        segment_path = plan.record_temp / segment["name"]
        if segment_path.exists():
            _verified_download(segment_path, expected, segment["name"])
        else:
            client.download_file(
                _download_url(
                    plan.depositions[storage_record],
                    plan.remotes[storage_record][segment["name"]],
                    segment["name"],
                ),
                segment_path,
                expected_size=segment["size_bytes"],
                expected_sha256=segment["sha256"],
            )
            _verified_download(segment_path, expected, segment["name"])
        _append_verified_segment(plan.part, segment_path, offset)
        if plan.part.stat().st_size != offset + segment["size_bytes"]:
            raise ZenodoReleaseError(
                f"Reassembled partial archive size differs: {relative}"
            )
        segment_path.unlink()
        offset += segment["size_bytes"]
        events.put(("segment", plan.segment_ids[index]))

    sha256, md5, size = hashes_file(plan.part)
    if (
        size != plan.archive["size_bytes"]
        or sha256 != plan.archive["sha256"]
        or md5 != plan.archive["md5"]
    ):
        raise ZenodoReleaseError(f"Reassembled canonical archive differs: {relative}")
    os.replace(plan.part, plan.destination)
    events.put(("archive", relative))


def _download_segmented_archives(
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    *,
    target: Path,
    transport: Mapping[str, Any],
    depositions: Mapping[str, Mapping[str, Any]],
    remotes: Mapping[str, Mapping[str, Mapping[str, Any]]],
    completed: set[str],
    download_workers: int,
    download_client_factory: Callable[[], Any] | None = None,
) -> tuple[set[str], int]:
    session = state["remote_download"]
    completed_segments = set(session.get("completed_segments", []))
    all_segment_ids = {
        f"{record}/{segment['name']}"
        for record in RECORDS
        for archive in transport["records"][record]["archives"]
        for segment in archive["segments"]
    }
    if not completed_segments <= all_segment_ids:
        raise ZenodoReleaseError("Roundtrip state contains unknown transport segments")
    session["completed_segments"] = sorted(completed_segments)

    expected_paths = {
        item["relative_local_path"]
        for record in RECORDS
        for item in state["records"][record]["files"]
    }
    expected_paths.add(TRANSPORT_MANIFEST_FILENAME)
    if not completed <= expected_paths:
        raise ZenodoReleaseError("Roundtrip state contains unknown completed files")

    remote_items = {record: expected_remote_items(state, record) for record in RECORDS}
    expected_items = {
        record: {item["name"]: item for item in remote_items[record]}
        for record in RECORDS
    }
    for item in remote_items["software_core"]:
        if item.get("transport_kind") not in {
            "canonical_manifest",
            "transport_manifest",
        }:
            continue
        _download_direct_item(
            state=state,
            client=client,
            persist=persist,
            target=target,
            deposition=depositions["software_core"],
            remote=remotes["software_core"],
            item=item,
            completed=completed,
        )

    temporary_root = target / ".transport-segments"
    if temporary_root.exists() and (
        temporary_root.is_symlink() or not temporary_root.is_dir()
    ):
        raise ZenodoReleaseError(
            f"Transport temporary path is unsafe: {temporary_root}"
        )
    temporary_root.mkdir(mode=0o700, exist_ok=True)
    for record in RECORDS:
        record_temp = temporary_root / record
        if record_temp.exists() and (
            record_temp.is_symlink() or not record_temp.is_dir()
        ):
            raise ZenodoReleaseError(f"Transport record cache is unsafe: {record_temp}")
        record_temp.mkdir(mode=0o700, exist_ok=True)

    plans: list[_SegmentedArchivePlan] = []
    preflight_changed = False
    archive_segments: dict[str, set[str]] = {}
    for record in RECORDS:
        record_temp = temporary_root / record
        for archive in transport["records"][record]["archives"]:
            relative = archive["relative_local_path"]
            destination = target.joinpath(*PurePosixPath(relative).parts)
            part = destination.with_name(f".{destination.name}.part")
            segments = tuple(archive["segments"])
            segment_ids = tuple(f"{record}/{item['name']}" for item in segments)
            archive_segments[relative] = set(segment_ids)
            prefix = 0
            while prefix < len(segments) and segment_ids[prefix] in completed_segments:
                prefix += 1
            if any(
                identifier in completed_segments for identifier in segment_ids[prefix:]
            ):
                raise ZenodoReleaseError(
                    f"Completed transport segments are not a prefix: {relative}"
                )

            if destination.exists():
                if part.exists():
                    raise ZenodoReleaseError(
                        f"Both final and partial archive exist: {relative}"
                    )
                _verified_download(destination, archive, relative)
                if hashes_file(destination)[1] != archive["md5"]:
                    raise ZenodoReleaseError(
                        f"Reassembled archive MD5 differs: {relative}"
                    )
                before_segments = len(completed_segments)
                before_files = len(completed)
                completed_segments.update(segment_ids)
                completed.add(relative)
                preflight_changed |= (
                    len(completed_segments) != before_segments
                    or len(completed) != before_files
                )
                _discard_completed_segment_cache(record_temp, segments)
                continue

            reconciled = _reconcile_partial_archive(part, segments, prefix)
            if reconciled > prefix:
                completed_segments.update(segment_ids[prefix:reconciled])
                preflight_changed = True
                prefix = reconciled
            _discard_completed_segment_cache(record_temp, segments[:prefix])
            plans.append(
                _SegmentedArchivePlan(
                    record=record,
                    archive=archive,
                    destination=destination,
                    part=part,
                    record_temp=record_temp,
                    segments=segments,
                    segment_ids=segment_ids,
                    prefix=prefix,
                    offset=sum(item["size_bytes"] for item in segments[:prefix]),
                    expected_items=expected_items,
                    depositions=depositions,
                    remotes=remotes,
                )
            )

    # Start with one archive from each record before taking a second archive
    # from any record.  This avoids concentrating all initial connections on a
    # single Zenodo bucket while preserving canonical order within each record.
    plans_by_record = {
        record: [plan for plan in plans if plan.record == record] for record in RECORDS
    }
    plans = [
        plans_by_record[record][index]
        for index in range(
            max((len(rows) for rows in plans_by_record.values()), default=0)
        )
        for record in RECORDS
        if index < len(plans_by_record[record])
    ]

    if preflight_changed:
        session["completed_segments"] = sorted(completed_segments)
        session["completed_files"] = sorted(completed)
        _touch(state)
        persist(state)

    factory = download_client_factory
    if factory is None:
        candidate = getattr(client, "new_download_worker", None)
        factory = candidate if callable(candidate) else lambda: client
    worker_local = threading.local()

    def worker_client() -> Any:
        value = getattr(worker_local, "client", None)
        if value is None:
            value = factory()
            if not callable(getattr(value, "download_file", None)):
                raise ZenodoReleaseError(
                    "Download-client factory returned an invalid client"
                )
            worker_local.client = value
        return value

    events: queue.Queue[tuple[str, str]] = queue.Queue()
    stop = threading.Event()

    def run_plan(plan: _SegmentedArchivePlan) -> None:
        _run_segmented_archive_plan(plan, worker_client(), events, stop)

    def apply_event(event: tuple[str, str]) -> None:
        kind, identifier = event
        if kind == "segment":
            if identifier in completed_segments:
                raise ZenodoReleaseError(
                    f"Duplicate transport-segment completion: {identifier}"
                )
            completed_segments.add(identifier)
        elif kind == "archive":
            required = archive_segments.get(identifier)
            if required is None or not required <= completed_segments:
                raise ZenodoReleaseError(
                    f"Canonical archive completed before its segments: {identifier}"
                )
            completed.add(identifier)
        else:
            raise ZenodoReleaseError(f"Unknown downloader event: {kind}")
        session["completed_segments"] = sorted(completed_segments)
        session["completed_files"] = sorted(completed)
        _touch(state)
        persist(state)

    errors: list[Exception] = []
    if plans:
        with ThreadPoolExecutor(
            max_workers=min(download_workers, len(plans)),
            thread_name_prefix="zenodo-roundtrip",
        ) as executor:
            plan_iterator = iter(plans)
            pending = {
                executor.submit(run_plan, next(plan_iterator))
                for _ in range(min(download_workers, len(plans)))
            }
            while pending:
                try:
                    apply_event(events.get(timeout=0.1))
                except queue.Empty:
                    pass
                while True:
                    try:
                        apply_event(events.get_nowait())
                    except queue.Empty:
                        break
                finished = {future for future in pending if future.done()}
                for future in finished:
                    pending.remove(future)
                    try:
                        future.result()
                    except Exception as exc:
                        errors.append(exc)
                        stop.set()
                if not errors:
                    for _ in finished:
                        try:
                            plan = next(plan_iterator)
                        except StopIteration:
                            break
                        pending.add(executor.submit(run_plan, plan))
            while True:
                try:
                    apply_event(events.get_nowait())
                except queue.Empty:
                    break
    if errors:
        first = errors[0]
        if isinstance(first, ZenodoReleaseError):
            raise first
        raise ZenodoReleaseError(
            f"Parallel archive download failed: {first}"
        ) from first

    for directory in sorted(
        (path for path in temporary_root.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    ):
        directory.rmdir()
    temporary_root.rmdir()
    session["completed_segments"] = sorted(completed_segments)
    session["completed_files"] = sorted(completed)
    verified_bytes = sum(
        item["size_bytes"] for record in RECORDS for item in remote_items[record]
    )
    return expected_paths, verified_bytes


def download_and_reconstruct(
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    *,
    download_root: Path,
    reconstructed_root: Path,
    download_workers: int = DEFAULT_DOWNLOAD_WORKERS,
    download_client_factory: Callable[[], Any] | None = None,
) -> None:
    """Download authenticated draft bytes, SHA-verify, then reconstruct exactly."""
    validate_state(state)
    if (
        not _plain_int(download_workers)
        or not 1 <= download_workers <= MAX_DOWNLOAD_WORKERS
    ):
        raise ZenodoReleaseError(
            f"Download workers must be between 1 and {MAX_DOWNLOAD_WORKERS}"
        )
    if any(not state["records"][name].get("remote_verified") for name in RECORDS):
        raise ZenodoReleaseError(
            "Every remote record must verify before roundtrip download"
        )
    target = download_root.expanduser().absolute()
    output = reconstructed_root.expanduser().absolute()
    session = state["remote_download"]
    recorded_root = session.get("root")
    if recorded_root is None:
        if target.exists() or target.is_symlink():
            raise ZenodoReleaseError(
                f"First roundtrip download requires an absent directory: {target}"
            )
        if not target.parent.is_dir() or target.parent.is_symlink():
            raise ZenodoReleaseError(
                f"Download parent must be a regular directory: {target.parent}"
            )
        session["root"] = str(target)
        session["status"] = "IN_PROGRESS"
        session["completed_files"] = []
        _touch(state)
        persist(state)
        target.mkdir(mode=0o700)
    elif Path(recorded_root) != target:
        raise ZenodoReleaseError(
            f"Roundtrip download is already bound to a different path: {recorded_root}"
        )
    elif (
        not target.exists()
        and session.get("status") == "IN_PROGRESS"
        and target.parent.is_dir()
        and not target.parent.is_symlink()
    ):
        target.mkdir(mode=0o700)
    elif not target.is_dir() or target.is_symlink():
        raise ZenodoReleaseError(
            f"Bound download directory is missing or linked: {target}"
        )

    completed = set(session.get("completed_files", []))
    transport = load_segmented_transport(state)
    expected_paths: set[str] = set()
    depositions: dict[str, Mapping[str, Any]] = {}
    remotes: dict[str, Mapping[str, Mapping[str, Any]]] = {}
    for record in RECORDS:
        row = state["records"][record]
        deposition = client.get_deposition(row["deposition_id"])
        if deposition.get("submitted") is not False:
            raise ZenodoReleaseError(
                f"Roundtrip verification expected an unpublished draft: {record}"
            )
        remote = normal_remote_files(deposition)
        remote_items = expected_remote_items(state, record)
        expected_names = {item["name"] for item in remote_items}
        if set(remote) != expected_names:
            raise ZenodoReleaseError(
                f"Remote inventory changed before download: {record}"
            )
        for item in remote_items:
            if not _remote_matches(remote[item["name"]], item):
                raise ZenodoReleaseError(
                    f"Remote size/MD5 changed before download: "
                    f"{record}/{item['name']}"
                )
        depositions[record] = deposition
        remotes[record] = remote
        if transport is not None:
            continue
        for item in remote_items:
            relative = _download_relative(item)
            expected_paths.add(relative)
            remote_item = remote[item["name"]]
            destination = target.joinpath(*PurePosixPath(relative).parts)
            if relative in completed:
                if (
                    destination.is_symlink()
                    or not destination.is_file()
                    or destination.stat().st_size != item["size_bytes"]
                    or sha256_file(destination) != item["sha256"]
                ):
                    raise ZenodoReleaseError(
                        f"Previously verified roundtrip file changed: {relative}"
                    )
                continue
            client.download_file(
                _download_url(deposition, remote_item, item["name"]),
                destination,
                expected_size=item["size_bytes"],
                expected_sha256=item["sha256"],
            )
            if (
                destination.stat().st_size != item["size_bytes"]
                or sha256_file(destination) != item["sha256"]
            ):
                raise ZenodoReleaseError(
                    f"Downloaded size/SHA-256 differs after client return: {relative}"
                )
            completed.add(relative)
            session["completed_files"] = sorted(completed)
            _touch(state)
            persist(state)
    if transport is not None:
        expected_paths, verified_bytes = _download_segmented_archives(
            state,
            client,
            persist,
            target=target,
            transport=transport,
            depositions=depositions,
            remotes=remotes,
            completed=completed,
            download_workers=download_workers,
            download_client_factory=download_client_factory,
        )
    else:
        verified_bytes = sum(
            item["size_bytes"]
            for record in RECORDS
            for item in state["records"][record]["files"]
        )
    _closed_download_inventory(target, expected_paths)
    downloaded_manifest = target / "manifest.json"
    if (
        downloaded_manifest.stat().st_size != state["manifest"]["size_bytes"]
        or sha256_file(downloaded_manifest) != state["manifest"]["sha256"]
    ):
        raise ZenodoReleaseError("Remotely downloaded manifest binding differs")

    recorded_output = session.get("reconstructed_root")
    if recorded_output is None:
        if output.exists() or output.is_symlink():
            raise ZenodoReleaseError(
                f"Roundtrip reconstruction requires an absent output: {output}"
            )
        reconstruction.reconstruct(
            downloaded_manifest,
            output,
            archive_root=target,
        )
        session["reconstructed_root"] = str(output)
    elif Path(recorded_output) != output or not output.is_dir() or output.is_symlink():
        raise ZenodoReleaseError(
            "Roundtrip reconstructed tree differs from the state-bound output"
        )
    session["status"] = "PASS"
    session["manifest_sha256"] = state["manifest"]["sha256"]
    session["tree_sha256"] = state["manifest"]["tree_sha256"]
    if transport is not None:
        session["transport_manifest_sha256"] = state["segmented_transport"]["sha256"]
    session["verified_bytes"] = verified_bytes
    _touch(state)
    persist(state)


def _default_runner(command: Sequence[str], *, cwd: Path, log: BinaryIO) -> int:
    environment = dict(os.environ)
    for name in (TOKEN_ENV, "ZENODO_SANDBOX_TOKEN", "GH_TOKEN", "GITHUB_TOKEN"):
        environment.pop(name, None)
    completed = subprocess.run(
        list(command),
        cwd=cwd,
        env=environment,
        stdin=subprocess.DEVNULL,
        stdout=log,
        stderr=subprocess.STDOUT,
        check=False,
    )
    return completed.returncode


def reproduce_roundtrip_figures(
    state: MutableMapping[str, Any],
    persist: Persist,
    *,
    log_path: Path,
    execute: bool,
    runner: Callable[..., int] = _default_runner,
) -> None:
    if not execute:
        raise ZenodoReleaseError(
            "Figure reproduction refused without --execute-reproduction"
        )
    session = state["remote_download"]
    binding = _segmented_binding(state)
    if (
        session.get("status") != "PASS"
        or session.get("manifest_sha256") != state["manifest"]["sha256"]
        or session.get("tree_sha256") != state["manifest"]["tree_sha256"]
        or (
            binding is not None
            and session.get("transport_manifest_sha256") != binding["sha256"]
        )
    ):
        raise ZenodoReleaseError(
            "Verified remote download/reconstruction gate has not passed"
        )
    root = Path(session["reconstructed_root"]).resolve()
    script = root / "reproduce_all_figures.sh"
    if (
        script.is_symlink()
        or not script.is_file()
        or not os.access(script, os.X_OK)
        or not script.resolve().is_relative_to(root)
    ):
        raise ZenodoReleaseError(f"Figure reproduction entry point is unsafe: {script}")
    target = log_path.expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise ZenodoReleaseError(f"Figure log output must be absent: {target}")
    if not target.parent.is_dir() or target.parent.is_symlink():
        raise ZenodoReleaseError(f"Figure log parent is invalid: {target.parent}")
    command = [str(script), "--output", str(root), "--force"]
    state["reproduction"] = {
        "status": "RUNNING",
        "log_path": str(target),
        "log_sha256": None,
        "command": command,
        "remote_manifest_sha256": state["manifest"]["sha256"],
        "remote_tree_sha256": state["manifest"]["tree_sha256"],
        "remote_transport_manifest_sha256": (
            binding["sha256"] if binding is not None else None
        ),
    }
    _touch(state)
    persist(state)
    with target.open("xb") as log:
        returncode = runner(command, cwd=root, log=log)
        log.flush()
        os.fsync(log.fileno())
    os.chmod(target, 0o644)
    state["reproduction"]["returncode"] = returncode
    state["reproduction"]["log_sha256"] = sha256_file(target)
    state["reproduction"]["status"] = "PASS" if returncode == 0 else "FAIL"
    _touch(state)
    persist(state)
    if returncode != 0:
        raise ZenodoReleaseError(
            f"All-figure reproduction failed with exit code {returncode}; see {target}"
        )


def _roundtrip_gate(state: Mapping[str, Any]) -> None:
    download = state.get("remote_download", {})
    reproduction = state.get("reproduction", {})
    log_value = reproduction.get("log_path")
    log_path = (
        Path(log_value).expanduser().absolute() if isinstance(log_value, str) else None
    )
    binding = _segmented_binding(state)
    if (
        download.get("status") != "PASS"
        or download.get("manifest_sha256") != state["manifest"]["sha256"]
        or download.get("tree_sha256") != state["manifest"]["tree_sha256"]
        or reproduction.get("status") != "PASS"
        or reproduction.get("returncode") != 0
        or reproduction.get("remote_manifest_sha256") != state["manifest"]["sha256"]
        or reproduction.get("remote_tree_sha256") != state["manifest"]["tree_sha256"]
        or (
            binding is not None
            and (
                download.get("transport_manifest_sha256") != binding["sha256"]
                or reproduction.get("remote_transport_manifest_sha256")
                != binding["sha256"]
            )
        )
        or not SHA256_RE.fullmatch(str(reproduction.get("log_sha256", "")))
        or log_path is None
        or log_path.is_symlink()
        or not log_path.is_file()
        or sha256_file(log_path) != reproduction.get("log_sha256")
    ):
        raise ZenodoReleaseError(
            "Publish blocked: remote downloads and all-figure reproduction are not PASS"
        )


def _is_published(deposition: Mapping[str, Any]) -> bool:
    return deposition.get("submitted") is True and deposition.get("state") == "done"


def _public_url(deposition: Mapping[str, Any]) -> str:
    value = deposition.get("record_url")
    if isinstance(value, str):
        return _safe_public_url(value)
    links = deposition.get("links", {})
    if isinstance(links, Mapping):
        for key in ("html", "record"):
            value = links.get(key)
            if isinstance(value, str) and "/deposit/" not in value:
                return _safe_public_url(value)
    value = deposition.get("doi_url")
    if isinstance(value, str):
        return _safe_public_url(value)
    raise ZenodoReleaseError("Published deposition has no public record URL")


def publish_all(
    state: MutableMapping[str, Any],
    client: Any,
    persist: Persist,
    *,
    execute: bool,
    production: bool,
    confirmation: str,
    expected_state_sha256: str,
) -> None:
    """Publish only after one full preflight; persist after each non-atomic action."""
    require_mutation_gate("publish", execute=execute, production=production)
    if confirmation != PUBLISH_CONFIRMATION:
        raise ZenodoReleaseError(
            f"Publish requires --confirm-irreversible {PUBLISH_CONFIRMATION}"
        )
    observed_state_sha = state_sha256(state)
    if expected_state_sha256 != observed_state_sha:
        raise ZenodoReleaseError(
            "Publish state SHA-256 differs; run status and review the current state"
        )
    validate_state(state)
    _roundtrip_gate(state)
    preflight: dict[str, Mapping[str, Any]] = {}
    for record in RECORDS:
        row = state["records"][record]
        _quota_ready(record, row)
        if not row.get("remote_verified"):
            raise ZenodoReleaseError(f"Remote verification is not PASS: {record}")
        deposition = client.get_deposition(row["deposition_id"])
        verify_metadata(
            deposition, effective_metadata(state, record), require_related=True
        )
        _verify_remote_inventory(state, record, deposition)
        if deposition.get("submitted") not in {False, True}:
            raise ZenodoReleaseError(f"Zenodo publication state is invalid: {record}")
        if deposition.get("submitted") is True and not _is_published(deposition):
            raise ZenodoReleaseError(
                f"Zenodo publication is incomplete/error: {record}"
            )
        if row.get("reserved_doi") != _reserved_doi(deposition):
            raise ZenodoReleaseError(f"Reserved DOI changed: {record}")
        preflight[record] = deposition

    # Zenodo has no multi-record transaction. A failure after one POST leaves a
    # partial public set; state is saved immediately so rerunning reconciles and
    # publishes only the remaining preflighted records.
    for record in RECORDS:
        row = state["records"][record]
        deposition = preflight[record]
        if not _is_published(deposition):
            try:
                deposition = client.publish(row["deposition_id"])
            except ZenodoReleaseError:
                # The POST may have succeeded before a transport failure. Reconcile
                # once rather than issuing a second irreversible request.
                deposition = client.get_deposition(row["deposition_id"])
            if not _is_published(deposition):
                # An ambiguous response is reconciled once; the POST is never retried.
                deposition = client.get_deposition(row["deposition_id"])
            if not _is_published(deposition):
                raise ZenodoReleaseError(
                    f"Zenodo did not confirm publication for {record}; reconcile manually"
                )
        if _reserved_doi(deposition) != row["reserved_doi"]:
            raise ZenodoReleaseError(f"Published DOI differs for {record}")
        row["published"] = True
        row["public_url"] = _public_url(deposition)
        row["record_id"] = int(deposition.get("record_id", row["deposition_id"]))
        _touch(state)
        persist(state)


def _public_file_url(record_id: int, filename: str) -> str:
    if record_id <= 0 or PurePosixPath(filename).name != filename:
        raise ZenodoReleaseError("Cannot construct a stable public Zenodo file URL")
    return _safe_public_url(
        f"https://zenodo.org/records/{record_id}/files/"
        + urllib.parse.quote(filename, safe="")
        + "?download=1"
    )


def build_public_url_map(
    state: MutableMapping[str, Any], client: Any
) -> dict[str, Any]:
    validate_state(state)
    if any(not state["records"][name].get("published") for name in RECORDS):
        raise ZenodoReleaseError("All four records must be published first")
    transport = load_segmented_transport(state)
    urls: dict[str, str] = {}
    segmented_urls: dict[str, Any] = {}
    for record in RECORDS:
        row = state["records"][record]
        public = client.get_public_record(row["record_id"])
        remote = normal_remote_files(public)
        expected = {item["name"]: item for item in expected_remote_items(state, record)}
        if set(remote) != set(expected):
            raise ZenodoReleaseError(f"Public file inventory differs: {record}")
        for filename, item in expected.items():
            if not _remote_matches(remote[filename], item):
                raise ZenodoReleaseError(
                    f"Public size/MD5 differs: {record}/{filename}"
                )
            if transport is not None or filename == "manifest.json":
                continue
            relative = item["relative_local_path"]
            urls[relative] = _public_file_url(row["record_id"], filename)
        if transport is not None:
            for archive in transport["records"][record]["archives"]:
                relative = archive["relative_local_path"]
                segmented_urls[relative] = {
                    "size_bytes": archive["size_bytes"],
                    "sha256": archive["sha256"],
                    "segments": [
                        {
                            "name": segment["name"],
                            "size_bytes": segment["size_bytes"],
                            "sha256": segment["sha256"],
                            "url": _public_file_url(
                                state["records"][segment.get("storage_record", record)][
                                    "record_id"
                                ],
                                segment["name"],
                            ),
                        }
                        for segment in archive["segments"]
                    ],
                }
    if transport is not None:
        expected_urls = {
            item["relative_local_path"]
            for record in RECORDS
            for item in state["records"][record]["files"]
            if item["name"] != "manifest.json"
        }
        if set(segmented_urls) != expected_urls:
            raise ZenodoReleaseError(
                "Segmented public URL map does not match the archive manifest"
            )
        return {
            "schema": reconstruction.URL_MAP_SCHEMA_V2,
            "manifest": {
                "sha256": state["manifest"]["sha256"],
                "size_bytes": state["manifest"]["size_bytes"],
            },
            "archives": dict(sorted(segmented_urls.items())),
        }
    expected_urls = {
        item["relative_local_path"]
        for record in RECORDS
        for item in state["records"][record]["files"]
        if item["name"] != "manifest.json"
    }
    if set(urls) != expected_urls:
        raise ZenodoReleaseError("Public URL map does not match the archive manifest")
    return {
        "schema": reconstruction.URL_MAP_SCHEMA,
        "archives": dict(sorted(urls.items())),
    }


def write_new_json(path: Path, document: Mapping[str, Any], *, mode: int) -> None:
    target = path.expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise ZenodoReleaseError(f"JSON output must be absent: {target}")
    if not target.parent.is_dir() or target.parent.is_symlink():
        raise ZenodoReleaseError(f"JSON output parent is invalid: {target.parent}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    descriptor = os.open(target, flags, mode)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(canonical_json(document))
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(target, mode)
    except Exception:
        try:
            target.unlink()
        except FileNotFoundError:
            pass
        raise


def status_summary(state: Mapping[str, Any]) -> dict[str, Any]:
    transport = _segmented_binding(state)
    return {
        "status": "READY",
        "environment": state["environment"],
        "api_base": state["api_base"],
        "release_id": state["release_id"],
        "state_sha256": state_sha256(state),
        "manifest_sha256": state["manifest"]["sha256"],
        "tree_sha256": state["manifest"]["tree_sha256"],
        "segmented_transport": (
            {
                "schema": transport["schema"],
                "sha256": transport["sha256"],
                "size_bytes": transport["size_bytes"],
            }
            if transport is not None
            else None
        ),
        "remote_download": state["remote_download"]["status"],
        "reproduction": state["reproduction"]["status"],
        "records": {
            name: {
                "deposition_id": state["records"][name]["deposition_id"],
                "reserved_doi": state["records"][name]["reserved_doi"],
                "license": state["records"][name]["license"],
                "metadata_synced": state["records"][name]["metadata_synced"],
                "remote_verified": state["records"][name]["remote_verified"],
                "published": state["records"][name]["published"],
                "public_url": state["records"][name]["public_url"],
                "files_verified": sum(
                    (
                        bool(state["records"][name]["remote_verified"])
                        if transport is not None
                        else bool(row["remote_verified"])
                    )
                    for row in (
                        expected_remote_items(state, name)
                        if transport is not None
                        else state["records"][name]["files"]
                    )
                ),
                "file_count": len(expected_remote_items(state, name)),
                "quota": state["records"][name]["quota"],
            }
            for name in RECORDS
        },
        "credentials_embedded": False,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--state", type=Path, required=True)
    result.add_argument(
        "--token-file",
        type=Path,
        help="mode-0600 file containing the production token; never pass tokens inline",
    )
    result.add_argument("--timeout-seconds", type=float, default=300.0)
    result.add_argument("--retries", type=int, default=4)
    sub = result.add_subparsers(dest="command", required=True)

    init = sub.add_parser("init", help="locally verify and bind an approved request")
    init.add_argument("--request", type=Path, required=True)
    init.add_argument("--archive-root", type=Path, required=True)
    init.add_argument(
        "--adopt-draft",
        action="append",
        default=[],
        metavar="RECORD=DEPOSITION_ID",
        help="seed one reviewed empty production draft by exact record and ID",
    )
    init.add_argument("--plan", type=Path, default=HERE / "release_plan.json")
    init.add_argument(
        "--author-audit", type=Path, default=PAPER / control.EXPECTED_AUTHOR_AUDIT
    )
    init.add_argument(
        "--contract-script",
        type=Path,
        default=PAPER / "scripts/utilities/01_validate_bundle.py",
    )
    sub.add_parser("status", help="show credential-free phase status and state digest")
    sub.add_parser(
        "effective-metadata",
        help="print DOI-cross-linked metadata after drafts have reserved all DOIs",
    )

    drafts = sub.add_parser(
        "drafts", help="create/adopt drafts and synchronize metadata"
    )
    drafts.add_argument("--execute-drafts", action="store_true")
    drafts.add_argument("--production", action="store_true")

    quota = sub.add_parser(
        "confirm-quota",
        help="bind a manually reviewed Manage-storage allocation to one draft ID",
    )
    quota.add_argument("--record", choices=RECORDS, required=True)
    quota.add_argument("--deposition-id", type=int, required=True)
    quota.add_argument("--allocated-record-bytes", type=int, required=True)
    quota.add_argument("--confirm-quota-allocation", action="store_true")

    upload = sub.add_parser(
        "upload", help="stream missing files and reconcile exact MD5/size"
    )
    upload.add_argument("--execute-uploads", action="store_true")
    upload.add_argument("--production", action="store_true")
    sub.add_parser(
        "verify-remote", help="read-only exact draft metadata/file reconciliation"
    )

    download = sub.add_parser(
        "download-roundtrip",
        help="authenticated resumable SHA-256 download and exact reconstruction",
    )
    download.add_argument("--download-root", type=Path, required=True)
    download.add_argument("--reconstructed-root", type=Path, required=True)
    download.add_argument(
        "--download-workers",
        type=int,
        choices=range(1, MAX_DOWNLOAD_WORKERS + 1),
        default=DEFAULT_DOWNLOAD_WORKERS,
        help="parallel canonical-archive downloads (default: 4; maximum: 8)",
    )

    reproduce = sub.add_parser(
        "reproduce-roundtrip",
        help="run all figures from the remotely reconstructed tree",
    )
    reproduce.add_argument("--log", type=Path, required=True)
    reproduce.add_argument("--execute-reproduction", action="store_true")

    publish = sub.add_parser(
        "publish", help="irreversibly publish all preflighted records"
    )
    publish.add_argument("--execute-publish", action="store_true")
    publish.add_argument("--production", action="store_true")
    publish.add_argument("--confirm-irreversible", required=True)
    publish.add_argument("--expected-state-sha256", required=True)

    urls = sub.add_parser(
        "public-url-map",
        help="verify public records and write reconstruct_paper-compatible URLs",
    )
    urls.add_argument("--output", type=Path, required=True)
    return result


def _needs_client(command: str) -> bool:
    return command in {
        "drafts",
        "upload",
        "verify-remote",
        "download-roundtrip",
        "publish",
        "public-url-map",
    }


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "init":
            target = args.state.expanduser().absolute()
            if target.exists() or target.is_symlink():
                raise ZenodoReleaseError(f"Initial state path must be absent: {target}")
            state = initialize_state(
                args.request,
                args.archive_root,
                plan_path=args.plan,
                audit_path=args.author_audit,
                contract_script=args.contract_script,
                adopt_drafts=parse_adopted_drafts(args.adopt_draft),
            )
            write_state(target, state)
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
            print("[LOCAL-ONLY] no Zenodo request was sent")
            return 0

        state = load_state(args.state)
        persist = lambda document: write_state(args.state, document)
        client = None
        if _needs_client(args.command):
            token = load_token(args.token_file)
            try:
                client = ZenodoClient(
                    token,
                    timeout_seconds=args.timeout_seconds,
                    retries=args.retries,
                )
            finally:
                token = ""

        if args.command == "status":
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "effective-metadata":
            document = {record: effective_metadata(state, record) for record in RECORDS}
            print(json.dumps(document, ensure_ascii=False, indent=2))
        elif args.command == "drafts":
            create_and_link_drafts(
                state,
                client,
                persist,
                execute=args.execute_drafts,
                production=args.production,
            )
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "confirm-quota":
            confirm_quota(
                state,
                args.record,
                deposition_id=args.deposition_id,
                allocated_record_bytes=args.allocated_record_bytes,
                confirmed=args.confirm_quota_allocation,
            )
            persist(state)
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "upload":
            upload_all(
                state,
                client,
                persist,
                execute=args.execute_uploads,
                production=args.production,
            )
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "verify-remote":
            verify_remote(state, client, persist)
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "download-roundtrip":
            download_and_reconstruct(
                state,
                client,
                persist,
                download_root=args.download_root,
                reconstructed_root=args.reconstructed_root,
                download_workers=args.download_workers,
            )
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "reproduce-roundtrip":
            reproduce_roundtrip_figures(
                state,
                persist,
                log_path=args.log,
                execute=args.execute_reproduction,
            )
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "publish":
            publish_all(
                state,
                client,
                persist,
                execute=args.execute_publish,
                production=args.production,
                confirmation=args.confirm_irreversible,
                expected_state_sha256=args.expected_state_sha256,
            )
            print(json.dumps(status_summary(state), ensure_ascii=False, indent=2))
        elif args.command == "public-url-map":
            document = build_public_url_map(state, client)
            write_new_json(args.output, document, mode=0o644)
            print(f"[PASS] public URL map: {args.output.expanduser().absolute()}")
        else:
            raise AssertionError(args.command)
        return 0
    except (
        ZenodoReleaseError,
        control.ReleaseControlError,
        archives.ArchiveError,
        reconstruction.ReconstructionError,
    ) as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
