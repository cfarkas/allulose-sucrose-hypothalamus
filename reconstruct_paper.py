#!/usr/bin/env python3
"""Download verified public Zenodo shards and reconstruct one exact Paper tree.

The local manifest is the authority for every path, size, SHA-256 digest and
POSIX mode.  Archives can come from a local archive-set directory or from an
exact HTTPS URL map.  URL-map v2 can describe ordered transport segments for
each canonical ZIP; segments are verified individually and concatenated back
to the manifest-bound ZIP before extraction.  Interrupted v2 downloads retain
only manifest/URL-map-bound, checksum-verified progress for a safe rerun.  The
destination is created transactionally and is never merged with or written
over an existing tree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import re
import shutil
import stat
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from contextlib import contextmanager
from pathlib import Path, PurePosixPath
from typing import Any, BinaryIO, Callable, Iterator, Mapping, Sequence

import build_archives as archive


URL_MAP_SCHEMA_V1 = "apotome-paper-archive-url-map-v1"
URL_MAP_SCHEMA_V2 = "apotome-paper-archive-url-map-v2"
# Keep the historical name as the v1 identifier for callers that still emit
# the original relative-path-to-URL mapping.  New producers should explicitly
# use URL_MAP_SCHEMA_V2.
URL_MAP_SCHEMA = URL_MAP_SCHEMA_V1
MAX_MANIFEST_BYTES = 1_000_000_000
MAX_URL_MAP_BYTES = 10_000_000
MAX_RESUME_RECEIPT_BYTES = 100_000
USER_AGENT = "Apotome-Paper-Reconstructor/1"
SEGMENT_DOWNLOAD_ATTEMPTS = 5
SEGMENT_RETRY_INITIAL_SECONDS = 2.0
SEGMENT_RETRY_MAX_SECONDS = 30.0
RESUME_SCHEMA = "apotome-paper-segmented-download-resume-v1"
RESUME_RECEIPT_NAME = ".apotome-download-resume.json"
PLACEHOLDER_RE = re.compile(
    r"(?i)(<[^>]+>|YOUR[_-]|REPLACE[_-]?ME|CHANGEME|TODO|TBD|example\.com)"
)
SECRET_QUERY_KEYS = {
    "access_token",
    "api_key",
    "apikey",
    "auth",
    "authorization",
    "key",
    "signature",
    "sig",
    "token",
}
Opener = Callable[..., Any]
LOG = logging.getLogger(__name__)
LOG.addHandler(logging.NullHandler())


def progress(label: str, count: int, total: int, last: float) -> float:
    """Emit periodic byte progress without changing verification."""
    now = time.monotonic()
    if LOG.isEnabledFor(logging.INFO) and now - last >= 10:
        LOG.info("%s: %s / %s bytes (%.1f%%)", label, f"{count:,}",
                 f"{total:,}", 100 * count / max(total, 1))
        return now
    return last


class ReconstructionError(archive.ArchiveError):
    """A manifest, download, archive or output violates the release contract."""


def _unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ReconstructionError(f"Duplicate JSON key: {key!r}")
        result[key] = value
    return result


def read_json_object(path: Path, *, label: str, max_bytes: int) -> dict[str, Any]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink():
        raise ReconstructionError(f"{label} may not be a symlink: {lexical}")
    try:
        observed = lexical.lstat()
    except FileNotFoundError as exc:
        raise ReconstructionError(f"{label} is missing: {lexical}") from exc
    if not stat.S_ISREG(observed.st_mode):
        raise ReconstructionError(f"{label} is not a regular file: {lexical}")
    if observed.st_size > max_bytes:
        raise ReconstructionError(f"{label} is unexpectedly large: {observed.st_size}")
    try:
        document = json.loads(
            lexical.read_text(encoding="utf-8"),
            object_pairs_hook=_unique_object,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReconstructionError(
            f"Cannot read valid {label}: {lexical}: {exc}"
        ) from exc
    if not isinstance(document, dict):
        raise ReconstructionError(f"{label} root must be a JSON object")
    return document


def _plain_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _manifest_shards(document: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    return [
        shard
        for record in archive.RECORDS
        for shard in document["records"][record]["shards"]
    ]


def validate_release_manifest(document: Mapping[str, Any]) -> None:
    """Apply builder validation plus publication and reconstruction invariants."""
    if not isinstance(document, dict):
        raise ReconstructionError("Manifest root must be a JSON object")
    try:
        archive.validate_manifest(document)
    except ReconstructionError:
        raise
    except Exception as exc:
        raise ReconstructionError(f"Invalid release manifest: {exc}") from exc

    if (
        document.get("paper_root_name") != "Paper"
        or document.get("archive_member_prefix") != "Paper/"
    ):
        raise ReconstructionError("Manifest Paper root/member prefix differs")
    compression = document.get("compression")
    if (
        not isinstance(compression, dict)
        or compression.get("method") != "ZIP_DEFLATED"
        or not _plain_int(compression.get("level"))
        or not 1 <= compression["level"] <= 9
    ):
        raise ReconstructionError("Manifest compression contract differs")

    contract = document.get("canonical_exclusion_contract")
    if not isinstance(contract, dict):
        raise ReconstructionError("Canonical exclusion contract is missing")
    if (
        contract.get("script") != "scripts/utilities/01_validate_bundle.py"
        or contract.get("local_only_directories") != sorted(archive.LOCAL_ONLY)
        or contract.get("local_only_root_prefixes")
        != list(archive.LOCAL_ONLY_ROOT_PREFIXES)
        or contract.get("local_only_root_files") != sorted(archive.LOCAL_ONLY_FILES)
        or contract.get("local_only_file_suffixes")
        != list(archive.LOCAL_ONLY_FILE_SUFFIXES)
        or contract.get("public_excluded_relative_prefixes")
        != list(archive.PUBLIC_EXCLUDED_RELATIVE_PREFIXES)
        or contract.get("redistribution_excluded_basenames")
        != sorted(archive.REDISTRIBUTION_EXCLUDED_BASENAMES)
        or not archive.SHA_RE.fullmatch(str(contract.get("script_sha256", "")))
    ):
        raise ReconstructionError("Canonical exclusion contract differs")

    limits = document["limits"]
    numeric_limit_keys = (
        "uncompressed_payload_per_shard",
        "compressed_bytes_per_shard",
        "calculated_additional_quota_bytes",
    )
    if any(
        not _plain_int(limits.get(key)) or limits[key] < 0 for key in numeric_limit_keys
    ):
        raise ReconstructionError("Manifest numeric limits are invalid")
    fixed_limits = {
        "zenodo_individual_file_bytes": archive.ZENODO_INDIVIDUAL_FILE_LIMIT,
        "zenodo_files_per_record": archive.ZENODO_MAX_FILES_PER_RECORD,
        "zenodo_quota_raised_record_bytes": archive.ZENODO_MAX_RECORD_BYTES,
        "zenodo_account_additional_allowance_bytes": (
            archive.ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE
        ),
    }
    if any(limits.get(key) != value for key, value in fixed_limits.items()):
        raise ReconstructionError("Manifest Zenodo limits differ")

    compressed_totals: dict[str, int] = {}
    for record in archive.RECORDS:
        row = document["records"][record]
        numeric_record_keys = (
            "shard_count",
            "zenodo_uploaded_shard_files",
            "tree_file_count",
            "uncompressed_bytes",
            "compressed_archive_bytes",
        )
        if any(
            not _plain_int(row.get(key)) or row[key] < 0 for key in numeric_record_keys
        ):
            raise ReconstructionError(f"Invalid numeric totals for {record}")
        shards = row["shards"]
        if not shards:
            raise ReconstructionError(f"Logical record has no archive shards: {record}")
        expected_names = [
            f"{archive.RECORD_PREFIX[record]}.part{index:03d}.zip"
            for index in range(1, len(shards) + 1)
        ]
        if [shard["name"] for shard in shards] != expected_names:
            raise ReconstructionError(f"Shard names/order differ for {record}")
        for shard in shards:
            numeric_shard_keys = (
                "member_count",
                "uncompressed_bytes",
                "compressed_bytes",
            )
            if any(
                not _plain_int(shard.get(key)) or shard[key] < 0
                for key in numeric_shard_keys
            ):
                raise ReconstructionError(f"Invalid numeric values for {shard['name']}")
            if shard["compressed_bytes"] > archive.ZENODO_INDIVIDUAL_FILE_LIMIT:
                raise ReconstructionError(
                    f"Shard exceeds Zenodo individual-file limit: {shard['name']}"
                )
        total = row["compressed_archive_bytes"]
        if not _plain_int(total) or total < 0:
            raise ReconstructionError(f"Invalid compressed total for {record}")
        compressed_totals[record] = total
    if compressed_totals["software_core"] >= archive.ZENODO_DEFAULT_RECORD_BYTES:
        raise ReconstructionError("software/core record is not below 50 GB")
    if any(
        total > archive.ZENODO_MAX_RECORD_BYTES for total in compressed_totals.values()
    ):
        raise ReconstructionError("A logical record exceeds the 200 GB ceiling")
    additional = sum(
        max(0, total - archive.ZENODO_DEFAULT_RECORD_BYTES)
        for total in compressed_totals.values()
    )
    if (
        additional > archive.ZENODO_ACCOUNT_ADDITIONAL_ALLOWANCE
        or limits.get("calculated_additional_quota_bytes") != additional
    ):
        raise ReconstructionError("Account-level additional quota calculation differs")

    tree = document["tree"]
    directories = tree["directories"]
    files = tree["files"]
    if [row["path"] for row in directories] != sorted(
        row["path"] for row in directories
    ):
        raise ReconstructionError("Directory manifest is not path-sorted")
    if [row["path"] for row in files] != sorted(row["path"] for row in files):
        raise ReconstructionError("File manifest is not path-sorted")

    directory_paths = {row["path"] for row in directories}
    file_paths = {row["path"] for row in files}
    for path in directory_paths | file_paths:
        parent = PurePosixPath(path).parent
        while parent != PurePosixPath("."):
            if parent.as_posix() not in directory_paths:
                raise ReconstructionError(
                    f"Manifest omits parent directory {parent} for {path}"
                )
            parent = parent.parent
    direct_file_parents = {PurePosixPath(path).parent.as_posix() for path in file_paths}
    direct_directory_parents = {
        PurePosixPath(path).parent.as_posix() for path in directory_paths
    }
    for row in directories:
        expected_empty = (
            row["path"] not in direct_file_parents
            and row["path"] not in direct_directory_parents
        )
        if row["empty"] is not expected_empty:
            raise ReconstructionError(
                f"Empty-directory marker differs for {row['path']}"
            )

    contract_path = contract["script"]
    contract_rows = [row for row in files if row["path"] == contract_path]
    if (
        len(contract_rows) != 1
        or contract_rows[0]["sha256"] != contract["script_sha256"]
    ):
        raise ReconstructionError(
            "Canonical validator file is absent or differs from its contract hash"
        )


def load_manifest(path: Path) -> dict[str, Any]:
    document = read_json_object(
        path,
        label="manifest",
        max_bytes=MAX_MANIFEST_BYTES,
    )
    validate_release_manifest(document)
    return document


def validate_public_https_url(url: Any, *, reject_secret_query: bool = True) -> str:
    if not isinstance(url, str) or url != url.strip() or PLACEHOLDER_RE.search(url):
        raise ReconstructionError(
            f"URL is missing, placeholder-like, or malformed: {url!r}"
        )
    parsed = urllib.parse.urlsplit(url)
    try:
        port = parsed.port
    except ValueError as exc:
        raise ReconstructionError(f"URL port is invalid: {url!r}") from exc
    del port
    if (
        parsed.scheme.lower() != "https"
        or not parsed.hostname
        or parsed.username is not None
        or parsed.password is not None
        or parsed.fragment
    ):
        raise ReconstructionError(
            f"Only public, credential-free HTTPS URLs are allowed: {url!r}"
        )
    if reject_secret_query:
        query_keys = {key.lower() for key, _ in urllib.parse.parse_qsl(parsed.query)}
        if query_keys & SECRET_QUERY_KEYS:
            raise ReconstructionError(
                f"Public URL must not embed access credentials: {url!r}"
            )
    return url


def _manifest_binding(document: Mapping[str, Any]) -> dict[str, Any]:
    payload = archive.canonical_json(document)
    return {
        "size_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
    }


def _valid_sha256(value: Any) -> bool:
    return isinstance(value, str) and archive.SHA_RE.fullmatch(value) is not None


def _validate_segment_name(value: Any) -> str:
    if (
        not isinstance(value, str)
        or value != value.strip()
        or not value
        or "\\" in value
        or "\x00" in value
        or PurePosixPath(value).name != value
        or value in {".", ".."}
    ):
        raise ReconstructionError(f"Unsafe transport-segment name: {value!r}")
    return value


def _load_url_map_v1(
    url_document: Mapping[str, Any], document: Mapping[str, Any]
) -> dict[str, str]:
    if set(url_document) != {"schema", "archives"}:
        raise ReconstructionError("URL-map v1 must contain only schema and archives")
    urls = url_document["archives"]
    if not isinstance(urls, dict):
        raise ReconstructionError("URL-map archives must be an object")
    expected = {shard["relative_archive_path"] for shard in _manifest_shards(document)}
    if set(urls) != expected:
        missing = sorted(expected - set(urls))
        extra = sorted(set(urls) - expected)
        raise ReconstructionError(
            f"URL map must match manifest shards exactly; missing={missing}, extra={extra}"
        )
    return {
        relative: validate_public_https_url(url)
        for relative, url in sorted(urls.items())
    }


def _load_url_map_v2(
    url_document: Mapping[str, Any], document: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    if set(url_document) != {"schema", "manifest", "archives"}:
        raise ReconstructionError(
            "URL-map v2 must contain only schema, manifest, and archives"
        )

    binding = url_document["manifest"]
    if not isinstance(binding, dict) or set(binding) != {"size_bytes", "sha256"}:
        raise ReconstructionError(
            "URL-map v2 manifest binding must contain size_bytes and sha256"
        )
    expected_binding = _manifest_binding(document)
    if (
        not _plain_int(binding.get("size_bytes"))
        or binding["size_bytes"] <= 0
        or not _valid_sha256(binding.get("sha256"))
        or binding != expected_binding
    ):
        raise ReconstructionError(
            "URL-map v2 does not bind to the exact canonical release manifest"
        )

    rows = url_document["archives"]
    if not isinstance(rows, dict):
        raise ReconstructionError("URL-map v2 archives must be an object")
    shards = {
        shard["relative_archive_path"]: shard for shard in _manifest_shards(document)
    }
    if set(rows) != set(shards):
        missing = sorted(set(shards) - set(rows))
        extra = sorted(set(rows) - set(shards))
        raise ReconstructionError(
            f"URL map must match manifest shards exactly; missing={missing}, extra={extra}"
        )

    normalized: dict[str, dict[str, Any]] = {}
    all_urls: set[str] = set()
    for relative in sorted(rows):
        row = rows[relative]
        shard = shards[relative]
        if not isinstance(row, dict) or set(row) != {
            "size_bytes",
            "sha256",
            "segments",
        }:
            raise ReconstructionError(
                f"URL-map v2 archive entry has invalid fields: {relative}"
            )
        if (
            not _plain_int(row.get("size_bytes"))
            or row["size_bytes"] != shard["compressed_bytes"]
            or not _valid_sha256(row.get("sha256"))
            or row["sha256"] != shard["sha256"]
        ):
            raise ReconstructionError(
                f"URL-map v2 canonical ZIP binding differs: {relative}"
            )
        segments = row["segments"]
        if not isinstance(segments, list) or not segments:
            raise ReconstructionError(
                f"URL-map v2 archive requires ordered segments: {relative}"
            )
        normalized_segments: list[dict[str, Any]] = []
        names: set[str] = set()
        segment_total = 0
        for index, segment in enumerate(segments, start=1):
            if not isinstance(segment, dict) or set(segment) != {
                "name",
                "size_bytes",
                "sha256",
                "url",
            }:
                raise ReconstructionError(
                    f"URL-map v2 segment has invalid fields: {relative} #{index}"
                )
            name = _validate_segment_name(segment.get("name"))
            size = segment.get("size_bytes")
            digest = segment.get("sha256")
            if not _plain_int(size) or size <= 0:
                raise ReconstructionError(
                    f"Invalid transport-segment size: {relative} #{index}"
                )
            if not _valid_sha256(digest):
                raise ReconstructionError(
                    f"Invalid transport-segment SHA-256: {relative} #{index}"
                )
            if name in names:
                raise ReconstructionError(
                    f"Duplicate transport-segment name for {relative}: {name}"
                )
            url = validate_public_https_url(segment.get("url"))
            if url in all_urls:
                raise ReconstructionError(f"Transport-segment URL is reused: {url!r}")
            names.add(name)
            all_urls.add(url)
            segment_total += size
            normalized_segments.append(
                {
                    "name": name,
                    "size_bytes": size,
                    "sha256": digest,
                    "url": url,
                }
            )
        if segment_total != row["size_bytes"]:
            raise ReconstructionError(
                f"Transport-segment byte total differs for {relative}: "
                f"{segment_total} != {row['size_bytes']}"
            )
        normalized[relative] = {
            "size_bytes": row["size_bytes"],
            "sha256": row["sha256"],
            "segments": normalized_segments,
        }
    return normalized


def load_url_map(
    path: Path, document: Mapping[str, Any]
) -> dict[str, str | dict[str, Any]]:
    url_document = read_json_object(
        path,
        label="URL map",
        max_bytes=MAX_URL_MAP_BYTES,
    )
    schema = url_document.get("schema")
    if schema == URL_MAP_SCHEMA_V1:
        return _load_url_map_v1(url_document, document)
    if schema == URL_MAP_SCHEMA_V2:
        return _load_url_map_v2(url_document, document)
    raise ReconstructionError("URL-map schema differs")


def _resolved_new_path(path: Path, *, label: str) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.exists() or lexical.is_symlink():
        raise ReconstructionError(f"{label} must not already exist: {lexical}")
    return lexical.parent.resolve() / lexical.name


def _resolved_existing_directory(path: Path, *, label: str) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink():
        raise ReconstructionError(f"{label} may not be a symlink: {lexical}")
    resolved = lexical.resolve()
    if not resolved.is_dir():
        raise ReconstructionError(f"{label} is missing or not a directory: {resolved}")
    return resolved


def _paths_overlap(first: Path, second: Path) -> bool:
    return (
        first == second or first.is_relative_to(second) or second.is_relative_to(first)
    )


def _regular_file_matches(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    label: str,
) -> bool:
    """Return whether a regular, non-linked file has the exact bound bytes."""
    try:
        observed = path.lstat()
    except FileNotFoundError:
        return False
    if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
        raise ReconstructionError(f"{label} is linked or not a regular file: {path}")
    if observed.st_size != expected_size:
        return False
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise ReconstructionError(
            f"Cannot open {label} safely: {path}: {exc}"
        ) from exc
    digest = hashlib.sha256()
    count = 0
    with os.fdopen(descriptor, "rb") as handle:
        observed_open = os.fstat(handle.fileno())
        if not stat.S_ISREG(observed_open.st_mode):
            raise ReconstructionError(f"{label} is not a regular file: {path}")
        if observed_open.st_size != expected_size:
            return False
        for block in iter(lambda: handle.read(archive.CHUNK), b""):
            digest.update(block)
            count += len(block)
    return count == expected_size and digest.hexdigest() == expected_sha256


def _require_exact_file(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    label: str,
) -> None:
    if not _regular_file_matches(
        path,
        expected_size=expected_size,
        expected_sha256=expected_sha256,
        label=label,
    ):
        raise ReconstructionError(f"{label} size/SHA-256 differs: {path}")


def _archive_path(root: Path, relative: str) -> Path:
    safe = archive.safe_relative(relative)
    candidate = root.joinpath(*PurePosixPath(safe).parts)
    current = root
    for part in PurePosixPath(safe).parts:
        current /= part
        try:
            observed = current.lstat()
        except FileNotFoundError as exc:
            raise ReconstructionError(f"Archive is missing: {candidate}") from exc
        if stat.S_ISLNK(observed.st_mode):
            raise ReconstructionError(f"Archive path contains a symlink: {current}")
    resolved = candidate.resolve()
    if not resolved.is_relative_to(root):
        raise ReconstructionError(f"Archive escapes its root: {relative}")
    if not resolved.is_file():
        raise ReconstructionError(f"Archive is not a regular file: {resolved}")
    return resolved


@contextmanager
def verified_archive_handle(path: Path, shard: Mapping[str, Any]) -> Iterator[BinaryIO]:
    LOG.info("VERIFY ZIP %s (%s bytes), expected SHA-256 %s", path,
             f"{shard['compressed_bytes']:,}", shard["sha256"])
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        raise ReconstructionError(f"Cannot open archive safely: {path}: {exc}") from exc
    handle = os.fdopen(descriptor, "rb")
    try:
        observed = os.fstat(handle.fileno())
        if not stat.S_ISREG(observed.st_mode):
            raise ReconstructionError(f"Archive is not a regular file: {path}")
        if observed.st_size != shard["compressed_bytes"]:
            raise ReconstructionError(
                f"Archive byte count differs: {path}: "
                f"{observed.st_size} != {shard['compressed_bytes']}"
            )
        digest = hashlib.sha256()
        count = 0
        last = time.monotonic()
        for block in iter(lambda: handle.read(archive.CHUNK), b""):
            digest.update(block)
            count += len(block)
            last = progress(f"HASH {path.name}", count, shard["compressed_bytes"], last)
        if count != shard["compressed_bytes"] or digest.hexdigest() != shard["sha256"]:
            raise ReconstructionError(f"Archive SHA-256 differs: {path}")
        handle.seek(0)
        LOG.info("VERIFIED ZIP %s", path.name)
        yield handle
    finally:
        handle.close()


def download_one(
    url: str,
    target: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    timeout_seconds: float,
    opener: Opener,
) -> None:
    validate_public_https_url(url)
    LOG.info("GET %s", url)
    LOG.info("EXPECT %s bytes; SHA-256 %s", f"{expected_size:,}", expected_sha256)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"},
        method="GET",
    )
    try:
        response_context = opener(request, timeout=timeout_seconds)
        with response_context as response:
            status = getattr(response, "status", None)
            if status is not None and status != 200:
                raise ReconstructionError(f"HTTP status {status} for archive download")
            final_url = response.geturl() if hasattr(response, "geturl") else url
            validate_public_https_url(final_url, reject_secret_query=False)
            encoding = response.headers.get("Content-Encoding", "identity").lower()
            if encoding not in {"", "identity"}:
                raise ReconstructionError(
                    f"Unexpected HTTP content encoding for ZIP bytes: {encoding}"
                )
            length_text = response.headers.get("Content-Length")
            if length_text is not None:
                try:
                    length = int(length_text)
                except ValueError as exc:
                    raise ReconstructionError(
                        f"Invalid HTTP Content-Length: {length_text!r}"
                    ) from exc
                if length != expected_size:
                    raise ReconstructionError(
                        f"HTTP Content-Length differs: {length} != {expected_size}"
                    )
            target.parent.mkdir(parents=True, exist_ok=True)
            digest = hashlib.sha256()
            count = 0
            last = time.monotonic()
            with target.open("xb") as destination:
                while True:
                    block = response.read(archive.CHUNK)
                    if not block:
                        break
                    count += len(block)
                    if count > expected_size:
                        raise ReconstructionError(
                            "Downloaded archive exceeds manifest size"
                        )
                    destination.write(block)
                    digest.update(block)
                    last = progress("DOWNLOAD", count, expected_size, last)
    except ReconstructionError:
        raise
    except Exception as exc:
        raise ReconstructionError(f"Archive download failed: {exc}") from exc
    if count != expected_size or digest.hexdigest() != expected_sha256:
        raise ReconstructionError(
            f"Downloaded archive size/SHA-256 differs: {target.name}"
        )
    os.chmod(target, 0o644)
    LOG.info("VERIFIED DOWNLOAD %s bytes; SHA-256 %s", f"{count:,}", expected_sha256)


def _verified_partial_prefix(
    partial: Path,
    segments: Sequence[Mapping[str, Any]],
    *,
    canonical_size: int,
) -> tuple[int, int, int]:
    """Verify complete segment boundaries in a resumable canonical-ZIP partial."""
    try:
        observed = partial.lstat()
    except FileNotFoundError:
        return 0, 0, 0
    if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
        raise ReconstructionError(
            f"Canonical archive partial is linked or not regular: {partial}"
        )
    observed_size = observed.st_size
    if observed_size > canonical_size:
        raise ReconstructionError(
            f"Canonical archive partial exceeds its bound size: {partial}"
        )

    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(partial, flags)
    completed = 0
    boundary = 0
    with os.fdopen(descriptor, "rb") as handle:
        opened = os.fstat(handle.fileno())
        if not stat.S_ISREG(opened.st_mode) or opened.st_size != observed_size:
            raise ReconstructionError(
                f"Canonical archive partial changed during verification: {partial}"
            )
        for segment in segments:
            segment_size = segment["size_bytes"]
            if observed_size - boundary < segment_size:
                break
            digest = hashlib.sha256()
            remaining = segment_size
            while remaining:
                block = handle.read(min(archive.CHUNK, remaining))
                if not block:
                    raise ReconstructionError(
                        f"Canonical archive partial ended unexpectedly: {partial}"
                    )
                digest.update(block)
                remaining -= len(block)
            if digest.hexdigest() != segment["sha256"]:
                raise ReconstructionError(
                    "Canonical archive partial differs at a completed "
                    f"segment boundary: {partial}"
                )
            boundary += segment_size
            completed += 1
    return completed, boundary, observed_size


def _truncate_partial(partial: Path, *, observed_size: int, boundary: int) -> None:
    flags = os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(partial, flags)
    with os.fdopen(descriptor, "r+b") as handle:
        opened = os.fstat(handle.fileno())
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_size != observed_size
            or not 0 <= boundary <= observed_size
        ):
            raise ReconstructionError(
                f"Canonical archive partial changed before rollback: {partial}"
            )
        os.ftruncate(handle.fileno(), boundary)
        os.fsync(handle.fileno())


def _append_verified_segment(
    partial: Path,
    segment_path: Path,
    segment: Mapping[str, Any],
    *,
    boundary: int,
) -> None:
    """Append one twice-verified segment and roll back an interrupted write."""
    _require_exact_file(
        segment_path,
        expected_size=segment["size_bytes"],
        expected_sha256=segment["sha256"],
        label="Downloaded transport segment",
    )
    source_flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    destination_flags = os.O_RDWR | getattr(os, "O_NOFOLLOW", 0)
    source_descriptor = os.open(segment_path, source_flags)
    try:
        destination_descriptor = os.open(partial, destination_flags)
    except Exception:
        os.close(source_descriptor)
        raise
    with (
        os.fdopen(source_descriptor, "rb") as source,
        os.fdopen(destination_descriptor, "r+b") as destination,
    ):
        source_stat = os.fstat(source.fileno())
        destination_stat = os.fstat(destination.fileno())
        if (
            not stat.S_ISREG(source_stat.st_mode)
            or source_stat.st_size != segment["size_bytes"]
            or not stat.S_ISREG(destination_stat.st_mode)
            or destination_stat.st_size != boundary
        ):
            raise ReconstructionError(
                f"Segment or canonical partial changed before append: {partial}"
            )
        destination.seek(0, os.SEEK_END)
        digest = hashlib.sha256()
        count = 0
        try:
            for block in iter(lambda: source.read(archive.CHUNK), b""):
                count += len(block)
                if count > segment["size_bytes"]:
                    raise ReconstructionError(
                        f"Transport segment exceeds its bound size: {segment_path}"
                    )
                destination.write(block)
                digest.update(block)
            if (
                count != segment["size_bytes"]
                or digest.hexdigest() != segment["sha256"]
            ):
                raise ReconstructionError(
                    f"Transport segment changed during append: {segment_path}"
                )
            destination.flush()
            os.fsync(destination.fileno())
        except Exception:
            try:
                os.ftruncate(destination.fileno(), boundary)
                os.fsync(destination.fileno())
            except OSError as rollback_error:
                raise ReconstructionError(
                    f"Cannot roll back interrupted segment append: {partial}"
                ) from rollback_error
            raise


def download_segmented_archive(
    entry: Mapping[str, Any],
    target: Path,
    *,
    timeout_seconds: float,
    opener: Opener,
) -> None:
    """Resume one canonical ZIP while retaining at most one transport segment."""
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_suffix(target.suffix + ".part")
    segment_target = target.with_suffix(target.suffix + ".segment")
    if target.exists() or target.is_symlink():
        raise ReconstructionError(
            f"Segmented archive target already exists unexpectedly: {target}"
        )

    if not partial.exists() and not partial.is_symlink():
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = os.open(partial, flags, 0o600)
        os.close(descriptor)

    completed, boundary, observed_size = _verified_partial_prefix(
        partial,
        entry["segments"],
        canonical_size=entry["size_bytes"],
    )
    LOG.info("RESUME %s: %s/%s chunks already verified (%s bytes)", target.name,
             completed, len(entry["segments"]), f"{boundary:,}")
    if observed_size != boundary:
        _truncate_partial(
            partial,
            observed_size=observed_size,
            boundary=boundary,
        )

    if completed == len(entry["segments"]):
        if segment_target.exists() or segment_target.is_symlink():
            observed = segment_target.lstat()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
                raise ReconstructionError(
                    f"Transport-segment cache is linked or not regular: {segment_target}"
                )
            segment_target.unlink()
        _require_exact_file(
            partial,
            expected_size=entry["size_bytes"],
            expected_sha256=entry["sha256"],
            label="Concatenated archive",
        )
        os.chmod(partial, 0o644)
        os.replace(partial, target)
        return

    for index, segment in enumerate(entry["segments"][completed:], start=completed + 1):
        LOG.info("CHUNK %s/%s for %s: %s", index, len(entry["segments"]),
                 target.name, segment["name"])
        if segment_target.exists() or segment_target.is_symlink():
            if not _regular_file_matches(
                segment_target,
                expected_size=segment["size_bytes"],
                expected_sha256=segment["sha256"],
                label="Cached transport segment",
            ):
                segment_target.unlink()
        if not segment_target.exists():
            for attempt in range(1, SEGMENT_DOWNLOAD_ATTEMPTS + 1):
                try:
                    download_one(
                        segment["url"],
                        segment_target,
                        expected_size=segment["size_bytes"],
                        expected_sha256=segment["sha256"],
                        timeout_seconds=timeout_seconds,
                        opener=opener,
                    )
                    break
                except ReconstructionError as exc:
                    try:
                        segment_target.unlink()
                    except FileNotFoundError:
                        pass
                    if attempt == SEGMENT_DOWNLOAD_ATTEMPTS:
                        raise
                    delay = min(
                        SEGMENT_RETRY_INITIAL_SECONDS * (2 ** (attempt - 1)),
                        SEGMENT_RETRY_MAX_SECONDS,
                    )
                    LOG.info("RETRY %s/%s in %.0f seconds: %s", attempt + 1,
                             SEGMENT_DOWNLOAD_ATTEMPTS, delay, exc)
                    time.sleep(delay)

        _append_verified_segment(
            partial,
            segment_target,
            segment,
            boundary=boundary,
        )
        boundary += segment["size_bytes"]
        segment_target.unlink()
        LOG.info("ASSEMBLED %s: %s / %s verified bytes", target.name,
                 f"{boundary:,}", f"{entry['size_bytes']:,}")

    LOG.info("VERIFY concatenated ZIP %s; expected SHA-256 %s", target.name, entry["sha256"])
    _require_exact_file(
        partial,
        expected_size=entry["size_bytes"],
        expected_sha256=entry["sha256"],
        label="Concatenated archive",
    )
    os.chmod(partial, 0o644)
    os.replace(partial, target)
    LOG.info("VERIFIED canonical ZIP %s", target)


def _expected_resume_receipt(
    document: Mapping[str, Any],
    urls: Mapping[str, str | Mapping[str, Any]],
) -> dict[str, Any]:
    url_payload = archive.canonical_json(urls)
    return {
        "schema": RESUME_SCHEMA,
        "manifest": _manifest_binding(document),
        "url_map": {
            "size_bytes": len(url_payload),
            "sha256": hashlib.sha256(url_payload).hexdigest(),
        },
    }


def _resume_inventory(
    document: Mapping[str, Any],
) -> tuple[set[str], set[str], list[str]]:
    directories: set[str] = set()
    allowed_files = {RESUME_RECEIPT_NAME}
    canonical_paths: list[str] = []
    for shard in _manifest_shards(document):
        relative = archive.safe_relative(shard["relative_archive_path"])
        canonical_paths.append(relative)
        pure = PurePosixPath(relative)
        parent = pure.parent
        while parent != PurePosixPath("."):
            directories.add(parent.as_posix())
            parent = parent.parent
        allowed_files.update({relative, relative + ".part", relative + ".segment"})
    return directories, allowed_files, canonical_paths


def _validate_resume_stage(
    stage: Path,
    document: Mapping[str, Any],
    receipt: Mapping[str, Any],
) -> None:
    try:
        observed_stage = stage.lstat()
    except FileNotFoundError as exc:
        raise ReconstructionError(f"Resume staging directory is missing: {stage}") from exc
    if stat.S_ISLNK(observed_stage.st_mode) or not stat.S_ISDIR(
        observed_stage.st_mode
    ):
        raise ReconstructionError(
            f"Resume staging path is linked or not a directory: {stage}"
        )
    observed_receipt = read_json_object(
        stage / RESUME_RECEIPT_NAME,
        label="download resume receipt",
        max_bytes=MAX_RESUME_RECEIPT_BYTES,
    )
    if observed_receipt != receipt:
        raise ReconstructionError(
            "Retained download is bound to a different manifest or URL map"
        )

    expected_directories, allowed_files, canonical_paths = _resume_inventory(document)
    observed_directories: set[str] = set()
    observed_files: set[str] = set()
    for directory, dirnames, filenames in os.walk(stage, followlinks=False):
        base = Path(directory)
        for name in dirnames:
            path = base / name
            relative = path.relative_to(stage).as_posix()
            observed = path.lstat()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISDIR(observed.st_mode):
                raise ReconstructionError(
                    f"Unsafe object in retained download: {relative}"
                )
            observed_directories.add(relative)
        for name in filenames:
            path = base / name
            relative = path.relative_to(stage).as_posix()
            observed = path.lstat()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
                raise ReconstructionError(
                    f"Unsafe object in retained download: {relative}"
                )
            observed_files.add(relative)

    if observed_directories != expected_directories:
        raise ReconstructionError(
            "Retained download directory inventory differs; "
            f"missing={sorted(expected_directories - observed_directories)}, "
            f"extra={sorted(observed_directories - expected_directories)}"
        )
    if not observed_files <= allowed_files:
        raise ReconstructionError(
            "Retained download contains unexpected files: "
            f"{sorted(observed_files - allowed_files)}"
        )

    first_incomplete: str | None = None
    active_partial: str | None = None
    for relative in canonical_paths:
        complete = relative in observed_files
        partial = relative + ".part" in observed_files
        segment = relative + ".segment" in observed_files
        if complete and (partial or segment):
            raise ReconstructionError(
                f"Retained canonical archive has conflicting staging files: {relative}"
            )
        if complete:
            if first_incomplete is not None:
                raise ReconstructionError(
                    f"Retained canonical archives are not a strict prefix: {relative}"
                )
            continue
        if first_incomplete is None:
            first_incomplete = relative
        if segment and not partial:
            raise ReconstructionError(
                f"Retained transport segment lacks its archive partial: {relative}"
            )
        if partial or segment:
            if relative != first_incomplete or active_partial is not None:
                raise ReconstructionError(
                    "Retained download has more than one active canonical archive"
                )
            active_partial = relative


def _initialize_resume_stage(
    stage: Path,
    document: Mapping[str, Any],
    receipt: Mapping[str, Any],
) -> None:
    try:
        stage.mkdir(mode=0o700)
        directories, _, _ = _resume_inventory(document)
        for relative in sorted(
            directories,
            key=lambda value: (value.count("/"), value),
        ):
            stage.joinpath(*PurePosixPath(relative).parts).mkdir(mode=0o755)
        receipt_path = stage / RESUME_RECEIPT_NAME
        payload = archive.canonical_json(receipt)
        with receipt_path.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(receipt_path, 0o600)
        _validate_resume_stage(stage, document, receipt)
    except Exception:
        if stage.exists() and not stage.is_symlink():
            shutil.rmtree(stage)
        raise


def download_archives(
    document: Mapping[str, Any],
    urls: Mapping[str, str | Mapping[str, Any]],
    destination: Path,
    *,
    timeout_seconds: float = 120.0,
    opener: Opener = urllib.request.urlopen,
) -> Path:
    if not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ReconstructionError("Download timeout must be positive")
    target = _resolved_new_path(destination, label="download directory")
    target.parent.mkdir(parents=True, exist_ok=True)

    sources = list(urls.values())
    segmented = bool(sources) and all(
        isinstance(source, Mapping) for source in sources
    )
    direct = bool(sources) and all(isinstance(source, str) for source in sources)
    if not (segmented or direct):
        raise ReconstructionError("Normalized URL map mixes archive source types")

    if segmented:
        stage = target.with_name(f".{target.name}.downloading")
        receipt = _expected_resume_receipt(document, urls)
    else:
        stage = target.with_name(f".{target.name}.downloading.{os.getpid()}")
        receipt = None

    LOG.info("DOWNLOAD STAGING %s; completed archives will be installed at %s", stage, target)

    try:
        if stage.exists() or stage.is_symlink():
            if not segmented:
                raise ReconstructionError(
                    f"Download staging path already exists: {stage}"
                )
            if receipt is None:
                raise AssertionError("Segmented download receipt is missing")
            _validate_resume_stage(stage, document, receipt)
        elif segmented:
            if receipt is None:
                raise AssertionError("Segmented download receipt is missing")
            _initialize_resume_stage(stage, document, receipt)
        else:
            stage.mkdir(mode=0o700)
            for record in archive.RECORDS:
                (stage / record).mkdir(mode=0o755)

        shards = _manifest_shards(document)
        for index, shard in enumerate(shards, start=1):
            relative = shard["relative_archive_path"]
            LOG.info("ARCHIVE %s/%s: %s", index, len(shards), relative)
            target_path = stage.joinpath(*PurePosixPath(relative).parts)
            source = urls[relative]
            if isinstance(source, str):
                download_one(
                    source,
                    target_path,
                    expected_size=shard["compressed_bytes"],
                    expected_sha256=shard["sha256"],
                    timeout_seconds=timeout_seconds,
                    opener=opener,
                )
            elif isinstance(source, Mapping):
                if target_path.exists() or target_path.is_symlink():
                    LOG.info("RECHECK retained ZIP %s", target_path)
                    _require_exact_file(
                        target_path,
                        expected_size=shard["compressed_bytes"],
                        expected_sha256=shard["sha256"],
                        label="Retained canonical archive",
                    )
                    LOG.info("REUSED verified ZIP %s", target_path)
                    continue
                download_segmented_archive(
                    source,
                    target_path,
                    timeout_seconds=timeout_seconds,
                    opener=opener,
                )
            else:
                raise ReconstructionError(
                    f"Invalid normalized URL-map entry: {relative}"
                )

        if segmented:
            if receipt is None:
                raise AssertionError("Segmented download receipt is missing")
            _validate_resume_stage(stage, document, receipt)
        LOG.info("All archives verified; installing download directory %s", target)
        os.chmod(stage, 0o755)
        os.replace(stage, target)
        if segmented:
            (target / RESUME_RECEIPT_NAME).unlink()
    except Exception:
        if not segmented and stage.exists() and not stage.is_symlink():
            shutil.rmtree(stage)
        raise
    return target


def _files_by_shard(
    document: Mapping[str, Any],
) -> dict[tuple[str, str], list[Mapping[str, Any]]]:
    grouped: dict[tuple[str, str], list[Mapping[str, Any]]] = {}
    for row in document["tree"]["files"]:
        grouped.setdefault((row["record"], row["shard"]), []).append(row)
    return grouped


def extract_shard(
    archive_path: Path,
    shard: Mapping[str, Any],
    file_rows: Sequence[Mapping[str, Any]],
    stage: Path,
) -> set[str]:
    expected_members = [row["member"] for row in file_rows]
    extracted: set[str] = set()
    with verified_archive_handle(archive_path, shard) as handle:
        try:
            with zipfile.ZipFile(handle, mode="r") as zipped:
                infos = zipped.infolist()
                if [info.filename for info in infos] != expected_members:
                    raise ReconstructionError(
                        f"ZIP member inventory/order differs: {archive_path.name}"
                    )
                if len(infos) != len(set(expected_members)):
                    raise ReconstructionError(
                        f"ZIP contains duplicate members: {archive_path.name}"
                    )
                for index, (info, row) in enumerate(zip(infos, file_rows), start=1):
                    expected_mode = int(row["mode"], 8)
                    encoded_mode = info.external_attr >> 16
                    if (
                        info.is_dir()
                        or info.flag_bits & 0x1
                        or info.compress_type != zipfile.ZIP_DEFLATED
                        or info.file_size != row["size_bytes"]
                        or stat.S_IFMT(encoded_mode) != stat.S_IFREG
                        or stat.S_IMODE(encoded_mode) != expected_mode
                    ):
                        raise ReconstructionError(
                            f"ZIP member metadata differs: {info.filename}"
                        )
                    target = stage.joinpath(*PurePosixPath(row["path"]).parts)
                    if not target.parent.is_dir() or target.parent.is_symlink():
                        raise ReconstructionError(
                            f"Manifest parent directory is unavailable: {row['path']}"
                        )
                    digest = hashlib.sha256()
                    count = 0
                    last = time.monotonic()
                    with zipped.open(info, mode="r") as source, target.open(
                        "xb"
                    ) as destination:
                        while True:
                            block = source.read(archive.CHUNK)
                            if not block:
                                break
                            count += len(block)
                            if count > row["size_bytes"]:
                                raise ReconstructionError(
                                    f"ZIP member exceeds manifest size: {info.filename}"
                                )
                            destination.write(block)
                            digest.update(block)
                            last = progress(f"EXTRACT {row['path']}", count, row["size_bytes"], last)
                    if (
                        count != row["size_bytes"]
                        or digest.hexdigest() != row["sha256"]
                    ):
                        raise ReconstructionError(
                            f"ZIP member size/SHA-256 differs: {info.filename}"
                        )
                    os.chmod(target, expected_mode)
                    extracted.add(row["path"])
                    LOG.info("VERIFIED MEMBER %s/%s in %s: %s (%s bytes)",
                             index, len(file_rows), archive_path.name, row["path"], f"{count:,}")
        except ReconstructionError:
            raise
        except Exception as exc:
            raise ReconstructionError(
                f"Cannot safely read ZIP {archive_path}: {exc}"
            ) from exc
    return extracted


def verify_reconstructed_tree(
    stage: Path,
    document: Mapping[str, Any],
    extracted: set[str],
) -> None:
    expected_directories = {
        row["path"]: int(row["mode"], 8) for row in document["tree"]["directories"]
    }
    expected_files = {
        row["path"]: (row["size_bytes"], int(row["mode"], 8))
        for row in document["tree"]["files"]
    }
    observed_directories: dict[str, int] = {}
    observed_files: dict[str, tuple[int, int]] = {}
    for directory, dirnames, filenames in os.walk(stage, followlinks=False):
        base = Path(directory)
        for name in dirnames:
            path = base / name
            observed = path.lstat()
            relative = path.relative_to(stage).as_posix()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISDIR(observed.st_mode):
                raise ReconstructionError(f"Unexpected directory entry: {relative}")
            observed_directories[relative] = stat.S_IMODE(observed.st_mode)
        for name in filenames:
            path = base / name
            observed = path.lstat()
            relative = path.relative_to(stage).as_posix()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
                raise ReconstructionError(f"Unexpected file entry: {relative}")
            observed_files[relative] = (
                observed.st_size,
                stat.S_IMODE(observed.st_mode),
            )
    if extracted != set(expected_files):
        raise ReconstructionError("Not every manifest file was extracted exactly once")
    if (
        observed_directories != expected_directories
        or observed_files != expected_files
        or stat.S_IMODE(stage.stat().st_mode) != int(document["tree"]["root_mode"], 8)
    ):
        raise ReconstructionError("Reconstructed tree paths, sizes, or modes differ")


def reconstruct_tree(
    document: Mapping[str, Any],
    archive_root: Path,
    output: Path,
) -> Path:
    source = _resolved_existing_directory(archive_root, label="archive root")
    target = _resolved_new_path(output, label="output Paper tree")
    if _paths_overlap(source, target):
        raise ReconstructionError(
            "Archive root and reconstructed tree must not overlap"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = target.with_name(f".{target.name}.reconstructing.{os.getpid()}")
    if stage.exists() or stage.is_symlink():
        raise ReconstructionError(f"Reconstruction staging path exists: {stage}")

    grouped = _files_by_shard(document)
    extracted: set[str] = set()
    LOG.info("RECONSTRUCT %s files, %s bytes into %s", len(document["tree"]["files"]),
             f"{document['tree']['payload_bytes']:,}", stage)
    try:
        stage.mkdir(mode=0o700)
        for row in document["tree"]["directories"]:
            path = stage.joinpath(*PurePosixPath(row["path"]).parts)
            path.mkdir(mode=0o700)
        for record in archive.RECORDS:
            for shard in document["records"][record]["shards"]:
                key = (record, shard["name"])
                rows = grouped.get(key, [])
                if len(rows) != shard["member_count"]:
                    raise ReconstructionError(
                        f"Manifest member count differs: {shard['name']}"
                    )
                path = _archive_path(source, shard["relative_archive_path"])
                LOG.info("EXTRACT ARCHIVE %s (%s members)", path.name, len(rows))
                new_paths = extract_shard(path, shard, rows, stage)
                if extracted & new_paths:
                    raise ReconstructionError("A file was extracted more than once")
                extracted |= new_paths
        for row in sorted(
            document["tree"]["directories"],
            key=lambda value: (value["path"].count("/"), value["path"]),
            reverse=True,
        ):
            os.chmod(
                stage.joinpath(*PurePosixPath(row["path"]).parts),
                int(row["mode"], 8),
            )
        os.chmod(stage, int(document["tree"]["root_mode"], 8))
        verify_reconstructed_tree(stage, document, extracted)
        LOG.info("VERIFIED complete tree inventory and modes; manifest tree SHA-256 %s",
                 document["tree"]["sha256"])
        os.replace(stage, target)
        LOG.info("INSTALLED exact Paper tree at %s", target)
    except Exception:
        if stage.exists() and not stage.is_symlink():
            shutil.rmtree(stage)
        raise
    return target


def reconstruct(
    manifest_path: Path,
    output: Path,
    *,
    archive_root: Path | None = None,
    url_map_path: Path | None = None,
    download_dir: Path | None = None,
    timeout_seconds: float = 120.0,
    opener: Opener = urllib.request.urlopen,
) -> tuple[Path, Path]:
    document = load_manifest(manifest_path)
    prospective_output = _resolved_new_path(
        output, label="reconstructed destination"
    )
    if (archive_root is None) == (url_map_path is None):
        raise ReconstructionError(
            "Choose exactly one archive source: local archive root or HTTPS URL map"
        )
    if archive_root is not None:
        if download_dir is not None:
            raise ReconstructionError(
                "Download directory is only valid with an HTTPS URL map"
            )
        source = _resolved_existing_directory(archive_root, label="archive root")
    else:
        if url_map_path is None or download_dir is None:
            raise ReconstructionError(
                "HTTPS URL-map mode requires --url-map and --download-dir"
            )
        urls = load_url_map(url_map_path, document)
        prospective_download = (
            download_dir.expanduser().absolute().parent.resolve()
            / download_dir.expanduser().absolute().name
        )
        if _paths_overlap(prospective_download, prospective_output):
            raise ReconstructionError(
                "Download directory and reconstructed tree must not overlap"
            )
        if all(isinstance(source, Mapping) for source in urls.values()):
            prospective_stage = prospective_download.with_name(
                f".{prospective_download.name}.downloading"
            )
            if _paths_overlap(prospective_stage, prospective_output):
                raise ReconstructionError(
                    "Resume staging directory and reconstructed tree "
                    "must not overlap"
                )
        source = download_archives(
            document,
            urls,
            download_dir,
            timeout_seconds=timeout_seconds,
            opener=opener,
        )
    result = reconstruct_tree(document, source, output)
    return source, result


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="New destination directory for the exact reconstructed Paper tree.",
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--archive-root",
        type=Path,
        help="Existing directory containing record/shard ZIP paths.",
    )
    source.add_argument(
        "--url-map",
        type=Path,
        help=(
            "JSON URL map using schema apotome-paper-archive-url-map-v1 for "
            "canonical ZIP URLs or v2 for ordered transport-segment URLs."
        ),
    )
    parser.add_argument(
        "--download-dir",
        type=Path,
        help=(
            "New persistent archive directory; required with --url-map. "
            "Interrupted URL-map-v2 downloads resume on an exact rerun."
        ),
    )
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument("--verbose", action="store_true",
                        help="Log every chunk, retry, checksum, archive and extracted member.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    try:
        source, output = reconstruct(
            args.manifest,
            args.output,
            archive_root=args.archive_root,
            url_map_path=args.url_map,
            download_dir=args.download_dir,
            timeout_seconds=args.timeout_seconds,
        )
    except ReconstructionError:
        raise
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise ReconstructionError(f"Reconstruction failed: {exc}") from exc
    print(f"[PASS] verified archive root: {source}")
    print(f"[PASS] exact reconstructed Paper tree: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ReconstructionError as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        raise SystemExit(2)
