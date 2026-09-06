#!/usr/bin/env python3
"""Build and validate the small public GitHub reconstruction repository.

This module deliberately has no GitHub or Zenodo mutation code. A final
repository can only be assembled from a valid release manifest, an exact map
of public Zenodo file URLs, and four resolved Zenodo record identities.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import urllib.parse
from pathlib import Path
from typing import Any, Mapping, Sequence

import reconstruct_paper


METADATA_SCHEMA = "apotome-github-publication-metadata-v1"
TRANSPORT_SCHEMA_V1 = "allulose-zenodo-segmented-transport-v1"
TRANSPORT_SCHEMA = "allulose-zenodo-segmented-transport-v2"
TRANSPORT_SCHEMAS = {TRANSPORT_SCHEMA_V1, TRANSPORT_SCHEMA}
TRANSPORT_REASSEMBLY = (
    "Concatenate each archive's segments in ascending number order; "
    "the resulting ZIP must match the recorded size and SHA-256."
)
REPOSITORY_URL = "https://github.com/cfarkas/allulose-sucrose-hypothalamus"
RELEASE_TAG = "v1.0.0"
RECORDS = (
    "software_core",
    "main_figure_data",
    "supplementary_figure_data",
    "figs3_wsi",
)
EXPECTED_SHARD_COUNTS = {
    "software_core": 1,
    "main_figure_data": 4,
    "supplementary_figure_data": 1,
    "figs3_wsi": 4,
}
PLACEHOLDER_RE = re.compile(
    r"(?:@@[A-Z0-9_]+@@|REPLACE(?:[_ -]?ME)?|CHANGEME|YOUR[_-]|"
    r"example\.com|\bTODO\b|\bTBD\b)",
    re.IGNORECASE,
)
RECORD_ID_RE = re.compile(r"^[1-9][0-9]*$")
RELEASE_ID_RE = re.compile(r"^apotome-[0-9a-f]{24}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")
LICENSE_SHA256 = {
    "LICENSE": "ab0fe511c596148fad61fd320619bf91c0b71b4083c5821948175dcfc87b36ed",
    "LICENSES/MIT.txt": "212090e1e8934e8a110cb6aefa0db81589ba4997cf4e79bd754be4d198bec6c0",
    "LICENSES/CC-BY-4.0.txt": "63cb76cfde49829d3e28968762415559967f6dcd70e6025ef3a079f5e4c67a58",
}
MAX_JSON_BYTES = 1_000_000_000


class PublicationError(RuntimeError):
    """A GitHub publication input or assembled repository is invalid."""


def _unique_object(pairs: Sequence[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise PublicationError(f"Duplicate JSON key: {key!r}")
        result[key] = value
    return result


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink():
        raise PublicationError(f"{label} may not be a symlink: {lexical}")
    try:
        observed = lexical.lstat()
    except FileNotFoundError as exc:
        raise PublicationError(f"{label} is missing: {lexical}") from exc
    if not stat.S_ISREG(observed.st_mode) or observed.st_size > MAX_JSON_BYTES:
        raise PublicationError(
            f"{label} is not a suitably sized regular file: {lexical}"
        )
    try:
        document = json.loads(
            lexical.read_text(encoding="utf-8"), object_pairs_hook=_unique_object
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PublicationError(f"Cannot read valid {label}: {lexical}: {exc}") from exc
    if not isinstance(document, dict):
        raise PublicationError(f"{label} root must be a JSON object")
    return document


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require_exact_keys(
    document: Mapping[str, Any], expected: set[str], *, label: str
) -> None:
    observed = set(document)
    if observed != expected:
        raise PublicationError(
            f"{label} keys differ; missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )


def _reject_placeholders(value: Any, *, label: str) -> None:
    if isinstance(value, str):
        if PLACEHOLDER_RE.search(value):
            raise PublicationError(f"{label} contains an unresolved placeholder")
    elif isinstance(value, Mapping):
        for key, child in value.items():
            _reject_placeholders(child, label=f"{label}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _reject_placeholders(child, label=f"{label}[{index}]")


def validate_metadata_document(document: Mapping[str, Any]) -> dict[str, Any]:
    _require_exact_keys(
        document,
        {"schema", "repository_url", "release_tag", "records"},
        label="publication metadata",
    )
    _reject_placeholders(document, label="publication metadata")
    if document["schema"] != METADATA_SCHEMA:
        raise PublicationError("Publication metadata schema differs")
    if document["repository_url"] != REPOSITORY_URL:
        raise PublicationError(f"Repository URL must be exactly {REPOSITORY_URL}")
    if document["release_tag"] != RELEASE_TAG:
        raise PublicationError(f"Release tag must be exactly {RELEASE_TAG}")
    records = document["records"]
    if not isinstance(records, dict):
        raise PublicationError("Publication metadata records must be an object")
    _require_exact_keys(records, set(RECORDS), label="publication metadata records")
    seen_ids: set[str] = set()
    seen_dois: set[str] = set()
    normalized: dict[str, dict[str, str]] = {}
    for record in RECORDS:
        row = records[record]
        if not isinstance(row, dict):
            raise PublicationError(f"Metadata row for {record} must be an object")
        _require_exact_keys(
            row,
            {"record_id", "record_url", "doi", "doi_url"},
            label=f"metadata row {record}",
        )
        if any(not isinstance(row[key], str) for key in row):
            raise PublicationError(f"Metadata values for {record} must be strings")
        record_id = row["record_id"]
        if not RECORD_ID_RE.fullmatch(record_id):
            raise PublicationError(f"Invalid Zenodo record ID for {record}")
        doi = f"10.5281/zenodo.{record_id}"
        expected = {
            "record_id": record_id,
            "record_url": f"https://zenodo.org/records/{record_id}",
            "doi": doi,
            "doi_url": f"https://doi.org/{doi}",
        }
        if row != expected:
            raise PublicationError(
                f"Zenodo DOI/URL values are inconsistent for {record}; expected {expected}"
            )
        if record_id in seen_ids or doi in seen_dois:
            raise PublicationError(
                "Each logical record must have a distinct Zenodo ID/DOI"
            )
        seen_ids.add(record_id)
        seen_dois.add(doi)
        normalized[record] = expected
    return {
        "schema": METADATA_SCHEMA,
        "repository_url": REPOSITORY_URL,
        "release_tag": RELEASE_TAG,
        "records": normalized,
    }


def load_metadata(path: Path) -> dict[str, Any]:
    return validate_metadata_document(load_json(path, label="publication metadata"))


def load_project_manifest(path: Path) -> dict[str, Any]:
    try:
        document = reconstruct_paper.load_manifest(path)
    except Exception as exc:
        raise PublicationError(f"Release manifest is invalid: {exc}") from exc
    counts = {record: len(document["records"][record]["shards"]) for record in RECORDS}
    if counts != EXPECTED_SHARD_COUNTS:
        raise PublicationError(
            f"Expected the reviewed 1+4+1+4 archive layout, found {counts}"
        )
    if sum(counts.values()) != 10:
        raise PublicationError("Exactly ten Zenodo archive files are required")
    return document


def _canonical_shards(
    manifest: Mapping[str, Any],
) -> dict[str, tuple[str, Mapping[str, Any]]]:
    return {
        shard["relative_archive_path"]: (record, shard)
        for record, shard in manifest_shards(manifest)
    }


def load_project_transport_manifest(
    path: Path, manifest_path: Path, manifest: Mapping[str, Any]
) -> dict[str, Any]:
    """Validate the sealed Zenodo chunk inventory against the canonical ZIPs."""
    document = load_json(path, label="transport manifest")
    if set(document) != {
        "schema",
        "release_id",
        "canonical_manifest",
        "records",
        "reassembly",
    }:
        raise PublicationError("Transport manifest fields differ")
    expected_binding = {
        "sha256": sha256_file(manifest_path),
        "size_bytes": manifest_path.stat().st_size,
        "tree_sha256": manifest["tree"]["sha256"],
        "tree_files": manifest["tree"]["file_count"],
    }
    schema = document.get("schema")
    if (
        schema not in TRANSPORT_SCHEMAS
        or not RELEASE_ID_RE.fullmatch(str(document.get("release_id", "")))
        or document.get("canonical_manifest") != expected_binding
        or document.get("reassembly") != TRANSPORT_REASSEMBLY
    ):
        raise PublicationError("Transport manifest binding differs")
    records = document.get("records")
    if not isinstance(records, dict):
        raise PublicationError("Transport-manifest records must be an object")
    _require_exact_keys(records, set(RECORDS), label="transport-manifest records")

    canonical = _canonical_shards(manifest)
    seen_segments: set[tuple[str, str]] = set()
    physical_names: dict[str, set[str]] = {record: set() for record in RECORDS}
    physical_names["software_core"].update({"manifest.json", "transport_manifest.json"})
    for record in RECORDS:
        record_row = records[record]
        if (
            not isinstance(record_row, dict)
            or set(record_row) != {"archives"}
            or not isinstance(record_row["archives"], list)
        ):
            raise PublicationError(f"Transport record is invalid: {record}")
        expected_relatives = [
            shard["relative_archive_path"]
            for current_record, shard in manifest_shards(manifest)
            if current_record == record
        ]
        archives = record_row["archives"]
        if (
            len(archives) != len(expected_relatives)
            or [
                row.get("relative_local_path")
                for row in archives
                if isinstance(row, dict)
            ]
            != expected_relatives
        ):
            raise PublicationError(f"Transport archive inventory differs: {record}")
        for archive_row, relative in zip(archives, expected_relatives):
            _, shard = canonical[relative]
            if not isinstance(archive_row, dict) or set(archive_row) != {
                "name",
                "relative_local_path",
                "size_bytes",
                "sha256",
                "md5",
                "segment_bytes",
                "segments",
            }:
                raise PublicationError(f"Transport archive fields differ: {relative}")
            segment_bytes = archive_row.get("segment_bytes")
            segments = archive_row.get("segments")
            if (
                archive_row.get("name") != shard["name"]
                or archive_row.get("size_bytes") != shard["compressed_bytes"]
                or archive_row.get("sha256") != shard["sha256"]
                or not MD5_RE.fullmatch(str(archive_row.get("md5", "")))
                or not isinstance(segment_bytes, int)
                or isinstance(segment_bytes, bool)
                or segment_bytes <= 0
                or not isinstance(segments, list)
                or not segments
                or len(segments) > 999
            ):
                raise PublicationError(f"Transport archive binding differs: {relative}")
            offset = 0
            count = len(segments)
            for number, segment in enumerate(segments, start=1):
                expected_name = (
                    f"{archive_row['name']}.chunk{number:03d}-of-{count:03d}"
                )
                expected_segment_keys = {
                    "name",
                    "number",
                    "offset_bytes",
                    "size_bytes",
                    "sha256",
                    "md5",
                }
                segment_keys = set(segment) if isinstance(segment, dict) else set()
                storage_record = (
                    segment.get("storage_record", record)
                    if isinstance(segment, dict)
                    else None
                )
                if (
                    not isinstance(segment, dict)
                    or (
                        segment_keys != expected_segment_keys
                        and not (
                            schema == TRANSPORT_SCHEMA
                            and segment_keys
                            == expected_segment_keys | {"storage_record"}
                        )
                    )
                    or storage_record not in RECORDS
                    or ("storage_record" in segment_keys and storage_record == record)
                ):
                    raise PublicationError(
                        f"Transport segment fields differ: {relative} #{number}"
                    )
                size = segment.get("size_bytes")
                identity = (record, str(segment.get("name")))
                if (
                    segment.get("name") != expected_name
                    or identity in seen_segments
                    or segment.get("number") != number
                    or segment.get("offset_bytes") != offset
                    or not isinstance(size, int)
                    or isinstance(size, bool)
                    or size <= 0
                    or size > segment_bytes
                    or (number < count and size != segment_bytes)
                    or not SHA256_RE.fullmatch(str(segment.get("sha256", "")))
                    or not MD5_RE.fullmatch(str(segment.get("md5", "")))
                ):
                    raise PublicationError(
                        f"Transport segment binding differs: {relative} #{number}"
                    )
                physical_name = str(segment["name"])
                if physical_name in physical_names[storage_record]:
                    raise PublicationError(
                        f"Transport physical filename collides: "
                        f"{storage_record}/{physical_name}"
                    )
                seen_segments.add(identity)
                physical_names[storage_record].add(physical_name)
                offset += size
            if offset != shard["compressed_bytes"]:
                raise PublicationError(
                    f"Transport segment byte total differs: {relative}"
                )
    for record in RECORDS:
        if len(physical_names[record]) > 100:
            raise PublicationError(
                f"Transport exceeds Zenodo's 100-file limit: {record}"
            )
    return document


def expected_zenodo_file_url(record_id: str, filename: str) -> str:
    quoted = urllib.parse.quote(filename, safe="._-")
    return f"https://zenodo.org/records/{record_id}/files/{quoted}?download=1"


def load_project_url_map(
    path: Path,
    manifest: Mapping[str, Any],
    metadata: Mapping[str, Any],
    transport: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    raw = load_json(path, label="Zenodo URL map")
    if raw.get("schema") != reconstruct_paper.URL_MAP_SCHEMA_V2:
        raise PublicationError("Public Zenodo URL map must use schema v2")
    try:
        urls = reconstruct_paper.load_url_map(path, manifest)
    except Exception as exc:
        raise PublicationError(f"Zenodo URL map is invalid: {exc}") from exc
    transport_archives = {
        archive["relative_local_path"]: (record, archive)
        for record in RECORDS
        for archive in transport["records"][record]["archives"]
    }
    for relative, row in urls.items():
        record, archive = transport_archives[relative]
        expected_segments = archive["segments"]
        if len(row["segments"]) != len(expected_segments):
            raise PublicationError(f"Transport segment count differs: {relative}")
        for observed, expected_segment in zip(row["segments"], expected_segments):
            storage_record = expected_segment.get("storage_record", record)
            record_id = metadata["records"][storage_record]["record_id"]
            expected_url = expected_zenodo_file_url(record_id, expected_segment["name"])
            if (
                observed["name"] != expected_segment["name"]
                or observed["size_bytes"] != expected_segment["size_bytes"]
                or observed["sha256"] != expected_segment["sha256"]
                or observed["url"] != expected_url
            ):
                raise PublicationError(
                    f"Zenodo transport URL/binding differs: {relative}/"
                    f"{expected_segment['name']}"
                )
    return urls


def human_size(size: int) -> str:
    return f"{size:,} bytes ({size / 1_000_000_000:.2f} GB; {size / 2**30:.2f} GiB)"


def manifest_shards(manifest: Mapping[str, Any]):
    for record in RECORDS:
        for shard in manifest["records"][record]["shards"]:
            yield record, shard


def transport_segments(transport: Mapping[str, Any]):
    for record in RECORDS:
        for archive in transport["records"][record]["archives"]:
            for segment in archive["segments"]:
                yield record, archive, segment


def render_readme(
    template: str,
    manifest: Mapping[str, Any],
    transport: Mapping[str, Any],
    metadata: Mapping[str, Any],
    urls: Mapping[str, Mapping[str, Any]],
    manifest_sha256: str,
    transport_manifest_sha256: str,
) -> str:
    archive_rows = []
    for record, shard in manifest_shards(manifest):
        relative = shard["relative_archive_path"]
        segment_count = len(urls[relative]["segments"])
        archive_rows.append(
            "| `{}` | `{}` | {} | `{}` | {} independently verified chunks |".format(
                record,
                shard["name"],
                human_size(shard["compressed_bytes"]),
                shard["sha256"],
                segment_count,
            )
        )
    download_bytes = sum(
        shard["compressed_bytes"] for _, shard in manifest_shards(manifest)
    )
    segment_count = sum(1 for _ in transport_segments(transport))
    transport_file_count = segment_count + 2
    tree_bytes = manifest["tree"]["payload_bytes"]
    replacements = {
        "@@ARCHIVE_TABLE@@": "\n".join(archive_rows),
        "@@ARCHIVE_COUNT@@": str(len(archive_rows)),
        "@@TRANSPORT_SEGMENT_COUNT@@": str(segment_count),
        "@@TRANSPORT_FILE_COUNT@@": str(transport_file_count),
        "@@DOWNLOAD_SIZE@@": human_size(download_bytes),
        "@@TREE_SIZE@@": human_size(tree_bytes),
        "@@COMBINED_SIZE@@": human_size(download_bytes + tree_bytes),
        "@@MANIFEST_SHA256@@": manifest_sha256,
        "@@TRANSPORT_MANIFEST_SHA256@@": transport_manifest_sha256,
        "@@TREE_SHA256@@": manifest["tree"]["sha256"],
        "@@REPOSITORY_URL@@": metadata["repository_url"],
        "@@RELEASE_TAG@@": metadata["release_tag"],
    }
    for record in RECORDS:
        upper = record.upper()
        for key in ("record_id", "record_url", "doi", "doi_url"):
            replacements[f"@@{upper}_{key.upper()}@@"] = metadata["records"][record][
                key
            ]
    rendered = template
    for token, value in replacements.items():
        if token not in rendered:
            raise PublicationError(f"README template omits required token {token}")
        rendered = rendered.replace(token, value)
    if PLACEHOLDER_RE.search(rendered):
        raise PublicationError("Rendered README retains an unresolved placeholder")
    return rendered


def _resolved_root(path: Path, *, label: str) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink():
        raise PublicationError(f"{label} may not be a symlink: {lexical}")
    resolved = lexical.resolve()
    if not resolved.is_dir():
        raise PublicationError(f"{label} is not a directory: {resolved}")
    return resolved


def _validate_license_hashes(root: Path) -> None:
    for relative, expected in LICENSE_SHA256.items():
        path = root / relative
        if not path.is_file() or path.is_symlink() or sha256_file(path) != expected:
            raise PublicationError(f"License file is absent or differs: {relative}")


def validate_repository(root_path: Path) -> dict[str, Any]:
    root = _resolved_root(root_path, label="repository root")
    required = {
        "README.md",
        "manifest.json",
        "transport_manifest.json",
        "zenodo-urls.json",
        "release-metadata.json",
        "release_plan.json",
        "LICENSE",
        "LICENSES/MIT.txt",
        "LICENSES/CC-BY-4.0.txt",
        "build_archives.py",
        "reconstruct_paper.py",
        "repository_tools.py",
        "validate_repository.py",
        ".github/workflows/validate-publication.yml",
    }
    for relative in sorted(required):
        path = root / relative
        if not path.is_file() or path.is_symlink():
            raise PublicationError(f"Required repository file is absent: {relative}")
    forbidden = {
        "README.template.md",
        "release-metadata.template.json",
        "zenodo-urls.template.json",
        "TEMPLATE_NOT_A_RELEASE.txt",
    }
    for relative in forbidden:
        if (root / relative).exists() or (root / relative).is_symlink():
            raise PublicationError(
                f"Template-only file leaked into final repository: {relative}"
            )
    _validate_license_hashes(root)
    metadata = load_metadata(root / "release-metadata.json")
    manifest = load_project_manifest(root / "manifest.json")
    transport = load_project_transport_manifest(
        root / "transport_manifest.json", root / "manifest.json", manifest
    )
    urls = load_project_url_map(
        root / "zenodo-urls.json", manifest, metadata, transport
    )
    readme = (root / "README.md").read_text(encoding="utf-8")
    documentation = readme
    if "python3 reproduce.py" in readme:
        if "(DETAILED_GUIDE.md)" not in readme or "(SETUP.md)" not in readme:
            raise PublicationError("Beginner README must link the detailed and setup guides")
        for relative in ("DETAILED_GUIDE.md", "SETUP.md", "reproduce.py",
                         "check_reproduction.py", "requirements-replay-compatibility.txt",
                         "reproduction-runtime.env"):
            path = root / relative
            if not path.is_file() or path.is_symlink():
                raise PublicationError(f"Required beginner-guide file is absent: {relative}")
        documentation += "\n" + (root / "DETAILED_GUIDE.md").read_text(encoding="utf-8")
    if PLACEHOLDER_RE.search(documentation):
        raise PublicationError("Published guides contain an unresolved publication placeholder")
    required_readme_values = {
        metadata["repository_url"],
        metadata["release_tag"],
        sha256_file(root / "manifest.json"),
        sha256_file(root / "transport_manifest.json"),
        manifest["tree"]["sha256"],
        "MIT",
        "CC BY 4.0",
        './reproduce_all_figures.sh --output "$PWD" --force',
    }
    for record in RECORDS:
        required_readme_values.update(metadata["records"][record].values())
    missing = sorted(value for value in required_readme_values if value not in documentation)
    if missing:
        raise PublicationError(f"Published guides omit required resolved values: {missing}")
    plan = load_json(root / "release_plan.json", label="release plan")
    if (
        plan.get("github", {}).get("repository_name") != "allulose-sucrose-hypothalamus"
        or plan.get("github", {}).get("release_tag") != RELEASE_TAG
        or plan.get("visibility", {}).get("github") != "public"
        or plan.get("licensing", {}).get("scheme")
        != "MIT-for-code-and-CC-BY-4.0-for-documents-and-data"
    ):
        raise PublicationError(
            "Copied release plan differs from the approved public plan"
        )
    return {
        "archive_count": 10,
        "transport_segment_count": sum(1 for _ in transport_segments(transport)),
        "transport_file_count": sum(1 for _ in transport_segments(transport)) + 2,
        "download_bytes": sum(
            shard["compressed_bytes"] for _, shard in manifest_shards(manifest)
        ),
        "tree_bytes": manifest["tree"]["payload_bytes"],
        "manifest_sha256": sha256_file(root / "manifest.json"),
        "tree_sha256": manifest["tree"]["sha256"],
    }


def assemble_repository(
    template_root_path: Path,
    manifest_path: Path,
    transport_manifest_path: Path,
    url_map_path: Path,
    metadata_path: Path,
    output_path: Path,
) -> Path:
    template_root = _resolved_root(template_root_path, label="template root")
    target_lexical = output_path.expanduser().absolute()
    if target_lexical.exists() or target_lexical.is_symlink():
        raise PublicationError(f"Output path must not exist: {target_lexical}")
    target = target_lexical.parent.resolve() / target_lexical.name
    if target == template_root or target.is_relative_to(template_root):
        raise PublicationError(
            "Final repository must be assembled outside the template"
        )
    metadata = load_metadata(metadata_path)
    manifest = load_project_manifest(manifest_path)
    transport = load_project_transport_manifest(
        transport_manifest_path, manifest_path, manifest
    )
    urls = load_project_url_map(url_map_path, manifest, metadata, transport)
    template = (template_root / "README.template.md").read_text(encoding="utf-8")
    manifest_digest = sha256_file(manifest_path)
    readme = render_readme(
        template,
        manifest,
        transport,
        metadata,
        urls,
        manifest_digest,
        sha256_file(transport_manifest_path),
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    stage = target.with_name(f".{target.name}.assembling.{os.getpid()}")
    if stage.exists() or stage.is_symlink():
        raise PublicationError(f"Assembly staging path already exists: {stage}")
    try:
        stage.mkdir(mode=0o700)
        (stage / "LICENSES").mkdir()
        (stage / ".github" / "workflows").mkdir(parents=True)
        copies = {
            "build_archives.py": "build_archives.py",
            "reconstruct_paper.py": "reconstruct_paper.py",
            "repository_tools.py": "repository_tools.py",
            "validate_repository.py": "validate_repository.py",
            "release_plan.json": "release_plan.json",
            "LICENSE": "LICENSE",
            "LICENSES/MIT.txt": "LICENSES/MIT.txt",
            "LICENSES/CC-BY-4.0.txt": "LICENSES/CC-BY-4.0.txt",
            "workflow/validate-publication.yml": ".github/workflows/validate-publication.yml",
        }
        for source_relative, destination_relative in copies.items():
            source = template_root / source_relative
            if not source.is_file() or source.is_symlink():
                raise PublicationError(f"Template asset is absent: {source_relative}")
            shutil.copyfile(source, stage / destination_relative)
        shutil.copyfile(manifest_path, stage / "manifest.json")
        shutil.copyfile(transport_manifest_path, stage / "transport_manifest.json")
        shutil.copyfile(url_map_path, stage / "zenodo-urls.json")
        shutil.copyfile(metadata_path, stage / "release-metadata.json")
        (stage / "README.md").write_text(readme, encoding="utf-8")
        for path in stage.rglob("*"):
            if path.is_file():
                os.chmod(path, 0o755 if path.suffix == ".py" else 0o644)
        validate_repository(stage)
        os.chmod(stage, 0o755)
        os.replace(stage, target)
    except Exception:
        if stage.exists() and not stage.is_symlink():
            shutil.rmtree(stage)
        raise
    return target
