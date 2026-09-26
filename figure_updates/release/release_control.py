#!/usr/bin/env python3
"""Validate and prepare the public release without contacting remote services.

This module deliberately contains no HTTP or GitHub client.  It binds release
metadata to the reviewed manuscript author audit, verifies every archive and
ZIP member, evaluates the required license/rights/credential gates, and can
write a credential-free request document for a separately reviewed executor.
Every remote-action command fails closed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import stat
import sys
import zipfile
from datetime import date
from pathlib import Path, PurePosixPath
from typing import Any, Mapping, MutableMapping, Sequence


HERE = Path(__file__).resolve().parent
PAPER = HERE.parent
sys.path.insert(0, str(HERE))

import build_archives as archives  # noqa: E402


PLAN_SCHEMA = "apotome-public-release-plan-v1"
APPROVALS_SCHEMA = "apotome-public-release-approvals-v1"
REQUEST_SCHEMA = "apotome-public-release-request-v1"
EXPECTED_AUTHOR_AUDIT = "manuscript_sources/PAPER_AUTHOR_PUBMED_AUDIT_20260830.json"
EXPECTED_ZENODO_NAMES = (
    "Segura, Nancy",
    "Azócar, Vinka",
    "Aguilera, Claudia",
    "Sanhueza, Constanza",
    "Aravena, Bastían",
    "Casas, Rocío",
    "Magdalena, Rocío",
    "Troncoso, Isidora",
    "Quiroz, Aracelly",
    "Mennickent Barros, Daniela",
    "Flores, Ricardo",
    "Espinoza, Francisca",
    "Recabal-Beyer, Antonia",
    "Elizondo, Roberto",
    "Dorfmann, Mauricio D",
    "García, María de los Ángeles",
    "González-Pecchi, Valentina",
    "Farkas, Carlos",
)
REQUIRED_VISIBILITY = {
    "zenodo": {"record": "public", "files": "public", "access_right": "open"},
    "github": "public",
}
EXPECTED_LICENSING = {
    "scheme": "MIT-for-code-and-CC-BY-4.0-for-documents-and-data",
    "root_notice": "LICENSE",
    "zenodo_license_ids": {
        "software_core": "other-open",
        "main_figure_data": "cc-by-4.0",
        "supplementary_figure_data": "cc-by-4.0",
        "figs3_wsi": "cc-by-4.0",
    },
    "software_core_license_matrix": {
        "code": "MIT",
        "documents_and_data": "CC-BY-4.0",
    },
}
PLACEHOLDER_VALUES = {
    "",
    "none",
    "null",
    "required",
    "tbd",
    "todo",
    "unknown",
    "choose-a-license",
}
TOKEN_ENVIRONMENT = {
    "zenodo_sandbox": ("ZENODO_SANDBOX_TOKEN",),
    "zenodo_production": ("ZENODO_TOKEN",),
    "github": ("GH_TOKEN", "GITHUB_TOKEN"),
}
OPERATION_POLICY = {
    "zenodo-sandbox-draft": {
        "credential": "zenodo_sandbox",
        "record_ids": None,
        "github_owner": False,
    },
    "zenodo-sandbox-upload": {
        "credential": "zenodo_sandbox",
        "record_ids": "zenodo_sandbox_record_ids",
        "github_owner": False,
    },
    "zenodo-production-draft": {
        "credential": "zenodo_production",
        "record_ids": None,
        "github_owner": False,
    },
    "zenodo-production-upload": {
        "credential": "zenodo_production",
        "record_ids": "zenodo_production_record_ids",
        "github_owner": False,
    },
    "github-public-draft": {
        "credential": "github",
        "record_ids": None,
        "github_owner": True,
    },
    "github-public-release": {
        "credential": "github",
        "record_ids": None,
        "github_owner": True,
    },
}
QUOTA_REQUIRED_OPERATIONS = frozenset(
    {
        "zenodo-sandbox-upload",
        "zenodo-production-upload",
        "github-public-release",
    }
)
CHUNK = 8 * 1024 * 1024


class ReleaseControlError(RuntimeError):
    """A local release contract failed closed."""


def canonical_json(document: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def regular_file(path: Path, label: str) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_file():
        raise ReleaseControlError(f"Missing or linked {label}: {lexical}")
    return lexical.resolve()


def load_json(path: Path, label: str) -> dict[str, Any]:
    source = regular_file(path, label)
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise ReleaseControlError(f"Invalid {label}: {source}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ReleaseControlError(f"{label} must be one JSON object: {source}")
    return payload


def valid_zenodo_name(value: Any) -> bool:
    """Accept only canonical ``Family name, Given names`` creator strings."""
    if not isinstance(value, str) or value.count(",") != 1:
        return False
    family, given = value.split(",", 1)
    if not given.startswith(" "):
        return False
    given = given[1:]
    if not family or not given:
        return False
    if family != " ".join(family.split()) or given != " ".join(given.split()):
        return False
    return not any(ord(character) < 32 or ord(character) == 127 for character in value)


def validate_release_plan(
    plan: Mapping[str, Any], audit: Mapping[str, Any], audit_path: Path
) -> dict[str, Any]:
    """Bind the plan exactly to the author-approved audit and public policy."""
    if plan.get("schema") != PLAN_SCHEMA:
        raise ReleaseControlError("Release-plan schema differs")
    if audit.get("status") != "reviewed":
        raise ReleaseControlError("Author audit is not reviewed")
    source = plan.get("source_author_audit")
    if not isinstance(source, dict):
        raise ReleaseControlError("Release plan has no source-author binding")
    if source.get("path") != EXPECTED_AUTHOR_AUDIT:
        raise ReleaseControlError("Release plan points to the wrong author audit")
    observed_audit_hash = sha256_file(regular_file(audit_path, "author audit"))
    if source.get("sha256") != observed_audit_hash:
        raise ReleaseControlError("Release-plan author-audit SHA-256 differs")

    manuscript = plan.get("manuscript")
    if not isinstance(manuscript, dict):
        raise ReleaseControlError("Release plan has no manuscript metadata")
    if manuscript.get("title") != audit.get("paper_title"):
        raise ReleaseControlError("Manuscript title differs from the reviewed audit")
    audit_authors = audit.get("authors")
    plan_authors = manuscript.get("authors")
    if not isinstance(audit_authors, list) or not isinstance(plan_authors, list):
        raise ReleaseControlError("Author metadata is not a list")
    exact_audit = [
        (row.get("order"), row.get("display"))
        for row in audit_authors
        if isinstance(row, dict)
    ]
    exact_plan = [
        (row.get("order"), row.get("display"))
        for row in plan_authors
        if isinstance(row, dict)
    ]
    if exact_plan != exact_audit or [row[0] for row in exact_plan] != list(
        range(1, len(exact_audit) + 1)
    ):
        raise ReleaseControlError("Author display list/order differs from the audit")
    for row in plan_authors:
        display = str(row.get("display", ""))
        corresponding = display.endswith("*")
        if (
            row.get("service_name") != display.removesuffix("*")
            or row.get("corresponding") is not corresponding
        ):
            raise ReleaseControlError(
                f"Author service-name/correspondence mapping differs: {display}"
            )
    zenodo_names = [row.get("zenodo_name") for row in plan_authors]
    if any(not valid_zenodo_name(name) for name in zenodo_names):
        raise ReleaseControlError(
            "Every Zenodo creator name must use canonical Family name, Given names form"
        )
    if zenodo_names != list(EXPECTED_ZENODO_NAMES):
        raise ReleaseControlError(
            "Zenodo creator names/order differ from the author-confirmed list"
        )
    if sum(bool(row.get("corresponding")) for row in plan_authors) != 1:
        raise ReleaseControlError(
            "Exactly one audited corresponding author is required"
        )

    if plan.get("visibility") != REQUIRED_VISIBILITY:
        raise ReleaseControlError("Release visibility must be public/open everywhere")
    if plan.get("licensing") != EXPECTED_LICENSING:
        raise ReleaseControlError("Approved mixed-license matrix differs")
    license_notice = regular_file(PAPER / "LICENSE", "root license notice")
    if license_notice.name != EXPECTED_LICENSING["root_notice"]:
        raise ReleaseControlError("Root license notice differs")
    records = plan.get("records")
    if not isinstance(records, dict) or set(records) != set(archives.RECORDS):
        raise ReleaseControlError("Release plan must contain exactly four records")
    for name in archives.RECORDS:
        row = records[name]
        if (
            not isinstance(row, dict)
            or row.get("archive_record") != name
            or row.get("license_input") != name
            or row.get("resource_type") not in {"software", "dataset"}
            or not str(row.get("title", "")).startswith(str(manuscript["title"]) + ": ")
            or not str(row.get("description", "")).strip()
        ):
            raise ReleaseControlError(f"Invalid logical-record metadata: {name}")
    if records["software_core"]["resource_type"] != "software" or any(
        records[name]["resource_type"] != "dataset"
        for name in archives.RECORDS
        if name != "software_core"
    ):
        raise ReleaseControlError("Four-record resource types differ")

    github = plan.get("github")
    if (
        not isinstance(github, dict)
        or github.get("visibility") != "public"
        or github.get("release_is_draft_until_clean_room_passes") is not True
        or not re.fullmatch(r"[A-Za-z0-9._-]+", str(github.get("repository_name", "")))
        or not re.fullmatch(
            r"v[0-9]+\.[0-9]+\.[0-9]+", str(github.get("release_tag", ""))
        )
    ):
        raise ReleaseControlError("GitHub public-draft policy differs")

    required = plan.get("required_inputs")
    if not isinstance(required, dict):
        raise ReleaseControlError("Required-input declaration is absent")
    if set(required.get("licenses", [])) != set(archives.RECORDS):
        raise ReleaseControlError("License gates do not cover all four records")
    if set(required.get("quota_confirmations", [])) != (
        set(archives.RECORDS) - {"software_core"}
    ):
        raise ReleaseControlError("Quota gates do not cover all three data records")
    rights = required.get("rights_confirmations")
    if (
        not isinstance(rights, list)
        or len(rights) < 6
        or len(rights) != len(set(rights))
    ):
        raise ReleaseControlError("Rights gates are incomplete or duplicated")
    if required.get("credential_environment") != {
        key: list(value) for key, value in TOKEN_ENVIRONMENT.items()
    }:
        raise ReleaseControlError("Credential environment declaration differs")
    policy = plan.get("remote_execution_policy")
    if not isinstance(policy, dict) or any(
        policy.get(key) is not False
        for key in (
            "network_enabled",
            "draft_creation_enabled",
            "upload_enabled",
            "publication_enabled",
        )
    ):
        raise ReleaseControlError("Local-only remote-execution policy differs")
    return {
        "status": "PASS",
        "title": manuscript["title"],
        "authors": len(plan_authors),
        "records": list(archives.RECORDS),
        "visibility": plan["visibility"],
        "author_audit_sha256": observed_audit_hash,
    }


def load_validated_plan(
    plan_path: Path, audit_path: Path
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = load_json(plan_path, "release plan")
    audit = load_json(audit_path, "author audit")
    return plan, validate_release_plan(plan, audit, audit_path)


def _contained(root: Path, relative: str) -> Path:
    archives.safe_relative(relative)
    target = root.joinpath(*PurePosixPath(relative).parts)
    resolved = target.resolve(strict=False)
    if resolved == root or not resolved.is_relative_to(root):
        raise ReleaseControlError(f"Archive path escapes release root: {relative}")
    return target


def verify_archive_set(
    archive_root: Path, *, contract_script: Path | None = None
) -> dict[str, Any]:
    """Verify archive files and every decompressed member against manifest SHA-256."""
    lexical = archive_root.expanduser().absolute()
    if lexical.is_symlink() or not lexical.is_dir():
        raise ReleaseControlError(f"Archive root is missing or linked: {lexical}")
    root = lexical.resolve()
    manifest_path = regular_file(root / "manifest.json", "archive manifest")
    manifest = load_json(manifest_path, "archive manifest")
    try:
        archives.validate_manifest(manifest)
    except archives.ArchiveError as exc:
        raise ReleaseControlError(f"Archive manifest failed: {exc}") from exc

    contract_bound = False
    recorded_contract = manifest.get("canonical_exclusion_contract", {}).get(
        "script_sha256"
    )
    if contract_script is not None:
        contract = regular_file(contract_script, "canonical bundle contract")
        if recorded_contract != sha256_file(contract):
            raise ReleaseControlError("Archive exclusion-contract SHA-256 differs")
        contract_bound = True
    elif isinstance(recorded_contract, str) and archives.SHA_RE.fullmatch(
        recorded_contract
    ):
        contract_bound = True

    shard_to_rows: dict[str, list[Mapping[str, Any]]] = {}
    for row in manifest["tree"]["files"]:
        shard_to_rows.setdefault(str(row["shard"]), []).append(row)
    expected_files = {"manifest.json"}
    records_summary: dict[str, Any] = {}
    verified_members = 0
    for record in archives.RECORDS:
        record_row = manifest["records"][record]
        summaries = []
        for shard in record_row["shards"]:
            relative = str(shard["relative_archive_path"])
            path = _contained(root, relative)
            expected_files.add(relative)
            if path.is_symlink() or not path.is_file():
                raise ReleaseControlError(f"Missing or linked ZIP shard: {relative}")
            if path.stat().st_size != shard["compressed_bytes"]:
                raise ReleaseControlError(f"ZIP shard size differs: {relative}")
            observed_sha = sha256_file(path)
            if observed_sha != shard["sha256"]:
                raise ReleaseControlError(f"ZIP shard SHA-256 differs: {relative}")
            rows = shard_to_rows.get(str(shard["name"]), [])
            expected_members = [str(row["member"]) for row in rows]
            try:
                with zipfile.ZipFile(path, "r") as zipped:
                    infos = zipped.infolist()
                    if [info.filename for info in infos] != expected_members or len(
                        infos
                    ) != len(set(expected_members)):
                        raise ReleaseControlError(
                            f"ZIP member inventory/order differs: {relative}"
                        )
                    for info, row in zip(infos, rows):
                        if (
                            info.date_time != archives.ZIP_DATETIME
                            or info.compress_type != zipfile.ZIP_DEFLATED
                            or info.file_size != row["size_bytes"]
                            or ((info.external_attr >> 16) & 0o7777)
                            != int(str(row["mode"]), 8)
                        ):
                            raise ReleaseControlError(
                                f"ZIP member metadata differs: {info.filename}"
                            )
                        digest = hashlib.sha256()
                        size = 0
                        with zipped.open(info, "r") as handle:
                            for block in iter(lambda: handle.read(CHUNK), b""):
                                digest.update(block)
                                size += len(block)
                        if (
                            size != row["size_bytes"]
                            or digest.hexdigest() != row["sha256"]
                        ):
                            raise ReleaseControlError(
                                f"ZIP member checksum differs: {info.filename}"
                            )
                        verified_members += 1
            except zipfile.BadZipFile as exc:
                raise ReleaseControlError(
                    f"Invalid ZIP shard: {relative}: {exc}"
                ) from exc
            summaries.append(
                {
                    "name": shard["name"],
                    "size_bytes": shard["compressed_bytes"],
                    "sha256": observed_sha,
                    "member_count": len(rows),
                }
            )
        records_summary[record] = {
            "compressed_archive_bytes": record_row["compressed_archive_bytes"],
            "shard_count": record_row["shard_count"],
            "shards": summaries,
        }

    actual_files: set[str] = set()
    for path in root.rglob("*"):
        if path.is_symlink():
            raise ReleaseControlError(f"Linked archive-set entry: {path}")
        if path.is_file():
            actual_files.add(path.relative_to(root).as_posix())
        elif not path.is_dir():
            raise ReleaseControlError(f"Non-regular archive-set entry: {path}")
    if actual_files != expected_files:
        raise ReleaseControlError(
            "Archive-set closed-world inventory differs: "
            f"missing={sorted(expected_files-actual_files)} "
            f"extra={sorted(actual_files-expected_files)}"
        )
    if verified_members != manifest["tree"]["file_count"]:
        raise ReleaseControlError("Verified member count differs from public tree")
    return {
        "status": "PASS",
        "archive_root": str(root),
        "manifest_sha256": sha256_file(manifest_path),
        "manifest_size_bytes": manifest_path.stat().st_size,
        "tree_sha256": manifest["tree"]["sha256"],
        "tree_files": manifest["tree"]["file_count"],
        "verified_zip_members": verified_members,
        "canonical_contract_bound": contract_bound,
        "records": records_summary,
    }


def _valid_license(value: Any) -> bool:
    if not isinstance(value, str) or value.strip().casefold() in PLACEHOLDER_VALUES:
        return False
    return bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.+-]{1,127}", value.strip()))


def _valid_release_date(value: Any) -> bool:
    if not isinstance(value, str):
        return False
    try:
        return date.fromisoformat(value).isoformat() == value
    except ValueError:
        return False


def evaluate_readiness(
    plan: Mapping[str, Any],
    approvals: Mapping[str, Any],
    archive_summary: Mapping[str, Any],
    operation: str,
    *,
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Return blockers without ever exposing credential values."""
    if operation not in OPERATION_POLICY:
        raise ReleaseControlError(f"Unknown remote operation: {operation}")
    blockers: list[str] = []
    if approvals.get("schema") != APPROVALS_SCHEMA:
        blockers.append("approvals schema is missing or invalid")
    licenses = approvals.get("licenses", {})
    planned_licenses = plan.get("licensing", {}).get("zenodo_license_ids", {})
    for record in archives.RECORDS:
        if not isinstance(licenses, dict) or not _valid_license(licenses.get(record)):
            blockers.append(f"license not approved: {record}")
        elif licenses.get(record) != planned_licenses.get(record):
            blockers.append(f"license differs from approved plan: {record}")
    rights = approvals.get("rights_confirmations", {})
    for key in plan["required_inputs"]["rights_confirmations"]:
        if not isinstance(rights, dict) or rights.get(key) is not True:
            blockers.append(f"rights confirmation is false: {key}")
    if approvals.get("public_visibility_acknowledged") is not True:
        blockers.append("public visibility is not acknowledged")
    if approvals.get("clean_room_reproduction_passed") is not True:
        blockers.append("clean-room reproduction is not confirmed")
    if not _valid_release_date(approvals.get("release_date")):
        blockers.append("release_date is missing or not ISO YYYY-MM-DD")
    if archive_summary.get("status") != "PASS":
        blockers.append("archive set is not verified")
    if archive_summary.get("canonical_contract_bound") is not True:
        blockers.append("archive set is not bound to the canonical exclusion contract")

    if operation in QUOTA_REQUIRED_OPERATIONS:
        quotas = approvals.get("quota_confirmations", {})
        for record, row in archive_summary.get("records", {}).items():
            if int(
                row.get("compressed_archive_bytes", 0)
            ) > archives.ZENODO_DEFAULT_RECORD_BYTES and (
                not isinstance(quotas, dict) or quotas.get(record) is not True
            ):
                blockers.append(f"Zenodo >50 GB quota not confirmed: {record}")

    policy = OPERATION_POLICY[operation]
    credential_group = str(policy["credential"])
    scopes = approvals.get("credential_scope_confirmations", {})
    if not isinstance(scopes, dict) or scopes.get(credential_group) is not True:
        blockers.append(f"credential scopes not confirmed: {credential_group}")
    environment = os.environ if environ is None else environ
    credential_names = TOKEN_ENVIRONMENT[credential_group]
    credential_present = any(
        bool(environment.get(name, "").strip()) for name in credential_names
    )
    if not credential_present:
        blockers.append(
            "credential is absent; set one of: " + ", ".join(credential_names)
        )

    targets = approvals.get("remote_targets", {})
    if not isinstance(targets, dict):
        targets = {}
    if policy["github_owner"]:
        owner = targets.get("github_owner")
        if not isinstance(owner, str) or not re.fullmatch(
            r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,38})", owner
        ):
            blockers.append("GitHub owner is missing or invalid")
        if targets.get("github_repository_name") != plan["github"]["repository_name"]:
            blockers.append("GitHub repository name differs from release plan")
    record_id_key = policy["record_ids"]
    if record_id_key:
        identifiers = targets.get(record_id_key)
        if not isinstance(identifiers, dict):
            identifiers = {}
        for record in archives.RECORDS:
            value = identifiers.get(record)
            if not isinstance(value, str) or not re.fullmatch(r"[0-9]+", value):
                blockers.append(f"remote record ID is missing: {record}")

    return {
        "status": "PASS" if not blockers else "BLOCKED",
        "operation": operation,
        "public_visibility": plan["visibility"],
        "blockers": blockers,
        "credential_environment_names": list(credential_names),
        "credential_present": credential_present,
        "credential_values_embedded": False,
        "remote_execution_enabled": False,
    }


def prepare_request_document(
    plan: Mapping[str, Any],
    approvals: Mapping[str, Any],
    archive_summary: Mapping[str, Any],
    operation: str,
    *,
    environ: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    readiness = evaluate_readiness(
        plan, approvals, archive_summary, operation, environ=environ
    )
    if readiness["status"] != "PASS":
        raise ReleaseControlError(
            "Release request is blocked: " + "; ".join(readiness["blockers"])
        )
    creators = [
        {
            "name": row["zenodo_name"],
            "service_name": row["service_name"],
            "manuscript_display": row["display"],
            "order": row["order"],
            "corresponding": row["corresponding"],
        }
        for row in plan["manuscript"]["authors"]
    ]
    zenodo_records: dict[str, Any] = {}
    for record in archives.RECORDS:
        metadata = plan["records"][record]
        zenodo_records[record] = {
            "metadata": {
                "title": metadata["title"],
                "description": metadata["description"],
                "resource_type": metadata["resource_type"],
                "publication_date": approvals["release_date"],
                "license": approvals["licenses"][record],
                "creators": creators,
            },
            "access": dict(plan["visibility"]["zenodo"]),
            "archive_shards": list(archive_summary["records"][record]["shards"]),
            "ancillary_files": (
                [
                    {
                        "name": "manifest.json",
                        "size_bytes": archive_summary["manifest_size_bytes"],
                        "sha256": archive_summary["manifest_sha256"],
                    }
                ]
                if record == "software_core"
                else []
            ),
            "related_logical_records": [
                other for other in archives.RECORDS if other != record
            ],
            "article_doi": approvals.get("article_doi"),
        }
    document = {
        "schema": REQUEST_SCHEMA,
        "operation": operation,
        "status": "READY_FOR_SEPARATELY_REVIEWED_EXECUTOR",
        "local_only": True,
        "remote_execution_enabled": False,
        "credentials_embedded": False,
        "credential_environment_names": readiness["credential_environment_names"],
        "manuscript": {
            "title": plan["manuscript"]["title"],
            "authors": creators,
        },
        "visibility": plan["visibility"],
        "zenodo_records": zenodo_records,
        "github": {
            "owner": approvals["remote_targets"].get("github_owner"),
            "repository_name": plan["github"]["repository_name"],
            "visibility": "public",
            "default_branch": plan["github"]["default_branch"],
            "release_tag": plan["github"]["release_tag"],
            "release_draft": True,
        },
        "archive_manifest": {
            "sha256": archive_summary["manifest_sha256"],
            "tree_sha256": archive_summary["tree_sha256"],
            "tree_files": archive_summary["tree_files"],
        },
        "approvals": {
            "rights_confirmations": approvals["rights_confirmations"],
            "public_visibility_acknowledged": approvals[
                "public_visibility_acknowledged"
            ],
            "clean_room_reproduction_passed": approvals[
                "clean_room_reproduction_passed"
            ],
            "quota_confirmations": approvals["quota_confirmations"],
        },
    }
    serialized = canonical_json(document)
    for name in sum((list(names) for names in TOKEN_ENVIRONMENT.values()), []):
        value = (os.environ if environ is None else environ).get(name, "")
        if value and value.encode("utf-8") in serialized:
            raise ReleaseControlError("Credential value leaked into request document")
    return document


def write_request(path: Path, document: Mapping[str, Any]) -> Path:
    target = path.expanduser().absolute()
    if target.exists() or target.is_symlink():
        raise ReleaseControlError(f"Request output must be absent: {target}")
    if not target.parent.is_dir() or target.parent.is_symlink():
        raise ReleaseControlError(
            f"Request output parent must be an existing regular directory: {target.parent}"
        )
    with target.open("xb") as handle:
        handle.write(canonical_json(document))
        handle.flush()
        os.fsync(handle.fileno())
    os.chmod(target, 0o600)
    return target


def remote_action(operation: str) -> None:
    if operation not in OPERATION_POLICY:
        raise ReleaseControlError(f"Unknown remote operation: {operation}")
    raise ReleaseControlError(
        "Remote draft/upload execution is intentionally disabled in this local-only "
        "tool. Review a prepared request and use a separately authorized executor."
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--plan", type=Path, default=HERE / "release_plan.json")
    result.add_argument(
        "--author-audit", type=Path, default=archives.local_input_path(PAPER, EXPECTED_AUTHOR_AUDIT)
    )
    sub = result.add_subparsers(dest="command", required=True)
    sub.add_parser("validate-plan")
    verify = sub.add_parser("verify-archives")
    verify.add_argument("--archive-root", type=Path, required=True)
    verify.add_argument(
        "--contract-script",
        type=Path,
        default=PAPER / "scripts/utilities/01_validate_bundle.py",
    )
    for name in ("readiness", "prepare-request"):
        item = sub.add_parser(name)
        item.add_argument("--archive-root", type=Path, required=True)
        item.add_argument("--approvals", type=Path, required=True)
        item.add_argument("--operation", choices=tuple(OPERATION_POLICY), required=True)
        item.add_argument(
            "--contract-script",
            type=Path,
            default=PAPER / "scripts/utilities/01_validate_bundle.py",
        )
        if name == "prepare-request":
            item.add_argument("--output", type=Path, required=True)
    remote = sub.add_parser("remote-action")
    remote.add_argument("--operation", choices=tuple(OPERATION_POLICY), required=True)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        plan, plan_summary = load_validated_plan(args.plan, args.author_audit)
        if args.command == "validate-plan":
            print(json.dumps(plan_summary, ensure_ascii=False, indent=2))
            return 0
        if args.command == "remote-action":
            remote_action(args.operation)
        archive_summary = verify_archive_set(
            args.archive_root, contract_script=args.contract_script
        )
        if args.command == "verify-archives":
            print(json.dumps(archive_summary, ensure_ascii=False, indent=2))
            return 0
        approvals = load_json(args.approvals, "release approvals")
        readiness = evaluate_readiness(plan, approvals, archive_summary, args.operation)
        if args.command == "readiness":
            print(json.dumps(readiness, ensure_ascii=False, indent=2))
            return 0 if readiness["status"] == "PASS" else 2
        document = prepare_request_document(
            plan, approvals, archive_summary, args.operation
        )
        target = write_request(args.output, document)
        print(f"[PASS] local request document: {target}")
        print("[LOCAL-ONLY] no Zenodo or GitHub request was sent")
        return 0
    except (ReleaseControlError, archives.ArchiveError) as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
