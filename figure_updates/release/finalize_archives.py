#!/usr/bin/env python3
"""Finalize a staged archive set while safely reusing verified ZIP shards.

The Paper tree is inventoried again. Existing shards are reused only when
their prior manifest settings, checksums, and members agree with the new plan.
The staged manifest may be hard-linked to a prior generation, so the new
manifest is always written to a new inode and atomically replaced.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import stat
import sys
import time
from pathlib import Path
from typing import Any, Mapping, Sequence

import build_archives as archives


class FinalizationError(archives.ArchiveError):
    """A staged archive set cannot be finalized safely."""


def _overlap(first: Path, second: Path) -> bool:
    return (
        first == second or first.is_relative_to(second) or second.is_relative_to(first)
    )


def _existing_directory(path: Path, label: str) -> Path:
    lexical = path.expanduser().absolute()
    if lexical.is_symlink():
        raise FinalizationError(f"{label} may not be a symlink: {lexical}")
    resolved = lexical.resolve()
    if not resolved.is_dir():
        raise FinalizationError(f"{label} is missing or not a directory: {resolved}")
    return resolved


def _validate_positive(value: int, label: str, maximum: int | None = None) -> None:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise FinalizationError(f"{label} must be a positive integer")
    if maximum is not None and value > maximum:
        raise FinalizationError(f"{label} must not exceed {maximum}")


def _audit_stage(stage: Path, expected: set[Path]) -> None:
    allowed_directories = set(archives.RECORDS)
    for child in stage.iterdir():
        observed = child.lstat()
        if stat.S_ISLNK(observed.st_mode):
            raise FinalizationError(f"staged entry may not be a symlink: {child}")
        if stat.S_ISDIR(observed.st_mode):
            if child.name not in allowed_directories:
                raise FinalizationError(f"unexpected staged directory: {child}")
        elif not stat.S_ISREG(observed.st_mode) or child.name != "manifest.json":
            raise FinalizationError(f"unexpected staged root file: {child}")
    for record in archives.RECORDS:
        directory = stage / record
        if directory.exists() or directory.is_symlink():
            observed = directory.lstat()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISDIR(observed.st_mode):
                raise FinalizationError(f"invalid staged record directory: {directory}")
        else:
            directory.mkdir(mode=0o755)
        for child in directory.iterdir():
            observed = child.lstat()
            relative = child.relative_to(stage)
            if (
                stat.S_ISLNK(observed.st_mode)
                or not stat.S_ISREG(observed.st_mode)
                or relative not in expected
            ):
                raise FinalizationError(f"unexpected staged shard entry: {relative}")


def _load_prior_manifest(
    path: Path, trusted_legacy_sha256: str | None = None
) -> Mapping[str, Any] | None:
    if not path.exists():
        return None
    observed = path.lstat()
    if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
        raise FinalizationError("staged manifest must be a regular file")
    document: Any = None
    try:
        document = json.loads(path.read_bytes())
        archives.validate_manifest(document)
    except Exception as exc:
        if trusted_legacy_sha256 is not None:
            if not archives.SHA_RE.fullmatch(trusted_legacy_sha256):
                raise FinalizationError(
                    "trusted legacy manifest SHA-256 is invalid"
                ) from exc
            if archives.sha256_file(path) != trusted_legacy_sha256:
                raise FinalizationError(
                    "trusted legacy manifest SHA-256 differs"
                ) from exc
            records = document.get("records") if isinstance(document, Mapping) else None
            if (
                not isinstance(document, Mapping)
                or document.get("schema") != archives.SCHEMA
                or not isinstance(records, Mapping)
                or set(records) != set(archives.RECORDS)
            ):
                raise FinalizationError(
                    "trusted legacy manifest schema or record inventory differs"
                ) from exc
            print(
                "[PRIOR] full legacy tree validation failed, but explicitly "
                f"SHA-256-bound shard rows are eligible: {type(exc).__name__}: {exc}",
                flush=True,
            )
            return document
        print(
            f"[PRIOR] manifest cannot authorize reuse: {type(exc).__name__}: {exc}",
            flush=True,
        )
        return None
    return document


def _reusable_rows(
    document: Mapping[str, Any] | None,
    payload_limit: int,
    max_archive_bytes: int,
    compression_level: int,
) -> dict[str, Mapping[str, Any]]:
    if document is None:
        return {}
    if (
        document.get("compression")
        != {"method": "ZIP_DEFLATED", "level": compression_level}
        or document.get("limits", {}).get("uncompressed_payload_per_shard")
        != payload_limit
        or document.get("limits", {}).get("compressed_bytes_per_shard")
        != max_archive_bytes
    ):
        print(
            "[PRIOR] build settings differ; staged shards will be rebuilt", flush=True
        )
        return {}
    return {
        shard["relative_archive_path"]: shard
        for record in archives.RECORDS
        for shard in document["records"][record]["shards"]
    }


def _replace_manifest(path: Path, payload: bytes) -> None:
    """Create a new inode before replacing a possibly hard-linked manifest."""
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if temporary.exists() or temporary.is_symlink():
        raise FinalizationError(f"temporary manifest path must be absent: {temporary}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(temporary, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def finalize_archives(
    paper_root: Path,
    stage_path: Path,
    output_path: Path,
    *,
    contract_script: Path | None = None,
    payload_limit: int = archives.DEFAULT_SHARD_PAYLOAD_BYTES,
    max_archive_bytes: int = archives.DEFAULT_MAX_ARCHIVE_BYTES,
    compression_level: int = 6,
    jobs: int = 10,
    trusted_legacy_manifest_sha256: str | None = None,
) -> Path:
    """Finalize stage_path into absent output_path and return the output."""
    _validate_positive(jobs, "jobs")
    _validate_positive(
        payload_limit, "payload limit", archives.DEFAULT_SHARD_PAYLOAD_BYTES
    )
    _validate_positive(
        max_archive_bytes, "archive-byte limit", archives.DEFAULT_MAX_ARCHIVE_BYTES
    )
    _validate_positive(compression_level, "compression level", 9)

    paper = _existing_directory(paper_root, "Paper root")
    stage = _existing_directory(stage_path, "staged archive set")
    output_lexical = output_path.expanduser().absolute()
    if output_lexical.exists() or output_lexical.is_symlink():
        raise FinalizationError(f"final output must be absent: {output_lexical}")
    output = output_lexical.parent.resolve() / output_lexical.name
    if stage.parent != output.parent:
        raise FinalizationError(
            "stage and output must share a parent for atomic rename"
        )
    if _overlap(paper, stage) or _overlap(paper, output):
        raise FinalizationError("Paper tree and archive-set paths must not overlap")
    contract = (
        contract_script.expanduser().resolve()
        if contract_script is not None
        else paper / "scripts" / "utilities" / "01_validate_bundle.py"
    )

    started = time.monotonic()
    print("[START] Re-inventorying the Paper tree", flush=True)
    inventory, contract_hash = archives.inventory_tree(paper, contract_script=contract)
    plans = archives.plan_shards(inventory.files, payload_limit=payload_limit)
    print(
        f"[INVENTORY] {len(inventory.files)} files, "
        f"{sum(item.size for item in inventory.files)} bytes, "
        f"{len(plans)} shards in {time.monotonic() - started:.1f}s",
        flush=True,
    )
    expected = {Path(plan.record) / plan.name for plan in plans}
    _audit_stage(stage, expected)
    reusable = _reusable_rows(
        _load_prior_manifest(
            stage / "manifest.json",
            trusted_legacy_sha256=trusted_legacy_manifest_sha256,
        ),
        payload_limit,
        max_archive_bytes,
        compression_level,
    )

    quarantine = stage.parent / f".{output.name}.incomplete-quarantine.{os.getpid()}"
    if quarantine.exists() or quarantine.is_symlink():
        raise FinalizationError(f"quarantine path must be absent: {quarantine}")
    quarantine.mkdir(mode=0o700)

    def ensure_shard(
        plan: archives.Shard,
    ) -> tuple[str, dict[str, Any], str, bool]:
        target = stage / plan.record / plan.name
        relative = f"{plan.record}/{plan.name}"
        if target.exists() or target.is_symlink():
            observed = target.lstat()
            if stat.S_ISLNK(observed.st_mode) or not stat.S_ISREG(observed.st_mode):
                raise FinalizationError(f"staged shard must be regular: {relative}")
            try:
                prior_row = reusable.get(relative)
                if prior_row is None:
                    raise FinalizationError("no trusted prior-manifest row")
                digest = archives.sha256_file(target)
                if (
                    prior_row.get("compressed_bytes") != observed.st_size
                    or prior_row.get("sha256") != digest
                ):
                    raise FinalizationError("prior-manifest size or SHA-256 differs")
                archives.verify_shard(target, plan)
                return (
                    plan.name,
                    {
                        "name": target.name,
                        "relative_archive_path": relative,
                        "member_count": len(plan.files),
                        "uncompressed_bytes": plan.payload_bytes,
                        "compressed_bytes": observed.st_size,
                        "sha256": digest,
                    },
                    "reused",
                    False,
                )
            except Exception as exc:
                destination = quarantine / plan.record / plan.name
                destination.parent.mkdir(parents=True, exist_ok=True)
                os.replace(target, destination)
                print(
                    f"[REPAIR] {relative}: {type(exc).__name__}: {exc}",
                    flush=True,
                )
        row = archives.write_shard(
            plan,
            target,
            compression_level=compression_level,
            max_archive_bytes=max_archive_bytes,
        )
        return plan.name, row, "rebuilt", True

    rows: dict[str, dict[str, Any]] = {}
    repairs = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(ensure_shard, plan): plan for plan in plans}
        for future in concurrent.futures.as_completed(futures):
            plan = futures[future]
            name, row, disposition, repaired = future.result()
            rows[name] = row
            repairs += int(repaired)
            print(
                f"[SHARD] {disposition}: {plan.record}/{name} "
                f"({row['compressed_bytes']} bytes, sha256={row['sha256']})",
                flush=True,
            )

    manifest = archives.make_manifest(
        inventory,
        plans,
        rows,
        contract_hash=contract_hash,
        payload_limit=payload_limit,
        max_archive_bytes=max_archive_bytes,
        compression_level=compression_level,
    )
    archives.validate_manifest(manifest)
    _replace_manifest(stage / "manifest.json", archives.canonical_json(manifest))
    os.chmod(stage, 0o755)
    os.replace(stage, output)
    print(f"[PASS] Final archive set: {output}", flush=True)
    if repairs:
        size = sum(
            path.stat().st_size for path in quarantine.glob("*/*.zip") if path.is_file()
        )
        print(f"[QUARANTINE] {quarantine} ({size} bytes)", flush=True)
    else:
        quarantine.rmdir()
        print("[QUARANTINE] none", flush=True)
    print(f"[ELAPSED] {time.monotonic() - started:.1f}s", flush=True)
    return output


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, required=True)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--contract-script", type=Path)
    parser.add_argument(
        "--trusted-legacy-manifest-sha256",
        help=(
            "Explicitly authorize shard rows from an otherwise legacy-invalid "
            "staged manifest with this exact SHA-256. Every shard still receives "
            "full current-plan verification before reuse."
        ),
    )
    parser.add_argument(
        "--compression-level", type=int, choices=range(1, 10), default=6
    )
    parser.add_argument("--jobs", type=int, default=10)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    finalize_archives(
        args.paper_root,
        args.stage,
        args.output,
        contract_script=args.contract_script,
        compression_level=args.compression_level,
        jobs=args.jobs,
        trusted_legacy_manifest_sha256=args.trusted_legacy_manifest_sha256,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (archives.ArchiveError, OSError) as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        raise SystemExit(2)
