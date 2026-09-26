#!/usr/bin/env python3
"""Transactionally remove review/private metadata from selected Office files.

DOCX processing removes Word comment and people parts, their relationships and
content-type overrides, comment range/reference markers, comment-only styles,
and identities stored only in those parts. XLSX processing removes the
OneDrive absolute-path extension, Office revision identifiers/UUID attributes,
and core lastModifiedBy metadata. All other package members are retained.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import stat
import sys
import zipfile
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Sequence

from lxml import etree


COMMENT_PART_RE = re.compile(
    r"^word/(?:comments(?:extended|extensible|ids)?[0-9]*|people)\.xml$",
    re.IGNORECASE,
)
EMAIL_RE = re.compile(r"[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}", re.IGNORECASE)
COMMENT_MARKERS = {"commentRangeStart", "commentRangeEnd", "commentReference"}
DELETED_REVISION_ELEMENTS = {"del", "moveFrom"}
PRESERVED_REVISION_WRAPPERS = {"ins", "moveTo"}
REVISION_MARKERS = {
    "delText",
    "moveFromRangeStart",
    "moveFromRangeEnd",
    "moveToRangeStart",
    "moveToRangeEnd",
    "customXmlInsRangeStart",
    "customXmlInsRangeEnd",
    "customXmlDelRangeStart",
    "customXmlDelRangeEnd",
    "numberingChange",
}
COMMENT_STYLE_RE = re.compile(
    r"^Comment(?:Text|Reference|Subject|TextChar|SubjectChar)[0-9]*$",
    re.IGNORECASE,
)
REVISION_ATTRIBUTE_NAMES = {
    "uid",
    "uidLastSave",
    "coauthVersionLast",
    "coauthVersionMax",
}
CHUNK = 8 * 1024 * 1024


class SanitizationError(RuntimeError):
    """An Office package cannot be sanitized without violating the contract."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def local_name(value: str) -> str:
    return value.rsplit("}", 1)[-1]


def namespace(value: str) -> str:
    return value[1:].split("}", 1)[0] if value.startswith("{") else ""


def parse_xml(payload: bytes, label: str) -> etree._Element:
    parser = etree.XMLParser(
        resolve_entities=False,
        no_network=True,
        recover=False,
        remove_blank_text=False,
        huge_tree=True,
    )
    try:
        return etree.fromstring(payload, parser=parser)
    except (etree.XMLSyntaxError, ValueError) as exc:
        raise SanitizationError(f"invalid XML in {label}: {exc}") from exc


def serialize_xml(root: etree._Element, original: bytes) -> bytes:
    declaration = original.lstrip().startswith(b"<?xml")
    encoding = root.getroottree().docinfo.encoding or "UTF-8"
    return etree.tostring(
        root.getroottree(),
        encoding=encoding,
        xml_declaration=declaration,
    )


def remove_preserving_tail(element: etree._Element) -> None:
    parent = element.getparent()
    if parent is None:
        raise SanitizationError("refusing to remove an XML document root")
    if element.tail:
        previous = element.getprevious()
        if previous is None:
            parent.text = (parent.text or "") + element.tail
        else:
            previous.tail = (previous.tail or "") + element.tail
    parent.remove(element)


def identity_candidates(archive: zipfile.ZipFile, removed_parts: set[str]) -> set[str]:
    candidates: set[str] = set()
    identity_attributes = {"author", "userName", "userId", "email"}
    for name in removed_parts:
        payload = archive.read(name)
        for match in EMAIL_RE.findall(payload.decode("utf-8", errors="ignore")):
            candidates.add(match.strip())
        root = parse_xml(payload, name)
        for element in root.iter():
            for key, value in element.attrib.items():
                if local_name(key) in identity_attributes and len(value.strip()) >= 4:
                    candidates.add(value.strip())
    return candidates


def sanitize_content_types(root: etree._Element, removed_parts: set[str]) -> int:
    removed = 0
    for element in list(root):
        part = element.get("PartName", "").lstrip("/")
        content_type = element.get("ContentType", "").lower()
        if part in removed_parts or (
            local_name(element.tag) == "Override"
            and ("comments" in content_type or content_type.endswith(".people+xml"))
        ):
            remove_preserving_tail(element)
            removed += 1
    return removed


def sanitize_relationships(root: etree._Element, removed_parts: set[str]) -> int:
    removed = 0
    removed_basenames = {PurePosixPath(name).name.lower() for name in removed_parts}
    for element in list(root):
        relation_type = element.get("Type", "").lower()
        target_name = PurePosixPath(element.get("Target", "")).name.lower()
        if (
            "comments" in relation_type
            or relation_type.endswith("/people")
            or target_name in removed_basenames
        ):
            remove_preserving_tail(element)
            removed += 1
    return removed


def remove_comment_markers(root: etree._Element) -> int:
    removed = 0
    references: list[etree._Element] = []
    for element in list(root.iter()):
        if local_name(element.tag) not in COMMENT_MARKERS:
            continue
        parent = element.getparent()
        is_reference = local_name(element.tag) == "commentReference"
        remove_preserving_tail(element)
        removed += 1
        if is_reference and parent is not None and local_name(parent.tag) == "r":
            references.append(parent)
    for run in references:
        parent = run.getparent()
        if parent is not None and all(local_name(child.tag) == "rPr" for child in run):
            remove_preserving_tail(run)
    return removed


def remove_comment_styles(root: etree._Element) -> int:
    removed = 0
    for element in list(root):
        if local_name(element.tag) != "style":
            continue
        style_id = next(
            (
                value
                for key, value in element.attrib.items()
                if local_name(key) == "styleId"
            ),
            "",
        )
        names = [
            next(
                (
                    value
                    for key, value in child.attrib.items()
                    if local_name(key) == "val"
                ),
                "",
            )
            for child in element
            if local_name(child.tag) == "name"
        ]
        if COMMENT_STYLE_RE.fullmatch(style_id) or any(
            re.fullmatch(r"Comment (?:Text|Reference|Subject)[0-9]*", name, re.I)
            for name in names
        ):
            remove_preserving_tail(element)
            removed += 1
    return removed


def unwrap_element(element: etree._Element) -> None:
    parent = element.getparent()
    if parent is None:
        raise SanitizationError("refusing to unwrap an XML document root")
    position = parent.index(element)
    if element.text:
        previous = element.getprevious()
        if previous is None:
            parent.text = (parent.text or "") + element.text
        else:
            previous.tail = (previous.tail or "") + element.text
    children = list(element)
    for child in children:
        element.remove(child)
        parent.insert(position, child)
        position += 1
    tail = element.tail
    parent.remove(element)
    if tail:
        previous = parent[position - 1] if position else None
        if previous is None:
            parent.text = (parent.text or "") + tail
        else:
            previous.tail = (previous.tail or "") + tail


def accept_tracked_revisions(root: etree._Element) -> dict[str, int]:
    changes = {
        "deleted_revisions": 0,
        "accepted_insertions": 0,
        "revision_markers": 0,
        "property_changes": 0,
        "revision_session_metadata": 0,
    }
    for element in list(root.iter()):
        if element.getparent() is None:
            continue
        lname = local_name(element.tag)
        if lname in DELETED_REVISION_ELEMENTS:
            remove_preserving_tail(element)
            changes["deleted_revisions"] += 1
    for element in list(root.iter()):
        if element.getparent() is None:
            continue
        lname = local_name(element.tag)
        if lname in PRESERVED_REVISION_WRAPPERS:
            unwrap_element(element)
            changes["accepted_insertions"] += 1
    for element in list(root.iter()):
        if element.getparent() is None:
            continue
        lname = local_name(element.tag)
        if lname in REVISION_MARKERS:
            remove_preserving_tail(element)
            changes["revision_markers"] += 1
        elif lname.endswith("PrChange"):
            remove_preserving_tail(element)
            changes["property_changes"] += 1
    for element in list(root.iter()):
        if element.getparent() is not None and local_name(element.tag) == "rsids":
            remove_preserving_tail(element)
            changes["revision_session_metadata"] += 1
            continue
        for key in list(element.attrib):
            if local_name(key).lower().startswith("rsid"):
                del element.attrib[key]
                changes["revision_session_metadata"] += 1
    return changes


def sanitize_docx_xml(
    name: str, payload: bytes, removed_parts: set[str]
) -> tuple[bytes, dict[str, int], set[str]]:
    root = parse_xml(payload, name)
    changes = {
        "content_types": 0,
        "relationships": 0,
        "markers": 0,
        "styles": 0,
        "last_modified_by": 0,
    }
    sensitive: set[str] = set()
    for element in root.iter():
        lname = local_name(element.tag)
        if (
            lname in DELETED_REVISION_ELEMENTS
            or lname in PRESERVED_REVISION_WRAPPERS
            or lname.endswith("PrChange")
        ):
            for key, value in element.attrib.items():
                if local_name(key) in {"author", "userName", "userId", "email"}:
                    sensitive.add(value)
    if name == "[Content_Types].xml":
        changes["content_types"] = sanitize_content_types(root, removed_parts)
    if name.endswith(".rels"):
        changes["relationships"] = sanitize_relationships(root, removed_parts)
    changes["markers"] = remove_comment_markers(root)
    if name.startswith("word/") and name.endswith(".xml"):
        for key, value in accept_tracked_revisions(root).items():
            changes[key] = changes.get(key, 0) + value
    if name == "word/styles.xml":
        changes["styles"] = remove_comment_styles(root)
    if name == "docProps/core.xml":
        for element in list(root.iter()):
            if local_name(element.tag) == "lastModifiedBy":
                if element.text:
                    sensitive.add(element.text.strip())
                remove_preserving_tail(element)
                changes["last_modified_by"] += 1
    if not any(changes.values()):
        return payload, changes, sensitive
    return serialize_xml(root, payload), changes, sensitive


def accepted_word_tokens(payload: bytes, label: str) -> tuple[str, ...]:
    root = parse_xml(payload, label)
    remove_comment_markers(root)
    accept_tracked_revisions(root)
    tokens: list[str] = []
    text_elements = {"t", "instrText"}
    structural_elements = {"tab", "br", "cr", "noBreakHyphen", "softHyphen"}
    for element in root.iter():
        lname = local_name(element.tag)
        if lname in text_elements:
            tokens.append(f"{lname}:{element.text or ''}")
        elif lname in structural_elements:
            tokens.append(f"<{lname}>")
    return tuple(tokens)


def sanitize_xlsx_xml(
    name: str, payload: bytes
) -> tuple[bytes, dict[str, int], set[str]]:
    root = parse_xml(payload, name)
    changes = {
        "absolute_paths": 0,
        "revision_nodes": 0,
        "revision_attributes": 0,
        "last_modified_by": 0,
    }
    sensitive: set[str] = set()
    for element in list(root.iter()):
        lname = local_name(element.tag)
        if name == "xl/workbook.xml" and lname in {"absPath", "revisionPtr"}:
            sensitive.update(value for value in element.attrib.values() if value)
            remove_preserving_tail(element)
            changes["absolute_paths" if lname == "absPath" else "revision_nodes"] += 1
            continue
        if name == "docProps/core.xml" and lname == "lastModifiedBy":
            if element.text:
                sensitive.add(element.text.strip())
            remove_preserving_tail(element)
            changes["last_modified_by"] += 1
            continue
        for key in list(element.attrib):
            if (
                local_name(key) in REVISION_ATTRIBUTE_NAMES
                and "revision" in namespace(key).lower()
            ):
                sensitive.add(element.attrib[key])
                del element.attrib[key]
                changes["revision_attributes"] += 1
    if not any(changes.values()):
        return payload, changes, sensitive
    return serialize_xml(root, payload), changes, sensitive


def clone_info(info: zipfile.ZipInfo) -> zipfile.ZipInfo:
    cloned = zipfile.ZipInfo(info.filename, info.date_time)
    for attribute in (
        "compress_type",
        "comment",
        "extra",
        "create_system",
        "create_version",
        "extract_version",
        "reserved",
        "flag_bits",
        "volume",
        "internal_attr",
        "external_attr",
    ):
        setattr(cloned, attribute, getattr(info, attribute))
    return cloned


def package_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == ".docx":
        return "docx"
    if suffix == ".xlsx":
        return "xlsx"
    raise SanitizationError(f"unsupported Office package type: {path}")


def _casefold_payloads(payloads: Iterable[bytes]) -> str:
    return "\n".join(
        payload.decode("utf-8", errors="ignore").casefold() for payload in payloads
    )


def sanitize_package(source: Path, output: Path) -> dict[str, Any]:
    kind = package_kind(source)
    if source.is_symlink() or not source.is_file():
        raise SanitizationError(f"source must be a regular file: {source}")
    if output.exists() or output.is_symlink():
        raise SanitizationError(f"output must be absent: {output}")
    transformed: dict[str, bytes] = {}
    removed_parts: set[str] = set()
    changes: dict[str, int] = {}
    sensitive: set[str] = set()
    with zipfile.ZipFile(source, "r") as archive:
        infos = archive.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise SanitizationError(f"duplicate ZIP member names: {source}")
        if any(info.flag_bits & 0x1 for info in infos):
            raise SanitizationError(f"encrypted ZIP member: {source}")
        bad = archive.testzip()
        if bad is not None:
            raise SanitizationError(f"corrupt ZIP member {bad}: {source}")
        if kind == "docx":
            removed_parts = {name for name in names if COMMENT_PART_RE.fullmatch(name)}
            sensitive |= identity_candidates(archive, removed_parts)
        for info in infos:
            name = info.filename
            if name in removed_parts or not name.endswith((".xml", ".rels")):
                continue
            payload = archive.read(info)
            if kind == "docx":
                updated, observed, found = sanitize_docx_xml(
                    name, payload, removed_parts
                )
                sensitive |= {value for value in found if value}
                for key, value in observed.items():
                    changes[key] = changes.get(key, 0) + value
            else:
                updated, observed, found = sanitize_xlsx_xml(name, payload)
                sensitive |= {value for value in found if value}
                for key, value in observed.items():
                    changes[key] = changes.get(key, 0) + value
            if updated != payload:
                transformed[name] = updated
        retained_text = _casefold_payloads(
            transformed.get(info.filename, archive.read(info))
            for info in infos
            if info.filename not in removed_parts
            and info.filename.endswith((".xml", ".rels"))
        )
        comment_only = {
            value
            for value in sensitive
            if len(value) >= 4 and value.casefold() not in retained_text
        }
        with zipfile.ZipFile(output, "x", allowZip64=True) as destination:
            destination.comment = archive.comment
            for info in infos:
                if info.filename in removed_parts:
                    continue
                copied = clone_info(info)
                if info.filename in transformed:
                    destination.writestr(copied, transformed[info.filename])
                else:
                    with archive.open(info, "r") as src, destination.open(
                        copied, "w", force_zip64=True
                    ) as dst:
                        shutil.copyfileobj(src, dst, CHUNK)
    os.chmod(output, stat.S_IMODE(source.stat().st_mode))
    validation = validate_sanitized(source, output, kind, removed_parts, comment_only)
    validation["changes"] = changes
    validation["removed_parts"] = sorted(removed_parts)
    validation["sensitive_value_sha256"] = sorted(
        hashlib.sha256(value.encode("utf-8")).hexdigest() for value in sensitive
    )
    return validation


def validate_sanitized(
    source: Path,
    output: Path,
    kind: str,
    removed_parts: set[str],
    sensitive: set[str],
) -> dict[str, Any]:
    with zipfile.ZipFile(source) as before, zipfile.ZipFile(output) as after:
        if after.testzip() is not None:
            raise SanitizationError(f"sanitized ZIP integrity failed: {output}")
        before_names = before.namelist()
        after_names = after.namelist()
        if after_names != [name for name in before_names if name not in removed_parts]:
            raise SanitizationError(
                "sanitized ZIP inventory differs beyond removed parts"
            )
        transformed_names: list[str] = []
        retained_xml: list[bytes] = []
        for name in after_names:
            new = after.read(name)
            old = before.read(name)
            if new != old:
                transformed_names.append(name)
            if name.endswith((".xml", ".rels")):
                retained_xml.append(new)
                root = parse_xml(new, name)
                if kind == "docx":
                    for element in root.iter():
                        lname = local_name(element.tag)
                        if lname in COMMENT_MARKERS:
                            raise SanitizationError(f"comment marker remains in {name}")
                        if (
                            lname in DELETED_REVISION_ELEMENTS
                            or lname in PRESERVED_REVISION_WRAPPERS
                            or lname in REVISION_MARKERS
                            or lname.endswith("PrChange")
                        ):
                            raise SanitizationError(
                                f"tracked-revision element remains in {name}: {lname}"
                            )
                        if any(
                            local_name(key).lower().startswith("rsid")
                            for key in element.attrib
                        ):
                            raise SanitizationError(
                                f"revision session metadata remains in {name}"
                            )
                        if name == "docProps/core.xml" and lname == "lastModifiedBy":
                            raise SanitizationError("lastModifiedBy remains")
                        if name.endswith(".rels"):
                            relation_type = element.get("Type", "").lower()
                            target = element.get("Target", "").lower()
                            if (
                                "comments" in relation_type
                                or relation_type.endswith("/people")
                                or COMMENT_PART_RE.search(target)
                            ):
                                raise SanitizationError(
                                    f"comment relationship remains in {name}"
                                )
                        if name == "[Content_Types].xml" and (
                            "comments" in element.get("ContentType", "").lower()
                            or element.get("PartName", "").lstrip("/") in removed_parts
                        ):
                            raise SanitizationError("comment content type remains")
                else:
                    for element in root.iter():
                        lname = local_name(element.tag)
                        if name == "xl/workbook.xml" and lname in {
                            "absPath",
                            "revisionPtr",
                        }:
                            raise SanitizationError(
                                f"private workbook node remains: {lname}"
                            )
                        if name == "docProps/core.xml" and lname == "lastModifiedBy":
                            raise SanitizationError("lastModifiedBy remains")
                        if any(
                            local_name(key) in REVISION_ATTRIBUTE_NAMES
                            and "revision" in namespace(key).lower()
                            for key in element.attrib
                        ):
                            raise SanitizationError(
                                f"revision UUID attribute remains in {name}"
                            )
        retained = _casefold_payloads(retained_xml)
        leaked = [
            hashlib.sha256(value.encode("utf-8")).hexdigest()
            for value in sensitive
            if value.casefold() in retained
        ]
        if leaked:
            raise SanitizationError(f"comment/private identity values remain: {leaked}")
        if kind == "docx":
            story = re.compile(
                r"^word/(?:document|header[0-9]*|footer[0-9]*|"
                r"footnotes|endnotes)\.xml$"
            )
            before_text = tuple(
                token
                for name in before_names
                if story.fullmatch(name)
                for token in accepted_word_tokens(before.read(name), name)
            )
            after_text = tuple(
                token
                for name in after_names
                if story.fullmatch(name)
                for token in accepted_word_tokens(after.read(name), name)
            )
            if before_text != after_text:
                raise SanitizationError("accepted-revisions visible Word text changed")
        else:
            for name in before_names:
                if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                    old = parse_xml(before.read(name), name)
                    new = parse_xml(after.read(name), name)
                    for root in (old, new):
                        for element in root.iter():
                            for key in list(element.attrib):
                                if (
                                    local_name(key) in REVISION_ATTRIBUTE_NAMES
                                    and "revision" in namespace(key).lower()
                                ):
                                    del element.attrib[key]
                    if etree.tostring(old, method="c14n") != etree.tostring(
                        new, method="c14n"
                    ):
                        raise SanitizationError(f"worksheet content changed: {name}")
    return {
        "kind": kind,
        "source_sha256": sha256_file(source),
        "sanitized_sha256": sha256_file(output),
        "transformed_members": transformed_names,
    }


def _safe_sources(root: Path, values: Sequence[Path], backup: Path) -> list[Path]:
    result: list[Path] = []
    for value in values:
        lexical = value if value.is_absolute() else root / value
        if lexical.is_symlink():
            raise SanitizationError(f"source may not be a symlink: {lexical}")
        source = lexical.resolve()
        if not source.is_file() or not source.is_relative_to(root):
            raise SanitizationError(f"source must be a file below {root}: {source}")
        if source.is_relative_to(backup):
            raise SanitizationError(f"backup content cannot be sanitized: {source}")
        package_kind(source)
        result.append(source)
    if len(result) != len(set(result)):
        raise SanitizationError("duplicate source path")
    return result


def sanitize_in_place(
    root_path: Path, source_values: Sequence[Path], backup_path: Path
) -> list[dict[str, Any]]:
    root = root_path.expanduser().absolute()
    if root.is_symlink():
        raise SanitizationError(f"root may not be a symlink: {root}")
    root = root.resolve()
    backup = backup_path if backup_path.is_absolute() else root / backup_path
    if backup.exists() or backup.is_symlink():
        raise SanitizationError(f"backup directory must be absent: {backup}")
    backup = backup.parent.resolve() / backup.name
    sources = _safe_sources(root, source_values, backup)
    backup.mkdir(parents=True, mode=0o700)
    staged: list[tuple[Path, Path, Path, dict[str, Any]]] = []
    try:
        for index, source in enumerate(sources):
            relative = source.relative_to(root)
            saved = backup / relative
            saved.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, saved)
            if sha256_file(saved) != sha256_file(source):
                raise SanitizationError(f"backup checksum differs: {relative}")
            temporary = source.with_name(
                f".{source.name}.sanitized.{os.getpid()}.{index}"
            )
            if temporary.exists() or temporary.is_symlink():
                raise SanitizationError(f"sanitization stage exists: {temporary}")
            report = sanitize_package(source, temporary)
            staged.append((source, saved, temporary, report))
        for source, _saved, temporary, _report in staged:
            os.replace(temporary, source)
        manifest = {
            "schema": "apotome-office-sanitization-backup-v1",
            "files": [
                {
                    "path": source.relative_to(root).as_posix(),
                    "backup": saved.relative_to(backup).as_posix(),
                    "original_size": saved.stat().st_size,
                    "original_sha256": sha256_file(saved),
                    "sanitized_size": source.stat().st_size,
                    "sanitized_sha256": sha256_file(source),
                    "report": report,
                }
                for source, saved, _temporary, report in staged
            ],
        }
        manifest_path = backup / "backup_manifest.json"
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.chmod(manifest_path, 0o600)
        return manifest["files"]
    except Exception:
        for _source, _saved, temporary, _report in staged:
            if temporary.exists() and not temporary.is_symlink():
                temporary.unlink()
        raise


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--backup-dir", type=Path, required=True)
    parser.add_argument("documents", type=Path, nargs="+")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    reports = sanitize_in_place(args.root, args.documents, args.backup_dir)
    for report in reports:
        print(
            f"[PASS] {report['path']}: {report['original_sha256']} -> "
            f"{report['sanitized_sha256']}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (SanitizationError, OSError, zipfile.BadZipFile) as exc:
        print(f"[FAIL-CLOSED] {exc}", file=sys.stderr)
        raise SystemExit(2)
