#!/usr/bin/env python3
"""Refresh only HIL terminology and canonical paper art in the final thesis DOCX."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import posixpath
import struct
import tempfile
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from lxml import etree


ROOT = Path(__file__).resolve().parents[2]
MANUSCRIPT_SOURCES = ROOT / "manuscript_sources"
if not MANUSCRIPT_SOURCES.is_dir():
    MANUSCRIPT_SOURCES = ROOT / "thesis_and_manuscript/insumos/manuscript_sources"
if not MANUSCRIPT_SOURCES.is_dir():
    MANUSCRIPT_SOURCES = ROOT / "revision_profesional_20260906/insumos/manuscript_sources"
MANIFEST = MANUSCRIPT_SOURCES / "TESIS_FINAL_NS4_spanish_panel_manifest_20260830.csv"
DOCUMENT = "word/document.xml"
RELATIONSHIPS = "word/_rels/document.xml.rels"
W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
P_NS = "http://schemas.openxmlformats.org/package/2006/relationships"
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"
R_EMBED = f"{{{R_NS}}}embed"
NS = {"w": W_NS, "a": A_NS, "p": P_NS}

MASTER_PATHS = (
    "Fig1/Figure_1_behavior_experiment_1.png",
    "Fig2/Figure_2.png",
    "Fig3/Figure_3_cFos_NPY.png",
    "Fig4/Figure_4.png",
    "Fig5/Figure_5.png",
    "FigS1/Figure_S1_single_bottle_male.png",
    "FigS2/Figure_S2_single_bottle_female.png",
    "FigS3/Figure_S3.png",
    "FigS4/Figure_S4_behavior_preference.png",
    "FigS5/Figure_S5.png",
)

TEXT_REPLACEMENTS = (
    ("participaron en la revisión anatómica y", "participaron en la revisión anatómica mediante HIL y"),
    ("Delimitación anatómica revisada por una persona", "Delimitación anatómica aceptada mediante HIL"),
    ("HIL: humano en el circuito", "HIL: human-in-the-loop"),
    ("delimitación anatómica revisada por una persona", "delimitación anatómica aceptada mediante HIL"),
    ("máscaras aceptadas por revisión humana", "máscaras aceptadas mediante HIL"),
    ("aceptadas por revisión humana", "aceptadas mediante HIL"),
    ("Las máscaras ARC/EM/VMN se revisaron en las coordenadas originales", "Las máscaras ARC/EM/VMN se aceptaron mediante HIL en las coordenadas originales"),
    ("Los límites ventriculares de los campos representativos se revisaron manualmente", "Los límites ventriculares de los campos representativos se aceptaron mediante HIL"),
    ("las 76 secciones revisadas", "las 76 secciones aceptadas mediante HIL"),
    ("máscaras ARC, EM y VMN revisadas imagen por imagen", "máscaras ARC, EM y VMN aceptadas mediante HIL imagen por imagen"),
    ("máscaras de ARC, EM y VMN revisadas imagen por imagen", "máscaras de ARC, EM y VMN aceptadas mediante HIL imagen por imagen"),
    ("la anatomía se revisó sin usar la intensidad de los marcadores", "la anatomía se aceptó mediante HIL sin usar la intensidad de los marcadores"),
    ("corrección humana de las regiones de interés", "corrección mediante HIL de las regiones de interés"),
    ("la revisión anatómica establecieron", "la revisión anatómica mediante HIL establecieron"),
    ("aceptada mediante revisión HIL", "aceptada mediante HIL"),
    ("máscaras anatómicas revisadas de ARC y EM", "máscaras anatómicas aceptadas mediante HIL de ARC y EM"),
    ("Los límites de ARC y EM proceden de revisión HIL sobre DAPI", "Los límites de ARC y EM proceden de máscaras aceptadas mediante HIL sobre DAPI"),
    ("máscaras aceptadas por revisión HIL", "máscaras aceptadas mediante HIL"),
    ("definidas por revisión HIL", "definidas mediante HIL"),
    ("delimitación anatómica humana documentada", "delimitación anatómica mediante HIL documentada"),
    ("Reviewer-defined anatomy", "HIL-defined anatomy"),
    ("reviewer-defined anatomy", "HIL-defined anatomy"),
    ("human review", "HIL"),
    ("Human input", "HIL input"),
    ("human-accepted", "HIL-accepted"),
    ("reviewer-accepted", "HIL-accepted"),
    ("61 of 61 sections reviewed image by image; automated masks that were not reviewed were not used", "61 of 61 sections were accepted via HIL image by image; automated masks without HIL acceptance were not used"),
)

FORBIDDEN_TEXT = (
    "humano en el circuito",
    "revisión humana",
    "revisada por una persona",
    "corrección humana",
    "revisión HIL",
    "delimitación anatómica humana",
    "human review",
    "reviewer-defined",
    "human-accepted",
    "human input",
    "reviewer-accepted",
    "secciones revisadas",
    "revisadas imagen por imagen",
    "se revisaron en las coordenadas originales",
    "se revisaron manualmente",
    "la anatomía se revisó",
    "sections reviewed image by image",
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_png(path: Path) -> tuple[int, int, str]:
    data = path.read_bytes()
    if len(data) < 24 or data[:8] != b"\x89PNG\r\n\x1a\n" or data[12:16] != b"IHDR":
        raise ValueError(f"Not a valid PNG: {path}")
    width, height = struct.unpack(">II", data[16:24])
    if width < 1 or height < 1:
        raise ValueError(f"Invalid PNG dimensions for {path}: {width}x{height}")
    return width, height, sha256_bytes(data)


def set_space(node: etree._Element) -> None:
    value = node.text or ""
    if value[:1].isspace() or value[-1:].isspace():
        node.set(XML_SPACE, "preserve")
    else:
        node.attrib.pop(XML_SPACE, None)


def replace_across_nodes(nodes: list[etree._Element], old: str, new: str) -> int:
    count = 0
    while True:
        values = [node.text or "" for node in nodes]
        full = "".join(values)
        start = full.find(old)
        if start < 0:
            return count
        end = start + len(old)
        first = last = None
        cursor = 0
        for index, value in enumerate(values):
            following = cursor + len(value)
            if first is None and cursor <= start < following:
                first = (index, start - cursor)
            if cursor < end <= following:
                last = (index, end - cursor)
                break
            cursor = following
        if first is None or last is None:
            raise RuntimeError(f"Could not map text replacement across w:t nodes: {old!r}")
        first_index, first_offset = first
        last_index, last_offset = last
        if first_index == last_index:
            nodes[first_index].text = values[first_index][:first_offset] + new + values[first_index][last_offset:]
            set_space(nodes[first_index])
        else:
            nodes[first_index].text = values[first_index][:first_offset] + new
            set_space(nodes[first_index])
            for index in range(first_index + 1, last_index):
                nodes[index].text = ""
                set_space(nodes[index])
            nodes[last_index].text = values[last_index][last_offset:]
            set_space(nodes[last_index])
        count += 1


def patch_document(data: bytes) -> tuple[bytes, dict[str, int], int]:
    parser = etree.XMLParser(resolve_entities=False, no_network=True, remove_blank_text=False)
    root = etree.fromstring(data, parser)
    counts = {old: 0 for old, _ in TEXT_REPLACEMENTS}
    for paragraph in root.xpath("//w:p", namespaces=NS):
        nodes = list(paragraph.xpath(".//w:t", namespaces=NS))
        if not nodes:
            continue
        for old, new in TEXT_REPLACEMENTS:
            counts[old] += replace_across_nodes(nodes, old, new)
    text = "\n".join(
        "".join(node.text or "" for node in paragraph.xpath(".//w:t", namespaces=NS))
        for paragraph in root.xpath("//w:p", namespaces=NS)
    )
    folded = text.casefold()
    remaining = [pattern for pattern in FORBIDDEN_TEXT if pattern.casefold() in folded]
    if remaining:
        raise ValueError(f"Deprecated HIL wording remains in document.xml: {remaining}")
    if "hil: human-in-the-loop" not in folded:
        raise ValueError("Required glossary definition 'HIL: human-in-the-loop' is missing")
    if folded.count("hil") < 15:
        raise ValueError("Unexpectedly few HIL references after patch")
    output = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)
    return output, counts, folded.count("hil")


def load_manifest() -> list[dict[str, object]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 63 or [int(row["index"]) for row in rows] != list(range(1, 64)):
        raise ValueError("Spanish panel manifest must contain exactly indices 1..63")
    loaded = []
    for row in rows:
        source = (ROOT / row["relative_path"]).resolve()
        if ROOT.resolve() not in source.parents or not source.is_file():
            raise ValueError(f"Invalid or missing canonical panel path: {source}")
        width, height, current_hash = validate_png(source)
        old_hash = row["sha256"].lower()
        if len(old_hash) != 64 or any(char not in "0123456789abcdef" for char in old_hash):
            raise ValueError(f"Invalid historical SHA-256 in manifest row {row['index']}")
        loaded.append({"source": source, "relative_path": row["relative_path"], "old_hash": old_hash, "sha256": current_hash, "width": width, "height": height})
    return loaded


def member_for_blip(blip: etree._Element, relationships: dict[str, etree._Element]) -> str:
    rel_id = blip.get(R_EMBED)
    if not rel_id or rel_id not in relationships:
        raise ValueError(f"Image blip has unresolved relationship: {rel_id}")
    rel = relationships[rel_id]
    if rel.get("TargetMode") == "External" or not (rel.get("Type") or "").endswith("/image"):
        raise ValueError(f"Relationship {rel_id} is not an internal image")
    member = posixpath.normpath(posixpath.join("word", rel.get("Target") or ""))
    if not member.startswith("word/media/") or member.startswith("../"):
        raise ValueError(f"Unsafe image target for {rel_id}: {member}")
    return member


def write_json_atomic(path: Path, payload: dict[str, object], force: bool) -> None:
    if path.exists() and not force:
        raise FileExistsError(f"Receipt already exists (use --force): {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def refresh(input_path: Path, output_path: Path, receipt_path: Path | None, force: bool) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if output_path.exists() and not force:
        raise FileExistsError(f"Output already exists (use --force): {output_path}")
    panel_rows = load_manifest()
    master_rows = []
    for relative_path in MASTER_PATHS:
        source = (ROOT / relative_path).resolve()
        width, height, current_hash = validate_png(source)
        master_rows.append({"source": source, "relative_path": relative_path, "sha256": current_hash, "width": width, "height": height})

    with zipfile.ZipFile(input_path, "r") as source_zip:
        infos = source_zip.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise ValueError("DOCX contains duplicate ZIP member names")
        if DOCUMENT not in names or RELATIONSHIPS not in names:
            raise ValueError("DOCX is missing document XML or its relationships")
        before_hashes = {name: sha256_bytes(source_zip.read(name)) for name in names}
        document_data, text_counts, hil_count = patch_document(source_zip.read(DOCUMENT))
        parser = etree.XMLParser(resolve_entities=False, no_network=True, remove_blank_text=False)
        document_root = etree.fromstring(source_zip.read(DOCUMENT), parser)
        rel_root = etree.fromstring(source_zip.read(RELATIONSHIPS), parser)
        relationships = {rel.get("Id"): rel for rel in rel_root.xpath("//p:Relationship", namespaces=NS)}
        blips = document_root.xpath("//a:blip", namespaces=NS)
        if len(blips) != 84:
            raise ValueError(f"Expected exactly 84 document image blips, found {len(blips)}")
        panel_members = [member_for_blip(blip, relationships) for blip in blips[11:74]]
        master_members = [member_for_blip(blip, relationships) for blip in blips[74:84]]
        expected_panels = [f"word/media/image{number}.png" for number in range(30, 34)] + ["word/media/image29.png"] + [f"word/media/image{number}.png" for number in range(34, 92)]
        expected_masters = [f"word/media/image{number}.png" for number in range(92, 102)]
        if panel_members != expected_panels or master_members != expected_masters:
            raise ValueError("Unexpected final-thesis panel/master image relationship contract")

        media_sources: dict[str, Path] = {}
        media_receipt = []
        for member, row in zip(panel_members, panel_rows, strict=True):
            before = before_hashes[member]
            if before not in {row["old_hash"], row["sha256"]}:
                raise ValueError(f"Embedded Spanish panel is neither historical nor current canonical art: {member}")
            media_sources[member] = row["source"]
            media_receipt.append({"kind": "spanish_panel", "member": member, "source": row["relative_path"], "before_sha256": before, "after_sha256": row["sha256"], "width": row["width"], "height": row["height"]})
        for member, row in zip(master_members, master_rows, strict=True):
            media_sources[member] = row["source"]
            media_receipt.append({"kind": "paper_master", "member": member, "source": row["relative_path"], "before_sha256": before_hashes[member], "after_sha256": row["sha256"], "width": row["width"], "height": row["height"]})

        output_path.parent.mkdir(parents=True, exist_ok=True)
        fd, temporary = tempfile.mkstemp(prefix=f".{output_path.name}.", suffix=".tmp", dir=output_path.parent)
        os.close(fd)
        try:
            with zipfile.ZipFile(temporary, "w") as output_zip:
                output_zip.comment = source_zip.comment
                for info in infos:
                    if info.filename == DOCUMENT:
                        data = document_data
                    elif info.filename in media_sources:
                        data = media_sources[info.filename].read_bytes()
                    else:
                        data = source_zip.read(info.filename)
                    output_zip.writestr(info, data)
            with open(temporary, "rb") as handle:
                os.fsync(handle.fileno())
            os.replace(temporary, output_path)
        except Exception:
            try:
                os.unlink(temporary)
            except FileNotFoundError:
                pass
            raise

    with zipfile.ZipFile(output_path, "r") as result_zip:
        result_names = [info.filename for info in result_zip.infolist()]
        if result_names != names or result_zip.testzip() is not None:
            raise ValueError("Output DOCX ZIP structure or CRC validation failed")
        changed = {DOCUMENT, *media_sources}
        for name in names:
            result_hash = sha256_bytes(result_zip.read(name))
            if name not in changed and result_hash != before_hashes[name]:
                raise ValueError(f"Unrelated DOCX member changed: {name}")
        for member, source in media_sources.items():
            if sha256_bytes(result_zip.read(member)) != sha256_file(source):
                raise ValueError(f"Canonical media replacement hash mismatch: {member}")
        patch_document(result_zip.read(DOCUMENT))

    payload = {
        "status": "complete",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(input_path.resolve()),
        "output": str(output_path.resolve()),
        "input_sha256": sha256_file(input_path),
        "output_sha256": sha256_file(output_path),
        "zip_members": len(names),
        "image_blips": 84,
        "spanish_panels_refreshed": len(panel_rows),
        "paper_masters_refreshed": len(master_rows),
        "unrelated_members_verified_unchanged": len(names) - len({DOCUMENT, *media_sources}),
        "text_replacements_applied": {old: count for old, count in text_counts.items() if count},
        "hil_mentions": hil_count,
        "media": media_receipt,
    }
    if receipt_path is not None:
        write_json_atomic(receipt_path, payload, force)
        payload["receipt"] = str(receipt_path.resolve())
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    result = refresh(args.input.resolve(), args.output.resolve(), args.receipt.resolve() if args.receipt else None, args.force)
    print(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
