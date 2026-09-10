#!/usr/bin/env python3
"""Update only single-bottle terminology and art in the final thesis DOCX.

The source DOCX is never modified. The output changes word/document.xml and
the nine image members that display Figure 1, S1, or S2 consumption wording.
Every other ZIP member is copied and verified byte-for-byte.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import struct
import sys
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
sys.path.insert(0, str(MANUSCRIPT_SOURCES))

from caps_en import CAPS_EN  # noqa: E402
from thesis_figure_revision_content import (  # noqa: E402
    FIGURE_CAPTION_OVERRIDES,
    FIGURE_RESULT_NARRATIVES,
)


W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
XMLSPACE = "{http://www.w3.org/XML/1998/namespace}space"
DOCUMENT_XML = "word/document.xml"
CHUNK = 8 * 1024 * 1024

MEDIA_SOURCES = {
    "word/media/image30.png": ROOT
    / "Fig1/panels/Figure_1_behavior_experiment_1_Panel_A_spanish.png",
    "word/media/image31.png": ROOT
    / "Fig1/panels/Figure_1_behavior_experiment_1_Panel_B_spanish.png",
    "word/media/image63.png": ROOT
    / "FigS1/panels/Figure_S1_single_bottle_male_Panel_A_spanish.png",
    "word/media/image64.png": ROOT
    / "FigS1/panels/Figure_S1_single_bottle_male_Panel_B_spanish.png",
    "word/media/image67.png": ROOT
    / "FigS2/panels/Figure_S2_single_bottle_female_Panel_A_spanish.png",
    "word/media/image68.png": ROOT
    / "FigS2/panels/Figure_S2_single_bottle_female_Panel_B_spanish.png",
    "word/media/image92.png": ROOT / "Fig1/Figure_1_behavior_experiment_1.png",
    "word/media/image97.png": ROOT / "FigS1/Figure_S1_single_bottle_male.png",
    "word/media/image98.png": ROOT / "FigS2/Figure_S2_single_bottle_female.png",
}

# Unique prefixes are accepted in either old or corrected form. This makes
# reruns idempotent while still failing closed if document structure drifts.
PARAGRAPH_TARGETS = (
    (
        (
            "Las trayectorias por jaula mostraron una extracción acumulada",
            "Las trayectorias por jaula mostraron un consumo acumulado",
        ),
        FIGURE_RESULT_NARRATIVES[9],
    ),
    (
        (
            "En el día 6 se detectó una diferencia global en el volumen retirado",
            "En el día 6 se detectó una diferencia global en el consumo",
        ),
        FIGURE_RESULT_NARRATIVES[10],
    ),
    (
        ("Figura 9. Volumen retirado", "Figura 9. Consumo"),
        FIGURE_CAPTION_OVERRIDES[9],
    ),
    (
        ("Figura 10. Volumen retirado", "Figura 10. Consumo"),
        FIGURE_CAPTION_OVERRIDES[10],
    ),
    (
        (
            "En machos, las trayectorias de volumen",
            "En machos, las trayectorias de consumo",
        ),
        FIGURE_RESULT_NARRATIVES[43],
    ),
    (
        (
            "El volumen retirado por machos en el día 6",
            "El consumo por machos en el día 6",
            "El consumo de los machos en el día 6",
        ),
        FIGURE_RESULT_NARRATIVES[44],
    ),
    (
        ("Figura 43. Volumen retirado", "Figura 43. Consumo"),
        FIGURE_CAPTION_OVERRIDES[43],
    ),
    (
        ("Figura 44. Volumen retirado", "Figura 44. Consumo"),
        FIGURE_CAPTION_OVERRIDES[44],
    ),
    (
        (
            "En hembras, la trayectoria de volumen",
            "En hembras, la trayectoria de consumo",
        ),
        FIGURE_RESULT_NARRATIVES[47],
    ),
    (
        (
            "El volumen del día 6 en hembras",
            "El consumo del día 6 en hembras",
        ),
        FIGURE_RESULT_NARRATIVES[48],
    ),
    (
        ("Figura 47. Volumen retirado", "Figura 47. Consumo"),
        FIGURE_CAPTION_OVERRIDES[47],
    ),
    (
        ("Figura 48. Volumen retirado", "Figura 48. Consumo"),
        FIGURE_CAPTION_OVERRIDES[48],
    ),
    (
        (
            "(A) Cage-level volume removed from the bottle between day 1",
            "(A) Cage-level consumption from day 1",
        ),
        CAPS_EN[628],
    ),
    (
        (
            "(A) Volume removed from the bottle per male cage between day 1",
            "(A) Consumption from day 1 through day 6 per male cage",
        ),
        CAPS_EN[660],
    ),
    (
        (
            "(A) Volume removed from the bottle per female cage between day 1",
            "(A) Consumption from day 1 through day 6 per female cage",
        ),
        CAPS_EN[666],
    ),
)

TEXT_REPLACEMENTS = (
    (
        "El volumen retirado se calculó a partir del cambio de masa de la botella "
        "bajo la aproximación de 1 g/mL; por ello se denomina volumen retirado y "
        "no ingesta individual.",
        "El consumo se calculó a partir del cambio de masa de la botella bajo la "
        "aproximación de 1 g/mL y se informó por jaula, no como ingesta individual.",
    ),
    (
        "Figura 9. Volumen retirado por jaula durante la adaptación",
        "Figura 9. Consumo por jaula desde el día 1",
    ),
    (
        "Figura 10. Volumen retirado por jaula en el día 6",
        "Figura 10. Consumo por jaula en el día 6",
    ),
    (
        "Figura 43. Volumen retirado por jaula en machos",
        "Figura 43. Consumo por jaula desde el día 1 en machos",
    ),
    (
        "Figura 44. Volumen retirado por jaula de machos en el día 6",
        "Figura 44. Consumo por jaula de machos en el día 6",
    ),
    (
        "Figura 47. Volumen retirado por jaula en hembras",
        "Figura 47. Consumo por jaula desde el día 1 en hembras",
    ),
    (
        "Figura 48. Volumen retirado por jaula de hembras en el día 6",
        "Figura 48. Consumo por jaula de hembras en el día 6",
    ),
    (
        "Evaluar el volumen retirado de las soluciones",
        "Evaluar el consumo de las soluciones desde el día 1",
    ),
    (
        "Se registraron el peso corporal y el volumen retirado de la botella.",
        "Se registraron el peso corporal y el consumo desde el día 1.",
    ),
    (
        "registros diarios de volumen retirado",
        "registros de consumo desde el día 1",
    ),
    (
        "El volumen retirado y el cambio de peso",
        "El consumo desde el día 1 y el cambio de peso",
    ),
    ("mayor volumen del grupo Sacarosa", "mayor consumo del grupo Sacarosa"),
    ("mayor extracción de sacarosa", "mayor consumo de sacarosa"),
    ("mayor extracción de Sacarosa", "mayor consumo de Sacarosa"),
    ("alta extracción de sacarosa", "alto consumo de sacarosa"),
    ("alta extracción de Sacarosa", "alto consumo de Sacarosa"),
    (
        "Este ANOVA nominal positivo de volumen",
        "Este ANOVA nominal positivo de consumo",
    ),
    (
        "diferencia global positiva de volumen retirado",
        "diferencia global positiva de consumo",
    ),
    ("Day-6 volume removed", "Day-6 consumption"),
    (
        "experimental unit for intake and body weight",
        "experimental unit for consumption and body weight",
    ),
    ("Volume removed", "Consumption"),
    ("volume removed", "consumption"),
    ("Volume removal", "Consumption"),
    ("volume removal", "consumption"),
    ("Volumen retirado", "Consumo"),
    ("volumen retirado", "consumo"),
)

FORBIDDEN_TEXT = (
    "volume removed",
    "volume removal",
    "cumulative bottle-volume removal",
    "solution-removal",
    "volumen retirado",
    "volumen retirado acumulado",
    "eliminación de solución",
    "extracción acumulada",
    "mayor extracción",
    "alta extracción",
    "trayectorias de volumen",
    "trayectoria de volumen",
)


class UpdateError(RuntimeError):
    """The targeted DOCX update could not be completed safely."""


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def paragraph_text(paragraph: etree._Element) -> str:
    return "".join(node.text or "" for node in paragraph.iter(W + "t"))


def set_space(node: etree._Element) -> None:
    value = node.text or ""
    if value[:1].isspace() or value[-1:].isspace():
        node.set(XMLSPACE, "preserve")
    else:
        node.attrib.pop(XMLSPACE, None)


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
            raise UpdateError(f"Could not map replacement across text nodes: {old!r}")
        first_index, first_offset = first
        last_index, last_offset = last
        if first_index == last_index:
            nodes[first_index].text = (
                values[first_index][:first_offset]
                + new
                + values[first_index][last_offset:]
            )
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


def replace_paragraph_text(paragraph: etree._Element, desired: str) -> bool:
    nodes = list(paragraph.iter(W + "t"))
    if not nodes:
        raise UpdateError("Target paragraph contains no text nodes")
    current = "".join(node.text or "" for node in nodes)
    if current == desired:
        return False
    nodes[0].text = desired
    set_space(nodes[0])
    for node in nodes[1:]:
        node.text = ""
        set_space(node)
    return True


def patch_document(data: bytes) -> tuple[bytes, dict[str, object]]:
    parser = etree.XMLParser(
        resolve_entities=False, no_network=True, remove_blank_text=False, huge_tree=True
    )
    root = etree.fromstring(data, parser)
    paragraphs = list(root.iter(W + "p"))
    target_changes = 0
    for prefixes, desired in PARAGRAPH_TARGETS:
        matches = [
            paragraph
            for paragraph in paragraphs
            if paragraph_text(paragraph).startswith(prefixes)
        ]
        if len(matches) != 1:
            raise UpdateError(
                f"Expected one paragraph beginning with {prefixes!r}; found {len(matches)}"
            )
        target_changes += int(replace_paragraph_text(matches[0], desired))

    replacement_counts: dict[str, int] = {}
    for old, new in TEXT_REPLACEMENTS:
        count = 0
        for paragraph in paragraphs:
            nodes = list(paragraph.iter(W + "t"))
            if nodes:
                count += replace_across_nodes(nodes, old, new)
        replacement_counts[old] = count

    full_text = "\n".join(paragraph_text(paragraph) for paragraph in paragraphs)
    folded = full_text.casefold()
    remaining = [value for value in FORBIDDEN_TEXT if value.casefold() in folded]
    if remaining:
        raise UpdateError(f"Deprecated single-bottle wording remains: {remaining}")
    required = (
        "consumption from day 1",
        "consumo desde el día 1",
        "consumo por jaula desde el día 1",
    )
    missing = [value for value in required if value.casefold() not in folded]
    if missing:
        raise UpdateError(f"Required corrected wording is missing: {missing}")

    revision_tags = (
        "commentRangeStart",
        "commentReference",
        "ins",
        "del",
        "moveFrom",
        "moveTo",
    )
    present_revision_tags = [
        name for name in revision_tags if root.find(f".//{W}{name}") is not None
    ]
    if present_revision_tags:
        raise UpdateError(
            f"Comments or tracked-revision markup present: {present_revision_tags}"
        )

    output = etree.tostring(
        root,
        xml_declaration=True,
        encoding=root.getroottree().docinfo.encoding or "UTF-8",
        standalone=True,
    )
    return output, {
        "paragraph_targets_changed": target_changes,
        "text_replacement_counts": replacement_counts,
        "forbidden_wording_absent": True,
        "comments_and_revision_markup_absent": True,
    }


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


def png_dimensions(payload: bytes) -> tuple[int, int]:
    if len(payload) < 24 or payload[:8] != b"\x89PNG\r\n\x1a\n":
        raise UpdateError("Expected a valid PNG payload")
    return struct.unpack(">II", payload[16:24])


def update_docx(source: Path, output: Path) -> dict[str, object]:
    if not source.is_file():
        raise UpdateError(f"Input DOCX not found: {source}")
    if output.exists() or output.is_symlink():
        raise UpdateError(f"Output must be absent: {output}")
    for member, image_path in MEDIA_SOURCES.items():
        if not image_path.is_file():
            raise UpdateError(f"Canonical image is missing for {member}: {image_path}")

    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    os.close(fd)
    temporary = Path(temporary_name)
    temporary.unlink()
    try:
        with zipfile.ZipFile(source, "r") as archive:
            if archive.testzip() is not None:
                raise UpdateError("Input DOCX ZIP integrity check failed")
            infos = archive.infolist()
            names = [info.filename for info in infos]
            if len(names) != len(set(names)):
                raise UpdateError("Input DOCX has duplicate ZIP member names")
            if DOCUMENT_XML not in names:
                raise UpdateError(f"Input DOCX lacks {DOCUMENT_XML}")
            if any(member not in names for member in MEDIA_SOURCES):
                raise UpdateError("One or more expected Figure 1/S1/S2 image members are missing")
            if any(
                name.startswith("word/comments")
                or name.startswith("word/people")
                for name in names
            ):
                raise UpdateError("Input DOCX contains comments or people parts")

            before_hashes = {
                info.filename: sha256_bytes(archive.read(info.filename))
                for info in infos
            }
            document_payload, text_report = patch_document(
                archive.read(DOCUMENT_XML)
            )
            media_report = []
            replacement_payloads: dict[str, bytes] = {}
            for member, image_path in MEDIA_SOURCES.items():
                before = archive.read(member)
                after = image_path.read_bytes()
                before_dimensions = png_dimensions(before)
                after_dimensions = png_dimensions(after)
                if before_dimensions != after_dimensions:
                    raise UpdateError(
                        f"Image dimensions changed for {member}: "
                        f"{before_dimensions} != {after_dimensions}"
                    )
                replacement_payloads[member] = after
                media_report.append(
                    {
                        "member": member,
                        "source": str(image_path.relative_to(ROOT)),
                        "dimensions_px": list(after_dimensions),
                        "before_sha256": sha256_bytes(before),
                        "after_sha256": sha256_bytes(after),
                    }
                )

            with zipfile.ZipFile(temporary, "x", allowZip64=True) as destination:
                destination.comment = archive.comment
                for info in infos:
                    copied = clone_info(info)
                    if info.filename == DOCUMENT_XML:
                        destination.writestr(copied, document_payload)
                    elif info.filename in replacement_payloads:
                        destination.writestr(copied, replacement_payloads[info.filename])
                    else:
                        with archive.open(info, "r") as src, destination.open(
                            copied, "w", force_zip64=True
                        ) as dst:
                            shutil.copyfileobj(src, dst, CHUNK)

        with zipfile.ZipFile(temporary, "r") as candidate:
            if candidate.testzip() is not None:
                raise UpdateError("Output DOCX ZIP integrity check failed")
            output_names = [info.filename for info in candidate.infolist()]
            if output_names != names:
                raise UpdateError("DOCX ZIP member inventory or order changed")
            changed_members = {DOCUMENT_XML, *MEDIA_SOURCES}
            for name in names:
                after_hash = sha256_bytes(candidate.read(name))
                if name not in changed_members and after_hash != before_hashes[name]:
                    raise UpdateError(f"Unrelated DOCX member changed: {name}")
            for member, payload in replacement_payloads.items():
                if candidate.read(member) != payload:
                    raise UpdateError(f"Embedded image differs from canonical source: {member}")
            patch_document(candidate.read(DOCUMENT_XML))

        os.chmod(temporary, source.stat().st_mode & 0o777)
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    return {
        "schema": "apotome-single-bottle-consumption-docx-update-v1",
        "recorded_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {
            "path": str(source),
            "size_bytes": source.stat().st_size,
            "sha256": sha256_file(source),
        },
        "output": {
            "path": str(output),
            "size_bytes": output.stat().st_size,
            "sha256": sha256_file(output),
        },
        "changed_zip_members": [DOCUMENT_XML, *MEDIA_SOURCES],
        "unchanged_zip_members_verified": True,
        "text": text_report,
        "media": media_report,
    }


def write_receipt(path: Path, report: dict[str, object]) -> None:
    if path.exists() or path.is_symlink():
        raise UpdateError(f"Receipt must be absent: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = update_docx(args.input.resolve(), args.output.resolve())
    if args.receipt:
        write_receipt(args.receipt.resolve(), report)
    print(
        f"[PASS] updated {args.output}: {len(report['changed_zip_members'])} "
        "targeted DOCX members; every other member verified unchanged"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
