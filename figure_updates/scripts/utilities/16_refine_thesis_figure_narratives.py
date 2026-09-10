#!/usr/bin/env python3
"""Refine thesis result figures without rebuilding or recompressing their media.

The patch is deliberately semantic: it identifies the Spanish thesis captions
``Figura 9.`` through ``Figura 71.`` by their text and verifies that each one
is immediately preceded by an image paragraph.  It then inserts one marked
results paragraph immediately before that image, replaces the caption with the
reviewed text, removes the obsolete group-note paragraphs, and normalizes every
paragraph using the Word ``Caption`` style (including the paper appendix).

Only ``word/document.xml`` and ``word/styles.xml`` are changed.  Every ZIP
member is otherwise copied with its original metadata, and the SHA-256 digest
of every ``word/media/`` member is checked after writing.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import json
import os
import re
import sys
import tempfile
import zipfile
from pathlib import Path

from lxml import etree


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MANUSCRIPT_SOURCES = ROOT / "manuscript_sources"
if not MANUSCRIPT_SOURCES.is_dir():
    MANUSCRIPT_SOURCES = ROOT / "thesis_and_manuscript/insumos/manuscript_sources"
if not MANUSCRIPT_SOURCES.is_dir():
    MANUSCRIPT_SOURCES = ROOT / "revision_profesional_20260906/insumos/manuscript_sources"
sys.path.insert(0, str(MANUSCRIPT_SOURCES))

from thesis_figure_revision_content import (  # noqa: E402
    FIGURE_CAPTION_OVERRIDES,
    FIGURE_RESULT_NARRATIVES,
    INTRO_FIGURE_CAPTION_OVERRIDES,
    INTRO_FIGURE_RESULT_NARRATIVES,
)


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = f"{{{W_NS}}}"
NS = {"w": W_NS}
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"

DOCUMENT_MEMBER = "word/document.xml"
STYLES_MEMBER = "word/styles.xml"
CHANGED_MEMBERS = {DOCUMENT_MEMBER, STYLES_MEMBER}
EXPECTED_FIGURES = set(range(9, 72))
EXPECTED_INTRO_FIGURES = set(range(1, 9))
EXPECTED_INTRO_NARRATIVE_FIGURES = set(range(5, 9))
INTRO_CONTINUATION_FIGURES = {2, 3, 4, 7, 8}
CAPTION_STYLE_ID = "Caption"
FONT_NAME = "Times New Roman"
FONT_SIZE_HALF_POINTS = "24"  # 12 pt
LINE_SPACING_TWIPS = "360"  # 1.5 lines with lineRule=auto
NARRATIVE_BOOKMARK_PREFIX = "ThesisFigureNarrative_"

FIGURE_CAPTION_RE = re.compile(r"^Figura\s+(\d+)\.\s*(.*)$", re.DOTALL)
TABLE_CAPTION_RE = re.compile(r"^Tabla\s+(\d+)\.", re.IGNORECASE)
ANY_FIGURE_START_RE = re.compile(r"^Figura\s+\d+\.", re.IGNORECASE)
GROUP_NOTE_RE = re.compile(r"^Nota de las Figuras\b", re.IGNORECASE)


def w(tag: str) -> str:
    return W + tag


def normalized_text(value: object, *, label: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} debe ser str; se recibió {type(value).__name__}")
    cleaned = re.sub(r"\s+", " ", value).strip()
    if not cleaned:
        raise ValueError(f"{label} está vacío")
    return cleaned


def validate_content() -> tuple[dict[int, str], dict[int, str]]:
    """Validate and normalize the two externally maintained content maps."""
    narrative_keys = set(FIGURE_RESULT_NARRATIVES)
    caption_keys = set(FIGURE_CAPTION_OVERRIDES)
    if narrative_keys != EXPECTED_FIGURES:
        raise ValueError(
            "FIGURE_RESULT_NARRATIVES debe contener exactamente las claves "
            f"9..71; faltan={sorted(EXPECTED_FIGURES - narrative_keys)}, "
            f"sobran={sorted(narrative_keys - EXPECTED_FIGURES, key=str)}"
        )
    if caption_keys != EXPECTED_FIGURES:
        raise ValueError(
            "FIGURE_CAPTION_OVERRIDES debe contener exactamente las claves "
            f"9..71; faltan={sorted(EXPECTED_FIGURES - caption_keys)}, "
            f"sobran={sorted(caption_keys - EXPECTED_FIGURES, key=str)}"
        )
    if any(type(key) is not int for key in narrative_keys | caption_keys):
        raise TypeError("Las claves de ambos diccionarios deben ser enteros")

    narratives: dict[int, str] = {}
    captions: dict[int, str] = {}
    for number in sorted(EXPECTED_FIGURES):
        narrative = normalized_text(
            FIGURE_RESULT_NARRATIVES[number],
            label=f"FIGURE_RESULT_NARRATIVES[{number}]",
        )
        if ANY_FIGURE_START_RE.match(narrative):
            raise ValueError(
                f"La narrativa {number} no debe comenzar con 'Figura N.': "
                "ese patrón está reservado para el pie y el mapa de páginas"
            )
        caption = normalized_text(
            FIGURE_CAPTION_OVERRIDES[number],
            label=f"FIGURE_CAPTION_OVERRIDES[{number}]",
        )
        match = FIGURE_CAPTION_RE.fullmatch(caption)
        if match is None or int(match.group(1)) != number or not match.group(2).strip():
            raise ValueError(
                f"El pie {number} debe comenzar exactamente con 'Figura {number}.' "
                "y contener una descripción"
            )
        narratives[number] = narrative
        captions[number] = caption
    return narratives, captions
def validate_intro_content() -> tuple[dict[int, str], dict[int, str]]:
    """Validate and normalize the Figure 1--8 caption and narrative maps."""
    narrative_keys = set(INTRO_FIGURE_RESULT_NARRATIVES)
    caption_keys = set(INTRO_FIGURE_CAPTION_OVERRIDES)
    if narrative_keys != EXPECTED_INTRO_NARRATIVE_FIGURES:
        raise ValueError(
            "INTRO_FIGURE_RESULT_NARRATIVES debe contener exactamente 5..8; "
            f"faltan={sorted(EXPECTED_INTRO_NARRATIVE_FIGURES - narrative_keys)}, "
            f"sobran={sorted(narrative_keys - EXPECTED_INTRO_NARRATIVE_FIGURES, key=str)}"
        )
    if caption_keys != EXPECTED_INTRO_FIGURES:
        raise ValueError(
            "INTRO_FIGURE_CAPTION_OVERRIDES debe contener exactamente 1..8; "
            f"faltan={sorted(EXPECTED_INTRO_FIGURES - caption_keys)}, "
            f"sobran={sorted(caption_keys - EXPECTED_INTRO_FIGURES, key=str)}"
        )
    if any(type(key) is not int for key in narrative_keys | caption_keys):
        raise TypeError(
            "Las claves de los diccionarios de figuras introductorias deben ser enteros"
        )

    narratives: dict[int, str] = {}
    captions: dict[int, str] = {}
    for number in sorted(EXPECTED_INTRO_NARRATIVE_FIGURES):
        narrative = normalized_text(
            INTRO_FIGURE_RESULT_NARRATIVES[number],
            label=f"INTRO_FIGURE_RESULT_NARRATIVES[{number}]",
        )
        if ANY_FIGURE_START_RE.match(narrative):
            raise ValueError(
                f"La narrativa introductoria {number} no debe comenzar con 'Figura N.': "
                "ese patrón está reservado para el pie y el mapa de páginas"
            )
        narratives[number] = narrative
    for number in sorted(EXPECTED_INTRO_FIGURES):
        caption = normalized_text(
            INTRO_FIGURE_CAPTION_OVERRIDES[number],
            label=f"INTRO_FIGURE_CAPTION_OVERRIDES[{number}]",
        )
        match = FIGURE_CAPTION_RE.fullmatch(caption)
        if match is None or int(match.group(1)) != number or not match.group(2).strip():
            raise ValueError(
                f"El pie introductorio {number} debe comenzar exactamente con "
                f"'Figura {number}.' y contener una descripción"
            )
        captions[number] = caption
    return narratives, captions




def paragraph_text(paragraph: etree._Element) -> str:
    return "".join(node.text or "" for node in paragraph.iter(w("t"))).strip()


def paragraph_style(paragraph: etree._Element) -> str | None:
    node = paragraph.find("./w:pPr/w:pStyle", namespaces=NS)
    return node.get(w("val")) if node is not None else None


def get_or_add(parent: etree._Element, tag: str) -> etree._Element:
    child = parent.find(w(tag))
    if child is None:
        child = etree.SubElement(parent, w(tag))
    return child


def get_or_add_ppr(paragraph: etree._Element) -> etree._Element:
    ppr = paragraph.find(w("pPr"))
    if ppr is None:
        ppr = etree.Element(w("pPr"))
        paragraph.insert(0, ppr)
    return ppr


def apply_paragraph_format(
    paragraph: etree._Element,
    *,
    style_id: str,
    justify: bool,
    keep_next: bool,
    keep_lines: bool,
    first_line_twips: str | None,
    after_twips: str,
) -> None:
    ppr = get_or_add_ppr(paragraph)
    pstyle = get_or_add(ppr, "pStyle")
    pstyle.set(w("val"), style_id)

    spacing = get_or_add(ppr, "spacing")
    spacing.set(w("line"), LINE_SPACING_TWIPS)
    spacing.set(w("lineRule"), "auto")
    spacing.set(w("after"), after_twips)

    if justify:
        jc = get_or_add(ppr, "jc")
        jc.set(w("val"), "both")

    for tag, enabled in (("keepNext", keep_next), ("keepLines", keep_lines)):
        existing = ppr.find(w(tag))
        if enabled:
            node = existing if existing is not None else etree.SubElement(ppr, w(tag))
            node.set(w("val"), "1")
        elif existing is not None:
            ppr.remove(existing)

    ind = ppr.find(w("ind"))
    if first_line_twips is not None:
        ind = ind if ind is not None else etree.SubElement(ppr, w("ind"))
        ind.set(w("firstLine"), first_line_twips)
    elif ind is not None:
        ind.attrib.pop(w("firstLine"), None)


def apply_run_font(run: etree._Element, *, bold: bool | None = None) -> None:
    rpr = run.find(w("rPr"))
    if rpr is None:
        rpr = etree.Element(w("rPr"))
        run.insert(0, rpr)

    rfonts = get_or_add(rpr, "rFonts")
    for attribute in ("ascii", "hAnsi", "eastAsia", "cs"):
        rfonts.set(w(attribute), FONT_NAME)

    size = get_or_add(rpr, "sz")
    size.set(w("val"), FONT_SIZE_HALF_POINTS)
    size_cs = get_or_add(rpr, "szCs")
    size_cs.set(w("val"), FONT_SIZE_HALF_POINTS)

    if bold is not None:
        for tag in ("b", "bCs"):
            node = rpr.find(w(tag))
            if bold:
                node = node if node is not None else etree.SubElement(rpr, w(tag))
                node.set(w("val"), "1")
            elif node is not None:
                rpr.remove(node)


def append_text_run(
    paragraph: etree._Element,
    text: str,
    *,
    bold: bool = False,
) -> etree._Element:
    run = etree.SubElement(paragraph, w("r"))
    apply_run_font(run, bold=bold)
    text_node = etree.SubElement(run, w("t"))
    text_node.set(XML_SPACE, "preserve")
    text_node.text = text
    return run


def clear_paragraph_content(paragraph: etree._Element) -> None:
    """Keep paragraph properties but remove runs, hyperlinks and old fields."""
    for child in list(paragraph):
        if child.tag != w("pPr"):
            paragraph.remove(child)


def replace_caption(paragraph: etree._Element, number: int, text: str) -> None:
    match = FIGURE_CAPTION_RE.fullmatch(text)
    if match is None or int(match.group(1)) != number:
        raise AssertionError(f"Pie inválido para Figura {number}: {text[:80]}")
    body = match.group(2).strip()
    clear_paragraph_content(paragraph)
    apply_paragraph_format(
        paragraph,
        style_id=CAPTION_STYLE_ID,
        justify=True,
        keep_next=False,
        keep_lines=True,
        first_line_twips=None,
        after_twips="160",  # 8 pt
    )
    append_text_run(paragraph, f"Figura {number}. ", bold=True)
    append_text_run(paragraph, body, bold=False)


def bookmark_id_map(root: etree._Element) -> tuple[dict[str, str], int]:
    mapping: dict[str, str] = {}
    maximum = 0
    for node in root.xpath(".//w:bookmarkStart", namespaces=NS):
        name = node.get(w("name"))
        identifier = node.get(w("id"))
        if name and identifier is not None:
            mapping[name] = identifier
        try:
            maximum = max(maximum, int(identifier or "0"))
        except ValueError:
            continue
    return mapping, maximum


def has_bookmark(paragraph: etree._Element, name: str) -> bool:
    return bool(
        paragraph.xpath(
            ".//w:bookmarkStart[@w:name=$name]",
            namespaces=NS,
            name=name,
        )
    )


def set_narrative_paragraph(
    paragraph: etree._Element,
    *,
    number: int,
    text: str,
    bookmark_id: str,
) -> None:
    clear_paragraph_content(paragraph)
    apply_paragraph_format(
        paragraph,
        style_id="Normal",
        justify=True,
        keep_next=True,
        keep_lines=False,
        first_line_twips="504",  # 0.35 in
        after_twips="120",  # 6 pt
    )
    bookmark_name = f"{NARRATIVE_BOOKMARK_PREFIX}{number}"
    start = etree.SubElement(paragraph, w("bookmarkStart"))
    start.set(w("id"), bookmark_id)
    start.set(w("name"), bookmark_name)
    append_text_run(paragraph, text)
    end = etree.SubElement(paragraph, w("bookmarkEnd"))
    end.set(w("id"), bookmark_id)


def is_image_paragraph(paragraph: etree._Element | None) -> bool:
    if paragraph is None or paragraph.tag != w("p"):
        return False
    return bool(
        paragraph.xpath(".//w:drawing | .//w:pict", namespaces=NS)
    )

def image_before_caption(caption: etree._Element, number: int) -> etree._Element:
    image = caption.getprevious()
    if not is_image_paragraph(image):
        previous_text = paragraph_text(image)[:100] if image is not None else ""
        raise RuntimeError(
            f"Figura {number}: el párrafo inmediatamente anterior al pie "
            f"no contiene una imagen; texto previo={previous_text!r}"
        )
    if image.getparent() is not caption.getparent():
        raise RuntimeError(f"Figura {number}: imagen y pie no son hermanos XML")
    return image


def intro_image_before_caption(
    caption: etree._Element, number: int
) -> etree._Element:
    """Locate an intro image and remove only intervening empty paragraphs."""
    candidate = caption.getprevious()
    empty_gap: list[etree._Element] = []
    while (
        candidate is not None
        and candidate.tag == w("p")
        and not is_image_paragraph(candidate)
        and not paragraph_text(candidate)
    ):
        empty_gap.append(candidate)
        candidate = candidate.getprevious()
    if not is_image_paragraph(candidate):
        previous_text = paragraph_text(candidate)[:100] if candidate is not None else ""
        raise RuntimeError(
            f"Figura introductoria {number}: no se encontró una imagen antes del pie; "
            f"texto previo={previous_text!r}"
        )
    parent = caption.getparent()
    if parent is None or candidate.getparent() is not parent:
        raise RuntimeError(f"Figura introductoria {number}: imagen y pie no son hermanos XML")
    for paragraph in empty_gap:
        parent.remove(paragraph)
    return candidate


def intro_figure2_required_after(image: etree._Element) -> int:
    """Return the paragraph spacing needed below the floating Figure 2 group."""
    anchors = image.xpath(".//*[local-name()=\"anchor\"]")
    group_anchors = [
        anchor
        for anchor in anchors
        if len(anchor.xpath(".//*[local-name()=\"blip\"]")) == 2
    ]
    if len(group_anchors) != 1:
        raise RuntimeError(
            "Figura introductoria 2: se esperaba un grupo flotante con dos imágenes; "
            f"se encontraron {len(group_anchors)}"
        )
    anchor = group_anchors[0]
    pos_offsets = anchor.xpath(
        "./*[local-name()=\"positionV\"]/*[local-name()=\"posOffset\"]/text()"
    )
    extent_ys = anchor.xpath("./*[local-name()=\"extent\"]/@cy")
    effect_bottoms = anchor.xpath("./*[local-name()=\"effectExtent\"]/@b")
    if len(pos_offsets) != 1 or len(extent_ys) != 1 or len(effect_bottoms) != 1:
        raise RuntimeError(
            "Figura introductoria 2: geometría flotante incompleta o ambigua"
        )
    pos_offset = int(pos_offsets[0])
    extent_y = int(extent_ys[0])
    effect_bottom = int(effect_bottoms[0])
    distance_bottom = int(anchor.get("distB", "0"))
    total_twips = math.ceil(
        (pos_offset + extent_y + effect_bottom + distance_bottom) / 635
    )
    required_after = max(
        0, total_twips - int(LINE_SPACING_TWIPS) + 160
    )
    return int(math.ceil(required_after / 100) * 100)


def reserve_intro_figure2_floating_height(image: etree._Element) -> int:
    required_after = intro_figure2_required_after(image)
    ppr = get_or_add_ppr(image)
    spacing = get_or_add(ppr, "spacing")
    spacing.set(w("line"), LINE_SPACING_TWIPS)
    spacing.set(w("lineRule"), "auto")
    spacing.set(w("after"), str(required_after))
    keep_next = get_or_add(ppr, "keepNext")
    keep_next.set(w("val"), "1")
    return required_after


def validate_intro_figure2_floating_height(image: etree._Element) -> None:
    required_after = intro_figure2_required_after(image)
    spacing = image.find("./w:pPr/w:spacing", namespaces=NS)
    keep_next = image.find("./w:pPr/w:keepNext", namespaces=NS)
    if spacing is None or int(spacing.get(w("after"), "0")) < required_after:
        raise RuntimeError(
            "Figura introductoria 2: no se reservó la altura del grupo flotante"
        )
    if keep_next is None or keep_next.get(w("val"), "1") in ("0", "false", "off"):
        raise RuntimeError(
            "Figura introductoria 2: el párrafo de imagen no conserva el pie"
        )


def upsert_marked_narrative(
    image: etree._Element,
    *,
    number: int,
    text: str,
    bookmark_ids: dict[str, str],
    next_bookmark_id: int,
) -> tuple[int, bool]:
    bookmark_name = f"{NARRATIVE_BOOKMARK_PREFIX}{number}"
    previous = image.getprevious()
    if (
        previous is not None
        and previous.tag == w("p")
        and has_bookmark(previous, bookmark_name)
    ):
        narrative = previous
        bookmark_id = bookmark_ids.get(bookmark_name)
        if bookmark_id is None:
            raise RuntimeError(f"Figura {number}: marcador narrativo sin id")
        inserted = False
    else:
        next_bookmark_id += 1
        bookmark_id = str(next_bookmark_id)
        narrative = etree.Element(w("p"))
        image.addprevious(narrative)
        bookmark_ids[bookmark_name] = bookmark_id
        inserted = True
    set_narrative_paragraph(
        narrative,
        number=number,
        text=text,
        bookmark_id=bookmark_id,
    )
    return next_bookmark_id, inserted


def locate_thesis_captions(root: etree._Element) -> dict[int, etree._Element]:
    found: dict[int, etree._Element] = {}
    duplicates: list[int] = []
    for paragraph in root.xpath(".//w:p", namespaces=NS):
        if paragraph_style(paragraph) != CAPTION_STYLE_ID:
            continue
        match = FIGURE_CAPTION_RE.fullmatch(paragraph_text(paragraph))
        if match is None:
            continue
        number = int(match.group(1))
        if number not in EXPECTED_FIGURES:
            continue
        if number in found:
            duplicates.append(number)
        else:
            found[number] = paragraph
    if duplicates:
        raise RuntimeError(f"Pies duplicados para Figuras: {sorted(set(duplicates))}")
    if set(found) != EXPECTED_FIGURES:
        raise RuntimeError(
            "No se localizaron exactamente los pies Figura 9..71; "
            f"faltan={sorted(EXPECTED_FIGURES - set(found))}, "
            f"sobran={sorted(set(found) - EXPECTED_FIGURES)}"
        )
    return found
def locate_intro_captions(root: etree._Element) -> dict[int, etree._Element]:
    found: dict[int, etree._Element] = {}
    duplicates: list[int] = []
    for paragraph in root.xpath(".//w:p", namespaces=NS):
        if paragraph_style(paragraph) != CAPTION_STYLE_ID:
            continue
        match = FIGURE_CAPTION_RE.fullmatch(paragraph_text(paragraph))
        if match is None:
            continue
        number = int(match.group(1))
        if number not in EXPECTED_INTRO_FIGURES:
            continue
        if number in found:
            duplicates.append(number)
        else:
            found[number] = paragraph
    if duplicates:
        raise RuntimeError(
            f"Pies duplicados para Figuras 1..8: {sorted(set(duplicates))}"
        )
    if set(found) != EXPECTED_INTRO_FIGURES:
        raise RuntimeError(
            "No se localizaron exactamente los pies Figura 1..8; "
            f"faltan={sorted(EXPECTED_INTRO_FIGURES - set(found))}, "
            f"sobran={sorted(set(found) - EXPECTED_INTRO_FIGURES)}"
        )
    return found

def intro_caption_continuation(
    caption: etree._Element, number: int
) -> etree._Element | None:
    candidate = caption.getnext()
    if (
        candidate is None
        or candidate.tag != w("p")
        or paragraph_style(candidate) != CAPTION_STYLE_ID
    ):
        return None
    text = paragraph_text(candidate)
    if ANY_FIGURE_START_RE.match(text) or TABLE_CAPTION_RE.match(text):
        return None
    if not text:
        raise RuntimeError(f"Figura {number}: continuación Caption inmediata vacía")
    if is_image_paragraph(candidate):
        raise RuntimeError(f"Figura {number}: la continuación Caption contiene una imagen")
    return candidate


def remove_intro_caption_continuations(
    captions: dict[int, etree._Element],
) -> int:
    removed = 0
    for number in sorted(INTRO_CONTINUATION_FIGURES):
        continuation = intro_caption_continuation(captions[number], number)
        if continuation is None:
            continue
        parent = continuation.getparent()
        if parent is None:
            raise RuntimeError(f"Figura {number}: continuación Caption sin padre XML")
        parent.remove(continuation)
        removed += 1
    return removed


def remove_group_notes(root: etree._Element) -> int:
    notes = [
        paragraph
        for paragraph in root.xpath(".//w:p", namespaces=NS)
        if GROUP_NOTE_RE.match(paragraph_text(paragraph))
    ]
    for paragraph in notes:
        parent = paragraph.getparent()
        if parent is None:
            raise RuntimeError("Se encontró una nota sin elemento padre")
        parent.remove(paragraph)
    return len(notes)


def normalize_caption_paragraph(paragraph: etree._Element) -> None:
    apply_paragraph_format(
        paragraph,
        style_id=CAPTION_STYLE_ID,
        justify=True,
        keep_next=False,
        keep_lines=True,
        first_line_twips=None,
        after_twips="160",
    )
    for run in paragraph.xpath(".//w:r", namespaces=NS):
        apply_run_font(run)


def normalize_all_captions(root: etree._Element) -> int:
    captions = [
        paragraph
        for paragraph in root.xpath(".//w:p", namespaces=NS)
        if paragraph_style(paragraph) == CAPTION_STYLE_ID
    ]
    for paragraph in captions:
        normalize_caption_paragraph(paragraph)
    return len(captions)


def normalize_caption_style(styles_root: etree._Element) -> None:
    matches = styles_root.xpath(
        ".//w:style[@w:type='paragraph' and @w:styleId=$style_id]",
        namespaces=NS,
        style_id=CAPTION_STYLE_ID,
    )
    if len(matches) != 1:
        raise RuntimeError(
            f"Se esperaba un estilo Caption; se encontraron {len(matches)}"
        )
    style = matches[0]
    ppr = style.find(w("pPr"))
    if ppr is None:
        ppr = etree.SubElement(style, w("pPr"))
    spacing = get_or_add(ppr, "spacing")
    spacing.set(w("line"), LINE_SPACING_TWIPS)
    spacing.set(w("lineRule"), "auto")
    spacing.set(w("after"), "160")
    jc = get_or_add(ppr, "jc")
    jc.set(w("val"), "both")

    rpr = style.find(w("rPr"))
    if rpr is None:
        rpr = etree.SubElement(style, w("rPr"))
    synthetic_run = etree.Element(w("r"))
    synthetic_run.append(rpr)
    apply_run_font(synthetic_run)
    # apply_run_font reuses the same rPr element; return it to the style.
    style.append(synthetic_run.remove(rpr) if False else rpr)


def parse_xml(data: bytes, *, member: str) -> etree._Element:
    parser = etree.XMLParser(
        remove_blank_text=False,
        resolve_entities=False,
        no_network=True,
        huge_tree=True,
    )
    try:
        return etree.fromstring(data, parser=parser)
    except etree.XMLSyntaxError as error:
        raise RuntimeError(f"XML inválido en {member}: {error}") from error


def serialize_xml(root: etree._Element) -> bytes:
    return etree.tostring(
        root,
        encoding="UTF-8",
        xml_declaration=True,
        standalone=True,
    )


def patch_document(
    document_data: bytes,
    styles_data: bytes,
    narratives: dict[int, str],
    captions: dict[int, str],
    intro_narratives: dict[int, str],
    intro_captions: dict[int, str],
) -> tuple[bytes, bytes, dict[str, int]]:
    root = parse_xml(document_data, member=DOCUMENT_MEMBER)
    styles_root = parse_xml(styles_data, member=STYLES_MEMBER)
    located = locate_thesis_captions(root)
    intro_located = locate_intro_captions(root)
    continuations_removed = remove_intro_caption_continuations(intro_located)
    bookmark_ids, next_bookmark_id = bookmark_id_map(root)

    intro_inserted = 0
    intro_updated = 0
    intro_figure2_reserved_twips = 0
    for number in sorted(EXPECTED_INTRO_FIGURES):
        caption = intro_located[number]
        image = intro_image_before_caption(caption, number)
        if number == 2:
            intro_figure2_reserved_twips = reserve_intro_figure2_floating_height(image)
        if number in intro_narratives:
            next_bookmark_id, was_inserted = upsert_marked_narrative(
                image,
                number=number,
                text=intro_narratives[number],
                bookmark_ids=bookmark_ids,
                next_bookmark_id=next_bookmark_id,
            )
            if was_inserted:
                intro_inserted += 1
            else:
                intro_updated += 1
        replace_caption(caption, number, intro_captions[number])

    inserted = 0
    updated = 0
    for number in sorted(EXPECTED_FIGURES):
        caption = located[number]
        image = image_before_caption(caption, number)
        next_bookmark_id, was_inserted = upsert_marked_narrative(
            image,
            number=number,
            text=narratives[number],
            bookmark_ids=bookmark_ids,
            next_bookmark_id=next_bookmark_id,
        )
        if was_inserted:
            inserted += 1
        else:
            updated += 1
        replace_caption(caption, number, captions[number])

    notes_removed = remove_group_notes(root)
    captions_normalized = normalize_all_captions(root)
    normalize_caption_style(styles_root)
    validate_patched_structure(
        root,
        narratives,
        captions,
        intro_narratives,
        intro_captions,
    )
    return (
        serialize_xml(root),
        serialize_xml(styles_root),
        {
            "intro_narratives_inserted": intro_inserted,
            "intro_narratives_updated": intro_updated,
            "intro_figure_captions_replaced": len(intro_located),
            "intro_figure2_reserved_twips": intro_figure2_reserved_twips,
            "intro_caption_continuations_removed": continuations_removed,
            "narratives_inserted": inserted,
            "narratives_updated": updated,
            "figure_captions_replaced": len(located),
            "group_notes_removed": notes_removed,
            "caption_paragraphs_normalized": captions_normalized,
        },
    )




def validate_run_font(run: etree._Element, *, context: str) -> None:
    rpr = run.find(w("rPr"))
    if rpr is None:
        raise RuntimeError(f"{context}: run sin rPr")
    rfonts = rpr.find(w("rFonts"))
    if rfonts is None or any(
        rfonts.get(w(attribute)) != FONT_NAME
        for attribute in ("ascii", "hAnsi", "eastAsia", "cs")
    ):
        raise RuntimeError(f"{context}: fuente distinta de {FONT_NAME}")
    for tag in ("sz", "szCs"):
        size = rpr.find(w(tag))
        if size is None or size.get(w("val")) != FONT_SIZE_HALF_POINTS:
            raise RuntimeError(f"{context}: tamaño distinto de 12 pt")


def validate_spacing(paragraph: etree._Element, *, context: str) -> None:
    spacing = paragraph.find("./w:pPr/w:spacing", namespaces=NS)
    if (
        spacing is None
        or spacing.get(w("line")) != LINE_SPACING_TWIPS
        or spacing.get(w("lineRule")) != "auto"
    ):
        raise RuntimeError(f"{context}: interlineado distinto de 1,5")


def validate_marked_narrative(
    image: etree._Element,
    *,
    number: int,
    expected_text: str,
) -> None:
    narrative = image.getprevious()
    bookmark_name = f"{NARRATIVE_BOOKMARK_PREFIX}{number}"
    if (
        narrative is None
        or narrative.tag != w("p")
        or not has_bookmark(narrative, bookmark_name)
    ):
        raise RuntimeError(f"Figura {number}: falta narrativa marcada antes de la imagen")
    if paragraph_text(narrative) != expected_text:
        raise RuntimeError(f"Figura {number}: narrativa distinta del registro")
    if paragraph_style(narrative) != "Normal":
        raise RuntimeError(f"Figura {number}: narrativa sin estilo Normal")
    keep_next = narrative.find("./w:pPr/w:keepNext", namespaces=NS)
    jc = narrative.find("./w:pPr/w:jc", namespaces=NS)
    if keep_next is None or keep_next.get(w("val"), "1") in ("0", "false", "off"):
        raise RuntimeError(f"Figura {number}: narrativa sin keepNext")
    if jc is None or jc.get(w("val")) != "both":
        raise RuntimeError(f"Figura {number}: narrativa sin justificación")
    validate_spacing(narrative, context=f"Figura {number}, narrativa")
    for run in narrative.xpath(".//w:r[w:t]", namespaces=NS):
        validate_run_font(run, context=f"Figura {number}, narrativa")

def validate_patched_structure(
    root: etree._Element,
    narratives: dict[int, str],
    captions: dict[int, str],
    intro_narratives: dict[int, str],
    intro_captions: dict[int, str],
) -> None:
    located = locate_thesis_captions(root)
    intro_located = locate_intro_captions(root)
    for number in sorted(EXPECTED_INTRO_FIGURES):
        caption = intro_located[number]
        if paragraph_text(caption) != intro_captions[number]:
            raise RuntimeError(
                f"Figura introductoria {number}: pie distinto del registro"
            )
        image = intro_image_before_caption(caption, number)
        if number == 2:
            validate_intro_figure2_floating_height(image)
        if number in intro_narratives:
            validate_marked_narrative(
                image,
                number=number,
                expected_text=intro_narratives[number],
            )
    for number in sorted(INTRO_CONTINUATION_FIGURES):
        if intro_caption_continuation(intro_located[number], number) is not None:
            raise RuntimeError(f"Figura introductoria {number}: persiste continuación")
    residual_notes = [
        paragraph_text(paragraph)
        for paragraph in root.xpath(".//w:p", namespaces=NS)
        if GROUP_NOTE_RE.match(paragraph_text(paragraph))
    ]
    if residual_notes:
        raise RuntimeError(f"Persisten notas grupales: {residual_notes[:3]}")

    for number in sorted(EXPECTED_FIGURES):
        caption = located[number]
        if paragraph_text(caption) != captions[number]:
            raise RuntimeError(f"Figura {number}: el pie final no coincide con el registro")
        image = caption.getprevious()
        if not is_image_paragraph(image):
            raise RuntimeError(f"Figura {number}: falta imagen inmediatamente antes del pie")
        narrative = image.getprevious()
        bookmark_name = f"{NARRATIVE_BOOKMARK_PREFIX}{number}"
        if (
            narrative is None
            or narrative.tag != w("p")
            or not has_bookmark(narrative, bookmark_name)
        ):
            raise RuntimeError(f"Figura {number}: falta narrativa marcada antes de la imagen")
        if paragraph_text(narrative) != narratives[number]:
            raise RuntimeError(f"Figura {number}: narrativa distinta del registro")
        if paragraph_style(narrative) != "Normal":
            raise RuntimeError(f"Figura {number}: narrativa sin estilo Normal")
        keep_next = narrative.find("./w:pPr/w:keepNext", namespaces=NS)
        jc = narrative.find("./w:pPr/w:jc", namespaces=NS)
        if keep_next is None or keep_next.get(w("val"), "1") in ("0", "false", "off"):
            raise RuntimeError(f"Figura {number}: narrativa sin keepNext")
        if jc is None or jc.get(w("val")) != "both":
            raise RuntimeError(f"Figura {number}: narrativa sin justificación")
        validate_spacing(narrative, context=f"Figura {number}, narrativa")
        for run in narrative.xpath(".//w:r[w:t]", namespaces=NS):
            validate_run_font(run, context=f"Figura {number}, narrativa")

    for position, paragraph in enumerate(
        paragraph
        for paragraph in root.xpath(".//w:p", namespaces=NS)
        if paragraph_style(paragraph) == CAPTION_STYLE_ID
    ):
        validate_spacing(paragraph, context=f"Caption {position + 1}")
        for run in paragraph.xpath(".//w:r[w:t]", namespaces=NS):
            validate_run_font(run, context=f"Caption {position + 1}")


def member_sha256(archive: zipfile.ZipFile, name: str) -> str:
    digest = hashlib.sha256()
    with archive.open(name) as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_patched_package(
    source: Path,
    temporary: Path,
    document_data: bytes,
    styles_data: bytes,
) -> tuple[int, int]:
    with zipfile.ZipFile(source, "r") as input_zip:
        infos = input_zip.infolist()
        names = [info.filename for info in infos]
        if len(names) != len(set(names)):
            raise RuntimeError("El DOCX contiene nombres ZIP duplicados")
        for required in CHANGED_MEMBERS:
            if required not in names:
                raise RuntimeError(f"Falta el miembro DOCX requerido: {required}")
        archive_comment = input_zip.comment
        media_names = [name for name in names if name.startswith("word/media/")]
        media_before = {
            name: member_sha256(input_zip, name) for name in media_names
        }

        with zipfile.ZipFile(
            temporary, "w", allowZip64=True, strict_timestamps=False
        ) as output_zip:
            output_zip.comment = archive_comment
            for info in infos:
                if info.filename == DOCUMENT_MEMBER:
                    data = document_data
                elif info.filename == STYLES_MEMBER:
                    data = styles_data
                else:
                    data = input_zip.read(info.filename)
                output_zip.writestr(info, data)

    with zipfile.ZipFile(temporary, "r") as output_zip:
        bad_member = output_zip.testzip()
        if bad_member is not None:
            raise RuntimeError(f"Miembro dañado después del parche: {bad_member}")
        if output_zip.namelist() != names:
            raise RuntimeError("La lista u orden de miembros ZIP cambió")
        media_after = {
            name: member_sha256(output_zip, name) for name in media_names
        }
        if media_after != media_before:
            changed = sorted(
                name for name in media_names if media_after.get(name) != media_before.get(name)
            )
            raise RuntimeError(f"Se modificaron miembros de word/media: {changed}")
        # Reparse the two edited members from the package actually written.
        final_root = parse_xml(output_zip.read(DOCUMENT_MEMBER), member=DOCUMENT_MEMBER)
        final_styles = parse_xml(output_zip.read(STYLES_MEMBER), member=STYLES_MEMBER)
        if final_root.tag != w("document") or final_styles.tag != w("styles"):
            raise RuntimeError("Los XML OOXML editados no conservaron sus raíces esperadas")
    return len(names), len(media_names)


def refine(source: Path, output: Path, *, force: bool) -> dict[str, object]:
    source = source.resolve()
    output = output.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    same_path = source == output
    if output.exists() and not force:
        raise FileExistsError(
            f"La salida ya existe: {output}. Use --force para reemplazarla de forma atómica."
        )

    narratives, captions = validate_content()
    intro_narratives, intro_captions = validate_intro_content()
    with zipfile.ZipFile(source, "r") as archive:
        bad_member = archive.testzip()
        if bad_member is not None:
            raise RuntimeError(f"DOCX de entrada dañado: {bad_member}")
        document_data = archive.read(DOCUMENT_MEMBER)
        styles_data = archive.read(STYLES_MEMBER)

    patched_document, patched_styles, audit = patch_document(
        document_data,
        styles_data,
        narratives,
        captions,
        intro_narratives,
        intro_captions,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=output.name + ".",
        suffix=".building",
        dir=output.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        members, media_members = write_patched_package(
            source, temporary, patched_document, patched_styles
        )
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    result: dict[str, object] = {
        "status": "complete",
        "source": str(source),
        "output": str(output),
        "in_place": same_path,
        "sha256": file_sha256(output),
        "bytes": output.stat().st_size,
        "zip_members": members,
        "media_members_verified_unchanged": media_members,
        **audit,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="DOCX de entrada")
    parser.add_argument("output", type=Path, help="DOCX de salida")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Permitir reemplazo atómico si la salida ya existe, incluida salida in-place",
    )
    args = parser.parse_args()
    result = refine(args.input, args.output, force=args.force)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
