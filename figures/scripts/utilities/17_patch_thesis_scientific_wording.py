#!/usr/bin/env python3
"""Make nominal-versus-adjusted statistical wording consistent in the thesis."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
import zipfile
from pathlib import Path

from lxml import etree


W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
W = f"{{{W_NS}}}"
NS = {"w": W_NS}
XML_SPACE = "{http://www.w3.org/XML/1998/namespace}space"
DOCUMENT_MEMBER = "word/document.xml"
FONT = "Times New Roman"


SUBSTITUTIONS = [
    (
        "Ninguna comparación por pares ajustada por Holm alcanzó p<0,05.",
        "Los contrastes exactos Agua–Sacarosa (p nominal no ajustado=0,0357) y "
        "Sacarosa–Alulosa (p nominal no ajustado=0,0238) fueron positivos antes "
        "del ajuste, pero ambos quedaron en p=0,0714 tras Holm; Agua–Alulosa fue "
        "p=0,1429 antes y después del ajuste.",
    ),
    (
        "el contraste exacto Sacarosa–Alulosa en ARC fue significativo "
        "(p=0,0411), con una proporción mayor en Alulosa; el ANOVA regional fue "
        "p=0,133 y el contraste no se mantuvo tras la corrección conjunta por "
        "multiplicidad.",
        "el contraste exacto Sacarosa–Alulosa en ARC fue nominalmente "
        "significativo y no ajustado (p=0,0411), con una proporción mayor en "
        "Alulosa. El ANOVA de ARC fue p=0,133 y su valor p ajustado por Holm dentro de "
        "la familia de cuatro ANOVA regionales fue 0,531; este ajuste corresponde "
        "al ANOVA regional y no al contraste por pares.",
    ),
    (
        "la clasificación morfométrica microglial identificó comparaciones "
        "positivas en ARC (Agua–Sacarosa p=0,0476; Sacarosa–Alulosa p=0,0260) y "
        "en eminencia media (Agua–Sacarosa p=0,0357)",
        "la clasificación morfométrica microglial identificó comparaciones "
        "nominales no ajustadas en ARC (Agua–Sacarosa p=0,0476; "
        "Sacarosa–Alulosa p=0,0260) y en eminencia media (Agua–Sacarosa "
        "p=0,0357)",
    ),
    (
        "No Holm-adjusted pairwise comparison reached p<0.05",
        "The exact water–sucrose (nominal unadjusted p=0.0357) and "
        "sucrose–allulose (nominal unadjusted p=0.0238) contrasts were positive "
        "before correction but both became p=0.0714 after Holm adjustment; the "
        "water–allulose contrast was p=0.1429 before and after adjustment",
    ),
    (
        "the exact sucrose–allulose c-FOS/DAPI contrast was significant "
        "(p=0.0411), with a higher proportion under allulose; the regional ANOVA "
        "was p=0.133 and the contrast did not remain significant after joint "
        "multiplicity correction.",
        "the exact sucrose–allulose c-FOS/DAPI contrast was nominally significant "
        "and unadjusted (p=0.0411), with a higher proportion under allulose. The "
        "ARC ANOVA was p=0.133 and its Holm-adjusted p value within the four "
        "regional ANOVAs was 0.531; this correction applies to the regional ANOVA, "
        "not to the pairwise contrast.",
    ),
    (
        "yielded positive pairwise comparisons in ARC (water–sucrose p=0.0476; "
        "sucrose–allulose p=0.0260) and median eminence (water–sucrose p=0.0357)",
        "yielded nominal unadjusted pairwise signals in ARC (water–sucrose "
        "p=0.0476; sucrose–allulose p=0.0260) and median eminence "
        "(water–sucrose p=0.0357)",
    ),
    (
        "Las comparaciones exactas ajustadas por Holm fueron Agua–Sacarosa "
        "p=0,0714, Agua–Alulosa p=0,1429 y Sacarosa–Alulosa p=0,0714.",
        "Los contrastes exactos nominales no ajustados fueron Agua–Sacarosa "
        "p=0,0357, Agua–Alulosa p=0,1429 y Sacarosa–Alulosa p=0,0238; después "
        "de Holm fueron p=0,0714, p=0,1429 y p=0,0714, respectivamente.",
    ),
    (
        "la comparación exacta Sacarosa–Alulosa fue significativa (p=0,0411)",
        "la comparación exacta Sacarosa–Alulosa fue nominalmente significativa "
        "y no ajustada (p=0,0411)",
    ),
    (
        "El ANOVA global del ARC fue p=0,133 y, al considerar conjuntamente las "
        "cuatro regiones, el ajuste de Holm produjo p=0,531",
        "El ANOVA global del ARC fue p=0,133 y su valor p ajustado por Holm "
        "dentro de la familia de cuatro ANOVA regionales fue 0,531",
    ),
    (
        "produjo comparaciones positivas. En ARC, Agua fue mayor que Sacarosa "
        "(p=0,0476) y Alulosa fue mayor que Sacarosa (p=0,0260)",
        "produjo comparaciones nominales no ajustadas. En ARC, Agua fue mayor "
        "que Sacarosa (p nominal no ajustado=0,0476) y Alulosa fue mayor que "
        "Sacarosa (p nominal no ajustado=0,0260)",
    ),
    (
        "En EM también se observó una diferencia Agua–Sacarosa (p=0,0357)",
        "En EM también se observó una diferencia nominal no ajustada "
        "Agua–Sacarosa (p=0,0357)",
    ),
    (
        "el cambio de peso presentó un ANOVA positivo (p=0,049)",
        "el cambio de peso presentó un ANOVA nominalmente significativo y no "
        "ajustado (p=0,0487)",
    ),
    (
        "el resultado se conserva como una señal global positiva de esta subcohorte",
        "el resultado se conserva como una señal global nominal positiva de esta "
        "subcohorte",
    ),
    (
        "En hembras, el volumen retirado del día 6 difirió globalmente "
        "(ANOVA p=0,032)",
        "En hembras, el consumo del día 6 mostró un ANOVA nominalmente "
        "significativo y no ajustado (p=0,0315)",
    ),
    (
        "Este ANOVA positivo de volumen",
        "Este ANOVA nominal positivo de consumo",
    ),
    (
        "Alulosa presentó mayor c-FOS/DAPI que Sacarosa en la comparación exacta "
        "(p=0,0411), y el perfil espacial",
        "Alulosa presentó mayor c-FOS/DAPI que Sacarosa en el contraste exacto nominal no ajustado (p=0,0411), y el perfil espacial",
    ),
    (
        "también mostró comparaciones positivas en ARC y EM",
        "también mostró comparaciones nominales no ajustadas en ARC y EM",
    ),
    (
        "fueron significativas; en EM también lo fue Agua–Sacarosa (p=0,0357)",
        "fueron nominalmente significativas sin ajuste; en EM, Agua–Sacarosa "
        "(p=0,0357) también fue nominalmente significativa sin ajuste",
    ),
    (
        "En machos, el cambio de peso presentó ANOVA p=0,049",
        "En machos, el cambio de peso presentó un ANOVA nominal no ajustado "
        "p=0,0487",
    ),
    (
        "en hembras, el volumen retirado presentó ANOVA p=0,032",
        "en hembras, el consumo presentó un ANOVA nominal no ajustado "
        "p=0,0315",
    ),
    (
        "en la comparación exacta (p=0,0411).",
        "en el contraste exacto nominal no ajustado (p=0,0411).",
    ),

    (
        "la fracción c-FOS+/NPY+ entre núcleos NPY+ fue mayor con alulosa que "
        "con sacarosa",
        "la media de la fracción c-FOS+/NPY+ entre núcleos NPY+ fue "
        "descriptivamente mayor con alulosa que con sacarosa",
    ),
    (
        "the c-FOS+/NPY+ fraction was higher under allulose than sucrose",
        "the mean c-FOS+/NPY+ fraction was descriptively higher under allulose "
        "than under sucrose",
    ),
    (
        "El experimento de dos botellas también mostró ANOVA positivos para "
        "consumo de solución, agua, fracción de solución y alimento (todos "
        "p<0,004)",
        "El experimento de dos botellas también mostró cuatro ANOVA globales "
        "nominalmente significativos, con valores p no ajustados por multiplicidad "
        "entre desenlaces, para consumo de solución, agua, fracción de solución y "
        "alimento (todos p<0,004)",
    ),
    (
        "The two-bottle experiment also yielded positive ANOVAs for solution, "
        "water, solution share, and food intake (all p<0.004), driven by sucrose",
        "The two-bottle experiment also yielded four nominally significant global "
        "ANOVAs, with p values unadjusted for multiplicity across endpoints, for "
        "solution, water, solution share, and food intake (all p<0.004), driven by "
        "sucrose",
    ),
    (
        "En el experimento de dos botellas se obtuvieron cuatro resultados globales "
        "positivos.",
        "En el experimento de dos botellas se obtuvieron cuatro ANOVA globales "
        "nominalmente significativos, con valores p no ajustados por multiplicidad "
        "entre desenlaces.",
    ),
    (
        "los ANOVA se informan como señales globales positivas de un experimento "
        "exploratorio",
        "los ANOVA se informan como señales globales nominales no ajustadas de un "
        "experimento exploratorio",
    ),
    (
        "Los análisis por sexo y de dos botellas generaron resultados globales "
        "positivos.",
        "Los análisis por sexo y de dos botellas generaron resultados globales "
        "nominales no ajustados.",
    ),
    (
        "los cuatro ANOVA fueron positivos (todos p<0,004)",
        "los cuatro ANOVA fueron nominalmente significativos, con valores p no "
        "ajustados por multiplicidad entre desenlaces (todos p<0,004)",
    ),
    (
        "presentaron ANOVA globales positivos (todos p<0,004)",
        "presentaron cuatro ANOVA globales nominalmente significativos, con valores "
        "p no ajustados por multiplicidad entre desenlaces (todos p<0,004)",
    ),
]


GLOBAL_REPLACEMENTS = [
    ("núcleos doble positivos", "núcleos doblemente positivos"),
    ("células doble positivas", "células doblemente positivas"),
    (" y doble positivos", " y doblemente positivos"),
]


def w(tag: str) -> str:
    return W + tag


def paragraph_text(paragraph: etree._Element) -> str:
    return "".join(node.text or "" for node in paragraph.iter(w("t")))


def get_or_add(parent: etree._Element, tag: str) -> etree._Element:
    child = parent.find(w(tag))
    if child is None:
        child = etree.SubElement(parent, w(tag))
    return child


def replace_paragraph_text(paragraph: etree._Element, text: str) -> None:
    ppr = paragraph.find(w("pPr"))
    for child in list(paragraph):
        if child is not ppr:
            paragraph.remove(child)
    if ppr is None:
        ppr = etree.Element(w("pPr"))
        paragraph.insert(0, ppr)
    spacing = get_or_add(ppr, "spacing")
    spacing.set(w("line"), "360")
    spacing.set(w("lineRule"), "auto")
    jc = get_or_add(ppr, "jc")
    jc.set(w("val"), "both")
    run = etree.SubElement(paragraph, w("r"))
    rpr = etree.SubElement(run, w("rPr"))
    fonts = etree.SubElement(rpr, w("rFonts"))
    for attribute in ("ascii", "hAnsi", "eastAsia", "cs"):
        fonts.set(w(attribute), FONT)
    for tag in ("sz", "szCs"):
        node = etree.SubElement(rpr, w(tag))
        node.set(w("val"), "24")
    text_node = etree.SubElement(run, w("t"))
    text_node.set(XML_SPACE, "preserve")
    text_node.text = text


def patch_document(data: bytes) -> tuple[bytes, dict[str, int]]:
    parser = etree.XMLParser(remove_blank_text=False, no_network=True, huge_tree=True)
    root = etree.fromstring(data, parser)
    paragraphs = root.xpath(".//w:p", namespaces=NS)
    applied = 0
    already = 0
    for old, new in SUBSTITUTIONS:
        old_matches = [p for p in paragraphs if old in paragraph_text(p)]
        if len(old_matches) == 1:
            paragraph = old_matches[0]
            replace_paragraph_text(paragraph, paragraph_text(paragraph).replace(old, new))
            applied += 1
            continue
        if old_matches:
            raise RuntimeError(
                f"Sustitución ambigua ({len(old_matches)} coincidencias): {old[:80]}"
            )
        new_matches = [p for p in paragraphs if new in paragraph_text(p)]
        if len(new_matches) == 1:
            already += 1
            continue
        raise RuntimeError(
            f"No se encontró texto original ni reemplazo único: {old[:100]}"
        )
    global_applied = 0
    for paragraph in paragraphs:
        original = paragraph_text(paragraph)
        updated = original
        for old, new in GLOBAL_REPLACEMENTS:
            global_applied += updated.count(old)
            updated = updated.replace(old, new)
        if updated != original:
            replace_paragraph_text(paragraph, updated)
    rendered = etree.tostring(
        root, encoding="UTF-8", xml_declaration=True, standalone=True
    )
    return rendered, {
        "substitutions_applied": applied,
        "global_replacements_applied": global_applied,
        "already_applied": already,
    }


def media_hashes(archive: zipfile.ZipFile) -> dict[str, str]:
    return {
        info.filename: hashlib.sha256(archive.read(info)).hexdigest()
        for info in archive.infolist()
        if info.filename.startswith("word/media/") and not info.is_dir()
    }


def patch_docx(source: Path, output: Path, force: bool) -> dict[str, object]:
    source = source.resolve()
    output = output.resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if output.exists() and not force:
        raise FileExistsError(f"La salida existe: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source) as archive:
        if archive.testzip() is not None:
            raise RuntimeError("DOCX de entrada dañado")
        infos = archive.infolist()
        names = [info.filename for info in infos]
        before = media_hashes(archive)
        document_xml, audit = patch_document(archive.read(DOCUMENT_MEMBER))
        comment = archive.comment
        payloads = {
            info.filename: (
                document_xml
                if info.filename == DOCUMENT_MEMBER
                else archive.read(info.filename)
            )
            for info in infos
        }
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output.parent, prefix=f".{output.name}.", suffix=".building"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        with zipfile.ZipFile(temporary, "w", allowZip64=True) as archive:
            archive.comment = comment
            for info in infos:
                archive.writestr(info, payloads[info.filename])
        with zipfile.ZipFile(temporary) as archive:
            bad = archive.testzip()
            if bad is not None:
                raise RuntimeError(f"DOCX de salida dañado: {bad}")
            if archive.namelist() != names:
                raise RuntimeError("Cambió la lista u orden de miembros ZIP")
            after = media_hashes(archive)
        if after != before:
            raise RuntimeError("Cambió al menos un SHA-256 de word/media")
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    return {
        "status": "complete",
        "source": str(source),
        "output": str(output),
        "media_verified_unchanged": len(before),
        **audit,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    result = patch_docx(args.input, args.output, args.force)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
