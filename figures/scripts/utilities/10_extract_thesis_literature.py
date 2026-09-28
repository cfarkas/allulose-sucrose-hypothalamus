#!/usr/bin/env python3
"""Extract the complete Zotero CSL corpus embedded in the thesis DOCX.

The script is deliberately read-only with respect to the Word document.  It
writes a deterministic literature manifest that can be used by the thesis
builder and by the web-literature collector.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from zipfile import ZipFile

from lxml import etree


WORD_NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
ZOTERO_PREFIX = "ADDIN ZOTERO_ITEM CSL_CITATION "


def citation_year(item: dict[str, object]) -> str:
    parts = item.get("issued", {})
    if isinstance(parts, dict):
        date_parts = parts.get("date-parts", [])
        if isinstance(date_parts, list) and date_parts and date_parts[0]:
            return str(date_parts[0][0])
    return "s.f."


def author_label(item: dict[str, object]) -> str:
    authors = item.get("author", [])
    if not isinstance(authors, list) or not authors:
        return str(item.get("publisher", "Autor no indicado"))
    labels: list[str] = []
    for author in authors:
        if not isinstance(author, dict):
            continue
        labels.append(str(author.get("family") or author.get("literal") or ""))
    labels = [label for label in labels if label]
    if not labels:
        return "Autor no indicado"
    if len(labels) == 1:
        return labels[0]
    if len(labels) == 2:
        return f"{labels[0]} y {labels[1]}"
    return f"{labels[0]} et al."


def stable_key(item: dict[str, object]) -> str:
    doi = str(item.get("DOI", "")).strip().lower()
    if doi:
        return f"doi:{doi}"
    url = str(item.get("URL", "")).strip().lower()
    if url:
        return f"url:{url}"
    title = re.sub(r"\s+", " ", str(item.get("title", "")).strip().lower())
    return f"title:{title}"


def extract(docx: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    with ZipFile(docx) as archive:
        root = etree.fromstring(archive.read("word/document.xml"))

    citations: list[dict[str, object]] = []
    unique: dict[str, dict[str, object]] = {}
    for instruction in root.xpath("//w:instrText/text()", namespaces=WORD_NS):
        if ZOTERO_PREFIX not in instruction:
            continue
        payload = instruction.split(ZOTERO_PREFIX, 1)[1].strip()
        citation = json.loads(payload)
        citations.append(citation)
        for cited in citation.get("citationItems", []):
            item = cited.get("itemData", {})
            if not isinstance(item, dict):
                continue
            key = stable_key(item)
            previous = unique.get(key)
            if previous is None or len(json.dumps(item, ensure_ascii=False)) > len(
                json.dumps(previous, ensure_ascii=False)
            ):
                unique[key] = item

    items = sorted(
        unique.values(),
        key=lambda item: (
            author_label(item).casefold(),
            citation_year(item),
            str(item.get("title", "")).casefold(),
        ),
    )
    return citations, items


def write_outputs(
    output_dir: Path,
    source: Path,
    citations: list[dict[str, object]],
    items: list[dict[str, object]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "zotero_citations.json").write_text(
        json.dumps(citations, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "zotero_items.json").write_text(
        json.dumps(items, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    columns = [
        "author",
        "year",
        "title",
        "container_title",
        "type",
        "doi",
        "url",
        "abstract",
    ]
    with (output_dir / "zotero_items.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for item in items:
            writer.writerow(
                {
                    "author": author_label(item),
                    "year": citation_year(item),
                    "title": item.get("title", ""),
                    "container_title": item.get("container-title", ""),
                    "type": item.get("type", ""),
                    "doi": item.get("DOI", ""),
                    "url": item.get("URL", ""),
                    "abstract": item.get("abstract", ""),
                }
            )

    item_types = Counter(str(item.get("type", "unknown")) for item in items)
    lines = [
        "# Corpus bibliográfico extraído de la tesis",
        "",
        f"Documento fuente: `{source.name}`",
        "",
        f"Campos de cita Zotero: **{len(citations)}**",
        "",
        f"Fuentes únicas: **{len(items)}**",
        "",
        "Tipos: " + ", ".join(f"{key}={value}" for key, value in sorted(item_types.items())),
        "",
        "## Fuentes",
        "",
    ]
    for index, item in enumerate(items, 1):
        doi = str(item.get("DOI", "")).strip()
        url = str(item.get("URL", "")).strip()
        locator = f"https://doi.org/{doi}" if doi else url
        suffix = f" — {locator}" if locator else ""
        lines.append(
            f"{index}. {author_label(item)} ({citation_year(item)}). "
            f"{str(item.get('title', '')).strip()}.{suffix}"
        )
    (output_dir / "CITATION_MANIFEST.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--docx", type=Path, default=Path("TESIS_FINAL_NS(4).docx"))
    parser.add_argument(
        "--output-dir", type=Path, default=Path("literature_webscrap_30_08_2026")
    )
    args = parser.parse_args()
    if not args.docx.is_file():
        raise SystemExit(f"Thesis DOCX not found: {args.docx}")
    citations, items = extract(args.docx)
    if not citations or not items:
        raise SystemExit("No embedded Zotero citations were found")
    write_outputs(args.output_dir, args.docx, citations, items)
    print(
        json.dumps(
            {
                "status": "PASS",
                "docx": str(args.docx),
                "zotero_fields": len(citations),
                "unique_sources": len(items),
                "output_dir": str(args.output_dir),
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
