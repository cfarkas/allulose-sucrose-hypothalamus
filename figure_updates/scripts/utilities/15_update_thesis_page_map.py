#!/usr/bin/env python3
"""Derive the thesis' logical page map from its rendered PDF."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path


def normalized_lines(page: str) -> list[str]:
    return [re.sub(r"\s+", " ", line).strip()
            for line in page.splitlines() if line.strip()]


def find_page(pages: list[list[str]], pattern: str, *, start: int = 1) -> int:
    expression = re.compile(pattern, re.IGNORECASE)
    for physical_page in range(max(1, start), len(pages) + 1):
        if any(expression.fullmatch(line) for line in pages[physical_page - 1]):
            return physical_page
    raise RuntimeError(f"No se encontró en el PDF: {pattern}")


def derive(pdf: Path) -> dict[str, int]:
    result = subprocess.run(
        ["pdftotext", "-layout", str(pdf), "-"], check=True,
        stdout=subprocess.PIPE, text=True)
    raw_pages = result.stdout.split("\f")
    if raw_pages and not raw_pages[-1].strip():
        raw_pages.pop()
    pages = [normalized_lines(page) for page in raw_pages]

    front_patterns = {
        "acknowledgements": r"AGRADECIMIENTOS",
        "summary": r"RESUMEN",
        "abstract": r"ABSTRACT",
        "figure_index": r"ÍNDICE DE FIGURAS",
        "table_index": r"ÍNDICE DE TABLAS",
        "abbreviations": r"ABREVIATURAS",
    }
    page_map = {key: find_page(pages, pattern)
                for key, pattern in front_patterns.items()}

    introduction_physical = find_page(
        pages, r"(?:I\.\s*)?INTRODUCCIÓN",
        start=page_map["abbreviations"])
    main_offset = introduction_physical - 1

    main_patterns = {
        "introduction": r"(?:I\.\s*)?INTRODUCCIÓN",
        "hypothesis": r"(?:II\.\s*)?HIPÓTESIS Y OBJETIVOS",
        "methodology": r"(?:III\.\s*)?METODOLOGÍA",
        "results": r"(?:IV\.\s*)?RESULTADOS",
        "discussion": r"(?:V\.\s*)?DISCUSIÓN",
        "conclusions": r"(?:VI\.\s*)?CONCLUSIONES",
        "references": r"(?:VII\.\s*)?REFERENCIAS",
        "paper_appendix": r"ANEXO\. MANUSCRITO EN PREPARACIÓN PARA ENVÍO",
    }
    for key, pattern in main_patterns.items():
        physical = find_page(pages, pattern, start=introduction_physical)
        page_map[key] = physical - main_offset

    for number in range(1, 72):
        physical = find_page(
            pages, rf"Figura {number}\..*", start=introduction_physical)
        page_map[f"fig_{number}"] = physical - main_offset

    methodology_physical = page_map["methodology"] + main_offset
    for number in (1, 2):
        physical = find_page(
            pages, rf"Tabla {number}\..*", start=methodology_physical)
        page_map[f"table_{number}"] = physical - main_offset

    appendix_physical = page_map["paper_appendix"] + main_offset
    for label in ("1", "2", "3", "4", "5", "S1", "S2", "S3", "S4", "S5"):
        physical = find_page(
            pages, rf"Artículo\s*[—-]\s*Figure {label}",
            start=appendix_physical)
        page_map[f"paper_fig_{label.casefold()}"] = physical - main_offset

    expected_keys = 6 + 8 + 71 + 2 + 10
    if len(page_map) != expected_keys:
        raise RuntimeError(
            f"Mapa incompleto: {len(page_map)} claves; se esperaban {expected_keys}")
    page_map["_physical_pdf_pages"] = len(pages)
    page_map["_main_numbering_physical_offset"] = main_offset
    return page_map


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    page_map = derive(args.pdf.resolve())
    rendered = json.dumps(page_map, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
