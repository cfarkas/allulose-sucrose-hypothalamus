#!/usr/bin/env python3
"""Collect legal web records and open-access full text for thesis references.

Metadata are requested from Crossref, OpenAlex and Europe PMC.  A PDF is
downloaded only when one of those scholarly indexes explicitly exposes an
open-access PDF URL.  Paywalls and access controls are never bypassed.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode, urlparse
from urllib.request import Request, urlopen


USER_AGENT = "Apotome-thesis-literature-audit/2026-08-30"


def text(value: object) -> str:
    return str(value or "").strip()


def year(item: dict[str, Any]) -> str:
    parts = item.get("issued", {}).get("date-parts", [])
    return str(parts[0][0]) if parts and parts[0] else "sf"


def first_author(item: dict[str, Any]) -> str:
    authors = item.get("author") or []
    if authors:
        return text(authors[0].get("family") or authors[0].get("literal")) or "autor"
    return "autor"


def slug(value: str, limit: int = 70) -> str:
    value = value.casefold()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")[:limit] or "fuente"


def request_json(url: str, timeout: float) -> tuple[dict[str, Any] | None, str]:
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    try:
        with urlopen(request, timeout=timeout) as response:
            payload = response.read()
        return json.loads(payload.decode("utf-8")), "ok"
    except HTTPError as exc:
        return None, f"HTTP {exc.code}"
    except (URLError, TimeoutError, json.JSONDecodeError, UnicodeDecodeError) as exc:
        return None, f"{type(exc).__name__}: {exc}"


def write_json(path: Path, payload: dict[str, Any] | list[Any]) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def candidate_pdf_urls(openalex: dict[str, Any] | None, europe_pmc: dict[str, Any] | None) -> list[str]:
    urls: list[str] = []
    if openalex:
        for field in ("best_oa_location", "primary_location"):
            location = openalex.get(field) or {}
            candidate = text(location.get("pdf_url"))
            if candidate:
                urls.append(candidate)
        for location in openalex.get("locations") or []:
            candidate = text((location or {}).get("pdf_url"))
            if candidate:
                urls.append(candidate)
    if europe_pmc:
        results = europe_pmc.get("resultList", {}).get("result", [])
        for record in results:
            pmcid = text(record.get("pmcid"))
            if pmcid and text(record.get("isOpenAccess")).upper() == "Y":
                urls.append(
                    f"https://www.ebi.ac.uk/europepmc/webservices/rest/{quote(pmcid)}/fullTextPDF"
                )
    seen: set[str] = set()
    return [url for url in urls if not (url in seen or seen.add(url))]


def download_pdf(url: str, destination: Path, timeout: float) -> tuple[bool, str]:
    if urlparse(url).scheme not in {"http", "https"}:
        return False, "unsupported URL scheme"
    request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/pdf"})
    temporary = destination.with_suffix(".downloading")
    try:
        with urlopen(request, timeout=timeout) as response, temporary.open("wb") as handle:
            first = response.read(8192)
            if not first.startswith(b"%PDF"):
                return False, "response is not a PDF"
            handle.write(first)
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                handle.write(block)
        temporary.replace(destination)
        return True, "ok"
    except HTTPError as exc:
        return False, f"HTTP {exc.code}"
    except (URLError, TimeoutError, OSError) as exc:
        return False, f"{type(exc).__name__}: {exc}"
    finally:
        if temporary.exists():
            temporary.unlink()


def collect(items_path: Path, output_dir: Path, timeout: float, pause: float) -> None:
    items = json.loads(items_path.read_text(encoding="utf-8"))
    if not isinstance(items, list):
        raise SystemExit(f"Expected a JSON list: {items_path}")
    sources_dir = output_dir / "sources"
    sources_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    total = len(items)
    for index, item in enumerate(items, 1):
        if not isinstance(item, dict):
            continue
        doi = text(item.get("DOI")).lower()
        title = text(item.get("title"))
        folder = sources_dir / (
            f"{index:03d}_{slug(first_author(item), 24)}_{year(item)}_{slug(title, 48)}"
        )
        folder.mkdir(parents=True, exist_ok=True)
        write_json(folder / "embedded_zotero_metadata.json", item)

        errors: list[str] = []
        crossref: dict[str, Any] | None = None
        openalex: dict[str, Any] | None = None
        europe_pmc: dict[str, Any] | None = None
        if doi:
            crossref_url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
            crossref_raw, status = request_json(crossref_url, timeout)
            if crossref_raw:
                crossref = crossref_raw.get("message", crossref_raw)
                write_json(folder / "crossref_record.json", crossref)
            else:
                errors.append(f"Crossref: {status}")

            openalex_url = f"https://api.openalex.org/works/https://doi.org/{quote(doi, safe='/')}"
            openalex, status = request_json(openalex_url, timeout)
            if openalex:
                write_json(folder / "openalex_record.json", openalex)
            else:
                errors.append(f"OpenAlex: {status}")

            query = urlencode({"query": f'DOI:"{doi}"', "format": "json", "pageSize": 5})
            europe_pmc, status = request_json(
                f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?{query}", timeout
            )
            if europe_pmc:
                write_json(folder / "europe_pmc_record.json", europe_pmc)
            else:
                errors.append(f"Europe PMC: {status}")
        elif title:
            query = urlencode({"query.bibliographic": title, "rows": 3})
            candidates, status = request_json(f"https://api.crossref.org/works?{query}", timeout)
            if candidates:
                write_json(folder / "crossref_title_candidates.json", candidates)
            else:
                errors.append(f"Crossref title search: {status}")

        pdf_url = ""
        pdf_status = "not available from scholarly OA indexes"
        for candidate in candidate_pdf_urls(openalex, europe_pmc):
            ok, status = download_pdf(candidate, folder / "full_text_open_access.pdf", timeout)
            if ok:
                pdf_url = candidate
                pdf_status = "downloaded and PDF signature verified"
                break
            errors.append(f"PDF {candidate}: {status}")

        provenance = {
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "title": title,
            "doi": doi,
            "doi_url": f"https://doi.org/{doi}" if doi else "",
            "embedded_url": text(item.get("URL")),
            "pdf_url": pdf_url,
            "pdf_status": pdf_status,
            "errors": errors,
            "policy": (
                "Metadata from scholarly APIs; full text only when an index supplied an "
                "explicit open-access PDF URL. No paywall or access control was bypassed."
            ),
        }
        write_json(folder / "retrieval_provenance.json", provenance)
        rows.append(
            {
                "index": index,
                "author": first_author(item),
                "year": year(item),
                "title": title,
                "doi": doi,
                "crossref": bool(crossref),
                "openalex": bool(openalex),
                "europe_pmc": bool(europe_pmc),
                "open_access_pdf": bool(pdf_url),
                "pdf_url": pdf_url,
                "errors": " | ".join(errors),
                "folder": str(folder.relative_to(output_dir)),
            }
        )
        print(
            f"[{index:02d}/{total:02d}] {first_author(item)} {year(item)} | "
            f"DOI={'yes' if doi else 'no'} | OA PDF={'yes' if pdf_url else 'no'}",
            flush=True,
        )
        if pause:
            time.sleep(pause)

    fieldnames = list(rows[0]) if rows else []
    with (output_dir / "COLLECTION_REPORT.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    downloaded = sum(bool(row["open_access_pdf"]) for row in rows)
    readme = [
        "# Literatura de la tesis recopilada el 30-08-2026",
        "",
        f"Fuentes Zotero procesadas: **{len(rows)}**.",
        "",
        f"Textos completos de acceso abierto descargados: **{downloaded}**.",
        "",
        "Cada carpeta conserva los metadatos embebidos en Word, los registros web disponibles,",
        "la procedencia de la descarga y, cuando el índice académico publicó una URL PDF de",
        "acceso abierto, una copia local cuya firma `%PDF` fue verificada.",
        "",
        "La ausencia de PDF no implica ausencia de la fuente: el registro bibliográfico y el enlace",
        "DOI se conservan. No se eludieron muros de pago ni controles de acceso.",
        "",
        "Archivos principales:",
        "",
        "- `CITATION_MANIFEST.md`: inventario legible.",
        "- `zotero_items.json` y `.csv`: metadatos extraídos del DOCX.",
        "- `COLLECTION_REPORT.csv`: estado de verificación y descarga por fuente.",
        "- `sources/`: evidencia por fuente.",
    ]
    (output_dir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "sources": len(rows), "oa_pdfs": downloaded}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--items",
        type=Path,
        default=Path("literature_webscrap_30_08_2026/zotero_items.json"),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("literature_webscrap_30_08_2026")
    )
    parser.add_argument("--timeout", type=float, default=45.0)
    parser.add_argument("--pause", type=float, default=0.15)
    args = parser.parse_args()
    collect(args.items, args.output_dir, args.timeout, args.pause)


if __name__ == "__main__":
    main()
