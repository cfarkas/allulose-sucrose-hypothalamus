#!/usr/bin/env python3
"""Verify the ten current Figure 3 PDFs, their hashes and PDF text bounds."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import fitz


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def expected_pdfs():
    paths = {f"quantitative_inputs/Figure_3_cFos_NPY_Panel_{letter}{language}.pdf"
             for letter in ("D", "E") for language in ("", "_spanish")}
    paths.update(f"spatial_panels/Fig3/Spatial_shape_definition_NPY{language}.pdf"
                 for language in ("", "_spanish"))
    paths.update(f"spatial_panels/Fig3/{language}/Figure3_Spatial_{endpoint}_occurrence.pdf"
                 for language in ("en", "es") for endpoint in ("cFOS", "cFOS_NPY"))
    return paths


def validate(directory):
    directory = Path(directory).resolve()
    receipt_path = directory / "panel_render_receipt.json"
    receipt = json.loads(receipt_path.read_text())
    current_pdfs = {name for name in receipt["files"] if Path(name).suffix == ".pdf"}
    expected = expected_pdfs()
    if current_pdfs != expected:
        raise ValueError(f"Expected precisely ten current panel PDFs; missing={sorted(expected-current_pdfs)}, unexpected={sorted(current_pdfs-expected)}")
    panels = []
    for name in sorted(expected):
        path = directory / name
        actual_hash = digest(path)
        if actual_hash != receipt["files"][name]:
            raise ValueError(f"PDF changed after rendering: {name}")
        with fitz.open(path) as document:
            if len(document) != 1:
                raise ValueError(f"Panel must contain one page: {name}")
            page = document[0]
            bounds = page.rect
            tolerance = 0.5  # PDF font metrics can differ by a fraction of a point.
            outside = []
            text_spans = 0
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line["spans"]:
                        if not span["text"].strip():
                            continue
                        text_spans += 1
                        x0, y0, x1, y1 = span["bbox"]
                        if x0 < bounds.x0-tolerance or y0 < bounds.y0-tolerance or x1 > bounds.x1+tolerance or y1 > bounds.y1+tolerance:
                            outside.append({"text": span["text"], "bbox": list(span["bbox"])})
            if not text_spans:
                raise ValueError(f"Panel has no extractable text to verify: {name}")
            panels.append({"path": name, "sha256": actual_hash,
                           "page_size_pt": [bounds.width, bounds.height],
                           "text_spans_checked": text_spans, "text_outside_canvas": outside})
    passed = all(not panel["text_outside_canvas"] for panel in panels)
    result = {"status": "PASS" if passed else "FAIL", "pdf_count": len(panels),
              "all_text_inside_canvas": passed, "text_bound_tolerance_pt": 0.5,
              "render_receipt_sha256": digest(receipt_path), "validator_sha256": digest(__file__),
              "panels": panels}
    (directory / "panel_render_validation.json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    if not passed:
        raise ValueError("Text extends outside one or more current Figure 3 panels; see panel_render_validation.json")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", required=True, type=Path)
    args = parser.parse_args()
    result = validate(args.analysis_dir)
    print(f"PASS: {result['pdf_count']} current Figure 3 PDFs; all text inside canvas.")


if __name__ == "__main__":
    main()
