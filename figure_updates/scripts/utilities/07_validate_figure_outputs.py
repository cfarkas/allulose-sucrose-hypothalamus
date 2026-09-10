#!/usr/bin/env python3
"""Validate complete bilingual publication output for every live paper figure."""

from __future__ import annotations

import argparse
import filecmp
import json
import hashlib
import re
from dataclasses import dataclass, replace
from pathlib import Path

from PIL import Image


@dataclass(frozen=True)
class Contract:
    masters: tuple[str, ...]
    panel_count: int
    legend_english: str
    legend_spanish: str


CONTRACTS = {
    "FigS7": Contract(("Figure_S7",), 3, "Figure_S7_LEGEND.txt", "Figure_S7_LEGEND_spanish.txt"),
    "FigS8": Contract(("Figure_S8",), 2, "Figure_S8_LEGEND.txt", "Figure_S8_LEGEND_spanish.txt"),
    "FigS6": Contract(
        ("Figure_S6",), 4,
        "Figure_S6_LEGEND.txt", "Figure_S6_LEGEND_spanish.txt",
    ),
    "Fig1": Contract(
        ("Figure_1_behavior_experiment_1",), 4,
        "Figure_1_behavior_experiment_1_LEGEND.txt",
        "Figure_1_behavior_experiment_1_LEGEND_spanish.txt",
    ),
    "Fig2": Contract(
        ("Figure_2",), 4,
        "Figure_2_legend.txt", "Figure_2_legend_spanish.txt",
    ),
    "Fig3": Contract(
        ("Figure_3_cFos_NPY",), 7,
        "Figure_3_cFos_NPY_LEGEND.txt",
        "Figure_3_cFos_NPY_LEGEND_spanish.txt",
    ),
    "Fig4": Contract(
        ("Figure_4", "Figure_ARC_ME_multipanel"), 10,
        "Figure_4_LEGEND.txt", "Figure_4_LEGEND_spanish.txt",
    ),
    "Fig5": Contract(
        ("Figure_5", "Figure_GFAP_Iba1_microglia_ARC_ME"), 9,
        "Figure_GFAP_Iba1_microglia_ARC_ME_LEGEND.txt",
        "Figure_GFAP_Iba1_microglia_ARC_ME_LEGEND_spanish.txt",
    ),
    "FigS1": Contract(
        ("Figure_S1_single_bottle_male",), 4,
        "Figure_S1_single_bottle_male_LEGEND.txt",
        "Figure_S1_single_bottle_male_LEGEND_spanish.txt",
    ),
    "FigS2": Contract(
        ("Figure_S2_single_bottle_female",), 4,
        "Figure_S2_single_bottle_female_LEGEND.txt",
        "Figure_S2_single_bottle_female_LEGEND_spanish.txt",
    ),
    "FigS3": Contract(
        ("Figure_S3",), 6,
        "Figure_S3_LEGEND.txt", "Figure_S3_LEGEND_spanish.txt",
    ),
    "FigS4": Contract(
        ("Figure_S4_behavior_preference",), 5,
        "Figure_S4_behavior_preference_LEGEND.txt",
        "Figure_S4_behavior_preference_LEGEND_spanish.txt",
    ),
    "FigS5": Contract(
        ("Figure_S5",), 10,
        "Figure_S5_LEGEND.txt", "Figure_S5_LEGEND_spanish.txt",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def output_layout(root: Path) -> tuple[Path, Path, Path]:
    """Return master, panels and legends locations for normal/Figure-2 layouts."""
    if (root / "outputs").is_dir():
        return root / "outputs", root / "outputs" / "panels", root / "legends"
    return root, root / "panels", root / "legends"


def master_stem(directory: Path, contract: Contract) -> str:
    for stem in contract.masters:
        if (directory / f"{stem}.png").is_file() and (
            directory / f"{stem}.pdf"
        ).is_file():
            return stem
    raise RuntimeError(
        f"No complete PNG/PDF master in {directory}; expected one of "
        + ", ".join(contract.masters)
    )


def image_dpi(path: Path) -> tuple[float, float]:
    with Image.open(path) as image:
        dpi = image.info.get("dpi", (0.0, 0.0))
        require(image.width > 0 and image.height > 0, f"Empty raster: {path}")
    if isinstance(dpi, (int, float)):
        return float(dpi), float(dpi)
    return float(dpi[0]), float(dpi[1])


def panel_letter(path: Path) -> str:
    match = re.search(r"(?:Panel|panel)_([A-Z])(?:_|$)", path.stem)
    require(match is not None, f"Cannot identify panel letter: {path.name}")
    return str(match.group(1))


def validate(figure: str, root: Path) -> dict[str, object]:
    contract = CONTRACTS[figure]
    complete=root/'provenance/complete_conditions_20260908.json'
    if figure=='Fig3' and complete.exists():
        receipt=json.loads(complete.read_text())
        require(receipt['schema']=='figure3_complete_three_conditions_v1','Unknown complete Figure 3 provenance')
        require(set(receipt['panel_map'])==set('ABCDEFGH'),'Incomplete three-condition panel map')
        require(receipt.get('publication_language')=='en','Complete Figure 3 must have an English master')
        master=root/'Figure_3_cFos_NPY.pdf'
        require(hashlib.sha256(master.read_bytes()).hexdigest()==next(x['sha256'] for x in receipt['figures'] if x['language']=='en'),'Master does not match complete Figure 3 receipt')
        require(not list(root.glob('Figure_3_cFos_NPY_spanish.*')),'Unexpected complete Spanish master; retain only bilingual subpanels and legends')
        contract=replace(contract,panel_count=8)
    master_dir, panels_dir, legends_dir = output_layout(root)
    stem = master_stem(master_dir, contract)
    master_png = master_dir / f"{stem}.png"
    master_pdf = master_dir / f"{stem}.pdf"
    require(master_pdf.read_bytes()[:4] == b"%PDF", f"Invalid PDF: {master_pdf}")
    require(min(image_dpi(master_png)) >= 590.0, f"Master is below 600 dpi: {master_png}")

    require(panels_dir.is_dir(), f"Panel directory missing: {panels_dir}")
    per_format: dict[str, dict[str, int]] = {}
    expected_letters: set[str] | None = None
    for suffix in (".png", ".pdf"):
        all_files = sorted(panels_dir.glob(f"*{suffix}"))
        english = [path for path in all_files if not path.stem.endswith("_spanish")]
        spanish = [path for path in all_files if path.stem.endswith("_spanish")]
        require(
            len(english) == contract.panel_count,
            f"{figure}: expected {contract.panel_count} English {suffix} panels, got {len(english)}",
        )
        require(
            len(spanish) == contract.panel_count,
            f"{figure}: expected {contract.panel_count} Spanish {suffix} panels, got {len(spanish)}",
        )
        english_letters = {panel_letter(path) for path in english}
        spanish_letters = {panel_letter(path) for path in spanish}
        require(
            len(english_letters) == contract.panel_count,
            f"{figure}: duplicate/missing English panel letters: {sorted(english_letters)}",
        )
        require(
            spanish_letters == english_letters,
            f"{figure}: English/Spanish panel letters differ",
        )
        if expected_letters is None:
            expected_letters = english_letters
        else:
            require(english_letters == expected_letters, f"{figure}: PNG/PDF letters differ")
        for path in english:
            partner = path.with_name(f"{path.stem}_spanish{path.suffix}")
            require(partner.is_file(), f"Missing Spanish pair for {path.name}")
            if suffix == ".png":
                require(
                    not filecmp.cmp(path, partner, shallow=False),
                    f"{figure}: English/Spanish PNG panels are byte-identical: {path.name}",
                )
        for path in all_files:
            require(path.stat().st_size > 100, f"Empty panel: {path}")
            if suffix == ".pdf":
                require(path.read_bytes()[:4] == b"%PDF", f"Invalid panel PDF: {path}")
            else:
                require(min(image_dpi(path)) >= 590.0, f"Panel is below 600 dpi: {path}")
        per_format[suffix[1:]] = {"english": len(english), "spanish": len(spanish)}

    english_legend = legends_dir / contract.legend_english
    spanish_legend = legends_dir / contract.legend_spanish
    legend_bodies: dict[str, str] = {}
    for path, language in ((english_legend, "English"), (spanish_legend, "Spanish")):
        require(path.is_file(), f"{figure}: {language} figure legend missing: {path}")
        body = path.read_text(encoding="utf-8")
        require(len(body.strip()) >= 80, f"{figure}: {language} figure legend is too short")
        legend_bodies[language] = body
    require(
        legend_bodies["English"] != legend_bodies["Spanish"],
        f"{figure}: English and Spanish figure legends are identical",
    )
    require(
        re.search(r"\bFigure\b", legend_bodies["English"], re.IGNORECASE) is not None,
        f"{figure}: English legend lacks an English Figure heading",
    )
    require(
        re.search(r"\bFigura\b", legend_bodies["Spanish"], re.IGNORECASE) is not None,
        f"{figure}: Spanish legend lacks a Spanish Figura heading",
    )

    return {
        "status": "PASS",
        "figure": figure,
        "root": str(root),
        "master": stem,
        "master_language": "English",
        "master_dpi": image_dpi(master_png),
        "panel_letters": sorted(expected_letters or ()),
        "panels": per_format,
        "figure_legends": {
            "english": str(english_legend),
            "spanish": str(spanish_legend),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure", choices=sorted(CONTRACTS))
    parser.add_argument("--root", type=Path)
    parser.add_argument(
        "--all-live", action="store_true",
        help="Validate every promoted FigN package under the Paper root.",
    )
    args = parser.parse_args()
    paper = Path(__file__).resolve().parents[2]
    if args.all_live:
        require(args.figure is None and args.root is None, "--all-live takes no figure/root")
        results = [validate(figure, paper / figure) for figure in CONTRACTS]
        print(json.dumps({"status": "PASS", "figures": results}, indent=2))
        return 0
    require(args.figure is not None and args.root is not None, "Use --all-live or --figure with --root")
    result = validate(str(args.figure), args.root.expanduser().resolve())
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(json.dumps({"status": "FAIL", "error": str(error)}, indent=2))
        raise SystemExit(1)
