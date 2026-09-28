#!/usr/bin/env python3
"""Write a fresh Figure 2 animal-level analysis stage from raw_data only."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve()
# The paper root is found by walking up to the directory that holds README.txt
# and scripts/setup, so this works from Fig2/ and from its scripts/Fig2 mirror
# alike. Deriving it as HERE.parent.parent silently resolved to Paper/scripts
# when run from the mirror, and every raw path under it then pointed at
# scripts/Fig2/raw_data, which does not exist.
PAPER = next((path for path in HERE.parents
              if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                  and (path / "README.txt").is_file()))), None)
if PAPER is None:  # pragma: no cover
    raise RuntimeError(f"Could not locate the Paper directory above {HERE}")
FIGURE = PAPER / "Fig2"
RAW = FIGURE / "raw_data" / "Apotome"
ANALYZER = FIGURE / "01_analyze_global_cfos.py"
INPUT = RAW / "combined_human_region_plots" / "combined_raw_per_animal_human_region_summary.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_analyzer():
    spec = importlib.util.spec_from_file_location("fig2_staged_analyzer", ANALYZER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {ANALYZER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", type=Path, help="Fresh absolute stage directory")
    args = parser.parse_args()
    stage = args.stage.expanduser()
    if not stage.is_absolute():
        raise ValueError("Stage path must be absolute")
    if stage.exists():
        raise FileExistsError(f"Fresh-only stage already exists: {stage}")
    analysis = stage / "analysis"
    analysis.mkdir(parents=True, exist_ok=False)

    module = load_analyzer()
    module.CONDITION_EVIDENCE = {
        "behavior_weights": RAW / "Comportamiento/raw/experiment_1_single_bottle/weights.csv",
        "binning_metadata": RAW / "24_06_2025/analysis/bins_heuristic_stats_all.csv",
        "assay_samplesheet": RAW / "24_06_2025/analysis/apotome_samplesheet_template.csv",
        "image_discovery": RAW / "24_06_2025/analysis/human_regions_v9/discovered_images.csv",
        "annotation_status": RAW / "24_06_2025/analysis/human_regions_v9/annotation_status.csv",
    }
    previous_argv = sys.argv
    try:
        sys.argv = [
            str(ANALYZER), "--input", str(INPUT), "--outdir", str(analysis),
            "--strict-paper-validation",
        ]
        result = int(module.main())
    finally:
        sys.argv = previous_argv
    if result != 0:
        return result

    outputs = []
    for path in sorted(analysis.iterdir()):
        if path.is_file():
            outputs.append(
                {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path)}
            )
    receipt = {
        "status": "passed",
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "fresh_stage": str(stage),
        "input": str(INPUT),
        "input_sha256": sha256(INPUT),
        "analyzer": str(ANALYZER),
        "analyzer_sha256": sha256(ANALYZER),
        "outputs": outputs,
        "next": "run stage_native_from_raw_bundle.py with this same stage",
    }
    (stage / "STAGE_INPUT_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"[OK] Fresh Figure 2 analysis stage: {stage}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
