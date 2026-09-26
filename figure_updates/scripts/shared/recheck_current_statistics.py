#!/usr/bin/env python3
"""Recompute the current Fig. 3–4 statistics from portable cell/geometry inputs.

This starts after reviewed segmentation and cell classification; it does not
retrain models or resegment microscopy. Run from Paper or figure_updates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import tempfile

import pandas as pd
from pandas.testing import assert_frame_equal

import fig3_scale_resolved as spatial
import fig3_robust_abundance as abundance


def check(paper: Path, output: Path) -> dict:
    receipt = {}
    for cohort, figure in (("NPY", "Fig3"), ("POMC", "Fig4")):
        source = paper / figure / "source_data" / "current_statistics"
        target = output / cohort
        # An existing destination is rejected rather than mixed with old results.
        shutil.copytree(source, target)
        spatial.run(target, cohort=cohort)
        names = ["scale_resolved_curves.csv", "scale_resolved_strata.csv.gz",
                 "scale_resolved_test.csv", "scale_resolved_envelope.csv"]
        if cohort == "NPY":
            abundance.run(target)
            names.append("abundance_robust_statistics.csv")
        for name in names:
            expected = pd.read_csv(source / name)
            actual = pd.read_csv(target / name)
            assert_frame_equal(actual, expected, check_dtype=False,
                               check_exact=False, rtol=1e-10, atol=1e-10,
                               obj=f"{cohort}/{name}")
        results = pd.read_csv(target / "scale_resolved_test.csv")
        receipt[cohort] = {"status": "PASS", "tables_recomputed": names,
                          "results": results[["endpoint", "unit", "n_water",
                                              "n_sucrose", "n_allulose", "p_value",
                                              "p_holm_two_endpoint_family"]].to_dict("records")}
        print(f"PASS: {cohort}, {len(names)} complete numerical tables", flush=True)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paper-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path,
                        help="Optional new directory in which to retain recomputed files")
    args = parser.parse_args()
    if args.output_dir is None:
        with tempfile.TemporaryDirectory(prefix="paper-statistics-") as directory:
            receipt = check(args.paper_root.resolve(), Path(directory))
    else:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        receipt = check(args.paper_root.resolve(), args.output_dir.resolve())
        (args.output_dir / "recheck_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
