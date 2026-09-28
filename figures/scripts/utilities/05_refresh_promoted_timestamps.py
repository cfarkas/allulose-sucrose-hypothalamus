#!/usr/bin/env python3
"""Rewrite each promoted FigN/ file from a fresh render, without changing bytes.

04_promote_figure.py now atomically rewrites all files when the canonical
launcher passes --replace-current, including files whose bytes already match.
This compatibility helper remains for older staged run-books and targeted
manual refreshes. It makes a figure directory reflect the current successful
run while retaining exact bytes and uses the same canonical mapping and atomic
installer as the main promotion path.

This rewrites those files from the render that produced them, so FigN/ holds the
current rebuild's output and carries its time. It reuses 04_promote_figure.py's
own planned_copies(), so the render-to-promoted name mapping cannot drift from
the promoter's.

It is deliberately not a promotion. It rewrites only files whose bytes are
already identical, and preserves each destination's permissions. Anything whose
bytes differ, and anything not already in FigN/, is reported and left alone: a
real content change has to go through 04_promote_figure.py, so that the
superseded copy is archived to FigN/legacy/superseded_promotion_<stamp>/ first
and a receipt is written.

Figure 3 is handled conservatively. Its promoted bytes are asserted from inside
its own frozen reproduction bundle, and its PDFs are not byte-stable across runs
while its PNGs are pixel-bound, so its byte-identical files are refreshed and
everything else is left exactly as it is.

Run it after the README.txt or README2.txt run-book, against the same $OUT:

    python scripts/utilities/05_refresh_promoted_timestamps.py --out "$OUT"
    python scripts/utilities/05_refresh_promoted_timestamps.py --out "$OUT" --apply

Without --apply it reports and writes nothing.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve()
PROMOTER = HERE.parent / "04_promote_figure.py"

_spec = importlib.util.spec_from_file_location("promote_figure", PROMOTER)
_promote = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_promote)

PAPER = _promote.PAPER
# Where the run-book leaves each figure's finished render, relative to $OUT.
RENDER_SUBDIRS = {
    "Fig1": "Fig1/figure",
    "Fig2": "Fig2/render/outputs",
    "Fig3": "Fig3/out/final",
    "Fig4": "Fig4/figure",
    "Fig5": "Fig5/figure",
    "FigS5": "FigS5/figure",
    "FigS1": "FigS1/figure",
    "FigS2": "FigS2/figure",
    "FigS3": "FigS3/figure",
    "FigS4": "FigS4/figure",
    "FigS6": "FigS6/figure",
    "FigS7": "FigS7/figure",
    "FigS8": "FigS8/figure",
}
# Promoted bytes asserted by a frozen reproduction bundle: refresh only what is
# byte-identical, never skip the whole figure because the rest differs.
CONSERVATIVE = {"Fig3"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, default=Path("/tmp/apotome_rebuild"),
                        help="The run-book's $OUT, holding this rebuild's renders.")
    parser.add_argument("--figure", choices=sorted(RENDER_SUBDIRS), action="append",
                        help="Limit to one figure; repeatable. Default is all thirteen.")
    parser.add_argument("--apply", action="store_true",
                        help="Rewrite the files. Without it nothing is written.")
    args = parser.parse_args()

    out = args.out.expanduser().resolve()
    figures = args.figure or sorted(RENDER_SUBDIRS)
    identical = differ = absent = written = 0

    for figure in figures:
        source = out / RENDER_SUBDIRS[figure]
        target = PAPER / figure
        if not source.is_dir():
            print(f"{figure}: no render at {source}\n")
            continue
        if not target.is_dir():
            print(f"{figure}: no figure directory at {target}\n")
            continue

        same, different, missing = [], [], []
        for src, dst in _promote.planned_copies(source, target, figure):
            if not dst.exists():
                missing.append(dst)
            elif _promote.sha256(src) == _promote.sha256(dst):
                same.append((src, dst))
            else:
                different.append(dst)

        tag = " (conservative)" if figure in CONSERVATIVE else ""
        print(f"{figure}{tag}: {len(same)} byte-identical, {len(different)} differ, "
              f"{len(missing)} not in {figure}/")
        for path in different:
            print(f"    left alone, bytes differ  {path.relative_to(PAPER)}")
        for path in missing:
            print(f"    left alone, not promoted  {path.relative_to(PAPER)}")
        identical += len(same); differ += len(different); absent += len(missing)

        if args.apply:
            if (different or missing) and figure not in CONSERVATIVE:
                print(f"    -> skipped {figure}: content changed, "
                      f"promote it with 04_promote_figure.py instead")
                print()
                continue
            for src, dst in same:
                _promote.atomic_install(src, dst, stamp_current=True)
                if _promote.sha256(src) != _promote.sha256(dst):
                    raise RuntimeError(f"Byte mismatch after refreshing {dst}")
            written += len(same)
            print(f"    -> rewrote {len(same)} files from {source}")
        print()

    print(f"{identical} byte-identical, {differ} differ, {absent} not promoted, "
          f"{written} rewritten")
    if not args.apply:
        print("Nothing was written. Pass --apply to rewrite.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
