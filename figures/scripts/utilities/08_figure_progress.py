#!/usr/bin/env python3
"""Show validated figure progress for ``reproduce_all_figures.sh``.

The launcher writes tab-separated ``START``, ``DONE`` and ``FAIL`` events on
standard input. One tqdm unit represents one complete figure package and is
credited only after that package passes its output validator.
"""

from __future__ import annotations

import argparse
import sys

from tqdm import tqdm


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--total", type=int, required=True)
    args = parser.parse_args()
    if args.total < 1:
        parser.error("--total must be positive")

    completed = 0
    current = ""
    failed = False
    bar = tqdm(
        total=args.total,
        desc="Preparing figures",
        unit="figure",
        dynamic_ncols=True,
        file=sys.stderr,
    )
    try:
        for raw_line in sys.stdin:
            line = raw_line.rstrip("\n")
            try:
                event, label = line.split("\t", 1)
            except ValueError:
                print(f"Invalid figure-progress event: {line!r}", file=sys.stderr)
                return 2

            if event == "START":
                current = label
                bar.set_description_str(f"Rebuilding {label}", refresh=True)
                bar.set_postfix_str("running", refresh=True)
            elif event == "DONE":
                if completed >= args.total:
                    print("Too many completed figure-progress events", file=sys.stderr)
                    return 2
                completed += 1
                bar.update(1)
                bar.set_description_str(f"Validated {label}", refresh=True)
                bar.set_postfix_str("complete", refresh=True)
                current = ""
            elif event == "FAIL":
                current = label
                failed = True
                bar.set_description_str(f"Failed: {label}", refresh=True)
                bar.set_postfix_str("stopped", refresh=True)
            else:
                print(f"Unknown figure-progress event: {event!r}", file=sys.stderr)
                return 2
    finally:
        if completed == args.total:
            bar.set_description_str("All figures validated", refresh=True)
            bar.set_postfix_str("complete", refresh=True)
        elif current:
            bar.set_description_str(f"Stopped: {current}", refresh=True)
        bar.close()

    if completed != args.total and not failed:
        print(
            f"Figure progress ended at {completed}/{args.total} without a failure event",
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
