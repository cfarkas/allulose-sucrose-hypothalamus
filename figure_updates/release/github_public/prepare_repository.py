#!/usr/bin/env python3
"""Assemble an absent, publication-ready GitHub repository without publishing it."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

from repository_tools import PublicationError, assemble_repository


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--transport-manifest", type=Path, required=True)
    parser.add_argument("--url-map", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Absent destination for the final small GitHub repository.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output = assemble_repository(
        Path(__file__).resolve().parent,
        args.manifest,
        args.transport_manifest,
        args.url_map,
        args.metadata,
        args.output,
    )
    print(f"[PASS] publication-ready repository assembled: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PublicationError as exc:
        print(f"[BLOCKED] {exc}", file=sys.stderr)
        raise SystemExit(2)
