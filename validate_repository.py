#!/usr/bin/env python3
"""Fail-closed validation for an assembled public GitHub repository."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from repository_tools import PublicationError, validate_repository


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    result = validate_repository(parse_args(argv).root)
    print("[PASS] resolved DOI/URL, manifest, archive, README, and license gates")
    print(json.dumps(result, indent=2, sort_keys=True))
    print("[STORAGE WARNING] Full release: "
          f"{result['download_bytes'] / 1e9:.2f} GB archives + "
          f"{result['tree_bytes'] / 1e9:.2f} GB extracted = "
          f"{(result['download_bytes'] + result['tree_bytes']) / 1e9:.2f} GB before environments and work space.")
    print("[NEXT] Plan for 700 GB free on the clone filesystem + 100 GB in /tmp "
          "(800 GB if shared); these are planning allowances.")
    print("[NEXT] Check this machine: python3 reproduce.py --check-only")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PublicationError as exc:
        print(f"[BLOCKED] {exc}", file=sys.stderr)
        raise SystemExit(2)
