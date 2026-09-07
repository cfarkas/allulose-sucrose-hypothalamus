#!/usr/bin/env python3
"""Portability and reproducibility self-check for the whole Paper tree.

Run this first on any machine, and before publishing the tree, to see which
figures can be rebuilt here. It reads and parses; it never executes a figure
step and never writes outside its own report.

Checks, per figure:
  * every front NN_*.py parses and every front NN_*.sh passes bash -n
  * no script resolves the Paper root by fixed depth, which breaks as soon as a
    script is moved or the tree is relocated
  * no script resolves it by the directory name alone either, which breaks as
    soon as the archive is unpacked under any name other than Paper
  * no script opens a hard-coded absolute path outside the tree
  * front scripts are byte-identical to their scripts/<Fig> mirrors
  * every sibling script a runner invokes actually exists
  * the third-party modules the figure needs are importable
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = next(path for path in HERE.parents if (path.name == "Paper" or ((path / "scripts" / "setup").is_dir()
                              and (path / "README.txt").is_file())))
FIGURES = ("Fig1", "Fig2", "Fig3", "Fig4", "Fig5", "FigS1", "FigS2", "FigS3", "FigS4", "FigS5", "FigS6")

FIXED_DEPTH = re.compile(r"\.parents\[\d+\]")
# Chained .parent walks are the same defect wearing a different hat, and the
# .parents[N] rule never saw them. Fig2's staging scripts derived the paper root
# as HERE.parent.parent, which silently resolved to Paper/scripts when they were
# run from their scripts/Fig2 mirror, and every raw path under it then pointed
# somewhere that does not exist.
CHAINED_PARENT = re.compile(r"=\s*\w+\.parent\.parent\b")
DERIVED_PARENT = re.compile(r"^\s*(\w+)\s*=\s*(\w+)\.parent\s*$")
SHELL_UP = re.compile(r"\$\{?SCRIPT_DIR\}?/\.\./\.\.|\$\{?script_dir\}?/\.\./\.\.")
ABSOLUTE = re.compile(r"[\"'](/(?:media|home|mnt|data)/[^\"']+)[\"']")
INSIDE_TREE = re.compile(r"/Apotome/Paper(?:/|\"|')")

MODULES = {
    "FigS6": ("numpy", "pandas", "matplotlib", "PIL"),
    "Fig1": ("numpy", "pandas", "matplotlib", "scipy"),
    "Fig2": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile", "aicspylibczi", "cv2", "skimage"),
    "Fig3": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile"),
    "Fig4": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile"),
    "Fig5": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile", "skimage"),
    "FigS5": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile", "aicspylibczi", "skimage"),
    "FigS1": ("numpy", "pandas", "matplotlib", "scipy"),
    "FigS2": ("numpy", "pandas", "matplotlib", "scipy"),
    "FigS4": ("numpy", "pandas", "matplotlib", "scipy"),
    "FigS3": ("numpy", "pandas", "matplotlib", "scipy", "PIL", "tifffile", "skimage", "sklearn", "statsmodels"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def front_scripts(figure: str) -> list[Path]:
    directory = PAPER / figure
    return sorted(p for p in directory.glob("[0-9]*") if p.suffix in {".py", ".sh"} and p.is_file())


def check_figure(figure: str) -> dict:
    problems: list[str] = []
    notes: list[str] = []
    scripts = front_scripts(figure)
    if not scripts:
        problems.append("no numbered scripts at the front of the figure directory")

    for path in scripts:
        relative = path.relative_to(PAPER).as_posix()
        text = path.read_text(encoding="utf-8", errors="replace")

        if path.suffix == ".py":
            try:
                ast.parse(text)
            except SyntaxError as error:
                problems.append(f"{relative}: syntax error line {error.lineno}")
                continue
        else:
            result = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
            if result.returncode != 0:
                problems.append(f"{relative}: bash -n failed: {result.stderr.strip()}")

        derived_from: dict[str, str] = {}
        for line_number, line in enumerate(text.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if FIXED_DEPTH.search(line) and "parents" in line and "Paper" not in line:
                problems.append(f"{relative}:{line_number}: resolves a root by fixed depth")
            if CHAINED_PARENT.search(line):
                problems.append(
                    f"{relative}:{line_number}: resolves a root by chained .parent.parent, "
                    "which breaks as soon as the script is run from its mirror"
                )
            match = DERIVED_PARENT.match(line)
            if match:
                name, source = match.group(1), match.group(2)
                if derived_from.get(source) == "file":
                    problems.append(
                        f"{relative}:{line_number}: {name} is two .parent steps from the "
                        "script, so it resolves to Paper/scripts when run from the mirror"
                    )
                derived_from[name] = "file" if source in {"HERE", "_HERE"} else derived_from.get(source, "")
            if SHELL_UP.search(line):
                problems.append(f"{relative}:{line_number}: shell root derived by fixed ../..")
            for match in ABSOLUTE.finditer(line):
                candidate = match.group(1)
                if INSIDE_TREE.search(candidate + "/"):
                    notes.append(f"{relative}:{line_number}: absolute in-tree path {candidate}")
                else:
                    notes.append(f"{relative}:{line_number}: external path {candidate}")

        # A tree unpacked as anything but "Paper" must still be found. The name is
        # tried first and is still the usual answer; the marker directory that
        # holds README.txt and scripts/setup is what makes the name optional.
        if path.suffix == ".py" and 'name == "Paper"' in text and '"scripts" / "setup"' not in text:
            problems.append(
                f"{relative}: locates the tree only by the directory name Paper, so a "
                "downloaded archive unpacked under any other name cannot run it"
            )

        # every numbered script a runner names must resolve somewhere in the tree
        for reference in sorted(set(re.findall(r"[\"'\$\{/]([0-9]{2}[a-z]?_[A-Za-z0-9_]+\.(?:py|sh))", text))):
            if reference == path.name:
                continue
            if (path.parent / reference).is_file():
                continue
            elsewhere = [
                candidate.relative_to(PAPER).as_posix()
                for candidate in PAPER.rglob(reference)
                if candidate.is_file()
            ]
            if elsewhere:
                notes.append(f"{relative}: names {reference}, resolved at {elsewhere[0]}")
            else:
                problems.append(f"{relative}: names {reference}, which exists nowhere in the tree")

    # front and mirror must agree
    mirror = PAPER / "scripts" / figure
    if mirror.is_dir():
        for path in sorted(mirror.rglob("*")):
            if not path.is_file() or "__pycache__" in path.parts:
                continue
            front = PAPER / figure / path.relative_to(mirror)
            if not front.is_file():
                problems.append(f"scripts/{figure}/{path.relative_to(mirror)} has no front copy")
            elif sha256(front) != sha256(path):
                problems.append(f"{figure}/{path.relative_to(mirror)} differs from its scripts/ mirror")
    else:
        notes.append(f"no scripts/{figure} mirror; the front copies are canonical")

    missing = []
    for module in MODULES[figure]:
        try:
            __import__(module)
        except Exception:
            missing.append(module)
    if missing:
        notes.append(f"modules not importable with this interpreter: {', '.join(missing)}")

    return {
        "figure": figure,
        "scripts": [p.name for p in scripts],
        "status": "FAIL" if problems else "PASS",
        "problems": problems,
        "notes": notes,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="also write the report as JSON")
    args = parser.parse_args()

    report = [check_figure(figure) for figure in FIGURES]
    failed = [entry for entry in report if entry["status"] == "FAIL"]

    print(f"Paper root: {PAPER}")
    print(f"Interpreter: {sys.executable}\n")
    for entry in report:
        print(f"{entry['status']:4}  {entry['figure']:6}  {len(entry['scripts'])} numbered scripts")
        for problem in entry["problems"]:
            print(f"        PROBLEM  {problem}")
        for note in entry["notes"]:
            print(f"        note     {note}")
    print()
    print("OVERALL:", "PASS" if not failed else f"FAIL ({', '.join(e['figure'] for e in failed)})")

    if args.json:
        args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
