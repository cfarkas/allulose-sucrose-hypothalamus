#!/usr/bin/env python3
"""Check source hashes, script syntax and figure/script mirror agreement."""
import ast
from pathlib import Path
import subprocess
from figure_updates import load_updates, PAYLOAD


def main():
    root = Path(__file__).resolve().parent
    update = load_updates(root)
    payload = root / PAYLOAD
    scripts = 0
    mirrors = 0
    for row in update["files"]:
        path = payload / row["path"]
        if path.suffix == ".py":
            ast.parse(path.read_text(), filename=row["path"])
            scripts += 1
        elif path.suffix == ".sh":
            subprocess.run(["bash", "-n", path], check=True)
            scripts += 1
    for figure in ("Fig1", "Fig4", "Fig5", "FigS6", "FigS7", "FigS8"):
        for path in (payload / figure).iterdir():
            if path.suffix not in {".py", ".sh"}:
                continue
            mirror = (payload / "scripts/setup" / path.name if path.name == "00_create_conda_envs.sh"
                      else payload / "scripts" / figure / path.name)
            if not mirror.is_file() or path.read_bytes() != mirror.read_bytes():
                raise ValueError(f"Missing or different script mirror: {figure}/{path.name}")
            mirrors += 1
    print(f"PASS: {len(update['files'])} manifest files, {scripts} script syntax checks, {mirrors} exact mirrors.")


if __name__ == "__main__":
    main()
