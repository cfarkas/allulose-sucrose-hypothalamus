#!/usr/bin/env python3
"""Fetch the exact public DINOv2 source/weights used in the S6 benchmark.

Only an absent destination is created. Use --verify-only to inspect an existing
asset directory without network requests or changes. This downloads assets, not
human labels or microscopy, and does not train a classifier.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tempfile
import urllib.request
import zipfile

COMMIT = "7764ea0f912e53c92e82eb78a2a1631e92725fc8"
SOURCE_URL = f"https://codeload.github.com/facebookresearch/dinov2/zip/{COMMIT}"
SOURCE_SHA256 = "04276715cddb29d45d05bff3a6fc132224dc27749b279ac98ad2ce4620e20d48"
WEIGHTS_URL = "https://dl.fbaipublicfiles.com/dinov2/dinov2_vits14/dinov2_vits14_pretrain.pth"
WEIGHTS_SHA256 = "b938bf1bc15cd2ec0feacfe3a1bb553fe8ea9ca46a7e1d8d00217f29aef60cd9"
PAPER = next(p for p in Path(__file__).resolve().parents if (p / "Fig5/09_improve_microglia_classifier.py").is_file())
DEFAULT_OUTPUT = PAPER / "analyses/Fig5/model_assets/dinov2"


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def verify(directory):
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError("Asset directory is absent or linked")
    metadata = json.loads((directory / "asset_manifest.json").read_text())
    if metadata["source_commit"] != COMMIT or metadata["weights_sha256"] != WEIGHTS_SHA256:
        raise ValueError("Asset metadata does not describe the frozen benchmark encoder")
    weights = Path(metadata["weights_path"])
    source = Path(metadata["source_dir"])
    if weights.is_symlink() or not weights.resolve().is_relative_to(directory.resolve()) or digest(weights) != WEIGHTS_SHA256:
        raise ValueError("Encoder weights differ or escape the asset directory")
    archive = directory / "source_7764ea0f.zip"
    if archive.is_symlink() or digest(archive) != SOURCE_SHA256:
        raise ValueError("Encoder source archive differs")
    if source.is_symlink() or source.resolve() != (directory / f"dinov2-{COMMIT}").resolve():
        raise ValueError("Encoder source directory differs")
    with zipfile.ZipFile(archive) as zipped:
        for member in zipped.infolist():
            if member.is_dir():
                continue
            target = directory / member.filename
            if target.is_symlink() or not target.resolve().is_relative_to(directory.resolve()) or not target.is_file() or target.read_bytes() != zipped.read(member):
                raise ValueError(f"Encoder source differs: {member.filename}")
    print("PASS: exact DINOv2 source and pretrained weights verified.")


def prepare(destination):
    if destination.exists() or destination.is_symlink():
        raise ValueError("Destination already exists; use --verify-only or choose an absent directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".dinov2-download-", dir=destination.parent) as temporary:
        stage = Path(temporary)
        for url, filename, expected in [(SOURCE_URL, "source_7764ea0f.zip", SOURCE_SHA256),
                                         (WEIGHTS_URL, "dinov2_vits14_pretrain.pth", WEIGHTS_SHA256)]:
            target = stage / filename
            print(f"Downloading {filename}", flush=True)
            with urllib.request.urlopen(url, timeout=90) as incoming, target.open("wb") as output:
                shutil.copyfileobj(incoming, output, length=1024 * 1024)
            if digest(target) != expected:
                raise ValueError(f"Downloaded SHA-256 differs: {filename}")
        with zipfile.ZipFile(stage / "source_7764ea0f.zip") as zipped:
            for member in zipped.infolist():
                path = PurePosixPath(member.filename)
                if (path.is_absolute() or ".." in path.parts or "\\" in member.filename
                        or not path.parts or path.parts[0] != f"dinov2-{COMMIT}"
                        or (member.external_attr >> 16) & 0o170000 == 0o120000):
                    raise ValueError("Unsafe source archive member")
            zipped.extractall(stage)
        metadata = {"model": "dinov2_vits14", "source_commit": COMMIT,
                    "source_url": f"https://github.com/facebookresearch/dinov2/tree/{COMMIT}",
                    "source_dir": str(destination / f"dinov2-{COMMIT}"),
                    "source_archive_sha256": SOURCE_SHA256, "weights_url": WEIGHTS_URL,
                    "weights_path": str(destination / "dinov2_vits14_pretrain.pth"),
                    "weights_sha256": WEIGHTS_SHA256, "license": "Apache-2.0"}
        (stage / "asset_manifest.json").write_text(json.dumps(metadata, indent=2) + "\n")
        # Refuse replacing a destination created during download.
        destination.mkdir()
        for path in stage.iterdir():
            shutil.move(str(path), destination / path.name)
    verify(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    destination = args.output_dir.expanduser().absolute()
    if args.verify_only:
        verify(destination)
    else:
        prepare(destination)


if __name__ == "__main__":
    main()
