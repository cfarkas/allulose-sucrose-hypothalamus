#!/usr/bin/env python3
"""List or safely download model-release assets (Python 3.9+, POSIX).

Examples:
  python3 download_models.py
  python3 download_models.py --asset cellpose__cfos_channel3 --destination Paper
  python3 download_models.py --all --destination Paper

Paths come from the adjacent manifest.json. Downloads require an explicit selection;
existing files must already match. The helper never loads serialized model objects.
"""
import argparse
import hashlib
from http.client import HTTPException
import json
import os
from pathlib import Path, PureWindowsPath
import re
import secrets
import stat
import sys
from urllib.parse import urlparse
from urllib.request import Request, urlopen

CHUNK = 1024 * 1024


def relative_parts(value):
    if not isinstance(value, str) or not value or "\\" in value or "\x00" in value:
        raise ValueError(f"Invalid relative path: {value!r}")
    parts = value.split("/")
    if any(part in ("", ".", "..") for part in parts) or PureWindowsPath(value).drive:
        raise ValueError(f"Unsafe relative path: {value!r}")
    return parts


def read_manifest(path):
    with open(path, encoding="utf-8") as stream:
        manifest = json.load(stream)
    if manifest.get("schema_version") != 1 or not isinstance(manifest.get("assets"), list):
        raise ValueError("Expected schema_version 1 and an assets list")
    names, paths = {}, set()
    for asset in manifest["assets"]:
        if not isinstance(asset, dict):
            raise ValueError("Each asset must be an object")
        name = asset.get("asset_name")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError(f"Missing or duplicate asset name: {name!r}")
        size = asset.get("size_bytes")
        digest = asset.get("sha256", "")
        if type(size) is not int or size < 0 or not isinstance(digest, str) or not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(f"Invalid size or SHA256 for {name}")
        url = urlparse(asset.get("download_url", ""))
        if url.scheme not in ("https", "http") or not url.hostname or url.username or url.password:
            raise ValueError(f"Invalid HTTP(S) download URL for {name}")
        aliases = asset.get("aliases", [])
        if not isinstance(aliases, list):
            raise ValueError(f"Invalid aliases for {name}")
        names[name] = asset
        for target in [asset] + aliases:
            if not isinstance(target, dict):
                raise ValueError(f"Invalid alias object for {name}")
            relative_parts(target.get("restore_path"))
            path_value = target["restore_path"]
            if path_value in paths:
                raise ValueError(f"Duplicate restore path: {path_value}")
            paths.add(path_value)
            if target is not asset:
                alias_name = target.get("asset_name")
                if not isinstance(alias_name, str) or not alias_name or alias_name in names:
                    raise ValueError(f"Missing or duplicate alias name: {alias_name!r}")
                names[alias_name] = asset
    return manifest, names


def parent_directory(root_fd, relative_path):
    """Walk beneath an open root; reject symlinks, including concurrent swaps."""
    parts = relative_parts(relative_path)
    current = os.dup(root_fd)
    try:
        for component in parts[:-1]:
            try:
                os.mkdir(component, mode=0o755, dir_fd=current)
            except FileExistsError:
                pass
            child = os.open(component, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=current)
            os.close(current)
            current = child
        return current, parts[-1]
    except BaseException:
        os.close(current)
        raise


def verify_stream(stream, asset):
    digest, total = hashlib.sha256(), 0
    for block in iter(lambda: stream.read(CHUNK), b""):
        total += len(block)
        if total > asset["size_bytes"]:
            raise ValueError(f"Size mismatch for {asset['asset_name']}")
        digest.update(block)
    if total != asset["size_bytes"] or digest.hexdigest() != asset["sha256"]:
        raise ValueError(f"Size/SHA256 mismatch for {asset['asset_name']}")


def existing_valid(directory, name, asset):
    try:
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory)
    except FileNotFoundError:
        return False
    with os.fdopen(fd, "rb") as stream:
        if not stat.S_ISREG(os.fstat(stream.fileno()).st_mode):
            raise ValueError(f"Existing target is not a regular file: {name}")
        try:
            verify_stream(stream, asset)
        except ValueError as exc:
            raise ValueError(f"Refusing to overwrite different existing file: {name}") from exc
    return True


def install(directory, name, asset, source):
    """Write a verified temporary file and atomically link without overwriting."""
    temporary = f".download-{secrets.token_hex(16)}.part"
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o644, dir_fd=directory)
    try:
        with os.fdopen(fd, "wb") as output:
            digest, total = hashlib.sha256(), 0
            for block in iter(lambda: source.read(CHUNK), b""):
                total += len(block)
                if total > asset["size_bytes"]:
                    raise ValueError(f"Downloaded size exceeds manifest for {asset['asset_name']}")
                output.write(block)
                digest.update(block)
            if total != asset["size_bytes"] or digest.hexdigest() != asset["sha256"]:
                raise ValueError(f"Downloaded size/SHA256 mismatch for {asset['asset_name']}")
            output.flush()
            os.fsync(output.fileno())
        try:
            os.link(temporary, name, src_dir_fd=directory, dst_dir_fd=directory, follow_symlinks=False)
        except FileExistsError:
            # Another downloader may have completed this same asset meanwhile.
            if not existing_valid(directory, name, asset):
                raise ValueError(f"Target disappeared during installation: {name}; retry")
    finally:
        os.unlink(temporary, dir_fd=directory)


def download_asset(root_fd, asset):
    parent, name = parent_directory(root_fd, asset["restore_path"])
    try:
        if existing_valid(parent, name, asset):
            print(f"Verified existing: {asset['restore_path']}")
        else:
            print(f"Downloading: {asset['asset_name']} ({asset['size_bytes']} bytes)", flush=True)
            request = Request(asset["download_url"], headers={"User-Agent": "ScientificReports-model-download/1.0"})
            with urlopen(request, timeout=120) as source:
                install(parent, name, asset, source)
            print(f"Installed: {asset['restore_path']}")
        for alias in asset.get("aliases", []):
            alias_parent, alias_name = parent_directory(root_fd, alias["restore_path"])
            try:
                if not existing_valid(alias_parent, alias_name, asset):
                    source_fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=parent)
                    with os.fdopen(source_fd, "rb") as source:
                        install(alias_parent, alias_name, asset, source)
                print(f"Verified alias: {alias['restore_path']}")
            finally:
                os.close(alias_parent)
    finally:
        os.close(parent)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).resolve().with_name("manifest.json"))
    parser.add_argument("--destination", type=Path, default=Path.cwd(), help="Root for manifest restore paths (default: current directory)")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--list", action="store_true", help="List assets and aliases; download nothing (default)")
    action.add_argument("--all", action="store_true", help="Download every manifest asset")
    action.add_argument("--asset", action="append", metavar="NAME", help="Download named asset or alias; repeat to select several")
    args = parser.parse_args(argv)
    try:
        manifest, names = read_manifest(args.manifest)
        if not args.all and not args.asset:
            for asset in manifest["assets"]:
                print(f"{asset['asset_name']}\t{asset['size_bytes']} bytes\t{asset.get('role', '')}\t{asset['restore_path']}")
                for alias in asset.get("aliases", []):
                    print(f"  alias: {alias['asset_name']} -> {asset['asset_name']}")
            print(f"{len(manifest['assets'])} downloadable assets; {sum(a['size_bytes'] for a in manifest['assets'])} bytes total")
            return 0
        selected = manifest["assets"] if args.all else []
        for requested in args.asset or []:
            if requested not in names:
                raise ValueError(f"Unknown asset: {requested}; use --list")
            if names[requested] not in selected:
                selected.append(names[requested])
        if not all(hasattr(os, name) for name in ("O_NOFOLLOW", "O_DIRECTORY")) or os.name != "posix":
            raise ValueError("Safe downloading requires POSIX directory-descriptor support; listing works on other systems")
        args.destination.mkdir(parents=True, exist_ok=True)
        root = args.destination.resolve(strict=True)
        root_fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            # Check existing targets, including aliases, before transferring any bytes.
            for asset in selected:
                for target in [asset] + asset.get("aliases", []):
                    parent, name = parent_directory(root_fd, target["restore_path"])
                    try:
                        existing_valid(parent, name, asset)
                    finally:
                        os.close(parent)
            for asset in selected:
                download_asset(root_fd, asset)
        finally:
            os.close(root_fd)
        return 0
    except (OSError, ValueError, KeyError, TypeError, HTTPException) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
