#!/usr/bin/env python3
"""Reproduce all ten paper figures. Run with --check-only to check your computer."""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import time

if sys.version_info < (3, 10):
    print("Please install Python 3.10 or newer. See SETUP.md.", file=sys.stderr)
    raise SystemExit(2)

from check_reproduction import Report, check_paper_sources, find_conda
from reconstruct_paper import load_manifest

SUCCESS_MARKER = "All ten figures were installed, including certified Figure 3 (--force)."
SETUP_FILES = (
    "scripts/setup/environment.yml",
    "scripts/setup/environment_s3_conda_explicit_linux-64.txt",
    "scripts/setup/requirements_s3_pip_lock.txt",
)


class ReproductionError(RuntimeError):
    pass


class Runner:
    def __init__(self, root, archive_root=None, environ=None):
        self.root = root.resolve()
        self.paper = self.root / "Paper"
        self.archives = archive_root
        self.env = dict(os.environ if environ is None else environ)
        self.env.update({name: "64" for name in
                         ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")})
        self.env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1", CONDA_ALWAYS_YES="true")
        self.logs = self.root / "logs"
        self.local = self.root / ".paper-envs"
        self.conda = None

    def command(self, title, command, log_name, cwd=None):
        """Keep full diagnostics on disk and show a heartbeat during long steps."""
        self.logs.mkdir(exist_ok=True)
        log = self.logs / log_name
        print(f"\n{title}\n  Details: {log}", flush=True)
        started = time.monotonic()
        # Append so an interrupted attempt's diagnostics remain available.
        with log.open("a", encoding="utf-8") as output:
            output.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n{shlex.join(map(str, command))}\n")
            output.flush()
            offset = output.tell()
            child = subprocess.Popen(list(map(str, command)), cwd=cwd or self.root,
                                     env=self.env, stdout=output, stderr=subprocess.STDOUT,
                                     stdin=subprocess.DEVNULL, start_new_session=True)
            try:
                while True:
                    try:
                        code = child.wait(timeout=30)
                        break
                    except subprocess.TimeoutExpired:
                        minutes = (time.monotonic() - started) / 60
                        print(f"  Still working ({minutes:.0f} min). Full progress is in the log.", flush=True)
            except BaseException:
                with contextlib.suppress(ProcessLookupError):
                    os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    with contextlib.suppress(ProcessLookupError):
                        os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                raise
        with log.open(encoding="utf-8", errors="replace") as source:
            source.seek(offset)
            current_output = source.read()
        if code:
            lines = current_output.splitlines()
            problems = [line for line in lines if line.startswith(("[FAIL]", "[WARN]"))]
            for line in (problems or lines[-8:])[-12:]:
                print(f"  {line}", flush=True)
            raise ReproductionError(f"This step stopped (exit {code}). See {log}.\nFix the reported problem, then run the same command again.")
        print("  Done.", flush=True)
        return current_output

    def check_computer(self):
        self.conda = find_conda()
        if self.conda:
            self.env["CONDA_EXE"] = self.conda
        main = self.env.get("PAPER_PYTHON")
        s3 = self.env.get("PAPER_S3_PYTHON") or self.env.get("PAPER_REVIEW_PYTHON")
        for label, value in (("PAPER_PYTHON", main), ("PAPER_S3_PYTHON", s3)):
            if value:
                executable = shutil.which(value, path=self.env.get("PATH"))
                if not executable:
                    raise ReproductionError(f"{label} points to a missing program: {value}. Correct or unset it and retry.")
                # Keep a virtual environment's bin/python symlink; resolving its
                # target can silently select the base interpreter instead.
                self.env[label] = os.path.abspath(executable)
        if not self.conda and not (main and s3):
            raise ReproductionError("Conda is missing. Follow SETUP.md to install Miniforge, then run this command again. No paper data were downloaded.")
        if self.archives is None and (self.root / "zenodo-archives").exists():
            self.archives = self.root / "zenodo-archives"
        stage = "reproduce" if self.paper.exists() else "download"
        command = [sys.executable, self.root / "check_reproduction.py", "--root", self.root,
                   "--stage", stage, "--system-only", "--strict"]
        if stage == "download" and self.archives is not None:
            command += ["--archive-root", self.archives]
        output = self.command("1/5 Checking your computer and free space", command, "01-computer.log")
        for line in output.splitlines():
            if "Host:" in line:
                print("  Computer: " + line.split("Host: ", 1)[1].split(";", 1)[0], flush=True)
            elif "logical CPU capacity:" in line:
                count = line.split("capacity: ", 1)[1].split(";", 1)[0]
                print(f"  Available CPU capacity: {count}. A GPU is not needed.", flush=True)
            elif "RAM:" in line:
                memory = re.search(r"([0-9.]+) GiB available", line)
                if memory:
                    print(f"  Available memory: {memory[1]} GiB (32 GiB recommended).", flush=True)
            elif "Reproduction filesystem" in line or "Launcher scratch filesystem" in line:
                space = re.search(r"\(([0-9.]+) GB;", line)
                if space:
                    location = "Project drive" if "Reproduction filesystem" in line else "Temporary work drive (/tmp)"
                    print(f"  {location}: {space[1]} GB free.", flush=True)

    def reconstruct(self):
        if self.paper.exists():
            print("\n2/5 Paper data are already present; checking them before use.", flush=True)
        else:
            command = [sys.executable, self.root / "reconstruct_paper.py",
                       "--manifest", self.root / "manifest.json", "--output", self.paper, "--verbose"]
            if self.archives is not None:
                command += ["--archive-root", self.archives]
            else:
                command += ["--url-map", self.root / "zenodo-urls.json",
                            "--download-dir", self.root / "zenodo-archives"]
            self.command("2/5 Downloading and unpacking the paper data (about 240 GB for a new download)",
                         command, "02-download.log")
        # Authenticate original programs and installer specifications before use.
        self.logs.mkdir(exist_ok=True)
        log = self.logs / "02-source-check.log"
        manifest = load_manifest(self.root / "manifest.json")
        with log.open("w") as output, contextlib.redirect_stdout(output):
            report = Report()
            valid = check_paper_sources(report, self.root, manifest)
        if not valid:
            raise ReproductionError(f"Some original paper programs are missing or changed. See {log} and DETAILED_GUIDE.md.")
        rows = {row["path"]: row for row in manifest["tree"]["files"]}
        for relative in SETUP_FILES:
            path = self.paper / relative
            if (relative not in rows or path.is_symlink() or not path.is_file()
                    or not path.resolve().is_relative_to(self.paper)
                    or hashlib.sha256(path.read_bytes()).hexdigest() != rows[relative]["sha256"]):
                raise ReproductionError(f"Original software specification is missing or changed: {path}. See DETAILED_GUIDE.md.")

    def install_environment(self, name, specifications):
        prefix = self.local / name
        python = prefix / "bin/python"
        receipt = self.local / f"{name}-prepared.json"
        signature = {str(path.relative_to(self.root)): hashlib.sha256(path.read_bytes()).hexdigest()
                     for path in specifications}
        if python.is_file() and receipt.is_file():
            try:
                if json.loads(receipt.read_text()) == signature:
                    return python
            except (OSError, ValueError):
                pass
        if not self.conda:
            raise ReproductionError("Conda is needed to install the project software. See SETUP.md.")
        self.local.mkdir(exist_ok=True)
        self.env["CONDA_PKGS_DIRS"] = str(self.local / "cache/conda")
        self.env["PIP_CACHE_DIR"] = str(self.local / "cache/pip")
        exists = (prefix / "conda-meta/history").is_file()
        if name == "main":
            command = [self.conda, "env", "update" if exists else "create", "--prefix", prefix,
                       "--file", specifications[0]]
        else:
            command = [self.conda, "install" if exists else "create", "--prefix", prefix,
                       "--file", specifications[0], "--yes"]
        self.command(f"3/5 Installing {name} software (first run can take a while)",
                     command, f"03-software-{name}.log")
        self.command(f"3/5 Applying the tested {name} package versions",
                     [python, "-m", "pip", "install", "--no-deps", "-r", specifications[1]],
                     f"03-software-{name}.log")
        receipt.write_text(json.dumps(signature, indent=2) + "\n")
        return python

    def setup_software(self):
        main = self.env.get("PAPER_PYTHON")
        s3 = self.env.get("PAPER_S3_PYTHON") or self.env.get("PAPER_REVIEW_PYTHON")
        if not main:
            main = self.install_environment("main", [self.paper / SETUP_FILES[0],
                                                   self.root / "requirements-replay-compatibility.txt"])
        if not s3:
            s3 = self.install_environment("s3", [self.paper / SETUP_FILES[1], self.paper / SETUP_FILES[2]])
        self.env.update(PAPER_PYTHON=str(main), PAPER_PYTHON_CZI=str(main), FIG3_PYTHON=str(main),
                        PAPER_S3_PYTHON=str(s3), PAPER_REVIEW_PYTHON=str(s3))
        print("\n3/5 Project software is installed. Checking it next.", flush=True)

    def readiness(self):
        self.command("4/5 Checking the data, software, and figure settings",
                     [sys.executable, self.root / "check_reproduction.py", "--root", self.root,
                      "--stage", "reproduce", "--strict"], "04-readiness.log")

    def figures(self):
        output = self.command("5/5 Reproducing and verifying all 10 figures",
                              ["bash", "./reproduce_all_figures.sh", "--output", self.paper, "--force"],
                              "05-figures.log", cwd=self.paper)
        if SUCCESS_MARKER not in output.splitlines():
            raise ReproductionError(f"The figure program did not confirm completion. Inspect {self.logs / '05-figures.log'} before retrying.")
        print(f"\nSUCCESS: All 10 figures reproduced and verified.\nFind the results in {self.paper}/Fig1 through Fig5 and FigS1 through FigS5.", flush=True)

    def execute(self, check_only=False, prepare_only=False):
        self.check_computer()
        if check_only:
            print("\nComputer check passed. No data were downloaded or software installed.\nRun python3 reproduce.py to reproduce the paper.", flush=True)
            return
        self.reconstruct()
        self.setup_software()
        self.readiness()
        if prepare_only:
            print("\nPreparation passed. Figures have not been rerun.\nRun python3 reproduce.py to reproduce the paper.", flush=True)
            return
        self.figures()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--check-only", action="store_true", help="Check the computer and space without downloading or installing anything.")
    modes.add_argument("--prepare-only", action="store_true", help="Download, install software, and check readiness; stop before running figures.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent, help="Repository folder (defaults to this script's folder).")
    parser.add_argument("--archive-root", type=Path, help="Advanced: reuse an existing complete set of downloaded ZIP archives.")
    args = parser.parse_args(argv)
    try:
        if not sys.platform.startswith("linux"):
            raise ReproductionError("This workflow needs a 64-bit Linux computer. See SETUP.md.")
        import fcntl
        with (args.root / ".reproduce.lock").open("a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise ReproductionError("Another reproduction is already running in this folder. Use its terminal to follow progress.") from None
            print("Paper reproduction: all 10 figures\nA new download is about 240 GB; plan for 800 GB free space including work files.", flush=True)
            runner = Runner(args.root, args.archive_root.expanduser().absolute() if args.archive_root else None)
            runner.execute(check_only=args.check_only, prepare_only=args.prepare_only)
    except KeyboardInterrupt:
        print("\nStopped. Run the same command again to continue with completed downloads.", file=sys.stderr)
        return 130
    except (ReproductionError, OSError, ValueError) as exc:
        print(f"\nSTOPPED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
