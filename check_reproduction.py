#!/usr/bin/env python3
"""Check this clone's storage and system before download or figure reproduction.

Uses only the Python standard library. No data are downloaded or environments
installed. Temporary filesystem probes are removed immediately. --check-zenodo
makes four small HEAD requests; --stage reproduce runs the shipped readiness
checks, which may create review preparation/cache files but do not render figures.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import urllib.request

if sys.version_info < (3, 10):
    print("[FAIL] This checker requires Python 3.10 or newer.", file=sys.stderr)
    raise SystemExit(2)

from repository_tools import PublicationError, human_size, validate_repository
from reconstruct_paper import load_manifest, validate_public_https_url

GB = 1_000_000_000
GIB = 2**30
RECOMMENDED_CLONE_BYTES = 700 * GB
RECOMMENDED_TMP_BYTES = 100 * GB
RECOMMENDED_RAM_BYTES = 32 * GIB


class Report:
    def __init__(self):
        self.failures = 0
        self.warnings = 0

    def emit(self, status, message):
        self.failures += status == "FAIL"
        self.warnings += status == "WARN"
        print(f"[{status}] {message}", flush=True)


def run(command, **kwargs):
    return subprocess.run(command, text=True, capture_output=True, timeout=60,
                          check=False, **kwargs)


def filesystem_probe(directory):
    """Prove writes, case distinction and executable POSIX modes in our own probe."""
    with tempfile.TemporaryDirectory(prefix=".paper-preflight-", dir=directory) as tmp:
        probe = Path(tmp) / "CaseSensitive"
        probe.write_text("#!/bin/sh\nexit 0\n")
        if (Path(tmp) / "casesensitive").exists():
            raise ValueError("a case-sensitive filesystem is required")
        probe.chmod(0o751)
        if stat.S_IMODE(probe.stat().st_mode) != 0o751:
            raise ValueError("filesystem must preserve POSIX executable modes")
        result = run([str(probe)])
        if result.returncode:
            raise ValueError("filesystem must allow executing scripts (check noexec mounts)")


def memory_available():
    """Account for the host and, where exposed, container memory limits."""
    try:
        fields = dict(line.split(":", 1) for line in Path("/proc/meminfo").read_text().splitlines())
        total = int(fields["MemTotal"].split()[0]) * 1024
        available = int(fields["MemAvailable"].split()[0]) * 1024
        for limit_path, used_path in (
            ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
            ("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes"),
        ):
            try:
                limit = int(Path(limit_path).read_text())
                used = int(Path(used_path).read_text())
                total = min(total, limit)
                available = min(available, max(0, limit - used))
            except (OSError, ValueError):
                pass
        return total, available
    except (OSError, ValueError, KeyError):
        return None


def cpu_available():
    count = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count() or 1
    try:
        quota, period = Path("/sys/fs/cgroup/cpu.max").read_text().split()
        if quota != "max":
            count = min(count, max(0.01, int(quota) / int(period)))
    except (OSError, ValueError):
        pass
    return count


def find_conda():
    explicit = os.environ.get("CONDA_EXE")
    if explicit:
        return shutil.which(explicit)
    candidates = [shutil.which("conda")]
    candidates += [str(Path.home() / name / "bin/conda")
                   for name in ("miniforge3", "miniconda3", "anaconda3", "mambaforge")]
    candidates += ["/opt/conda/bin/conda", "/usr/local/bin/conda"]
    return next((path for path in candidates if path and os.access(path, os.X_OK)), None)


def find_environments(conda):
    prefixes = []
    if conda:
        result = run([conda, "env", "list", "--json"])
        if result.returncode:
            raise ValueError(f"Conda environment discovery failed: {result.stderr.strip()}")
        prefixes = json.loads(result.stdout).get("envs", [])
    found = {}
    for variable, name_variable, default in (
        ("PAPER_PYTHON", "PAPER_ENV", "paper_apotome_repro"),
        ("PAPER_S3_PYTHON", "PAPER_S3_ENV", "paper_apotome_s3_repro"),
    ):
        explicit = os.environ.get(variable)
        if variable == "PAPER_S3_PYTHON":
            explicit = explicit or os.environ.get("PAPER_REVIEW_PYTHON")
        if explicit:
            found[variable] = shutil.which(explicit)
        else:
            name = os.environ.get(name_variable, default)
            prefix = next((value for value in prefixes if Path(value).name == name), None)
            candidate = Path(prefix) / "bin/python" if prefix else None
            found[variable] = str(candidate) if candidate and os.access(candidate, os.X_OK) else None
    return found


def storage_budget(archive_bytes, tree_bytes, stage, have_archives, same_disk):
    """Additional payload bytes and the recommended free-space planning budget."""
    remaining = 0 if stage == "reproduce" else tree_bytes
    if stage == "download" and not have_archives:
        remaining += archive_bytes
    headroom = max(0, RECOMMENDED_CLONE_BYTES - archive_bytes - tree_bytes)
    recommended = remaining + headroom
    if same_disk:
        recommended += RECOMMENDED_TMP_BYTES
    return remaining, recommended


def check_storage(report, root, tmp, manifest, stage, archive_root):
    archive_bytes = sum(shard["compressed_bytes"] for record in manifest["records"].values()
                        for shard in record["shards"])
    tree_bytes = manifest["tree"]["payload_bytes"]
    report.emit("INFO", f"Archive transfer: {human_size(archive_bytes)}")
    report.emit("INFO", f"Extracted Paper payload: {human_size(tree_bytes)}")
    report.emit("INFO", f"Archives + payload: {human_size(archive_bytes + tree_bytes)}; environments and work space are additional.")
    if stage == "download":
        if archive_root is None:
            if (root / ".zenodo-archives.downloading").exists():
                report.emit("INFO", "Resumable downloads found. Space estimate conservatively includes a full download; the reconstructor alone validates retained bytes.")
        else:
            source_path = archive_root.expanduser().absolute()
            if source_path.is_symlink() or not source_path.is_dir():
                report.emit("FAIL", f"Archive root must be an existing non-symlink directory: {source_path}")
            source = source_path.resolve()
            if source == root / "Paper" or (root / "Paper").is_relative_to(source) or source.is_relative_to(root / "Paper"):
                report.emit("FAIL", "Archive root and reconstructed Paper directory must not overlap.")
            for record in manifest["records"].values():
                for shard in record["shards"]:
                    path = source / shard["relative_archive_path"]
                    if path.is_symlink() or path.resolve() != path or not path.is_file() or path.stat().st_size != shard["compressed_bytes"]:
                        report.emit("FAIL", f"Missing/nonregular/wrong-size cached archive: {path}")
            report.emit("INFO", f"Using existing archives at {source}; sizes checked here, SHA-256 checked during reconstruction.")
    work_root = root / "Paper" if stage == "reproduce" and (root / "Paper").is_dir() else root
    same_disk = work_root.stat().st_dev == tmp.stat().st_dev
    remaining, recommended = storage_budget(archive_bytes, tree_bytes, stage,
                                            archive_root is not None, same_disk)
    free = shutil.disk_usage(work_root).free
    report.emit("INFO", f"Reproduction filesystem {work_root}: {human_size(free)} free")
    if free <= 0 or free < remaining:
        report.emit("FAIL", f"Insufficient reproduction space: at least {human_size(remaining)} additional payload space is required.")
    elif free < recommended:
        report.emit("WARN", f"Payload fits, but recommended free space for this stage is {human_size(recommended)} including work reserves.")
    else:
        report.emit("PASS", f"Reproduction storage meets this stage's recommendation: {human_size(recommended)} free.")
    tmp_free = shutil.disk_usage(tmp).free
    report.emit("INFO", f"Launcher scratch filesystem {tmp}: {human_size(tmp_free)} free; {'shares the reproduction filesystem (reserves added together)' if same_disk else 'separate filesystem'}.")
    if tmp_free <= 0:
        report.emit("FAIL", "No free space in /tmp for required figure stages.")
    elif not same_disk and tmp_free < RECOMMENDED_TMP_BYTES:
        report.emit("WARN", "Recommend at least 100 GB free in /tmp for figure staging; measured peak scratch use is not established.")
    else:
        report.emit("PASS", "/tmp storage reserve considered with the reproduction reserve.")
    if hasattr(os, "statvfs") and stage == "download":
        vfs = os.statvfs(root)
        needed = len(manifest["tree"]["files"]) + len(manifest["tree"]["directories"]) + 100
        if vfs.f_files and vfs.f_favail < needed:
            report.emit("FAIL", f"Too few free filesystem inodes: {vfs.f_favail}; need at least {needed} for reconstruction.")
    report.emit("INFO", "Storage/RAM recommendations are planning allowances, not measured minima or a guarantee; account separately for Conda caches, quotas and concurrent jobs.")


def probe_zenodo(segment):
    url = segment["url"]
    validate_public_https_url(url)
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Paper-Reproduction-Preflight/1", "Accept-Encoding": "identity"})
    try:
        with urllib.request.urlopen(request, timeout=20) as response:
            validate_public_https_url(response.geturl())
            if response.status != 200:
                return "FAIL", f"Zenodo HTTP {response.status}: {url}"
            length = response.headers.get("Content-Length")
            if length is None:
                return "WARN", f"Zenodo reachable but HEAD has no Content-Length: {url}"
            if int(length) != segment["size_bytes"]:
                return "FAIL", f"Zenodo Content-Length differs from manifest: {url}"
        return "PASS", f"Public Zenodo file reachable with expected size (HEAD only): {url}"
    except Exception as exc:
        return "FAIL", f"Zenodo access failed: {url}: {exc}; check HTTPS/proxy access and retry."


def check_zenodo(report, root):
    urls = json.loads((root / "zenodo-urls.json").read_text())
    representatives = {}
    for entry in urls["archives"].values():
        for segment in entry["segments"]:
            record = segment["url"].split("/records/", 1)[1].split("/", 1)[0]
            representatives.setdefault(record, segment)
    report.emit("INFO", f"Testing one public file from each of {len(representatives)} physical Zenodo records. This does not verify all 244 chunks or their hashes.")
    with ThreadPoolExecutor(max_workers=4) as pool:
        for status, message in pool.map(probe_zenodo, representatives.values()):
            report.emit(status, message)


def check_replay_threads(report, root):
    """Require the explicit numerical profile used by the frozen CPU CNN."""
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        value = os.environ.get(name)
        report.emit("PASS" if value == "64" else "FAIL",
                    f"{name}={value or 'unset'}; the frozen Figure 5 CPU CNN requires the tested 64-thread profile.")
    report.emit("INFO", f'Set the profile in your reproduction shell: source "{root / "reproduction-runtime.env"}". This is a numerical thread setting, not a 64-core hardware minimum.')


def check_scientific_environment(report, interpreter, spec, modules, overrides=None):
    """Exercise imports and report drift from the released pip/Python pins."""
    probe = """
import importlib, importlib.metadata, json, platform, sys
errors = {}
for module in json.loads(sys.argv[1]):
    try:
        importlib.import_module(module)
    except Exception as exc:
        errors[module] = str(exc)
versions = {}
for package in json.loads(sys.argv[2]):
    try:
        versions[package] = importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        versions[package] = None
threads = sys.modules['torch'].get_num_threads() if 'torch' in sys.modules and 'torch' not in errors else None
print(json.dumps({'python': platform.python_version(), 'errors': errors, 'versions': versions, 'torch_threads': threads}))
"""
    specification = spec.read_text()
    pins = dict(re.findall(r"^\s*-\s+([A-Za-z0-9_.-]+)==([^\s]+)\s*$", specification, re.MULTILINE))
    for name, version in (overrides or {}).items():
        report.emit("INFO", f"Documented replay compatibility pin: {name}=={version} (archived {spec.name}: {pins.get(name)}); see requirements-replay-compatibility.txt.")
        pins[name] = version
    result = run([interpreter, "-c", probe, json.dumps(modules), json.dumps(list(pins))])
    if result.returncode:
        report.emit("FAIL", f"Scientific interpreter probe failed: {interpreter}: {result.stderr.strip()}")
        return
    data = json.loads(result.stdout.strip().splitlines()[-1])
    if "torch" in modules and "torch" not in data["errors"]:
        report.emit("PASS" if data["torch_threads"] == 64 else "FAIL",
                    f"Actual PyTorch numerical threads: {data['torch_threads']}; need 64 for the frozen Figure 5 CNN.")
    for module, message in data["errors"].items():
        report.emit("FAIL", f"Cannot import {module} with {interpreter}: {message}")
    if not data["errors"]:
        report.emit("PASS", f"All {len(modules)} scientific imports succeeded: {interpreter}")
    python_pin = re.search(r"^\s*-\s+python=([^\s]+)\s*$", specification, re.MULTILINE)
    if python_pin and not (data["python"] == python_pin[1] or data["python"].startswith(python_pin[1] + ".")):
        report.emit("WARN", f"Python version {data['python']} differs from {spec.name} pin {python_pin[1]}.")
    drift = [f"{name}: installed {data['versions'][name]}, pinned {version}"
             for name, version in pins.items() if data["versions"][name] != version]
    if drift:
        for message in drift:
            report.emit("WARN", f"Environment version drift ({spec.name}): {message}")
    else:
        report.emit("PASS", f"All {len(pins)} effective pip pins match {spec.name} plus any documented replay correction; Python {data['python']}.")


def replay_overrides(root):
    path = root / "requirements-replay-compatibility.txt"
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"Required replay compatibility specification missing or linked: {path}")
    rows = [line.strip() for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if rows != ["opencv-python-headless==4.10.0.84"]:
        raise ValueError("Replay compatibility specification differs from the validated v1.0 OpenCV correction.")
    return {"opencv-python-headless": "4.10.0.84"}


def check_hil_polygons(report, interpreter, paper):
    """Exercise the renderer's polygon operation, not just the import/version."""
    probe = """
import cv2, json, numpy as np, sys
from pathlib import Path
root = Path(sys.argv[1]) / 'Fig3/hil_review/ventricle_cartoon_20260827_v1'
rows = []
for name in ('water', 'allulose'):
    receipt = json.loads((root / (name + '_ventricle_receipt.json')).read_text())
    image = cv2.imread(str(root / (name + '_ventricle_mask.png')), cv2.IMREAD_UNCHANGED)
    if image is None or image.ndim != 2:
        raise ValueError('Missing or invalid accepted mask: ' + name)
    vertices = np.rint(np.asarray(receipt['polygon_vertices_xy'], dtype=np.float64)).astype(np.int32)
    raster = np.zeros(image.shape, dtype=np.uint8)
    cv2.fillPoly(raster, [vertices.reshape(-1, 1, 2)], 1)
    rows.append({'condition': name, 'differing_pixels': int(np.count_nonzero((raster > 0) != (image > 0)))})
print(json.dumps({'opencv': cv2.__version__, 'masks': rows}))
"""
    result = run([interpreter, "-c", probe, str(paper)])
    if result.returncode:
        report.emit("FAIL", f"Figure 3 HIL runtime check failed with {interpreter}: {result.stderr.strip()}")
        return
    data = json.loads(result.stdout.strip().splitlines()[-1])
    for mask in data["masks"]:
        count = mask["differing_pixels"]
        report.emit("FAIL" if count else "PASS", f"Figure 3 {mask['condition']} HIL polygon versus accepted mask with OpenCV {data['opencv']}: {count} differing pixels.")
    if any(mask["differing_pixels"] for mask in data["masks"]):
        report.emit("INFO", 'Apply the documented main-environment correction: "$PAPER_PYTHON" -m pip install --no-deps -r "$REPRO_REPO/requirements-replay-compatibility.txt". Preserve the accepted masks and receipts.')


def check_paper_sources(report, root, manifest):
    paper = root / "Paper"
    if paper.is_symlink() or not paper.is_dir():
        report.emit("FAIL", f"Reconstructed Paper directory is missing or linked: {paper}; complete reconstruction first.")
        return False
    from figure_updates import effective_manifest
    try:
        manifest = effective_manifest(root, manifest)
    except (OSError, ValueError, KeyError) as exc:
        report.emit("FAIL", f"Figure update verification failed: {exc}")
        return False
    # Authenticate executable source before invoking the shipped checks. Generated
    # publication outputs can legitimately differ after a successful rerender.
    checked = 0
    import hashlib
    for row in manifest["tree"]["files"]:
        if Path(row["path"]).suffix not in {".py", ".sh"}:
            continue
        path = paper / row["path"]
        if (path.is_symlink() or not path.is_file()
                or not path.resolve().is_relative_to(paper)
                or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]
                or stat.S_IMODE(path.stat().st_mode) != int(row["mode"], 8)):
            report.emit("FAIL", f"Manifest-bound script missing, changed, linked or wrong mode: {path}")
        checked += 1
    report.emit("INFO", f"Checked SHA-256 and modes of {checked} manifest-bound scripts.")
    return not report.failures


def check_paper(report, root, manifest, envs):
    if not check_paper_sources(report, root, manifest):
        return
    paper = root / "Paper"
    figure3_python = os.environ.get("FIG3_PYTHON") or envs["PAPER_PYTHON"]
    check_hil_polygons(report, figure3_python, paper)
    if report.failures:
        return
    child_env = os.environ.copy()
    child_env.update({key: value for key, value in envs.items() if value})
    child_env.setdefault("PAPER_PYTHON_CZI", child_env["PAPER_PYTHON"])
    child_env["PAPER_REVIEW_PYTHON"] = child_env["PAPER_S3_PYTHON"]
    child_env["PYTHONDONTWRITEBYTECODE"] = "1"
    command = ["bash", "./reproduce_all_figures.sh", "--check-only", "--output", str(paper)]
    report.emit("INFO", f"Running the shipped dependency, bundle, raw-data, Figure S3 and accepted-HIL checks in {paper}. Output follows; this may take several minutes.")
    result = subprocess.run(command, cwd=paper, env=child_env, check=False)
    report.emit("PASS" if result.returncode == 0 else "FAIL", f"Shipped reproduction readiness check exited {result.returncode}.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent, help="Git clone directory (default: directory containing this checker).")
    parser.add_argument("--stage", choices=("download", "reproduce"), default="download")
    parser.add_argument("--archive-root", type=Path, help="Existing canonical archive set; skip budgeting new downloads (still verified during reconstruction).")
    parser.add_argument("--check-zenodo", action="store_true", help="Probe one file per physical public Zenodo record, without downloading data.")
    parser.add_argument("--system-only", action="store_true", help="Check computer/storage before software setup; does not certify scientific readiness.")
    parser.add_argument("--strict", action="store_true", help="Treat recommendations/warnings as blockers too.")
    args = parser.parse_args(argv)
    report = Report()
    try:
        validate_repository(args.root)
        root = args.root.resolve()
        manifest = load_manifest(root / "manifest.json")
        compatibility = replay_overrides(root)
        report.emit("PASS", f"Repository release metadata and manifest are valid: {root}")
        report.emit("PASS" if sys.version_info >= (3, 10) else "FAIL", f"Bootstrap Python {platform.python_version()}; need Python >= 3.10.")
        report.emit("PASS" if platform.system() == "Linux" and struct.calcsize("P") == 8 else "FAIL", f"Host: {platform.system()} {platform.machine()}, {struct.calcsize('P') * 8}-bit Python; full reproduction requires 64-bit Linux.")
        if platform.machine().lower() not in {"x86_64", "amd64"}:
            report.emit("WARN", "Pinned explicit S3 locks target Linux x86-64. Other Linux architectures use YAML and have not been validated here.")
        cpus = cpu_available()
        report.emit("WARN" if cpus < 4 else "INFO", f"Available logical CPU capacity: {cpus:g}; 16 cores is the reference run, fewer cores take longer. A GPU is not required for portable replay.")
        memory = memory_available()
        if memory:
            total, available = memory
            report.emit("WARN" if available < RECOMMENDED_RAM_BYTES else "PASS", f"RAM: {total / GIB:.1f} GiB effective total, {available / GIB:.1f} GiB available; recommend >= 32 GiB available (not a measured minimum).")
        else:
            report.emit("WARN", "Cannot determine available RAM; budget >= 32 GiB for planning and check your scheduler/container memory limit.")
        for command in ("bash", "readlink", "mktemp", "find", "wc", "cp", "mv", "rm",
                        "date", "sha256sum", "tee", "cat", "dirname", "mkdir", "sed", "uname", "cmp", "chmod"):
            if not shutil.which(command):
                report.emit("FAIL", f"Required Linux command missing from PATH: {command}; install Bash/GNU coreutils/findutils/diffutils.")
        if shutil.which("cp"):
            cp_version = run(["cp", "--version"])
            if cp_version.returncode or "coreutils" not in cp_version.stdout:
                report.emit("FAIL", "GNU coreutils cp is required: installation uses cp -a --reflink=never.")
            else:
                report.emit("PASS", cp_version.stdout.splitlines()[0])
        if shutil.which("bash"):
            bash_version = run(["bash", "-c", 'printf "%s\\n" "$BASH_VERSION"; ((BASH_VERSINFO[0] > 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] >= 1)))'])
            report.emit("PASS" if bash_version.returncode == 0 else "FAIL",
                        f"Bash {bash_version.stdout.strip()}; launcher requires Bash >= 4.1 for coprocesses/dynamic descriptors.")
        tmp = Path("/tmp").resolve()
        probe_directories = [root, tmp]
        if args.stage == "reproduce" and (root / "Paper").is_dir() and not (root / "Paper").is_symlink():
            probe_directories.append(root / "Paper")
        for directory in probe_directories:
            try:
                filesystem_probe(directory)
                report.emit("PASS", f"Writable, case-sensitive filesystem with executable POSIX scripts: {directory}")
            except (OSError, ValueError) as exc:
                report.emit("FAIL", f"Filesystem probe failed for {directory}: {exc}")
        if args.stage == "download" and (root / "Paper").exists():
            report.emit("FAIL", "Paper/ already exists; reconstruction refuses overwrite. Use --stage reproduce to check the existing tree.")
        if (root / "Paper").is_symlink():
            report.emit("FAIL", "Paper/ must not be a symlink.")
        archives = args.archive_root
        if args.stage == "download" and archives is None and ((root / "zenodo-archives").exists() or (root / "zenodo-archives").is_symlink()):
            archives = root / "zenodo-archives"
            report.emit("INFO", "Canonical archive directory already exists: use reconstruct_paper.py --archive-root zenodo-archives instead of --url-map.")
        check_storage(report, root, tmp, manifest, args.stage, archives)
        conda = find_conda()
        envs = find_environments(conda)
        if conda:
            report.emit("PASS", f"Conda: {conda}")
        elif envs and all(envs.values()):
            report.emit("INFO", "Conda not found; using explicitly supplied scientific interpreters.")
        else:
            report.emit("WARN", "Conda not found. Install Miniforge/Miniconda and expose conda on PATH or set CONDA_EXE before creating the scientific environments.")
        for variable, interpreter in envs.items():
            if interpreter:
                report.emit("PASS", f"{variable}: {interpreter}")
            else:
                report.emit("FAIL" if args.stage == "reproduce" and not args.system_only else "INFO", f"{variable} not ready. After reconstruction run Paper/Fig1/00_create_conda_envs.sh and Paper/scripts/setup/01_create_s3_environment.sh, or set explicit interpreter paths.")
        if args.check_zenodo:
            check_zenodo(report, root)
        else:
            report.emit("INFO", "Zenodo connectivity not tested; add --check-zenodo for four public-file probes.")
        if args.stage == "reproduce" and not args.system_only:
            check_replay_threads(report, root)
            paper = root / "Paper"
            if not report.failures and paper.is_dir() and not paper.is_symlink():
                check_scientific_environment(report, envs["PAPER_PYTHON"],
                    paper / "scripts/setup/environment.yml",
                    ["numpy", "pandas", "scipy", "matplotlib", "PIL", "tifffile",
                     "imagecodecs", "skimage", "cv2", "aicspylibczi", "pypdf", "sklearn",
                     "torch", "cellpose", "openpyxl", "flask", "click", "tqdm"], compatibility)
                check_scientific_environment(report, envs["PAPER_S3_PYTHON"],
                    paper / "scripts/setup/environment_s3.yml",
                    ["cv2", "matplotlib", "numpy", "pandas", "pyarrow", "pyvips", "scipy",
                     "seaborn", "shapely", "sklearn", "statsmodels.api", "tifffile", "PIL", "pypdf"])
            check_paper(report, root, manifest, envs)
    except (PublicationError, OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        report.emit("FAIL", str(exc))
    if report.failures or (args.strict and report.warnings):
        report.emit("FAIL", f"Preflight blocked: {report.failures} failed checks, {report.warnings} warnings. Resolve these and rerun; no figures were reproduced.")
        return 2
    if args.system_only:
        report.emit("PASS", "Computer/storage check passed. Software setup and the final scientific readiness check are still required.")
    elif args.stage == "download":
        report.emit("PASS", "READY TO DOWNLOAD/RECONSTRUCT. This does not certify scientific environments or figure reproduction; run --stage reproduce after reconstruction and environment setup.")
    else:
        report.emit("PASS", "READY TO REPRODUCE in this clone's Paper/. This was a readiness check; figure reproduction success requires running the actual launcher.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
