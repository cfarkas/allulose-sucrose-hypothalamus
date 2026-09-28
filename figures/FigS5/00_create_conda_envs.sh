#!/usr/bin/env bash
# Step 00 for every figure: create the conda environment the workflows need.
#
#   ./Fig1/00_create_conda_envs.sh            create it if missing
#   ./Fig1/00_create_conda_envs.sh --force    remove and recreate it
#   ./Fig1/00_create_conda_envs.sh --check    only report what is present
#
# One environment, paper_apotome_repro, satisfies every figure including the
# Figure 2 CZI and Cellpose stages. The pinned specification is
# scripts/setup/environment.yml. Identical copies of this script sit at the
# front of each figure directory, so it can be run from wherever you start.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
# Locate the Paper root by walking up to the directory that holds scripts/shared,
# so this works from a figure folder, from scripts/setup, or from a relocated tree.
PAPER_ROOT="$(cd -- "$SCRIPT_DIR" && while [[ ! -d scripts/shared && "$PWD" != / ]]; do cd ..; done; pwd -P)"
[[ -d "$PAPER_ROOT/scripts/shared" ]] || {
  printf 'Cannot locate the Paper root above %s\n' "$SCRIPT_DIR" >&2
  exit 2
}

ENV_FILE="${PAPER_ROOT}/scripts/setup/environment.yml"
ENV_NAME="${PAPER_ENV:-paper_apotome_repro}"
FORCE=0
CHECK_ONLY=0

while (($#)); do
  case "$1" in
    --force) FORCE=1 ;;
    --check) CHECK_ONLY=1 ;;
    -h|--help) sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift
done

[[ -f "$ENV_FILE" ]] || { printf 'Missing environment specification: %s\n' "$ENV_FILE" >&2; exit 2; }

# Find conda: an explicit CONDA_EXE, then PATH, then the usual install roots.
CONDA_BIN=""
if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
  CONDA_BIN="$CONDA_EXE"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BIN="$(command -v conda)"
else
  for candidate in \
    "$HOME/anaconda3/bin/conda" "$HOME/miniconda3/bin/conda" \
    "$HOME/mambaforge/bin/conda" "$HOME/miniforge3/bin/conda" \
    /opt/conda/bin/conda /usr/local/bin/conda /home/server/anaconda3/bin/conda; do
    [[ -x "$candidate" ]] && { CONDA_BIN="$candidate"; break; }
  done
fi
[[ -n "$CONDA_BIN" ]] || {
  printf 'No conda found. Install Miniforge or Miniconda, or set CONDA_EXE.\n' >&2
  exit 2
}

ENV_PREFIX="$("$CONDA_BIN" env list --json 2>/dev/null \
  | "${PAPER_PYTHON:-python3}" -c '
import json, os, sys
name = sys.argv[1]
try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(0)
for prefix in data.get("envs", []):
    if os.path.basename(prefix) == name:
        print(prefix)
        break
' "$ENV_NAME" || true)"

printf 'Paper root : %s\n' "$PAPER_ROOT"
printf 'conda      : %s\n' "$CONDA_BIN"
printf 'environment: %s\n' "$ENV_NAME"

if [[ -n "$ENV_PREFIX" && "$FORCE" -eq 0 ]]; then
  printf 'status     : already present at %s\n' "$ENV_PREFIX"
elif [[ "$CHECK_ONLY" -eq 1 ]]; then
  printf 'status     : absent (run without --check to create it)\n'
  exit 1
else
  if [[ -n "$ENV_PREFIX" ]]; then
    printf 'status     : removing existing environment before recreating it\n'
    "$CONDA_BIN" env remove --name "$ENV_NAME" --yes
  fi
  printf 'status     : creating from %s\n' "$ENV_FILE"
  "$CONDA_BIN" env create --name "$ENV_NAME" --file "$ENV_FILE"
  ENV_PREFIX="$("$CONDA_BIN" env list --json | "${PAPER_PYTHON:-python3}" -c '
import json, os, sys
for prefix in json.load(sys.stdin).get("envs", []):
    if os.path.basename(prefix) == sys.argv[1]:
        print(prefix)
        break
' "$ENV_NAME")"
fi

PYTHON_BIN="${ENV_PREFIX}/bin/python"
[[ -x "$PYTHON_BIN" ]] || { printf 'Environment has no interpreter: %s\n' "$PYTHON_BIN" >&2; exit 1; }

"$PYTHON_BIN" - <<'PY'
import importlib, sys
required = ["numpy", "pandas", "scipy", "matplotlib", "PIL", "tifffile",
            "skimage", "cv2", "aicspylibczi", "pypdf"]
missing = []
for module in required:
    try:
        importlib.import_module(module)
    except Exception:
        missing.append(module)
if missing:
    print("MISSING modules:", ", ".join(missing))
    sys.exit(1)
print("all required modules import cleanly")
PY

cat <<EOF

Use this environment for every figure:

  export PAPER_PYTHON=${PYTHON_BIN}
  export PAPER_PYTHON_CZI=${PYTHON_BIN}

Then follow README.txt from the Paper root.
EOF
