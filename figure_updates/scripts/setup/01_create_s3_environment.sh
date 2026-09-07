#!/usr/bin/env bash
# Create the separately pinned Figure S3 replay/review environment.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PAPER_ROOT="$(cd -- "$SCRIPT_DIR/../.." && pwd -P)"
ENV_FILE="$SCRIPT_DIR/environment_s3.yml"
CONDA_LOCK="$SCRIPT_DIR/environment_s3_conda_explicit_linux-64.txt"
PIP_LOCK="$SCRIPT_DIR/requirements_s3_pip_lock.txt"
ENV_NAME="${PAPER_S3_ENV:-paper_apotome_s3_repro}"
FORCE=0
CHECK_ONLY=0

while (($#)); do
  case "$1" in
    --force) FORCE=1 ;;
    --check) CHECK_ONLY=1 ;;
    -h|--help)
      printf '%s\n' 'Usage: scripts/setup/01_create_s3_environment.sh [--check|--force]'
      exit 0
      ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; exit 2 ;;
  esac
  shift
done

[[ -f "$ENV_FILE" && -d "$PAPER_ROOT/FigS3" ]] || {
  printf 'Cannot locate Figure S3 environment inputs below %s\n' "$PAPER_ROOT" >&2
  exit 2
}

CONDA_BIN=""
if [[ -n "${CONDA_EXE:-}" && -x "$CONDA_EXE" ]]; then
  CONDA_BIN="$CONDA_EXE"
elif command -v conda >/dev/null 2>&1; then
  CONDA_BIN="$(command -v conda)"
else
  for candidate in \
    /opt/conda/bin/conda /usr/local/bin/conda \
    /home/server/anaconda3/bin/conda; do
    if [[ -x "$candidate" ]]; then
      CONDA_BIN="$candidate"
      break
    fi
  done
fi
[[ -n "$CONDA_BIN" ]] || {
  printf '%s\n' 'No conda found; set CONDA_EXE to Miniforge/Miniconda conda.' >&2
  exit 2
}

JSON_PYTHON="${PAPER_PYTHON:-python3}"
ENV_PREFIX="$("$CONDA_BIN" env list --json 2>/dev/null | "$JSON_PYTHON" -c '
import json, os, sys
name = sys.argv[1]
for prefix in json.load(sys.stdin).get("envs", []):
    if os.path.basename(prefix) == name:
        print(prefix)
        break
' "$ENV_NAME" || true)"

if [[ -n "$ENV_PREFIX" && "$FORCE" -eq 0 ]]; then
  printf 'Figure S3 environment already exists: %s\n' "$ENV_PREFIX"
elif [[ "$CHECK_ONLY" -eq 1 ]]; then
  printf 'Figure S3 environment is absent: %s\n' "$ENV_NAME" >&2
  exit 1
else
  if [[ -n "$ENV_PREFIX" ]]; then
    "$CONDA_BIN" env remove --name "$ENV_NAME" --yes
  fi
  if [[ "$(uname -s)" == "Linux" && "$(uname -m)" == "x86_64"         && -f "$CONDA_LOCK" && -f "$PIP_LOCK" ]]; then
    "$CONDA_BIN" create --name "$ENV_NAME" --file "$CONDA_LOCK" --yes
    "$CONDA_BIN" run --name "$ENV_NAME"       python -m pip install --no-deps --requirement "$PIP_LOCK"
  else
    printf '%s\n'       'No platform lock is provided for this host; using the pinned cross-platform YAML.'
    "$CONDA_BIN" env create --name "$ENV_NAME" --file "$ENV_FILE"
  fi
  ENV_PREFIX="$("$CONDA_BIN" env list --json | "$JSON_PYTHON" -c '
import json, os, sys
name = sys.argv[1]
for prefix in json.load(sys.stdin).get("envs", []):
    if os.path.basename(prefix) == name:
        print(prefix)
        break
' "$ENV_NAME")"
fi

PYTHON_BIN="$ENV_PREFIX/bin/python"
[[ -x "$PYTHON_BIN" ]] || { printf 'Missing interpreter: %s\n' "$PYTHON_BIN" >&2; exit 1; }
PATH="$ENV_PREFIX/bin:$PATH" PAPER_S3_PYTHON="$PYTHON_BIN" \
  "$PAPER_ROOT/FigS3/00_check_environment.sh" --replay

cat <<EOF

Use the main and Figure S3 environments for the complete rebuild:

  export PAPER_PYTHON=<paper_apotome_repro>/bin/python
  export PAPER_PYTHON_CZI=\$PAPER_PYTHON
  export PAPER_S3_PYTHON=${PYTHON_BIN}
  export PAPER_REVIEW_PYTHON=${PYTHON_BIN}

Then run ./reproduce_all_figures.sh --check-only.
The root launcher scopes the S3 environment to the commands that require it;
do not prepend the S3 environment to PATH globally.
EOF
