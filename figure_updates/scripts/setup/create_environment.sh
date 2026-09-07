#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="${PAPER_ENV:-paper_apotome_repro}"
CONDA_EXE="${CONDA_EXE:-/home/server/anaconda3/bin/conda}"
FORCE=0

usage() {
  cat <<'EOF'
Usage: ./scripts/setup/create_environment.sh [--force]

Create the single conda environment used by all Paper figure workflows.
--force removes an existing environment with the same name before recreating
it from the pinned environment.yml file.
EOF
}

while (($#)); do
  case "$1" in
    --force) FORCE=1 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

[[ -x "$CONDA_EXE" ]] || {
  printf 'Conda executable not found: %s\n' "$CONDA_EXE" >&2
  exit 1
}
[[ -s "$ENV_FILE" ]] || {
  printf 'Environment file not found: %s\n' "$ENV_FILE" >&2
  exit 1
}

if "$CONDA_EXE" env list --json | \
    "$CONDA_EXE" run -n base python -c \
      'import json,sys; n=sys.argv[1]; raise SystemExit(0 if any(p.rstrip("/").split("/")[-1] == n for p in json.load(sys.stdin)["envs"]) else 1)' \
      "$ENV_NAME"; then
  if [[ "$FORCE" != 1 ]]; then
    printf 'Environment %s already exists. Re-run with --force to recreate it.\n' "$ENV_NAME" >&2
    exit 2
  fi
  "$CONDA_EXE" env remove --name "$ENV_NAME" --yes
fi

"$CONDA_EXE" env create --name "$ENV_NAME" --file "$ENV_FILE"
PYTHONNOUSERSITE=1 "$CONDA_EXE" run --no-capture-output --name "$ENV_NAME" \
  python -c 'import aicspylibczi, cv2, flask, imagecodecs, matplotlib, numpy, openpyxl, pandas, scipy, sklearn, skimage, tifffile, torch, tqdm; print("[OK] Paper environment import check passed")'

printf 'Created environment: %s\n' "$ENV_NAME"
