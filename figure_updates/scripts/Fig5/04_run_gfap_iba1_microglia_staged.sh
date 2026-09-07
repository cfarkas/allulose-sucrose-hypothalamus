#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Locate the Paper root by walking up to the directory that holds scripts/shared,
# so this works from the figure folder or any relocated copy of the tree.
PAPER_ROOT="$(cd -- "$SCRIPT_DIR" && while [[ ! -d scripts/shared && "$PWD" != / ]]; do cd ..; done; pwd -P)"
[[ -d "$PAPER_ROOT/scripts/shared" ]] || { printf 'Cannot locate the Paper root above %s\n' "$SCRIPT_DIR" >&2; exit 2; }
DRIVER="${PAPER_ROOT}/scripts/shared/03_safe_stage_fig45.py"
PYTHON_BIN="${PAPER_PYTHON:-/home/server/anaconda3/envs/paper_apotome_repro/bin/python}"

[[ -f "${DRIVER}" ]] || {
  printf 'Missing safe Figure 4/5 staging driver: %s\n' "${DRIVER}" >&2
  exit 2
}
[[ -x "${PYTHON_BIN}" ]] || {
  printf 'Missing executable Python environment: %s\n' "${PYTHON_BIN}" >&2
  exit 2
}

exec "${PYTHON_BIN}" "${DRIVER}" --figure fig5 "$@"
