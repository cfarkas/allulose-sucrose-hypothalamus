#!/usr/bin/env bash
# Reproduce Figure 2 end to end from Fig2/raw_data.
#
#   Fig2/07_run_all.sh /tmp/fig2_stage /tmp/fig2_render
#
# Both paths must be absolute, under /tmp, and must not already exist.
set -euo pipefail

if [[ $# -ne 2 ]]; then
  printf 'Usage: %s /absolute/fresh/stage /absolute/fresh/render\n' "$0" >&2
  exit 2
fi

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
STAGE="$1"
RENDER="$2"
for path in "$STAGE" "$RENDER"; do
  case "$path" in
    /tmp/*) ;;
    *) printf 'Path must be absolute and under /tmp: %s\n' "$path" >&2; exit 2 ;;
  esac
  if [[ -e "$path" || -L "$path" ]]; then
    printf 'Fresh-only path already exists: %s\n' "$path" >&2
    exit 1
  fi
done

# Override on another machine: PAPER_PYTHON and PAPER_PYTHON_CZI.
# One environment, paper_apotome_repro, satisfies every step. Run
# Fig2/00_create_conda_envs.sh first; it prints the two exports.
DEFAULT_PY=/home/server/anaconda3/envs/paper_apotome_repro/bin/python
PY_BASE="${PAPER_PYTHON:-$DEFAULT_PY}"
PY_CZI="${PAPER_PYTHON_CZI:-$PY_BASE}"
for interpreter in "$PY_BASE" "$PY_CZI"; do
  command -v "$interpreter" >/dev/null 2>&1 || [[ -x "$interpreter" ]] || {
    printf 'Interpreter not found: %s\nSet PAPER_PYTHON and PAPER_PYTHON_CZI.\n' "$interpreter" >&2
    exit 2
  }
done
export PYTHONDONTWRITEBYTECODE=1

# A requested reproduction is a foreground task by default. The previous
# unconditional idle I/O and nice-19 wrapper could starve CZI reads for hours
# whenever another process was using STORAGE. Set PAPER_BACKGROUND_PRIORITY=1
# explicitly when a low-impact background run is preferred.
case "${PAPER_BACKGROUND_PRIORITY:-0}" in
  0) RUN=() ;;
  1) RUN=(ionice -c3 nice -n19) ;;
  *)
    printf 'PAPER_BACKGROUND_PRIORITY must be 0 or 1, found: %s\n' \
      "$PAPER_BACKGROUND_PRIORITY" >&2
    exit 2
    ;;
esac

run_step() {
  local label="$1"
  local started elapsed
  shift
  started=$SECONDS
  printf '[Fig2] START %s\n' "$label"
  "${RUN[@]}" "$@"
  elapsed=$((SECONDS - started))
  printf '[Fig2] DONE  %s (%dm %02ds)\n' \
    "$label" "$((elapsed / 60))" "$((elapsed % 60))"
}

# 00  validate the figure-local raw bundle (high I/O)
run_step 'validate raw bundle' \
  "$PY_BASE" "${HERE}/00_validate_raw_bundle.py"
# 02  animal-level analysis stage (01 is the analysis engine it imports)
run_step 'stage animal-level analysis' \
  "$PY_BASE" "${HERE}/02_stage_analysis.py" "$STAGE"
# 04  panel C native RAW-CZI display assets (03 is the CZI engine it imports)
run_step 'reconstruct panel C native CZI assets' \
  "$PY_CZI" "${HERE}/04_stage_microscopy.py" "${STAGE}/native_czi"
# 05  panel B accepted HIL anatomy assets
run_step 'stage panel B HIL anatomy assets' \
  "$PY_CZI" "${HERE}/05_stage_human_regions.py" "${STAGE}/human_regions"
# 06  master figure plus isolated English/Spanish subpanels
run_step 'render master and bilingual panels' \
  "$PY_BASE" "${HERE}/06_make_figure_2.py" \
  --stage-root "$STAGE" --outdir "$RENDER"

printf '[OK] Figure 2 rebuilt: %s/outputs\n' "$RENDER"
