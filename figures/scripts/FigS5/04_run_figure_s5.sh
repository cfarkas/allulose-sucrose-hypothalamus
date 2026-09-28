#!/usr/bin/env bash
set -euo pipefail

FIGS5_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER_DIR="$FIGS5_DIR"
while [[ ! -d "$PAPER_DIR/scripts/setup" || ! -d "$PAPER_DIR/FigS5" ]]; do
  [[ "$PAPER_DIR" != / ]] || { printf 'Paper root not found\n' >&2; exit 2; }
  PAPER_DIR="$(dirname "$PAPER_DIR")"
done
SOURCE_DIR="$PAPER_DIR/FigS5"
PY="${PAPER_PYTHON:-/home/server/anaconda3/envs/paper_apotome_repro/bin/python}"
WORK="${FIGS5_WORKDIR:-$PAPER_DIR/analyses/FigS5/results/including_NPY_M_20260908}"
EXPORT_DIR="$SOURCE_DIR/channel_tiffs"
PIPELINE="$SOURCE_DIR/02a_july_cfos_acth_clip_human_in_loop_s5.py"
RENDERER="$FIGS5_DIR/03_make_figure_s5_acth_clip_cfos.py"
IMPORTER="$PAPER_DIR/scripts/utilities/06_import_figure_s5_hil.py"
SPATIAL_REVIEWER="$SOURCE_DIR/05_review_spatial_tissue_hil.py"
SPATIAL_HIL_DIR="$SOURCE_DIR/hil_review/spatial_3v_including_NPY_M_20260908"
ANATOMY_HIL_DIR="$SOURCE_DIR/hil_review/human_final_including_NPY_M_20260908"
OVERLAY_DIR="$SOURCE_DIR/provenance/inclusion_NPY_M_20260908"
ACTION="${1:-reproduce}"

if [[ ! -x "$PY" ]]; then
  printf 'Python is not executable: %s\nSet PAPER_PYTHON to the paper environment.\n' "$PY" >&2
  exit 2
fi

case "$ACTION" in
  prepare)
    "$PY" "$SOURCE_DIR/06_include_npy_m.py" --root "$PAPER_DIR"
    "$PY" "$PIPELINE" --self-test
    "$PY" "$PIPELINE" \
      --prepare \
      --sample-manifest "$OVERLAY_DIR/sample_manifest.csv" \
      --channel-map "$OVERLAY_DIR/channel_map.csv" \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    "$PY" "$PIPELINE" \
      --audit-segmentations \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    ;;
  review)
    PORT="${2:-34249}"
    "$PY" "$PIPELINE" \
      --annotate \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK" \
      --host 0.0.0.0 \
      --port "$PORT"
    ;;
  spatial-review)
    PORT="${2:-34250}"
    "$PY" "$SPATIAL_REVIEWER" \
      --workdir "$WORK" \
      --review-dir "$SPATIAL_HIL_DIR" \
      --host 0.0.0.0 \
      --port "$PORT"
    ;;
  spatial-status)
    "$PY" "$SPATIAL_REVIEWER" \
      --workdir "$WORK" \
      --review-dir "$SPATIAL_HIL_DIR" \
      --status
    ;;
  status)
    "$PY" "$PIPELINE" \
      --status \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    ;;
  render)
    "$PY" "$RENDERER" --workdir "$WORK" --spatial-hil-dir "$SPATIAL_HIL_DIR" \
      --output-dir "${FIGS5_OUTPUT_DIR:?Set an absent FIGS5_OUTPUT_DIR}" --dpi 600
    ;;
  build)
    if [[ ! -f "$SPATIAL_HIL_DIR/SPATIAL_TISSUE_HIL_RECEIPT.json" ]]; then
      printf 'Finalized DAPI-only 3V exclusion HIL for I/J is missing: %s\n' "$SPATIAL_HIL_DIR" >&2
      printf 'Complete: %s spatial-review 34250\n' "$0" >&2
      exit 7
    fi
    RESULT_SENTINEL="$WORK/results/per_analysis_unit_region_summary.csv"
    if [[ -e "$RESULT_SENTINEL" ]]; then
      printf 'Refusing to replace an existing Figure S5 result set: %s\n' "$RESULT_SENTINEL" >&2
      printf 'Set FIGS5_WORKDIR to a fresh directory and repeat prepare/review/build.\n' >&2
      exit 3
    fi
    "$PY" "$PIPELINE" \
      --audit-segmentations \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    "$PY" "$PIPELINE" \
      --quantify \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK" \
      --qc-format both \
      --qc-dpi 600 \
      --plot-dpi 600
    if [[ -n "${FIGS5_OUTPUT_DIR:-}" ]]; then
      OUT="$FIGS5_OUTPUT_DIR"
    else
      BUILD_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
      OUT="$PAPER_DIR/analyses/FigS5/figure/acth_clip_candidate_$BUILD_STAMP"
    fi
    if [[ -e "$OUT" ]]; then
      printf 'Choose an absent FIGS5_OUTPUT_DIR; refusing to replace %s\n' "$OUT" >&2
      exit 4
    fi
    "$PY" "$RENDERER" \
      --workdir "$WORK" \
      --output-dir "$OUT" \
      --spatial-hil-dir "$SPATIAL_HIL_DIR" \
      --dpi 600
    printf 'Figure S5 candidate complete: %s\n' "$OUT"
    ;;
  reproduce)
    if [[ ! -f "$SPATIAL_HIL_DIR/SPATIAL_TISSUE_HIL_RECEIPT.json" ]]; then
      printf 'Finalized DAPI-only 3V exclusion HIL for I/J is missing: %s\n' "$SPATIAL_HIL_DIR" >&2
      printf 'Complete: %s spatial-review 34250\n' "$0" >&2
      exit 7
    fi
    "$PY" "$RENDERER" \
      --spatial-hil-dir "$SPATIAL_HIL_DIR" \
      --validate-spatial-hil-source-only
    if [[ -n "${FIGS5_WORKDIR:-}" ]]; then
      if [[ -e "$WORK" || -L "$WORK" ]]; then
        printf 'Choose an absent FIGS5_WORKDIR; refusing to reuse %s\n' "$WORK" >&2
        exit 5
      fi
      mkdir -p -m 700 -- "$WORK"
    else
      WORK="$(mktemp -d /tmp/figs5_rebuild_XXXXXXXX)"
    fi
    "$PY" "$SOURCE_DIR/06_include_npy_m.py" --root "$PAPER_DIR"
    "$PY" "$PIPELINE" --self-test
    "$PY" "$PIPELINE" \
      --prepare \
      --sample-manifest "$OVERLAY_DIR/sample_manifest.csv" \
      --channel-map "$OVERLAY_DIR/channel_map.csv" \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    "$PY" "$PIPELINE" \
      --audit-segmentations \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK"
    "$PY" "$IMPORTER" --from "$ANATOMY_HIL_DIR" --workdir "$WORK"
    "$PY" "$PIPELINE" \
      --quantify \
      --export-dir "$EXPORT_DIR" \
      --workdir "$WORK" \
      --qc-format both \
      --qc-dpi 600 \
      --plot-dpi 600
    if [[ -n "${FIGS5_OUTPUT_DIR:-}" ]]; then
      OUT="$FIGS5_OUTPUT_DIR"
    else
      OUT="$WORK/figure_candidate"
    fi
    if [[ -e "$OUT" || -L "$OUT" ]]; then
      printf 'Choose an absent FIGS5_OUTPUT_DIR; refusing to replace %s\n' "$OUT" >&2
      exit 6
    fi
    "$PY" "$RENDERER" \
      --workdir "$WORK" \
      --output-dir "$OUT" \
      --spatial-hil-dir "$SPATIAL_HIL_DIR" \
      --dpi 600
    printf 'Figure S5 reproduced from accepted HIL.\nWork: %s\nCandidate: %s\n' \
      "$WORK" "$OUT"
    ;;
  *)
    printf 'Usage: %s {prepare|review [port]|spatial-review [port]|spatial-status|status|build|render|reproduce}\n' "$0" >&2
    exit 2
    ;;
esac
