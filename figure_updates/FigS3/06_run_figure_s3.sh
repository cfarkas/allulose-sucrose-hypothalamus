#!/usr/bin/env bash
set -euo pipefail

FIGS3_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PAPER_DIR="$(cd "$FIGS3_DIR/.." && pwd)"
SOURCE_ROOT="${FIGS3_SOURCE_ROOT:-/media/server/STORAGE/Motic_AnatomiaPatologica_2025/HE_Liver_Spleen_Kidney_Nancy_2026}"
SYSTEM_PYTHON="${FIGS3_SYSTEM_PYTHON:-/home/server/anaconda3/bin/python}"
LAZYSLIDE_PYTHON="${FIGS3_LAZYSLIDE_PYTHON:-/home/server/anaconda3/envs/lazyslide311/bin/python}"
WEIGHT_FILE="${HISTOPLUS_WEIGHT_FILE:-/home/server/.cache/histoplus/histoplus_cellvit_segmentor_20x.pt}"
EXPORT_WORKERS="${FIGS3_EXPORT_WORKERS:-2}"
GPU_BATCH_SIZE="${FIGS3_GPU_BATCH_SIZE:-6}"
GPU_DATA_WORKERS="${FIGS3_GPU_DATA_WORKERS:-6}"

cd "$PAPER_DIR"
mkdir -p "$FIGS3_DIR/logs"
RUN_LOG="$FIGS3_DIR/logs/run_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "$RUN_LOG") 2>&1

echo "Figure S3 run started: $(date -u --iso-8601=seconds)"
"$FIGS3_DIR/00_check_environment.sh"
"$SYSTEM_PYTHON" "$FIGS3_DIR/01_prepare_metadata.py" --source-root "$SOURCE_ROOT"
"$SYSTEM_PYTHON" "$FIGS3_DIR/02_export_mds_levels.py" --workers "$EXPORT_WORKERS" --resume
"$SYSTEM_PYTHON" "$FIGS3_DIR/03_segment_organs.py"
"$LAZYSLIDE_PYTHON" "$FIGS3_DIR/04_run_histoplus_gpu.py" \
  --weight-file "$WEIGHT_FILE" \
  --tiles-per-organ 24 \
  --batch-size "$GPU_BATCH_SIZE" \
  --num-workers "$GPU_DATA_WORKERS" \
  --resume
"$SYSTEM_PYTHON" "$FIGS3_DIR/05_analyze_make_figure.py" --permutations 99999 --dpi 600
"$SYSTEM_PYTHON" "$FIGS3_DIR/07_validate_figure_s3.py"
echo "Figure S3 run completed: $(date -u --iso-8601=seconds)"
