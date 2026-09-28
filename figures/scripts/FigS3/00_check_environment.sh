#!/usr/bin/env bash
set -euo pipefail

FIGS3_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ "$(basename "$(dirname "$FIGS3_DIR")")" == scripts ]]; then
  FIGS3_DIR="$(cd "$FIGS3_DIR/../.." && pwd -P)/FigS3"
fi
MODE="${1:---replay}"

case "$MODE" in
  --replay)
    REPLAY_PYTHON="${PAPER_S3_PYTHON:-${PAPER_REVIEW_PYTHON:-${PAPER_PYTHON:-}}}"
    if [[ -z "$REPLAY_PYTHON" ]]; then
      for candidate in \
        /home/server/anaconda3/envs/paper_apotome_s3_repro/bin/python \
        /opt/conda/envs/paper_apotome_s3_repro/bin/python \
        /home/server/anaconda3/bin/python; do
        if [[ -x "$candidate" ]]; then
          REPLAY_PYTHON="$candidate"
          break
        fi
      done
    fi
    [[ -n "$REPLAY_PYTHON" && -x "$REPLAY_PYTHON" ]] || {
      printf '%s\n' \
        'No Figure S3 replay environment found.' \
        'Create it with: ./scripts/setup/01_create_s3_environment.sh' >&2
      exit 2
    }
    for command_name in pdfinfo pdffonts pdftoppm pdftotext; do
      command -v "$command_name" >/dev/null || {
        printf 'Missing Figure S3 replay command: %s\n' "$command_name" >&2
        exit 2
      }
    done
    test -f "$FIGS3_DIR/exports/Figure_S3_L0_L2_export_manifest.csv"
    test -f "$FIGS3_DIR/exports/Figure_S3_L0_L2_export_summary.json"
    test -f "$FIGS3_DIR/histoplus/inference_summary.json"
    test -f "$FIGS3_DIR/organ_segmentation/segmentation_summary.json"
    test ! -e "$FIGS3_DIR/registration"
    test "$(find "$FIGS3_DIR/exports" -mindepth 2 -maxdepth 2 -name '1_L0_rgb.tif' -type f | wc -l)" -eq 24
    test "$(find "$FIGS3_DIR/exports" -mindepth 2 -maxdepth 2 -name '1_L2_rgb.tif' -type f | wc -l)" -eq 24
    "$REPLAY_PYTHON" -c 'import cv2, matplotlib, numpy, pandas, pyarrow, pyvips, scipy, seaborn, shapely, sklearn, statsmodels.api, tifffile; from PIL import Image; from pypdf import PdfReader; print("Figure S3 portable replay environment: OK")'
    "$REPLAY_PYTHON" -c 'import pandas as pd, sys; frame = pd.read_parquet(sys.argv[1]); assert not frame.empty; print("Figure S3 Parquet replay: OK")' "$FIGS3_DIR/histoplus/m26-001/cells.parquet"
    printf 'Figure S3 replay Python: %s\n' "$REPLAY_PYTHON"
    printf '%s\n' 'Figure S3 portable replay check: PASS'
    ;;
  --full)
    SOURCE_ROOT="${FIGS3_SOURCE_ROOT:-/media/server/STORAGE/Motic_AnatomiaPatologica_2025/HE_Liver_Spleen_Kidney_Nancy_2026}"
    SYSTEM_PYTHON="${FIGS3_SYSTEM_PYTHON:-/home/server/anaconda3/bin/python}"
    LAZYSLIDE_PYTHON="${FIGS3_LAZYSLIDE_PYTHON:-/home/server/anaconda3/envs/lazyslide311/bin/python}"
    WEIGHT_FILE="${HISTOPLUS_WEIGHT_FILE:-/home/server/.cache/histoplus/histoplus_cellvit_segmentor_20x.pt}"
    test -d "$SOURCE_ROOT"
    test -f "$SOURCE_ROOT/EXCEL_ORGANOS.xlsx"
    test -x "$SYSTEM_PYTHON"
    test -x "$LAZYSLIDE_PYTHON"
    test -f "$WEIGHT_FILE"
    test "$(find "$SOURCE_ROOT" -mindepth 2 -maxdepth 2 -name '1.mds' -type f | wc -l)" -eq 24
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader
    "$SYSTEM_PYTHON" -c 'import cv2, olefile, openpyxl, pandas, scipy, skimage, sklearn, statsmodels, tifffile; print("system analysis environment: OK")'
    "$LAZYSLIDE_PYTHON" -c 'import lazyslide, torch, wsidata; assert torch.cuda.is_available(); print(f"LazySlide {lazyslide.__version__}; torch {torch.__version__}; CUDA device {torch.cuda.get_device_name(0)}")'
    sha256sum "$WEIGHT_FILE"
    du -sh "$SOURCE_ROOT" "$FIGS3_DIR/exports" 2>/dev/null || true
    printf '%s\n' 'Figure S3 full raw/GPU environment check: PASS'
    ;;
  -h|--help)
    printf '%s\n' \
      'Usage: FigS3/00_check_environment.sh [--replay|--full]' \
      '  --replay  Check the portable shipped-export renderer (default).' \
      '  --full    Check external MDS, HistoPLUS weight, LazySlide, and CUDA.'
    ;;
  *)
    printf 'Unknown Figure S3 environment mode: %s\n' "$MODE" >&2
    exit 2
    ;;
esac
