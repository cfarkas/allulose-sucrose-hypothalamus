#!/usr/bin/env bash
# Rebuild every established paper figure into the canonical Paper tree.
#
# Simplest use from the Paper root:
#
#   ./reproduce_all_figures.sh
#
# Every figure is first rendered and validated in isolated /tmp storage, then
# installed into the existing FigN/ and analyses/ layout only after all thirteen
# packages pass. Raw data and active scripts remain immutable package inputs.
# Supplementary Figure S5 reuses the accepted DAPI-only anatomy HIL and
# hash-audits the separate 3V-exclusion review.

# Stable Python hashing prevents set/dict iteration from changing raster output.
export PYTHONHASHSEED=0
export PYTHONUTF8=1
# Matplotlib otherwise embeds the wall-clock save time in every PDF. A fixed
# source epoch makes independent reproductions byte-identical as well as
# visually identical.
export SOURCE_DATE_EPOCH=0
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PAPER_ROOT="$SCRIPT_DIR"
FIGS5_SPATIAL_HIL="$PAPER_ROOT/FigS5/hil_review/spatial_3v_including_NPY_M_20260908"
ENV_NAME="${PAPER_ENV:-paper_apotome_repro}"
OUT_ARG=""
PROMOTE_COMPAT=0
FORCE_FIG3=0
CHECK_ONLY=0
PLAN_ONLY=0
DO_MICROGLIAL_CHOICES=0
BACKGROUND_PRIORITY="${PAPER_BACKGROUND_PRIORITY:-0}"
ACTIVE_INTERNAL_STAGE=""
BUILD_ROOT=""
RAW_INPUT_SNAPSHOT=""
SCRIPT_INPUT_SNAPSHOT=""
INSTALL_RECOVERY_ROOT=""

usage() {
  cat <<'EOF'
Usage: ./reproduce_all_figures.sh [options]

Rebuild Figures 1-5 plus S1-S8, including S3 from its completed WSI receipts
and S5 from accepted HIL anatomy. All figure graphics are rendered at 600 dpi.

Options:
  --check-only       Validate dependencies, workflows, and accepted HIL inputs.
  --plan             Print the steps without creating output or running them.
  --do_microglial_choices
                     Optional: open a blank 300-cell microglial review and
                     retrain from the new choices. Default reuses the current
                     203 choices and saved candidate without a review prompt.
  --background       Run Figure 2 at idle CPU/I/O priority instead of the
                     default foreground/full-speed priority.
  --output PATH      Canonical output root. PATH must be this Paper directory.
  --out PATH         Backward-compatible alias for --output.
  --promote          Accepted for compatibility; canonical installation is now
                     automatic after all thirteen staged builds pass.
  --force            Also install the validated certified Figure 3 rebuild.
                     Without it, Figure 3 is verified but remains unchanged.
  -h, --help         Show this help.

The canonical output is this existing Paper tree. Every analysis and render is
first produced in fresh /tmp staging. Only after all thirteen figures pass validation
are current generated analyses and publication artifacts installed into the
canonical analyses/ and FigN/ paths. Active scripts, raw, raw_data, accepted HIL
inputs and figure READMEs remain in place and are audited before and after the
run. No top-level apotome_rebuild_* output is created. Figures 2 and 3 retain
their additional mandatory isolated computation stages. Certified Figure 3 is
installed only when --force is passed.
EOF
}

while (($#)); do
  case "$1" in
    --check-only)
      CHECK_ONLY=1
      ;;
    --plan)
      PLAN_ONLY=1
      ;;
    --do_microglial_choices|--do-microglial-choices)
      DO_MICROGLIAL_CHOICES=1
      ;;
    --background)
      BACKGROUND_PRIORITY=1
      ;;
    --promote)
      PROMOTE_COMPAT=1
      ;;
    --force)
      FORCE_FIG3=1
      ;;
    --output|--out)
      shift
      (($#)) || { printf '%s\n' 'Missing path after --output/--out.' >&2; usage >&2; exit 2; }
      OUT_ARG="$1"
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

case "$BACKGROUND_PRIORITY" in
  0|1) ;;
  *)
    printf 'PAPER_BACKGROUND_PRIORITY must be 0 or 1, found: %s\n' \
      "$BACKGROUND_PRIORITY" >&2
    exit 2
    ;;
esac
if ((BACKGROUND_PRIORITY)); then
  PRIORITY_LABEL='background (idle I/O, nice 19)'
else
  PRIORITY_LABEL='foreground (normal CPU/I/O)'
fi

if ((CHECK_ONLY && PLAN_ONLY)); then
  printf '%s\n' '--check-only and --plan cannot be combined.' >&2
  exit 2
fi

validate_output_argument() {
  local resolved
  if [[ -z "$OUT_ARG" ]]; then
    OUT_ARG="$PAPER_ROOT"
  fi
  if [[ "$OUT_ARG" != /* ]]; then
    printf 'The canonical --output path must be absolute: %s\n' "$OUT_ARG" >&2
    return 2
  fi
  if ! resolved="$(readlink -f -- "$OUT_ARG")"; then
    printf 'Cannot resolve --output path: %s\n' "$OUT_ARG" >&2
    return 2
  fi
  if [[ "$resolved" != "$PAPER_ROOT" ]]; then
    printf 'The complete canonical rebuild must output directly to %s, not %s\n' \
      "$PAPER_ROOT" "$OUT_ARG" >&2
    return 2
  fi
  OUT_ARG="$resolved"
}
validate_output_argument

if ((PROMOTE_COMPAT)); then
  printf '%s\n' 'Note: --promote is now implicit for a validated canonical rebuild.'
fi


find_conda_base() {
  local conda_bin=""
  if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    conda_bin="$CONDA_EXE"
  elif command -v conda >/dev/null 2>&1; then
    conda_bin="$(command -v conda)"
  elif [[ -x /home/server/anaconda3/bin/conda ]]; then
    conda_bin=/home/server/anaconda3/bin/conda
  fi
  if [[ -n "$conda_bin" ]]; then
    "$conda_bin" info --base 2>/dev/null || true
  fi
}

select_scientific_python() {
  local conda_base="$1"
  local candidate
  local -a candidates=()
  [[ -n "${PAPER_PYTHON:-}" ]] && candidates+=("$PAPER_PYTHON")
  [[ -n "$conda_base" ]] && candidates+=("$conda_base/envs/$ENV_NAME/bin/python")
  candidates+=(
    "/home/server/anaconda3/envs/$ENV_NAME/bin/python"
    "/opt/conda/envs/$ENV_NAME/bin/python"
  )
  for candidate in "${candidates[@]}"; do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

can_review_fig3() {
  local candidate="$1"
  [[ -x "$candidate" ]] || return 1
  "$candidate" -c 'import cv2, matplotlib, numpy, pandas, pyarrow, pyvips, scipy, seaborn, shapely, sklearn, statsmodels.api, tifffile; from PIL import Image; from pypdf import PdfReader' \
    >/dev/null 2>&1
}

select_review_python() {
  local scientific="$1"
  local conda_base="$2"
  local candidate
  local -a candidates=()
  [[ -n "${PAPER_S3_PYTHON:-}" ]] && candidates+=("$PAPER_S3_PYTHON")
  [[ -n "${PAPER_REVIEW_PYTHON:-}" ]] && candidates+=("$PAPER_REVIEW_PYTHON")
  [[ -n "$conda_base" ]] && candidates+=("$conda_base/envs/paper_apotome_s3_repro/bin/python")
  candidates+=(
    /home/server/anaconda3/envs/paper_apotome_s3_repro/bin/python
    /opt/conda/envs/paper_apotome_s3_repro/bin/python
    "$scientific"
  )
  [[ -n "$conda_base" ]] && candidates+=("$conda_base/bin/python")
  candidates+=(/home/server/anaconda3/bin/python)
  for candidate in "${candidates[@]}"; do
    if can_review_fig3 "$candidate"; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done
  return 1
}

CONDA_BASE="$(find_conda_base)"
if ! PY="$(select_scientific_python "$CONDA_BASE")"; then
  printf '%s\n' \
    "Cannot find the $ENV_NAME Python environment." \
    'Create it once with: ./Fig1/00_create_conda_envs.sh' \
    'Or set PAPER_PYTHON=/absolute/path/to/python.' >&2
  exit 2
fi
PY_CZI="${PAPER_PYTHON_CZI:-$PY}"
if [[ ! -x "$PY_CZI" ]]; then
  printf 'PAPER_PYTHON_CZI is not executable: %s\n' "$PY_CZI" >&2
  exit 2
fi
if ! REVIEW_PY="$(select_review_python "$PY" "$CONDA_BASE")"; then
  printf '%s\n' \
    'No Python interpreter can import the Figure 3 review dependencies.' \
    'Recreate the paper environment with ./Fig1/00_create_conda_envs.sh.' >&2
  exit 2
fi

# The complete Figure 3 compositor uses PyMuPDF; it may live in a separate
# existing interpreter from the pinned image-review environment.
FIG3_PDF_BIN=""
for candidate in "${FIG3_PDF_PYTHON:-}" "$REVIEW_PY" "$PY" "$CONDA_BASE/envs/lazyslide311/bin/python" "$CONDA_BASE/bin/python"; do
  if [[ -n "$candidate" && -x "$candidate" ]] && "$candidate" -c 'import fitz, matplotlib, numpy, pandas, scipy, skimage; from PIL import Image' >/dev/null 2>&1; then
    FIG3_PDF_BIN="$candidate"
    break
  fi
done
if [[ -z "$FIG3_PDF_BIN" ]]; then
  printf '%s\n' 'Figure 3 completion requires PyMuPDF, Pillow, Matplotlib, NumPy, pandas, SciPy and scikit-image. Set FIG3_PDF_PYTHON to an interpreter containing these packages.' >&2
  exit 2
fi

export PAPER="$PAPER_ROOT"
export FIG2_PAPER_ROOT="$PAPER_ROOT"
export FIG4_PAPER_ROOT="$PAPER_ROOT"
export PAPER_PYTHON="$PY"
export PAPER_PYTHON_CZI="$PY_CZI"
export PAPER_REVIEW_PYTHON="$REVIEW_PY"
export PAPER_S3_PYTHON="$REVIEW_PY"
export PAPER_BACKGROUND_PRIORITY="$BACKGROUND_PRIORITY"
export PYTHONDONTWRITEBYTECODE=1
REVIEW_BIN="$(dirname -- "$REVIEW_PY")"
# Scope the pinned review/Poppler bin directory to commands that need it.
# Exporting it globally can make `conda run -n ... python` resolve the review
# interpreter instead of the requested environment.

printf 'Paper root       : %s\n' "$PAPER_ROOT"
printf 'Scientific Python: %s\n' "$PY"
printf 'CZI Python       : %s\n' "$PY_CZI"
printf 'Review Python    : %s\n' "$REVIEW_PY"
printf 'Resource priority: %s\n' "$PRIORITY_LABEL"

if ((PLAN_ONLY)); then
  if ((DO_MICROGLIAL_CHOICES)); then
    printf '%s\n' 'Microglial choices: independent blank review (optional flag enabled). No reviewer is opened by --plan.'
  else
    printf '%s\n' 'Microglial choices: reuse current 203 saved choices; no review prompt.'
  fi
  cat <<'EOF'

Plan (no commands were run):
  0. Validate the repository, completed S3 receipts and accepted Figure 3/S5 HIL inputs.
  1. Rebuild Figure 1 at 600 dpi.
  2. Rebuild Figure 2 at 600 dpi.
  3. Rebuild the certified English Figure 3 master plus bilingual isolated panels at 600 dpi.
  4. Rebuild the English Figure 4 master plus bilingual isolated/spatial panels at 600 dpi.
  5. Rebuild the English S5 master plus bilingual isolated panels at 600 dpi; accepted per-acquisition 3V masks exclude nuclei/density while outer tissue geometry remains anatomy-HIL-only.
  6. Load the current 203 microglial choices (or review 300 cells with --do_microglial_choices), retain the candidate results, and rebuild Figure 5 at 600 dpi.
  7. Rebuild Figure S1 at 600 dpi.
  8. Rebuild Figure S2 at 600 dpi.
  9. Recompute S3 statistics and rebuild its bilingual A-F panels from the
     completed, fully validated WSI/HistoPLUS receipts.
 10. Rebuild Figure S4 at 600 dpi.
 11. Recompute and render Figure S6 from frozen human-reference/prediction tables at 600 dpi; no classifier retraining.
 12. Render Figure S7 from verified training crops, masks and recorded losses.
 13. Render Figure S8 from the pilot-based Monte Carlo power curves.
     Every rebuild step ends by validating its master, bilingual panel pairs and bilingual figure legends.

Run ./reproduce_all_figures.sh to execute this plan.
EOF
  printf 'Canonical output: %s\n' "$PAPER_ROOT"
  printf '%s\n' 'All computation is staged under /tmp; no nested rebuild is delivered.'
  if ((FORCE_FIG3)); then
    printf '%s\n' 'Certified Figure 3 will be installed after frozen-reference validation (--force).'
  else
    printf '%s\n' 'Certified Figure 3 will be verified but its promoted bytes will remain unchanged.'
  fi
  exit 0
fi

cd "$PAPER_ROOT"

printf '%s\n' '' '== Step 0: validate repository and dependencies =='
"$PY" -c 'import cv2, matplotlib, numpy, pandas, scipy, skimage, tifffile, tqdm; from PIL import Image'
"$PY" scripts/utilities/03_check_reproducibility.py
"$PY" scripts/utilities/01_validate_bundle.py
"$PY" Fig2/00_validate_raw_bundle.py
PATH="$REVIEW_BIN:$PATH" FigS3/00_check_environment.sh --replay
PATH="$REVIEW_BIN:$PATH" "$REVIEW_PY" FigS3/07_validate_figure_s3.py --check-only
"$REVIEW_PY" Fig3/02a_review_cfos_npy_ventricles.py --prepare-only
"$PY" scripts/utilities/06_import_figure_s5_hil.py --from FigS5/hil_review/human_final_including_NPY_M_20260908 --source-only
"$PY" FigS5/03_make_figure_s5_acth_clip_cfos.py \
  --spatial-hil-dir "$FIGS5_SPATIAL_HIL" \
  --validate-spatial-hil-source-only

"$PY" FigS6/01_make_figure_s6_microglia_classifier.py --validate-only
"$PY" Fig5/12_microglial_choices.py --validate-only

if ((CHECK_ONLY)); then
  printf '%s\n' '' '[PASS] Ready to reproduce all figures.'
  exit 0
fi

umask 077
BUILD_ROOT="$(mktemp -d /tmp/apotome_full_rebuild_XXXXXXXX)"
OUT="$BUILD_ROOT"
export OUT
RAW_INPUT_SNAPSHOT="$BUILD_ROOT/canonical_inputs_before.json"
INSTALL_RECOVERY_ROOT="/tmp/apotome_canonical_superseded_$(date -u +%Y%m%dT%H%M%SZ)"
printf 'Temporary validation root: %s\n' "$BUILD_ROOT"
printf 'Canonical final root     : %s\n' "$PAPER_ROOT"
printf '%s\n' 'Raw data and scripts are immutable and must match the pre-run snapshot.'
"$PY" scripts/utilities/09_audit_canonical_inputs.py \
  --root "$PAPER_ROOT" --snapshot "$RAW_INPUT_SNAPSHOT"

start_internal_stage() {
  local figure_slug="$1"
  if [[ -n "$ACTIVE_INTERNAL_STAGE" ]]; then
    printf 'Internal stage already active: %s\n' "$ACTIVE_INTERNAL_STAGE" >&2
    return 2
  fi
  case "$figure_slug" in
    fig2|fig3) ;;
    *) printf 'Invalid internal stage label: %s\n' "$figure_slug" >&2; return 2 ;;
  esac
  ACTIVE_INTERNAL_STAGE="$(mktemp -d "/tmp/apotome_${figure_slug}_internal_XXXXXXXX")"
  printf 'Isolated %s computation stage: %s\n' "$figure_slug" "$ACTIVE_INTERNAL_STAGE"
}

remove_internal_stage() {
  local completed_stage="$ACTIVE_INTERNAL_STAGE"
  [[ -n "$completed_stage" ]] || return 0
  case "$completed_stage" in
    /tmp/apotome_fig2_internal_*|/tmp/apotome_fig3_internal_*) ;;
    *) printf 'Refusing to remove unexpected internal stage: %s\n' "$completed_stage" >&2; return 2 ;;
  esac
  [[ -d "$completed_stage" && ! -L "$completed_stage" ]] || {
    printf 'Internal stage is missing or linked: %s\n' "$completed_stage" >&2
    return 2
  }
  rm -rf -- "$completed_stage"
  ACTIVE_INTERNAL_STAGE=""
}

promote_canonical() {
  local figure="$1"
  local source="$2"
  shift 2
  "$PY" scripts/utilities/04_promote_figure.py \
    --figure "$figure" --from "$source" --replace-current "$@"
}
verify_canonical_bytes() {
  local figure="$1"
  local source="$2"
  "$PY" scripts/utilities/04_promote_figure.py \
    --figure "$figure" --from "$source" --verify-only
}
verify_reference_images() {
  local figure="$1"
  local source="$2"
  "$PY" scripts/utilities/04_promote_figure.py \
    --figure "$figure" --from "$source" --verify-reference-images
}
validate_render() {
  local figure="$1"
  local source="$2"
  PATH="$REVIEW_BIN:$PATH" "$PY" scripts/utilities/07_validate_figure_outputs.py \
    --figure "$figure" --root "$source"
}

replace_canonical_analysis() {
  local source="$1"
  local target="$2"
  local relative_target
  local backup

  case "$source" in
    "$BUILD_ROOT"/*) ;;
    *) printf 'Refusing non-stage analysis source: %s\n' "$source" >&2; return 2 ;;
  esac
  case "$target" in
    "$PAPER_ROOT/analyses/"*) ;;
    *) printf 'Refusing non-canonical analysis target: %s\n' "$target" >&2; return 2 ;;
  esac
  [[ -d "$source" && ! -L "$source" ]] || {
    printf 'Validated analysis source is missing or linked: %s\n' "$source" >&2
    return 2
  }
  relative_target="${target#"$PAPER_ROOT/"}"
  backup="$INSTALL_RECOVERY_ROOT/$relative_target"
  mkdir -p -- "$(dirname -- "$target")"
  if [[ -e "$target" || -L "$target" ]]; then
    mkdir -p -- "$(dirname -- "$backup")"
    mv -- "$target" "$backup"
  fi
  if ! cp -a --reflink=never -- "$source" "$target"; then
    printf 'Canonical analysis install failed: %s\n' "$target" >&2
    rm -rf -- "$target"
    if [[ -e "$backup" || -L "$backup" ]]; then
      mv -- "$backup" "$target"
    fi
    return 1
  fi
  printf '[OK] canonical analysis installed: %s\n' "$relative_target"

}
# Keep one truthful terminal bar across the thirteen figure packages. A unit advances
# only after that figure's renderer and strict output validator both succeed.
FIGURE_TOTAL=13
FIGURE_PROGRESS_FD=""
FIGURE_PROGRESS_PID_VALUE=""
CURRENT_FIGURE=""

start_figure_progress() {
  local progress_read_fd
  coproc FIGURE_PROGRESS_PROCESS {
    "$PY" scripts/utilities/08_figure_progress.py \
      --total "$FIGURE_TOTAL" >/dev/null
  }
  FIGURE_PROGRESS_FD="${FIGURE_PROGRESS_PROCESS[1]}"
  FIGURE_PROGRESS_PID_VALUE="$FIGURE_PROGRESS_PROCESS_PID"
  progress_read_fd="${FIGURE_PROGRESS_PROCESS[0]}"
  exec {progress_read_fd}<&-
}

figure_progress_event() {
  printf '%s\t%s\n' "$1" "$2" >&"$FIGURE_PROGRESS_FD"
}

figure_progress_start() {
  CURRENT_FIGURE="$1"
  figure_progress_event START "$CURRENT_FIGURE"
}

figure_progress_done() {
  figure_progress_event DONE "$CURRENT_FIGURE"
  CURRENT_FIGURE=""
}

finish_figure_progress() {
  local helper_status=0
  if [[ -n "$FIGURE_PROGRESS_FD" ]]; then
    exec {FIGURE_PROGRESS_FD}>&-
    FIGURE_PROGRESS_FD=""
  fi
  if [[ -n "$FIGURE_PROGRESS_PID_VALUE" ]]; then
    wait "$FIGURE_PROGRESS_PID_VALUE" || helper_status=$?
    FIGURE_PROGRESS_PID_VALUE=""
  fi
  return "$helper_status"
}

cleanup_figure_progress() {
  local workflow_status=$?
  trap - EXIT
  if [[ -n "${FIGURE_PROGRESS_FD:-}" ]]; then
    if ((workflow_status != 0)); then
      figure_progress_event FAIL "${CURRENT_FIGURE:-reproduction workflow}" || true
    fi
    exec {FIGURE_PROGRESS_FD}>&- || true
  fi
  if [[ -n "${FIGURE_PROGRESS_PID_VALUE:-}" ]]; then
    wait "$FIGURE_PROGRESS_PID_VALUE" || true
  fi
  if [[ -n "${ACTIVE_INTERNAL_STAGE:-}" ]]; then
    printf 'Failed isolated stage preserved for diagnosis: %s\n' \
      "$ACTIVE_INTERNAL_STAGE" >&2
  fi
  if [[ -n "${BUILD_ROOT:-}" && -d "${BUILD_ROOT:-}" ]]; then
    printf 'Failed full rebuild stage preserved for diagnosis: %s\n' \
      "$BUILD_ROOT" >&2
  fi
  exit "$workflow_status"
}

trap cleanup_figure_progress EXIT
start_figure_progress

printf '%s\n' '' '== Step 1: Figure 1 =='
figure_progress_start 'Figure 1'
"$PY" Fig1/01_analyze_behavior_experiment1.py \
  --consumption "$PAPER_ROOT/Fig1/raw_data/consumption.csv" \
  --weights "$PAPER_ROOT/Fig1/raw_data/weights.csv" \
  --outdir "$OUT/Fig1/analysis"
"$PY" Fig1/02_make_figure_1_behavior.py \
  --analysis-dir "$OUT/Fig1/analysis" \
  --output-dir "$OUT/Fig1/figure" \
  --dpi 600
validate_render Fig1 "$OUT/Fig1/figure"
figure_progress_done

printf '%s\n' '' '== Step 2: Figure 2 =='
figure_progress_start 'Figure 2'
start_internal_stage fig2
FIG2_INTERNAL_STAGE="$ACTIVE_INTERNAL_STAGE"
./Fig2/07_run_all.sh \
  "$FIG2_INTERNAL_STAGE/stage" \
  "$FIG2_INTERNAL_STAGE/render"
validate_render Fig2 "$FIG2_INTERNAL_STAGE/render"
mkdir -m 700 -- "$OUT/Fig2"
cp -a --reflink=never -- "$FIG2_INTERNAL_STAGE/stage" "$OUT/Fig2/stage"
cp -a --reflink=never -- "$FIG2_INTERNAL_STAGE/render" "$OUT/Fig2/render"
validate_render Fig2 "$OUT/Fig2/render"
remove_internal_stage
figure_progress_done

printf '%s\n' '' '== Step 3: Figure 3 (certified; promotion requires --force) =='
figure_progress_start 'Figure 3'
start_internal_stage fig3
FIG3_INTERNAL_STAGE="$ACTIVE_INTERNAL_STAGE"
mkdir -p "$FIG3_INTERNAL_STAGE/work"
cp -r --reflink=never \
  Fig3/raw_data/legacy_experiment_2025_07_28 \
  "$FIG3_INTERNAL_STAGE/work/raw"
"$PY" Fig3/reproduce_final_20260823/prepare_raw_project.py \
  --project-root "$FIG3_INTERNAL_STAGE/work/project" \
  --copy-from "$FIG3_INTERNAL_STAGE/work/raw" \
  --receipt "$FIG3_INTERNAL_STAGE/work/prepare_receipt.json"
FIG3_UTILITY_PYTHON="$REVIEW_PY" PATH="$REVIEW_BIN:$PATH" \
  ./Fig3/reproduce_final_20260823/RUN_REPRODUCE_FINAL.sh \
  "$FIG3_INTERNAL_STAGE/work/project" \
  "$FIG3_INTERNAL_STAGE/work/out" \
  "$PAPER_ROOT/Fig3/hil_review/ventricle_cartoon_20260827_v1"
validate_render Fig3 "$FIG3_INTERNAL_STAGE/work/out/final"
mkdir -p -m 700 -- "$OUT/Fig3/out"
cp -a --reflink=never -- \
  "$FIG3_INTERNAL_STAGE/work/out/final" "$OUT/Fig3/out/final"
cp -a --reflink=never -- \
  "$FIG3_INTERNAL_STAGE/work/out/accepted_ventricle_hil" \
  "$OUT/Fig3/out/accepted_ventricle_hil"
cp --reflink=never -- \
  "$FIG3_INTERNAL_STAGE/work/out/RAW_PROJECT_VALIDATION.json" \
  "$OUT/Fig3/out/RAW_PROJECT_VALIDATION.json"
cp --reflink=never -- \
  "$FIG3_INTERNAL_STAGE/work/prepare_receipt.json" \
  "$OUT/Fig3/prepare_receipt.json"
mkdir -p -m 700 -- "$OUT/Fig3"
cp -a --reflink=never -- \
  "$FIG3_INTERNAL_STAGE/work/out/Paper/analyses/Fig3/results" \
  "$OUT/Fig3/analysis"
validate_render Fig3 "$OUT/Fig3/out/final"
remove_internal_stage
figure_progress_done

printf '%s\n' '' '== Step 4: Figure 4 =='
figure_progress_start 'Figure 4'
"$PY" Fig4/02_analyze_pomc_cfos.py \
  --output "$OUT/Fig4/analysis" --no-auto-conda </dev/null
"$PY" Fig4/07_audit_pomc_size_sensitivity.py --analysis-dir "$OUT/Fig4/analysis"
"$PY" Fig4/05_analyze_spatial_distributions.py \
  --analysis-dir "$OUT/Fig4/analysis" \
  --output-dir "$OUT/Fig4/analysis/spatial" \
  --dpi 600
"$PY" Fig4/03_make_figure_4_pomc_cfos.py \
  --analysis-dir "$OUT/Fig4/analysis" \
  --outdir "$OUT/Fig4/figure" \
  --figure-name Figure_4 \
  --dpi 600
validate_render Fig4 "$OUT/Fig4/figure"
figure_progress_done

printf '%s\n' '' '== Step 5: S5 from anatomy HIL plus reviewed per-acquisition 3V exclusion =='
figure_progress_start 'Figure S5'
"$PY" FigS5/06_include_npy_m.py --root "$PAPER_ROOT"
"$PY" FigS5/02a_july_cfos_acth_clip_human_in_loop_s5.py \
  --prepare \
  --sample-manifest FigS5/provenance/inclusion_NPY_M_20260908/sample_manifest.csv \
  --channel-map FigS5/provenance/inclusion_NPY_M_20260908/channel_map.csv \
  --export-dir FigS5/channel_tiffs \
  --workdir "$OUT/FigS5/work"
"$PY" FigS5/02a_july_cfos_acth_clip_human_in_loop_s5.py \
  --audit-segmentations \
  --export-dir FigS5/channel_tiffs \
  --workdir "$OUT/FigS5/work"
"$PY" scripts/utilities/06_import_figure_s5_hil.py \
  --from FigS5/hil_review/human_final_including_NPY_M_20260908 \
  --workdir "$OUT/FigS5/work"
"$PY" FigS5/02a_july_cfos_acth_clip_human_in_loop_s5.py \
  --quantify \
  --export-dir FigS5/channel_tiffs \
  --workdir "$OUT/FigS5/work" \
  --qc-format both \
  --qc-dpi 600 \
  --plot-dpi 600
"$PY" FigS5/03_make_figure_s5_acth_clip_cfos.py \
  --workdir "$OUT/FigS5/work" \
  --spatial-hil-dir "$FIGS5_SPATIAL_HIL" \
  --output-dir "$OUT/FigS5/figure" \
  --dpi 600
validate_render FigS5 "$OUT/FigS5/figure"
figure_progress_done

printf '%s\n' '' '== Step 6: Figure 5 =='
figure_progress_start 'Figure 5'
MICROGLIAL_CHOICE_ARGS=()
if ((DO_MICROGLIAL_CHOICES)); then
  MICROGLIAL_CHOICE_ARGS+=(--do_microglial_choices)
fi
"$PY" Fig5/12_microglial_choices.py \
  --output-dir "$OUT/Fig5/microglial_choices" "${MICROGLIAL_CHOICE_ARGS[@]}"
"$PY" Fig5/01_analyze_gfap_iba1_microglia.py \
  --manual-region-dir "$PAPER_ROOT/Fig5/hil_review/human_final_20260823_v1" \
  --reviewed-dapi-root "$PAPER_ROOT/Fig5/hil_review/dapi_final_20260830_v1" \
  --reviewed-iba1-root "$PAPER_ROOT/Fig5/hil_review/iba1_final_20260830_v1" \
  --cnn-device cpu \
  --output "$OUT/Fig5/analysis" --no-auto-conda </dev/null
"$PY" Fig5/02_make_figure_5_gfap_iba1_microglia.py \
  --analysis-dir "$OUT/Fig5/analysis" \
  --output-dir "$OUT/Fig5/figure" \
  --dpi 600
validate_render Fig5 "$OUT/Fig5/figure"
figure_progress_done

printf '%s\n' '' '== Step 7: Figure S1 =='
figure_progress_start 'Figure S1'
"$PY" FigS1/01_analyze_single_bottle_male.py \
  --consumption "$PAPER_ROOT/FigS1/raw_data/consumption.csv" \
  --weights "$PAPER_ROOT/FigS1/raw_data/weights.csv" \
  --figure1-reference-dir "$OUT/Fig1/analysis" \
  --outdir "$OUT/FigS1/analysis"
"$PY" FigS1/02_make_figure_s1_male.py \
  --analysis-dir "$OUT/FigS1/analysis" \
  --output-dir "$OUT/FigS1/figure" \
  --dpi 600
validate_render FigS1 "$OUT/FigS1/figure"
figure_progress_done

printf '%s\n' '' '== Step 8: Figure S2 =='
figure_progress_start 'Figure S2'
"$PY" FigS2/01_analyze_single_bottle_female.py \
  --consumption "$PAPER_ROOT/FigS2/raw_data/consumption.csv" \
  --weights "$PAPER_ROOT/FigS2/raw_data/weights.csv" \
  --figure1-reference-dir "$OUT/Fig1/analysis" \
  --outdir "$OUT/FigS2/analysis"
"$PY" FigS2/02_make_figure_s2_female.py \
  --analysis-dir "$OUT/FigS2/analysis" \
  --output-dir "$OUT/FigS2/figure" \
  --dpi 600
validate_render FigS2 "$OUT/FigS2/figure"
figure_progress_done

printf '%s\n' '' '== Step 9: Figure S3 =='
figure_progress_start 'Figure S3'
PATH="$REVIEW_BIN:$PATH" "$REVIEW_PY" FigS3/05_analyze_make_figure.py \
  --figs3-root "$PAPER_ROOT/FigS3" \
  --output-root "$OUT/FigS3/figure" \
  --permutations 99999 \
  --dpi 600
validate_render FigS3 "$OUT/FigS3/figure"
figure_progress_done

printf '%s\n' '' '== Step 10: Figure S4 =='
figure_progress_start 'Figure S4'
"$PY" FigS4/01_analyze_behavior_experiment2.py \
  --output-dir "$OUT/FigS4/analysis"
"$PY" FigS4/02_make_figure_s4_behavior_preference.py \
  --analysis-dir "$OUT/FigS4/analysis" \
  --output-dir "$OUT/FigS4/figure" \
  --dpi 600
validate_render FigS4 "$OUT/FigS4/figure"
figure_progress_done

printf '%s\n' '' '== Step 11: Figure S6 =='
figure_progress_start 'Figure S6'
"$PY" FigS6/01_make_figure_s6_microglia_classifier.py \
  --output-dir "$OUT/FigS6/figure" \
  --dpi 600
validate_render FigS6 "$OUT/FigS6/figure"
figure_progress_done

printf '%s\n' '' '== Step 12: Figure S7 =='
figure_progress_start 'Figure S7'
"$PY" FigS7/01_make_figure_s7_training.py --output-dir "$OUT/FigS7/figure" --dpi 600
validate_render FigS7 "$OUT/FigS7/figure"
figure_progress_done

printf '%s\n' '' '== Step 13: Figure S8 =='
figure_progress_start 'Figure S8'
"$PY" FigS8/01_make_figure_s8_power.py --output-dir "$OUT/FigS8/figure" --dpi 600
validate_render FigS8 "$OUT/FigS8/figure"
figure_progress_done

finish_figure_progress

printf '%s\n' '' '== Byte-for-byte reproduction against shipped raster images =='
verify_reference_images Fig1 "$OUT/Fig1/figure"
verify_reference_images Fig2 "$OUT/Fig2/render/outputs"
verify_reference_images Fig3 "$OUT/Fig3/out/final"
verify_reference_images Fig4 "$OUT/Fig4/figure"
verify_reference_images FigS5 "$OUT/FigS5/figure"
verify_reference_images Fig5 "$OUT/Fig5/figure"
verify_reference_images FigS1 "$OUT/FigS1/figure"
verify_reference_images FigS2 "$OUT/FigS2/figure"
verify_reference_images FigS3 "$OUT/FigS3/figure"
verify_reference_images FigS4 "$OUT/FigS4/figure"
verify_reference_images FigS6 "$OUT/FigS6/figure"
verify_reference_images FigS7 "$OUT/FigS7/figure"
verify_reference_images FigS8 "$OUT/FigS8/figure"

printf '%s\n' '' '== All rebuilds passed: install canonical figure packages =='
promote_canonical Fig1 "$OUT/Fig1/figure"
promote_canonical Fig2 "$OUT/Fig2/render/outputs"
if ((FORCE_FIG3)); then
  promote_canonical Fig3 "$OUT/Fig3/out/final" --force
fi
promote_canonical Fig4 "$OUT/Fig4/figure"
promote_canonical FigS5 "$OUT/FigS5/figure"
promote_canonical Fig5 "$OUT/Fig5/figure"
promote_canonical FigS1 "$OUT/FigS1/figure"
promote_canonical FigS2 "$OUT/FigS2/figure"
promote_canonical FigS3 "$OUT/FigS3/figure"
promote_canonical FigS4 "$OUT/FigS4/figure"
promote_canonical FigS6 "$OUT/FigS6/figure"
promote_canonical FigS7 "$OUT/FigS7/figure"
promote_canonical FigS8 "$OUT/FigS8/figure"

printf '%s\n' '' '== Byte-for-byte staged-to-canonical publication audit =='
verify_canonical_bytes Fig1 "$OUT/Fig1/figure"
verify_canonical_bytes Fig2 "$OUT/Fig2/render/outputs"
if ((FORCE_FIG3)); then
  verify_canonical_bytes Fig3 "$OUT/Fig3/out/final"
else
  printf '%s\n' 'Fig3 byte audit skipped because certified Figure 3 was not installed; pass --force to install and audit it.'
fi
verify_canonical_bytes Fig4 "$OUT/Fig4/figure"
verify_canonical_bytes FigS5 "$OUT/FigS5/figure"
verify_canonical_bytes Fig5 "$OUT/Fig5/figure"
verify_canonical_bytes FigS1 "$OUT/FigS1/figure"
verify_canonical_bytes FigS2 "$OUT/FigS2/figure"
verify_canonical_bytes FigS3 "$OUT/FigS3/figure"
verify_canonical_bytes FigS4 "$OUT/FigS4/figure"
verify_canonical_bytes FigS6 "$OUT/FigS6/figure"
verify_canonical_bytes FigS7 "$OUT/FigS7/figure"
verify_canonical_bytes FigS8 "$OUT/FigS8/figure"

printf '%s\n' '' '== Install fresh generated analyses into the canonical tree =='
replace_canonical_analysis "$OUT/Fig1/analysis" "$PAPER_ROOT/analyses/Fig1/results"
replace_canonical_analysis "$OUT/Fig2/stage/analysis" "$PAPER_ROOT/analyses/Fig2/results"
replace_canonical_analysis "$OUT/Fig3/analysis" "$PAPER_ROOT/analyses/Fig3/results"
replace_canonical_analysis "$OUT/Fig4/analysis" "$PAPER_ROOT/analyses/Fig4/results/human_final_run"
replace_canonical_analysis "$OUT/Fig4/figure" "$PAPER_ROOT/analyses/Fig4/figure/human_final_render"
replace_canonical_analysis "$OUT/FigS5/work" "$PAPER_ROOT/analyses/FigS5/results/including_NPY_M_20260908"
replace_canonical_analysis "$OUT/Fig5/analysis" "$PAPER_ROOT/analyses/Fig5/results/human_final_run"
replace_canonical_analysis "$OUT/Fig5/figure" "$PAPER_ROOT/analyses/Fig5/figure/human_final_render"
replace_canonical_analysis "$OUT/Fig5/microglial_choices" "$PAPER_ROOT/analyses/Fig5/results/current_microglial_choices"
replace_canonical_analysis "$OUT/FigS1/analysis" "$PAPER_ROOT/analyses/FigS1/results"
replace_canonical_analysis "$OUT/FigS2/analysis" "$PAPER_ROOT/analyses/FigS2/results"
replace_canonical_analysis "$OUT/FigS4/analysis" "$PAPER_ROOT/analyses/FigS3/results"
replace_canonical_analysis "$OUT/FigS6/figure" "$PAPER_ROOT/analyses/FigS6/figure/current"
replace_canonical_analysis "$OUT/FigS7/figure" "$PAPER_ROOT/analyses/FigS7/figure/current"
replace_canonical_analysis "$OUT/FigS8/figure" "$PAPER_ROOT/analyses/FigS8/figure/current"

printf '%s\n' '' '== Complete Figure 3 with accepted sucrose and bilingual subpanels =='
FIG3_RENDER_PYTHON="$PY" FIG3_PDF_PYTHON="$FIG3_PDF_BIN" \
  ./Fig3/07_rebuild_complete_figure3.sh --base-panels-dir "$OUT/Fig3/out/final/panels"

printf '%s\n' '' '== Final canonical tree and immutable-input validation =='
"$PY" scripts/utilities/09_audit_canonical_inputs.py \
  --root "$PAPER_ROOT" --verify "$RAW_INPUT_SNAPSHOT"
"$PY" scripts/utilities/07_validate_figure_outputs.py --all-live
"$PY" scripts/utilities/03_check_reproducibility.py
"$PY" scripts/utilities/01_validate_bundle.py

case "$BUILD_ROOT" in
  /tmp/apotome_full_rebuild_*) rm -rf -- "$BUILD_ROOT" ;;
  *) printf 'Refusing to remove unexpected build stage: %s\n' "$BUILD_ROOT" >&2; exit 2 ;;
esac
BUILD_ROOT=""
trap - EXIT
printf '%s\n' '' '============================================================'
printf '[PASS] Reproduction workflow completed.\n'
printf 'Canonical output: %s\n' "$PAPER_ROOT"
printf '%s\n' 'Raw data and active scripts are present and unchanged.'
printf 'Previous generated analyses are recoverable at: %s\n' "$INSTALL_RECOVERY_ROOT"
if ((FORCE_FIG3)); then
  printf '%s\n' 'All thirteen figures were installed, including certified Figure 3 (--force).'
else
  printf '%s\n' 'Ten figures were installed; certified Figure 3 measurements were retained and its three-condition presentation was rebuilt.'
fi
