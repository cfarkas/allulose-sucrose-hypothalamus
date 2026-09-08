#!/bin/bash -p
set -euo pipefail
IFS=$'\n\t'
umask 077
unset PYTHONPATH PYTHONHOME BASH_ENV ENV CDPATH GLOBIGNORE
export PYTHONDONTWRITEBYTECODE=1
export MPLBACKEND=Agg

die() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 64
}

[[ "$#" -eq 3 ]] || die "usage: $0 /ABSOLUTE/RAW_PROJECT/Paper /tmp/ABSENT_OUTPUT /ABSOLUTE/ACCEPTED_VENTRICLE_HIL"
RAW_PROJECT="$1"
OUTPUT_ROOT="$2"
VENTRICLE_HIL_SOURCE="$3"
[[ "$RAW_PROJECT" = /* ]] || die "raw project must be absolute"
[[ "$OUTPUT_ROOT" = /tmp/* ]] || die "output must be an absolute /tmp child"
[[ "$VENTRICLE_HIL_SOURCE" = /* ]] || die "accepted ventricle HIL directory must be absolute"
[[ -d "$RAW_PROJECT" && ! -L "$RAW_PROJECT" ]] || die "raw project missing or linked"
[[ ! -e "$OUTPUT_ROOT" && ! -L "$OUTPUT_ROOT" ]] || die "output already exists"
[[ -d "$VENTRICLE_HIL_SOURCE" && ! -L "$VENTRICLE_HIL_SOURCE" ]] ||
  die "accepted ventricle HIL directory missing or linked"

SCRIPT_DIR="$(cd -- "$(dirname -- "$0")" && pwd -P)"
if [[ -n "${FIG3_PYTHON:-}" ]]; then
  PYTHON="$FIG3_PYTHON"
elif [[ -n "${PAPER_PYTHON:-}" ]]; then
  PYTHON="$PAPER_PYTHON"
else
  PYTHON="/home/server/anaconda3/envs/paper_apotome_repro/bin/python"
fi
[[ -x "$PYTHON" ]] || die "Python is not executable: $PYTHON"
if [[ -n "${FIG3_UTILITY_PYTHON:-}" ]]; then
  UTILITY_PYTHON="$FIG3_UTILITY_PYTHON"
elif "$PYTHON" -c 'import pypdf' >/dev/null 2>&1; then
  UTILITY_PYTHON="$PYTHON"
else
  UTILITY_PYTHON="/home/server/anaconda3/bin/python"
fi
[[ -x "$UTILITY_PYTHON" ]] || die "Utility Python is not executable: $UTILITY_PYTHON"
"$UTILITY_PYTHON" -c 'import pypdf' >/dev/null 2>&1 ||
  die "Utility Python must import pypdf; set FIG3_UTILITY_PYTHON"

(cd -- "$SCRIPT_DIR" && sha256sum -c SHA256_MANIFEST.tsv)
mkdir -m 700 -- "$OUTPUT_ROOT"
VENTRICLE_HIL_STAGE="$OUTPUT_ROOT/accepted_ventricle_hil"
mkdir -m 700 -- "$VENTRICLE_HIL_STAGE"
for name in \
  water_ventricle_mask.png water_ventricle_receipt.json \
  allulose_ventricle_mask.png allulose_ventricle_receipt.json
do
  source_file="$VENTRICLE_HIL_SOURCE/$name"
  staged_file="$VENTRICLE_HIL_STAGE/$name"
  [[ -f "$source_file" && ! -L "$source_file" ]] || die "missing or linked accepted HIL file: $source_file"
  cp --reflink=never -- "$source_file" "$staged_file"
  cmp --silent -- "$source_file" "$staged_file" || die "accepted HIL copy mismatch: $name"
  chmod 400 -- "$staged_file"
done
STAGE_PAPER="$OUTPUT_ROOT/Paper"
SHADOW_SCRIPTS="$STAGE_PAPER/scripts/Fig3"
mkdir -p -m 700 -- "$SHADOW_SCRIPTS"
for name in   00_extract_water3_neun_display.py   01_analyze_cfos_npy.py   02_make_cfos_npy_cartoons.py   05_analyze_spatial_distributions.py   03_make_figure_3_cfos_npy.py   02_make_cfos_npy_cartoons_spanish.py   05_analyze_spatial_distributions_spanish.py   03_make_figure_3_cfos_npy_spanish.py
do
  cp --reflink=never -- "$SCRIPT_DIR/scripts/$name" "$SHADOW_SCRIPTS/$name"
  cmp --silent -- "$SCRIPT_DIR/scripts/$name" "$SHADOW_SCRIPTS/$name" ||
    die "shadow script copy mismatch: $name"
done

"$PYTHON" "$SCRIPT_DIR/prepare_raw_project.py" \
  --project-root "$RAW_PROJECT" \
  --receipt "$OUTPUT_ROOT/RAW_PROJECT_VALIDATION.json" \
  --validate-only

RAW_ROOT="$RAW_PROJECT/Fig3/raw/legacy_experiment_2025_07_28"
ANALYSIS="$STAGE_PAPER/analyses/Fig3/results/final_run"
SPATIAL_EN="$STAGE_PAPER/analyses/Fig3/results/spatial_english"
PUBLICATION_EN="$STAGE_PAPER/Fig3/english_publication"
ANALYSIS_ES="$STAGE_PAPER/analyses/Fig3/results/spanish_analysis_work"
PANELS_ES="$ANALYSIS_ES/panels"
SPATIAL_ES="$STAGE_PAPER/analyses/Fig3/results/spatial_spanish"
PUBLICATION_ES="$STAGE_PAPER/Fig3/spanish_publication_work"
FINAL="$OUTPUT_ROOT/final"

"$PYTHON" "$SHADOW_SCRIPTS/01_analyze_cfos_npy.py" \
  --raw-root "$RAW_ROOT" --output-dir "$ANALYSIS"
"$PYTHON" "$SHADOW_SCRIPTS/02_make_cfos_npy_cartoons.py" \
  --root "$RAW_PROJECT" --output-dir "$ANALYSIS/panels" --dpi 600 \
  --ventricle-hil-dir "$VENTRICLE_HIL_STAGE"
"$PYTHON" "$SHADOW_SCRIPTS/05_analyze_spatial_distributions.py" \
  --raw-root "$RAW_ROOT" --output-dir "$SPATIAL_EN" --dpi 600
"$PYTHON" "$SHADOW_SCRIPTS/03_make_figure_3_cfos_npy.py" \
  --analysis-dir "$ANALYSIS" --output-dir "$PUBLICATION_EN" \
  --raw-root "$RAW_ROOT" --spatial-dir "$SPATIAL_EN" \
  --per-animal "$ANALYSIS/per_animal_cfos_npy.csv" \
  --statistics "$ANALYSIS/anova_mwu_statistics.csv" \
  --panel-b "$ANALYSIS/panels/Figure3_Panel_B_Water_marker_positive_nuclei.png" \
  --panel-c "$ANALYSIS/panels/Figure3_Panel_C_Allulose_marker_positive_nuclei.png" \
  --panel-f "$SPATIAL_EN/Figure3_Spatial_cFOS_occurrence.png" \
  --panel-g "$SPATIAL_EN/Figure3_Spatial_cFOS_NPY_occurrence.png" \
  --figure-name Figure_3_cFos_NPY --dpi 600 --panel-a-dpi 600 --seed 31

"$PYTHON" "$SHADOW_SCRIPTS/02_make_cfos_npy_cartoons_spanish.py" \
  --root "$RAW_PROJECT" --output-dir "$PANELS_ES" --dpi 600 \
  --ventricle-hil-dir "$VENTRICLE_HIL_STAGE"
"$PYTHON" "$SHADOW_SCRIPTS/05_analyze_spatial_distributions_spanish.py" \
  --raw-root "$RAW_ROOT" --output-dir "$SPATIAL_ES" --dpi 600
"$PYTHON" "$SHADOW_SCRIPTS/03_make_figure_3_cfos_npy_spanish.py" \
  --analysis-dir "$ANALYSIS_ES" --output-dir "$PUBLICATION_ES" \
  --raw-root "$RAW_ROOT" --spatial-dir "$SPATIAL_ES" \
  --per-animal "$ANALYSIS/per_animal_cfos_npy.csv" \
  --statistics "$ANALYSIS/anova_mwu_statistics.csv" \
  --panel-b "$PANELS_ES/Figure3_Panel_B_Water_marker_positive_nuclei.png" \
  --panel-c "$PANELS_ES/Figure3_Panel_C_Allulose_marker_positive_nuclei.png" \
  --panel-f "$SPATIAL_ES/Figure3_Spatial_cFOS_occurrence.png" \
  --panel-g "$SPATIAL_ES/Figure3_Spatial_cFOS_NPY_occurrence.png" \
  --figure-name Figure_3_cFos_NPY_spanish --dpi 600 --panel-a-dpi 600 --seed 31

"$UTILITY_PYTHON" "$SCRIPT_DIR/assemble_final.py" \
  --english-publication "$PUBLICATION_EN" \
  --spanish-publication "$PUBLICATION_ES" \
  --output "$FINAL" \
  --work-dir "$OUTPUT_ROOT/normalization_work"
"$UTILITY_PYTHON" "$SCRIPT_DIR/validate_final.py" \
  --final "$FINAL" \
  --raw-validation "$OUTPUT_ROOT/RAW_PROJECT_VALIDATION.json"

printf '[PASS] Figure 3 accepted bilingual reproduction: %s\n' "$FINAL"
