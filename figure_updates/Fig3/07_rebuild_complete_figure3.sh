#!/usr/bin/env bash
# Finish the certified quantitative build with the accepted third-condition cartoon.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
paper_root="$(cd -- "$script_dir" && while [[ ! -d scripts/shared && "$PWD" != / ]]; do cd ..; done; pwd -P)"
[[ -d "$paper_root/scripts/shared" ]] || { printf 'Cannot locate the Paper root\n' >&2; exit 2; }
render_python="${FIG3_RENDER_PYTHON:-python3}"
pdf_python="${FIG3_PDF_PYTHON:-$render_python}"
"$render_python" "$paper_root/Fig3/02c_make_sucrose_cartoon.py" --root "$paper_root"
"$pdf_python" "$paper_root/Fig3/04_complete_figure3_with_sucrose.py" --root "$paper_root" "$@"
"$pdf_python" "$paper_root/scripts/utilities/07_validate_figure_outputs.py" --figure Fig3 --root "$paper_root/Fig3"
