#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
paper_root="$(cd -- "$script_dir" && while [[ ! -d scripts/shared && "$PWD" != / ]]; do cd ..; done; pwd -P)"
[[ -d "$paper_root/scripts/shared" ]] || { printf 'Cannot locate the Paper root above %s\n' "$script_dir" >&2; exit 2; }
driver="${paper_root}/scripts/shared/04_safe_stage_fig3.py"
python_bin="${FIG3_STAGE_PYTHON:-/home/server/anaconda3/envs/paper_apotome_repro/bin/python}"
default_ventricle_hil="${paper_root}/Fig3/hil_review/ventricle_cartoon_20260827_v1"

if [[ ! -f "${driver}" ]]; then
    echo "REFUSED: missing staged Figure 3 driver: ${driver}" >&2
    exit 2
fi
if [[ ! -x "${python_bin}" ]]; then
    echo "REFUSED: FIG3 staging Python is not executable: ${python_bin}" >&2
    exit 2
fi

hil_option_present=false
for argument in "$@"; do
    if [[ "${argument}" == "--ventricle-hil-dir" || "${argument}" == --ventricle-hil-dir=* ]]; then
        hil_option_present=true
        break
    fi
done

if [[ "${hil_option_present}" == true ]]; then
    exec "${python_bin}" "${driver}" --paper-root "${paper_root}" "$@"
fi
exec "${python_bin}" "${driver}" --paper-root "${paper_root}" \
    --ventricle-hil-dir "${default_ventricle_hil}" "$@"
