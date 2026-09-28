#!/usr/bin/env bash
# The supported Figure 3 rebuild uses the eleven-animal cohort: the original
# nine animals plus the two 2026 Zeiss LSM 780 Water controls.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
paper_root="$(cd -- "$script_dir" && while [[ ! -d scripts/shared && "$PWD" != / ]]; do cd ..; done; pwd -P)"
[[ -d "$paper_root/scripts/shared" ]] || { printf 'Cannot locate the Paper root\n' >&2; exit 2; }
exec bash "$paper_root/scripts/run_fig3_new_water_cohort.sh" "$@"
