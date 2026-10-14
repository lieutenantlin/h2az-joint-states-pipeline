#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
source config/tools.sh
export THREADS

"$H2AZ_RSCRIPT" tests/selftest.R
"$H2AZ_RSCRIPT" scripts/01_validate_inventory.R
JOBS=results/01_inventory/peak_jobs.tsv
if [[ -f "$JOBS" ]]; then
  bash scripts/02_call_peaks.sh "$JOBS"
fi
"$H2AZ_RSCRIPT" scripts/03_build_master_regions.R
for mark in H2AZ H2AZac H2AZub; do
  "$H2AZ_RSCRIPT" scripts/04_diffbind_per_mark.R "$mark"
done
"$H2AZ_RSCRIPT" scripts/05_joint_states.R
"$H2AZ_RSCRIPT" scripts/06_plot_and_export.R
"$H2AZ_RSCRIPT" scripts/07_review_loci.R
