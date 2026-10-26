#!/usr/bin/env bash
# Run the validated H2A.Z analysis stages in dependency order from the project root.
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
source config/tools.sh
export THREADS

# Verify shared helper invariants before accepting the sample inventory.
"$H2AZ_RSCRIPT" tests/selftest.R
"$H2AZ_RSCRIPT" tests/inventory_selftest.R
"$H2AZ_RSCRIPT" scripts/01_validate_inventory.R
# Call peaks only for included samples that lack a supplied peak path.
JOBS=results/01_inventory/peak_jobs.tsv
if [[ -f "$JOBS" ]]; then
  bash scripts/02_call_peaks.sh "$JOBS"
fi
# Build one fixed region set, then analyze each mark against that same set.
"$H2AZ_RSCRIPT" scripts/03_build_master_regions.R
for mark in H2AZ H2AZac H2AZub; do
  "$H2AZ_RSCRIPT" scripts/04_diffbind_per_mark.R "$mark"
done
# Join mark-level calls, export figures/BEDs, and prepare manual-review candidates.
"$H2AZ_RSCRIPT" scripts/05_joint_states.R
"$H2AZ_RSCRIPT" scripts/06_plot_and_export.R
"$H2AZ_RSCRIPT" scripts/07_review_loci.R
