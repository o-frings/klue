#!/bin/zsh
# Full empirical re-run for version11: 5 adapters (LCMNL C=1..6 + independent
# MMNL) -> extend-C search -> correlated MMNL, all on the availability-masked /
# canonical samples. v10 outputs are backed up first for diffing; the resumable
# CSVs (corr-MMNL, extend-C) are cleared so they regenerate on the new samples.
set -e
cd "$(dirname "$0")/.."   # repo root, wherever the repo lives
export KLUE_CORES=4

# First run only: back up the v10 outputs and clear the resumable CSVs so they
# regenerate. On a re-launch (backup already populated) skip both, so the v10
# reference is never clobbered and the resumable steps continue where they left
# off.
if [ ! -d output/v10_backup ] || [ -z "$(ls -A output/v10_backup 2>/dev/null)" ]; then
  mkdir -p output/v10_backup
  cp output/*_lcmnl_results.csv output/*_model_comparison.csv output/*_class_betas.csv \
     output/empirical_corr_mmnl_3000draws.csv output/emp_cmax_extend_*.csv \
     output/v10_backup/ 2>/dev/null || true
  rm -f output/empirical_corr_mmnl_3000draws.csv output/emp_cmax_extend_*.csv
  echo "backed up v10 outputs; cleared resumable CSVs"
else
  echo "v10_backup exists; re-launch mode (resuming, not re-backing-up)"
fi

echo "===== ADAPTERS (LCMNL + independent MMNL) ====="
for s in application mode_choice swiss_route electricity swissmetro; do
  echo "----- adapter: $s ($(date +%H:%M:%S)) -----"
  Rscript R/empirical_$s.R || echo "!!! adapter $s exited non-zero"
done

echo "===== EXTEND-C ($(date +%H:%M:%S)) ====="
Rscript dev/run_emp_cmax_extend.R || echo "!!! extend-C exited non-zero"

echo "===== CORRELATED MMNL ($(date +%H:%M:%S)) ====="
Rscript dev/empirical_corr_mmnl.R || echo "!!! corr-MMNL exited non-zero"

echo "===== RE-RUN COMPLETE ($(date +%H:%M:%S)) ====="
