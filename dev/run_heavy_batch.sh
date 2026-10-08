#!/usr/bin/env bash
# =============================================================================
# Sequential, low-heat batch of the parked revision reruns.
#
# Thermal policy: single-threaded BLAS, Apollo on one core, nice -19, and a
# cooling pause between steps -- only one core is ever warm, so the laptop stays
# cool at the cost of wall-clock time. Every step is resumable, so re-running
# continues where it stopped (and a failure in one step does not abort the rest).
#
# Excludes the gmnl head-to-head (Minor 6 stress ladder, Gate 1 Apollo) because
# gmnl needs a pinned mlogit (<1.1) side-library; run those separately once
# KLUE_OLD_MLOGIT_LIB is set up.
#
# Launch (detached):
#   nohup bash dev/run_heavy_batch.sh > /tmp/heavy_batch.log 2>&1 &
# Watch:  tail -f /tmp/heavy_batch.log
# Cooling between steps defaults to 120s; override with COOL=<seconds>.
# =============================================================================
set -u
cd "$(dirname "$0")/.."

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export KLUE_CORES=1
COOL=${COOL:-120}

RUN() {
  echo "[$(date +%H:%M:%S)] >>> $*"
  nice -n 19 "$@"
  echo "[$(date +%H:%M:%S)] <<< exit $? ; cooling ${COOL}s"
  sleep "$COOL"
}

echo "==== HEAVY BATCH START $(date) ===="

# --- Gate 4: correlated MMNL, one dataset per call (resumable; cool between) ---
for D in Mode SwissRoute Electricity Vittel Swissmetro; do
  RUN Rscript dev/empirical_corr_mmnl.R 3000 "$D"
done

# --- Gate 5: extend C until BIC turns; 20s cooling between C-fits, cool between datasets ---
export KLUE_SLEEP=20
for D in Electricity Vittel Swissmetro; do
  RUN Rscript dev/run_emp_cmax_extend.R "$D" 12
done

# --- Minor 7: full K*=1 blocked grid (its own cooling via KLUE_SLEEP) ---
KLUE_SLEEP=30 RUN Rscript dev/run_k1_mmnl_blocked.R

echo "==== HEAVY BATCH DONE $(date) ===="
