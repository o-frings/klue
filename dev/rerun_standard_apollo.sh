#!/usr/bin/env bash
# =============================================================================
# Re-run every Apollo-dependent result of the paper with klue >= 0.10.0, whose
# Apollo code follows Apollo's example scripts (see klue/NEWS.md 0.10.0 and the
# "Apollo latent-class model code" block in dev/compare_packages.R).
#
# 1. Archives the superseded outputs to output/superseded_klue095/ (moved when a
#    driver would otherwise resume from them, copied when it tops up in place).
# 2. Runs every step in one single-core lane (default), or in two lanes with
#    --two-lanes. Each driver is seeded and resumable, so the command can be
#    rerun after an interruption (step 1 skips files already archived; the
#    MMNL benchmark skips fitted stages).
#
#    Memory: Apollo's analytic gradient keeps about one observations-by-draws
#    array per parameter; the Vittel correlated MMNL at 3000 draws took about
#    12 GB, and two lanes together overflowed the 32 GB machine on 2026-09-30.
#    Every step therefore runs with Apollo's memorySaver
#    (KLUE_MMNL_MEMORY_SAVER=TRUE), which gives bitwise-identical estimates
#    with about half the peak memory (klue/NEWS.md, 0.10.0).
#
#   Steps: MMNL benchmark Vittel; MMNL benchmark Swissmetro, Mode, SwissRoute,
#   Electricity; second start where the benchmark could decide the reading
#   (dev/run_mmnl_second_start.R) and the refreshed benchmark summary; H4 K*=1
#   blocked; H4 discrete blocked; H4 misspecification blocked; blocked stress
#   ladder; randomised stress ladder; randomised MMNL studies (tab:mmnl,
#   tab:mmnl_corr). The empirical benchmark runs before the H4 steps (order
#   changed 2026-10-01): its results move the paper's readings, and the order
#   does not change any number, since every step is seeded and runs alone.
#   A relaunch resumes: finished steps skip their fitted stages, and an
#   existing driver log is kept as <log>.<timestamp> instead of overwritten.
#
# Run from the repository root:  bash dev/rerun_standard_apollo.sh [--two-lanes]
# Logs: output/rerun_std_lane*.log (+ one log per driver, listed there).
# =============================================================================
set -u
cd "$(dirname "$0")/.."
ARCH=output/superseded_klue095
mkdir -p "$ARCH"
move() { for f in "$@"; do [ -e "$f" ] && [ ! -e "$ARCH/$(basename "$f")" ] && mv "$f" "$ARCH/"; done; }
copy() { for f in "$@"; do [ -e "$f" ] && [ ! -e "$ARCH/$(basename "$f")" ] && cp -p "$f" "$ARCH/"; done; }
# Drivers that resume from their CSV: move, so they recompute every cell.
# (h4_misspec_blocked.csv is new, with no superseded version: never moved, so a
# relaunch cannot archive the rerun's own cells.)
move output/k1_mmnl_blocked_full.csv output/discrete_mmnl_blocked.csv
# Drivers that overwrite or top up in place: keep a copy of the published file.
copy output/mmnl_results.csv output/mmnl_correlated_results.csv output/mmnl_benchmark_v12.csv \
     output/stress_replicate_blocked.rds output/stress_replicate.rds

ENVS="KLUE_CORES=1 KLUE_MMNL_MEMORY_SAVER=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1"
# run <driver log> <script> [args]: the lane log records START/EXIT, the driver
# log the R output.
run() {
  local log=$1; shift
  [ -e "$log" ] && mv "$log" "$log.$(date -r "$log" '+%m%d%H%M')"
  echo "[$(date '+%m-%d %H:%M')] START $* (log: $log)"
  env $ENVS nice -n 15 Rscript "$@" > "$log" 2>&1
  echo "[$(date '+%m-%d %H:%M')] EXIT $? $1"
}

laneA() {
  run output/mmnl_bench_std_vittel.log     dev/rerun_mmnl_benchmark.R Vittel
  run output/mmnl_bench_std_others.log     dev/rerun_mmnl_benchmark.R Swissmetro Mode SwissRoute Electricity
  run output/mmnl_second_start.log         dev/run_mmnl_second_start.R
  run output/mmnl_bench_std_summary.log    dev/rerun_mmnl_benchmark.R summary
  run output/k1_mmnl_blocked_std.log       dev/run_k1_mmnl_blocked.R
  run output/discrete_mmnl_blocked_std.log dev/run_discrete_mmnl_blocked.R
  run output/h4_misspec_blocked.log        dev/run_h4_misspec_blocked.R
}
laneB() {
  run output/stress_replicate_blocked_std.log dev/stress_replicate_blocked.R
  run output/stress_replicate_std.log         dev/stress_replicate.R
  run output/mmnl_studies_std.log             dev/run_mmnl_studies.R mmnl mmnl_corr
}

if [ "${1:-}" = "--two-lanes" ]; then
  laneA > output/rerun_std_laneA.log 2>&1 &
  PID_A=$!
  laneB > output/rerun_std_laneB.log 2>&1 &
  PID_B=$!
  echo "lanes started: A (pid $PID_A), B (pid $PID_B); see output/rerun_std_lane{A,B}.log"
else
  ( laneA; laneB ) >> output/rerun_std_lane.log 2>&1 &
  echo "single lane started (pid $!); see output/rerun_std_lane.log"
fi
