#!/usr/bin/env Rscript
# Move non-converged MMNL benchmark stages aside (output/mmnl_bench_std/failed/)
# so that dev/rerun_mmnl_benchmark.R refits them. Only useful after a code
# change: the fits are deterministic, so an unchanged refit fails again. A
# stage that Apollo reports as unsuccessful is a result, not a crash.
# Run: Rscript dev/requeue_failed_mmnl.R Swissmetro Mode SwissRoute Electricity
ds <- commandArgs(trailingOnly = TRUE)
dir.create("output/mmnl_bench_std/failed", showWarnings = FALSE)
for (f in list.files("output/mmnl_bench_std", pattern = sprintf("^(%s)_M[1-4]\\.rds$",
                     paste(ds, collapse = "|")), full.names = TRUE)) {
  if (!isTRUE(readRDS(f)$converged)) {
    file.rename(f, file.path("output/mmnl_bench_std/failed", basename(f)))
    cat("requeued", basename(f), "\n")
  }
}
