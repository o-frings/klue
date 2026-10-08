#!/usr/bin/env Rscript
# Three-arm initialisation ablation (RP / one-hot / random) on the H1 design,
# randomised reference design; the blocked-design run is stage 3 of
# dev/run_blocked_baseline.R. Feeds the H1a numbers in the manuscript.
# Output: output/init_ablation.rds + output/init_ablation_summary.txt
#         (+ the marker output/init_ablation.done)
# Seeds: condition data under .cond_seed(K, kappa, sigma, rep) =
#   1000*K + 100*(kappa*100) + 10*(sigma*100) + rep (studies/klue_studies.R),
#   which klue_simulate() passes to set.seed(); random start r under
#   set.seed(seed*1000 + r); klue's clustering starts (RP and one-hot arms)
#   under seed 123. The output is the same at any KLUE_CORES.
# Run from the repository root (KLUE_CORES = conditions fitted in parallel):
#   KLUE_CORES=4 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/run_init_ablation.R
suppressMessages(pkgload::load_all("klue", quiet = TRUE))
source("studies/klue_studies.R")   # the klue_study_* drivers (moved out of the package 2026-09-17)
n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else klue:::.klue_cores()
options(klue.cores = n_cores)      # condition-level parallelism of the drivers
options(klue.screen = FALSE)    # full-precision multistart, as when the published results were
                                # produced (before klue 0.9.1); README_REPRODUCE.md, "Package history"

cat("=== INIT ABLATION runner ===\n")
cat("R:", R.version.string, "\n")
cat("Cores:", parallel::detectCores(), "  klue.cores:", n_cores, "\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

t0 <- Sys.time()
df <- klue_study_initialisation(n_random = 50, n_cond = 40, verbose = TRUE)
t1 <- Sys.time()

cat("\nEnd:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Elapsed:", format(t1 - t0), "\n")

dir.create("output", showWarnings = FALSE)
saveRDS(df, "output/init_ablation.rds")

sink("output/init_ablation_summary.txt")
cat("=== INITIALISATION ABLATION RESULTS ===\n\n")
cat("Conditions evaluated:", nrow(df), "\n\n")
cat("Overall:\n")
cat(sprintf("  RP contrasts reach global:    %.1f%% (%d/%d)\n",
            100 * mean(df$rp_at_global), sum(df$rp_at_global), nrow(df)))
cat(sprintf("  One-hot indicators at global: %.1f%% (%d/%d)\n",
            100 * mean(df$oh_at_global), sum(df$oh_at_global), nrow(df)))
cat(sprintf("  Random starts at global:      %.1f%% per start\n",
            100 * mean(df$random_pct_global)))
cat("\nBy true K:\n")
print(aggregate(cbind(rp_at_global, oh_at_global, random_pct_global) ~ K,
                data = df, FUN = mean), row.names = FALSE)
cat("\nLog-likelihood gap (RP minus one-hot), per condition:\n")
print(summary(df$gap_rp_oh))
sink()

cat("\nWritten output/init_ablation.rds and output/init_ablation_summary.txt\n")
file.create("output/init_ablation.done")
