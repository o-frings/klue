# dev/run_design_discrete.R
# Re-run the robustness experiment "Discrete versus random attribute levels"
# (paper Appendix, Table tab:design) with klue_study_design(), whose discrete
# arm now draws each attribute independently from five levels
# (klue_simulate_discrete(), 2026-10-19). The earlier discrete arm, the cyclic
# grid of klue_simulate_deff(), made the price level collinear with the
# quality levels before its jitter.
# Run from the repository root (KLUE_CORES = conditions fitted at once):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   KLUE_CORES=2 nice -n 10 Rscript dev/run_design_discrete.R
# Writes: output/design_results.csv (one row per condition) and
#         output/design_results.rds.
# Grid: K* in {2, 3} x kappa in {0.75, 1.0} x 5 replications, sigma = 0.2,
# N = 150 per class, T = 20; LCMNLs with C = 1..5, the true C chosen by BIC.
# Seeds: data seed per condition as.integer(4000 + 100*K + 10*(kappa*100) + rep),
# shared by the continuous and discrete arms, so each pair has the same
# respondents and coefficients; the discrete levels continue the same RNG
# stream. klue's clustering starts use seed 123. Same output at any KLUE_CORES.
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("studies/klue_studies.R")
OUT <- "output"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else klue:::.klue_cores()
options(klue.cores = n_cores)
options(klue.screen = FALSE)    # full-precision multistart, as for the other randomised-design
                                # studies (README_REPRODUCE.md, "Package history")
t0 <- Sys.time()
res <- klue_study_design(verbose = TRUE, n_cores = n_cores)
saveRDS(res, file.path(OUT, "design_results.rds"))
write.csv(res, file.path(OUT, "design_results.csv"), row.names = FALSE)
print(res)
cat(sprintf("Done in %.1f min, cores = %d\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), n_cores))
