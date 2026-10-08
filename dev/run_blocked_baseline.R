# =============================================================================
# Blocked-baseline LCMNL study (NO MMNL).
# Re-runs the headline enumeration + H1/H1a/H1b initialisation benchmarks under
# the harder realistic-DCE baseline (shared blocked D-efficient design, T=12).
# Resumable: each stage skips if its .rds output already exists (delete to redo).
# Launch from the repository root (KLUE_CORES = conditions fitted at once):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   KLUE_CORES=8 nohup nice -n 10 Rscript dev/run_blocked_baseline.R \
#     > /tmp/blocked_run.log 2>&1 &
# Writes: output/{main,convergence,init_ablation,estimator}_blocked.{rds,csv}.
#   Stage 4 is the computation dev/run_estimator_blocked.R runs in chunks;
#   whichever runs first writes output/estimator_blocked.{rds,csv}, and the
#   other then skips it.
# Seeds: per-condition data seed .cond_seed() from studies/klue_studies.R,
#   as.integer(1000*K + 100*(kappa*100) + 10*(sigma*100) + rep); the shared
#   design klue_design(n_cards = 48, n_blocks = 4) under set.seed(20240601);
#   klue's clustering starts under seed 123; diffuse random start r (stages
#   2-4) under set.seed(seed*1000 + r); in stage 4 also perturbation start s
#   under set.seed(seed*100000 + s) and random partition s under
#   set.seed(seed*100000 + 50000 + s). Same output at any KLUE_CORES.
# =============================================================================
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("studies/klue_studies.R")   # the klue_study_* drivers (moved out of the package 2026-09-17)

OUT <- "output"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else klue:::.klue_cores()
options(klue.cores = n_cores)   # the drivers read .klue_cores(); the old N_CORES_LCMNL namespace poke was a silent no-op
options(klue.screen = FALSE)    # full-precision multistart, as when the published results were
                                # produced (before klue 0.9.1); README_REPRODUCE.md, "Package history"

stamp <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }
save_stage <- function(obj, name) {
  saveRDS(obj, file.path(OUT, paste0(name, ".rds")))
  if (is.data.frame(obj))
    tryCatch(write.csv(obj, file.path(OUT, paste0(name, ".csv")), row.names = FALSE),
             error = function(e) NULL)
}
run_stage <- function(name, fn) {
  rds <- file.path(OUT, paste0(name, ".rds"))
  if (file.exists(rds)) { stamp(sprintf("SKIP %s (exists)", name)); return(invisible()) }
  stamp(sprintf("START %s", name))
  t <- system.time(obj <- fn())[3]
  save_stage(obj, name)
  stamp(sprintf("DONE  %s (%.1f min) -> %s.rds", name, t / 60, name))
}

stamp(sprintf("Blocked-baseline LCMNL study (no MMNL), cores=%d", n_cores))

# 1. Headline: enumeration accuracy under the blocked baseline (n_reps = 20).
run_stage("main_blocked",
          function() klue_study_main(blocked = TRUE, n_reps = 20, verbose = TRUE))

# 2. H1: clustering vs uninformed-random convergence, blocked baseline.
run_stage("convergence_blocked",
          function() klue_study_convergence(blocked = TRUE, verbose = TRUE))

# 3. H1a: RP-contrast vs one-hot feature ablation, blocked baseline.
run_stage("init_ablation_blocked",
          function() klue_study_initialisation(blocked = TRUE, verbose = TRUE))

# 4. H1b: {ML,EM} x {random, perturbation, random partition, clustering},
#    blocked baseline. Heaviest stage; starts reduced 50 -> 20 (per-start rates
#    already clean at 20).
run_stage("estimator_blocked",
          function() klue_study_estimator(blocked = TRUE, n_random = 20L,
                                          n_perturb = 20L, verbose = TRUE))

stamp("ALL BLOCKED LCMNL STAGES COMPLETE")
