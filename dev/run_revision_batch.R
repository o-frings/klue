# dev/run_revision_batch.R
# One-shot heavy batch for the v9 revision. Runs every compute-bearing piece
# sequentially, saves each result to output/ as .rds, and logs progress. Each
# stage is wrapped so one failure does not abort the rest.
#
# Launch (background, overnight) from repo root:
#   nohup Rscript dev/run_revision_batch.R > output/revision_batch.log 2>&1 &
#
# Optional faster pass (fewer starts / draws) for a same-day sanity batch:
#   Rscript dev/run_revision_batch.R quick

# From the repo, load the dev source (latest code); else the installed package.
if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
source("studies/klue_studies.R")   # the klue_study_* drivers (moved out of the package 2026-09-17)
} else {
  suppressWarnings(suppressMessages(library(klue)))
}
if (!dir.exists("output")) dir.create("output")

QUICK <- "quick" %in% commandArgs(trailingOnly = TRUE)
N_START <- if (QUICK) 20L else 30L      # diffuse-random / perturbation budget
N_DRAWS <- if (QUICK) 500L else 3000L   # MMNL MLHS draws
N_COND  <- 36L                          # full convergence design (3 K x 3 kappa x 2 sigma x 2 rep)

# Thermal cap: run "gentle" on a laptop in active use. Set BOTH the LCMNL
# condition-level parallelism and Apollo's within-model parallelism low so the
# machine stays cool and usable. Override with env var KLUE_CORES.
n_use <- as.integer(Sys.getenv("KLUE_CORES", unset = "2"))
n_use <- max(1L, min(n_use, parallel::detectCores() - 1L))
options(klue.cores = n_use)   # the drivers read .klue_cores(); the old N_CORES_LCMNL namespace poke was a silent no-op
options(klue.mmnl.n_cores = n_use)   # cap Apollo (MMNL) cores too

stamp <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   sprintf(...)))
stage <- function(label, file, expr) {
  # Resumable: skip a stage whose result already exists (delete the .rds to redo).
  if (file.exists(file.path("output", file))) {
    stamp("SKIP   %s (output/%s already exists)", label, file)
    return(readRDS(file.path("output", file)))
  }
  stamp("START  %s", label)
  t0 <- Sys.time()
  res <- tryCatch(force(expr), error = function(e) {
    stamp("FAILED %s : %s", label, conditionMessage(e)); NULL
  })
  if (!is.null(res)) {
    saveRDS(res, file.path("output", file))
    stamp("DONE   %s -> output/%s (%.1f min)", label, file,
          as.numeric(difftime(Sys.time(), t0, units = "mins")))
  }
  res
}

stamp("revision batch start  (quick=%s, n_start=%d, n_draws=%d, cores=%d)",
      QUICK, N_START, N_DRAWS, n_use)

# --- WS1: estimator x start-strategy, orthogonal design ----------------------
stage("estimator sweep [orthogonal]", "estimator_orthogonal.rds",
      klue_study_estimator(n_random = N_START, n_perturb = N_START,
                           n_cond = N_COND, verbose = TRUE))

# --- WS1+WS2: same sweep under correlated attributes -------------------------
stage("estimator sweep [attr_corr=0.6]", "estimator_corr.rds",
      klue_study_estimator(n_random = N_START, n_perturb = N_START,
                           n_cond = N_COND, attr_corr = 0.6, verbose = TRUE))

# --- WS2: contrasts-vs-one-hot ablation, orthogonal vs correlated ------------
stage("init ablation [orthogonal]", "init_ablation_orthogonal.rds",
      klue_study_initialisation(n_random = N_START, n_cond = N_COND,
                                verbose = TRUE))
stage("init ablation [attr_corr=0.6]", "init_ablation_corr.rds",
      klue_study_initialisation(n_random = N_START, n_cond = N_COND,
                                attr_corr = 0.6, verbose = TRUE))

# --- WS3: correlated vs independent MMNL on the 5 public datasets ------------
# (sourced runner returns its data.frame; see dev/empirical_corr_mmnl.R)
stage("empirical corr-MMNL [5 datasets]", "empirical_corr_mmnl.rds", {
  sys.source("dev/empirical_corr_mmnl.R", envir = globalenv())  # defines run_one, SPECS
  do.call(rbind, lapply(names(SPECS), run_one, n_draws = N_DRAWS))
})

stamp("revision batch complete")
