# dev/run_robustness_blocked.R
# Re-run the robustness analyses under the blocked D-efficient baseline
# (present-both with the orthogonal versions). LCMNL-only, multicore-safe.
# Arms: unbalanced class sizes, concomitant (covariate membership), sample
# size / panel length. The random-vs-D-efficient arm is NOT re-run (the blocked
# baseline IS a D-efficient design; it is reframed in the manuscript instead).
# Resumable: skips any stage whose .rds exists. Run from the repository root
# (KLUE_CORES = conditions fitted at once):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   KLUE_CORES=8 nohup nice -n 10 Rscript dev/run_robustness_blocked.R \
#     > /tmp/robustness_blocked.log 2>&1 &
# Writes: output/{unbalanced,concomitant,sample}_blocked.{rds,csv}.
# Seeds: these drivers (studies/klue_studies.R) do not use .cond_seed(). Data
# seed per condition, with ci = 1-4 for the balanced, mild, moderate and severe
# shares, cs the covariate strength and N the respondents per class:
#   unbalanced   as.integer(3000 + 1000*ci + 100*(kappa*100) + rep)
#   concomitant  as.integer(5000 + 100*K + 10*(kappa*100) + rep + 1000*cs)
#   sample       as.integer(7000 + 100*K + 10*(kappa*100) + rep + 100*T + N)
# The concomitant and sample formulas (from klue 0.6.x) give 12 and 24 pairs of
# conditions the same seed; the paired conditions differ, so their data do too.
# The design is klue_design(n_cards = 48, n_blocks = 4) under
# set.seed(20240601); the sample arm builds one per panel length (32, 48 or 80
# cards in 4 blocks), each under that seed. klue's clustering starts use seed
# 123. Same output at any KLUE_CORES.
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("studies/klue_studies.R")   # the klue_study_* drivers (moved out of the package 2026-09-17)
OUT <- "output"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else klue:::.klue_cores()
options(klue.cores = n_cores)   # the drivers read .klue_cores(); the old N_CORES_LCMNL namespace poke was a silent no-op
options(klue.screen = FALSE)    # full-precision multistart, as when the published results were
                                # produced (before klue 0.9.1); README_REPRODUCE.md, "Package history"
stamp <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }

run_stage <- function(name, fn) {
  rds <- file.path(OUT, paste0(name, ".rds"))
  if (file.exists(rds)) { stamp(sprintf("SKIP %s (exists)", name)); return(invisible()) }
  stamp(sprintf("START %s", name))
  t <- system.time(obj <- fn())[3]
  saveRDS(obj, rds)
  if (is.data.frame(obj)) tryCatch(write.csv(obj, file.path(OUT, paste0(name, ".csv")),
                                             row.names = FALSE), error = function(e) NULL)
  stamp(sprintf("DONE  %s (%.1f min)", name, t / 60)); print(obj)
}

stamp(sprintf("Robustness under blocked baseline, cores=%d", n_cores))
run_stage("unbalanced_blocked",  function() klue_study_unbalanced(blocked = TRUE, verbose = TRUE))
run_stage("concomitant_blocked", function() klue_study_concomitant(blocked = TRUE, verbose = TRUE))
run_stage("sample_blocked",      function() klue_study_sample(blocked = TRUE, verbose = TRUE))
stamp("ALL ROBUSTNESS-BLOCKED STAGES DONE")
