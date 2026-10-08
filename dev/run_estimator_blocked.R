# =============================================================================
# H1b: {ML,EM} x {clustering, perturbation, random partition, random} start
# factorial under the BLOCKED realistic-DCE baseline.  CHECKPOINTED in chunks of
# 6 conditions so an interruption costs at most one chunk.  Needs NO internet
# (klue ML/EM + local idefix design).  Resumable: skips any chunk whose
# part-file exists, then combines into output/estimator_blocked.{rds,csv}.
# Writes: output/estimator_blocked_chunk01.rds to _chunk06.rds (the part-files)
#   and output/estimator_blocked.{rds,csv}. Quits at once if
#   output/estimator_blocked.rds exists (stage 4 of dev/run_blocked_baseline.R
#   writes the same file).
# Seeds: per-condition data seed .cond_seed() from studies/klue_studies.R,
#   as.integer(1000*K + 100*(kappa*100) + 10*(sigma*100) + rep); the shared
#   design klue_design(n_cards = 48, n_blocks = 4) under set.seed(20240601);
#   klue's clustering starts under seed 123; the other three strategies get 20
#   starts each: perturbation start s under set.seed(seed*100000 + s), random
#   partition s under set.seed(seed*100000 + 50000 + s), diffuse random start r
#   under set.seed(seed*1000 + r). Same output at any KLUE_CORES or chunking.
# Launch from the repository root (detached, survives disconnection;
# KLUE_CORES = conditions fitted at once):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   KLUE_CORES=8 nohup nice -n 10 Rscript dev/run_estimator_blocked.R \
#     > /tmp/estimator_blocked.log 2>&1 &
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

final <- file.path(OUT, "estimator_blocked.rds")
if (file.exists(final)) { stamp("estimator_blocked.rds already exists - nothing to do"); quit(save = "no") }

chunks <- split(1:36, ceiling((1:36) / 6))   # 6 chunks x 6 conditions
stamp(sprintf("H1b estimator (blocked), %d chunks of 6 conditions, cores=%d", length(chunks), n_cores))

for (ci in seq_along(chunks)) {
  pf <- file.path(OUT, sprintf("estimator_blocked_chunk%02d.rds", ci))
  if (file.exists(pf)) { stamp(sprintf("SKIP chunk %d (exists)", ci)); next }
  stamp(sprintf("START chunk %d/%d (conditions %s)", ci, length(chunks),
                paste(range(chunks[[ci]]), collapse = "-")))
  t <- system.time(
    d <- klue_study_estimator(blocked = TRUE, n_random = 20L, n_perturb = 20L,
                              cond_idx = chunks[[ci]], verbose = TRUE)
  )[3]
  saveRDS(d, pf)
  stamp(sprintf("DONE chunk %d (%.1f min) -> %s", ci, t / 60, basename(pf)))
}

# Combine all chunks
parts <- lapply(seq_along(chunks),
                function(ci) readRDS(file.path(OUT, sprintf("estimator_blocked_chunk%02d.rds", ci))))
df <- do.call(rbind, parts[!vapply(parts, is.null, logical(1))])
saveRDS(df, final)
write.csv(df, file.path(OUT, "estimator_blocked.csv"), row.names = FALSE)
stamp(sprintf("ALL DONE -> estimator_blocked.{rds,csv}  (%d rows)", nrow(df)))

# Headline summary: per-start at-global rate by estimator x strategy
agg <- aggregate(cbind(at_global_rate, reached_global) ~ estimator + strategy, data = df, FUN = mean)
agg <- agg[order(agg$estimator, -agg$at_global_rate), ]
cat("\n=== H1b BLOCKED: per-start at-global / reached-global ===\n")
for (k in seq_len(nrow(agg)))
  cat(sprintf("  %-3s x %-12s  per-start %5.1f%%   reached %5.1f%%\n",
              agg$estimator[k], agg$strategy[k],
              100 * agg$at_global_rate[k], 100 * agg$reached_global[k]))
