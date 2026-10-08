# dev/stress_replicate_blocked.R
# BLOCKED-baseline version of the package stress ladder (twin of the randomised
# dev/stress_replicate.R). Method runners from dev/compare_packages.R on
# progressively less-separable data drawn from ONE shared blocked D-efficient
# design (klue_design, internally seeded -> reproducible). The correlated rung
# is omitted: a fixed D-efficient design cannot carry an imposed attr_corr.
#
# METHODS (registry below; one entry per ladder arm):
#   klue_ml            clustering starts + direct ML          (proposed workflow)
#   klue_em            clustering starts + EM
#   klue_em_rp         random-partition starts + EM           (Stata lclogit /
#                      Latent GOLD default; added 2026-07-20 to test whether
#                      EM's start-robustness alone survives the hard rungs)
#   apollo_searchStart Apollo's apollo_searchStart at its 0.3.5 DEFAULT settings
#                      (100 candidates within +-0.1 of the start, smartStart =
#                      FALSE, 5 stages, gTest = 1e-3), then apollo_estimate
#                      (BGW), from the examples-style start (pooled MNL, class c
#                      scaled by 1/c); model code as in LC_no_covariates.r
#   apollo_lcEM        Apollo's apollo_lcEM as in EM_LC_no_covariates.r, same
#                      start
#   apollo_clust       Apollo's estimator from klue's six clustering starts,
#                      best kept: same starts as klue_ml, other optimiser
#   gmnl_lc            gmnl latent-class estimator
#   apollo_ss_published the NON-default arm published in version12
#                      (30 candidates, 3 stages, gTest = 10, candidates in
#                      [-3, 3], smartStart = FALSE, generic start; written by
#                      rerun_apollo_blocked.R, whose $secs were not updated),
#                      kept so the published numbers can be reproduced
#
# apollo_lcEM runs in Apollo 0.3.5 when written as Apollo's EM example writes
# it. The earlier "parameters do not influence the log-likelihood" failure came
# from combining the length(pi_values) loop with apollo_classAlloc, which in
# 0.3.5 gives zero analytic gradients for the class-specific parameters.
#
# RESUMABLE + METHOD-AWARE: saves the rds after each (rung,seed) cell. A cell
# already on disk is topped up with any methods missing from it (the data are
# regenerated bit-exactly from the stored seed formula and the seeded design),
# existing per-method LLs are KEPT, and best/gap are recomputed over the union
# of methods -- so published numbers can only change if a new method finds a
# strictly better optimum, which the log then flags.
# LOW-INTENSITY: single core (options mc.cores=1, apollo nCores=1); run with
# single-thread BLAS + nice. Apollo is the bottleneck.
#
# klue_ml is fitted with options(klue.screen = FALSE) (SCREEN_OFF and
# run_method() below). Its stored LLs come from klue's full-precision
# multistart before 0.9.1, which that option restores: with it every stored
# klue_ml LL reproduces exactly, while at klue's default they agree only to
# within 4e-5 (dev/check_ladder_carryover.R). The other arms keep the default.
#
# Seeds: the shared design klue_design(n_cards = 48, n_blocks = 4) runs under
# set.seed(20240601); cell data seed as.integer(1000*K + 100*(kap*100) + seed),
# seeds 1-10 per rung (moderate 11500 + seed, hard 10000 + seed, very_hard
# 8000 + seed; klue_simulate calls set.seed with it); klue's clustering starts
# use seed 123 (klue_ml, klue_em, apollo_clust); klue_em_rp draws partition p
# (1-6) under set.seed(cell_seed*100000 + 50000 + p); apollo_searchStart draws
# its candidates after set.seed(17) (Apollo's default seed 13, plus 4), in both
# Apollo search arms. klue's EM, apollo_lcEM, apollo_estimate (BGW) and gmnl
# use no random numbers that affect the estimates. LLs reproduce; seconds do
# not.
#
# Run from the repository root (gentle, detached; one core; as
# dev/rerun_standard_apollo.sh runs it):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   nohup nice -n 15 Rscript dev/stress_replicate_blocked.R \
#     > output/stress_replicate_blocked_std.log 2>&1 &
# Writes output/stress_replicate_blocked.rds (topped up in place; a backup of
# the published 5-method file is at output/stress_replicate_blocked_5methods.rds).
# dev/ladder_statistics.R reads the log for apollo_lcEM's iteration-cap
# messages.

sys.source("dev/compare_packages.R", envir = globalenv())  # functions only (guarded)
options(mc.cores = 1L)

DGP    <- klue_dgp(n_generic = 4, n_alternatives = 3)
DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)  # one shared design
SEEDS  <- 1:10   # 10 seeds/rung for tighter stress estimates (editor flagged thin replication)
RUNGS  <- list(
  list(label = "moderate",  K = 4, kap = 0.75, sig = 0.25),
  list(label = "hard",      K = 5, kap = 0.50, sig = 0.30),
  list(label = "very_hard", K = 5, kap = 0.30, sig = 0.35)
)
N_PER_CLASS <- 60L
OUT <- "output/stress_replicate_blocked.rds"

# Ladder arms: name -> function(db, K, cell_seed) -> list(LL, converged, seconds).
METHODS <- list(
  klue_ml             = function(db, K, cs) run_klue(db, K, DGP, "ml"),
  klue_em             = function(db, K, cs) run_klue(db, K, DGP, "em"),
  klue_em_rp          = function(db, K, cs) run_klue_em_rp(db, K, DGP, seed = cs),
  apollo_searchStart  = function(db, K, cs) run_apollo_searchStart(db, K, DGP),
  apollo_lcEM         = function(db, K, cs) run_apollo_lcEM(db, K, DGP),
  apollo_clust        = function(db, K, cs) run_apollo_from_clustering(db, K, DGP),
  gmnl_lc             = function(db, K, cs) run_gmnl(db, K, DGP),
  apollo_ss_published = function(db, K, cs) run_apollo_ss_published_blocked(db, K, DGP)
)
# klue_ml's stored LLs come from klue's pre-0.9.1 full-precision multistart, so
# the loop fits that arm with options(klue.screen = FALSE); with it every stored
# klue_ml LL reproduces exactly (dev/check_ladder_carryover.R). The other arms
# keep klue's default: klue_em reproduces under it, and klue_em_rp and the
# clustering starts of apollo_clust were fitted with it.
SCREEN_OFF <- "klue_ml"
run_method <- function(m, ...) {
  if (m %in% SCREEN_OFF) { op <- options(klue.screen = FALSE); on.exit(options(op)) }
  METHODS[[m]](...)
}
# Apollo arms follow Apollo's example scripts since 2026-09-29 (see the
# "Apollo latent-class model code" block in dev/compare_packages.R). Cells
# written before hold the published non-default Apollo arm under the name
# apollo_searchStart: rename it apollo_ss_published, so the default-settings
# arm is computed fresh, and drop arms that are no longer run.
migrate_cell <- function(cell) {
  if (is.null(cell) || identical(cell$harness, "apollo_std")) return(cell)
  rn <- function(v) {
    if (!is.null(v)) names(v)[names(v) == "apollo_searchStart"] <- "apollo_ss_published"
    v
  }
  LL <- rn(cell$LL); secs <- rn(cell$secs)
  keep <- names(LL) %in% names(METHODS) &
          (is.finite(LL) | names(LL) == "apollo_ss_published")  # redo errored arms only
  cell$LL   <- LL[keep]
  cell$secs <- secs[names(cell$LL)]
  # best and gap follow the kept arms; the "apollo_std" tag is set only once
  # the cell has been topped up with the new arms (below), so an interrupted
  # run never leaves a gap under a name whose arm has not been computed.
  cell$best <- suppressWarnings(max(cell$LL[is.finite(cell$LL)]))
  cell$gap  <- cell$LL - cell$best
  cell$harness <- NULL
  cell
}


stamp <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   sprintf(...)))
`%||%` <- function(a, b) if (is.null(a)) b else a

results <- if (file.exists(OUT)) lapply(readRDS(OUT), migrate_cell) else list()
if (!dir.exists("output")) dir.create("output")

# The seed loop variable is not called `s`: a global `s` only makes Apollo warn
# that a global object is used inside apollo_probabilities (the class loop
# index), but the warnings clutter the log.
for (r in RUNGS) {
  for (sd_i in SEEDS) {
    key  <- paste0(r$label, "|seed", sd_i)
    cell <- results[[key]]
    todo <- setdiff(names(METHODS), names(cell$LL))
    if (!length(todo)) { stamp("SKIP %s (all %d methods done)", key,
                               length(METHODS)); next }
    stamp("RUN  %s  (K=%d kap=%.2f sig=%.2f, blocked)  methods: %s",
          key, r$K, r$kap, r$sig, paste(todo, collapse = ", "))
    seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + sd_i)
    d <- klue_simulate(N_per_class = N_PER_CLASS, true_K = r$K,
                       separation = r$kap, heterogeneity = r$sig,
                       seed = seed_i, dgp = DGP, design = DESIGN)
    db <- d$database

    fits <- lapply(setNames(todo, todo), run_method, db, r$K, seed_i)
    LL   <- c(cell$LL,   vapply(fits, function(x) as.numeric(x$LL), numeric(1)))
    secs <- c(cell$secs, vapply(fits, function(x) as.numeric(x$seconds), numeric(1)))
    LL   <- LL[names(METHODS)]; secs <- secs[names(METHODS)]

    old_best <- cell$best
    best <- suppressWarnings(max(LL[is.finite(LL)]))
    if (is.finite(old_best %||% NA_real_) && best > old_best + 1e-6)
      stamp("  NOTE %s: new method improved the best LL by %+.3f -- existing gaps shift",
            key, best - old_best)

    upd <- if (is.null(cell)) list(rung = r$label, seed = sd_i, K = r$K,
                                   kap = r$kap, sig = r$sig, corr = NULL) else cell
    upd$best <- best; upd$LL <- LL; upd$gap <- LL - best; upd$secs <- secs
    upd$harness <- "apollo_std"
    results[[key]] <- upd
    saveRDS(results, OUT)
    stamp("  done %s : gaps  %s", key,
          paste(sprintf("%s=%+.2f", names(METHODS), upd$gap), collapse = "  "))
  }
}

# ---- aggregate: mean gap and "stuck" rate (gap < -0.5) per rung x method ------
stamp("ALL CELLS DONE - aggregating")
meths <- names(METHODS)
cat("\n===== BLOCKED: MEAN gap_to_best by rung x method (n seeds per cell) =====\n")
cat(sprintf("%-10s", "rung"))
for (m in meths) cat(sprintf(" | %18s", m))
cat("\n")
for (r in RUNGS) {
  cells <- results[grepl(paste0("^", r$label, "\\|"), names(results))]
  cat(sprintf("%-10s", r$label))
  for (m in meths) {
    g <- sapply(cells, function(c) c$gap[m])
    g <- g[is.finite(g)]
    cat(sprintf(" | %6.2f (stuck %d/%d)", if (length(g)) mean(g) else NA_real_,
                sum(g < -0.5), length(g)))
  }
  cat("\n")
}
cat("\n(mean LL below best; 'stuck k/n' = cells >0.5 LL below best of n seeds)\n")
