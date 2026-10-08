# dev/stress_replicate.R
# Replicated stress test: the single-draw stress ladder showed klue (clustering
# + EM) reaching a better LL than Apollo searchStart / gmnl on less-clear data.
# This replicates each rung over several seeds to show the pattern is SYSTEMATIC,
# not draw-specific. Reuses the validated method runners in compare_packages.R.
#
# RESUMABLE + METHOD-AWARE: results are saved to output/stress_replicate.rds
# after every (rung, seed) cell; a rerun tops up cells with missing arms.
# Gentle: single core. Apollo searchStart is the bottleneck; it runs at its
# 0.3.5 default settings (100 candidates).
#
# klue_ml is fitted with options(klue.screen = FALSE) (SCREEN_OFF and
# run_method() below). Its stored LLs come from klue's full-precision
# multistart before 0.9.1, which that option restores: with it every stored
# klue_ml LL reproduces exactly, while at klue's default 5 of the 20 cells end
# elsewhere (dev/check_ladder_carryover.R). The other arms keep the default.
#
# Seeds: cell data seed as.integer(1000*K + 100*(kap*100) + seed), seeds 1-5
# per rung (moderate 11500 + seed, hard and hard_corr 10000 + seed, very_hard
# 8000 + seed); klue_simulate calls set.seed with it, so a hard and a hard_corr
# cell with the same seed share their respondents' tastes and differ in the
# attributes and choices. klue's clustering starts use seed 123;
# apollo_searchStart draws its candidates after set.seed(17) (Apollo's default
# seed 13, plus 4), in both Apollo search arms. klue's EM, apollo_lcEM,
# apollo_estimate (BGW) and gmnl use no random numbers that affect the
# estimates. LLs reproduce; seconds do not.
#
# Run from the repository root, one core (as dev/rerun_standard_apollo.sh runs
# it):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   nohup nice -n 15 Rscript dev/stress_replicate.R \
#     > output/stress_replicate_std.log 2>&1 &
# dev/ladder_statistics.R reads that log for apollo_lcEM's iteration-cap
# messages.

sys.source("dev/compare_packages.R", envir = globalenv())  # functions only (guarded)
options(mc.cores = 1L)

DGP   <- klue_dgp(n_generic = 4, n_alternatives = 3)
SEEDS <- 1:5   # 5 seeds, as reported in the paper (Appendix tab:stress); the
                # blocked ladder (tab:stress_blocked) uses 10 seeds in its own
                # script. Extend to 1:10 only for a rerun: the script resumes
                # per (rung, seed) cell, so existing cells are kept.
RUNGS <- list(
  list(label = "moderate",  K = 4, kap = 0.75, sig = 0.25, corr = NULL),
  list(label = "hard",      K = 5, kap = 0.50, sig = 0.30, corr = NULL),
  list(label = "very_hard", K = 5, kap = 0.30, sig = 0.35, corr = NULL),
  list(label = "hard_corr", K = 5, kap = 0.50, sig = 0.30, corr = 0.6)
)
N_PER_CLASS <- 60L   # smaller N => harder + faster than the single-draw ladder
OUT <- "output/stress_replicate.rds"

stamp <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                   sprintf(...)))
`%||%` <- function(a, b) if (is.null(a)) b else a

# Ladder arms: name -> function(db, K) -> list(LL, converged, seconds). The
# Apollo arms follow Apollo's example scripts (dev/compare_packages.R):
# apollo_searchStart at its 0.3.5 default settings and apollo_lcEM as in
# EM_LC_no_covariates.r, both from the examples-style start (pooled MNL, class
# c scaled by 1/c). apollo_ss_published is the NON-default arm published in
# version12 (30 candidates, 3 stages, gTest = 10, smartStart = TRUE, generic
# start), kept so the published numbers can be reproduced.
METHODS <- list(
  klue_ml             = function(db, K) run_klue(db, K, DGP, "ml"),
  klue_em             = function(db, K) run_klue(db, K, DGP, "em"),
  apollo_searchStart  = function(db, K) run_apollo_searchStart(db, K, DGP),
  apollo_lcEM         = function(db, K) run_apollo_lcEM(db, K, DGP),
  gmnl_lc             = function(db, K) run_gmnl(db, K, DGP),
  apollo_ss_published = function(db, K) run_apollo_ss_published_randomised(db, K, DGP)
)
# klue_ml's stored LLs come from klue's pre-0.9.1 full-precision multistart, so
# the loop fits that arm with options(klue.screen = FALSE); with it every stored
# klue_ml LL reproduces exactly (dev/check_ladder_carryover.R). The other arms
# keep klue's default, under which klue_em reproduces.
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

results <- if (file.exists(OUT)) lapply(readRDS(OUT), migrate_cell) else list()
if (!dir.exists("output")) dir.create("output")

# Method-aware resume: a cell on disk is topped up with any arm missing from it
# (the data are regenerated bit-exactly from the stored seed formula); existing
# per-arm LLs are kept, and best/gap are recomputed over the union of arms.
for (r in RUNGS) {
  for (sd_i in SEEDS) {
    key  <- paste0(r$label, "|seed", sd_i)
    cell <- results[[key]]
    todo <- setdiff(names(METHODS), names(cell$LL))
    if (!length(todo)) { stamp("SKIP %s (all %d arms done)", key, length(METHODS)); next }
    stamp("RUN  %s  (K=%d kap=%.2f sig=%.2f corr=%s)  arms: %s", key, r$K, r$kap, r$sig,
          if (is.null(r$corr)) "none" else as.character(r$corr), paste(todo, collapse = ", "))
    seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + sd_i)
    d <- klue_simulate(N_per_class = N_PER_CLASS, T_tasks = 12, true_K = r$K,
                       separation = r$kap, heterogeneity = r$sig,
                       seed = seed_i, dgp = DGP, attr_corr = r$corr)
    db <- d$database
    fits <- lapply(setNames(todo, todo), run_method, db, r$K)
    lls  <- c(cell$LL,   vapply(fits, function(x) as.numeric(x$LL), numeric(1)))
    secs <- c(cell$secs, vapply(fits, function(x) as.numeric(x$seconds), numeric(1)))
    lls  <- lls[names(METHODS)]; secs <- secs[names(METHODS)]
    best <- suppressWarnings(max(lls[is.finite(lls)]))
    if (is.finite(cell$best %||% NA_real_) && best > cell$best + 1e-6)
      stamp("  NOTE %s: a new arm improved the best LL by %+.3f -- existing gaps shift",
            key, best - cell$best)
    results[[key]] <- list(rung = r$label, seed = sd_i, K = r$K, kap = r$kap,
                           sig = r$sig, corr = r$corr, best = best,
                           LL = lls, gap = lls - best, secs = secs,
                           harness = "apollo_std")
    saveRDS(results, OUT)
    stamp("  done %s : gaps  %s", key,
          paste(sprintf("%s=%+.2f", names(METHODS), lls - best), collapse = "  "))
  }
}

# ---- aggregate: mean gap and "stuck" rate (gap < -0.5) per rung x method ------
stamp("ALL CELLS DONE - aggregating")
meths <- names(METHODS)
cat("\n===== MEAN gap_to_best by rung x method (n seeds per cell) =====\n")
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
