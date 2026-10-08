#!/usr/bin/env Rscript
# =============================================================================
# Vittel re-run under the IDENTIFIED specification (2026-09-21).
#
# The published Vittel numbers entered both levels of the two-level forest
# attribute alongside the status-quo ASC, so Forest_For_Water +
# Forest_For_Biodiv + ASC_sq was constant across alternatives and the model was
# rank-deficient by one parameter per class (see the identification note in
# R/empirical_application.R). The likelihood is unaffected, but k counted a
# parameter the design cannot identify, over-penalising BIC/AIC/ICL by
# C*log(N), and the three affected coefficients were identified only up to a
# constant. This script re-runs the whole Vittel arm with the forest attribute
# coded against a reference level (forest_ref = "biodiv"), so LCMNL and both
# MMNL benchmarks are penalised consistently.
#
# Stages, each resumable (a stage whose .rds exists is skipped; delete to redo):
#   1 lcmnl  C = 1..C_MAX, extended until BIC turns  -> vittel_respec_lcmnl.rds/.csv
#   2 mmnl   independent normals, lognormal price    -> vittel_respec_mmnl.rds
#   3 corr   full-Cholesky correlated MMNL           -> vittel_respec_mmnl_corr.rds
#   4 summary: model comparison + old-vs-new table   -> vittel_respec_comparison.csv
#
# Deterministic: the LCMNL multistart is seeded per C and the MMNL draws are
# fixed MLHS. Writes only vittel_respec_* files; the published artifacts are
# left untouched for comparison.
#
# Launch (2 cores, the machine's cap):
#   KLUE_CORES=2 nohup nice -n 10 Rscript dev/rerun_vittel_respec.R \
#     > output/rerun_vittel_respec.log 2>&1 &
# =============================================================================

n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else 2L
# Apollo's own cluster is kept to a single core: it can deadlock when detached
# at > 1 core (README_REPRODUCE.md, "Notes on reproducibility"), and a 2-core
# detached run of the correlated fit died here with a bus error ('invalid
# alignment') on 2026-09-21. The LCMNL multistart still uses n_cores.
options(klue.cores = n_cores, klue.mmnl.n_cores = 1L, mc.cores = n_cores,
        empirical.skip = TRUE)                     # load the adapter, skip its autorun
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(source("R/empirical_application.R")))

OUT    <- "output"
C_MAX  <- 12L          # safety cap; the search stops as soon as BIC turns
N_DRAWS <- 3000L
# MMNL starting values. "zero" is the conventional neutral start (attribute
# means at 0, lognormal price coefficient at -0.01, spreads at 0.1); "sign"
# starts the means at +0.01 instead; "informed" is klue's warm start. The first
# pass used "informed", whose correlated arm fell back to an arbitrary start
# (all means +0.5, b_price = -1) when the internal warm-start fit failed, so the
# neutral start is the one to trust. Stage files are suffixed by the scheme, so
# runs under different schemes sit side by side. Since klue 0.10.0 (standard
# Apollo estimation: sigma unconstrained, one estimation, correlated model
# started from the independent estimates) the MMNL stage files also carry
# "_std", so the superseded 0.9.x fits are never reused.
START <- Sys.getenv("KLUE_MMNL_START", "zero")
sfx <- if (identical(START, "informed")) "" else paste0("_", START)
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                                     paste0(...))); flush.console() }
stage <- function(name, fn) {
  f <- file.path(OUT, sprintf("vittel_respec_%s.rds", name))
  if (file.exists(f)) { stamp("SKIP ", name, " (exists)"); return(readRDS(f)) }
  stamp("START ", name)
  t <- system.time(obj <- fn())[3]
  saveRDS(obj, f)
  stamp(sprintf("DONE  %s (%.1f min)", name, t / 60))
  obj
}

db <- load_empirical_database()                    # forest_ref = "biodiv"
stamp(sprintf("Vittel: N=%d, T=%d, %d generic attributes + price, cores=%d, MMNL start=%s",
              length(unique(db$ID)), max(db$TASK), attr(db, "n_generic"), n_cores, START))

# ---- 1. LCMNL ladder, extended until BIC turns ------------------------------
lcmnl <- stage("lcmnl", function() {
  rows <- list(); fits <- list(); best <- Inf
  for (C in 1:C_MAX) {
    t0 <- Sys.time()
    m  <- klue_lcmnl(db, C, dgp = EMP_DGP)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    fits[[as.character(C)]] <- m
    rows[[length(rows) + 1]] <- data.frame(
      C = C, converged = m$converged, LL = m$LL, k = m$k, BIC = m$BIC,
      AIC = m$AIC, ICL = m$ICL, ICL_BIC = m$ICL_BIC,
      min_class_prob = min(m$class_probs), singular = isTRUE(m$singular),
      best_method = if (is.null(m$best_method)) "pooled" else m$best_method,
      secs = round(dt, 1), stringsAsFactors = FALSE)
    stamp(sprintf("  C=%2d  LL=%.2f  k=%d  BIC=%.1f%s  singular=%s  (%.0fs)",
                  C, m$LL, m$k, m$BIC, if (m$BIC < best) " *new min*" else "",
                  isTRUE(m$singular), dt))
    if (m$BIC < best) best <- m$BIC else { stamp("  BIC turned; stopping the ladder"); break }
  }
  list(table = do.call(rbind, rows), fits = fits)
})
write.csv(lcmnl$table, file.path(OUT, "vittel_respec_lcmnl.csv"), row.names = FALSE)

# ---- 2 & 3. MMNL benchmarks --------------------------------------------------
mmnl <- stage(paste0("mmnl", sfx, "_std"), function()
  klue_mmnl(db, n_draws = N_DRAWS, start = START, dgp = EMP_DGP))
stamp(sprintf("  MMNL  LL=%.2f k=%d BIC=%.1f conv=%s", mmnl$LL, mmnl$k, mmnl$BIC, mmnl$converged))

corr <- stage(paste0("mmnl_corr", sfx, "_std"), function()
  klue_mmnl(db, correlation = TRUE, n_draws = N_DRAWS, start = START, dgp = EMP_DGP,
            warm_start = if (isTRUE(mmnl$converged)) mmnl))
stamp(sprintf("  MMNL-corr  LL=%.2f k=%d BIC=%.1f conv=%s", corr$LL, corr$k, corr$BIC, corr$converged))

# ---- 4. Comparison ----------------------------------------------------------
tab   <- lcmnl$table
bestC <- tab$C[which.min(tab$BIC)]
cmp <- data.frame(
  model = c("MNL (C=1)", sprintf("LCMNL (C=%d, BIC-best)", bestC),
            "MMNL (independent normals)", "MMNL (correlated, full Cholesky)"),
  LL  = c(tab$LL[tab$C == 1],  tab$LL[tab$C == bestC],  mmnl$LL,  corr$LL),
  k   = c(tab$k[tab$C == 1],   tab$k[tab$C == bestC],   mmnl$k,   corr$k),
  BIC = c(tab$BIC[tab$C == 1], tab$BIC[tab$C == bestC], mmnl$BIC, corr$BIC),
  stringsAsFactors = FALSE)
cmp$dBIC <- cmp$BIC - min(cmp$BIC)
write.csv(cmp, file.path(OUT, sprintf("vittel_respec_comparison%s.csv", sfx)), row.names = FALSE)

cat("\n===== VITTEL, IDENTIFIED SPECIFICATION =====\n"); print(tab, row.names = FALSE)
cat("\n"); print(cmp, row.names = FALSE)
cat(sprintf("\nBIC-best: %s  (verdict: %s)\n",
            cmp$model[which.min(cmp$BIC)],
            if (grepl("LCMNL", cmp$model[which.min(cmp$BIC)])) "discrete" else "continuous"))

# Side-by-side against the published (rank-deficient) numbers, where available.
pub <- file.path(OUT, "empirical_results.csv")
if (file.exists(pub)) {
  o <- read.csv(pub)
  m <- merge(o[, c("C", "LL", "k", "BIC", "ICL")], tab[, c("C", "LL", "k", "BIC", "ICL")],
             by = "C", suffixes = c("_published", "_identified"))
  m$dLL <- m$LL_identified - m$LL_published
  m$dBIC <- m$BIC_identified - m$BIC_published
  cat("\n--- published vs identified (dLL should be ~0 at every C) ---\n")
  print(m[, c("C", "LL_published", "LL_identified", "dLL",
              "k_published", "k_identified", "BIC_published", "BIC_identified", "dBIC")],
        row.names = FALSE)
  write.csv(m, file.path(OUT, sprintf("vittel_respec_vs_published%s.csv", sfx)), row.names = FALSE)
}
stamp("ALL STAGES COMPLETE")
