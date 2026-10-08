# dev/rerun_apollo_blocked.R
# SUPERSEDED (2026-09-29): wrote the non-default Apollo arm published in
# version12, now kept in the ladder as apollo_ss_published; the ladder's
# apollo_searchStart arm runs Apollo's defaults (dev/stress_replicate_blocked.R).
# Re-run ONLY the Apollo searchStart arm of the blocked stress ladder with the
# robust setting (smartStart=FALSE + bounded candidates), because the default
# smartStart's Hessian step produced NA probabilities on the low-diversity
# blocked cards. klue_ml/klue_em/gmnl results are kept as-is; for each cell we
# regenerate the identical data (same seed), re-estimate Apollo, overwrite
# cell$LL["apollo_searchStart"], and recompute cell$gap = LL - max(finite LL).
# Resumable: cells flagged cell$apollo_robust=TRUE are skipped.
# Low-intensity: single core; run with single-thread BLAS + nice.

stop("dev/rerun_apollo_blocked.R is superseded: the published non-default Apollo arm is ",
     "apollo_ss_published in dev/stress_replicate_blocked.R", call. = FALSE)
sys.source("dev/compare_packages.R", envir = globalenv())
options(mc.cores = 1L)
DGP    <- klue_dgp(n_generic = 4, n_alternatives = 3)
DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)
OUT    <- "output/stress_replicate_blocked.rds"
stamp  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

results <- readRDS(OUT)
keys <- names(results)
stamp("Apollo robust re-run: %d cells", length(keys))

for (k in keys) {
  cell <- results[[k]]
  if (isTRUE(cell$apollo_robust)) { stamp("SKIP %s (done)", k); next }
  seed_i <- as.integer(1000 * cell$K + 100 * (cell$kap * 100) + cell$seed)
  d <- klue_simulate(N_per_class = 60, true_K = cell$K, separation = cell$kap,
                     heterogeneity = cell$sig, seed = seed_i, dgp = DGP, design = DESIGN)
  stamp("RUN  %s (K=%d kap=%.2f sig=%.2f)", k, cell$K, cell$kap, cell$sig)
  ap <- run_apollo_ss_published_blocked(d$database, cell$K, DGP)
  cell$LL[["apollo_searchStart"]] <- ap$LL
  best <- suppressWarnings(max(cell$LL[is.finite(cell$LL)]))
  cell$gap <- cell$LL - best
  cell$apollo_robust <- TRUE
  results[[k]] <- cell
  saveRDS(results, OUT)
  stamp("  %s apollo_LL=%s gap=%+.2f", k, format(ap$LL), cell$gap[["apollo_searchStart"]])
}

stamp("ALL APOLLO RE-RUN DONE")
# headline: mean gap + stuck per rung x method (now with finite Apollo)
meths <- c("klue_ml","klue_em","apollo_searchStart","gmnl_lc")
for (rg in c("moderate","hard","very_hard")) {
  cells <- results[sapply(results, function(c) c$rung == rg)]
  cat(sprintf("%-10s", rg))
  for (m in meths) {
    g <- sapply(cells, function(c) c$gap[[m]]); g <- g[is.finite(g)]
    cat(sprintf(" | %s %+5.2f (stuck %d/%d)", m, if (length(g)) mean(g) else NA, sum(g < -0.5), length(g)))
  }
  cat("\n")
}
