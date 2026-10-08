suppressMessages(pkgload::load_all("klue", quiet = TRUE)); source("studies/klue_studies.R")   # the klue_study_* drivers (moved out of the package 2026-09-17)
# ============================================================================
# Diagnose the standard-harness FAILs (2026-07-08, klue 0.9.2 vs klue 0.9.0):
#   study_recovery      rmse/bias rel diff ~1e-5/2e-4
#   study_concomitant   mean ARI rel diff ~4e-3
#   study_estimator     best_LL diff up to 1.24 LL on the random_partition arm
# while LL/BIC matched at 1e-5 in every direct-estimation component.
#
# Candidate 0.9.1/0.9.2 deltas:
#   (S) loose-tolerance multistart screen in klue_lcmnl (winner-tie reorder,
#       and method_results now hold screen-level LLs, which the estimator
#       study reads as the clustering arm and folds into global_LL)
#   (F) cluster-MNL fallback rule: non-converged finite-LL fits are kept
#       (affects clustering starts AND random-partition/perturbation starts)
#   (R) .with_seed RNG hygiene: clustering no longer clobbers the caller's
#       RNG stream
#
# Phase A (fast): generate the six estimator-study conditions in both engines
# and compare every start object bit-wise -> isolates (F)/(R) in seconds.
# Phase B (slow): re-run the new side of the three failed studies with
# options(klue.screen = FALSE) and compare to the old-side results stored in
# output/compare_klue_versions_standard.rds -> whatever still differs is not
# the screen.
#
# Deterministic: study drivers seed per condition; no ambient RNG use here.
# Output: verdicts to stdout; objects to output/diagnose_screen_effect.rds
# Usage:  Rscript dev/diagnose_screen_effect.R
# ============================================================================
suppressMessages({library(mclust); library(cluster)})
root <- normalizePath(".")
oldE <- new.env(parent = globalenv())
for (f in sort(list.files(file.path(root, "klue/R"), full.names = TRUE)))
  sys.source(f, envir = oldE)
newE <- new.env(parent = globalenv())
for (f in sort(list.files(file.path(root, "klue/R"), full.names = TRUE)))
  sys.source(f, envir = newE)

chk <- readRDS("output/compare_klue_versions_standard.rds")
out <- list()

# ---- Phase A: start-generation parity on the estimator-study conditions -----
cat("== Phase A: start objects, old vs new, estimator-study conditions 1:6 ==\n")
conds <- newE$.h1_conditions(NULL, 1:6)
outA <- list()
for (i in seq_len(nrow(conds))) {
  tK <- conds$true_K[i]; kap <- conds$kappa[i]
  sig <- conds$sigma[i]; rp <- conds$rep[i]
  seed <- newE$.cond_seed(tK, kap, sig, rp)
  dgp <- newE$klue_dgp()
  db_o <- oldE$klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                             separation = kap, heterogeneity = sig,
                             seed = seed)$database
  db_n <- newE$klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                             separation = kap, heterogeneity = sig,
                             seed = seed)$database
  cmp <- c(
    data      = identical(db_o, db_n),
    clustering = isTRUE(all.equal(oldE$get_all_starts(db_o, tK, dgp = dgp),
                                  newE$get_all_starts(db_n, tK, dgp = dgp),
                                  check.attributes = FALSE)),
    perturb   = isTRUE(all.equal(
      oldE$get_mnl_perturbation_starts(db_o, tK, n_starts = 10, seed = seed, dgp = dgp),
      newE$get_mnl_perturbation_starts(db_n, tK, n_starts = 10, seed = seed, dgp = dgp),
      check.attributes = FALSE)),
    partition = isTRUE(all.equal(
      oldE$get_random_partition_starts(db_o, tK, n_starts = 10, seed = seed, dgp = dgp),
      newE$get_random_partition_starts(db_n, tK, n_starts = 10, seed = seed, dgp = dgp),
      check.attributes = FALSE)),
    random    = isTRUE(all.equal(
      lapply(1:10, function(r) oldE$.random_start(seed, r, tK, dgp)),
      lapply(1:10, function(r) newE$.random_start(seed, r, tK, dgp)),
      check.attributes = FALSE)))
  cat(sprintf("  cond %d (K=%d kap=%.2f sig=%.2f rep=%d): %s\n", i, tK, kap,
              sig, rp, paste(names(cmp), ifelse(cmp, "OK", "DIFF"),
                             collapse = "  ")))
  outA[[i]] <- cmp
}
out$phaseA <- outA

# ---- Phase B: screen-off re-runs vs stored old results ----------------------
cat("\n== Phase B: new side re-run with klue.screen = FALSE ==\n")
options(klue.screen = FALSE)
runs <- list(
  study_recovery    = function() newE$klue_study_recovery(verbose = FALSE),
  study_concomitant = function() newE$klue_study_concomitant(verbose = FALSE),
  study_estimator_cond1to6 = function()
    newE$klue_study_estimator(n_random = 10L, n_perturb = 10L,
                              cond_idx = 1:6, verbose = FALSE))
for (comp in names(runs)) {
  cat("\n--", comp, "--\n")
  t0 <- proc.time()[3]
  rn <- runs[[comp]]()
  cat(sprintf("  (%.0fs)\n", proc.time()[3] - t0))
  ro <- chk[[comp]]$old
  if (is.data.frame(ro)) {
    key <- intersect(c("K","kappa","sigma","rep","estimator","strategy"), names(ro))
    if (length(key)) {   # without keys, ro[key] is 0-col and order() empties the frame
      ro <- ro[do.call(order, ro[key]), ]; rn <- rn[do.call(order, rn[key]), ]
    }
    for (col in setdiff(names(ro), key)) {
      if (!is.numeric(ro[[col]])) next
      d <- abs(ro[[col]] - rn[[col]])
      cat(sprintf("  %-16s max |diff| %.3g  -> %s\n", col,
                  suppressWarnings(max(d, na.rm = TRUE)),
                  if (isTRUE(all(d < 1e-9, na.rm = TRUE) &&
                             identical(is.na(ro[[col]]), is.na(rn[[col]]))))
                    "MATCH" else "DIFF"))
    }
  } else {
    for (nm in names(ro)) {
      eq <- all.equal(ro[[nm]], rn[[nm]], tolerance = 1e-9,
                      check.attributes = FALSE)
      cat(sprintf("  %-10s vs old (tol 1e-9): %s\n", nm,
                  paste(format(eq), collapse = "; ")))
    }
  }
  out[[comp]] <- list(old = ro, new_noscreen = rn, new_screen = chk[[comp]]$new)
}
saveRDS(out, "output/diagnose_screen_effect.rds")
cat("\nSaved output/diagnose_screen_effect.rds\n")
