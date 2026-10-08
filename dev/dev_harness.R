# dev/dev_harness.R
# Reusable iteration loop for the v9 revision work (EM estimator, MNL-perturbation
# starts, correlated-attribute DGP). Loads klue in-place via pkgload (no install).
# Run from repo root:  Rscript dev/dev_harness.R
#
# This file lives OUTSIDE the package on purpose: it is scratch tooling, not part
# of the klue build, so it never touches R CMD check.

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))

# --- tiny fixture: fast enough to iterate on (sub-second per estimate) ---------
dev_data <- function(true_K = 2, sep = 1.0, sigma = 0.2, Npc = 120, T = 12, seed = 42) {
  klue_simulate(N_per_class = Npc, T_tasks = T, true_K = true_K,
                separation = sep, heterogeneity = sigma, seed = seed)
}

# --- base-engine smoke test: confirms the known-good baseline ------------------
smoke_base <- function() {
  d <- dev_data()
  t1 <- Sys.time()
  m <- klue_lcmnl(d$database, C = 2)
  cat(sprintf("[base] klue_lcmnl C=2: %.2fs conv=%s LL=%.1f BIC=%.1f best=%s\n",
              as.numeric(Sys.time() - t1), m$converged, m$LL, m$BIC, m$best_method))
  invisible(m)
}

# --- EM acceptance test: from IDENTICAL starts, EM LL must match direct-ML -----
# Both estimators are MLEs of the same LCMNL likelihood, so from the same
# starting values they should converge to the same optimum.
smoke_em <- function(true_K = 3, C = 3, sep = 1.0, sigma = 0.2, tol = 1e-3) {
  d <- dev_data(true_K = true_K, sep = sep, sigma = sigma)
  starts <- get_all_starts(d$database, C)          # 6 clustering-derived starts
  cat(sprintf("[em] K*=%d C=%d  | comparing EM vs direct-ML from identical starts\n",
              true_K, C))
  ok <- TRUE
  for (nm in names(starts)) {
    s <- starts[[nm]]
    if (is.null(s)) next
    ml <- estimate_lcmnl(d$database, C, start_betas = s$betas, start_shares = s$shares)
    em <- estimate_lcmnl_em(d$database, C, start_betas = s$betas, start_shares = s$shares)
    dLL <- em$LL - ml$LL
    # Correctness: EM must never be meaningfully WORSE than direct-ML from the
    # same start (both maximise the same likelihood). EM > ML is fine and means
    # EM escaped a local optimum that trapped BFGS.
    flag <- if (abs(dLL) <= tol) "match" else if (dLL > tol) "EM better" else "FAIL: EM worse"
    if (dLL < -tol) ok <- FALSE
    cat(sprintf("    %-12s ML=%.4f  EM=%.4f  dLL=%+.5f  iters=%3d  [%s]\n",
                nm, ml$LL, em$LL, dLL, em$em_iters, flag))
  }
  cat(sprintf("[em] result: %s\n",
              if (ok) "PASS (EM never worse than direct-ML)" else "FAIL (EM worse on some start)"))
  invisible(ok)
}

# --- start-strategy comparison: clustering vs MNL-perturbation vs diffuse ------
# Sanity check that MNL-perturbation is a genuine middle-ground baseline:
# better than diffuse random, the realistic competitor clustering must beat.
smoke_starts <- function(true_K = 4, sep = 0.75, sigma = 0.2, n = 20) {
  d <- dev_data(true_K = true_K, sep = sep, sigma = sigma, Npc = 150, T = 20)
  db <- d$database
  clust <- klue_lcmnl(db, true_K)            # 6 clustering starts, keep best
  perturb <- get_mnl_perturbation_starts(db, true_K, n_starts = n, seed = 1)
  pLL <- vapply(perturb, function(s) {
    r <- estimate_lcmnl(db, true_K, start_betas = s$betas, start_shares = s$shares)
    if (r$converged) r$LL else -Inf
  }, numeric(1))
  rLL <- vapply(seq_len(n), function(r) {
    set.seed(99000 + r)
    rb <- cbind(matrix(rnorm(true_K * 4, 0, 2), true_K, 4), -exp(rnorm(true_K, 0, 1)))
    rr <- estimate_lcmnl(db, true_K, start_betas = rb)
    if (rr$converged) rr$LL else -Inf
  }, numeric(1))
  best <- max(c(clust$LL, pLL, rLL))
  hit <- function(x) 100 * mean(x >= best - 0.1)
  cat(sprintf("[starts] K*=%d sep=%.2f  global LL=%.1f\n", true_K, sep, best))
  cat(sprintf("    clustering (6 starts): best LL=%.1f  at-global=%s\n",
              clust$LL, ifelse(clust$LL >= best - 0.1, "YES", "no")))
  cat(sprintf("    MNL-perturbation (%d): at-global=%.0f%%\n", n, hit(pLL)))
  cat(sprintf("    diffuse random  (%d): at-global=%.0f%%\n", n, hit(rLL)))
}

if (sys.nframe() == 0L) {
  smoke_base()
  cat("\n")
  smoke_em()
  cat("\n")
  smoke_starts()
}
