# dev/validate_lcmnl_apollo.R
# ------------------------------------------------------------------------------
# DEFINITIVE point-estimate cross-check: klue's hand-rolled LCMNL estimator vs
# Apollo's latent-class estimator, compared COEFFICIENT BY COEFFICIENT (per-class
# betas, ASCs, and class shares) for the genuine multi-class case C > 1.
#
# The existing dev harnesses (compare_packages.R, stress_replicate*.R) compare
# only the log-likelihood reached. This script closes the remaining gap: it
# checks that klue and Apollo agree on the actual PARAMETER VALUES, not just the
# LL. Three independent checks per C:
#   (A) anchored LL   : Apollo's likelihood evaluated AT klue's converged params
#                       must equal klue$LL  -> proves both code the same model.
#   (B) no-move       : Apollo estimation STARTED at klue's params must not move
#                       -> proves klue's estimate is Apollo's optimum too.
#   (C) free estimate : Apollo estimated from a neutral start must land on the
#                       same LL and (after label alignment) the same betas.
# Plus recovery of the known true betas by both estimators.
#
# The Apollo spec (apollo_mnl / apollo_lc / apollo_classAlloc) is built here,
# independent of klue's estimator, so this is a genuine external check.
#
# Seeded and self-contained. Run:  Rscript dev/validate_lcmnl_apollo.R
# Output: console table + output/validate_lcmnl_apollo.rds
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

DGP  <- klue_dgp(n_generic = 4, n_alternatives = 3)   # default per-alt ASC map
J    <- DGP$n_alternatives
NG   <- DGP$n_generic
NB   <- DGP$n_beta        # generic + price = 5
NASC <- DGP$n_asc         # 2 (alts 1,2; alt 3 reference)

# ---- Apollo LC spec builders (standard Apollo primitives; no klue estimator) --
# String-built to dodge Apollo's checkIndices rejection of for-loops in the
# utility function. Matches klue's utility exactly: V_j = X_j' beta_c + asc_j.
make_lcPars <- function(C) {
  pnames <- c(paste0("asc_alt", 1:(J - 1)), paste0("b_x", 1:NG), "b_price")
  lines <- c("function(apollo_beta, apollo_inputs) {", "  lcpars <- list()")
  for (p in pnames)
    lines <- c(lines, sprintf('  lcpars[["%s"]] <- list(%s)', p,
                              paste(paste0(p, "_", 1:C), collapse = ", ")))
  lines <- c(lines,
    sprintf('  s <- list(classes = c(%s), utilities = list(%s))',
            paste(sprintf("class_%d = %d", 1:C, 1:C), collapse = ", "),
            paste(sprintf("class_%d = delta_%d", 1:C, 1:C), collapse = ", ")),
    '  lcpars[["pi_values"]] <- apollo_classAlloc(s)',
    '  return(lcpars) }')
  fn <- eval(parse(text = paste(lines, collapse = "\n"))); environment(fn) <- globalenv(); fn
}
make_probs <- function(C) {
  alt  <- paste(sprintf("alt%d = %d", 1:J, 1:J), collapse = ", ")
  avl  <- paste(sprintf("alt%d = 1", 1:J), collapse = ", ")
  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    sprintf('  ms <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE)', alt, avl))
  for (s in 1:C) {
    lines <- c(lines, '  V <- list()')
    for (j in 1:J) {
      terms <- c()
      if (j < J) terms <- c(terms, sprintf('asc_alt%d[[%d]]', j, s))
      for (a in 1:NG) terms <- c(terms, sprintf('b_x%d[[%d]] * x%d_%d', a, s, a, j))
      terms <- c(terms, sprintf('b_price[[%d]] * price_%d', s, j))
      lines <- c(lines, sprintf('  V[["alt%d"]] <- %s', j, paste(terms, collapse = " + ")))
    }
    lines <- c(lines, '  ms$utilities <- V',
      sprintf('  P[["class_%d"]] <- apollo_mnl(ms, functionality)', s),
      sprintf('  P[["class_%d"]] <- apollo_panelProd(P[["class_%d"]], apollo_inputs, functionality)', s, s))
  }
  lines <- c(lines,
    sprintf('  lcs <- list(inClassProb = list(%s), classProb = pi_values)',
            paste(sprintf('class_%d = P[["class_%d"]]', 1:C, 1:C), collapse = ", ")),
    '  P[["model"]] <- apollo_lc(lcs, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)',
    '  return(P) }')
  fn <- eval(parse(text = paste(lines, collapse = "\n"))); environment(fn) <- globalenv(); fn
}

# klue par (per class [b_x1..b_xNG, b_price, asc1..ascNASC]) + shares -> apollo_beta.
# klue uses last class as share reference; Apollo fixes delta_1. Both are the same
# softmax, so map deltas to Apollo's convention: delta_c = log(pi_c) - log(pi_1).
klue_to_apollo_beta <- function(fit, C) {
  npc <- NB + NASC
  b <- c()
  for (cc in 1:C) {
    par_c <- fit$par[((cc - 1) * npc + 1):(cc * npc)]
    for (a in 1:NG) b[paste0("b_x", a, "_", cc)] <- par_c[a]
    b[paste0("b_price_", cc)] <- par_c[NB]
    for (jj in 1:NASC) b[paste0("asc_alt", jj, "_", cc)] <- par_c[NB + jj]
  }
  lp <- log(fit$class_probs)
  for (cc in 1:C) b[paste0("delta_", cc)] <- lp[cc] - lp[1]
  b
}
apollo_beta_to_mat <- function(est, C) {
  betas <- matrix(0, C, NB); ascs <- matrix(0, C, NASC)
  for (cc in 1:C) {
    for (a in 1:NG) betas[cc, a] <- est[paste0("b_x", a, "_", cc)]
    betas[cc, NB] <- est[paste0("b_price_", cc)]
    for (jj in 1:NASC) ascs[cc, jj] <- est[paste0("asc_alt", jj, "_", cc)]
  }
  d <- vapply(1:C, function(cc) {
    nm <- paste0("delta_", cc); if (nm %in% names(est)) est[[nm]] else 0
  }, numeric(1))
  pi <- exp(d - max(d)); pi <- pi / sum(pi)
  list(betas = betas, ascs = ascs, shares = as.numeric(pi))
}

# Best class permutation aligning matrix B (rows) to reference A, min sum-sq.
perms <- function(n) if (n == 1) matrix(1) else
  do.call(rbind, lapply(1:n, function(i)
    cbind(i, matrix(perms(n - 1), ncol = n - 1) +
              (perms(n - 1) >= i))[, , drop = FALSE]))
best_align <- function(A, B) {
  C <- nrow(A); P <- perms(C); best <- NULL; bd <- Inf
  for (r in 1:nrow(P)) {
    p <- P[r, ]; d <- sum((A - B[p, , drop = FALSE])^2)
    if (d < bd) { bd <- d; best <- p }
  }
  best
}

# ---- Apollo LL evaluated at a fixed beta (no estimation) ---------------------
apollo_LL_at <- function(beta, C, db, probs_fn) {
  cleanup_apollo()
  apollo_control      <<- list(modelName = "chk", modelDescr = "chk",
                               indivID = "ID", nCores = 1L, outputDirectory = tempdir())
  apollo_beta         <<- beta
  apollo_fixed        <<- "delta_1"
  apollo_lcPars       <<- make_lcPars(C)
  apollo_probabilities <<- probs_fn
  apollo_inputs       <<- apollo_validateInputs(apollo_beta = apollo_beta,
                            apollo_fixed = apollo_fixed, database = db,
                            apollo_control = apollo_control, silent = TRUE)
  P <- apollo_probabilities(beta, apollo_inputs, "estimate")
  lik <- if (is.list(P)) P[["model"]] else P   # prepareProb strips to a vector
  sum(log(lik))
}

# ---- Apollo free/anchored estimation -----------------------------------------
apollo_fit <- function(start, C, db, probs_fn) {
  cleanup_apollo()
  apollo_control      <<- list(modelName = "fit", modelDescr = "fit",
                               indivID = "ID", nCores = 1L, outputDirectory = tempdir())
  apollo_beta         <<- start
  apollo_fixed        <<- "delta_1"
  apollo_lcPars       <<- make_lcPars(C)
  apollo_probabilities <<- probs_fn
  apollo_inputs       <<- apollo_validateInputs(apollo_beta = apollo_beta,
                            apollo_fixed = apollo_fixed, database = db,
                            apollo_control = apollo_control, silent = TRUE)
  m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                       estimate_settings = list(writeIter = FALSE, silent = TRUE,
                                                hessianRoutine = "none"))
  list(est = m$estimate, LL = as.numeric(m$maximum))
}
neutral_start <- function(C) {
  b <- c()
  for (cc in 1:C) {
    for (jj in 1:NASC) b[paste0("asc_alt", jj, "_", cc)] <- 0
    for (a in 1:NG)    b[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)  # break symmetry
    b[paste0("b_price_", cc)] <- -0.5
  }
  for (cc in 1:C) b[paste0("delta_", cc)] <- 0
  b
}

# ------------------------------------------------------------------------------
CASES <- list(list(C = 2, seed = 101), list(C = 3, seed = 202))
results <- list()

for (case in CASES) {
  C <- case$C; seed <- case$seed
  cat(sprintf("\n================  C = %d  (seed %d)  ================\n", C, seed))

  # Well-separated data so the global optimum is unambiguous for both estimators.
  d  <- klue_simulate(N_per_class = 200, T_tasks = 12, true_K = C,
                      separation = 2.0, heterogeneity = 0.1, seed = seed, dgp = DGP)
  db <- d$database

  probs_fn <- make_probs(C)

  # (1) klue reference fit
  k <- klue_lcmnl(db, C, dgp = DGP)
  cat(sprintf("[klue]   LL = %.4f\n", k$LL))
  kb <- klue_to_apollo_beta(k, C)

  # (A) Apollo LL AT klue's params
  ll_at <- apollo_LL_at(kb, C, db, probs_fn)
  cat(sprintf("[A] Apollo LL @ klue params = %.4f   dLL = %+.2e\n", ll_at, ll_at - k$LL))

  # (B) Apollo estimation STARTED at klue's params (should not move)
  fit_anchor <- apollo_fit(kb, C, db, probs_fn)
  a_anchor   <- apollo_beta_to_mat(fit_anchor$est, C)
  cat(sprintf("[B] Apollo from klue start: LL = %.4f   dLL = %+.2e\n",
              fit_anchor$LL, fit_anchor$LL - k$LL))

  # (C) Apollo estimation from a neutral start (independent)
  fit_free <- apollo_fit(neutral_start(C), C, db, probs_fn)
  a_free   <- apollo_beta_to_mat(fit_free$est, C)
  cat(sprintf("[C] Apollo from neutral start: LL = %.4f   dLL = %+.2e\n",
              fit_free$LL, fit_free$LL - k$LL))

  # klue's par is [class1(npc), .., classC(npc), delta_1..delta_{C-1}]; drop the
  # trailing deltas before reshaping into per-class [betas(NB), ascs(NASC)] rows.
  npc <- NB + NASC
  k_ascs <- matrix(k$par[1:(C * npc)], ncol = npc, byrow = TRUE)[, NB + 1:NASC, drop = FALSE]

  # Align Apollo(free) classes to klue classes, then compare coefficients
  kfull <- cbind(k$betas, k_ascs)
  afull_free <- cbind(a_free$betas, a_free$ascs)
  p <- best_align(kfull, afull_free)
  a_b <- a_free$betas[p, , drop = FALSE]
  a_a <- a_free$ascs[p, , drop = FALSE]
  a_s <- a_free$shares[p]

  d_beta  <- max(abs(k$betas - a_b))
  d_asc   <- max(abs(k_ascs - a_a))
  d_share <- max(abs(sort(k$class_probs) - sort(a_free$shares)))

  cat(sprintf("[C] aligned coeff diffs: max|dbeta| = %.2e  max|dasc| = %.2e  max|dshare| = %.2e\n",
              d_beta, d_asc, d_share))

  # Recovery of the known truth (align both to truth by betas)
  tp <- best_align(d$true_betas, k$betas)
  rec_klue   <- max(abs(d$true_betas - k$betas[tp, , drop = FALSE]))
  tp2 <- best_align(d$true_betas, a_free$betas)
  rec_apollo <- max(abs(d$true_betas - a_free$betas[tp2, , drop = FALSE]))
  cat(sprintf("    recovery vs truth: klue max|d| = %.3f   apollo max|d| = %.3f\n",
              rec_klue, rec_apollo))

  results[[paste0("C", C)]] <- list(
    C = C, seed = seed, LL_klue = k$LL, LL_apollo_at_klue = ll_at,
    LL_apollo_anchor = fit_anchor$LL, LL_apollo_free = fit_free$LL,
    d_beta = d_beta, d_asc = d_asc, d_share = d_share,
    rec_klue = rec_klue, rec_apollo = rec_apollo,
    klue_betas = k$betas, apollo_betas_aligned = a_b,
    klue_shares = k$class_probs, apollo_shares_aligned = a_s,
    true_betas = d$true_betas)
}

if (!dir.exists("output")) dir.create("output")
saveRDS(results, "output/validate_lcmnl_apollo.rds")

cat("\n================  SUMMARY  ================\n")
for (r in results) {
  cat(sprintf("C=%d | dLL(@klue)=%+.1e  dLL(free)=%+.1e | max|dbeta|=%.1e  max|dasc|=%.1e  max|dshare|=%.1e\n",
              r$C, r$LL_apollo_at_klue - r$LL_klue, r$LL_apollo_free - r$LL_klue,
              r$d_beta, r$d_asc, r$d_share))
}
cat("Saved: output/validate_lcmnl_apollo.rds\n")
