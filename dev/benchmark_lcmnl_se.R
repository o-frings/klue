# dev/benchmark_lcmnl_se.R
# ------------------------------------------------------------------------------
# Benchmark klue's NEW multi-class (C > 1) LCMNL covariance against Apollo.
# klue now returns, for any C:
#   fit$vcov         = inverse observed information (classical)
#   fit$robust_vcov  = sandwich clustered on respondent (robust)
# Apollo reports the same two: m$varcov (classical) and m$robvarcov (robust,
# clustered on indivID for panel data). This script fits the SAME LCMNL with
# both and compares standard errors:
#   - taste params (per-class betas + ASCs): SE by SE, after label alignment
#   - class shares pi_c: SE via the delta method (reference-invariant)
# for classical and robust variants.
#
# Seeds: klue_simulate seeds 101 (C = 2) and 202 (C = 3); klue's clustering
# starts use seed 123; Apollo starts from a fixed vector, and its BGW fit and
# covariances use no random numbers that affect the results. Self-contained.
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/benchmark_lcmnl_se.R
# Output: console table + output/benchmark_lcmnl_se.rds
#
# Superseded for the paper by dev/check_lcmnl_standard_errors.R, which starts
# Apollo at klue's solution. Here Apollo converges from its own start, so the
# two optima differ slightly: at C = 3 the classical SEs differ by 1.1e-5 in
# relative terms (also under klue 0.10.0), against 3e-8 at the same point.
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

DGP  <- klue_dgp(n_generic = 4, n_alternatives = 3)
J    <- DGP$n_alternatives; NG <- DGP$n_generic
NB   <- DGP$n_beta; NASC <- DGP$n_asc

# ---- Apollo LC spec (independent of klue's estimator) ------------------------
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
    '  lcpars[["pi_values"]] <- apollo_classAlloc(s)', '  return(lcpars) }')
  fn <- eval(parse(text = paste(lines, collapse = "\n"))); environment(fn) <- globalenv(); fn
}
make_probs <- function(C) {
  alt <- paste(sprintf("alt%d = %d", 1:J, 1:J), collapse = ", ")
  avl <- paste(sprintf("alt%d = 1", 1:J), collapse = ", ")
  lines <- c('function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))', '  P <- list()',
    sprintf('  ms <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE)', alt, avl))
  for (s in 1:C) {
    lines <- c(lines, '  V <- list()')
    for (j in 1:J) {
      terms <- c(); if (j < J) terms <- c(terms, sprintf('asc_alt%d[[%d]]', j, s))
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
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)', '  return(P) }')
  fn <- eval(parse(text = paste(lines, collapse = "\n"))); environment(fn) <- globalenv(); fn
}
neutral_start <- function(C) {
  b <- c()
  for (cc in 1:C) {
    for (jj in 1:NASC) b[paste0("asc_alt", jj, "_", cc)] <- 0
    for (a in 1:NG)    b[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)
    b[paste0("b_price_", cc)] <- -0.5
  }
  for (cc in 1:C) b[paste0("delta_", cc)] <- 0
  b
}

# label-switching alignment: permutation of B's rows minimising SSE to A
perms <- function(n) if (n == 1) matrix(1) else
  do.call(rbind, lapply(1:n, function(i)
    cbind(i, matrix(perms(n - 1), ncol = n - 1) + (perms(n - 1) >= i))))
best_align <- function(A, B) {
  P <- perms(nrow(A)); bd <- Inf; best <- NULL
  for (r in 1:nrow(P)) { d <- sum((A - B[P[r, ], , drop = FALSE])^2)
    if (d < bd) { bd <- d; best <- P[r, ] } }
  best
}

# delta-method SE of class shares pi = softmax(eta), eta has a 0 at fixed_idx,
# free entries = deltas; Sig = covariance of the free deltas (|free| x |free|).
share_se <- function(deltas, fixed_idx, Sig, C) {
  eta <- numeric(C); free <- setdiff(1:C, fixed_idx); eta[free] <- deltas
  pi <- exp(eta - max(eta)); pi <- pi / sum(pi)
  Jf <- outer(1:C, free, function(a, b) pi[a] * ((a == b) - pi[b]))  # C x |free|
  list(pi = pi, se = sqrt(pmax(diag(Jf %*% Sig %*% t(Jf)), 0)))
}

CASES <- list(list(C = 2, seed = 101), list(C = 3, seed = 202))
results <- list()

for (case in CASES) {
  C <- case$C; seed <- case$seed
  cat(sprintf("\n================  C = %d  (seed %d)  ================\n", C, seed))
  d  <- klue_simulate(N_per_class = 200, T_tasks = 12, true_K = C,
                      separation = 2.0, heterogeneity = 0.1, seed = seed, dgp = DGP)
  db <- d$database; npc <- NB + NASC

  # ---- klue fit + SEs -------------------------------------------------------
  k <- klue_lcmnl(db, C, dgp = DGP)
  k_se  <- sqrt(diag(k$vcov)); k_rse <- sqrt(diag(k$robust_vcov))

  # ---- Apollo fit + SEs -----------------------------------------------------
  cleanup_apollo()
  apollo_control <<- list(modelName = "bm", modelDescr = "bm", indivID = "ID",
                          nCores = 1L, outputDirectory = tempdir())
  apollo_beta          <<- neutral_start(C)
  apollo_fixed         <<- "delta_1"
  apollo_lcPars        <<- make_lcPars(C)
  apollo_probabilities <<- make_probs(C)
  apollo_inputs        <<- apollo_validateInputs(apollo_beta = apollo_beta,
                             apollo_fixed = apollo_fixed, database = db,
                             apollo_control = apollo_control, silent = TRUE)
  m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                       estimate_settings = list(writeIter = FALSE, silent = TRUE,
                                                bootstrapSE = 0))
  a_vc  <- m$varcov; a_rvc <- m$robvarcov
  a_se  <- sqrt(diag(a_vc)); a_rse <- sqrt(diag(a_rvc))

  # ---- align Apollo classes to klue classes (by betas) ----------------------
  kb <- k$betas
  ab <- t(vapply(1:C, function(cc)
           c(vapply(1:NG, function(a) m$estimate[[paste0("b_x", a, "_", cc)]], 0),
             m$estimate[[paste0("b_price_", cc)]]), numeric(NB)))
  p <- best_align(kb, ab)   # apollo class p[c] corresponds to klue class c

  # ---- taste-parameter SE comparison, class by class ------------------------
  rows <- list()
  for (cc in 1:C) {
    ac <- p[cc]
    kn <- c(paste0("b", 1:NB, "_class", cc), paste0("asc", 1:NASC, "_class", cc))
    an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac),
            paste0("asc_alt", 1:NASC, "_", ac))
    lbl <- c(paste0("b", 1:NB), paste0("asc", 1:NASC))
    for (i in seq_along(kn))
      rows[[length(rows) + 1]] <- data.frame(
        class = cc, param = lbl[i],
        klue_se = k_se[[kn[i]]],   apollo_se = a_se[[an[i]]],
        klue_rse = k_rse[[kn[i]]], apollo_rse = a_rse[[an[i]]])
  }
  tab <- do.call(rbind, rows)
  tab$rel_se  <- abs(tab$klue_se  - tab$apollo_se)  / tab$apollo_se
  tab$rel_rse <- abs(tab$klue_rse - tab$apollo_rse) / tab$apollo_rse

  # ---- class-share SEs via delta method (reference-invariant) ---------------
  kd <- k$par[paste0("delta", 1:(C - 1))]                       # klue: class C fixed
  kS <- share_se(kd, C, k$vcov[paste0("delta", 1:(C-1)), paste0("delta", 1:(C-1)), drop=FALSE], C)
  ad <- vapply(2:C, function(cc) m$estimate[[paste0("delta_", cc)]], 0)  # apollo: class 1 fixed
  aSig <- a_vc[paste0("delta_", 2:C), paste0("delta_", 2:C), drop = FALSE]
  aS <- share_se(ad, 1, aSig, C)
  # align apollo share SEs to klue class order
  share_cmp <- data.frame(class = 1:C,
                          klue_share = kS$pi, apollo_share = aS$pi[p],
                          klue_share_se = kS$se, apollo_share_se = aS$se[p])

  cat(sprintf("[LL] klue = %.4f   apollo = %.4f\n", k$LL, m$maximum))
  cat("--- taste-parameter standard errors (klue vs apollo) ---\n")
  disp <- tab[, c("class","param","klue_se","apollo_se","klue_rse","apollo_rse")]
  numc <- vapply(disp, is.numeric, logical(1)); disp[numc] <- round(disp[numc], 4)
  print(disp, row.names = FALSE)
  cat(sprintf("max rel diff  classical SE = %.2e   robust SE = %.2e\n",
              max(tab$rel_se), max(tab$rel_rse)))
  cat("--- class-share SEs (delta method) ---\n")
  print(round(share_cmp, 4), row.names = FALSE)

  results[[paste0("C", C)]] <- list(C = C, seed = seed, LL_klue = k$LL,
    LL_apollo = m$maximum, taste = tab, shares = share_cmp,
    max_rel_se = max(tab$rel_se), max_rel_rse = max(tab$rel_rse))
}

if (!dir.exists("output")) dir.create("output")
saveRDS(results, "output/benchmark_lcmnl_se.rds")

cat("\n================  SUMMARY  ================\n")
for (r in results)
  cat(sprintf("C=%d | max rel diff: classical SE = %.2e  robust SE = %.2e | share SE max rel = %.2e\n",
              r$C, r$max_rel_se, r$max_rel_rse,
              max(abs(r$shares$klue_share_se - r$shares$apollo_share_se) /
                  r$shares$apollo_share_se)))
cat("Saved: output/benchmark_lcmnl_se.rds\n")
