# dev/benchmark_concomitant.R
# ------------------------------------------------------------------------------
# Benchmark klue's NEW concomitant LCMNL (covariate-driven class membership,
# estimate_lcmnl_cov) against Apollo. Class membership: pi_c = softmax(Z'gamma_c)
# with Z = (intercept, Z1, Z2) and the last class as reference. Fit the SAME
# model in both and compare, for C = 2:
#   - taste params (per-class betas + ASCs): estimate + classical/robust SE
#   - membership params (gamma: intercept, Z1, Z2): estimate + classical/robust SE
#   - per-respondent class shares
# For C = 2 there is one non-reference class, so label alignment is a possible
# swap; the membership contrast then flips sign (SE unchanged).
#
# Seeded, self-contained. Run:  Rscript dev/benchmark_concomitant.R
# Output: console tables + output/benchmark_concomitant.rds
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

C <- 2; seed <- 7
DGP <- klue_dgp(n_generic = 4, n_alternatives = 3)
J <- DGP$n_alternatives; NG <- DGP$n_generic; NB <- DGP$n_beta; NASC <- DGP$n_asc
COV <- c("Z1", "Z2")

d  <- klue_simulate(N_per_class = 300, T_tasks = 12, true_K = C, separation = 2.0,
                    heterogeneity = 0.1, seed = seed, dgp = DGP,
                    covariates = TRUE, covariate_strength = 1.5)
db <- d$database

# ---- Apollo concomitant LC spec ----------------------------------------------
make_lcPars <- function() {
  tn <- c(paste0("asc_alt", 1:(J - 1)), paste0("b_x", 1:NG), "b_price")
  L <- c("function(apollo_beta, apollo_inputs) {", "  lcpars <- list()")
  for (p in tn) L <- c(L, sprintf('  lcpars[["%s"]] <- list(%s)', p,
                                   paste(paste0(p, "_", 1:C), collapse = ", ")))
  # membership utilities: class 1 = g0 + gZ1*Z1 + gZ2*Z2 ; class 2 = reference 0
  util <- c("class_1 = g0_1 + gZ1_1*Z1 + gZ2_1*Z2", "class_2 = 0")
  L <- c(L, sprintf('  s <- list(classes = c(class_1 = 1, class_2 = 2), utilities = list(%s))',
                    paste(util, collapse = ", ")),
         '  lcpars[["pi_values"]] <- apollo_classAlloc(s)', '  return(lcpars) }')
  fn <- eval(parse(text = paste(L, collapse = "\n"))); environment(fn) <- globalenv(); fn
}
make_probs <- function() {
  alt <- paste(sprintf("alt%d = %d", 1:J, 1:J), collapse = ", ")
  avl <- paste(sprintf("alt%d = 1", 1:J), collapse = ", ")
  L <- c('function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))', '  P <- list()',
    sprintf('  ms <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE)', alt, avl))
  for (s in 1:C) {
    L <- c(L, '  V <- list()')
    for (j in 1:J) {
      t <- c(); if (j < J) t <- c(t, sprintf('asc_alt%d[[%d]]', j, s))
      for (a in 1:NG) t <- c(t, sprintf('b_x%d[[%d]] * x%d_%d', a, s, a, j))
      t <- c(t, sprintf('b_price[[%d]] * price_%d', s, j))
      L <- c(L, sprintf('  V[["alt%d"]] <- %s', j, paste(t, collapse = " + ")))
    }
    L <- c(L, '  ms$utilities <- V',
      sprintf('  P[["class_%d"]] <- apollo_mnl(ms, functionality)', s),
      sprintf('  P[["class_%d"]] <- apollo_panelProd(P[["class_%d"]], apollo_inputs, functionality)', s, s))
  }
  L <- c(L, sprintf('  lcs <- list(inClassProb = list(%s), classProb = pi_values)',
                    paste(sprintf('class_%d = P[["class_%d"]]', 1:C, 1:C), collapse = ", ")),
    '  P[["model"]] <- apollo_lc(lcs, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)', '  return(P) }')
  fn <- eval(parse(text = paste(L, collapse = "\n"))); environment(fn) <- globalenv(); fn
}

# ---- klue fit ----------------------------------------------------------------
k <- estimate_lcmnl_cov(db, C, membership = COV, dgp = DGP)
k_se <- sqrt(diag(k$vcov)); k_rse <- sqrt(diag(k$robust_vcov))

# ---- Apollo fit --------------------------------------------------------------
start <- c()
for (cc in 1:C) { for (jj in 1:NASC) start[paste0("asc_alt", jj, "_", cc)] <- 0
  for (a in 1:NG) start[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)
  start[paste0("b_price_", cc)] <- -0.5 }
start[c("g0_1", "gZ1_1", "gZ2_1")] <- 0
cleanup_apollo()
apollo_control <<- list(modelName = "conc", modelDescr = "conc", indivID = "ID",
                        nCores = 1L, outputDirectory = tempdir())
apollo_beta <<- start; apollo_fixed <<- character(0)
apollo_lcPars <<- make_lcPars(); apollo_probabilities <<- make_probs()
apollo_inputs <<- apollo_validateInputs(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                    database = db, apollo_control = apollo_control, silent = TRUE)
m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                     estimate_settings = list(writeIter = FALSE, silent = TRUE, bootstrapSE = 0))
a_se <- sqrt(diag(m$varcov)); a_rse <- sqrt(diag(m$robvarcov))

# ---- align classes by taste betas (C = 2: identity or swap) ------------------
kb <- k$betas
ab <- t(vapply(1:C, function(cc)
         c(vapply(1:NG, function(a) m$estimate[[paste0("b_x", a, "_", cc)]], 0),
           m$estimate[[paste0("b_price_", cc)]]), numeric(NB)))
swap <- sum((kb - ab)^2) > sum((kb - ab[c(2, 1), , drop = FALSE])^2)
p <- if (swap) c(2, 1) else c(1, 2)

cat(sprintf("\nLL: klue = %.4f   apollo = %.4f   dLL = %+.2e   (class alignment: %s)\n",
            k$LL, m$maximum, m$maximum - k$LL, if (swap) "SWAP" else "identity"))

# ---- taste SE comparison -----------------------------------------------------
lab <- c(paste0("x", 1:NG), "price", paste0("asc", 1:NASC))
trows <- list()
for (cc in 1:C) { ac <- p[cc]
  kn <- c(paste0("b", 1:NB, "_class", cc), paste0("asc", 1:NASC, "_class", cc))
  an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac), paste0("asc_alt", 1:NASC, "_", ac))
  for (i in seq_along(kn)) trows[[length(trows) + 1]] <- data.frame(
    class = cc, term = lab[i],
    klue_est = k$par[[kn[i]]], apollo_est = m$estimate[[an[i]]],
    klue_se = k_se[[kn[i]]], apollo_se = a_se[[an[i]]],
    klue_rse = k_rse[[kn[i]]], apollo_rse = a_rse[[an[i]]]) }
taste <- do.call(rbind, trows)

# ---- membership (gamma) comparison; flip sign if classes are swapped ---------
sgn <- if (swap) -1 else 1
kg  <- c("(Intercept)_class1", "Z1_class1", "Z2_class1")
ag  <- c("g0_1", "gZ1_1", "gZ2_1")
memb <- data.frame(
  term = c("(Intercept)", "Z1", "Z2"),
  klue_est   = sapply(kg, function(z) k$par[[z]]),
  apollo_est = sgn * sapply(ag, function(z) m$estimate[[z]]),
  klue_se    = sapply(kg, function(z) k_se[[z]]),
  apollo_se  = sapply(ag, function(z) a_se[[z]]),
  klue_rse   = sapply(kg, function(z) k_rse[[z]]),
  apollo_rse = sapply(ag, function(z) a_rse[[z]]), row.names = NULL)

rnd <- function(df) { nc <- vapply(df, is.numeric, logical(1)); df[nc] <- round(df[nc], 4); df }
cat("\n--- taste-parameter estimates + SEs ---\n"); print(rnd(taste), row.names = FALSE)
cat("\n--- MEMBERSHIP (gamma) estimates + SEs ---\n"); print(rnd(memb), row.names = FALSE)

reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
cat(sprintf("\ntaste : max|dEst|=%.1e  classical SE reldiff=%.1e  robust SE reldiff=%.1e\n",
            max(abs(taste$klue_est - taste$apollo_est)),
            reldiff(taste$klue_se, taste$apollo_se), reldiff(taste$klue_rse, taste$apollo_rse)))
cat(sprintf("gamma : max|dEst|=%.1e  classical SE reldiff=%.1e  robust SE reldiff=%.1e\n",
            max(abs(memb$klue_est - memb$apollo_est)),
            reldiff(memb$klue_se, memb$apollo_se), reldiff(memb$klue_rse, memb$apollo_rse)))

if (!dir.exists("output")) dir.create("output")
saveRDS(list(taste = taste, memb = memb, swap = swap,
             LL_klue = k$LL, LL_apollo = m$maximum), "output/benchmark_concomitant.rds")
cat("Saved: output/benchmark_concomitant.rds\n")
