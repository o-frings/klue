# dev/benchmark_interactions.R
# ------------------------------------------------------------------------------
# Show that attribute INTERACTIONS work in klue with correct standard errors.
# klue's utility is linear over the columns x1_j..x{ng}_j, price_j, so an
# interaction is entered by building it as an extra attribute column and
# raising n_generic. Here we take C = 2 simulated data (4 attributes + price)
# and append two genuine product columns:
#   x5_j = x1_j * x2_j      (attribute x attribute)
#   x6_j = x3_j * price_j   (attribute x price)
# giving a 6-attribute model, then check klue's point estimates AND classical /
# robust SEs against Apollo for every coefficient, the interactions included.
# (The choices were simulated without these interactions, so their true betas
# are ~0; the test is that klue and Apollo AGREE on the estimate and its SE.)
#
# Seeded, self-contained. Run:  Rscript dev/benchmark_interactions.R
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

C <- 2; seed <- 101
DGP0 <- klue_dgp(n_generic = 4, n_alternatives = 3)      # simulation DGP
DGP  <- klue_dgp(n_generic = 6, n_alternatives = 3)      # estimation DGP (+2 interactions)
J <- DGP$n_alternatives; NG <- DGP$n_generic; NB <- DGP$n_beta; NASC <- DGP$n_asc

d  <- klue_simulate(N_per_class = 250, T_tasks = 12, true_K = C,
                    separation = 2.0, heterogeneity = 0.1, seed = seed, dgp = DGP0)
db <- d$database
for (j in 1:J) {                                         # append interaction columns
  db[[paste0("x5_", j)]] <- db[[paste0("x1_", j)]] * db[[paste0("x2_", j)]]
  db[[paste0("x6_", j)]] <- db[[paste0("x3_", j)]] * db[[paste0("price_", j)]]
}

# ---- Apollo LC spec (NG-generic; identical structure to the estimation DGP) --
make_lcPars <- function(C) {
  pn <- c(paste0("asc_alt", 1:(J - 1)), paste0("b_x", 1:NG), "b_price")
  L <- c("function(apollo_beta, apollo_inputs) {", "  lcpars <- list()")
  for (p in pn) L <- c(L, sprintf('  lcpars[["%s"]] <- list(%s)', p,
                                   paste(paste0(p, "_", 1:C), collapse = ", ")))
  L <- c(L, sprintf('  s <- list(classes = c(%s), utilities = list(%s))',
                    paste(sprintf("class_%d = %d", 1:C, 1:C), collapse = ", "),
                    paste(sprintf("class_%d = delta_%d", 1:C, 1:C), collapse = ", ")),
         '  lcpars[["pi_values"]] <- apollo_classAlloc(s)', '  return(lcpars) }')
  fn <- eval(parse(text = paste(L, collapse = "\n"))); environment(fn) <- globalenv(); fn
}
make_probs <- function(C) {
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
perms <- function(n) if (n == 1) matrix(1) else
  do.call(rbind, lapply(1:n, function(i)
    cbind(i, matrix(perms(n - 1), ncol = n - 1) + (perms(n - 1) >= i))))
best_align <- function(A, B) { P <- perms(nrow(A)); bd <- Inf; best <- NULL
  for (r in 1:nrow(P)) { dd <- sum((A - B[P[r, ], , drop = FALSE])^2)
    if (dd < bd) { bd <- dd; best <- P[r, ] } }; best }

# ---- klue fit + SEs (uses the 6-attribute estimation DGP) --------------------
k <- klue_lcmnl(db, C, dgp = DGP)
k_se  <- sqrt(diag(k$vcov)); k_rse <- sqrt(diag(k$robust_vcov))

# ---- Apollo fit + SEs --------------------------------------------------------
start <- c()
for (cc in 1:C) { for (jj in 1:NASC) start[paste0("asc_alt", jj, "_", cc)] <- 0
  for (a in 1:NG) start[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)
  start[paste0("b_price_", cc)] <- -0.3 }
for (cc in 1:C) start[paste0("delta_", cc)] <- 0
cleanup_apollo()
apollo_control <<- list(modelName = "int", modelDescr = "int", indivID = "ID",
                        nCores = 1L, outputDirectory = tempdir())
apollo_beta <<- start; apollo_fixed <<- "delta_1"
apollo_lcPars <<- make_lcPars(C); apollo_probabilities <<- make_probs(C)
apollo_inputs <<- apollo_validateInputs(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                    database = db, apollo_control = apollo_control, silent = TRUE)
m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                     estimate_settings = list(writeIter = FALSE, silent = TRUE, bootstrapSE = 0))
a_se <- sqrt(diag(m$varcov)); a_rse <- sqrt(diag(m$robvarcov))

# ---- align classes and compare -----------------------------------------------
kb <- k$betas
ab <- t(vapply(1:C, function(cc)
         c(vapply(1:NG, function(a) m$estimate[[paste0("b_x", a, "_", cc)]], 0),
           m$estimate[[paste0("b_price_", cc)]]), numeric(NB)))
p <- best_align(kb, ab)

lbl <- c(paste0("b", 1:NB), paste0("asc", 1:NASC))
lab <- c("x1","x2","x3","x4","x1:x2 (int)","x3:price (int)","price","asc1","asc2")
rows <- list()
for (cc in 1:C) { ac <- p[cc]
  kn <- c(paste0("b", 1:NB, "_class", cc), paste0("asc", 1:NASC, "_class", cc))
  an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac), paste0("asc_alt", 1:NASC, "_", ac))
  for (i in seq_along(kn)) rows[[length(rows)+1]] <- data.frame(
    class = cc, term = lab[i],
    klue_se = k_se[[kn[i]]], apollo_se = a_se[[an[i]]],
    klue_rse = k_rse[[kn[i]]], apollo_rse = a_rse[[an[i]]]) }
tab <- do.call(rbind, rows)
tab$klue_est <- vapply(seq_len(nrow(tab)), function(r) k$par[[
  c(paste0("b", 1:NB, "_class", tab$class[r]), paste0("asc", 1:NASC, "_class", tab$class[r]))[
    match(tab$term[r], lab)]]], 0)
tab$apollo_est <- vapply(seq_len(nrow(tab)), function(r) {
  ac <- p[tab$class[r]]
  an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac), paste0("asc_alt", 1:NASC, "_", ac))
  m$estimate[[an[match(tab$term[r], lab)]]] }, 0)
tab$rel_se  <- abs(tab$klue_se  - tab$apollo_se)  / tab$apollo_se
tab$rel_rse <- abs(tab$klue_rse - tab$apollo_rse) / tab$apollo_rse

cat(sprintf("\nLL: klue = %.4f   apollo = %.4f   dLL = %+.2e\n", k$LL, m$maximum, m$maximum - k$LL))
disp <- tab[, c("class","term","klue_est","apollo_est","klue_se","apollo_se","klue_rse","apollo_rse")]
nc <- vapply(disp, is.numeric, logical(1)); disp[nc] <- round(disp[nc], 4)
print(disp, row.names = FALSE)
cat(sprintf("\nmax |dEst| = %.2e   max rel-diff classical SE = %.2e   robust SE = %.2e\n",
            max(abs(tab$klue_est - tab$apollo_est)), max(tab$rel_se), max(tab$rel_rse)))
cat("Interaction rows are 'x1:x2 (int)' and 'x3:price (int)'.\n")
