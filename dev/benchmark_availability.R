# dev/benchmark_availability.R
# ------------------------------------------------------------------------------
# Benchmark klue's NEW alternative-availability handling against Apollo. klue
# masks unavailable alternatives in the per-task softmax (wide columns
# av_1..av_J, 1/0) instead of dropping tasks. To make an honest test, choices
# are RE-DRAWN from the availability-masked MNL using the true per-respondent
# betas (so unavailable alternatives are never chosen and the truth is known),
# with alternative 1 made unavailable for ~40% of respondents (a person-level
# "no car"-style restriction). Both estimators fit the SAME data with the SAME
# availability and are compared on betas, classical/robust SEs, and recovery.
#
# Seeded, self-contained. Run:  Rscript dev/benchmark_availability.R
# Output: console tables + output/benchmark_availability.rds
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

C <- 2; seed <- 101
DGP <- klue_dgp(n_generic = 4, n_alternatives = 3)
J <- DGP$n_alternatives; NG <- DGP$n_generic; NB <- DGP$n_beta; NASC <- DGP$n_asc

d  <- klue_simulate(N_per_class = 300, T_tasks = 14, true_K = C, separation = 1.6,
                    heterogeneity = 0.2, seed = seed, dgp = DGP)
db <- d$database
N <- d$N; Tt <- d$T
resp <- match(db$ID, unique(db$ID))            # 1..N respondent index per row
ib <- d$individual_betas                       # N x NB true per-respondent betas

# ---- assign person-level availability: alt 1 unavailable for ~40% of people --
set.seed(seed + 1)
no_alt1 <- which(runif(N) < 0.4)
av <- matrix(1L, nrow(db), J)
av[resp %in% no_alt1, 1] <- 0L

# ---- re-draw choices from the availability-masked MNL (true betas) -----------
V <- matrix(0, nrow(db), J)
for (j in 1:J) {
  Xj <- as.matrix(db[, c(paste0("x", 1:NG, "_", j), paste0("price_", j))])
  V[, j] <- rowSums(Xj * ib[resp, ])
}
V[av == 0] <- -Inf
P <- exp(V - apply(V, 1, max)); P <- P / rowSums(P)
set.seed(seed + 2)
db$CHOICE <- apply(P, 1, function(pr) sample.int(J, 1, prob = pr))
for (j in 1:J) db[[paste0("av_", j)]] <- av[, j]
cat(sprintf("Restricted choice sets: %.0f%% of respondents lack alt 1 (%.0f%% of tasks)\n",
            100 * length(no_alt1) / N, 100 * mean(av[, 1] == 0)))

# ---- Apollo LC spec with availability ----------------------------------------
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
  avl <- paste(sprintf("alt%d = av_%d", 1:J, 1:J), collapse = ", ")   # real availability
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

# ---- klue fit ----------------------------------------------------------------
k <- klue_lcmnl(db, C, dgp = DGP)
k_se <- sqrt(diag(k$vcov)); k_rse <- sqrt(diag(k$robust_vcov))

# ---- Apollo fit --------------------------------------------------------------
start <- c()
for (cc in 1:C) { for (jj in 1:NASC) start[paste0("asc_alt", jj, "_", cc)] <- 0
  for (a in 1:NG) start[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)
  start[paste0("b_price_", cc)] <- -0.5 }
for (cc in 1:C) start[paste0("delta_", cc)] <- 0
cleanup_apollo()
apollo_control <<- list(modelName = "avail", modelDescr = "avail", indivID = "ID",
                        nCores = 1L, outputDirectory = tempdir())
apollo_beta <<- start; apollo_fixed <<- "delta_1"
apollo_lcPars <<- make_lcPars(C); apollo_probabilities <<- make_probs(C)
apollo_inputs <<- apollo_validateInputs(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                    database = db, apollo_control = apollo_control, silent = TRUE)
m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                     estimate_settings = list(writeIter = FALSE, silent = TRUE, bootstrapSE = 0))
a_se <- sqrt(diag(m$varcov)); a_rse <- sqrt(diag(m$robvarcov))

# ---- align + compare ---------------------------------------------------------
kb <- k$betas
ab <- t(vapply(1:C, function(cc)
         c(vapply(1:NG, function(a) m$estimate[[paste0("b_x", a, "_", cc)]], 0),
           m$estimate[[paste0("b_price_", cc)]]), numeric(NB)))
p <- best_align(kb, ab)
tp <- best_align(d$true_betas, kb)   # klue -> truth alignment for recovery

lab <- c(paste0("x", 1:NG), "price", paste0("asc", 1:NASC))
rows <- list()
for (cc in 1:C) { ac <- p[cc]
  kn <- c(paste0("b", 1:NB, "_class", cc), paste0("asc", 1:NASC, "_class", cc))
  an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac), paste0("asc_alt", 1:NASC, "_", ac))
  for (i in seq_along(kn)) rows[[length(rows) + 1]] <- data.frame(
    class = cc, term = lab[i],
    klue_est = k$par[[kn[i]]], apollo_est = m$estimate[[an[i]]],
    klue_se = k_se[[kn[i]]], apollo_se = a_se[[an[i]]],
    klue_rse = k_rse[[kn[i]]], apollo_rse = a_rse[[an[i]]]) }
tab <- do.call(rbind, rows)

rnd <- function(df) { nc <- vapply(df, is.numeric, logical(1)); df[nc] <- round(df[nc], 4); df }
cat(sprintf("\nLL: klue = %.4f   apollo = %.4f   dLL = %+.2e\n", k$LL, m$maximum, m$maximum - k$LL))
print(rnd(tab), row.names = FALSE)
reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))
# recovery: klue class betas vs true (taste only)
rec <- max(abs(kb[tp, , drop = FALSE] - d$true_betas))
cat(sprintf("\nklue vs Apollo : max|dEst|=%.1e  classical SE reldiff=%.1e  robust SE reldiff=%.1e\n",
            max(abs(tab$klue_est - tab$apollo_est)),
            reldiff(tab$klue_se, tab$apollo_se), reldiff(tab$klue_rse, tab$apollo_rse)))
cat(sprintf("recovery of true class betas (taste): max|d| = %.3f\n", rec))

if (!dir.exists("output")) dir.create("output")
saveRDS(list(tab = tab, LL_klue = k$LL, LL_apollo = m$maximum, recovery = rec),
        "output/benchmark_availability.rds")
cat("Saved: output/benchmark_availability.rds\n")
