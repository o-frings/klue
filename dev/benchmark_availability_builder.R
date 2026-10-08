# dev/benchmark_availability_builder.R
# ------------------------------------------------------------------------------
# End-to-end benchmark of the availability pipeline: raw long data with an
# availability column -> klue_database(keep_unavailable = TRUE) (which now emits
# av_1..av_J instead of dropping tasks) -> klue_lcmnl -> compare to Apollo on the
# same built database. Availability-consistent choices are re-drawn from the
# masked MNL with known true betas (alt 1 unavailable for ~40% of respondents).
#
# Also shows the data retained vs the old drop-everything-unavailable behaviour.
#
# Seeded, self-contained. Run:  Rscript dev/benchmark_availability_builder.R
# Output: console + output/benchmark_availability_builder.rds
# ------------------------------------------------------------------------------

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))

C <- 2; seed <- 101
DGP <- klue_dgp(n_generic = 4, n_alternatives = 3)
J <- DGP$n_alternatives; NG <- DGP$n_generic; NB <- DGP$n_beta; NASC <- DGP$n_asc

d  <- klue_simulate(N_per_class = 300, T_tasks = 14, true_K = C, separation = 1.6,
                    heterogeneity = 0.2, seed = seed, dgp = DGP)
wide <- d$database
N <- d$N; Tt <- d$T
wide$TASK <- ave(seq_len(nrow(wide)), wide$ID, FUN = seq_along)
resp <- match(wide$ID, unique(wide$ID))
ib <- d$individual_betas

# person-level availability + availability-consistent re-drawn choices
set.seed(seed + 1); no_alt1 <- which(runif(N) < 0.4)
av <- matrix(1L, nrow(wide), J); av[resp %in% no_alt1, 1] <- 0L
V <- matrix(0, nrow(wide), J)
for (j in 1:J) V[, j] <- rowSums(as.matrix(wide[, c(paste0("x", 1:NG, "_", j),
                                                    paste0("price_", j))]) * ib[resp, ])
V[av == 0] <- -Inf
P <- exp(V - apply(V, 1, max)); P <- P / rowSums(P)
set.seed(seed + 2)
wide$CHOICE <- apply(P, 1, function(pr) sample.int(J, 1, prob = pr))

# ---- melt to RAW LONG (id, task, alt, choice indicator, attrs, price, avail) --
long <- do.call(rbind, lapply(1:J, function(j) data.frame(
  id = wide$ID, task = wide$TASK, alt = j,
  choice = as.integer(wide$CHOICE == j),
  X1 = wide[[paste0("x1_", j)]], X2 = wide[[paste0("x2_", j)]],
  X3 = wide[[paste0("x3_", j)]], X4 = wide[[paste0("x4_", j)]],
  price = wide[[paste0("price_", j)]], avail = av[, j])))
long <- long[order(long$id, long$task, long$alt), ]

# ---- build canonical db two ways ---------------------------------------------
common <- list(data = long, format = "long", id_col = "id", task_col = "task",
               alt_col = "alt", choice_col = "choice", choice_format = "indicator",
               attribute_cols = c("X1", "X2", "X3", "X4"), price_col = "price",
               avail_col = "avail", verbose = FALSE)
db_drop <- do.call(klue_database, common)                                 # old behaviour
db_keep <- do.call(klue_database, c(common, list(keep_unavailable = TRUE)))
cat(sprintf("Tasks retained: drop-unavailable = %d  vs  keep+mask = %d  (%.0f%% more)\n",
            nrow(db_drop), nrow(db_keep), 100 * (nrow(db_keep) / nrow(db_drop) - 1)))
cat(sprintf("Emitted availability columns: %s ; restricted tasks = %.0f%%\n",
            paste(grep("^av_", names(db_keep), value = TRUE), collapse = ","),
            100 * mean(db_keep$av_1 == 0 | db_keep$av_2 == 0 | db_keep$av_3 == 0)))

# ---- Apollo spec (reads av_j from the built db) ------------------------------
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
  avl <- paste(sprintf("alt%d = av_%d", 1:J, 1:J), collapse = ", ")
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

# ---- fit klue on the BUILT db, and Apollo on the same db ---------------------
k <- klue_lcmnl(db_keep, C, dgp = DGP)
k_se <- sqrt(diag(k$vcov)); k_rse <- sqrt(diag(k$robust_vcov))

start <- c()
for (cc in 1:C) { for (jj in 1:NASC) start[paste0("asc_alt", jj, "_", cc)] <- 0
  for (a in 1:NG) start[paste0("b_x", a, "_", cc)] <- 0.1 * (cc - 1)
  start[paste0("b_price_", cc)] <- -0.5 }
for (cc in 1:C) start[paste0("delta_", cc)] <- 0
cleanup_apollo()
apollo_control <<- list(modelName = "availb", modelDescr = "availb", indivID = "ID",
                        nCores = 1L, outputDirectory = tempdir())
apollo_beta <<- start; apollo_fixed <<- "delta_1"
apollo_lcPars <<- make_lcPars(C); apollo_probabilities <<- make_probs(C)
apollo_inputs <<- apollo_validateInputs(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                    database = db_keep, apollo_control = apollo_control, silent = TRUE)
m <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                     estimate_settings = list(writeIter = FALSE, silent = TRUE, bootstrapSE = 0))
a_se <- sqrt(diag(m$varcov)); a_rse <- sqrt(diag(m$robvarcov))

kb <- k$betas
ab <- t(vapply(1:C, function(cc)
         c(vapply(1:NG, function(a) m$estimate[[paste0("b_x", a, "_", cc)]], 0),
           m$estimate[[paste0("b_price_", cc)]]), numeric(NB)))
p <- best_align(kb, ab)
lab <- c(paste0("x", 1:NG), "price", paste0("asc", 1:NASC))
rows <- list()
for (cc in 1:C) { ac <- p[cc]
  kn <- c(paste0("b", 1:NB, "_class", cc), paste0("asc", 1:NASC, "_class", cc))
  an <- c(paste0("b_x", 1:NG, "_", ac), paste0("b_price_", ac), paste0("asc_alt", 1:NASC, "_", ac))
  for (i in seq_along(kn)) rows[[length(rows) + 1]] <- data.frame(
    class = cc, term = lab[i], klue_est = k$par[[kn[i]]], apollo_est = m$estimate[[an[i]]],
    klue_se = k_se[[kn[i]]], apollo_se = a_se[[an[i]]],
    klue_rse = k_rse[[kn[i]]], apollo_rse = a_rse[[an[i]]]) }
tab <- do.call(rbind, rows)
rnd <- function(df) { nc <- vapply(df, is.numeric, logical(1)); df[nc] <- round(df[nc], 4); df }
reldiff <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-8))

cat(sprintf("\nLL: klue = %.4f   apollo = %.4f   dLL = %+.2e\n", k$LL, m$maximum, m$maximum - k$LL))
print(rnd(tab), row.names = FALSE)
cat(sprintf("\nklue vs Apollo on the built db: max|dEst|=%.1e  classical SE reldiff=%.1e  robust SE reldiff=%.1e\n",
            max(abs(tab$klue_est - tab$apollo_est)),
            reldiff(tab$klue_se, tab$apollo_se), reldiff(tab$klue_rse, tab$apollo_rse)))

if (!dir.exists("output")) dir.create("output")
saveRDS(list(tab = tab, n_drop = nrow(db_drop), n_keep = nrow(db_keep),
             LL_klue = k$LL, LL_apollo = m$maximum), "output/benchmark_availability_builder.rds")
cat("Saved: output/benchmark_availability_builder.rds\n")
