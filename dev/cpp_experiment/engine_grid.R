# Exhaustive cpp-vs-r engine equivalence across DGP shapes, C, panel balance,
# and estimator. Run with the SDK libc++ on CPATH so klue's src/ compiles:
#   CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1 \
#     Rscript dev/engine_grid.R
suppressMessages(devtools::load_all("klue", quiet = TRUE))
stopifnot(klue:::.cpp_ok())

fails <- 0L; n <- 0L; worst <- 0
# The log-likelihood (the optimised objective) is the identifiable invariant
# and must match for BOTH estimators in every config. Class coefficients and
# posteriors are identified only when every class is populated, i.e. C <=
# true_K; on over-specified C the empty classes have flat/unbounded
# coefficient directions that differ between engines (and between any two
# 1e-9-perturbed R runs -- verified) without changing the likelihood. So
# betas/posteriors are asserted only for the deterministic ML path at
# identified C; otherwise LL alone gates.
chk <- function(lab, est, identified, dLL, dbet, dpost) {
  n <<- n + 1L; worst <<- max(worst, dLL)
  tight <- est == "ml" && identified
  ok <- dLL < 1e-4 && (!tight || (dbet < 1e-4 && dpost < 1e-4))
  if (!ok) { fails <<- fails + 1L
    cat(sprintf("  [FAIL] %-42s dLL=%.2e dbetas=%.2e dpost=%.2e\n", lab, dLL, dbet, dpost)) }
}

grid <- expand.grid(n_generic = c(2L, 4L), J = c(2L, 3L, 4L),
                    true_K = c(2L, 3L), est = c("ml", "em"),
                    stringsAsFactors = FALSE)
for (gi in seq_len(nrow(grid))) {
  ng <- grid$n_generic[gi]; J <- grid$J[gi]; tK <- grid$true_K[gi]; est <- grid$est[gi]
  dgp <- klue_dgp(n_generic = ng, n_alternatives = J)
  sim <- klue_simulate(N_per_class = 80, T_tasks = 10, true_K = tK,
                       separation = 1.0, heterogeneity = 0.2, seed = 1000 + gi, dgp = dgp)
  db <- sim$database
  fit <- if (est == "em") estimate_lcmnl_em else estimate_lcmnl
  for (C in 1:4) {
    st <- klue_starts(db, C, "kmeans", dgp = dgp)
    r <- fit(db, C, start_betas = st$betas, start_shares = st$shares, dgp = dgp, engine = "r")
    c <- fit(db, C, start_betas = st$betas, start_shares = st$shares, dgp = dgp, engine = "cpp")
    chk(sprintf("ng=%d J=%d K=%d %s C=%d", ng, J, tK, est, C), est, C <= tK,
        abs(r$LL - c$LL), max(abs(r$betas - c$betas)), max(abs(r$posteriors - c$posteriors)))
  }
}

# Unbalanced-class DGP + full klue_lcmnl multistart equivalence.
for (props in list(c(0.6, 0.25, 0.15), c(0.7, 0.3))) {
  tK <- length(props)
  sim <- klue_simulate(N_per_class = 120, T_tasks = 12, true_K = tK, separation = 1.0,
                       heterogeneity = 0.25, seed = 4242, class_proportions = props)
  for (C in 2:4) {
    mr <- klue_lcmnl(sim$database, C, engine = "r")
    mc <- klue_lcmnl(sim$database, C, engine = "cpp")
    # LL-only: best-of-6 multistart can tie at equal LL via different methods.
    chk(sprintf("unbal K=%d multistart C=%d", tK, C), "multistart", FALSE,
        abs(mr$LL - mc$LL), max(abs(mr$betas - mc$betas)), max(abs(mr$posteriors - mc$posteriors)))
  }
}

cat(sprintf("\n%s  %d configs, %d failures, worst dLL=%.2e\n",
            if (fails == 0L) "ALL ENGINE-EQUIVALENT" else "ENGINE MISMATCH", n, fails, worst))
if (fails > 0L) quit(status = 1L)
