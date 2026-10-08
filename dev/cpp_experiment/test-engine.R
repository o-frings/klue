# The optional RcppArmadillo kernel (engine = "cpp") must reproduce the pure-R
# reference kernel (engine = "r"). Skips when the compiled code is not loaded
# (e.g. sourced standalone, or no compiler available). See dev/engine_grid.R
# for the exhaustive cross-DGP sweep.

# Numerical equivalence is judged by ABSOLUTE max difference: gradient
# components span many magnitudes and pass through zero, where expect_equal's
# partly-relative tolerance is spuriously harsh even at ~1e-13 absolute error.
expect_close <- function(a, b, tol = 1e-8) {
  testthat::expect_lt(max(abs(as.numeric(a) - as.numeric(b))), tol)
}

test_that("C++ and R kernels agree at fixed parameter points (1e-8)", {
  skip_if_not(klue:::.cpp_ok(), "compiled kernel not available")
  set.seed(1)
  sim <- klue_simulate(N_per_class = 80, T_tasks = 10, true_K = 3,
                       separation = 1, heterogeneity = 0.2, seed = 42)
  db <- sim$database; dgp <- DGP_DEFAULT
  ctx_r <- klue:::.lcmnl_context(db, dgp, engine = "r")
  ctx_c <- klue:::.lcmnl_context(db, dgp, engine = "cpp")
  npc <- ctx_r$npc

  # Single-class nll/grad, unweighted and weighted.
  for (rep in 1:10) {
    par <- rnorm(npc, 0, 1.2)
    e <- klue:::.mnl_eval(ctx_r, par, probs = TRUE)
    s <- klue:::.mnl_score(ctx_r, e$probs)
    expect_close(-sum(e$tll), cpp_mnl_nll(ctx_c$dptr, par))
    expect_close(c(-colSums(s$beta), -colSums(s$asc)),
                 cpp_mnl_grad(ctx_c$dptr, par))
    rw <- runif(ctx_r$T_total, 0.1, 2)
    expect_close(-sum(rw * e$tll), cpp_mnl_nll(ctx_c$dptr, par, rw))
    expect_close(c(-colSums(rw * s$beta), -colSums(rw * s$asc)),
                 cpp_mnl_grad(ctx_c$dptr, par, rw))
  }

  # LCMNL nll for C = 2, 3, 4.
  panel <- function(par, C) {
    m <- matrix(0, ctx_r$N, C)
    for (ci in 1:C) m[, ci] <- klue:::.panel_sum(
      klue:::.mnl_eval(ctx_r, par[(ci - 1L) * npc + 1:npc])$tll, ctx_r)
    m
  }
  for (C in 2:4) {
    par <- c(rnorm(C * npc), rnorm(C - 1))
    lp <- panel(par, C)
    d <- c(par[(C * npc + 1):(C * npc + C - 1)], 0); dm <- max(d)
    logpi <- d - dm - log(sum(exp(d - dm)))
    lj <- sweep(lp, 2, logpi, "+"); lmx <- apply(lj, 1, max)
    r_nll <- -sum(lmx + log(rowSums(exp(lj - lmx))))
    expect_close(r_nll, cpp_lcmnl_nll(ctx_c$dptr, par, C))
  }
})

test_that("end-to-end ML fit matches across engines at identified C", {
  skip_if_not(klue:::.cpp_ok(), "compiled kernel not available")
  sim <- klue_simulate(N_per_class = 100, T_tasks = 12, true_K = 3,
                       separation = 1, heterogeneity = 0.2, seed = 7)
  db <- sim$database
  for (C in 1:3) {                       # C <= true_K: all classes identified
    st <- klue_starts(db, C, "kmeans")
    r <- estimate_lcmnl(db, C, start_betas = st$betas, start_shares = st$shares,
                        engine = "r")
    cc <- estimate_lcmnl(db, C, start_betas = st$betas, start_shares = st$shares,
                         engine = "cpp")
    expect_close(r$LL, cc$LL, tol = 1e-6)
    expect_close(r$betas, cc$betas, tol = 1e-4)
    expect_close(r$posteriors, cc$posteriors, tol = 1e-4)
  }
})
