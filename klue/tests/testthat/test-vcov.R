# Covariance of the LCMNL fit: shape, positive-definiteness, the vcov = FALSE
# switch, and an anti-drift guard tying the observed-information likelihood to
# the estimation likelihood. External agreement (vs Apollo) is checked in
# dev/benchmark_lcmnl_se.R.

test_that("C >= 2 fit returns a well-formed covariance", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 150, T_tasks = 10, true_K = 3,
                     separation = 2.0, heterogeneity = 0.1, seed = 202, dgp = dgp)
  f <- klue_lcmnl(d$database, 3, dgp = dgp)
  n_free <- 3 * dgp$npc + 2

  expect_false(is.null(f$vcov))
  expect_false(is.null(f$robust_vcov))
  expect_equal(dim(f$vcov), c(n_free, n_free))
  expect_equal(dim(f$robust_vcov), c(n_free, n_free))
  expect_equal(rownames(f$vcov), names(f$par))
  expect_true(isSymmetric(unname(round(f$vcov, 10))))
  expect_true(all(eigen(f$vcov, only.values = TRUE)$values > 0))  # PD
  se <- sqrt(diag(f$vcov))
  expect_true(all(is.finite(se) & se > 0))
})

test_that("vcov = FALSE skips the covariance", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 120, T_tasks = 8, true_K = 2,
                     separation = 2.0, heterogeneity = 0.1, seed = 11, dgp = dgp)
  st <- klue_starts(d$database, 2, "kmeans", dgp = dgp)
  f0 <- estimate_lcmnl(d$database, 2, start_betas = st$betas, dgp = dgp,
                       vcov = FALSE)
  f1 <- estimate_lcmnl(d$database, 2, start_betas = st$betas, dgp = dgp,
                       vcov = TRUE)
  expect_null(f0$vcov)
  expect_null(f0$robust_vcov)
  expect_false(is.null(f1$vcov))
  expect_equal(f0$LL, f1$LL)          # the flag must not touch the optimum
  expect_equal(unname(f0$par), unname(f1$par))
})

test_that("observed information uses the estimation likelihood (no drift)", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 120, T_tasks = 8, true_K = 2,
                     separation = 2.0, heterogeneity = 0.1, seed = 7, dgp = dgp)
  f <- klue_lcmnl(d$database, 2, dgp = dgp)
  ctx <- .lcmnl_context(d$database, dgp)
  obj <- .lcmnl_objective(ctx, 2)
  # the objective that formed the Hessian must reproduce the reported LL...
  expect_equal(-obj$neg_ll(f$par), f$LL, tolerance = 1e-6)
  # ...and the fit must sit at a stationary point (gradient ~ 0)
  expect_lt(max(abs(obj$grad_ll(f$par))), 1e-3)
})

test_that("C = 1 still returns a covariance", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 150, T_tasks = 10, true_K = 1,
                     separation = 1.0, heterogeneity = 0.1, seed = 3, dgp = dgp)
  f <- klue_lcmnl(d$database, 1, dgp = dgp)
  expect_false(is.null(f$vcov))
  expect_equal(dim(f$vcov), c(dgp$npc, dgp$npc))
  expect_true(all(sqrt(diag(f$vcov)) > 0))
})
