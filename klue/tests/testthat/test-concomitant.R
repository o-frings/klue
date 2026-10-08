# Concomitant (covariate-driven) class membership: reduction to the constant-
# share model, well-formed covariance, and the anti-drift guard tying the
# observed information to the estimation likelihood. External agreement (betas,
# gamma, and SEs vs Apollo) is checked in dev/benchmark_concomitant.R.

test_that("no covariates reproduces the constant-share LCMNL exactly", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 200, T_tasks = 10, true_K = 2,
                     separation = 2.0, heterogeneity = 0.1, seed = 101, dgp = dgp)
  st <- klue_starts(d$database, 2, "kmeans", dgp = dgp)
  f_const <- estimate_lcmnl(d$database, 2, start_betas = st$betas, dgp = dgp)
  f_cov0  <- estimate_lcmnl_cov(d$database, 2, membership = character(0),
                                start_betas = st$betas, dgp = dgp)
  expect_equal(f_cov0$LL, f_const$LL, tolerance = 1e-8)
  expect_equal(unname(f_cov0$betas), unname(f_const$betas), tolerance = 1e-6)
  expect_equal(unname(f_cov0$class_probs), unname(f_const$class_probs),
               tolerance = 1e-6)
})

test_that("concomitant fit returns gamma and a well-formed covariance", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 250, T_tasks = 10, true_K = 2, separation = 2.0,
                     heterogeneity = 0.1, seed = 7, dgp = dgp,
                     covariates = TRUE, covariate_strength = 1.5)
  f <- estimate_lcmnl_cov(d$database, 2, membership = c("Z1", "Z2"), dgp = dgp)
  n_free <- 2 * dgp$npc + 1 * 3          # 2 classes taste + 1 non-ref x (intercept+2)

  expect_true(f$converged)
  expect_equal(f$k, n_free)
  expect_equal(dim(f$gamma), c(3L, 1L))                 # (intercept,Z1,Z2) x (C-1)
  expect_equal(rownames(f$gamma), c("(Intercept)", "Z1", "Z2"))
  expect_equal(dim(f$vcov), c(n_free, n_free))
  expect_equal(dim(f$shares), c(length(unique(d$database$ID)), 2L))
  expect_true(all(eigen(f$vcov, only.values = TRUE)$values > 0))
  expect_true(all(sqrt(diag(f$vcov)) > 0))
  expect_equal(unname(rowSums(f$shares)), rep(1, nrow(f$shares)), tolerance = 1e-8)
})

test_that("concomitant observed information matches the estimation likelihood", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 200, T_tasks = 10, true_K = 2, separation = 2.0,
                     heterogeneity = 0.1, seed = 3, dgp = dgp,
                     covariates = TRUE, covariate_strength = 1.2)
  f <- estimate_lcmnl_cov(d$database, 2, membership = c("Z1", "Z2"), dgp = dgp)
  ctx <- .lcmnl_context(d$database, dgp)
  Z <- .lcmnl_membership_design(d$database, ctx, c("Z1", "Z2"))
  obj <- .lcmnl_cov_objective(ctx, 2, Z)
  expect_equal(-obj$neg_ll(f$par), f$LL, tolerance = 1e-6)
  expect_lt(max(abs(obj$grad_ll(f$par))), 1e-3)   # stationary point
})
