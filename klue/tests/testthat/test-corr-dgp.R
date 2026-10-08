# Correlated-attribute DGP (WS2): attr_corr induces cross-attribute correlation
# via a Gaussian copula while keeping marginals symmetric on [-1, 1]. The
# default (attr_corr = NULL) path must be unchanged.

test_that("attr_corr = NULL leaves the DGP byte-identical", {
  a <- klue_simulate(N_per_class = 80, T_tasks = 10, true_K = 2, seed = 123)
  b <- klue_simulate(N_per_class = 80, T_tasks = 10, true_K = 2, seed = 123,
                     attr_corr = NULL)
  expect_identical(a$database, b$database)
})

test_that("attr_corr induces the requested cross-attribute correlation", {
  d <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 3, seed = 7,
                     attr_corr = 0.6)
  X <- as.matrix(d$database[, c("x1_1", "x2_1", "x3_1", "x4_1")])
  mean_off <- mean(cor(X)[upper.tri(cor(X))])
  # Gaussian copula attenuates slightly; expect comfortably positive, near 0.6
  expect_gt(mean_off, 0.45)
  expect_lt(mean_off, 0.75)
  # marginals stay symmetric on [-1, 1] (means ~ 0)
  expect_lt(max(abs(colMeans(X))), 0.05)
  expect_lte(max(X), 1)
  expect_gte(min(X), -1)
})

test_that("independent DGP is uncorrelated across attributes", {
  d <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 3, seed = 7)
  X <- as.matrix(d$database[, c("x1_1", "x2_1", "x3_1", "x4_1")])
  expect_lt(abs(mean(cor(X)[upper.tri(cor(X))])), 0.05)
})

test_that("estimation converges under the correlated DGP", {
  d <- klue_simulate(N_per_class = 120, T_tasks = 12, true_K = 3, seed = 7,
                     attr_corr = 0.6)
  fit <- klue_lcmnl(d$database, 3)
  expect_true(fit$converged)
})
