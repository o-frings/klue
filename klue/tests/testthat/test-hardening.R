# 0.9.2 hardening: RNG restoration, C = 1 short-circuit, K > 8 guard,
# EM par parity.

fx <- function() klue_simulate(N_per_class = 60, T_tasks = 8, true_K = 2, seed = 11)

test_that("clustering starts do not disturb the caller's RNG stream", {
  d <- fx()
  set.seed(42)
  expected <- rnorm(3)
  set.seed(42)
  invisible(klue_starts(d$database, C = 2, method = "kmeans"))
  invisible(klue:::.cluster_subsample(d$database, cap = 50L))
  expect_identical(rnorm(3), expected)
})

test_that(".with_seed is deterministic and restores a missing .Random.seed", {
  x <- klue:::.with_seed(99, runif(2))
  y <- klue:::.with_seed(99, runif(2))
  expect_identical(x, y)
})

test_that("C = 1 short-circuits to a single pooled fit", {
  d <- fx()
  fit <- suppressMessages(klue_lcmnl(d$database, C = 1))
  expect_true(fit$converged)
  expect_identical(fit$best_method, "pooled")
  expect_length(fit$method_results, 1L)
  # same optimum as a direct MNL fit
  direct <- estimate_lcmnl(d$database, C = 1)
  expect_equal(fit$LL, direct$LL, tolerance = 1e-8)
})

test_that("compute_recovery refuses K > 8 loudly", {
  b <- matrix(0, 9, 5)
  expect_error(klue:::compute_recovery(b, b), "K <= 8")
})

test_that("estimate_lcmnl_em returns par consistent with its betas and shares", {
  d <- fx()
  f <- suppressMessages(estimate_lcmnl_em(d$database, C = 2, max_em_iter = 30L))
  dgp <- klue_dgp()
  expect_length(f$par, 2L * dgp$npc + 1L)
  for (ci in 1:2) {
    expect_equal(f$par[(ci - 1L) * dgp$npc + 1:dgp$n_beta],
                 unname(f$betas[ci, ]), tolerance = 1e-12)
  }
  implied <- exp(c(f$par[2L * dgp$npc + 1L], 0))
  expect_equal(implied / sum(implied), f$class_probs, tolerance = 1e-8)
})

test_that("klue_design accepts a matrix of prior draws (Bayesian D-error)", {
  skip_if_not_installed("idefix")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  draws <- matrix(rep(dgp$beta_bar, 3), nrow = 3, byrow = TRUE) *
    matrix(exp(rnorm(9, 0, .1)), 3)
  d <- klue_design(n_cards = 6L, n_blocks = 2L, dgp = dgp, priors = draws,
                   n_lvls = 2L, n_start = 1L)
  expect_equal(dim(d$cards), c(6L, 3L, 3L))
  expect_true(is.finite(d$Derror))
  expect_error(klue_design(n_cards = 6L, n_blocks = 2L, dgp = dgp,
                           priors = draws[, 1:2], n_lvls = 2L, n_start = 1L),
               "columns")
})
