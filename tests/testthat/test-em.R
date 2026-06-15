# EM estimator: maximises the same LCMNL likelihood as the direct-ML estimator,
# so from identical starting values it must never reach a meaningfully lower
# log-likelihood (it may reach a higher one by escaping a local optimum).

make_fixture <- function(true_K = 2, sigma = 0.2, sep = 1.0) {
  klue_simulate(N_per_class = 60, T_tasks = 8, true_K = true_K,
                separation = sep, heterogeneity = sigma, seed = 7)
}

test_that("EM is never worse than direct-ML from identical clustering starts", {
  d <- make_fixture(true_K = 2)
  starts <- klue:::get_all_starts(d$database, C = 2)
  tol <- 1e-3
  any_checked <- FALSE
  for (nm in names(starts)) {
    s <- starts[[nm]]
    if (is.null(s)) next
    ml <- klue:::estimate_lcmnl(d$database, 2,
                                start_betas = s$betas, start_shares = s$shares)
    em <- klue:::estimate_lcmnl_em(d$database, 2,
                                   start_betas = s$betas, start_shares = s$shares)
    expect_true(em$converged)
    # EM LL >= ML LL - tol  (EM allowed to be better, not meaningfully worse)
    expect_gte(em$LL, ml$LL - tol)
    any_checked <- TRUE
  }
  expect_true(any_checked)
})

test_that("EM returns a well-formed result", {
  d <- make_fixture(true_K = 3)
  em <- klue:::estimate_lcmnl_em(d$database, 3)
  dgp <- klue_dgp()
  N <- length(unique(d$database$ID))
  expect_equal(dim(em$betas), c(3L, dgp$n_beta))
  expect_equal(dim(em$posteriors), c(N, 3L))
  expect_equal(length(em$class_probs), 3L)
  expect_equal(sum(em$class_probs), 1, tolerance = 1e-8)
  # posterior rows are proper distributions
  expect_equal(unname(rowSums(em$posteriors)), rep(1, N), tolerance = 1e-6)
  expect_identical(em$estimator, "em")
  expect_gt(em$em_iters, 0L)
})

test_that("EM equals direct-ML at C = 1 (no mixture)", {
  d <- make_fixture(true_K = 2)
  ml <- klue:::estimate_lcmnl(d$database, 1)
  em <- klue:::estimate_lcmnl_em(d$database, 1)
  expect_equal(em$LL, ml$LL, tolerance = 1e-6)
  expect_identical(em$em_iters, 0L)
})

test_that("klue_lcmnl dispatches to the EM estimator", {
  d <- make_fixture(true_K = 2)
  fit <- klue_lcmnl(d$database, 2, estimator = "em")
  expect_true(fit$converged)
  expect_identical(fit$estimator, "em")
  expect_false(is.na(fit$best_method))
})
