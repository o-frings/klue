# 0.9.1 large-N behaviour: clustering subsample cap, two-stage multistart
# (screen keeps maxit-hit fits, winner polished from start_par), fit$par.

fx <- function() klue_simulate(N_per_class = 60, T_tasks = 8, true_K = 2, seed = 11)

test_that(".cluster_subsample caps respondents and keeps panels balanced", {
  d <- fx()
  sub <- klue:::.cluster_subsample(d$database, cap = 50L)
  ids <- unique(sub$ID)
  expect_length(ids, 50L)
  expect_true(all(table(sub$ID) == nrow(d$database) / length(unique(d$database$ID))))
  # deterministic
  sub2 <- klue:::.cluster_subsample(d$database, cap = 50L)
  expect_identical(sub, sub2)
  # no-op below the cap
  expect_identical(klue:::.cluster_subsample(d$database, cap = 1000L), d$database)
})

test_that("estimate_lcmnl returns par and warm-starts from it", {
  d <- fx()
  f <- estimate_lcmnl(d$database, C = 2)
  expect_true(f$converged)
  expect_length(f$par, 2L * klue_dgp()$npc + 1L)
  # warm start at the optimum reproduces it (and converges immediately)
  f2 <- estimate_lcmnl(d$database, C = 2, start_par = f$par)
  expect_true(f2$converged)
  expect_equal(f2$LL, f$LL, tolerance = 1e-6)
})

test_that("maxit-hit fits keep a finite LL with converged = FALSE", {
  d <- fx()
  f <- estimate_lcmnl(d$database, C = 2, maxit = 1L)
  expect_false(f$converged)
  expect_true(is.finite(f$LL))
  expect_length(f$par, 2L * klue_dgp()$npc + 1L)
})

test_that("two-stage klue_lcmnl polishes to the single-stage optimum", {
  d <- fx()
  fits <- suppressMessages(klue_lcmnl(d$database, C = 2))
  expect_true(fits$converged)
  expect_true(fits$best_method %in% klue:::KLUE_CLUSTER_METHODS)
  # the polished result matches a direct full-precision fit from its winner
  direct <- estimate_lcmnl(d$database, C = 2, start_par = fits$par)
  expect_equal(fits$LL, direct$LL, tolerance = 1e-6)
})
