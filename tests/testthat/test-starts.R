# Starting-value generators for the estimator x start-strategy benchmark:
# MNL-perturbation (gmnl-style) and random-partition (lclogit/Latent GOLD-style).

fx <- function() klue_simulate(N_per_class = 60, T_tasks = 8, true_K = 3, seed = 7)

test_that("get_random_partition_starts returns well-formed, deterministic starts", {
  d <- fx(); dgp <- klue_dgp()
  st <- klue:::get_random_partition_starts(d$database, C = 3, n_starts = 4, seed = 1)
  expect_length(st, 4)
  for (s in st) {
    expect_equal(dim(s$betas), c(3L, dgp$n_beta))
    expect_equal(length(s$shares), 3L)
    expect_equal(sum(s$shares), 1, tolerance = 1e-8)
  }
  # deterministic given the seed
  st2 <- klue:::get_random_partition_starts(d$database, C = 3, n_starts = 4, seed = 1)
  expect_equal(st[[1]]$betas, st2[[1]]$betas)
  # different seeds give different partitions
  st3 <- klue:::get_random_partition_starts(d$database, C = 3, n_starts = 4, seed = 2)
  expect_false(isTRUE(all.equal(st[[1]]$betas, st3[[1]]$betas)))
})

test_that("get_mnl_perturbation_starts keeps the price coefficient negative", {
  d <- fx(); dgp <- klue_dgp()
  st <- klue:::get_mnl_perturbation_starts(d$database, C = 3, n_starts = 4, seed = 1)
  expect_length(st, 4)
  for (s in st) {
    expect_equal(dim(s$betas), c(3L, dgp$n_beta))
    expect_true(all(s$betas[, dgp$n_beta] < 0))  # lognormal jitter => price < 0
  }
})
