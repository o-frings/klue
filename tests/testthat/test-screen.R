# Multistart screen opt-out: screen = FALSE must reproduce the pre-0.9.1
# path exactly (every start fit at full precision, winner taken as-is, no
# polish), and options(klue.screen = FALSE) must set the default.

fx <- function() klue_simulate(N_per_class = 50, T_tasks = 8, true_K = 2, seed = 11)

test_that("screen = FALSE reproduces the manual full-precision multistart exactly", {
  d <- fx(); dgp <- klue_dgp()

  fit <- suppressMessages(klue_lcmnl(d$database, C = 2, dgp = dgp, screen = FALSE))
  expect_true(fit$converged)

  # the pre-0.9.1 contract, spelled out: tight cluster-MNL start fits, each
  # start estimated at the tight default tolerance, highest finite LL wins,
  # first winner on ties
  starts <- klue:::get_all_starts(d$database, C = 2, dgp = dgp, screen = FALSE)
  best_LL <- -Inf; best_par <- NULL; best_nm <- NA_character_
  for (nm in names(starts)) {
    if (is.null(starts[[nm]])) next
    res <- tryCatch(
      suppressMessages(estimate_lcmnl(d$database, 2L,
                                      start_betas = starts[[nm]]$betas,
                                      start_shares = starts[[nm]]$shares,
                                      dgp = dgp)),
      error = function(e) NULL)
    if (is.null(res) || !is.finite(res$LL)) next
    if (res$LL > best_LL) { best_LL <- res$LL; best_par <- res$par; best_nm <- nm }
  }

  expect_identical(fit$LL, best_LL)
  expect_identical(fit$par, best_par)
  expect_identical(fit$best_method, best_nm)
})

test_that("screened and unscreened paths agree on the optimum here", {
  d <- fx(); dgp <- klue_dgp()
  f_scr <- suppressMessages(klue_lcmnl(d$database, C = 2, dgp = dgp, screen = TRUE))
  f_full <- suppressMessages(klue_lcmnl(d$database, C = 2, dgp = dgp, screen = FALSE))
  expect_true(f_scr$converged)
  expect_equal(f_scr$LL, f_full$LL, tolerance = 1e-5)
  expect_equal(f_scr$BIC, f_full$BIC, tolerance = 1e-4)
})

test_that("screen also gates the cluster-MNL start fits", {
  d <- fx(); dgp <- klue_dgp()
  s_loose <- klue:::get_all_starts(d$database, C = 2, dgp = dgp, screen = TRUE)
  s_tight <- klue:::get_all_starts(d$database, C = 2, dgp = dgp, screen = FALSE)
  # same clustering labels either way, so shares agree ...
  expect_identical(lapply(s_loose, `[[`, "shares"),
                   lapply(s_tight, `[[`, "shares"))
  # ... but the loose fits stop earlier, so the coefficients differ
  expect_false(isTRUE(all.equal(lapply(s_loose, `[[`, "betas"),
                                lapply(s_tight, `[[`, "betas"),
                                tolerance = 1e-9)))
})

test_that("options(klue.screen = FALSE) sets the default", {
  d <- fx(); dgp <- klue_dgp()
  old <- options(klue.screen = FALSE)
  on.exit(options(old), add = TRUE)
  f_opt <- suppressMessages(klue_lcmnl(d$database, C = 2, dgp = dgp))
  f_arg <- suppressMessages(klue_lcmnl(d$database, C = 2, dgp = dgp, screen = FALSE))
  expect_identical(f_opt$LL, f_arg$LL)
  expect_identical(f_opt$par, f_arg$par)
})
