# Rank-deficient designs: the covariance must report NA and name the aliased
# parameters rather than inverting a singular information matrix into negative
# variances. The trap reproduced here is the one found in a real DCE: a
# two-level attribute dummy-coded without a reference level, so that the two
# indicators sum to 1 on every non-status-quo alternative and are therefore
# collinear with the status-quo constant.

make_trap_db <- function(seed = 7) {
  dgp <- klue_dgp(n_generic = 3, n_alternatives = 3)
  d <- klue_simulate(N_per_class = 120, T_tasks = 8, true_K = 2,
                     separation = 1.2, heterogeneity = 0.2, seed = seed, dgp = dgp)
  db <- d$database
  # x2, x3 become the two levels of one attribute; alternative 3 is the
  # status quo (all attributes zero), so x2 + x3 + SQ is constant.
  for (j in 1:2) {
    lvl <- as.numeric(db[[paste0("x2_", j)]] > 0)
    db[[paste0("x2_", j)]] <- lvl
    db[[paste0("x3_", j)]] <- 1 - lvl
  }
  for (a in 1:3) db[[paste0("x", a, "_3")]] <- 0
  db[["price_3"]] <- 0
  list(db = db, dgp = dgp)
}

test_that("an exactly collinear design is reported, not silently inverted", {
  tr <- make_trap_db()
  m <- suppressMessages(estimate_lcmnl(tr$db, 1L, dgp = tr$dgp, vcov = TRUE))
  expect_true(m$singular)
  expect_true(all(is.na(m$vcov)))
  expect_true(all(is.na(m$robust_vcov)))
  expect_true(all(c("b2", "b3", "asc1") %in% m$aliased))
  expect_true(is.finite(m$LL))          # the fit itself is unaffected
})

test_that("the alias is reported per class for C >= 2", {
  tr <- make_trap_db()
  f <- suppressMessages(klue_lcmnl(tr$db, C = 2, dgp = tr$dgp))
  expect_true(f$singular)
  expect_true(all(is.na(f$vcov)))
  expect_true(all(c("b2_class1", "b3_class1", "asc1_class1",
                    "b2_class2", "b3_class2", "asc1_class2") %in% f$aliased))
})

test_that("a well-identified design is not flagged and keeps a PD covariance", {
  dgp <- klue_dgp(4, 3)
  d <- klue_simulate(N_per_class = 120, T_tasks = 8, true_K = 2,
                     separation = 1.5, heterogeneity = 0.15, seed = 11, dgp = dgp)
  f <- klue_lcmnl(d$database, C = 2, dgp = dgp)
  expect_false(isTRUE(f$singular))
  expect_length(f$aliased, 0)
  expect_true(all(eigen(f$vcov, only.values = TRUE)$values > 0))
  expect_true(all(is.finite(sqrt(diag(f$vcov)))))
})
