# klue_refit, klue_draw_betas, klue_predict, klue_simulate_choices with constants,
# and klue_derror (0.11.0).

.rf_data <- function(seed = 9, N = 150) {
  d <- klue_simulate(N_per_class = N, T_tasks = 8, true_K = 2, separation = 2, heterogeneity = 0.1,
                     seed = seed, dgp = klue_dgp(2, 3))
  d$database$x3_1 <- d$database$x1_1 * d$database$x2_1          # an extra attribute for added-term refits
  d$database$x3_2 <- d$database$x1_2 * d$database$x2_2
  d$database$x3_3 <- d$database$x1_3 * d$database$x2_3
  list(db = d$database, dgp = klue_dgp(2, 3, labels = c("a", "b", "cost")),
       dgp3 = klue_dgp(3, 3, labels = c("a", "b", "ab", "cost")))
}
.sb <- rbind(c(1, 1, -1), c(0, 0, -2))

test_that("a refit on the same data reproduces the fit", {
  d <- .rf_data()
  m <- estimate_lcmnl(d$db, 1, dgp = d$dgp)
  expect_equal(klue_refit(m, d$db, d$dgp)$LL, m$LL, tolerance = 1e-8)
  f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  r <- klue_refit(f, d$db, d$dgp)
  expect_equal(r$LL, f$LL, tolerance = 1e-6)
  expect_false(r$failed)
})

test_that("refits take subsets, added and dropped terms, and start overrides", {
  d <- .rf_data(); f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  sub <- d$db[d$db$ID <= 200, ]
  expect_false(klue_refit(f, sub, d$dgp)$failed)
  s3 <- klue:::.refit_start(f, d$dgp3, 2L)
  expect_equal(s3[c(3, 9)], c(0, 0))                                   # the new term starts at 0 (npc = 6)
  expect_equal(s3[c(1, 2, 4)], unname(f$par[c("a_class1", "b_class1", "cost_class1")]))
  g <- klue_refit(f, d$db, d$dgp3)
  expect_false(g$failed); expect_true("ab_class2" %in% names(g$par))
  dropped <- klue_refit(g, d$db, d$dgp)
  expect_equal(dropped$LL, f$LL, tolerance = 1e-4)
  s <- klue:::.refit_start(f, d$dgp, 2L, start = c(cost = -5, a_class2 = 9))
  expect_equal(s[c(3, 8, 6)], c(-5, -5, 9))
  expect_error(klue:::.refit_start(f, d$dgp, 2L, start = c(zz = 1)), "unknown")
})

test_that("a refit with common coefficients starts them at the class mean", {
  d <- .rf_data(); f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  g <- klue_refit(f, d$db, klue:::.with_common(d$dgp, "cost"))
  expect_equal(g$k, f$k - 1L)
  expect_equal(unname(g$par["cost_class1"]), unname(g$par["cost_class2"]))
  expect_false(g$failed)
})

test_that("a membership refit equals the concomitant model started by hand", {
  d <- .rf_data(); d$db$z <- rep(rnorm(300), each = 8)
  f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  r <- klue_refit(f, d$db, d$dgp, membership = "z")
  h <- estimate_lcmnl_cov(d$db, 2, membership = "z", dgp = d$dgp,
                          start_par = c(unname(f$par[1:10]), unname(f$par["delta1"]), 0))
  expect_equal(r$LL, h$LL, tolerance = 1e-8)
  expect_true("z_class1" %in% names(r$par))
})

test_that("multistart keeps the better clean fit", {
  d <- .rf_data(); f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  m <- suppressMessages(klue_refit(f, d$db, d$dgp, multistart = TRUE))
  expect_gte(m$LL, f$LL - 1e-6)
  expect_false(m$failed)
})

test_that("draws and predictions follow the fit; constants enter through the ASC block", {
  d <- .rf_data(); f <- estimate_lcmnl(d$db, 2, start_betas = .sb, dgp = d$dgp)
  set.seed(4); B <- klue_draw_betas(f, 20000)
  pi <- klue:::.klue_shares(f$par, grep("^delta", names(f$par)))
  expect_equal(as.numeric(table(attr(B, "class_draw")) / 20000), pi, tolerance = 0.02)
  expect_equal(unname(B[attr(B, "class_draw") == 2, , drop = FALSE][1, ]), unname(f$par[6:10]))
  P <- klue_predict(f, d$db, d$dgp)
  expect_equal(rowSums(P), rep(1, nrow(d$db)))
  m <- estimate_lcmnl(d$db, 1, dgp = d$dgp)
  V <- klue:::.row_utilities(klue:::.design_full(d$db, d$dgp), matrix(m$par, nrow(d$db), 5, byrow = TRUE))
  expect_equal(klue_predict(m, d$db, d$dgp), exp(V) / rowSums(exp(V)))
  # simulated shares from the class draws match the predicted mixture
  N <- max(d$db$ID); set.seed(5)
  sh <- rowMeans(vapply(1:40, function(r) tabulate(klue_simulate_choices(d$db, klue_draw_betas(f, N), d$dgp)$CHOICE, 3), numeric(3)))
  expect_equal(sh / nrow(d$db), colMeans(P), tolerance = 0.01)
  # the ASC block equals the same utilities computed by hand
  b <- c(0.5, -0.3, -1, 0.8, -0.4); set.seed(11)
  s1 <- klue_simulate_choices(d$db, b, d$dgp)$CHOICE
  set.seed(11); Vb <- klue:::.row_utilities(klue:::.design_full(d$db, d$dgp), matrix(b, nrow(d$db), 5, byrow = TRUE))
  expect_equal(s1, max.col(Vb - log(-log(matrix(runif(length(Vb)), nrow(Vb))))))
})

test_that("the D-error matches its closed form, averages over draws, and flags aliasing", {
  dgp <- klue_dgp(1, 3, asc_map = c(1, 0, 0))
  db <- data.frame(ID = 1:4, CHOICE = 1, x1_1 = c(1, 0, 1, 0), x1_2 = c(0, 1, 1, 0), x1_3 = 0,
                   price_1 = c(1, 2, 3, 1), price_2 = c(2, 1, 1, 3), price_3 = 0)
  Xf <- klue:::.design_full(db, dgp)
  I <- Reduce(`+`, lapply(1:4, function(t) { X <- t(sapply(Xf, function(Xj) Xj[t, ])); Z <- sweep(X, 2, colMeans(X)); crossprod(Z) / 3 }))
  r <- klue_derror(db, dgp)
  expect_equal(r$derror, det(I)^(-1 / 3))
  expect_true(r$full_rank); expect_equal(r$K, 3L)
  dr <- rbind(c(0.2, -0.5, 0.1), c(-0.1, -1, 0.3))
  expect_equal(klue_derror(db, dgp, draws = dr)$derror,
               mean(c(klue_derror(db, dgp, par = dr[1, ])$derror, klue_derror(db, dgp, par = dr[2, ])$derror)))
  db$x1_2 <- db$x1_1                                   # attribute 1 now has no variation within a task across A, B
  db$x1_3 <- db$x1_1
  expect_false(klue_derror(db, dgp)$full_rank)
})
