# Coefficients common to all latent classes (0.11.0): labels, the linear map
# par = M theta, the admissibility check and the fit-status fields.

# Two programmes and an opt-out (reference). Constants are coded as attributes:
# "prog" (1 on both programmes) and "lambda" (1 on the opt-out in arm 1).
# Class 2 never opts out in arm 1 when its lambda is very negative.
.arm_data <- function(lambda2, N = 600L, T = 8L, seed = 1L) {
  set.seed(seed)
  arm <- rep(rep(0:1, length.out = N), each = T)
  db <- data.frame(ID = rep(seq_len(N), each = T), TASK = rep(seq_len(T), N))
  for (j in 1:2) {
    db[[paste0("x1_", j)]] <- sample(0:1, N * T, TRUE)
    db[[paste0("x2_", j)]] <- 1
    db[[paste0("x3_", j)]] <- 0
    db[[paste0("price_", j)]] <- sample(1:5, N * T, TRUE) / 2
  }
  db$x1_3 <- 0; db$x2_3 <- 0; db$x3_3 <- arm; db$price_3 <- 0
  dgp <- klue_dgp(3, 3, asc_map = c(0, 0, 0), labels = c("x", "prog", "lambda", "price"))
  cls <- rep(1:2, each = N / 2)                 # both classes in both arms
  B <- rbind(c(1, 0.5, -0.5, -1), c(-1, 2, lambda2, -0.8))[cls, ]
  db$CHOICE <- 1L
  db$CHOICE <- klue_simulate_choices(db, B, dgp)$CHOICE
  list(db = db, dgp = dgp)
}

test_that("labels rename the parameters and leave the estimates unchanged", {
  d <- klue_simulate(N_per_class = 100, T_tasks = 8, true_K = 2, separation = 2, heterogeneity = 0.1,
                     seed = 3, dgp = klue_dgp(2, 3))
  sb <- rbind(c(1, 1, -1), c(0, 0, -2))
  f0 <- estimate_lcmnl(d$database, 2, start_betas = sb, dgp = klue_dgp(2, 3))
  f1 <- estimate_lcmnl(d$database, 2, start_betas = sb,
                       dgp = klue_dgp(2, 3, labels = c("a", "b", "cost"), asc_labels = c("k1", "k2")))
  expect_equal(unname(f1$par), unname(f0$par), tolerance = 1e-12)
  expect_equal(unname(f1$robust_vcov), unname(f0$robust_vcov), tolerance = 1e-12)
  expect_equal(names(f1$par)[1:5], c("a_class1", "b_class1", "cost_class1", "k1_class1", "k2_class1"))
  m <- estimate_lcmnl(d$database, 1, dgp = klue_dgp(2, 3, labels = c("a", "b", "cost")))
  expect_equal(names(m$par), c("a", "b", "cost", "asc1", "asc2"))
  expect_equal(klue_par(m, "a"), "a")
  expect_equal(klue_par(f1, "a", 2), "a_class2")
  expect_equal(klue_lincom(f1, setNames(1, klue_par(f1, "b", 1))), klue_lincom(f0, c(b2_class1 = 1)))
  expect_error(klue_par(f1, "a"), "class")
})

test_that("a dgp without constants has no constant names", {
  dgp <- klue_dgp(2, 3, asc_map = c(0, 0, 0))
  expect_equal(dgp$par_labels, c("b1", "b2", "b3"))
  d <- .arm_data(-0.5, N = 200)
  f <- estimate_lcmnl(d$db, 1, dgp = d$dgp)
  expect_equal(names(f$par), c("x", "prog", "lambda", "price"))
})

test_that("the analytic gradient with common coefficients matches finite differences", {
  d <- .arm_data(-0.5, N = 200)
  for (C in 2:3) {
    dgp <- klue:::.with_common(d$dgp, c("lambda", "price"))
    ctx <- klue:::.lcmnl_context(d$db, dgp); obj <- klue:::.lcmnl_objective(ctx, C)
    cm <- klue:::.common_map(dgp, C, klue:::.lab("delta", C - 1L))
    set.seed(C); th <- rnorm(ncol(cm$M), sd = 0.3)
    fn <- function(t) obj$neg_ll(drop(cm$M %*% t))
    gr <- drop(crossprod(cm$M, obj$grad_ll(drop(cm$M %*% th))))
    num <- vapply(seq_along(th), function(i) { h <- 1e-6; e <- replace(numeric(length(th)), i, h)
      (fn(th + e) - fn(th - e)) / (2 * h) }, 0)
    expect_equal(gr, num, tolerance = 1e-5)
  }
})

test_that("a common coefficient takes one value, and k falls by (C - 1) per common term", {
  d <- .arm_data(-0.5, N = 400)
  sb <- rbind(c(1, 0.5, -0.5, -1), c(-1, 2, -0.5, -0.8))
  f <- estimate_lcmnl(d$db, 2, start_betas = sb, dgp = d$dgp)
  g <- estimate_lcmnl(d$db, 2, start_betas = sb, dgp = klue_dgp(3, 3, asc_map = c(0, 0, 0),
                      labels = c("x", "prog", "lambda", "price"), common = c("lambda", "price")))
  expect_equal(g$k, f$k - 2L)
  expect_equal(g$common, c("lambda", "price"))
  expect_equal(unname(g$par["lambda_class1"]), unname(g$par["lambda_class2"]))
  expect_equal(unname(g$robust_vcov["price_class1", ]), unname(g$robust_vcov["price_class2", ]))
  expect_false(g$failed)
  expect_equal(klue_par(g, "lambda"), "lambda_class1")
  expect_equal(g$BIC, -2 * g$LL + g$k * log(400))
})

test_that("the admissibility check is one-sided on robust standard errors", {
  fit <- list(failed = FALSE, LL = -1, price = c("p_class1", "p_class2"),
              par = c(p_class1 = -3, p_class2 = -0.5),
              robust_vcov = matrix(c(1, 0, 0, 1), 2, dimnames = rep(list(c("p_class1", "p_class2")), 2)))
  expect_true(klue:::.lc_ok(fit))
  expect_false(klue:::.lc_ok(fit, alpha_cost = 0.05))
  fit$par["p_class2"] <- -1.7
  expect_true(klue:::.lc_ok(fit, alpha_cost = 0.05))
  fit$failed <- TRUE
  expect_false(klue:::.lc_ok(fit))
})

test_that("clean mode keeps a start that has not failed and reports the fit status", {
  d <- klue_simulate(N_per_class = 120, T_tasks = 8, true_K = 2, separation = 2, heterogeneity = 0.1,
                     seed = 5, dgp = klue_dgp(4, 3))
  f <- suppressMessages(klue_lcmnl(d$database, 2, dgp = klue_dgp(4, 3), clean = TRUE))
  expect_false(f$failed)
  expect_gte(f$n_same_optimum, 1L)
  expect_equal(f$entropy_norm, 1 - (f$ICL - f$BIC) / (2 * 240 * log(2)), tolerance = 1e-10)
  expect_equal(f$max_abs_par, max(abs(f$par)))
  expect_true(is.na(estimate_lcmnl(d$database, 1, dgp = klue_dgp(4, 3))$entropy_norm))
})

test_that("the concomitant model takes common coefficients; EM refuses them", {
  d <- .arm_data(-0.5, N = 400)
  d$db$z <- rep(rnorm(400), each = 8)
  sb <- rbind(c(1, 0.5, -0.5, -1), c(-1, 2, -0.5, -0.8))
  f <- estimate_lcmnl_cov(d$db, 2, membership = "z", start_betas = sb, dgp = d$dgp)
  g <- estimate_lcmnl_cov(d$db, 2, membership = "z", start_betas = sb,
                          dgp = klue:::.with_common(d$dgp, "lambda"))
  expect_equal(g$k, f$k - 1L)
  expect_equal(unname(g$par["lambda_class1"]), unname(g$par["lambda_class2"]))
  expect_true(all(is.finite(sqrt(diag(g$robust_vcov)))))
  expect_error(estimate_lcmnl_em(d$db, 2, dgp = klue:::.with_common(d$dgp, "lambda")), "common")
})
