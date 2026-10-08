# Delta-method inference on functions of the parameters (0.11.0): klue_nlcom,
# klue_wald(fn = ), klue_fieller and klue_tost with functions, and the population
# means of klue_mean_fn and klue_wtp_fn.

.lab_fits <- function() {
  dgp <- klue_dgp(2, 3, labels = c("a", "b", "cost"), asc_labels = c("k1", "k2"))
  d <- klue_simulate(N_per_class = 150, T_tasks = 8, true_K = 2, separation = 2, heterogeneity = 0.1,
                     seed = 9, dgp = klue_dgp(2, 3))
  list(mnl = estimate_lcmnl(d$database, 1, dgp = dgp),
       lc = estimate_lcmnl(d$database, 2, start_betas = rbind(c(1, 1, -1), c(0, 0, -2)), dgp = dgp),
       lcc = estimate_lcmnl(d$database, 2, start_betas = rbind(c(1, 1, -1), c(0, 0, -2)),
                            dgp = klue:::.with_common(dgp, "cost")))
}

test_that("a linear function reproduces the linear tests", {
  f <- .lab_fits()$lc
  w <- c(a_class1 = 1, b_class2 = -0.5)
  lin <- function(p) sum(w * p[names(w)])
  expect_equal(klue_nlcom(f, lin)[c("estimate", "se", "z", "p")], klue_lincom(f, w), tolerance = 1e-8)
  expect_equal(klue_tost(f, lin, 0.5), klue_tost(f, w, 0.5), tolerance = 1e-8)
  expect_equal(klue_wald(f, fn = function(p) p[c("a_class1", "b_class1")]),
               klue_wald(f, c("a_class1", "b_class1")), tolerance = 1e-7)
  num <- c(a_class1 = 1); den <- c(cost_class1 = 1)
  expect_equal(klue_fieller(f, function(p) p[["a_class1"]], function(p) p[["cost_class1"]]),
               klue_fieller(f, num, den), tolerance = 1e-7)
  expect_equal(klue_fieller(f, num, function(p) p[["cost_class1"]]), klue_fieller(f, num, den), tolerance = 1e-7)
})

test_that("a nonlinear function matches its closed-form delta method", {
  f <- .lab_fits()$mnl
  r <- klue_nlcom(f, function(p) exp(p[["a"]]))
  expect_equal(unname(r["estimate"]), exp(unname(f$par["a"])))
  expect_equal(unname(r["se"]), exp(unname(f$par["a"])) * sqrt(f$robust_vcov["a", "a"]), tolerance = 1e-6)
  expect_error(klue_nlcom(f, function(p) p[1:2]), "one number")
})

test_that("a Wald test on a redundant vector has the covariance's rank as df", {
  f <- .lab_fits()$mnl
  w3 <- klue_wald(f, fn = function(p) c(p[["a"]], p[["a"]], p[["b"]]))
  w2 <- klue_wald(f, c("a", "b"))
  expect_equal(unname(w3["df"]), 2)
  expect_equal(unname(w3["chisq"]), unname(w2["chisq"]), tolerance = 1e-6)
})

test_that("population means: shares times class values, and the MNL as the C = 1 case", {
  F <- .lab_fits(); f <- F$lc; m <- F$mnl
  p <- f$par; pi <- klue:::.klue_shares(p, grep("^delta", names(p)))
  expect_equal(klue_mean_fn(f, "a")(p), sum(pi * p[c("a_class1", "a_class2")]))
  expect_equal(klue_mean_fn(m, "a")(m$par), unname(m$par["a"]))
  ratio <- -100 * p[c("a_class1", "a_class2")] / p[c("cost_class1", "cost_class2")]
  expect_equal(klue_wtp_fn(f, c(a = 1), scale = 100)(p), sum(pi * ratio))
  expect_equal(klue_wtp_fn(f, c(a = 1), scale = 100, class = 2)(p), unname(ratio[2]))
  expect_equal(klue_nlcom(m, klue_wtp_fn(m, c(a = 1)))[c("estimate", "se")],
               setNames(klue_wtp(m, c(a = 1), c(cost = 1))[c("wtp", "se")], c("estimate", "se")), tolerance = 1e-7)
  s2 <- klue_wtp_fn(f, c(k2 = 1, k1 = -1))(p)
  expect_equal(s2, sum(pi * -(p[c("k2_class1", "k2_class2")] - p[c("k1_class1", "k1_class2")]) /
                        p[c("cost_class1", "cost_class2")]))
})

test_that("a common price enters every class with its one value; a ratio of means takes Fieller", {
  f <- .lab_fits()$lcc; p <- f$par
  expect_equal(unname(p["cost_class1"]), unname(p["cost_class2"]))
  pi <- klue:::.klue_shares(p, grep("^delta", names(p)))
  expect_equal(klue_wtp_fn(f, c(b = 1))(p), sum(pi * -p[c("b_class1", "b_class2")] / p[["cost_class1"]]))
  r <- klue_fieller(f, klue_mean_fn(f, "a"), klue_mean_fn(f, "cost"))
  expect_equal(unname(r["ratio"]), klue_mean_fn(f, "a")(p) / klue_mean_fn(f, "cost")(p))
  expect_true(all(is.finite(r)) && r["lo"] < r["ratio"] && r["ratio"] < r["hi"])
})
