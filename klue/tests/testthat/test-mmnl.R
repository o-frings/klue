# klue_mmnl specification builders (no estimation) and small estimation checks
# (Apollo, few draws). The generated code must be the model of Apollo's MMNL
# examples (mu + sigma * draws, sigma unconstrained; -exp() log-normal price;
# written-out Cholesky); random constants follow asc_map; a correlated fit
# started from an independent one cannot end below it; convergence is Apollo's.

gen_code <- function(fn) paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n")

test_that("default specification generates the Apollo example parameterisation", {
  dgp  <- klue_dgp(n_generic = 2, n_alternatives = 3)
  spec <- klue:::.mmnl_spec(dgp)
  expect_identical(spec$random, c("x1", "x2", "price"))
  expect_true(spec$lognormal)
  expect_identical(spec$asc_fixed, c("asc_alt1", "asc_alt2"))
  rc <- gen_code(klue:::.make_apollo_randCoeff(dgp, FALSE, spec))
  expect_match(rc, "mu_x1 + sigma_x1 * draws_x1", fixed = TRUE)
  expect_match(rc, "-exp(mu_price + sigma_price * draws_price)", fixed = TRUE)
  expect_false(grepl("exp(sigma", rc, fixed = TRUE))
  expect_identical(environment(klue:::.make_apollo_randCoeff(dgp, FALSE, spec)), globalenv())
  rcc <- gen_code(klue:::.make_apollo_randCoeff(dgp, TRUE, spec))
  expect_match(rcc, "mu_x2 + s_x2_1 * draws_x1 + s_x2_2 * draws_x2", fixed = TRUE)
  expect_match(rcc, "-exp(mu_price + s_pr_1 * draws_x1 + s_pr_2 * draws_x2 + s_pr_3 * draws_price)",
               fixed = TRUE)
  pr <- gen_code(klue:::.make_apollo_prob_mmnl(dgp, FALSE, spec))
  expect_match(pr, "asc_alt1 + b_x1 * x1_1", fixed = TRUE)
  expect_match(pr, 'V[["alt3"]] <- b_x1 * x1_3', fixed = TRUE)
  expect_false(grepl("Vc", pr, fixed = TRUE))            # no availability: untouched
  pa <- gen_code(klue:::.make_apollo_prob_mmnl(dgp, TRUE, spec))
  expect_match(pa, 'Vc <- (CHOICE == 1) * V[["alt1"]] + (CHOICE == 2) * V[["alt2"]]', fixed = TRUE)
  expect_match(pa, 'V[["alt3"]] <- V[["alt3"]] - Vc', fixed = TRUE)
  expect_match(pa, "alt3 = av_3", fixed = TRUE)
})

test_that("shared and random constants follow asc_map; fixed coefficients get no spread", {
  dgp  <- klue_dgp(n_generic = 3, n_alternatives = 4, asc_map = c(1, 1, 2, 0))
  spec <- klue:::.mmnl_spec(dgp, random = c("asc1", "x1", "x3"), price = "fixed")
  expect_identical(spec$random, c("x1", "x3", "asc1"))   # canonical order
  expect_false(spec$lognormal)
  expect_identical(spec$asc_fixed, c("asc1", "asc2"))     # non-default map
  pr <- gen_code(klue:::.make_apollo_prob_mmnl(dgp, FALSE, spec))
  expect_match(pr, 'V[["alt1"]] <- b_asc1 + b_x1', fixed = TRUE)
  expect_match(pr, 'V[["alt2"]] <- b_asc1 + b_x1', fixed = TRUE)   # shared constant
  expect_match(pr, 'V[["alt3"]] <- asc2 + b_x1', fixed = TRUE)     # fixed constant
  expect_match(pr, 'V[["alt4"]] <- b_x1 * x1_4', fixed = TRUE)     # reference
  rc <- gen_code(klue:::.make_apollo_randCoeff(dgp, FALSE, spec))
  expect_false(grepl("sigma_x2|sigma_price|b_price", rc))          # fixed: no draws
  expect_match(rc, "mu_asc1 + sigma_asc1 * draws_asc1", fixed = TRUE)
  expect_error(klue:::.mmnl_spec(dgp, random = c("price"), price = "fixed"), "conflicts")
  expect_error(klue:::.mmnl_spec(dgp, random = c("x9")), "unknown")
})

test_that("independent-to-correlated mapping copies the signed standard deviations onto the diagonal", {
  dgp  <- klue_dgp(n_generic = 2, n_alternatives = 3)
  spec <- klue:::.mmnl_spec(dgp, random = c("x1", "x2", "price", "asc1"))
  ind <- c(asc_alt2 = 0.3, mu_x1 = 1, mu_x2 = -0.5, mu_price = 0.2, mu_asc1 = 0.7,
           sigma_x1 = 0.4, sigma_x2 = -0.2, sigma_price = 0.3, sigma_asc1 = 1.5)
  b0 <- c(asc_alt2 = 0, mu_x1 = 0, mu_x2 = 0, mu_price = 0, mu_asc1 = 0)
  m <- klue:::.mmnl_indep_to_corr(ind, b0, spec)
  expect_equal(unname(m[c("s_x1_1", "s_x2_2", "s_pr_3", "s_asc1_4")]), c(0.4, -0.2, 0.3, 1.5))
  expect_equal(unname(m[c("s_x2_1", "s_pr_1", "s_asc1_3")]), c(0, 0, 0))
  expect_equal(unname(m[c("asc_alt2", "mu_x1", "mu_asc1")]), c(0.3, 1, 0.7))
})

test_that("random-constant and warm-started correlated fits estimate (Apollo)", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 60, T_tasks = 6, true_K = 1, heterogeneity = 0.4,
                       seed = 7, dgp = dgp)$database
  common <- list(database = db, n_draws = 60, n_cores = 1, quiet = TRUE, dgp = dgp)
  ind <- do.call(klue_mmnl, c(common, list(random = c("x1", "x2", "price", "asc1"))))
  expect_true(ind$converged)
  expect_true(ind$apollo_status$successfulEstimation)   # convergence is Apollo's
  expect_identical(ind$sd_scale, "linear")
  expect_identical(ind$start_used, "pooled_mnl")
  expect_equal(unname(ind$sigma), unname(abs(ind$par[paste0("sigma_", ind$random)])))
  expect_equal(unname(ind$start_values[paste0("sigma_", ind$random)]), rep(0.01, 4))
  expect_true(all(c("mu_asc1", "sigma_asc1", "asc_alt2") %in% names(ind$par)))
  expect_identical(rownames(ind$robust_vcov), names(ind$par))
  # Exact nesting: with the same draws, the correlated model at the mapped
  # independent point reproduces the independent log-likelihood.
  spec <- klue:::.mmnl_spec(dgp, ind$random)
  b0 <- klue:::.mmnl_indep_to_corr(ind$par, ind$par[!grepl("^sigma_", names(ind$par))], spec)
  invisible(utils::capture.output(ap <- suppressMessages(
    klue:::.run_apollo_mmnl(db, 60, b0, dgp = dgp, n_cores = 1, correlation = TRUE,
                            spec = spec))))
  expect_equal(unname(ap$LLStart), unname(ind$LL), tolerance = 1e-8)
  lc <- klue_lincom(ind, c(mu_x1 = 1, mu_x2 = -1))
  expect_true(is.finite(lc[["estimate"]]))
  cor <- do.call(klue_mmnl, c(common, list(random = ind$random, correlation = TRUE,
                                           warm_start = ind)))
  expect_true(cor$converged)
  expect_identical(cor$start_used, "warm_start")
  expect_equal(cor$k, 1 + 4 + 10)                   # asc_alt2, 4 means, 10 Cholesky
  expect_true(isTRUE(cor$nested_ok))
  expect_gte(cor$LL, ind$LL - 1e-6 * abs(ind$LL))
  fx <- do.call(klue_mmnl, c(common, list(random = c("x1", "asc1"), price = "fixed")))
  expect_true(fx$converged)
  expect_true(all(c("b_x2", "b_price", "mu_x1", "mu_asc1") %in% names(fx$par)))
  expect_false(any(c("sigma_x2", "sigma_price") %in% names(fx$par)))
})

test_that("deprecated two-stage and bounds arguments warn; old-scale warm starts are refused", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 40, T_tasks = 5, true_K = 1, heterogeneity = 0.4,
                       seed = 3, dgp = dgp)$database
  expect_warning(f <- klue_mmnl(db, n_draws = 30, n_cores = 1, quiet = TRUE, dgp = dgp,
                                n_draws_stage1 = 10),
                 "deprecated and ignored")
  expect_false("n_draws_stage1" %in% names(f$settings))
  expect_warning(klue_mmnl(db, n_draws = 30, n_cores = 1, quiet = TRUE, dgp = dgp,
                           mu_price_bounds = c(-5, 3)),
                 "deprecated and ignored")
  old <- f; old$sd_scale <- NULL                       # a klue < 0.10.0 fit
  expect_error(klue_mmnl(db, correlation = TRUE, n_draws = 30, n_cores = 1,
                         quiet = TRUE, dgp = dgp, warm_start = old),
               "klue < 0.10.0")
})

test_that("the correlated informed start is an independent fit with the same draws", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 40, T_tasks = 5, true_K = 1, heterogeneity = 0.4,
                       seed = 5, dgp = dgp)$database
  cor <- klue_mmnl(db, correlation = TRUE, n_draws = 40, n_cores = 1, quiet = TRUE,
                   dgp = dgp)
  expect_true(cor$converged)
  expect_identical(cor$start_used, "independent_fit")
  expect_identical(cor$independent_fit$settings$n_draws, cor$settings$n_draws)
  expect_true(isTRUE(cor$nested_ok))
  expect_gte(cor$LL, cor$independent_fit$LL - 1e-6 * abs(cor$independent_fit$LL))
})

test_that("warm starts must match the price distribution and specification", {
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 20, T_tasks = 4, true_K = 1, heterogeneity = 0.4,
                       seed = 9, dgp = dgp)$database
  fake <- list(converged = TRUE, sd_scale = "linear", correlation = FALSE,
               random = c("x1", "x2", "price"), price = "lognormal",
               par = c(asc_alt1 = 0, asc_alt2 = 0, mu_x1 = 1, mu_x2 = 1, mu_price = 0,
                       sigma_x1 = 0.1, sigma_x2 = 0.1, sigma_price = 0.1),
               settings = list(n_draws = 30, draws_type = "mlhs"))
  skip_if_not_installed("apollo")
  expect_error(klue_mmnl(db, correlation = TRUE, price = "normal", n_draws = 30,
                         quiet = TRUE, dgp = dgp, warm_start = fake), "price")
  bad <- fake; bad$par <- bad$par[names(bad$par) != "asc_alt2"]
  expect_error(klue_mmnl(db, correlation = TRUE, n_draws = 30, quiet = TRUE,
                         dgp = dgp, warm_start = bad), "lacks parameter")
})

test_that("normal price with a neutral start, and a shared random constant, estimate (Apollo)", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 60, T_tasks = 6, true_K = 1, heterogeneity = 0.4,
                       seed = 11, dgp = dgp)$database
  f <- klue_mmnl(db, n_draws = 40, n_cores = 1, quiet = TRUE, dgp = dgp,
                 price = "normal", start = "zero")
  expect_true(f$converged)
  expect_identical(f$start_used, "zero")
  expect_equal(unname(f$start_values[c("mu_x1", "mu_x2", "mu_price")]), c(0, 0, 0))
  expect_true(is.finite(f$BIC_apollo))
  dgp2 <- klue_dgp(n_generic = 2, n_alternatives = 3, asc_map = c(1, 1, 0))
  db2  <- klue_simulate(N_per_class = 60, T_tasks = 6, true_K = 1, heterogeneity = 0.4,
                        seed = 12, dgp = dgp2)$database
  g <- klue_mmnl(db2, n_draws = 40, n_cores = 1, quiet = TRUE, dgp = dgp2,
                 random = c("x1", "price", "asc1"))
  expect_true(g$converged)
  expect_true(all(c("mu_asc1", "sigma_asc1", "b_x2") %in% names(g$par)))
})

test_that("failures report Apollo's reason and keep the independent fit; caller state is restored", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 40, T_tasks = 5, true_K = 1, heterogeneity = 0.4,
                       seed = 13, dgp = dgp)$database
  had <- exists("apollo_beta", envir = globalenv(), inherits = FALSE)
  assign("apollo_beta", c(user = 42), envir = globalenv())
  on.exit(if (!had) rm("apollo_beta", envir = globalenv()), add = TRUE)
  set.seed(99); before <- .Random.seed
  ok <- klue_mmnl(db, n_draws = 30, n_cores = 1, quiet = TRUE, dgp = dgp)
  expect_identical(get("apollo_beta", envir = globalenv()), c(user = 42))
  expect_identical(.Random.seed, before)
  expect_true(ok$converged)
  bad <- klue_mmnl(db, correlation = TRUE, n_draws = 30, n_cores = 1, quiet = TRUE,
                   dgp = dgp, estimation_routine = "no_such_routine")
  expect_false(bad$converged)
  expect_match(bad$reason, "^independent_start_fit_failed")
  expect_false(is.null(bad$independent_fit))
  expect_identical(get("apollo_beta", envir = globalenv()), c(user = 42))
})

test_that("klue() leaves a neutral start of the correlated fit alone", {
  skip_on_cran()
  skip_if_not_installed("apollo")
  dgp <- klue_dgp(n_generic = 2, n_alternatives = 3)
  db  <- klue_simulate(N_per_class = 50, T_tasks = 6, true_K = 1, heterogeneity = 0.4,
                       seed = 14, dgp = dgp)$database
  res <- suppressWarnings(klue(database = db, C_cands = 1:2, run_mmnl = TRUE,
                               run_mmnl_corr = TRUE, verbose = FALSE, write_csv = FALSE,
                               mmnl_opts = list(n_draws = 30, n_cores = 1, start = "zero")))
  expect_identical(res$mmnl_corr$start_used, "zero")
  res2 <- suppressWarnings(klue(database = db, C_cands = 1:2, run_mmnl = TRUE,
                                run_mmnl_corr = TRUE, verbose = FALSE, write_csv = FALSE,
                                 mmnl_opts = list(n_draws = 30, n_cores = 1)))
  expect_identical(res2$mmnl_corr$start_used, "warm_start")
  expect_true(isTRUE(res2$mmnl_corr$nested_ok))
})
