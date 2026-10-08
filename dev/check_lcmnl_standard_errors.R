#!/usr/bin/env Rscript
# =============================================================================
# dev/check_lcmnl_standard_errors.R
#
# Why. The software paragraph of version12.tex (Methods) states, for latent
# class fits with C >= 2: klue takes the standard errors from the inverse of
# the Hessian at the solution, formed by finite differences of the analytic
# gradient with a fixed step, so attributes should be scaled to order one;
# robust (sandwich) standard errors are computed alongside; the classical
# standard errors match those of Apollo 0.3.5 and of an independently coded
# numerical Hessian to within 1e-5 in relative terms; and the robust ones
# differ from Apollo's only by Apollo's factor sqrt(N/(N-1)). This script
# checks each statement. It replaces two scratch checks run on 2026-10-05
# (se_check.R and se_check_scale.R, kept outside git in
# code_backups/scratch_2026-10-06/session_ed878b37/).
#
# Design. Three simulated LCMNLs on klue_dgp(4, 3) (four attributes and a
# price, three alternatives, two constants), T = 10 tasks, attributes of order
# one (uniform on [-1, 1], price uniform on [0.1, 0.9]):
#   A  C = 2, N = 300, separation 2.0, heterogeneity 0.10, data seed 11
#   B  C = 2, N = 300, separation 1.0, heterogeneity 0.25, data seed 12
#   C  C = 3, N = 450, separation 1.5, heterogeneity 0.25, data seed 13
# each fitted by klue_lcmnl() at its defaults (six clustering starts, direct
# ML). Per case:
#   (1) Gradient. klue's analytic gradient against numDeriv::grad of klue's
#       negative log-likelihood and of an independent one (per-class MNL,
#       product over a respondent's tasks, softmax shares with the last class
#       as reference, coded below from the model, not from klue's kernel), at
#       a point perturbed away from the optimum. Max-norm relative difference.
#   (2) Classical SEs. klue's (stats::optimHess on the analytic gradient,
#       absolute step 1e-3, then solve) against the inverse of
#       numDeriv::hessian (Richardson) of klue's and of the independent
#       log-likelihood, and against Apollo 0.3.5's apollo_estimate started at
#       klue's optimum. Measure: the largest |SE_klue - SE_ref| / SE_ref over
#       the parameters; the claim's tolerance is 1e-5. Apollo's model code is
#       that of the ladders' Apollo arms (dev/compare_packages.R builders:
#       apollo_classAlloc, the examples' literal class loop) with parameters
#       listed parameter by parameter and the last class as the share
#       reference, as in klue, so the parameters map one to one. Apollo moves
#       a little from klue's optimum before it computes its covariance, so
#       klue's covariance is also recomputed at Apollo's estimates: that row
#       compares the two methods at the same point.
#   (3) Robust SEs. klue's sandwich against one built here from the
#       independent Hessian and numDeriv per-respondent scores of the
#       independent log-likelihood, and against Apollo's robust SEs before
#       and after multiplying klue's by sqrt(N/(N-1)), N = respondents.
#   (4) Attribute scaling (case A). Price multiplied by s and the optimum
#       mapped exactly (b_price / s), so SE(b_price) must scale by 1/s and
#       every other SE stay put. For each s: the SE error of klue's covariance
#       routine (.lcmnl_vcov, what klue_lcmnl runs at its solution) against
#       its own SEs at s = 1 and against the Richardson reference at s = 1;
#       the same error for a relative-step Hessian (numDeriv::jacobian of the
#       analytic gradient, Apollo 0.3.5's default method) and for
#       numDeriv::hessian; the smallest-to-largest eigenvalue ratio of klue's
#       Hessian; and klue's rank flag (an eigenvalue at or below
#       sqrt(.Machine$double.eps) times the largest sets the covariance to NA,
#       which klue_lcmnl reports as fit$singular = TRUE).
#
# For comparison only: dev/benchmark_lcmnl_se.R (July 2026, an older klue;
# its output is not published) fitted Apollo from its own start and stored a
# largest relative SE difference of 1.1e-5 at C = 3, where the two packages'
# estimates also differ. Here Apollo starts at klue's optimum.
#
# Seeds. Data: klue_simulate seeds 11, 12 and 13 (it calls set.seed itself).
# The perturbed point of (1) is drawn after set.seed(1) (normal, sd 0.05) in
# every case. klue_lcmnl seeds its clustering starts internally; numDeriv,
# optimHess and Apollo's BGW optimiser are deterministic. Same numbers on
# every run; only the run times quoted in the CSV's note column differ.
#
# Run from the repository root, on one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/check_lcmnl_standard_errors.R
# About 35 s on one core (measured 2026-10-06), peak memory about 0.45 GB:
# three small fits, three Apollo estimates, numDeriv Hessians. Stops unless
# the library matches dev/pinned_versions.R (checked when
# dev/compare_packages.R is sourced).
# Writes: output/lcmnl_se_check.csv (long format: one row per part, case and
# measure, with the claim's tolerance where it has one). Apollo writes its
# files to R's tempdir(), removed when R exits.
# =============================================================================

if (!file.exists("dev/compare_packages.R")) stop("run from the repository root")
options(klue.cores = 1L, klue.mmnl.n_cores = 1L, mc.cores = 1L, width = 120)
OUT_CSV <- "output/lcmnl_se_check.csv"
TOL     <- 1e-5                                    # "to within 10^-5 in relative terms"
SCALES  <- c(0.01, 0.1, 10, 100, 1000, 1e4, 3e4, 1e5)   # price multipliers, part (4)

T0 <- proc.time()[["elapsed"]]
stamp <- function(...) { cat(sprintf("[%6.1f s] %s\n", proc.time()[["elapsed"]] - T0,
                                     sprintf(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a)) b else a
fmt_scale <- function(s) vapply(s, format, character(1), scientific = FALSE)

# klue (pkgload::load_all), apollo, the Apollo LC builders and the version
# check; the harness's autorun is guarded, so sourcing runs nothing.
sys.source("dev/compare_packages.R", envir = globalenv())
options(mc.cores = 1L)
suppressMessages(library(numDeriv))
stopifnot(exists(".apollo_lc_setup"), exists(".apollo_est"), exists("cleanup_apollo"))
.lcmnl_context   <- get(".lcmnl_context",   asNamespace("klue"))
.lcmnl_objective <- get(".lcmnl_objective", asNamespace("klue"))
.lcmnl_vcov      <- get(".lcmnl_vcov",      asNamespace("klue"))
VERSIONS <- c(klue = as.character(getNamespaceVersion("klue")),
              apollo = as.character(packageVersion("apollo")),
              numDeriv = as.character(packageVersion("numDeriv")),
              R = paste(R.version$major, R.version$minor, sep = "."))
stamp("klue %s | apollo %s | numDeriv %s | R %s", VERSIONS[["klue"]], VERSIONS[["apollo"]],
      VERSIONS[["numDeriv"]], VERSIONS[["R"]])

DGP <- klue_dgp(4, 3)
NG <- DGP$n_generic; NB <- DGP$n_beta; J <- DGP$n_alternatives
NASC <- DGP$n_asc; NPC <- DGP$npc; ASC_MAP <- DGP$asc_map

CASES <- list(
  A = list(C = 2L, true_K = 2, separation = 2.0, heterogeneity = 0.10, seed = 11,
           what = "C = 2, well separated"),
  B = list(C = 2L, true_K = 2, separation = 1.0, heterogeneity = 0.25, seed = 12,
           what = "C = 2, overlapping"),
  C = list(C = 3L, true_K = 3, separation = 1.5, heterogeneity = 0.25, seed = 13,
           what = "C = 3"))

# ---- helpers -----------------------------------------------------------------
ROWS <- list()
add_row <- function(part, case, measure, value, tolerance = NA_real_, note = "",
                    price_scale = 1) {
  cs <- CASES[[case]]
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    part = part, case = case, C = cs$C, N = cs$N %||% NA_integer_, data_seed = cs$seed,
    price_scale = price_scale, measure = measure, value = as.numeric(value),
    tolerance = tolerance,
    within_tolerance = if (is.na(tolerance) || is.na(value)) NA else value <= tolerance,
    note = note, stringsAsFactors = FALSE)
}
relmax <- function(a, b) max(abs(a - b) / abs(b))

# Independent per-respondent LCMNL log-likelihood, written from the model
# (per-class MNL, product over a respondent's tasks, softmax shares with the
# last class as reference), not from klue's kernel. par layout as klue's:
# class 1 (b1..b5, asc1, asc2), ..., class C, then delta1..delta_{C-1}.
ll_resp <- function(par, db, C) {
  ids <- db$ID; ch <- db$CHOICE; n_obs <- nrow(db)
  logpan <- NULL
  for (cc in 1:C) {
    th <- par[(cc - 1) * NPC + 1:NPC]
    b <- th[1:NB]; a <- if (NASC) th[NB + 1:NASC] else numeric(0)
    U <- matrix(0, n_obs, J)
    for (j in 1:J) {
      u <- b[NB] * db[[paste0("price_", j)]]
      for (k in 1:NG) u <- u + b[k] * db[[paste0("x", k, "_", j)]]
      if (ASC_MAP[j] > 0) u <- u + a[ASC_MAP[j]]
      U[, j] <- u
    }
    m  <- do.call(pmax, as.data.frame(U))
    lp <- U[cbind(seq_len(n_obs), ch)] - (m + log(rowSums(exp(U - m))))
    logpan <- cbind(logpan, as.numeric(rowsum(lp, ids, reorder = FALSE)))
  }
  d   <- c(par[C * NPC + seq_len(C - 1)], 0)
  lpi <- d - max(d); lpi <- lpi - log(sum(exp(lpi)))
  lj  <- sweep(logpan, 2, lpi, "+"); mm <- apply(lj, 1, max)
  mm + log(rowSums(exp(lj - mm)))
}

# klue -> Apollo parameter names (dev/compare_packages.R naming; the last
# class is the share reference in both, the reference alternative's constant
# asc_alt3_c is fixed at 0 in Apollo and absent in klue).
klue_to_apollo <- function(C) {
  k2a <- character(0)
  for (cc in 1:C) {
    for (a in 1:NG) k2a[paste0("b", a, "_class", cc)] <- paste0("b_x", a, "_", cc)
    k2a[paste0("b", NB, "_class", cc)] <- paste0("b_price_", cc)
    for (g in seq_len(NASC)) k2a[paste0("asc", g, "_class", cc)] <- paste0("asc_alt", g, "_", cc)
  }
  for (cc in seq_len(C - 1)) k2a[paste0("delta", cc)] <- paste0("delta_", cc)
  k2a
}

# ---- (1)-(3): one case ----------------------------------------------------------
run_case <- function(lab) {
  cs <- CASES[[lab]]; C <- cs$C
  db <- klue_simulate(N_per_class = 150, T_tasks = 10, true_K = cs$true_K,
                      separation = cs$separation, heterogeneity = cs$heterogeneity,
                      seed = cs$seed, dgp = DGP)$database
  N <- length(unique(db$ID)); CASES[[lab]]$N <<- N
  stamp("case %s: %s (separation %.1f, heterogeneity %.2f, seed %d), N = %d, %d choices",
        lab, cs$what, cs$separation, cs$heterogeneity, cs$seed, N, nrow(db))
  t0 <- proc.time()[["elapsed"]]
  f <- suppressMessages(klue_lcmnl(db, C, dgp = DGP))
  fit_secs <- proc.time()[["elapsed"]] - t0
  if (isTRUE(f$singular) || anyNA(f$vcov)) stop("case ", lab, ": klue's covariance is singular")
  p <- f$par
  ctx <- .lcmnl_context(db, DGP); obj <- .lcmnl_objective(ctx, C)
  ll_indep <- function(x) sum(ll_resp(x, db, C))
  add_row("fit", lab, "klue_LL", f$LL, note = sprintf("klue_lcmnl defaults, %.1f s", fit_secs))
  add_row("fit", lab, "independent_LL_minus_klue_LL", ll_indep(p) - f$LL,
          note = "the independent log-likelihood at klue's optimum")
  add_row("fit", lab, "max_abs_gradient_at_optimum", max(abs(obj$grad_ll(unname(p)))),
          note = "klue's analytic gradient of the negative log-likelihood")

  # (1) analytic gradient at a perturbed point
  set.seed(1); pp <- unname(p) + rnorm(length(p), sd = 0.05)
  g_an  <- obj$grad_ll(pp)
  g_nd  <- numDeriv::grad(obj$neg_ll, pp)
  g_ind <- -numDeriv::grad(ll_indep, pp)
  add_row("gradient", lab, "grad_rel_vs_numDeriv_klue_LL", max(abs(g_an - g_nd)) / max(abs(g_nd)),
          note = "max|g_analytic - g_numDeriv| / max|g_numDeriv| at optimum + N(0, 0.05^2) noise")
  add_row("gradient", lab, "grad_rel_vs_numDeriv_independent_LL",
          max(abs(g_an - g_ind)) / max(abs(g_ind)), note = "as above, independent log-likelihood")

  # (2) classical SEs against Richardson Hessians
  se_k <- sqrt(diag(f$vcov)); rse_k <- sqrt(diag(f$robust_vcov))
  H_k <- numDeriv::hessian(obj$neg_ll, unname(p))
  H_i <- numDeriv::hessian(function(x) -ll_indep(x), unname(p))
  V_i <- solve(H_i); se_i <- sqrt(diag(V_i))
  add_row("classical_se", lab, "se_rel_vs_numDeriv_hessian_klue_LL", relmax(se_k, sqrt(diag(solve(H_k)))),
          TOL, "inverse of numDeriv::hessian (Richardson) of klue's log-likelihood")
  add_row("classical_se", lab, "se_rel_vs_independent_numerical_hessian", relmax(se_k, se_i),
          TOL, "inverse of numDeriv::hessian (Richardson) of the independent log-likelihood")
  Dn <- diag(1 / se_i)
  add_row("classical_se", lab, "vcov_maxdiff_vs_independent_corr_units",
          max(abs(Dn %*% (f$vcov - V_i) %*% Dn)),
          note = "whole matrix, off-diagonals included, scaled by the independent SEs")
  add_row("classical_se", lab, "min_eigenvalue_independent_hessian",
          min(eigen((H_i + t(H_i)) / 2, symmetric = TRUE, only.values = TRUE)$values),
          note = "positive: the optimum is a maximum")

  # (3) robust SEs against an independent sandwich
  S <- numDeriv::jacobian(function(x) ll_resp(x, db, C), unname(p))   # N x n_free scores
  rse_i <- sqrt(diag(V_i %*% crossprod(S) %*% V_i))
  add_row("robust_se", lab, "rse_rel_vs_independent_sandwich", relmax(rse_k, rse_i),
          note = "independent Hessian and numDeriv::jacobian per-respondent scores")

  # (2)-(3) Apollo 0.3.5 from klue's optimum
  ap <- apollo_from_klue(f, db, C, lab)
  m <- ap$model; k2a <- ap$k2a; an <- unname(k2a)
  a_se  <- sqrt(diag(m$varcov))[an]; a_rse <- sqrt(diag(m$robvarcov))[an]
  ks  <- se_k[names(k2a)]; krs <- rse_k[names(k2a)]
  p_ap <- p; p_ap[names(k2a)] <- m$estimate[an]
  cv_ap <- .lcmnl_vcov(ctx, p_ap, C)
  if (anyNA(cv_ap$vcov)) stop("case ", lab, ": klue's covariance at Apollo's estimates is singular")
  fac <- sqrt(N / (N - 1))
  add_row("apollo_fit", lab, "apollo_successful_estimation", as.numeric(isTRUE(m$successfulEstimation)),
          note = "Apollo's own convergence flag (1 = TRUE)")
  add_row("apollo_fit", lab, "apollo_LL_minus_klue_LL", m$maximum - f$LL,
          note = sprintf("apollo_estimate from klue's optimum, %.1f s", ap$secs))
  add_row("apollo_fit", lab, "apollo_max_abs_parameter_move",
          max(abs(m$estimate[an] - p[names(k2a)])),
          note = paste("Apollo's Hessian:", m$hessianMethodUsed %||% "not recorded"))
  add_row("classical_se", lab, "se_rel_vs_apollo", relmax(ks, a_se), TOL,
          "SEs as each package reports them (klue at its optimum, Apollo at its own)")
  add_row("classical_se", lab, "se_rel_vs_apollo_same_point",
          relmax(sqrt(diag(cv_ap$vcov))[names(k2a)], a_se), TOL,
          "klue's covariance recomputed at Apollo's estimates")
  add_row("robust_se", lab, "rse_rel_vs_apollo", relmax(krs, a_rse),
          note = "before Apollo's factor")
  add_row("robust_se", lab, "sqrt_N_over_N_minus_1_minus_1", fac - 1,
          note = sprintf("N = %d respondents", N))
  add_row("robust_se", lab, "rse_rel_vs_apollo_after_factor", relmax(krs * fac, a_rse), TOL,
          "klue's robust SEs times sqrt(N/(N-1))")
  add_row("robust_se", lab, "rse_rel_vs_apollo_after_factor_same_point",
          relmax(sqrt(diag(cv_ap$robust_vcov))[names(k2a)] * fac, a_rse), TOL,
          "klue's sandwich recomputed at Apollo's estimates, times sqrt(N/(N-1))")

  tab <- data.frame(param = names(k2a), apollo = an, klue_est = unname(p[names(k2a)]),
                    apollo_est = unname(m$estimate[an]), se_klue = unname(ks),
                    se_indep = unname(se_i[match(names(k2a), names(p))]),
                    se_apollo = unname(a_se), rse_klue = unname(krs),
                    rse_indep = unname(rse_i[match(names(k2a), names(p))]),
                    rse_apollo = unname(a_rse))
  num <- vapply(tab, is.numeric, logical(1)); tab[num] <- signif(tab[num], 6)
  print(tab, row.names = FALSE)
  list(fit = f, db = db, C = C, se_indep = setNames(se_i, names(p)))
}

# Apollo 0.3.5's apollo_estimate (BGW, Apollo's default settings apart from
# writeIter = FALSE and silent = TRUE, as in the harness's .apollo_est)
# started at klue's optimum. .apollo_lc_setup builds the model and validates
# the inputs (start "generic" only fills apollo_beta's layout; every free
# entry is then replaced by klue's estimate, the fixed ones stay at 0).
apollo_from_klue <- function(fit, db, C, lab) {
  on.exit(cleanup_apollo(), add = TRUE)
  ok <- .apollo_lc_setup(db, C, DGP, paste0("lcmnl_se_check_", lab), start = "generic",
                         order = "parameter", ref = "last")
  if (!ok) stop("case ", lab, ": Apollo set-up failed")
  k2a <- klue_to_apollo(C)
  b <- apollo_beta
  stopifnot(all(unname(k2a) %in% names(b)),
            setequal(setdiff(names(b), unname(k2a)), apollo_fixed))
  b[unname(k2a)] <- unname(fit$par[names(k2a)])
  t0 <- proc.time()[["elapsed"]]
  m <- .apollo_est(b, paste("apollo_estimate, case", lab))
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(m) || is.null(m$varcov) || anyNA(m$varcov))
    stop("case ", lab, ": Apollo returned no covariance matrix")
  list(model = m, k2a = k2a, secs = secs)
}

# ---- (4): attribute scaling, case A ---------------------------------------------
scaling_check <- function(res, lab = "A") {
  f <- res$fit; db <- res$db; C <- res$C
  idx <- paste0("b", NB, "_class", 1:C)                 # the price coefficients
  se_jac <- function(d, p) {      # relative-step Hessian: numDeriv::jacobian of the analytic gradient
    ob <- .lcmnl_objective(.lcmnl_context(d, DGP), C)
    H <- numDeriv::jacobian(ob$grad_ll, unname(p)); H <- (H + t(H)) / 2
    tryCatch(setNames(sqrt(diag(solve(H))), names(p)),
             error = function(e) setNames(rep(NA_real_, length(p)), names(p)))
  }
  se_nd <- function(d, p) {       # Richardson Hessian of the negative log-likelihood
    ob <- .lcmnl_objective(.lcmnl_context(d, DGP), C)
    tryCatch(setNames(sqrt(diag(solve(numDeriv::hessian(ob$neg_ll, unname(p))))), names(p)),
             error = function(e) setNames(rep(NA_real_, length(p)), names(p)))
  }
  se1 <- sqrt(diag(f$vcov)); ref <- se_nd(db, f$par)
  add_row("scaling", lab, "klue_se_rel_vs_richardson_reference", relmax(se1, ref), TOL,
          "s = 1: klue's SEs against numDeriv::hessian of klue's log-likelihood")
  add_row("scaling", lab, "relstep_jacobian_se_rel_vs_richardson_reference",
          relmax(se_jac(db, f$par), ref), TOL, "s = 1")
  first_flag <- NA_real_
  for (s in SCALES) {
    d2 <- db
    for (j in 1:J) d2[[paste0("price_", j)]] <- d2[[paste0("price_", j)]] * s
    p2 <- f$par; p2[idx] <- p2[idx] / s
    cv2 <- suppressMessages(.lcmnl_vcov(.lcmnl_context(d2, DGP), p2, C))
    sing <- anyNA(cv2$vcov)
    back <- function(se) { se[idx] <- se[idx] * s; se }   # back to the s = 1 scale
    k_se <- if (sing) NULL else back(sqrt(diag(cv2$vcov)))
    ob2 <- .lcmnl_objective(.lcmnl_context(d2, DGP), C)
    Hk  <- stats::optimHess(unname(p2), ob2$neg_ll, ob2$grad_ll)
    ev  <- eigen((Hk + t(Hk)) / 2, symmetric = TRUE, only.values = TRUE)$values
    al  <- if (length(cv2$aliased)) paste(cv2$aliased, collapse = ",") else ""
    if (sing && is.na(first_flag)) first_flag <- s
    add_row("scaling", lab, "klue_singular_flag", as.numeric(sing), price_scale = s,
            note = if (sing) paste("covariance set to NA; aliased:", al) else "")
    add_row("scaling", lab, "klue_se_rel_err_vs_own_scale1",
            if (sing) NA_real_ else relmax(k_se, se1), price_scale = s,
            note = "klue's SEs mapped back to s = 1 against klue's SEs at s = 1")
    add_row("scaling", lab, "klue_se_rel_err_vs_richardson_reference",
            if (sing) NA_real_ else relmax(k_se, ref), TOL, price_scale = s,
            note = "against numDeriv::hessian at s = 1")
    add_row("scaling", lab, "relstep_jacobian_se_rel_err", relmax(back(se_jac(d2, p2)), ref), TOL,
            price_scale = s, note = "numDeriv::jacobian of the analytic gradient (Apollo 0.3.5's default)")
    add_row("scaling", lab, "richardson_hessian_se_rel_err", relmax(back(se_nd(d2, p2)), ref), TOL,
            price_scale = s, note = "numDeriv::hessian of the negative log-likelihood")
    add_row("scaling", lab, "klue_hessian_eigenvalue_ratio", min(ev) / max(ev), price_scale = s,
            note = sprintf("smallest / largest eigenvalue of optimHess; flag at <= %.2e",
                           sqrt(.Machine$double.eps)))
  }
  add_row("scaling", lab, "smallest_tested_scale_flagged_singular", first_flag,
          price_scale = NA_real_, note = paste("scales tested:", paste(fmt_scale(SCALES), collapse = ", ")))
}

# ---- run ---------------------------------------------------------------------
RES <- list()
for (lab in names(CASES)) RES[[lab]] <- run_case(lab)
stamp("part (4): price rescaling on case A")
scaling_check(RES$A)

tab <- do.call(rbind, ROWS)
tab$N <- vapply(tab$case, function(k) as.integer(CASES[[k]]$N), integer(1))
for (v in names(VERSIONS)) tab[[v]] <- VERSIONS[[v]]
dir.create("output", showWarnings = FALSE)
write.csv(tab, OUT_CSV, row.names = FALSE)

# ---- summary -------------------------------------------------------------------
pick <- function(part, measure) {
  r <- tab[tab$part == part & tab$measure == measure, , drop = FALSE]
  setNames(r$value, r$case)
}
by_case <- function(x) paste(sprintf("%s %.2e", names(x), x), collapse = ", ")
cat("\n================ SUMMARY (tolerance 1e-5) ================\n")
cat("(1) gradient vs numDeriv, independent LL:         ",
    by_case(pick("gradient", "grad_rel_vs_numDeriv_independent_LL")), "\n")
cat("(2) classical SE vs independent numerical Hessian:",
    by_case(pick("classical_se", "se_rel_vs_independent_numerical_hessian")), "\n")
cat("    classical SE vs Apollo 0.3.5:                 ",
    by_case(pick("classical_se", "se_rel_vs_apollo")), "\n")
cat("    same point (klue at Apollo's estimates):      ",
    by_case(pick("classical_se", "se_rel_vs_apollo_same_point")), "\n")
cat("(3) robust SE vs Apollo, before the factor:       ",
    by_case(pick("robust_se", "rse_rel_vs_apollo")), "\n")
cat("    sqrt(N/(N-1)) - 1:                            ",
    by_case(pick("robust_se", "sqrt_N_over_N_minus_1_minus_1")), "\n")
cat("    robust SE vs Apollo, after the factor:        ",
    by_case(pick("robust_se", "rse_rel_vs_apollo_after_factor")), "\n")
cat("    robust SE vs independent sandwich:            ",
    by_case(pick("robust_se", "rse_rel_vs_independent_sandwich")), "\n")
sc <- tab[tab$part == "scaling" & tab$measure == "klue_se_rel_err_vs_own_scale1", ]
cat("(4) klue SE error with price x s (case A):        ",
    paste(sprintf("x%s %s", fmt_scale(sc$price_scale),
                  ifelse(is.na(sc$value), "singular", sprintf("%.1e", sc$value))), collapse = ", "), "\n")
cl <- c(pick("classical_se", "se_rel_vs_independent_numerical_hessian"),
        pick("classical_se", "se_rel_vs_apollo"))
rb <- pick("robust_se", "rse_rel_vs_apollo_after_factor")
cat(sprintf("Classical SEs within 1e-5 of the independent Hessian and of Apollo in every case: %s (max %.2e)\n",
            all(cl <= TOL), max(cl)))
cat(sprintf("Robust SEs within 1e-5 of Apollo's after sqrt(N/(N-1)) in every case: %s (max %.2e)\n",
            all(rb <= TOL), max(rb)))
stamp("DONE: wrote %s (%d rows)", OUT_CSV, nrow(tab))
