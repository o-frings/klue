# MMNL benchmark via Apollo. One entry point, klue_mmnl(), covering the
# independent-normals specification (log-normal on price) and the correlated
# specification (full Cholesky covariance) via `correlation = TRUE`. Which
# parameters are random is set by `random` (default: every generic attribute
# and the price; constants fixed); constants follow dgp$asc_map and can be
# random too.
#
# The Apollo code follows Apollo's own example scripts
# (MMNL_preference_space.r and MMNL_preference_space_correlated.r): random
# coefficients mu + sigma * draws with sigma unconstrained, a negative
# log-normal price -exp(mu + sigma * draws), the correlated model written out
# as a lower-triangular Cholesky, apollo_estimate at Apollo's default
# estimation settings (writeIter = FALSE only skips the iterations file), the
# correlated model started from the independent estimates, and convergence as
# Apollo reports it (successfulEstimation). Where klue departs from the
# examples, and why:
#   - starting values: the examples hand-set them; klue starts from the pooled
#     MNL (a log-normal price at log(-b_MNL), or Apollo's recommended -3), with
#     standard deviations at the examples' 0.01;
#   - the correlated start sets the Cholesky off-diagonals to 0 (the current
#     example uses 0.01, an earlier revision 0), so that with the same draws
#     the start reproduces the independent log-likelihood exactly;
#   - 3000 MLHS draws (the examples use 500 Halton draws; Apollo's manual
#     advises against Halton beyond five random coefficients);
#   - BIC counts respondents, as for the LCMNL (Apollo's own model$BIC counts
#     observations and is returned as BIC_apollo);
#   - the reference alternative's constant is left out instead of fixed at 0
#     (the same model);
#   - on data with availability columns, utilities are measured from the
#     chosen alternative's (see .make_apollo_prob_mmnl), because Apollo 0.3.5
#     returns NaN otherwise; the probabilities are unchanged.
#
# Apollo's API requires apollo_* objects in the global environment. During an
# estimation klue assigns apollo_beta, apollo_fixed and apollo_probabilities
# there, as the examples' top-level script has them; cleanup_apollo() removes
# every .apollo_globals name afterwards, and klue_mmnl() then restores whatever
# the caller had under those names.

.apollo_globals <- c("apollo_control", "apollo_beta", "apollo_fixed",
                     "apollo_draws", "apollo_randCoeff",
                     "apollo_probabilities", "apollo_inputs", "apollo_lcPars")

#' Remove Apollo's global objects
#'
#' Deletes the \code{apollo_*} objects that Apollo requires in the global environment, so a session can mix klue and Apollo runs.
#' @export
cleanup_apollo <- function() {
  drop <- intersect(.apollo_globals, ls(envir = .GlobalEnv, all.names = TRUE))
  if (length(drop)) rm(list = drop, envir = .GlobalEnv)
}

# Parameter specification shared by the builders. The namespace is
# x1..xN (generic attributes), price, and asc1..ascK (one per asc_map group).
# `random` picks the random ones; everything else is a fixed coefficient.
# Random parameters keep the canonical order x's, price, asc's, which fixes the
# draw order and the Cholesky ordering. Parameter names in apollo_beta:
#   fixed attribute / price: b_x{a}, b_price
#   fixed constant:          asc_alt{g} under the default per-alternative
#                            asc_map (the historical name), else asc{g}
#   random p, independent:   mu_{p}, sigma_{p} (the standard deviation, sign
#                            arbitrary; for a log-normal price, of its log)
#   random p, correlated:    mu_{p}, Cholesky rows s_{p}_{1..i} (s_pr for price)
.mmnl_spec <- function(dgp = DGP_DEFAULT, random = NULL,
                       price = c("lognormal", "normal", "fixed")) {
  price <- match.arg(price)
  ng <- dgp$n_generic; K <- dgp$n_asc
  xs    <- paste0("x", seq_len(ng))
  ascs  <- if (K > 0) paste0("asc", seq_len(K)) else character(0)
  space <- c(xs, "price", ascs)
  given <- !is.null(random)
  if (!given) random <- c(xs, "price")
  random <- unique(as.character(random))
  bad <- setdiff(random, space)
  if (length(bad))
    stop("klue_mmnl: unknown random parameter(s) ", paste(bad, collapse = ", "),
         "; valid names are ", paste(space, collapse = ", "), call. = FALSE)
  if (identical(price, "fixed")) {
    if (given && "price" %in% random)
      stop("klue_mmnl: price = \"fixed\" conflicts with \"price\" in `random`",
           call. = FALSE)
    random <- setdiff(random, "price")
  }
  random <- space[space %in% random]
  if (!length(random))
    stop("klue_mmnl: no random parameters; use estimate_lcmnl(C = 1) for an MNL",
         call. = FALSE)
  default_map <- identical(as.integer(dgp$asc_map),
                           c(seq_len(dgp$n_alternatives - 1L), 0L))
  list(random    = random,
       price     = price,
       xs        = xs,
       ascs      = ascs,
       asc_fixed = if (default_map) paste0("asc_alt", seq_len(K)) else ascs,
       lognormal = "price" %in% random && identical(price, "lognormal"),
       chol      = ifelse(random == "price", "s_pr", paste0("s_", random)))
}

# apollo_randCoeff builder, as in Apollo's MMNL examples. Independent:
# b_p = mu_p + sigma_p * draw. Correlated: b = mu + L * draws with
# lower-triangular Cholesky L, every element unconstrained. A log-normal price
# wraps its normal core as -exp(.). The terms are written out as generated
# code, one line per random coefficient, as the examples write them.
.make_apollo_randCoeff <- function(dgp = DGP_DEFAULT, correlation = FALSE,
                                   spec = .mmnl_spec(dgp)) {
  r <- spec$random
  wrap <- function(p, core)
    if (p == "price" && spec$lognormal) sprintf("-exp(%s)", core) else core
  lines <- c("function(apollo_beta, apollo_inputs) {", "  randcoeff <- list()")
  for (i in seq_along(r)) {
    p <- r[i]
    core <- if (!correlation) {
      sprintf("mu_%s + sigma_%s * draws_%s", p, p, p)
    } else {
      sprintf("mu_%s + %s", p,
              paste(paste0(spec$chol[i], "_", seq_len(i), " * draws_", r[seq_len(i)]),
                    collapse = " + "))
    }
    lines <- c(lines, sprintf('  randcoeff[["b_%s"]] <- %s', p, wrap(p, core)))
  }
  lines <- c(lines, '  return(randcoeff)', '}')
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- globalenv()
  fn
}

# apollo_probabilities builder (same for both MMNL flavours). Each alternative
# gets its asc_map group's constant (random -> b_asc{g}, fixed -> the fixed
# name), so alternatives sharing a group share one constant.
.make_apollo_prob_mmnl <- function(dgp = DGP_DEFAULT, has_avail = FALSE,
                                   spec = .mmnl_spec(dgp)) {
  J <- dgp$n_alternatives; n_generic <- dgp$n_generic
  alt_entries   <- paste(sprintf('alt%d = %d', 1:J, 1:J), collapse = ", ")
  # Mask unavailable alternatives from the database's av_1..av_J columns when
  # present; otherwise every alternative is available (the default).
  avail_entries <- if (isTRUE(has_avail))
                     paste(sprintf('alt%d = av_%d', 1:J, 1:J), collapse = ", ")
                   else
                     paste(sprintf('alt%d = 1', 1:J), collapse = ", ")
  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    '  V <- list()'
  )
  for (j in 1:J) {
    g <- dgp$asc_map[j]
    asc_term <- if (g > 0L) {
      if (paste0("asc", g) %in% spec$random) paste0("b_asc", g) else spec$asc_fixed[g]
    }
    terms <- c(asc_term,
               sprintf('b_x%d * x%d_%d', 1:n_generic, 1:n_generic, j),
               sprintf('b_price * price_%d', j))
    lines <- c(lines, sprintf('  V[["alt%d"]] <- %s', j, paste(terms, collapse = " + ")))
  }
  if (isTRUE(has_avail)) {
    # Measure utilities from the chosen alternative's. The probabilities are
    # unchanged, but Apollo 0.3.5 (and 0.3.9) apollo_mnl sets unavailable
    # utilities to 0 and then multiplies exp(V_j - V_chosen) by availability,
    # which is Inf * 0 = NaN once a chosen utility falls below about -709
    # (extreme random-coefficient draws); no Apollo setting avoids it.
    # Relative to the chosen utility that term is exp(0) * 0 = 0.
    lines <- c(lines,
      sprintf('  Vc <- %s', paste(sprintf('(CHOICE == %d) * V[["alt%d"]]', 1:J, 1:J),
                                 collapse = " + ")),
      sprintf('  V[["alt%d"]] <- V[["alt%d"]] - Vc', 1:J, 1:J))
  }
  lines <- c(lines,
    sprintf('  mnl_settings <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE, utilities = V)',
            alt_entries, avail_entries),
    '  P[["model"]] <- apollo_mnl(mnl_settings, functionality)',
    '  P <- apollo_panelProd(P, apollo_inputs, functionality)',
    '  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)',
    '  return(P)',
    '}'
  )
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- globalenv()
  fn
}

# One Apollo MMNL estimation with given draws and starting values: a single
# apollo_estimate call at Apollo's default settings (estimationRoutine, BGW by
# default, is passed through; writeIter = FALSE only skips the iterations
# file). All components go to apollo_validateInputs as explicit arguments
# (officially supported); apollo_beta (the start), apollo_fixed and
# apollo_probabilities also sit in the global environment during estimation,
# as in the examples' top-level script (Apollo's apollo_expandLoop validates
# its rewrite at the global apollo_beta, and apollo_dVdB fetches the global
# apollo_probabilities for analytic gradients). cleanup_apollo() removes them
# afterwards and also clears stray apollo_* objects at entry (they would
# otherwise leak in through validateInputs' globalenv fallback).
# Returns Apollo's model, or list(klue_error = <message>) when validation or
# estimation stops with an error.
.run_apollo_mmnl <- function(database, n_draws, start_beta,
                             dgp                = DGP_DEFAULT,
                             n_cores            = NULL,
                             draws_type         = DRAWS_TYPE_MMNL,
                             estimation_routine = ESTIMATION_ROUTINE_MMNL,
                             correlation        = FALSE,
                             spec               = .mmnl_spec(dgp),
                             memory_saver       = FALSE) {
  # apollo_estimate expects the apollo package on the search path; attach it if
  # a caller only loaded the namespace (apollo is Suggests since 0.9.2).
  if (!"package:apollo" %in% search())
    try(suppressMessages(attachNamespace("apollo")), silent = TRUE)
  cleanup_apollo()
  if (is.null(n_cores)) n_cores <- getOption("klue.mmnl.n_cores", 1L)
  n_cores <- max(1L, as.integer(n_cores))

  control <- list(
    # unique per call: with nCores > 1 Apollo keeps a temporary inputs file
    # named after the model in tempdir(), shared by forked processes
    modelName       = paste0(if (correlation) "klue_mmnl_corr_" else "klue_mmnl_",
                             Sys.getpid(), "_", sub("^file", "", basename(tempfile()))),
    modelDescr      = if (correlation) "MMNL, correlated random coefficients"
                      else "MMNL, independent random coefficients",
    indivID         = "ID",
    nCores          = n_cores,
    memorySaver     = isTRUE(memory_saver),
    mixing          = TRUE,
    outputDirectory = tempdir()
  )
  draws <- list(
    interDrawsType = draws_type,
    # a double, as in the examples: with an integer, Apollo 0.3.5's core-count
    # hint (silent = FALSE, nCores = 1) computes nrow * draws^2 in integer
    # arithmetic and stops once that exceeds 2^31 - 1
    interNDraws    = as.numeric(n_draws),
    interUnifDraws = c(),
    interNormDraws = paste0("draws_", spec$random),
    intraDrawsType = draws_type,
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )
  has_avail <- all(paste0("av_", 1:dgp$n_alternatives) %in% names(database))
  probabilities <- .make_apollo_prob_mmnl(dgp, has_avail = has_avail, spec = spec)
  # register cleanup BEFORE the global assignment so an error anywhere below
  # cannot leak apollo_probabilities into the caller's global environment
  on.exit(cleanup_apollo(), add = TRUE)
  assign("apollo_beta",          start_beta,    envir = .GlobalEnv)
  assign("apollo_fixed",         c(),           envir = .GlobalEnv)
  assign("apollo_probabilities", probabilities, envir = .GlobalEnv)

  inputs <- tryCatch(
    apollo::apollo_validateInputs(
      apollo_beta      = start_beta,
      apollo_fixed     = c(),
      database         = database,
      apollo_control   = control,
      apollo_draws     = draws,
      apollo_randCoeff = .make_apollo_randCoeff(dgp, correlation, spec = spec)
    ),
    error = function(e) list(klue_error = paste0("apollo_validateInputs_error: ",
                                                   conditionMessage(e)))
  )
  if (!is.null(inputs$klue_error)) return(inputs)

  tryCatch(
    apollo::apollo_estimate(
      apollo_beta          = start_beta,
      apollo_fixed         = c(),
      apollo_probabilities = probabilities,
      apollo_inputs        = inputs,
      estimate_settings    = list(estimationRoutine = estimation_routine,
                                  writeIter         = FALSE)
    ),
    error = function(e) list(klue_error = paste0("apollo_estimate_error: ",
                                                 conditionMessage(e)))
  )
}

# Map an independent fit's estimates onto the correlated parameterisation at the
# point where the two models coincide, the apollo_readBeta step of Apollo's
# correlated example: constants, fixed coefficients and means carry over, the
# Cholesky diagonal takes the independent standard deviations sigma_p (sign
# kept), off-diagonals are 0. With the same draws, the correlated likelihood at
# this point equals the independent one, so a correlated fit started here
# cannot end below it.
.mmnl_indep_to_corr <- function(ind_par, beta0, spec) {
  r <- spec$random
  keep <- intersect(names(beta0), names(ind_par))
  beta0[keep] <- ind_par[keep]
  for (i in seq_along(r)) {
    for (l in seq_len(i)) {
      beta0[paste0(spec$chol[i], "_", l)] <-
        if (l == i) unname(ind_par[paste0("sigma_", r[i])]) else 0
    }
  }
  beta0
}

#' MMNL benchmark (independent or correlated normals)
#'
#' Estimates a mixed logit with Apollo the way Apollo's own example scripts
#' do: random coefficients \eqn{\mu + \sigma \xi} with \eqn{\sigma}
#' unconstrained, a negative log-normal price, `apollo_estimate` at Apollo's
#' default estimation settings, and convergence as Apollo reports it. The
#' departures from the examples (pooled-MNL starts, zero starting
#' correlations, 3000 MLHS draws, BIC on respondents, and utilities measured
#' from the chosen alternative on data with availability columns) are listed
#' in the source file's header.
#' By default every generic attribute is a random normal, the price a negative
#' log-normal, and the constants fixed; `random` and `price` change that.
#' Correlated flavour (`correlation = TRUE`): a full lower-triangular Cholesky
#' covariance over the random parameters, started from the independent
#' estimates (Apollo's `apollo_readBeta` step: the same means, the independent
#' standard deviations on the diagonal, zero correlations).
#'
#' @param database Data frame in klue's canonical wide format (one row per
#'   task) with an `ID` column and the attribute and choice columns expected by
#'   the data-generating process.
#' @param correlation Logical; if `TRUE`, estimate the correlated specification
#'   with a full lower-triangular Cholesky covariance instead of independent
#'   normals.
#' @param n_draws Number of inter-individual draws. `NULL` uses the package
#'   default.
#' @param draws_type Apollo inter-draws type (for example Halton or MLHS).
#'   `NULL` uses the package default.
#' @param estimation_routine Apollo estimation routine passed to
#'   `apollo_estimate`. `NULL` uses the package default (Apollo's default,
#'   BGW).
#' @param n_cores Number of cores for Apollo. `NULL` uses the package default
#'   (Apollo's default, 1).
#' @param quiet Logical; if `TRUE`, redirect Apollo output to a temporary log
#'   file whose tail is attached to a failing result. `NULL` uses the package
#'   default.
#' @param start Starting-value scheme. `"informed"` (default) starts the
#'   independent model from the pooled MNL estimates (means, fixed
#'   coefficients and constants; a log-normal price at
#'   \eqn{\mu = \log(-\beta_{MNL})}, or at Apollo's recommended \eqn{-3} when
#'   the pooled price coefficient is not negative) and the correlated model
#'   from an independent fit with the same draws. `"zero"` and `"sign"` are the
#'   conventional neutral starts: attribute means at `0` and at `+0.01`
#'   respectively, with a log-normal price coefficient started at `-0.01` in
#'   both (a log-normal cannot start at exactly zero, and `mu_price = 0` would
#'   mean `b_price = -1`, not a neutral start). Standard deviations start at
#'   `0.01` under every scheme, as in Apollo's examples. Use a neutral start
#'   to show that a fit does not depend on being warm-started. Ignored (with
#'   a warning) when `warm_start` is given.
#' @param random Character vector naming the random parameters among
#'   `x1..xN`, `price` and `asc1..ascK`, where `ascK` indexes the ASC groups of
#'   `dgp$asc_map` (with the default map, one per non-reference alternative).
#'   `NULL` (default) makes every generic attribute and the price random and
#'   the constants fixed. Parameters not named are fixed coefficients. Random
#'   constants are normal, with the reference alternative's constant fixed at
#'   zero, so an independent random-constant model depends on which
#'   alternative is the reference.
#' @param price Distribution of a random price coefficient: `"lognormal"`
#'   (default; negative log-normal), `"normal"`, or `"fixed"` (the price is not
#'   random, whatever `random` says by default).
#' @param warm_start Optional earlier `klue_mmnl` fit (klue >= 0.10.0) of the
#'   same data, `random` set and `price` distribution to start from. An
#'   independent fit given to a correlated call is mapped to the point where
#'   the two models coincide (same means, its standard deviations on the
#'   Cholesky diagonal, zero correlations), so with the same draws the
#'   correlated fit cannot end below it; the result's `nested_ok` reports
#'   whether it did not. A fit of the same flavour is copied by name.
#' @param dgp Data-generating-process specification giving the number of
#'   alternatives, generic attributes, parameters and the ASC groups.
#' @param memory_saver Logical; Apollo's `apollo_control$memorySaver`. `TRUE`
#'   computes the analytic gradient in chunks of about two respondents, which
#'   cuts its memory (about one observations-by-draws array per parameter)
#'   without changing the estimates. `NULL` uses the package default (Apollo's
#'   default, `FALSE`; set `options(klue.mmnl.memory_saver = TRUE)` or the
#'   environment variable `KLUE_MMNL_MEMORY_SAVER=TRUE` to change it).
#' @param n_draws_stage1,mu_price_bounds,sigma_price_bounds Deprecated and
#'   ignored (with a warning). Since klue 0.10.0 there is no coarse-draw first
#'   stage, and Apollo 0.3.5 has no `bounds` setting, so the price box of
#'   earlier versions was never applied.
#' @return A list describing the fit. `converged` is Apollo's
#'   `successfulEstimation`. On success: the log-likelihood `LL` (Apollo's
#'   final log-likelihood), information criteria `BIC` and `AIC` (with `N` =
#'   respondents, as for the LCMNL; Apollo's own BIC, on observations, is
#'   `BIC_apollo`), the number of free parameters `k`, `reason = "ok"`, the
#'   resolved `settings`,
#'   Apollo's status (`apollo_status`: `successfulEstimation`, `code`,
#'   `message`, `nIter`), the starting values `start_values` and their source
#'   `start_used`, the full estimate vector `par` with its classical and
#'   respondent-clustered covariances `vcov` and `robust_vcov` (Apollo's
#'   `varcov` and `robvarcov`; all `NA` when Apollo could not invert the
#'   Hessian), the resolved `random` set and `price` distribution,
#'   `sd_scale = "linear"`, and, for the independent specification, the means
#'   `mu` and standard deviations `sigma = |sigma_p|` of the random parameters
#'   (for a log-normal price, of the underlying normal). A correlated fit that
#'   fitted its own independent start also returns it as `independent_fit`.
#'   On failure: `converged = FALSE` with `LL`, `BIC`, `AIC`, and `k` set to
#'   non-informative values, a `reason` string, Apollo's status and the
#'   estimates where Apollo returned them (`LL_at_stop`, `par`), the
#'   `independent_fit` where there was one, and the Apollo log tail and path
#'   when `quiet` is `TRUE`.
#'
#'   Side effects: during estimation the Apollo objects live in the global
#'   environment, as Apollo requires; any `apollo_*` objects the caller had
#'   there are restored afterwards, and so is the random-number stream, which
#'   Apollo reseeds when it makes its draws.
#' @export
klue_mmnl <- function(database,
                      correlation        = FALSE,
                      n_draws            = NULL,
                      draws_type         = NULL,
                      estimation_routine = NULL,
                      n_cores            = NULL,
                      quiet              = NULL,
                      start              = c("informed", "zero", "sign"),
                      random             = NULL,
                      price              = c("lognormal", "normal", "fixed"),
                      warm_start         = NULL,
                      dgp                = DGP_DEFAULT,
                      memory_saver       = NULL,
                      n_draws_stage1     = NULL,
                      mu_price_bounds    = NULL,
                      sigma_price_bounds = NULL) {
  start <- match.arg(start)
  price <- match.arg(price)
  if (!requireNamespace("apollo", quietly = TRUE)) {
    stop("klue_mmnl requires the 'apollo' package (Suggests since 0.9.2); ",
         "install it with install.packages(\"apollo\")", call. = FALSE)
  }
  for (a in c("n_draws_stage1", "mu_price_bounds", "sigma_price_bounds"))
    if (!is.null(get(a)))
      warning("klue_mmnl: `", a, "` is deprecated and ignored (klue 0.10.0 ",
              "estimates once, as Apollo's examples do, and Apollo has no ",
              "bounds setting)", call. = FALSE)
  spec <- .mmnl_spec(dgp, random, price)
  r    <- spec$random
  # Restore on exit whatever the caller had under Apollo's global names and
  # the caller's random-number stream (Apollo reseeds it when making draws).
  # Registered first and run last (the handlers below are prepended).
  snap <- mget(intersect(.apollo_globals, ls(envir = .GlobalEnv, all.names = TRUE)),
               envir = .GlobalEnv)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    cleanup_apollo()
    for (nm in names(snap)) assign(nm, snap[[nm]], envir = .GlobalEnv)
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  if (!is.null(warm_start)) {
    if (!identical(warm_start$sd_scale, "linear"))
      stop("klue_mmnl: `warm_start` comes from klue < 0.10.0, whose sigma_* ",
           "were log standard deviations; refit it with this version",
           call. = FALSE)
    if (!isTRUE(warm_start$converged) || is.null(warm_start$par))
      stop("klue_mmnl: `warm_start` must be a converged klue_mmnl fit with $par",
           call. = FALSE)
    if (!identical(warm_start$random, r))
      stop("klue_mmnl: `warm_start` has random set (",
           paste(warm_start$random, collapse = ", "), "), this call (",
           paste(r, collapse = ", "), ")", call. = FALSE)
    if (!identical(warm_start$price, price))
      stop("klue_mmnl: `warm_start` has price = \"", warm_start$price,
           "\", this call price = \"", price, "\"", call. = FALSE)
    if (!identical(start, "informed"))
      warning("klue_mmnl: `start = \"", start, "\"` is ignored because a ",
              "`warm_start` is given", call. = FALSE)
  }
  # apollo_estimate expects the apollo package on the search path (it looks up
  # "package:apollo"); attach it if a caller only loaded the namespace.
  if (!"package:apollo" %in% search())
    try(suppressMessages(attachNamespace("apollo")), silent = TRUE)
  d <- klue_mmnl_defaults()
  if (is.null(n_draws))            n_draws            <- d$n_draws
  if (is.null(draws_type))         draws_type         <- d$draws_type
  if (is.null(estimation_routine)) estimation_routine <- d$estimation_routine
  if (is.null(n_cores))            n_cores            <- d$n_cores
  if (is.null(quiet))              quiet              <- d$quiet
  if (is.null(memory_saver))       memory_saver       <- d$memory_saver

  # The correlated "informed" start is an independent fit with the same draws;
  # fit it before this call's own log redirection so each fit keeps its log.
  indep_fit <- NULL
  if (correlation && is.null(warm_start) && identical(start, "informed")) {
    indep_fit <- klue_mmnl(database, correlation = FALSE, n_draws = n_draws,
                           draws_type = draws_type,
                           estimation_routine = estimation_routine,
                           n_cores = n_cores, quiet = quiet,
                           random = r, price = price, dgp = dgp,
                           memory_saver = memory_saver)
  }

  cleanup_apollo()

  # If `quiet`, redirect Apollo's output to a tempfile whose last 40 lines are
  # attached to a failing result (the file is kept for inspection).
  log_file <- tempfile(pattern = if (correlation) "klue_mmnl_corr_" else "klue_mmnl_",
                       fileext = ".log")
  if (isTRUE(quiet)) {
    old_sink_out <- sink.number()
    old_sink_msg <- sink.number(type = "message")
    msg_con <- file(log_file, open = "at")
    sink(log_file)
    sink(msg_con, type = "message")
    on.exit({
      while (sink.number() > old_sink_out) sink()
      while (sink.number(type = "message") > old_sink_msg) sink(type = "message")
      try(close(msg_con), silent = TRUE)
    }, add = TRUE, after = FALSE)
  }

  N <- length(unique(database$ID))
  n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  settings <- list(n_draws = n_draws, draws_type = draws_type,
                   estimation_routine = estimation_routine,
                   n_cores = n_cores, memory_saver = isTRUE(memory_saver),
                   N = N, n_rows = nrow(database))

  read_log_tail <- function() {
    if (!isTRUE(quiet) || !file.exists(log_file)) return(NULL)
    out <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) NULL)
    if (length(out) == 0) NULL else utils::tail(out, 40L)
  }
  fail_with <- function(reason, model = NULL, beta0 = NULL, start_used = NULL,
                        log_from = NULL) {
    res <- list(converged = FALSE, LL = -Inf, BIC = Inf, AIC = Inf, k = 0,
                reason          = reason,
                apollo_log_tail = if (!is.null(log_from)) log_from$apollo_log_tail
                                  else read_log_tail(),
                apollo_log_path = if (!is.null(log_from)) log_from$apollo_log_path
                                  else if (isTRUE(quiet) && file.exists(log_file)) log_file,
                settings        = settings,
                random          = r,
                price           = price,
                correlation     = isTRUE(correlation),
                sd_scale        = "linear",
                start_values    = beta0,
                start_used      = start_used)
    if (!is.null(model)) {
      res$apollo_status <- .apollo_status(model)
      res$LL_at_stop    <- if (is.numeric(model$maximum)) model$maximum else NA_real_
      res$par           <- model$estimate
    }
    if (!is.null(indep_fit)) res$independent_fit <- indep_fit
    if (!correlation) {
      res$mu <- rep(0, length(r)); res$sigma <- rep(0, length(r))
    }
    res
  }

  # ---- Starting values ------------------------------------------------------
  # Layout: fixed constants, then one location per attribute and price (mu_p if
  # random, b_p if fixed), then the means of random constants, then the spreads
  # (independent: sigma_p per random p; correlated: Cholesky rows). With the
  # default specification this is the historical asc_alt*, mu_*, sigma_* order.
  # "informed" is the apollo_readBeta-from-a-simpler-model start: every
  # location and constant takes its pooled-MNL estimate (klue's MNL matches
  # Apollo's to about 1e-5). Neutral starts: the log-normal price enters as
  # b_price = -exp(mu_price + sigma * draw), so a neutral price start is a small
  # NEGATIVE COEFFICIENT, mu_price = log(0.01), not mu_price = 0 (which would be
  # b_price = -1). Standard deviations start at 0.01 throughout, the small
  # positive spread of Apollo's current examples.
  NEUTRAL_MU    <- if (identical(start, "sign")) 0.01 else 0
  NEUTRAL_PRICE <- if (spec$lognormal) log(0.01) else -NEUTRAL_MU
  NEUTRAL_SD    <- 0.01
  loc_name <- function(p) if (p %in% r) paste0("mu_", p) else paste0("b_", p)
  beta0 <- c()
  for (g in seq_len(dgp$n_asc))
    if (!(spec$ascs[g] %in% r)) beta0[spec$asc_fixed[g]] <- 0
  # MNL-informed locations (MNL = LCMNL with C=1), unless a neutral start or a
  # correlated fit (which starts from an independent fit) is asked for.
  want_mnl <- !correlation && identical(start, "informed") && is.null(warm_start)
  mnl_fit  <- if (want_mnl)
    tryCatch(estimate_lcmnl(database, C = 1,
                            start_betas = matrix(0, nrow = 1, ncol = n_beta),
                            dgp = dgp, vcov = FALSE),
             error = function(e) NULL)
  # A pooled MNL that stopped at its iteration cap still gives a usable start.
  mnl_ok <- !want_mnl || (!is.null(mnl_fit$par) && all(is.finite(mnl_fit$par)))
  if (!mnl_ok)
    return(fail_with("pooled_mnl_start_failed", start_used = "pooled_mnl"))
  if (want_mnl && !isTRUE(mnl_fit$converged))
    message("klue_mmnl: the pooled MNL for the starting values stopped at its ",
            "iteration cap; using its estimates as the start")
  mnl_b <- if (want_mnl) mnl_fit$betas[1, ]
  mnl_a <- if (want_mnl && dgp$n_asc > 0) mnl_fit$par[paste0("asc", seq_len(dgp$n_asc))]
  for (g in seq_len(dgp$n_asc))
    if (!(spec$ascs[g] %in% r)) beta0[spec$asc_fixed[g]] <- if (want_mnl) unname(mnl_a[g]) else 0
  for (a in seq_len(n_generic))
    beta0[loc_name(spec$xs[a])] <- if (want_mnl) mnl_b[a] else NEUTRAL_MU
  if (want_mnl && spec$lognormal) {
    # Log-normal price: b_price = -exp(mu_price + sigma_price * draw) is always
    # negative, so a negative pooled-MNL price coefficient b gives
    # mu_price = log(-b), i.e. the starting median b_price = b. A price
    # coefficient that is not negative has no log-normal counterpart; Apollo's
    # FAQ then recommends starting the log-mean at -3.
    bp <- mnl_b[n_beta]
    if (is.finite(bp) && bp < 0) {
      beta0["mu_price"] <- log(-bp)
    } else {
      message("klue_mmnl: pooled-MNL price coefficient is ", signif(bp, 3),
              " (not negative); starting the log-normal price at mu_price = -3")
      beta0["mu_price"] <- -3
    }
  } else {
    beta0[loc_name("price")] <- if (want_mnl) mnl_b[n_beta] else NEUTRAL_PRICE
  }
  for (p in intersect(spec$ascs, r))
    beta0[paste0("mu_", p)] <- if (want_mnl) unname(mnl_a[match(p, spec$ascs)]) else 0
  if (!correlation) {
    for (p in r) beta0[paste0("sigma_", p)] <- NEUTRAL_SD
  } else {
    for (i in seq_along(r)) for (l in seq_len(i))
      beta0[paste0(spec$chol[i], "_", l)] <- if (l == i) NEUTRAL_SD else 0
  }
  start_used <- if (!want_mnl) start
                else if (isTRUE(mnl_fit$converged)) "pooled_mnl" else "pooled_mnl_maxit"

  # A fit to nest: the warm start, or the correlated model's own independent fit.
  nest_fit <- NULL
  if (!is.null(warm_start)) {
    # Same-flavour fits carry over by name; an independent fit maps onto the
    # correlated parameterisation exactly.
    wpar <- warm_start$par
    if (correlation && !isTRUE(warm_start$correlation)) {
      chol_nm <- unlist(lapply(seq_along(r), function(i) paste0(spec$chol[i], "_", seq_len(i))))
      miss <- setdiff(setdiff(names(beta0), chol_nm), names(wpar))
      if (length(miss))
        stop("klue_mmnl: `warm_start` lacks parameter(s) ", paste(miss, collapse = ", "),
             " of this specification (different data-generating process?)", call. = FALSE)
      beta0 <- .mmnl_indep_to_corr(wpar, beta0, spec)
    } else {
      if (!setequal(names(wpar), names(beta0)))
        stop("klue_mmnl: `warm_start` parameters do not match this specification",
             call. = FALSE)
      beta0[names(beta0)] <- wpar[names(beta0)]
    }
    nest_fit   <- warm_start
    start_used <- "warm_start"
  } else if (!is.null(indep_fit)) {
    if (!isTRUE(indep_fit$converged))
      return(fail_with(paste0("independent_start_fit_failed: ", indep_fit$reason),
                       start_used = "independent_fit", log_from = indep_fit))
    beta0      <- .mmnl_indep_to_corr(indep_fit$par, beta0, spec)
    nest_fit   <- indep_fit
    start_used <- "independent_fit"
  }

  # ---- One estimation at n_draws, as in Apollo's examples ---------------------
  model <- .run_apollo_mmnl(database, n_draws, beta0, dgp = dgp, n_cores = n_cores,
                            draws_type = draws_type,
                            estimation_routine = estimation_routine,
                            correlation = correlation, spec = spec,
                            memory_saver = memory_saver)
  if (!is.null(model$klue_error))
    return(fail_with(model$klue_error, beta0 = beta0, start_used = start_used))
  if (!isTRUE(model$successfulEstimation))
    return(fail_with(paste0("apollo_unsuccessful: ",
                            if (is.null(model$message)) "no message" else model$message),
                     model = model, beta0 = beta0, start_used = start_used))
  LL <- model$maximum                      # Apollo's LL(final)
  if (!is.numeric(LL) || length(LL) != 1L || !is.finite(LL))
    return(fail_with("apollo_non_finite_LL", model = model, beta0 = beta0,
                     start_used = start_used))

  est    <- model$estimate
  n_free <- length(est)
  res <- list(converged = TRUE, LL = LL,
              BIC = -2 * LL + n_free * log(N),
              AIC = -2 * LL + 2 * n_free, k = n_free,
              BIC_apollo = if (is.numeric(model$BIC)) model$BIC else NA_real_,
              reason = "ok", apollo_log_tail = NULL, apollo_log_path = NULL,
              settings = settings,
              apollo_status = .apollo_status(model),
              start_values = beta0,
              start_used   = start_used)
  if (!correlation) {
    res$mu    <- est[paste0("mu_", r)]
    res$sigma <- abs(est[paste0("sigma_", r)])
  }
  res$par         <- est
  res$vcov        <- model$varcov
  res$robust_vcov <- model$robvarcov
  res$random      <- r
  res$price       <- price
  res$correlation <- isTRUE(correlation)
  res$sd_scale    <- "linear"
  if (!is.null(nest_fit)) {
    # A start from a nested (or the same) model with the same draws cannot
    # legitimately end below it; if it does, the optimiser stopped early.
    same_draws <- isTRUE(all.equal(as.numeric(nest_fit$settings$n_draws),
                                   as.numeric(n_draws))) &&
                  identical(nest_fit$settings$draws_type, draws_type) &&
                  identical(nest_fit$settings$N, N) &&
                  identical(nest_fit$settings$n_rows, nrow(database))
    res$nested_ok <- if (same_draws) LL >= nest_fit$LL - 1e-6 * abs(nest_fit$LL) else NA
    if (isFALSE(res$nested_ok))
      warning(sprintf("klue_mmnl: LL %.3f is below the %.3f of the fit it starts from, which it nests",
                      LL, nest_fit$LL), call. = FALSE)
  }
  if (!is.null(indep_fit)) res$independent_fit <- indep_fit
  res
}

# Apollo's own convergence record for a fitted model.
.apollo_status <- function(model)
  list(successfulEstimation = isTRUE(model$successfulEstimation),
       code    = if (is.null(model$code)) NA_integer_ else model$code,
       message = if (is.null(model$message)) NA_character_ else model$message,
       nIter   = if (is.null(model$nIter)) NA_integer_ else model$nIter)

#' @rdname klue_mmnl
#' @param ... For \code{klue_mmnl_corr}, arguments passed to \code{klue_mmnl}.
#' @export
klue_mmnl_corr <- function(database, ...) klue_mmnl(database, correlation = TRUE, ...)
