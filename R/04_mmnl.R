# MMNL benchmark via Apollo. One entry point, klue_mmnl(), covering the
# independent-normals specification (log-normal on price) and the correlated
# specification (full Cholesky covariance) via `correlation = TRUE`.
#
# Apollo's API requires apollo_* objects in the global environment; we assign
# exactly the names in .apollo_globals and remove exactly those afterwards.

.apollo_globals <- c("apollo_control", "apollo_beta", "apollo_fixed",
                     "apollo_draws", "apollo_randCoeff",
                     "apollo_probabilities", "apollo_inputs", "apollo_lcPars")

cleanup_apollo <- function() {
  drop <- intersect(.apollo_globals, ls(envir = .GlobalEnv, all.names = TRUE))
  if (length(drop)) rm(list = drop, envir = .GlobalEnv)
}

# apollo_randCoeff builder. Independent: b_xa = mu + exp(sigma)*draw, price
# negative log-normal. Correlated: b = mu + L*draws with lower-triangular
# Cholesky L, price -exp(mu + L_row * draws). Built as generated code because
# Apollo's checkIndices rejects for-loops inside model functions.
.make_apollo_randCoeff <- function(dgp = DGP_DEFAULT, correlation = FALSE) {
  n_generic <- dgp$n_generic; n_beta <- dgp$n_beta
  lines <- c("function(apollo_beta, apollo_inputs) {", "  randcoeff <- list()")
  if (!correlation) {
    for (a in 1:n_generic) {
      lines <- c(lines, sprintf(
        '  randcoeff[["b_x%d"]] <- mu_x%d + exp(sigma_x%d) * draws_x%d', a, a, a, a))
    }
    lines <- c(lines, '  randcoeff[["b_price"]] <- -exp(mu_price + exp(sigma_price) * draws_price)')
  } else {
    attr_short  <- c(paste0("x", 1:n_generic), "price")
    draw_names  <- paste0("draws_", attr_short)
    chol_prefix <- c(paste0("s_x", 1:n_generic), "s_pr")
    for (a in 1:n_beta) {
      inner <- paste(paste0(chol_prefix[a], "_", 1:a, " * ", draw_names[1:a]),
                     collapse = " + ")
      if (a < n_beta) {
        lines <- c(lines, sprintf('  randcoeff[["b_%s"]] <- mu_%s + %s',
                                  attr_short[a], attr_short[a], inner))
      } else {
        lines <- c(lines, sprintf('  randcoeff[["b_price"]] <- -exp(mu_price + %s)', inner))
      }
    }
  }
  lines <- c(lines, '  return(randcoeff)', '}')
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

# apollo_probabilities builder (same for both MMNL flavours).
.make_apollo_prob_mmnl <- function(dgp = DGP_DEFAULT) {
  J <- dgp$n_alternatives; n_generic <- dgp$n_generic
  alt_entries   <- paste(sprintf('alt%d = %d', 1:J, 1:J), collapse = ", ")
  avail_entries <- paste(sprintf('alt%d = 1', 1:J), collapse = ", ")
  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    '  V <- list()'
  )
  for (j in 1:J) {
    terms <- c(if (j < J) sprintf('asc_alt%d', j),
               sprintf('b_x%d * x%d_%d', 1:n_generic, 1:n_generic, j),
               sprintf('b_price * price_%d', j))
    lines <- c(lines, sprintf('  V[["alt%d"]] <- %s', j, paste(terms, collapse = " + ")))
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
  environment(fn) <- asNamespace("apollo")
  fn
}

# One Apollo MMNL estimation with given draws and starting values.
# All components are passed to apollo_validateInputs as explicit arguments
# (officially supported; globalenv is only the fallback for missing ones).
# The single exception is apollo_probabilities: validateInputs' pre-processing
# looks it up in globalenv and otherwise fails for certain (J, n_generic)
# combinations, so it is assigned there BEFORE the call and removed afterwards
# by the caller's cleanup_apollo(). cleanup_apollo() at entry also clears any
# stray apollo_* objects from an interactive session, which would otherwise
# leak in through the globalenv fallback (e.g. a leftover apollo_lcPars).
.run_apollo_mmnl <- function(database, n_draws, start_beta,
                             dgp                = DGP_DEFAULT,
                             n_cores            = NULL,
                             draws_type         = DRAWS_TYPE_MMNL,
                             estimation_routine = ESTIMATION_ROUTINE_MMNL,
                             bounds             = NULL,
                             correlation        = FALSE) {
  cleanup_apollo()
  n_generic <- dgp$n_generic
  if (is.null(n_cores)) {
    n_cores <- getOption("klue.mmnl.n_cores", .klue_default_mmnl_cores())
  }
  n_cores <- max(1L, as.integer(n_cores))

  control <- list(
    modelName       = paste0(if (correlation) "MMNL_corr_" else "MMNL_sim_",
                             as.integer(Sys.time()) %% 100000,
                             "_", sample.int(10000, 1)),
    modelDescr      = if (correlation) "Correlated MMNL simulation" else "MMNL simulation",
    indivID         = "ID",
    nCores          = n_cores,
    mixing          = TRUE,
    outputDirectory = tempdir()
  )
  draws <- list(
    interDrawsType = draws_type,
    interNDraws    = as.integer(n_draws),
    interUnifDraws = c(),
    interNormDraws = c(paste0("draws_x", 1:n_generic), "draws_price"),
    intraDrawsType = draws_type,
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )
  probabilities <- .make_apollo_prob_mmnl(dgp)
  assign("apollo_probabilities", probabilities, envir = .GlobalEnv)

  inputs <- tryCatch(
    apollo_validateInputs(
      apollo_beta      = start_beta,
      apollo_fixed     = c(),
      database         = database,
      apollo_control   = control,
      apollo_draws     = draws,
      apollo_randCoeff = .make_apollo_randCoeff(dgp, correlation)
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl: apollo_validateInputs error: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(inputs)) return(NULL)

  est_settings <- list(
    estimationRoutine = estimation_routine,
    writeIter         = FALSE,
    silent            = TRUE
  )
  if (!is.null(bounds)) est_settings$bounds <- bounds

  tryCatch(
    apollo_estimate(
      apollo_beta          = start_beta,
      apollo_fixed         = c(),
      apollo_probabilities = probabilities,
      apollo_inputs        = inputs,
      estimate_settings    = est_settings
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl: apollo_estimate error: ",
              conditionMessage(e))
      NULL
    }
  )
}

#' MMNL benchmark (independent or correlated normals)
#'
#' Two-stage estimation: a cheap warm-start with `n_draws_stage1` draws, then
#' the main stage with `n_draws`. Independent flavour: normals on the generic
#' attributes, negative log-normal on price, optional box constraints on the
#' price parameters (pass NA to disable). Correlated flavour
#' (`correlation = TRUE`): full lower-triangular Cholesky covariance, warm-
#' started from an independent fit; the price bounds do not apply.
#'
#' @param database Data frame in long Apollo format with an `ID` column and the
#'   attribute and choice columns expected by the data-generating process.
#' @param correlation Logical; if `TRUE`, estimate the correlated specification
#'   with a full lower-triangular Cholesky covariance instead of independent
#'   normals.
#' @param n_draws Number of inter-individual draws for the main estimation
#'   stage. `NULL` uses the package default.
#' @param n_draws_stage1 Number of draws for the cheap warm-start stage. `NULL`
#'   uses the package default.
#' @param draws_type Apollo inter-draws type (for example Halton or MLHS).
#'   `NULL` uses the package default.
#' @param estimation_routine Apollo estimation routine passed to
#'   `apollo_estimate`. `NULL` uses the package default.
#' @param n_cores Number of cores for Apollo. `NULL` uses the package default.
#' @param quiet Logical; if `TRUE`, redirect Apollo output to a temporary log
#'   file whose tail is attached to a failing result. `NULL` uses the package
#'   default.
#' @param mu_price_bounds Length-2 numeric box constraint `c(lower, upper)` on
#'   the price mean parameter. `NA` disables the bound; `NULL` uses the package
#'   default. Ignored when `correlation = TRUE`.
#' @param sigma_price_bounds Length-2 numeric box constraint `c(lower, upper)`
#'   on the price spread parameter. `NA` disables the bound; `NULL` uses the
#'   package default. Ignored when `correlation = TRUE`.
#' @param dgp Data-generating-process specification giving the number of
#'   alternatives, generic attributes, and parameters.
#' @return A list describing the fit. On success: `converged = TRUE`, the
#'   log-likelihood `LL`, information criteria `BIC` and `AIC`, the number of
#'   free parameters `k`, `reason = "ok"`, the resolved `settings`, and, for the
#'   independent specification, the price-positive means `mu` and spreads
#'   `sigma`. On failure: `converged = FALSE` with `LL`, `BIC`, `AIC`, and `k`
#'   set to non-informative values, a `reason` string, and the Apollo log tail
#'   and path when `quiet` is `TRUE`.
#' @export
klue_mmnl <- function(database,
                      correlation        = FALSE,
                      n_draws            = NULL,
                      n_draws_stage1     = NULL,
                      draws_type         = NULL,
                      estimation_routine = NULL,
                      n_cores            = NULL,
                      quiet              = NULL,
                      mu_price_bounds    = NULL,
                      sigma_price_bounds = NULL,
                      dgp                = DGP_DEFAULT) {
  d <- klue_mmnl_defaults()
  if (is.null(n_draws))            n_draws            <- d$n_draws
  if (is.null(n_draws_stage1))     n_draws_stage1     <- d$n_draws_stage1
  if (is.null(draws_type))         draws_type         <- d$draws_type
  if (is.null(estimation_routine)) estimation_routine <- d$estimation_routine
  if (is.null(n_cores))            n_cores            <- d$n_cores
  if (is.null(quiet))              quiet              <- d$quiet
  if (!correlation) {
    # Bounds: NA disables; NULL means "use defaults".
    if (is.null(mu_price_bounds))    mu_price_bounds    <- d$mu_price_bounds
    if (is.null(sigma_price_bounds)) sigma_price_bounds <- d$sigma_price_bounds
    if (length(mu_price_bounds)    == 1 && is.na(mu_price_bounds))    mu_price_bounds    <- NULL
    if (length(sigma_price_bounds) == 1 && is.na(sigma_price_bounds)) sigma_price_bounds <- NULL
  } else {
    mu_price_bounds <- NULL; sigma_price_bounds <- NULL
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
      cleanup_apollo()
    }, add = TRUE)
  } else {
    on.exit(cleanup_apollo(), add = TRUE)
  }

  N <- length(unique(database$ID))
  n_beta <- dgp$n_beta; n_generic <- dgp$n_generic; J <- dgp$n_alternatives
  settings <- list(n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
                   draws_type = draws_type,
                   estimation_routine = estimation_routine,
                   n_cores = n_cores)
  if (!correlation) {
    settings$mu_price_bounds    <- mu_price_bounds
    settings$sigma_price_bounds <- sigma_price_bounds
  }

  read_log_tail <- function() {
    if (!isTRUE(quiet) || !file.exists(log_file)) return(NULL)
    out <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) NULL)
    if (length(out) == 0) NULL else utils::tail(out, 40L)
  }
  fail_with <- function(reason) {
    res <- list(converged = FALSE, LL = -Inf, BIC = Inf, AIC = Inf, k = 0,
                reason          = reason,
                apollo_log_tail = read_log_tail(),
                apollo_log_path = if (isTRUE(quiet) && file.exists(log_file)) log_file else NULL,
                settings        = settings)
    if (!correlation) {
      res$mu <- rep(0, n_beta); res$sigma <- rep(0, n_beta)
    }
    res
  }
  make_bounds <- function(beta_vec) {
    if (is.null(mu_price_bounds) && is.null(sigma_price_bounds)) return(NULL)
    lower <- rep(-Inf, length(beta_vec)); upper <- rep(Inf, length(beta_vec))
    names(lower) <- names(upper) <- names(beta_vec)
    if (!is.null(mu_price_bounds) && "mu_price" %in% names(beta_vec)) {
      lower["mu_price"] <- mu_price_bounds[1]; upper["mu_price"] <- mu_price_bounds[2]
    }
    if (!is.null(sigma_price_bounds) && "sigma_price" %in% names(beta_vec)) {
      lower["sigma_price"] <- sigma_price_bounds[1]; upper["sigma_price"] <- sigma_price_bounds[2]
    }
    list(lower = lower, upper = upper)
  }

  # ---- Starting values ------------------------------------------------------
  beta0 <- c()
  for (j in 1:(J - 1)) beta0[paste0("asc_alt", j)] <- 0
  if (!correlation) {
    # MNL-informed means (MNL = LCMNL with C=1)
    mnl_fit <- tryCatch(
      estimate_lcmnl(database, C = 1,
                     start_betas = matrix(0, nrow = 1, ncol = n_beta), dgp = dgp),
      error = function(e) NULL
    )
    if (!is.null(mnl_fit) && mnl_fit$converged) {
      mnl_b <- mnl_fit$betas[1, ]
      for (a in 1:n_generic) beta0[paste0("mu_x", a)] <- mnl_b[a]
      beta0["mu_price"] <- if (mnl_b[n_beta] < 0) log(-mnl_b[n_beta]) else 0.0
    } else {
      for (a in 1:n_generic) beta0[paste0("mu_x", a)] <- 0.5
      beta0["mu_price"] <- 0.0
    }
    for (a in 1:n_generic) beta0[paste0("sigma_x", a)] <- log(0.5)
    beta0["sigma_price"] <- log(0.3)
  } else {
    # Warm start from an independent fit at stage-1 draws.
    indep_fit <- klue_mmnl(database, correlation = FALSE,
                           n_draws = n_draws_stage1, n_draws_stage1 = 100L,
                           draws_type = draws_type,
                           estimation_routine = estimation_routine,
                           n_cores = n_cores, quiet = quiet, dgp = dgp)
    if (indep_fit$converged) {
      mu_starts  <- unname(indep_fit$mu)
      sig_starts <- log(unname(indep_fit$sigma))
    } else {
      mu_starts  <- c(rep(0.5, n_generic), 0.0)
      sig_starts <- rep(log(0.5), n_beta)
    }
    attr_short  <- c(paste0("x", 1:n_generic), "price")
    chol_prefix <- c(paste0("s_x", 1:n_generic), "s_pr")
    for (a in 1:n_beta) beta0[paste0("mu_", attr_short[a])] <- mu_starts[a]
    for (row in 1:n_beta) {
      for (col in 1:row) {
        beta0[paste0(chol_prefix[row], "_", col)] <-
          if (row == col) sig_starts[row] else 0
      }
    }
  }

  # ---- Stage 1: cheap warm start --------------------------------------------
  stage1 <- .run_apollo_mmnl(database, n_draws_stage1, beta0,
                             dgp = dgp, n_cores = n_cores,
                             draws_type = draws_type,
                             estimation_routine = estimation_routine,
                             bounds = make_bounds(beta0),
                             correlation = correlation)
  beta1 <- beta0
  if (!is.null(stage1) && !is.null(stage1$estimate)) {
    est_s1 <- stage1$estimate
    ok <- all(is.finite(est_s1))
    if (ok && !correlation) {
      # Reject stage-1 estimates pinned to a bound (failed convergence sign).
      bnd <- make_bounds(beta0)
      if (!is.null(bnd)) {
        rng <- bnd$upper - bnd$lower
        rng[!is.finite(rng)] <- 0
        slack <- 1e-3 * rng
        ok <- all(est_s1 >= bnd$lower + slack & est_s1 <= bnd$upper - slack,
                  na.rm = TRUE)
      }
    }
    if (ok) beta1 <- est_s1
  }

  # ---- Stage 2: main estimation ---------------------------------------------
  model <- .run_apollo_mmnl(database, n_draws, beta1,
                            dgp = dgp, n_cores = n_cores,
                            draws_type = draws_type,
                            estimation_routine = estimation_routine,
                            bounds = make_bounds(beta1),
                            correlation = correlation)
  if (is.null(model))                                return(fail_with("stage2_apollo_estimate_failed"))
  if (is.null(model$estimate) ||
      !all(is.finite(model$estimate)))               return(fail_with("stage2_non_finite_estimates"))
  if (is.null(model$LLout) || !is.finite(model$LLout[1]))
                                                     return(fail_with("stage2_non_finite_LL"))

  est    <- model$estimate
  LL     <- model$LLout[1]
  n_free <- length(est)
  res <- list(converged = TRUE, LL = LL,
              BIC = -2 * LL + n_free * log(N),
              AIC = -2 * LL + 2 * n_free, k = n_free,
              reason = "ok", apollo_log_tail = NULL, apollo_log_path = NULL,
              settings = settings)
  if (!correlation) {
    mu_names    <- c(paste0("mu_x", 1:n_generic), "mu_price")
    sigma_names <- c(paste0("sigma_x", 1:n_generic), "sigma_price")
    res$mu    <- est[mu_names]
    res$sigma <- exp(est[sigma_names])
  }
  res
}

#' @rdname klue_mmnl
#' @export
klue_mmnl_corr <- function(database,
                           n_draws            = NULL,
                           n_draws_stage1     = NULL,
                           draws_type         = NULL,
                           estimation_routine = NULL,
                           n_cores            = NULL,
                           quiet              = NULL,
                           dgp                = DGP_DEFAULT) {
  klue_mmnl(database, correlation = TRUE,
            n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
            draws_type = draws_type, estimation_routine = estimation_routine,
            n_cores = n_cores, quiet = quiet, dgp = dgp)
}
