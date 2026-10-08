#!/usr/bin/env Rscript
# =============================================================================
# dev/check_apollo_searchstart_gradients.R
#
# Why. version12.tex states, in the appendix on the Apollo arms, that in Apollo
# 0.3.5 apollo_searchStart builds the log-likelihood from the model code
# without the pre-processing that apollo_estimate applies to obtain analytic
# gradients, so that for models built on apollo_mnl and written as in Apollo's
# examples, a plain MNL included, its BFGS steps use central-difference
# gradients, each costing 2K + 1 log-likelihood evaluations for K parameters,
# and Apollo reports that it could not compute analytic gradients. The Speed
# paragraph (Results) and the notes to the blocked stress-ladder table rest on
# the same mechanism; the table notes say it of apollo_lcEM as well. This
# script checks these statements on tiny searches and records
# apollo_searchStart's default settings for README_REPRODUCE.md. It replaces
# two scratch checks run on 2026-10-05 (verify/indep_check.R and
# verify/mnl_searchstart_grad.R, kept outside git in
# code_backups/scratch_2026-10-06/session_ed878b37/).
#
# Mechanism, from the Apollo 0.3.5 source. apollo_estimate passes the model
# through apollo_modifyUserDefFunc, which wraps each utility in a function
# (and, for a latent class model, expands the class loop); Apollo's analytic
# gradient needs that form. apollo_searchStart calls apollo_makeLogLike and
# apollo_makeGrad on the function as written, so apollo_makeGrad returns NULL
# and prints "Apollo was not able to compute analytical gradients", and each
# candidate's BFGS run (maxLik::maxLik, method "bfgs") gets grad = NULL, i.e.
# maxLik::numericGradient: central differences, f at the point once and at
# point -/+ eps/2 for each of the K free parameters.
#
# Models, written below in the style of Apollo's example scripts
# (apollo_attach, a list V of utilities, apollo_mnl, apollo_panelProd,
# apollo_prepareProb; the latent class model with apollo_classAlloc and the
# literal class loop of LC_no_covariates.r; the EM form with an explicit-logit
# allocation, for (s in 1:length(pi_values)) and a weights column, as in
# EM_LC_no_covariates.r). The example scripts themselves are not copied.
#   MNL   simulated: 200 respondents x 6 tasks, 3 alternatives, 2 constants,
#         2 attributes and a price; K = 5 free parameters
#   MODE  Apollo's apollo_modeChoiceData, SP rows: 4 modes, constants (car
#         fixed), generic time and cost; K = 5
#   LC    simulated, 2 classes: 200 respondents x 8 tasks, class-specific
#         constants, attributes and price, delta_b fixed; K = 11
#   EM    the LC data and model in the EM example's form; K = 11
# Tests, one CSV row each (E1: one row per kind of maxLik call):
#   M1  apollo_searchStart on the MNL as written      expected: numerical, message
#   M2  M1 with apollo_control$analyticGrad = TRUE set by hand
#   M3  M1 with debug = TRUE: Apollo prints its component gradient table
#   M4  utilities hand-written as function() closures  analytic
#   M5  model pre-processed by apollo_modifyUserDefFunc first, as
#       apollo_estimate does                           analytic, same LLs as M1
#   M6  apollo_estimate on the MNL as written          analytic (BGW)
#   M7  apollo_estimate with noModification = TRUE      numerical, message
#   P1  apollo_searchStart on the mode-choice MNL      numerical, message
#   P2  P1 after pre-processing                        analytic, same LLs as P1
#   P3  apollo_estimate on the mode-choice MNL         analytic (BGW)
#   L0  per-component gradient flags of the LC log-likelihood that
#       apollo_searchStart builds, as written and pre-processed
#   L1  apollo_searchStart on the LC as written        numerical, message
#   L2a L1 after pre-processing, global apollo_beta = the start   analytic
#   L2b the same with the global apollo_beta moved off the start: Apollo
#       checks the expanded loop at the global apollo_beta and keeps the
#       unexpanded one                                 numerical
#   L3  apollo_estimate on the LC as written           analytic (BGW)
#   E1  apollo_lcEM on the EM form: its M-steps (allocation, K = 1; each
#       class, K = 5) and its closing maximum-likelihood step (K = 11)
#                                                      numerical
# plus one row per default setting of apollo_searchStart (Apollo 0.3.5 sets
# them in the function body; the formal default is searchStart_settings = NA).
#
# Instrumentation, by trace() in this session only (nothing is installed):
#   maxLik::maxLik            is grad a function on entry; the log-likelihood
#                             passed in is wrapped so every evaluation is counted
#   maxLik::numericGradient   evaluations made inside each call (inside maxLik)
#   bgw::bgw_mle              is calcJ, the analytic gradient, a function
# A row describes the optimiser runs that the tested function starts itself
# (its BFGS runs, else its BGW run); BFGS runs that other Apollo functions start
# inside it, such as apollo_mnl's LL(c) statistic on the mode-choice data, are
# listed in `note`. Columns: ll_calls_per_gradient = evaluations inside each
# numericGradient call (all values seen); matches_2K_plus_1 = every call cost
# exactly 2K + 1; ll_calls_outside_gradients = the remaining evaluations of
# the BFGS runs (line search, plus one evaluation per run that maxLik spends
# checking for a gradient attribute when no gradient function is given).
#
# Settings. Searches: nCandidates = 2, maxStages = 1, bfgsIter = 5 (all other
# settings at Apollo's defaults). apollo_estimate: writeIter = FALSE,
# hessianRoutine = "none", maxIterations = 100. apollo_lcEM:
# EMmaxIterations = 100 (the example's setting) and EMstoppingCriterion = 1,
# so EM stops after a few iterations and the closing step runs; closing step
# without a Hessian. Apollo's files go to a folder in R's tempdir().
#
# Seeds. Data: set.seed(20261005) before the MNL data, set.seed(20261006)
# before the LC data. apollo_control$seed = 13 (Apollo's default, set
# explicitly), so apollo_searchStart draws its candidates after set.seed(17)
# and the raw and pre-processed searches start from the same candidates.
# Deterministic: counts and LLs are the same on every run; seconds are not.
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/check_apollo_searchstart_gradients.R
# About 10 s on one core (measured 2026-10-06), well under 1 GB. Stops unless
# Apollo is the pinned 0.3.5 (dev/pinned_versions.R); the maxLik and bgw
# versions are recorded.
# Writes: output/apollo_searchstart_gradients.csv.
# =============================================================================

if (!file.exists("dev/pinned_versions.R")) stop("run from the repository root")
options(width = 160, mc.cores = 1L)
OUT_CSV <- "output/apollo_searchstart_gradients.csv"
OUTDIR  <- file.path(tempdir(), "apollo_searchstart_gradients")   # Apollo's own files
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
SS  <- list(nCandidates = 2, maxStages = 1, bfgsIter = 5)
EST <- list(writeIter = FALSE, hessianRoutine = "none", silent = FALSE, maxIterations = 100)

T0 <- proc.time()[["elapsed"]]
stamp <- function(...) { cat(sprintf("[%6.1f s] %s\n", proc.time()[["elapsed"]] - T0,
                                     sprintf(...))); flush.console() }

suppressPackageStartupMessages({ library(apollo); loadNamespace("maxLik"); invisible(loadNamespace("bgw")) })
source("dev/pinned_versions.R", local = TRUE)
check_pinned_versions("apollo")
VERSIONS <- c(apollo = as.character(packageVersion("apollo")),
              maxLik = as.character(packageVersion("maxLik")),
              bgw    = as.character(packageVersion("bgw")),
              R      = paste(R.version$major, R.version$minor, sep = "."))
stamp("apollo %s | maxLik %s | bgw %s | R %s", VERSIONS[["apollo"]], VERSIONS[["maxLik"]],
      VERSIONS[["bgw"]], VERSIONS[["R"]])

# ---- instrumentation -----------------------------------------------------------
# .TR is a dot-name, so Apollo's checks for global objects used inside the
# model functions do not see it.
.TR <- new.env()
.TR$reset <- function(test) {
  .TR$test <- test; .TR$calls <- list(); .TR$cur <- NULL
  .TR$in_ng <- FALSE; .TR$ng_n <- 0L; .TR$bgw <- logical(0)
}
.TR$reset("")
invisible(trace("maxLik", where = asNamespace("maxLik"), print = FALSE,
  tracer = quote({
    .TR$cur <- new.env()
    .TR$cur$K <- length(start); .TR$cur$analytic <- is.function(grad)
    .TR$cur$n_ll <- 0L; .TR$cur$ng <- integer(0)
    # the innermost Apollo function on the call stack: who ran this BFGS
    .tr_fn <- vapply(sys.calls(), function(cl) {
      f <- cl[[1]]
      if (is.name(f)) as.character(f)
      else if (is.call(f) && as.character(f[[1]]) %in% c("::", ":::")) as.character(f[[3]])
      else ""
    }, "")
    .tr_fn <- .tr_fn[grepl("^apollo_", .tr_fn) & .tr_fn != "apollo_logLike"]
    .TR$cur$caller <- if (length(.tr_fn)) .tr_fn[length(.tr_fn)] else "other"
    .tr_ll0 <- logLik
    logLik <- function(...) {
      .TR$cur$n_ll <- .TR$cur$n_ll + 1L
      if (.TR$in_ng) .TR$ng_n <- .TR$ng_n + 1L
      .tr_ll0(...)
    }
  }),
  exit = quote({
    .tr_rv <- tryCatch(returnValue(NULL), error = function(e) NULL)
    .TR$calls[[length(.TR$calls) + 1L]] <- list(
      K = .TR$cur$K, analytic = .TR$cur$analytic, n_ll = .TR$cur$n_ll, ng = .TR$cur$ng,
      caller = .TR$cur$caller,
      maximum = if (is.null(.tr_rv$maximum)) NA_real_ else as.numeric(.tr_rv$maximum))
    .TR$cur <- NULL
  })))
invisible(trace("numericGradient", where = asNamespace("maxLik"), print = FALSE,
  tracer = quote(if (!is.null(.TR$cur)) { .TR$in_ng <- TRUE; .TR$ng_n <- 0L }),
  exit = quote(if (!is.null(.TR$cur) && .TR$in_ng) {
    .TR$cur$ng <- c(.TR$cur$ng, .TR$ng_n); .TR$in_ng <- FALSE })))
invisible(trace("bgw_mle", where = asNamespace("bgw"), print = FALSE,
  tracer = quote(.TR$bgw <- c(.TR$bgw, is.function(calcJ)))))
stopifnot(inherits(get("maxLik", asNamespace("maxLik")), "functionWithTrace"),
          inherits(get("numericGradient", asNamespace("maxLik")), "functionWithTrace"),
          inherits(get("bgw_mle", asNamespace("bgw")), "functionWithTrace"))

# Run an Apollo call, keep its console output and messages, time it.
KEYPAT <- paste("not able to compute analytical", "numerical gradients will be used",
                "Analytical gradient is different", "ComponentName", "MNL +MNL",
                "analytic model derivatives", "EM algorithm stopped", sep = "|")
run_captured <- function(test, expr) {
  .TR$reset(test); msgs <- character(0)
  secs <- system.time(out <- utils::capture.output(
    val <- tryCatch(withCallingHandlers(expr,
             message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
             warning = function(w) { msgs <<- c(msgs, paste("WARNING:", conditionMessage(w)))
                                     invokeRestart("muffleWarning") }),
           error = function(e) structure(conditionMessage(e), class = "check_error")),
    type = "output"))[["elapsed"]]
  txt <- c(out, msgs)
  one <- gsub("[[:space:]]+", " ", paste(txt, collapse = " "))   # undo Apollo's line wrapping
  key <- unique(trimws(grep(KEYPAT, txt, value = TRUE)))
  stamp("%s (%.2f s)%s", test, secs,
        if (length(key)) paste0("\n    | ", paste(key, collapse = "\n    | ")) else "")
  if (inherits(val, "check_error")) stop(test, " failed: ", val)
  list(value = val, secs = secs, calls = .TR$calls, bgw = .TR$bgw, key = key,
       message = grepl("not able to compute analytical gradients", one, fixed = TRUE),
       numgrad_message = grepl("numerical gradients will be used", one, ignore.case = TRUE))
}

# One CSV row per test (what a row describes: see the header). LL: the end LL
# of each BFGS run (searches: one per candidate) or apollo_estimate's final LL.
ROWS <- list(); LLS <- list()
aux_note <- function(aux) {
  if (!length(aux)) return("")
  g <- split(aux, vapply(aux, function(cl) sprintf("%s|%d", cl$caller, cl$K), ""))
  paste0("other BFGS runs inside: ", paste(vapply(g, function(x) {
    ng <- unlist(lapply(x, `[[`, "ng"))
    sprintf("%d from %s (K = %d, grad %s, %s LL evaluations per gradient%s)", length(x),
            x[[1]]$caller, x[[1]]$K,
            if (all(vapply(x, `[[`, logical(1), "analytic"))) "analytic" else "NULL",
            if (length(ng)) paste(unique(ng), collapse = "/") else "no numerical",
            if (x[[1]]$caller == "apollo_mnl")
              "; apollo_mnl's LL(c), a constants-only fit run after estimation when availability varies"
            else "")
  }, ""), collapse = "; "))
}
add_test <- function(test, model, fun, K, r, note = "", calls = NULL) {
  if (is.null(calls)) {
    own   <- vapply(r$calls, `[[`, "", "caller") == fun
    calls <- r$calls[own]
    aux   <- aux_note(r$calls[!own])
    if (nzchar(aux)) note <- paste0(note, "; ", aux)
  }
  ml <- length(calls) > 0
  ng <- unlist(lapply(calls, `[[`, "ng"))
  ok <- if (length(ng)) all(unlist(lapply(calls, function(cl) cl$ng == 2L * cl$K + 1L))) else NA
  n_ll <- sum(vapply(calls, `[[`, integer(1), "n_ll"))
  an <- if (ml) unique(vapply(calls, `[[`, logical(1), "analytic")) else unique(r$bgw)
  LL <- if (ml) vapply(calls, `[[`, numeric(1), "maximum")
        else if (is.list(r$value) && !is.null(r$value$maximum)) r$value$maximum else NA_real_
  LLS[[test]] <<- LL
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    test = test, model = model, function_called = fun, K = K,
    optimizer = if (ml) "maxLik BFGS" else if (length(r$bgw)) "BGW" else "none",
    grad_function_passed = if (length(an) == 1L) an else NA,
    n_optimizer_calls = if (ml) length(calls) else length(r$bgw),
    n_numeric_gradients = length(ng),
    ll_calls_per_gradient = if (length(ng)) paste(sort(unique(ng)), collapse = "/") else "",
    matches_2K_plus_1 = ok,
    ll_calls_total = if (ml) n_ll else NA_integer_,
    ll_calls_outside_gradients = if (ml) n_ll - sum(ng) else NA_integer_,
    message_printed = r$message, numgrad_fallback_message = r$numgrad_message,
    LL = paste(sprintf("%.3f", LL), collapse = "/"), seconds = round(r$secs, 2),
    setting = "", default_value = "", note = note, stringsAsFactors = FALSE)
}

# Apollo reads its objects from the global environment, as in the examples;
# apollo_validateInputs also picks up a global apollo_lcPars, so the MNL
# set-ups remove it.
setup_apollo <- function(name, db, beta, fixed, extra_control = list(), lcPars = NULL) {
  apollo_control <<- c(list(modelName = name, modelDescr = name, indivID = "ID", nCores = 1,
                            outputDirectory = OUTDIR, seed = 13), extra_control)
  database     <<- db
  apollo_beta  <<- beta
  apollo_fixed <<- fixed
  if (is.null(lcPars)) {
    if (exists("apollo_lcPars", envir = globalenv(), inherits = FALSE))
      rm("apollo_lcPars", envir = globalenv())
  } else apollo_lcPars <<- lcPars
  invisible(capture.output(apollo_inputs <<- apollo_validateInputs(
    apollo_beta = apollo_beta, apollo_fixed = apollo_fixed, database = database,
    apollo_control = apollo_control, apollo_lcPars = if (is.null(lcPars)) NA else lcPars,
    silent = FALSE)))
}
# apollo_estimate's pre-processing (apollo_modifyUserDefFunc), applied before
# apollo_searchStart as in dev/run_benchmark_apollo_defaults.R: the
# pre-processed function and a copy of apollo_inputs with the pre-processed
# lcPars and scaling. Apollo checks the loop expansion at the global
# apollo_beta, so `beta` is passed explicitly and the caller decides what the
# global holds (test L2b).
preprocess <- function(beta, prob, inputs) {
  L <- NULL
  invisible(capture.output(L <- apollo_modifyUserDefFunc(beta, apollo_fixed, prob, inputs,
                                                         validate = TRUE, noModification = FALSE)))
  if (!isTRUE(L$success)) stop("apollo_modifyUserDefFunc failed")
  inputs$apollo_lcPars  <- L$apollo_lcPars
  inputs$apollo_scaling <- L$apollo_scaling; inputs$manualScaling <- L$manualScaling
  body_txt <- deparse(L$apollo_probabilities)
  list(prob = L$apollo_probabilities, inputs = inputs,
       wrapped = any(grepl("function()", body_txt, fixed = TRUE)),
       loop_left = any(grepl("for (s in", body_txt, fixed = TRUE)))
}

# ---- MNL, simulated data ---------------------------------------------------------
set.seed(20261005)
MNL_DB <- local({
  n_id <- 200; n_t <- 6
  d <- data.frame(ID = rep(1:n_id, each = n_t), task = rep(1:n_t, n_id)); n <- nrow(d)
  for (j in 1:3) {
    d[[paste0("x1_", j)]] <- rnorm(n); d[[paste0("x2_", j)]] <- rnorm(n)
    d[[paste0("p_", j)]] <- runif(n, 1, 3)
  }
  U <- sapply(1:3, function(j) c(0.3, -0.2, 0)[j] + 0.8 * d[[paste0("x1_", j)]] -
                0.5 * d[[paste0("x2_", j)]] - 0.7 * d[[paste0("p_", j)]] - log(-log(runif(n))))
  d$choice <- max.col(U)
  d
})

# As Apollo's MNL examples write it.
mnl_prob_raw <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  V <- list()
  V[["alt1"]] <- asc_1 + b_x1 * x1_1 + b_x2 * x2_1 + b_p * p_1
  V[["alt2"]] <- asc_2 + b_x1 * x1_2 + b_x2 * x2_2 + b_p * p_2
  V[["alt3"]] <-         b_x1 * x1_3 + b_x2 * x2_3 + b_p * p_3
  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
                       avail = list(alt1 = 1, alt2 = 1, alt3 = 1), choiceVar = choice,
                       utilities = V)
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}
# The same model with the utilities hand-written as function() closures, the
# form apollo_modifyUserDefFunc produces.
mnl_prob_fun <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  V <- list()
  V[["alt1"]] <- function() asc_1 + b_x1 * x1_1 + b_x2 * x2_1 + b_p * p_1
  V[["alt2"]] <- function() asc_2 + b_x1 * x1_2 + b_x2 * x2_2 + b_p * p_2
  V[["alt3"]] <- function()         b_x1 * x1_3 + b_x2 * x2_3 + b_p * p_3
  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
                       avail = list(alt1 = 1, alt2 = 1, alt3 = 1), choiceVar = choice,
                       utilities = V)
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}
MNL_BETA <- c(asc_1 = 0, asc_2 = 0, b_x1 = 0, b_x2 = 0, b_p = 0)
MNL_LAB  <- "MNL, simulated (200 x 6 choices, 3 alternatives)"

setup_apollo("chk_M1", MNL_DB, MNL_BETA, character(0))
r <- run_captured("M1", apollo_searchStart(apollo_beta, apollo_fixed, mnl_prob_raw, apollo_inputs, SS))
add_test("M1_searchStart_raw", MNL_LAB, "apollo_searchStart", 5L, r,
         "model as written; apollo_control$analyticGrad TRUE (Apollo's default)")

setup_apollo("chk_M2", MNL_DB, MNL_BETA, character(0), list(analyticGrad = TRUE))
r <- run_captured("M2", apollo_searchStart(apollo_beta, apollo_fixed, mnl_prob_raw, apollo_inputs, SS))
add_test("M2_searchStart_raw_analyticGrad_TRUE", MNL_LAB, "apollo_searchStart", 5L, r,
         sprintf("analyticGrad = TRUE set by hand (analyticGrad_manualSet %s)",
                 isTRUE(apollo_inputs$apollo_control$analyticGrad_manualSet)))

setup_apollo("chk_M3", MNL_DB, MNL_BETA, character(0), list(debug = TRUE))
r <- run_captured("M3", apollo_searchStart(apollo_beta, apollo_fixed, mnl_prob_raw, apollo_inputs, SS))
tab_line <- grep("^MNL +MNL", r$key, value = TRUE)
add_test("M3_searchStart_raw_debug", MNL_LAB, "apollo_searchStart", 5L, r,
         paste("debug = TRUE; Apollo's component table (name, type, gradient):",
               if (length(tab_line)) gsub(" +", " ", tab_line[1]) else "not printed"))

setup_apollo("chk_M4", MNL_DB, MNL_BETA, character(0))
r <- run_captured("M4", apollo_searchStart(apollo_beta, apollo_fixed, mnl_prob_fun, apollo_inputs, SS))
add_test("M4_searchStart_function_utilities", MNL_LAB, "apollo_searchStart", 5L, r,
         "utilities hand-written as function() closures")

setup_apollo("chk_M5", MNL_DB, MNL_BETA, character(0))
pp <- preprocess(apollo_beta, mnl_prob_raw, apollo_inputs)
r <- run_captured("M5", apollo_searchStart(apollo_beta, apollo_fixed, pp$prob, pp$inputs, SS))
add_test("M5_searchStart_preprocessed", MNL_LAB, "apollo_searchStart", 5L, r,
         sprintf("apollo_modifyUserDefFunc first (utilities wrapped as functions: %s)", pp$wrapped))

setup_apollo("chk_M6", MNL_DB, MNL_BETA, character(0))
r <- run_captured("M6", apollo_estimate(apollo_beta, apollo_fixed, mnl_prob_raw, apollo_inputs, EST))
add_test("M6_estimate_raw", MNL_LAB, "apollo_estimate", 5L, r, "model as written")

setup_apollo("chk_M7", MNL_DB, MNL_BETA, character(0), list(noModification = TRUE))
r <- run_captured("M7", apollo_estimate(apollo_beta, apollo_fixed, mnl_prob_raw, apollo_inputs, EST))
add_test("M7_estimate_noModification", MNL_LAB, "apollo_estimate", 5L, r,
         "noModification = TRUE: apollo_estimate skips the pre-processing")

# ---- MNL, Apollo's mode-choice data (SP rows) --------------------------------------
MODE_DB <- subset(apollo_modeChoiceData, apollo_modeChoiceData$SP == 1)
mode_prob_raw <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  V <- list()
  V[["car"]]  <- asc_car  + b_tt * time_car  + b_cost * cost_car
  V[["bus"]]  <- asc_bus  + b_tt * time_bus  + b_cost * cost_bus
  V[["air"]]  <- asc_air  + b_tt * time_air  + b_cost * cost_air
  V[["rail"]] <- asc_rail + b_tt * time_rail + b_cost * cost_rail
  mnl_settings <- list(alternatives = c(car = 1, bus = 2, air = 3, rail = 4),
                       avail = list(car = av_car, bus = av_bus, air = av_air, rail = av_rail),
                       choiceVar = choice, utilities = V)
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}
MODE_BETA <- c(asc_car = 0, asc_bus = 0, asc_air = 0, asc_rail = 0, b_tt = 0, b_cost = 0)
MODE_LAB  <- sprintf("MNL, apollo_modeChoiceData SP rows (%d choices, 4 modes)", nrow(MODE_DB))

setup_apollo("chk_P1", MODE_DB, MODE_BETA, "asc_car")
r <- run_captured("P1", apollo_searchStart(apollo_beta, apollo_fixed, mode_prob_raw, apollo_inputs, SS))
add_test("P1_searchStart_raw", MODE_LAB, "apollo_searchStart", 5L, r, "model as written")

setup_apollo("chk_P2", MODE_DB, MODE_BETA, "asc_car")
pp <- preprocess(apollo_beta, mode_prob_raw, apollo_inputs)
r <- run_captured("P2", apollo_searchStart(apollo_beta, apollo_fixed, pp$prob, pp$inputs, SS))
add_test("P2_searchStart_preprocessed", MODE_LAB, "apollo_searchStart", 5L, r,
         sprintf("apollo_modifyUserDefFunc first (utilities wrapped as functions: %s)", pp$wrapped))

setup_apollo("chk_P3", MODE_DB, MODE_BETA, "asc_car")
r <- run_captured("P3", apollo_estimate(apollo_beta, apollo_fixed, mode_prob_raw, apollo_inputs, EST))
add_test("P3_estimate_raw", MODE_LAB, "apollo_estimate", 5L, r, "model as written")

# ---- latent class, simulated data -------------------------------------------------
set.seed(20261006)
LC_DB <- local({
  n_id <- 200; n_t <- 8
  d <- data.frame(ID = rep(1:n_id, each = n_t)); n <- nrow(d)
  for (j in 1:3) {
    d[[paste0("x1_", j)]] <- rnorm(n); d[[paste0("x2_", j)]] <- rnorm(n)
    d[[paste0("p_", j)]] <- runif(n, 1, 3)
  }
  cls <- rep(ifelse(runif(n_id) < 0.6, 1L, 2L), each = n_t)
  B <- rbind(c(0.2, -0.1, 1.2, -0.2, -0.3), c(-0.3, 0.2, 0.1, -1.0, -1.5))
  U <- sapply(1:3, function(j) {
    b <- B[cls, , drop = FALSE]
    (if (j == 1) b[, 1] else if (j == 2) b[, 2] else 0) + b[, 3] * d[[paste0("x1_", j)]] +
      b[, 4] * d[[paste0("x2_", j)]] + b[, 5] * d[[paste0("p_", j)]] - log(-log(runif(n)))
  })
  d$choice <- max.col(U)
  d
})
LC_START <- c(asc_1_a = 0.1, asc_1_b = -0.1, asc_2_a = 0, asc_2_b = 0.1, b_x1_a = 0.8,
              b_x1_b = 0.2, b_x2_a = -0.1, b_x2_b = -0.8, b_p_a = -0.4, b_p_b = -1.2,
              delta_a = 0.3, delta_b = 0)
LC_LAB <- "latent class, 2 classes, simulated (200 x 8 choices)"

# As in Apollo's LC_no_covariates.r: apollo_classAlloc, literal class loop.
lc_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["asc_1"]] <- list(asc_1_a, asc_1_b)
  lcpars[["asc_2"]] <- list(asc_2_a, asc_2_b)
  lcpars[["b_x1"]]  <- list(b_x1_a, b_x1_b)
  lcpars[["b_x2"]]  <- list(b_x2_a, b_x2_b)
  lcpars[["b_p"]]   <- list(b_p_a, b_p_b)
  V <- list()
  V[["class_a"]] <- delta_a
  V[["class_b"]] <- delta_b
  classAlloc_settings <- list(classes = c(class_a = 1, class_b = 2), utilities = V)
  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)
  return(lcpars)
}
lc_prob_raw <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
                       avail = list(alt1 = 1, alt2 = 1, alt3 = 1), choiceVar = choice)
  for (s in 1:2) {
    V <- list()
    V[["alt1"]] <- asc_1[[s]] + b_x1[[s]] * x1_1 + b_x2[[s]] * x2_1 + b_p[[s]] * p_1
    V[["alt2"]] <- asc_2[[s]] + b_x1[[s]] * x1_2 + b_x2[[s]] * x2_2 + b_p[[s]] * p_2
    V[["alt3"]] <-              b_x1[[s]] * x1_3 + b_x2[[s]] * x2_3 + b_p[[s]] * p_3
    mnl_settings$utilities <- V
    P[[paste0("Class_", s)]] <- apollo_mnl(mnl_settings, functionality)
    P[[paste0("Class_", s)]] <- apollo_panelProd(P[[paste0("Class_", s)]], apollo_inputs, functionality)
  }
  lc_settings <- list(inClassProb = P, classProb = pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# L0: the per-component gradient flags that apollo_makeLogLike stores (the
# call apollo_searchStart makes), for the model as written and pre-processed.
component_flags <- function(ll) {
  ai <- environment(ll)$apollo_inputs
  nm <- grep("_settings$", names(ai), value = TRUE)
  f <- vapply(nm, function(n) isTRUE(ai[[n]]$gradient), logical(1))
  for (n in nm) if (is.list(ai[[n]]$classAlloc_settings) && !is.null(ai[[n]]$classAlloc_settings$gradient))
    f[paste0(n, "$classAlloc_settings")] <- isTRUE(ai[[n]]$classAlloc_settings$gradient)
  paste(sprintf("%s %s", sub("_settings", "", names(f)), ifelse(f, "analytic", "numeric")), collapse = ", ")
}
setup_apollo("chk_L0", LC_DB, LC_START, "delta_b", lcPars = lc_lcPars)
pp <- preprocess(apollo_beta, lc_prob_raw, apollo_inputs)
l0 <- tryCatch({
  ll_raw <- ll_mod <- NULL
  invisible(capture.output({
    ll_raw <- apollo_makeLogLike(apollo_beta, apollo_fixed, lc_prob_raw, apollo_inputs,
                                 list(estimationRoutine = "BHHH"))
    ll_mod <- apollo_makeLogLike(apollo_beta, apollo_fixed, pp$prob, pp$inputs,
                                 list(estimationRoutine = "BHHH"))
  }))
  b0 <- apollo_beta[!(names(apollo_beta) %in% apollo_fixed)]
  sprintf("as written: %s | pre-processed: %s | class loop expanded %s, utilities wrapped %s | LL at the start %.4f and %.4f",
          component_flags(ll_raw), component_flags(ll_mod), !pp$loop_left, pp$wrapped,
          ll_raw(b0, sumLL = TRUE), ll_mod(b0, sumLL = TRUE))
}, error = function(e) paste("flags not readable:", conditionMessage(e)))
stamp("L0 %s", l0)
ROWS[[length(ROWS) + 1L]] <- data.frame(
  test = "L0_component_gradient_flags", model = LC_LAB, function_called = "apollo_makeLogLike",
  K = 11L, optimizer = "none", grad_function_passed = NA, n_optimizer_calls = 0L,
  n_numeric_gradients = 0L, ll_calls_per_gradient = "", matches_2K_plus_1 = NA,
  ll_calls_total = NA_integer_, ll_calls_outside_gradients = NA_integer_,
  message_printed = NA, numgrad_fallback_message = NA, LL = "", seconds = NA_real_,
  setting = "", default_value = "", note = l0, stringsAsFactors = FALSE)

setup_apollo("chk_L1", LC_DB, LC_START, "delta_b", lcPars = lc_lcPars)
r <- run_captured("L1", apollo_searchStart(apollo_beta, apollo_fixed, lc_prob_raw, apollo_inputs, SS))
add_test("L1_searchStart_raw", LC_LAB, "apollo_searchStart", 11L, r, "model as written")

setup_apollo("chk_L2a", LC_DB, LC_START, "delta_b", lcPars = lc_lcPars)
pp <- preprocess(LC_START, lc_prob_raw, apollo_inputs)
r <- run_captured("L2a", apollo_searchStart(LC_START, apollo_fixed, pp$prob, pp$inputs, SS))
add_test("L2a_searchStart_preprocessed", LC_LAB, "apollo_searchStart", 11L, r,
         sprintf("apollo_modifyUserDefFunc first, global apollo_beta = the start (loop expanded %s)",
                 !pp$loop_left))

setup_apollo("chk_L2b", LC_DB, LC_START, "delta_b", lcPars = lc_lcPars)
apollo_beta <- LC_START + 0.37; apollo_beta["delta_b"] <- 0     # global moved off the start
pp <- preprocess(LC_START, lc_prob_raw, apollo_inputs)
r <- run_captured("L2b", apollo_searchStart(LC_START, apollo_fixed, pp$prob, pp$inputs, SS))
add_test("L2b_searchStart_preprocessed_global_beta_moved", LC_LAB, "apollo_searchStart", 11L, r,
         sprintf("as L2a with the global apollo_beta = start + 0.37 (loop expanded %s)", !pp$loop_left))

setup_apollo("chk_L3", LC_DB, LC_START, "delta_b", lcPars = lc_lcPars)
r <- run_captured("L3", apollo_estimate(apollo_beta, apollo_fixed, lc_prob_raw, apollo_inputs, EST))
add_test("L3_estimate_raw", LC_LAB, "apollo_estimate", 11L, r, "model as written")

# ---- apollo_lcEM, the EM example's form ---------------------------------------------
# Explicit-logit allocation and the loop for (s in 1:length(pi_values)), as in
# EM_LC_no_covariates.r; apollo_lcEM needs that loop form (dev/compare_packages.R).
lc_em_lcPars <- function(apollo_beta, apollo_inputs) {
  lcpars <- list()
  lcpars[["asc_1"]] <- list(asc_1_a, asc_1_b)
  lcpars[["asc_2"]] <- list(asc_2_a, asc_2_b)
  lcpars[["b_x1"]]  <- list(b_x1_a, b_x1_b)
  lcpars[["b_x2"]]  <- list(b_x2_a, b_x2_b)
  lcpars[["b_p"]]   <- list(b_p_a, b_p_b)
  lcpars[["pi_values"]] <- list(class_a = 1 / (1 + exp(delta_b - delta_a)),
                                class_b = 1 / (1 + exp(delta_a - delta_b)))
  return(lcpars)
}
lc_em_prob <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P <- list()
  mnl_settings <- list(alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
                       avail = list(alt1 = 1, alt2 = 1, alt3 = 1), choiceVar = choice)
  for (s in 1:length(pi_values)) {
    V <- list()
    V[["alt1"]] <- asc_1[[s]] + b_x1[[s]] * x1_1 + b_x2[[s]] * x2_1 + b_p[[s]] * p_1
    V[["alt2"]] <- asc_2[[s]] + b_x1[[s]] * x1_2 + b_x2[[s]] * x2_2 + b_p[[s]] * p_2
    V[["alt3"]] <-              b_x1[[s]] * x1_3 + b_x2[[s]] * x2_3 + b_p[[s]] * p_3
    mnl_settings$utilities <- V
    P[[paste0("Class_", s)]] <- apollo_mnl(mnl_settings, functionality)
    P[[paste0("Class_", s)]] <- apollo_panelProd(P[[paste0("Class_", s)]], apollo_inputs, functionality)
  }
  lc_settings <- list(inClassProb = P, classProb = pi_values)
  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}
EM_DB <- LC_DB; EM_DB$weights <- 1                 # the EM example's weights column
setup_apollo("chk_E1", EM_DB, LC_START, "delta_b",
             list(noValidation = TRUE, noDiagnostics = TRUE), lcPars = lc_em_lcPars)
r <- run_captured("E1", apollo_lcEM(apollo_beta, apollo_fixed, lc_em_prob, apollo_inputs,
                                    lcEM_settings = list(EMmaxIterations = 100, EMstoppingCriterion = 1),
                                    estimate_settings = list(writeIter = FALSE, silent = TRUE,
                                                             hessianRoutine = "none")))
em_note <- sprintf("%s; iterations: %s; its preliminary apollo_estimate ran %d BGW call(s), analytic gradient: %s",
                   if (any(grepl("improvements in LL smaller", r$key))) "EM stopped on the criterion, so the closing ML step ran"
                   else "EM stopped at the cap, no closing ML step",
                   if (is.list(r$value)) r$value$nIter else "?", length(r$bgw),
                   paste(unique(r$bgw), collapse = "/"))
# BFGS runs by who started them: apollo_lcEM's M-steps (class allocation,
# K = 1; one class at a time, K = 5) and the closing ML step, which
# apollo_lcEM hands to apollo_estimate with estimationRoutine "bfgs" (K = 11).
E1_GROUPS <- c(Mstep_allocation = "apollo_lcEM|1", Mstep_class = "apollo_lcEM|5",
               closing_ML = "apollo_estimate|11")
cl_key <- vapply(r$calls, function(cl) sprintf("%s|%d", cl$caller, cl$K), "")
other <- aux_note(r$calls[!(cl_key %in% E1_GROUPS)])
for (what in names(E1_GROUPS)) {
  sel <- cl_key == E1_GROUPS[[what]]
  if (!any(sel)) { stamp("E1: no BFGS run of kind %s", what); next }
  add_test(paste0("E1_lcEM_", what), paste(LC_LAB, "in the EM example's form"), "apollo_lcEM",
           as.integer(sub(".*\\|", "", E1_GROUPS[[what]])), r,
           note = paste0(sprintf("BFGS runs started by %s; ", sub("\\|.*", "", E1_GROUPS[[what]])),
                         em_note, if (nzchar(other)) paste0("; ", other) else ""),
           calls = r$calls[sel])
  # an M-step maximises a weighted class log-likelihood, not the model's: no LL
  if (what != "closing_ML") ROWS[[length(ROWS)]]$LL <- ""
}

# ---- apollo_searchStart's default settings ------------------------------------------
fm <- vapply(formals(apollo::apollo_searchStart), function(v) paste(deparse(v), collapse = " "), "")
fm <- fm[nzchar(fm)]                              # formals that have a default
def <- NULL
for (e in as.list(body(apollo::apollo_searchStart))[-1])
  if (is.call(e) && as.character(e[[1]]) %in% c("<-", "=") && identical(e[[2]], as.name("default"))) {
    def <- e[[3]]; break
  }
stopifnot(!is.null(def), identical(def[[1]], as.name("list")))
for (nm in names(as.list(def))[-1])
  ROWS[[length(ROWS) + 1L]] <- data.frame(
    test = "default_setting", model = "", function_called = "apollo_searchStart", K = NA_integer_,
    optimizer = "", grad_function_passed = NA, n_optimizer_calls = NA_integer_,
    n_numeric_gradients = NA_integer_, ll_calls_per_gradient = "", matches_2K_plus_1 = NA,
    ll_calls_total = NA_integer_, ll_calls_outside_gradients = NA_integer_,
    message_printed = NA, numgrad_fallback_message = NA, LL = "", seconds = NA_real_,
    setting = nm, default_value = paste(deparse(def[[nm]]), collapse = " "),
    note = sprintf("set in the function body (list 'default'); formal defaults: %s",
                   paste(sprintf("%s = %s", names(fm), fm), collapse = ", ")),
    stringsAsFactors = FALSE)

# ---- table -------------------------------------------------------------------------
tab <- do.call(rbind, ROWS)
same_end <- function(a, b) {
  x <- LLS[[a]]; y <- LLS[[b]]
  if (length(x) != length(y) || anyNA(c(x, y))) return("not comparable")
  sprintf("BFGS end LLs equal %s's to %.0e: %s", sub("_.*", "", b), 0.01, all(abs(x - y) < 0.01))
}
for (p in list(c("M2_searchStart_raw_analyticGrad_TRUE", "M1_searchStart_raw"),
               c("M4_searchStart_function_utilities", "M1_searchStart_raw"),
               c("M5_searchStart_preprocessed", "M1_searchStart_raw"),
               c("P2_searchStart_preprocessed", "P1_searchStart_raw"),
               c("L2a_searchStart_preprocessed", "L1_searchStart_raw"),
               c("L2b_searchStart_preprocessed_global_beta_moved", "L1_searchStart_raw"))) {
  i <- which(tab$test == p[1])
  tab$note[i] <- paste0(tab$note[i], "; ", same_end(p[1], p[2]))
}
for (v in names(VERSIONS)) tab[[v]] <- VERSIONS[[v]]
dir.create("output", showWarnings = FALSE)
write.csv(tab, OUT_CSV, row.names = FALSE)

# ---- summary -------------------------------------------------------------------------
tests <- tab[tab$test != "default_setting", ]
print(tests[, c("test", "function_called", "K", "optimizer", "grad_function_passed",
                "n_numeric_gradients", "ll_calls_per_gradient", "matches_2K_plus_1",
                "ll_calls_outside_gradients", "message_printed", "LL", "seconds")], row.names = FALSE)
cat("\napollo_searchStart defaults (Apollo", VERSIONS[["apollo"]], "):",
    paste(sprintf("%s = %s", tab$setting[tab$test == "default_setting"],
                  tab$default_value[tab$test == "default_setting"]), collapse = "; "), "\n")
raw_ss <- tests[tests$function_called == "apollo_searchStart" &
                  grepl("raw", tests$test) & !grepl("function_utilities", tests$test), ]
cat(sprintf("\nSearches on models as written (%s): grad passed = %s, every numerical gradient 2K+1 evaluations = %s, message printed = %s\n",
            paste(sub("_.*", "", raw_ss$test), collapse = ", "),
            paste(unique(raw_ss$grad_function_passed), collapse = "/"),
            all(raw_ss$matches_2K_plus_1), all(raw_ss$message_printed)))
stamp("DONE: wrote %s (%d rows)", OUT_CSV, nrow(tab))
