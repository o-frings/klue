#!/usr/bin/env Rscript
# =============================================================================
# dev/run_benchmark_apollo_defaults.R
#
# Why. version12.tex quotes one well-separated benchmark (randomised-ladder
# paragraph): klue's workflow, Apollo's apollo_searchStart + apollo_estimate and
# gmnl all reach LL -3021.2, klue in 1.5 s, Apollo in 173 s. That Apollo arm ran
# a reduced, non-default search (30 candidates, smartStart = TRUE, 3 stages,
# gTest = 10) from a hand-set generic start, kept in dev/compare_packages.R as
# run_apollo_ss_published_randomised(). This script times the same benchmark
# with Apollo's search as shipped, and with that search given analytic
# gradients. Apollo 0.3.5's apollo_searchStart passes the raw
# apollo_probabilities to apollo_makeLogLike and never calls
# apollo_modifyUserDefFunc, so maxLik's BFGS runs on numerical gradients;
# apollo_estimate pre-processes the same model and gets analytic ones. On a
# ladder cell (2026-10-05) pre-processing made each candidate-stage 8.4x faster,
# with the same candidate LLs (to 0.002).
#
# Data. Taken to be the validation gate in main() of dev/compare_packages.R:
#   klue_simulate(N_per_class = 100, T_tasks = 12, true_K = 3, separation = 1.5,
#                 heterogeneity = 0.2, seed = 42,
#                 dgp = klue_dgp(n_generic = 4, n_alternatives = 3))
# randomised attributes (no design), fitted at C = 3: 300 respondents, 3,600
# choices, 23 free parameters. NOT VERIFIED: no log or result file of the May
# 2026 run survives, and the harness entered git only in September (db3687b).
# The evidence: the gate is the repository's only one-dataset comparison of
# klue, Apollo and gmnl; output/compare_packages_stress.rds, dated 2026-05-30
# (the day before commit 98b00dc added the sentence), carries the same five
# method labels, so the same harness was in use; and klue_simulate reproduces
# the 0.6.x data bit-exactly (klue NEWS, 0.9.0). The run checks the
# identification: klue_ml's LL must round to -3021.2 (data_check), and
# apollo_ss_published, a deterministic rerun of the published arm, should match
# too (matches_published). If data_check is FALSE, this is not the published
# dataset.
#
# Arms, all on that dataset at C = 3:
#   klue_ml                      klue_lcmnl(estimator = "ml"), klue 0.10.0
#                                defaults: six clustering starts, direct ML
#   apollo_searchStart           run_apollo_searchStart(db, C, dgp): Apollo 0.3.5
#                                defaults as shipped (100 candidates within +-0.1
#                                of the start, smartStart FALSE, 5 stages, gTest
#                                1e-3, bfgsIter 20) from "mnl_scaled", then
#                                apollo_estimate (BGW)
#   apollo_searchStart_analytic  the same after apollo_modifyUserDefFunc
#                                pre-processing: same start, same candidates,
#                                same final apollo_estimate
#   apollo_ss_published          run_apollo_ss_published_randomised(), the
#                                published configuration, for reference
# Per arm: LL, seconds (set-up, search and estimate, as the harness times them),
# the gradient the search worked with (search_gradient: the two calls
# apollo_searchStart makes, repeated untimed after the arms) and the 1-minute
# load average at the start. Seconds depend on load; the published 1.5 s and
# 173 s come from the May run, load unknown. Apollo's changelog for 0.3.8 says
# apollo_searchStart now uses BGW, so these 0.3.5 timings need not carry over to
# current Apollo.
#
# Seeds. Data: klue_simulate seed 42. klue seeds its clustering starts
# internally. apollo_searchStart draws its candidates after set.seed(17)
# (apollo_control seed 13, Apollo's default, + 4); BGW is deterministic. LLs
# reproduce exactly; seconds do not.
#
# Run from the repository root, on one core (it counts toward the two-core cap):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 15 \
#     Rscript dev/run_benchmark_apollo_defaults.R >> output/benchmark_apollo_defaults.log 2>&1
# About 20-40 min (estimated, not measured), most of it the apollo_searchStart
# arm. Resumable: each arm is saved when it ends and skipped on a rerun with the
# same Apollo and klue versions and settings; delete the .rds to start over.
# BENCH_TEST=1: tiny searches (4 candidates, 1 stage, 5 BFGS iterations; the klue
# arm and the final apollo_estimate unchanged), about 1-3 min, nothing resumed,
# output only to BENCH_TEST_DIR (default: a folder in R's tempdir()). Point
# TMPDIR there as well: Apollo writes its search files to tempdir().
#
# Writes: output/benchmark_apollo_defaults.rds (per arm, for resuming) and
# output/benchmark_apollo_defaults.csv (one row per arm); Apollo writes
# <model>_searchStart.csv files to R's tempdir(), removed when R exits.
# =============================================================================

if (!file.exists("dev/compare_packages.R")) stop("run from the repository root")
TEST     <- identical(Sys.getenv("BENCH_TEST"), "1")
TEST_DIR <- Sys.getenv("BENCH_TEST_DIR")
if (!nzchar(TEST_DIR)) TEST_DIR <- file.path(tempdir(), "benchmark_apollo_defaults_test")
OUT_DIR  <- if (TEST) TEST_DIR else "output"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
OUT_RDS  <- file.path(OUT_DIR, "benchmark_apollo_defaults.rds")
OUT_CSV  <- file.path(OUT_DIR, "benchmark_apollo_defaults.csv")
PUBLISHED_LL <- -3021.2          # version12.tex: the LL every method reached

stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     sprintf(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a)) b else a

sys.source("dev/compare_packages.R", envir = globalenv())   # functions only (guarded)
options(mc.cores = 1L)
stopifnot(exists(".apollo_lc_setup"), exists(".apollo_est"), exists(".apollo_result"),
          exists("run_apollo_ss_published_randomised"), exists("SS_PUBLISHED_RANDOMISED"))
VERSIONS <- c(apollo = as.character(packageVersion("apollo")),
              klue   = as.character(packageVersion("klue")))

# ---- the benchmark dataset (dev/compare_packages.R, main()) ------------------
DGP     <- klue_dgp(n_generic = 4, n_alternatives = 3)
BENCH_C <- 3
DB <- klue_simulate(N_per_class = 100, T_tasks = 12, true_K = BENCH_C, separation = 1.5,
                    heterogeneity = 0.2, seed = 42, dgp = DGP)$database
stamp("data: %d respondents, %d choices, choice shares %s", length(unique(DB$ID)),
      nrow(DB), paste(sprintf("%.3f", as.numeric(prop.table(table(DB$CHOICE)))),
                      collapse = " / "))

# ---- search settings ---------------------------------------------------------
# NULL passes no searchStart_settings, i.e. Apollo's defaults as shipped.
SS_TINY      <- list(nCandidates = 4, maxStages = 1, bfgsIter = 5)
SS_DEFAULT   <- if (TEST) SS_TINY else NULL
SS_PUBLISHED <- if (TEST) modifyList(SS_PUBLISHED_RANDOMISED, SS_TINY) else SS_PUBLISHED_RANDOMISED
fmt_settings <- function(s) {
  if (is.null(s)) return(paste("Apollo 0.3.5 defaults: nCandidates 100, start +-0.1,",
                               "smartStart FALSE, maxStages 5, dTest 1, gTest 1e-3,",
                               "llTest 3, bfgsIter 20"))
  paste0(paste(names(s), vapply(s, function(v) paste(format(v), collapse = "/"), character(1)),
               collapse = ", "), "; others at Apollo defaults")
}

# ---- Apollo search on analytic gradients -------------------------------------
# apollo_estimate's pre-processing (apollo_modifyUserDefFunc: class loop
# expanded, utilities and lcPars wrapped as functions), applied to the globals
# that .apollo_lc_setup sets, before apollo_searchStart sees the model. The
# model is unchanged (same LL at the start). The global apollo_beta must be the
# start: apollo_expandLoop checks the expansion at the global value and keeps
# the unexpanded loop, whose class-specific gradients are zero, when the two
# differ. Scaling stays at 1. Returns the pre-processed apollo_probabilities and
# a copy of apollo_inputs carrying the pre-processed lcPars and scaling, or NULL.
.preprocess_for_search <- function() {
  ai <- apollo_inputs
  ai$silent <- FALSE   # read by apollo_modifyUserDefFunc (apollo_estimate sets it too)
  L <- tryCatch(apollo_modifyUserDefFunc(apollo_beta, apollo_fixed, apollo_probabilities, ai,
                                         validate = TRUE, noModification = FALSE),
                error = function(e) { message("apollo_modifyUserDefFunc error: ",
                                              conditionMessage(e)); NULL })
  if (is.null(L) || !isTRUE(L$success)) return(NULL)
  ai$apollo_lcPars  <- L$apollo_lcPars
  ai$apollo_scaling <- L$apollo_scaling
  ai$manualScaling  <- L$manualScaling
  list(probabilities = L$apollo_probabilities, inputs = ai)
}

# run_apollo_searchStart() with the pre-processing step added: same set-up, same
# timing, and the final apollo_estimate on the raw model (.apollo_est), which
# pre-processes it itself.
run_apollo_searchStart_analytic <- function(db, C, dgp, settings = NULL, start = "mnl_scaled",
                                            order = "parameter", ref = "last") {
  t0 <- Sys.time()
  if (!.apollo_lc_setup(db, C, dgp, "bench_search_an", start = start, order = order, ref = ref))
    return(.apollo_result(t0 = t0, note = "setup failed"))
  pp <- .preprocess_for_search()
  if (is.null(pp)) return(.apollo_result(t0 = t0, note = "pre-processing failed"))
  start_beta <- tryCatch(
    if (is.null(settings))
      apollo_searchStart(apollo_beta, apollo_fixed, pp$probabilities, pp$inputs)
    else
      apollo_searchStart(apollo_beta, apollo_fixed, pp$probabilities, pp$inputs,
                         searchStart_settings = settings),
    error = function(e) { message("apollo_searchStart error: ", conditionMessage(e)); NULL })
  m <- if (!is.null(start_beta)) .apollo_est(start_beta)
  .apollo_result(m, t0, if (is.null(start_beta)) "searchStart failed"
                        else if (is.null(m)) "estimate errored" else "")
}

# The gradient a search works with, from the calls apollo_searchStart makes:
# apollo_makeLogLike on the function it is given, then (with analyticGrad, the
# default) apollo_makeGrad(validateGrad = TRUE) at the start. NULL means maxLik's
# numerical gradient. Run untimed, after the arms.
search_gradient <- function(db, C, dgp, preprocess, start, order, ref) {
  on.exit(cleanup_apollo(), add = TRUE)
  if (!.apollo_lc_setup(db, C, dgp, "bench_gradcheck", start = start, order = order, ref = ref))
    return("setup failed")
  fn <- apollo_probabilities; ai <- apollo_inputs
  if (preprocess) {
    pp <- .preprocess_for_search()
    if (is.null(pp)) return("pre-processing failed")
    fn <- pp$probabilities; ai <- pp$inputs
  }
  ai$apollo_control$noDiagnostics <- TRUE   # as apollo_searchStart sets it
  tryCatch({
    ll <- apollo_makeLogLike(apollo_beta, apollo_fixed, fn, ai, list(estimationRoutine = "BHHH"))
    g  <- if (isTRUE(ai$apollo_control$analyticGrad))
            apollo_makeGrad(apollo_beta, apollo_fixed, ll, validateGrad = TRUE)
    if (is.null(g)) "numeric" else "analytic"
  }, error = function(e) { message("gradient check error: ", conditionMessage(e)); "check failed" })
}

# 1-minute load average (Linux /proc/loadavg, macOS sysctl), NA if unavailable.
load1 <- function() {
  x <- tryCatch(suppressWarnings(
         if (file.exists("/proc/loadavg")) readLines("/proc/loadavg", n = 1L)
         else system2("sysctl", c("-n", "vm.loadavg"), stdout = TRUE, stderr = FALSE)),
       error = function(e) "")
  v <- suppressWarnings(as.numeric(strsplit(gsub("[{}]", " ", paste(x, collapse = " ")),
                                            "[[:space:]]+")[[1]]))
  v <- v[is.finite(v)]
  if (length(v)) v[1] else NA_real_
}

# ---- arms --------------------------------------------------------------------
START_EX  <- "mnl_scaled (pooled MNL, class c at 1/c), parameter-major, last class reference"
START_PUB <- "generic (constants 0, b_x 0.1(c-1), b_price -0.5), class-major, first class reference"
GRAD_EX   <- list(start = "mnl_scaled", order = "parameter", ref = "last")
ARMS <- list(
  klue_ml = list(
    what = "klue_lcmnl, estimator ml: six clustering starts, direct ML",
    settings = "klue defaults", start = "clustering starts",
    published_LL = PUBLISHED_LL, published_seconds = 1.5, grad = NULL,
    run = function() run_klue(DB, BENCH_C, DGP, "ml")),
  apollo_searchStart = list(
    what = "apollo_searchStart as shipped (numerical gradients), then apollo_estimate",
    settings = fmt_settings(SS_DEFAULT), start = START_EX,
    published_LL = NA_real_, published_seconds = NA_real_,
    grad = c(list(preprocess = FALSE), GRAD_EX),
    run = function() run_apollo_searchStart(DB, BENCH_C, DGP, settings = SS_DEFAULT)),
  apollo_searchStart_analytic = list(
    what = paste("apollo_searchStart after apollo_modifyUserDefFunc pre-processing",
                 "(analytic gradients), then apollo_estimate"),
    settings = fmt_settings(SS_DEFAULT), start = START_EX,
    published_LL = NA_real_, published_seconds = NA_real_,
    grad = c(list(preprocess = TRUE), GRAD_EX),
    run = function() run_apollo_searchStart_analytic(DB, BENCH_C, DGP, settings = SS_DEFAULT)),
  apollo_ss_published = list(
    what = "reduced apollo_searchStart of the published 173 s arm, then apollo_estimate",
    settings = fmt_settings(SS_PUBLISHED), start = START_PUB,
    published_LL = PUBLISHED_LL, published_seconds = 173,
    grad = list(preprocess = FALSE, start = "generic", order = "class", ref = "first"),
    run = function()
      if (TEST) run_apollo_searchStart(DB, BENCH_C, DGP, settings = SS_PUBLISHED,
                                       start = "generic", order = "class", ref = "first")
      else run_apollo_ss_published_randomised(DB, BENCH_C, DGP))
)
RUN_ORDER <- c("klue_ml", "apollo_ss_published", "apollo_searchStart_analytic",
               "apollo_searchStart")   # cheap arms first
stopifnot(setequal(RUN_ORDER, names(ARMS)))

# ---- run (resumable per arm) -------------------------------------------------
results <- if (!TEST && file.exists(OUT_RDS)) readRDS(OUT_RDS) else list()
stamp("%s run: Apollo %s, klue %s -> %s", if (TEST) "TEST" else "full",
      VERSIONS[["apollo"]], VERSIONS[["klue"]], OUT_CSV)
for (arm in RUN_ORDER) {
  a <- ARMS[[arm]]; old <- results[[arm]]
  if (!is.null(old) && identical(old$versions, VERSIONS) && identical(old$settings, a$settings)) {
    stamp("SKIP %s (done %s: LL %.3f, %.1f s)", arm, old$started, old$LL, old$seconds)
    next
  }
  stamp("RUN  %s  [%s]", arm, a$settings)
  started <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ld  <- load1()
  fit <- a$run()
  results[[arm]] <- list(LL = as.numeric(fit$LL), converged = isTRUE(fit$converged),
                         seconds = as.numeric(fit$seconds), note = fit$note %||% "",
                         started = started, load1 = ld, versions = VERSIONS,
                         settings = a$settings)
  saveRDS(results, OUT_RDS)
  stamp("  done %s  LL=%.3f  %.1f s  converged=%s%s", arm, as.numeric(fit$LL),
        as.numeric(fit$seconds), isTRUE(fit$converged),
        if (nzchar(fit$note %||% "")) paste0("  note: ", fit$note) else "")
}

# ---- gradient checks (untimed) -----------------------------------------------
stamp("gradient checks: the gradient each Apollo search worked with")
grad_kind <- vapply(names(ARMS), function(arm) {
  g <- ARMS[[arm]]$grad
  if (is.null(g)) NA_character_
  else search_gradient(DB, BENCH_C, DGP, g$preprocess, g$start, g$order, g$ref)
}, character(1))

# ---- table -------------------------------------------------------------------
tab <- do.call(rbind, lapply(names(ARMS), function(arm) {
  a <- ARMS[[arm]]; r <- results[[arm]]
  data.frame(arm = arm, LL = r$LL, converged = r$converged, seconds = round(r$seconds, 1),
             search_gradient = grad_kind[[arm]], published_LL = a$published_LL,
             published_seconds = a$published_seconds, settings = a$settings,
             start = a$start, what = a$what, note = r$note, started = r$started,
             load1 = r$load1, apollo = r$versions[["apollo"]], klue = r$versions[["klue"]],
             mode = if (TEST) "test" else "full", stringsAsFactors = FALSE)
}))
fin  <- is.finite(tab$LL)
best <- if (any(fin)) max(tab$LL[fin]) else NA_real_
tab$gap_to_best <- round(tab$LL - best, 3)
tab$miss        <- tab$gap_to_best < -0.5       # the paper's miss: > 0.5 LL below the best
tab$times_klue  <- round(tab$seconds / tab$seconds[tab$arm == "klue_ml"], 1)
tab$matches_published <- ifelse(is.na(tab$published_LL), NA,
                                abs(tab$LL - tab$published_LL) <= 0.05 + 1e-9)
data_ok <- isTRUE(tab$matches_published[tab$arm == "klue_ml"])
tab$data_check <- data_ok
tab <- tab[, c("arm", "LL", "gap_to_best", "miss", "converged", "seconds", "times_klue",
               "search_gradient", "published_LL", "published_seconds", "matches_published",
               "data_check", "settings", "start", "what", "note", "started", "load1",
               "apollo", "klue", "mode")]
write.csv(tab, OUT_CSV, row.names = FALSE)

klue_LL <- tab$LL[tab$arm == "klue_ml"]
if (data_ok) {
  stamp("data check passed: klue_ml LL %.3f rounds to the published %.1f", klue_LL, PUBLISHED_LL)
} else {
  stamp("DATA CHECK FAILED: klue_ml LL %.3f, published %.1f; this is not the published dataset",
        klue_LL, PUBLISHED_LL)
}
if (!identical(grad_kind[["apollo_searchStart_analytic"]], "analytic"))
  stamp("WARNING: the pre-processed search did not get analytic gradients (%s)",
        grad_kind[["apollo_searchStart_analytic"]])
print(tab[, c("arm", "LL", "gap_to_best", "converged", "seconds", "times_klue",
              "search_gradient", "published_LL", "published_seconds", "matches_published",
              "load1")], row.names = FALSE)
stamp("DONE: wrote %s", OUT_CSV)
