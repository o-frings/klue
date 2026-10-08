#!/usr/bin/env Rscript
# =============================================================================
# dev/run_apollo_isolation_blocked.R
#
# Why. On the blocked stress ladder (dev/stress_replicate_blocked.R) Apollo's
# apollo_searchStart at its 0.3.5 defaults (100 candidates within +-0.1 of the
# start, 5 stages, gTest 1e-3) starts from a vector built by this project's
# harness: klue's pooled MNL with class c scaled by 1/c ("mnl_scaled"). The arm
# reported in version12 used reduced settings, a [-3, 3] window and a generic
# start. The two arms differ in start, search region, budget, pruning and
# parameter layout at once; arms A0-A4 separate those effects (audit of
# 2026-10-03). A5 and A6 (2026-10-05) separate two sources of Apollo's cost.
# Apollo 0.3.5's apollo_searchStart passes apollo_probabilities to
# apollo_makeLogLike without the pre-processing that apollo_estimate applies
# (apollo_modifyUserDefFunc), so the search's BFGS runs on numerical gradients
# (Apollo logs "not able to compute analytical gradients"). On hard|seed1 one
# candidate-stage (20 BFGS iterations) took 10.5 s that way and 1.25 s with
# the model pre-processed; the same candidates reached the same LLs (within
# 0.002). And the ladder's apollo_clust computes a covariance matrix after each
# of its six fits (on hard|seed1 about 11.5 s of the 13-15 s of a converged
# fit), where klue_ml computes one, at the best fit.
#
# Arms (parameter-major layout, last class as the allocation reference, as in the
# default arm):
#   A0  apollo_estimate (BGW) from mnl_scaled, no search
#   A4  apollo_estimate from the generic start (constants 0, b_x = 0.1(c-1),
#       b_price = -0.5), no search: how Apollo's course scripts run a latent
#       class model (hand-set start, no search)
#   A1  apollo_searchStart, reduced budget (30 candidates, 3 stages, gTest 10,
#       bfgsIter 20, smartStart FALSE) in the default +-0.1 box around mnl_scaled
#   A3  the reported blocked settings (same budget, [-3, 3] window) from mnl_scaled
#   A2  apollo_searchStart at Apollo's defaults from the generic start
#   A5  the default arm (apollo_searchStart at Apollo's defaults from
#       mnl_scaled, then apollo_estimate from its best candidate) with the model
#       pre-processed by apollo_modifyUserDefFunc before the search, so that
#       the search uses analytic gradients; same candidates, same final
#       apollo_estimate. Also records the search's stages and candidate-stages
#   A6  the ladder's apollo_clust (apollo_estimate from klue's six clustering
#       starts, best LL kept) with hessianRoutine = "none" in the six fits, then
#       one apollo_estimate at Apollo's defaults (analytic Hessian) from the
#       best fit: one covariance matrix, as in klue_ml. The clustering starts,
#       the six fits and the covariance are timed separately
#   chk klue's clustering + direct ML (as klue_ml on the ladder): its LL must
#       equal the ladder's klue_ml LL, which proves the data are the same
#
# Reading, fixed before running. "Default" = the ladder's apollo_searchStart arm
# (defaults, mnl_scaled start). A miss = more than 0.5 LL below the best LL found
# in that cell by any arm here or on the ladder.
#   - A1 hits where the default hits: the budget is not what matters.
#   - A2 misses where the default hits: the start does the work. Apollo's search
#     succeeds from an MNL-informed start; the clustering rule supplies starts
#     without such design choices.
#   - A2 hits wherever the default hits: Apollo's default search works from a
#     hand-set start too; the contrast with Apollo reduces to cost and to what
#     cheaper settings do.
#   - A3 misses where the reported arm missed: the reported misses come from the
#     wide window (or budget), not from the generic start.
#   - A0 and A4 show what one fit from each start does without any search.
#   - A5 misses only where the default arm missed (hard|seed5 among these cells)
#     and takes about 9-10 min per fit; any other hit pattern means the gradient
#     source changes the search, and both versions must be reported.
#   - A6 has no reading rule (cost only).
#
# Cells, fixed in advance: seeds 1-5 of the hard and very_hard rungs (10 cells;
# the reported arm missed hard 3, 4, 5 and very_hard 2, 3, 5 among them). The
# data are regenerated bit-exactly with the ladder's seed formula
# (1000 K + 100 (100 kappa) + seed) and its shared blocked design (klue_design,
# seeded internally). apollo_searchStart draws its candidates from its own
# fixed seed (17 under Apollo's default control seed), so A5 starts from the
# default arm's candidates; klue's clustering starts are deterministic.
#
# Usage (repository root; one core per process, Apollo nCores = 1; prefix each
# run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
# nice -n 15):
#   ARMS=chk,A0,A4 Rscript dev/run_apollo_isolation_blocked.R   # cheap arms
#   ARMS=A1,A3     Rscript dev/run_apollo_isolation_blocked.R   # about 15-30 min/cell
#   ARMS=A2        Rscript dev/run_apollo_isolation_blocked.R   # about 75 min/cell
#   ARMS=A6,A5     Rscript dev/run_apollo_isolation_blocked.R   # about 0.5 + 10 min/cell
#   Rscript dev/run_apollo_isolation_blocked.R summary           # merge and tabulate
# Each process first rebuilds the shared design (about 2.5 min).
# Writes. Each process writes output/apollo_isolation_blocked_<ARMS>.rds, saved
# after every (cell, arm) and resumable; "summary" merges all of them with the
# ladder's best LLs into output/apollo_isolation_blocked.csv (one row per cell
# and arm: LL, gap to the best, miss, seconds, note, then A6's starts_secs,
# fit_secs, cov_secs and LL_fits and A5's search_secs, stages and cand_stages,
# NA for the other arms) and prints the reading tables. The ladder's rds is
# only read (copied first, since its own job rewrites it after each cell).
#
# Test (about 4 min, most of it the design): ISO_TEST=1 runs the requested arms
# on hard|seed1 only, every search at nCandidates = 2, maxStages = 1,
# bfgsIter = 5 in its arm's box or window (A6 and the arms without a search run
# in full), and writes its rds and summary CSV under ISO_SMOKE (default
# apollo_isolation_test in the system temp directory) instead of output/.
# Apollo writes its search file to R's tempdir(), so set TMPDIR inside the
# smoke directory too:
#   mkdir -p "$SMOKE/rtmp"
#   ISO_TEST=1 ISO_SMOKE="$SMOKE" TMPDIR="$SMOKE/rtmp" ARMS=chk,A5,A6 \
#     OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/run_apollo_isolation_blocked.R
#   ISO_TEST=1 ISO_SMOKE="$SMOKE" TMPDIR="$SMOKE/rtmp" \
#     Rscript dev/run_apollo_isolation_blocked.R summary
# Expected: chk LL -3168.965 (the ladder's klue_ml), A5 analytic_grad=TRUE with
# stages=1 and cand_stages=2, A6's best of six equal to the ladder's apollo_clust.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
SUMMARY <- identical(args, "summary")
TEST <- identical(Sys.getenv("ISO_TEST"), "1")
ARMS_ALL <- c("chk", "A0", "A4", "A1", "A3", "A2", "A5", "A6")
LADDER <- "output/stress_replicate_blocked.rds"
RUNGS <- list(hard      = list(K = 5, kap = 0.50, sig = 0.30),
              very_hard = list(K = 5, kap = 0.30, sig = 0.35))
SEEDS <- 1:5
OUT_DIR <- "output"
SS_TEST <- list(nCandidates = 2, maxStages = 1, bfgsIter = 5)   # TEST: the budget of every search
stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     sprintf(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a)) b else a
read_ladder <- function() {   # copy first: the ladder job rewrites the file after each cell
  for (i in 1:5) {
    tmp <- tempfile(fileext = ".rds")
    ok <- file.copy(LADDER, tmp, overwrite = TRUE)
    x <- if (ok) tryCatch(readRDS(tmp), error = function(e) NULL)
    unlink(tmp)
    if (!is.null(x)) return(x)
    Sys.sleep(5)
  }
  stop("could not read ", LADDER)
}

if (TEST) {
  SMOKE <- Sys.getenv("ISO_SMOKE")
  if (!nzchar(SMOKE)) SMOKE <- file.path(dirname(tempdir()), "apollo_isolation_test")
  dir.create(SMOKE, recursive = TRUE, showWarnings = FALSE)
  RUNGS <- RUNGS["hard"]; SEEDS <- 1L; OUT_DIR <- SMOKE
  stamp("TEST: hard|seed1 only, searches at %s; files under %s; R tempdir %s",
        paste(names(SS_TEST), unlist(SS_TEST), sep = " = ", collapse = ", "), SMOKE, tempdir())
}

if (SUMMARY) {
  files <- list.files(OUT_DIR, pattern = "^apollo_isolation_blocked_.*\\.rds$", full.names = TRUE)
  iso <- list()
  for (f in files) for (key in names(r <- readRDS(f))) iso[[key]] <- c(iso[[key]], r[[key]])
  if (!length(iso)) { stamp("no results under %s", OUT_DIR); quit(save = "no", status = 0) }
  lad <- read_ladder()
  num <- function(v) if (is.null(v)) NA_real_ else as.numeric(v)
  rows <- list()
  for (key in names(iso)) {
    lc <- lad[[key]]
    best <- max(c(lc$best, vapply(iso[[key]], function(a) a$LL %||% NA_real_, numeric(1))), na.rm = TRUE)
    pick <- function(v, m) if (m %in% names(v)) unname(v[[m]]) else NA_real_
    # x: an arm's own fields (A6: time split and best of six; A5: the search)
    add <- function(arm, LL, secs, note = "", x = list())
      rows[[length(rows) + 1]] <<- data.frame(cell = key, arm = arm, LL = round(LL, 3),
                                              gap = round(LL - best, 3), miss = LL - best < -0.5,
                                              secs = round(secs, 1), note = note,
                                              starts_secs = round(num(x$starts_secs), 1),
                                              fit_secs = round(num(x$fit_secs), 1),
                                              cov_secs = round(num(x$cov_secs), 1),
                                              LL_fits = round(num(x$LL_fits), 3),
                                              search_secs = round(num(x$search_secs), 1),
                                              stages = num(x$stages), cand_stages = num(x$cand_stages))
    for (m in c("klue_ml", "apollo_searchStart", "apollo_clust", "apollo_ss_published"))
      add(paste0("ladder:", m), pick(lc$LL, m), pick(lc$secs, m))
    for (a in intersect(ARMS_ALL, names(iso[[key]]))) {
      z <- iso[[key]][[a]]
      add(a, z$LL %||% NA_real_, z$seconds %||% NA_real_, z$note %||% "", z)
    }
  }
  tab <- do.call(rbind, rows)
  write.csv(tab, file.path(OUT_DIR, "apollo_isolation_blocked.csv"), row.names = FALSE)
  chk <- tab[tab$arm == "chk", ]
  lml <- tab[tab$arm == "ladder:klue_ml", ]
  m <- merge(chk[, c("cell", "LL")], lml[, c("cell", "LL")], by = "cell", suffixes = c("_chk", "_ladder"))
  if (nrow(m)) stamp("data check: chk LL equals ladder klue_ml LL in %d of %d cells (|diff| < 1e-3)",
                     sum(abs(m$LL_chk - m$LL_ladder) < 1e-3), nrow(m))
  per_arm <- do.call(rbind, lapply(split(tab, tab$arm), function(d)
    data.frame(arm = d$arm[1], cells = sum(!is.na(d$LL)), misses = sum(d$miss %in% TRUE),
               missed_cells = paste(d$cell[d$miss %in% TRUE], collapse = " "),
               median_secs = median(d$secs, na.rm = TRUE))))
  print(per_arm[order(match(sub("^ladder:", "", per_arm$arm),
                            c("klue_ml", "apollo_searchStart", "apollo_clust",
                              "apollo_ss_published", ARMS_ALL))), ], row.names = FALSE)
  col_at <- function(arm, col, cells) { d <- tab[tab$arm == arm, ]; d[[col]][match(cells, d$cell)] }
  # A5 against the default arm, cell by cell (the reading in the header)
  a5 <- tab[tab$arm == "A5", ]
  if (nrow(a5)) {
    cat("\nA5 (default arm, search on analytic gradients) against the default arm:\n")
    print(data.frame(cell = a5$cell, gap_A5 = a5$gap,
                     gap_default = col_at("ladder:apollo_searchStart", "gap", a5$cell),
                     min_A5 = round(a5$secs / 60, 1),
                     min_default = round(col_at("ladder:apollo_searchStart", "secs", a5$cell) / 60, 1),
                     cand_stages = a5$cand_stages,
                     s_per_cand_stage = round(a5$search_secs / a5$cand_stages, 2),
                     note = a5$note), row.names = FALSE)
    run5 <- !is.na(a5$LL)
    miss5 <- a5$cell[run5 & a5$miss %in% TRUE]
    missd <- a5$cell[run5 & col_at("ladder:apollo_searchStart", "miss", a5$cell) %in% TRUE]
    stamp("A5 misses: %s; default misses (same %d cells): %s -> %s; A5 median %.1f min per cell",
          if (length(miss5)) paste(miss5, collapse = " ") else "none", sum(run5),
          if (length(missd)) paste(missd, collapse = " ") else "none",
          if (setequal(miss5, missd)) "same cells, as the reading expects"
          else "different cells: the gradient source changes the search, report both versions",
          median(a5$secs[run5]) / 60)
  }
  # A6: where Apollo's time goes once the covariance is computed once
  a6 <- tab[tab$arm == "A6", ]
  if (nrow(a6)) {
    sp <- data.frame(cell = a6$cell, klue_ml = col_at("ladder:klue_ml", "secs", a6$cell),
                     apollo_clust = col_at("ladder:apollo_clust", "secs", a6$cell),
                     A6_starts = a6$starts_secs, A6_fits = a6$fit_secs, A6_cov = a6$cov_secs,
                     A6_total = a6$secs)
    cat("\nA6 (apollo_clust, one covariance at the best fit), seconds:\n")
    print(sp, row.names = FALSE)
    md <- vapply(sp[-1], median, numeric(1), na.rm = TRUE)
    stamp("A6 medians: starts %.1f s, six fits %.1f s, covariance %.1f s, total %.1f s (ladder: apollo_clust %.1f s, klue_ml %.1f s)",
          md[["A6_starts"]], md[["A6_fits"]], md[["A6_cov"]], md[["A6_total"]],
          md[["apollo_clust"]], md[["klue_ml"]])
    stamp("A6 best of six LL equals the ladder's apollo_clust LL in %d of %d cells (|diff| < 1e-3)",
          sum(abs(a6$LL_fits - col_at("ladder:apollo_clust", "LL", a6$cell)) < 1e-3, na.rm = TRUE),
          nrow(a6))
  }
  quit(save = "no", status = 0)
}

arms <- strsplit(Sys.getenv("ARMS", "chk,A0,A4"), ",")[[1]]
stopifnot(length(arms) > 0, all(arms %in% ARMS_ALL))
OUT <- file.path(OUT_DIR, sprintf("apollo_isolation_blocked_%s.rds", paste(arms, collapse = "")))

sys.source("dev/compare_packages.R", envir = globalenv())   # functions only (guarded)
options(mc.cores = 1L)
DGP    <- klue_dgp(n_generic = 4, n_alternatives = 3)
DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)   # the ladder's shared design
N_PER_CLASS <- 60L

run_apollo_estimate_only <- function(db, C, dgp, start) {
  t0 <- Sys.time()
  if (!.apollo_lc_setup(db, C, dgp, "iso_est", start = start, order = "parameter", ref = "last"))
    return(.apollo_result(t0 = t0, note = "setup failed"))
  m <- .apollo_est(apollo_beta, "apollo_estimate (no search)")
  .apollo_result(m, t0, if (is.null(m)) "estimate errored" else "")
}
SS_REDUCED_BOX <- list(nCandidates = 30, smartStart = FALSE, dTest = 1, gTest = 10,
                       maxStages = 3, bfgsIter = 20)
ss <- function(s) if (TEST) modifyList(s %||% list(), SS_TEST) else s   # TEST: tiny budget, same box or window

# Stages and candidate-stages of one search, from the file apollo_searchStart
# keeps in its output directory: LL1 (initial candidates), then one LL column
# per stage, NA where a candidate was no longer active. Candidate-stages = the
# active candidates summed over stages, the count behind the ladder's timings
# ("Stage s, n active candidates" in its log).
ss_stages <- function(f) {
  x <- if (file.exists(f)) tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(x)) return(list(stages = NA_integer_, cand_stages = NA_integer_))
  ll <- grep("^LL[0-9]+$", names(x))
  list(stages = max(length(ll) - 1L, 0L),
       cand_stages = if (length(ll) > 1) sum(!is.na(as.matrix(x[, ll[-1], drop = FALSE]))) else 0L)
}

# A5. The default arm (run_apollo_searchStart with settings = NULL, start
# mnl_scaled) with the model pre-processed as apollo_estimate does before its
# optimiser (recipe of the gradient experiment of 2026-10-05), under two
# conditions:
#   - the global apollo_beta equals the start: Apollo checks its rewrites of
#     the model code at the global apollo_beta and, when it differs, keeps the
#     class loop unexpanded (see .apollo_est in dev/compare_packages.R);
#   - scaling is 1 for every parameter, so the pre-processed model has the raw
#     model's parameters and the search's best candidate goes to .apollo_est as is.
# apollo_modifyUserDefFunc reads apollo_inputs$silent (apollo_estimate sets it
# from its settings); FALSE is what the default arm's search sees. The search
# gets the pre-processed apollo_probabilities, apollo_lcPars, scaling and
# manualScaling; the final apollo_estimate (.apollo_est) gets the raw model, as
# in the default arm. Before the search, the gradient step it will take
# (apollo_makeLogLike, then apollo_makeGrad with validation) runs once on the
# same inputs: without an analytic gradient the search would fall back to
# numerical ones (the default arm again, at 75 min a cell), so the arm stops.
# seconds = the whole arm, these checks included (a few seconds).
run_apollo_searchStart_grad <- function(db, C, dgp, settings = NULL) {
  t0 <- Sys.time()
  name <- "iso_ss_grad"
  if (!.apollo_lc_setup(db, C, dgp, name, start = "mnl_scaled", order = "parameter", ref = "last"))
    return(.apollo_result(t0 = t0, note = "setup failed"))
  b0 <- apollo_beta                 # the start; the global apollo_beta stays equal to it
  ai <- apollo_inputs               # a copy: the global inputs stay raw for .apollo_est
  ai$silent <- FALSE
  L <- tryCatch(apollo_modifyUserDefFunc(b0, apollo_fixed, apollo_probabilities, ai,
                                         validate = TRUE, noModification = FALSE),
                error = function(e) { message("apollo_modifyUserDefFunc error: ",
                                              conditionMessage(e)); NULL })
  ok <- isTRUE(L$success) &&
        !any(grepl("for (s in", deparse(L$apollo_probabilities), fixed = TRUE)) &&   # loop expanded
        length(L$apollo_scaling) > 0 && isTRUE(all(L$apollo_scaling == 1))
  if (!ok) return(c(.apollo_result(t0 = t0, note = "pre-processing failed"),
                    list(analytic_grad = FALSE)))
  ai$apollo_lcPars  <- L$apollo_lcPars
  ai$apollo_scaling <- L$apollo_scaling
  ai$manualScaling  <- L$manualScaling
  ai_chk <- ai; ai_chk$apollo_control$noDiagnostics <- TRUE   # as apollo_searchStart sets it
  ll <- tryCatch(apollo_makeLogLike(b0, apollo_fixed, L$apollo_probabilities, ai_chk,
                                    list(estimationRoutine = "BHHH")), error = function(e) NULL)
  gr <- if (is.function(ll)) tryCatch(apollo_makeGrad(b0, apollo_fixed, ll, validateGrad = TRUE),
                                      error = function(e) NULL)
  analytic <- is.function(gr)
  rm(ll, gr, ai_chk)
  if (!analytic) return(c(.apollo_result(t0 = t0, note = "no analytic gradients after pre-processing"),
                          list(analytic_grad = FALSE)))
  csv <- file.path(sub("/+$", "", ai$apollo_control$outputDirectory), paste0(name, "_searchStart.csv"))
  unlink(csv)
  t1 <- Sys.time()
  start_beta <- tryCatch(
    if (is.null(settings))
      apollo_searchStart(b0, apollo_fixed, L$apollo_probabilities, ai)
    else
      apollo_searchStart(b0, apollo_fixed, L$apollo_probabilities, ai,
                         searchStart_settings = settings),
    error = function(e) { message("apollo_searchStart error: ", conditionMessage(e)); NULL })
  search_secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  st <- ss_stages(csv)
  m <- if (!is.null(start_beta)) .apollo_est(start_beta)
  c(.apollo_result(m, t0, if (is.null(start_beta)) "searchStart failed"
                          else if (is.null(m)) "estimate errored" else ""),
    list(analytic_grad = TRUE, search_secs = search_secs, stages = st$stages,
         cand_stages = st$cand_stages))
}

# A6. run_apollo_from_clustering (dev/compare_packages.R) with the covariance
# computed once: the six fits with hessianRoutine = "none", then
# apollo_estimate at Apollo's defaults (analytic Hessian) from the best of
# them, as klue_ml computes its covariance once, at the best fit. The
# covariance fit restarts BGW at that point (it stops at once at a converged
# fit) and gives the arm's LL; LL_fits is the best of the six, which should
# equal the ladder's apollo_clust LL. starts_secs = klue's clustering starts,
# fit_secs = the six fits, cov_secs = the covariance fit; seconds = the whole arm.
EST_NO_HESSIAN <- list(writeIter = FALSE, silent = TRUE, hessianRoutine = "none")
.apollo_est_with <- function(beta, settings, label) {   # .apollo_est with other settings
  apollo_beta <<- beta
  tryCatch(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                           estimate_settings = settings),
           error = function(e) { message(label, " error: ", conditionMessage(e)); NULL })
}
run_apollo_clust_one_cov <- function(db, C, dgp) {
  t0 <- Sys.time()
  starts <- tryCatch(Filter(Negate(is.null), get_all_starts(db, C, dgp = dgp)),
                     error = function(e) { message("get_all_starts error: ",
                                                   conditionMessage(e)); list() })
  starts_secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!length(starts) || !.apollo_lc_setup(db, C, dgp, "iso_clust_1cov", start = "generic"))
    return(.apollo_result(t0 = t0, note = "no starts or validateInputs failed"))
  ng <- dgp$n_generic; nb <- dgp$n_beta; J <- dgp$n_alternatives
  t1 <- Sys.time()
  fits <- lapply(seq_along(starts), function(i) {   # the start vectors of run_apollo_from_clustering
    st <- starts[[i]]; b <- apollo_beta; sh <- pmax(st$shares, 0.01)
    for (cc in 1:C) {
      for (j in 1:J) b[paste0("asc_alt", j, "_", cc)] <- 0
      for (a in 1:ng) b[paste0("b_x", a, "_", cc)] <- st$betas[cc, a]
      b[paste0("b_price_", cc)] <- st$betas[cc, nb]
      b[paste0("delta_", cc)]   <- log(sh[cc] / sh[C])
    }
    .apollo_est_with(b, EST_NO_HESSIAN, sprintf("apollo_estimate, no Hessian (start %d)", i))
  })
  fit_secs <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  lls <- vapply(fits, get_apollo_ll, numeric(1))
  if (all(!is.finite(lls))) return(.apollo_result(t0 = t0, note = "all starts errored"))
  win <- fits[[which.max(lls)]]
  b <- apollo_beta; k <- intersect(names(b), names(win$estimate)); b[k] <- win$estimate[k]
  t2 <- Sys.time()
  m <- .apollo_est(b, "apollo_estimate (covariance at the best start)")
  cov_secs <- as.numeric(difftime(Sys.time(), t2, units = "secs"))
  c(.apollo_result(if (is.null(m)) win else m, t0,
                   if (is.null(m)) "covariance fit errored" else ""),
    list(starts_secs = starts_secs, fit_secs = fit_secs, cov_secs = cov_secs,
         LL_fits = max(lls[is.finite(lls)])))
}

ARM_FUN <- list(
  chk = function(db, K) run_klue(db, K, DGP, "ml"),
  A0  = function(db, K) run_apollo_estimate_only(db, K, DGP, start = "mnl_scaled"),
  A4  = function(db, K) run_apollo_estimate_only(db, K, DGP, start = "generic"),
  A1  = function(db, K) run_apollo_searchStart(db, K, DGP, settings = ss(SS_REDUCED_BOX),
                                               start = "mnl_scaled"),
  A3  = function(db, K) run_apollo_searchStart(db, K, DGP, settings = ss(SS_PUBLISHED_BLOCKED),
                                               start = "mnl_scaled"),
  A2  = function(db, K) run_apollo_searchStart(db, K, DGP, settings = ss(NULL), start = "generic"),
  A5  = function(db, K) run_apollo_searchStart_grad(db, K, DGP, settings = ss(NULL)),
  A6  = function(db, K) run_apollo_clust_one_cov(db, K, DGP)
)

results <- if (file.exists(OUT)) readRDS(OUT) else list()
stamp("arms %s -> %s", paste(arms, collapse = ","), OUT)
for (rn in names(RUNGS)) {
  r <- RUNGS[[rn]]
  for (sd_i in SEEDS) {
    key  <- paste0(rn, "|seed", sd_i)
    todo <- setdiff(arms, names(results[[key]]))
    if (!length(todo)) next
    seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + sd_i)   # the ladder's formula
    db <- klue_simulate(N_per_class = N_PER_CLASS, true_K = r$K, separation = r$kap,
                        heterogeneity = r$sig, seed = seed_i, dgp = DGP, design = DESIGN)$database
    for (a in todo) {
      stamp("RUN  %s  %s", key, a)
      fit <- ARM_FUN[[a]](db, r$K)
      extra <- fit[setdiff(names(fit), c("LL", "converged", "seconds", "note"))]   # A5, A6 only
      results[[key]][[a]] <- c(list(LL = as.numeric(fit$LL), converged = fit$converged,
                                    seconds = as.numeric(fit$seconds), note = fit$note %||% ""),
                               extra)
      saveRDS(results, OUT)
      shown <- vapply(extra, function(v) format(if (is.numeric(v)) round(v, 3) else v), "")
      stamp("  done %s %s  LL=%.3f  (%.1f s)%s%s", key, a, as.numeric(fit$LL),
            as.numeric(fit$seconds), if (nzchar(fit$note %||% "")) paste0("  note: ", fit$note) else "",
            if (length(extra)) paste0("  ", paste(names(extra), shown, sep = "=", collapse = " ")) else "")
    }
  }
}
stamp("DONE %s", paste(arms, collapse = ","))
