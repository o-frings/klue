#!/usr/bin/env Rscript
# =============================================================================
# dev/run_mmnl_second_start.R
#
# Why. dev/rerun_mmnl_benchmark.R estimates each MMNL benchmark once, from the
# pooled MNL. On Swissmetro that start left M3 18.8 LL below a fit of the same
# model from another start (dev/run_swissmetro_m3_robustness.R). A lower MMNL
# maximum can only make the benchmark look worse, so it can only bias the
# LCMNL-versus-MMNL reading toward the LCMNL. This script therefore gives the
# benchmark a second start wherever that bias could decide the reading: where
# the best converged MMNL (M1-M4) is not at least MARGIN_BIC below the
# BIC-best LCMNL, i.e. the LCMNL wins or the MMNL wins by less than MARGIN_BIC.
# The rule is the same for all five datasets (added 2026-10-01, after the
# Swissmetro check). Until 2026-10-05 the second start went to M3 and M4 only,
# but the best MMNL can have fixed constants: on Mode, the one flagged dataset
# (margin -1.0 BIC), it is M1. Since 2026-10-05 the rule is
# symmetric, so the second start goes to the BIC-best MMNL whichever it is
# (the LCMNL side is dev/check_empirical_lcmnl_maxima.R).
#
# For each flagged dataset:
#   1. M3 (independent normals, random constants) from a neutral start
#      (start = "zero"), saved as output/mmnl_bench_std/robust/<ds>_M3_zero.rds.
#      Swissmetro's zero-start fit already exists from the Swissmetro check
#      and is reused.
#   2. If that fit reaches a higher LL than the benchmark M3 (by more than
#      MIN_GAIN), it becomes the benchmark M3 (the previous one is kept in
#      robust/), and M4 is refitted warm-started from it (the previous M4 is
#      kept in robust/), so M4 still nests the best M3.
#   3. If the best MMNL after step 2 has fixed constants (M1 or M2), steps 1
#      and 2 for M1 (independent normals, fixed constants; random = generic
#      attributes and price, as in dev/rerun_mmnl_benchmark.R) and M2: M1 from
#      start = "zero" (robust/<ds>_M1_zero.rds); if it beats the benchmark M1
#      by more than MIN_GAIN it replaces it, and M2 is refitted warm-started
#      from it; the previous fits are kept in robust/.
# After step 3 the best MMNL is always a model that had a second start: M3 or
# M1 directly, M4 or M2 through the M3 or M1 it is warm-started from. (Without
# step 3, M1 and M2 are unchanged and the best stays M3 or M4.)
#
# LCMNL ladders: the same files as dev/collect_standard_reruns.R.
# Each fit is one apollo_estimate at 3000 MLHS draws (klue >= 0.10.0) with
# Apollo's default draw seed, so the script is deterministic. Resumable: an
# existing zero-start fit is reused; a benchmark fit replaced by its zero
# start, and the correlated model refitted from it, carry second_start = TRUE
# and are not refitted again; a refit that an interrupted run left pending is
# finished even when the better fit took the dataset out of the flagged set.
#
# Run in the single lane of dev/rerun_standard_apollo.sh, after
# dev/rerun_mmnl_benchmark.R, or alone; never beside another Apollo lane (two
# lanes overflowed the 32 GB machine on 2026-09-30). From the repository root,
# on one core, then refresh the benchmark table:
#   KLUE_CORES=1 KLUE_MMNL_MEMORY_SAVER=TRUE OMP_NUM_THREADS=1 \
#     OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/run_mmnl_second_start.R
#   Rscript dev/rerun_mmnl_benchmark.R summary
# Decisions and the fits a run would make, no estimation (base R, safe beside
# a running lane):
#   KLUE_SECOND_START_PLAN=1 Rscript dev/run_mmnl_second_start.R
# Test (a few minutes): KLUE_SECOND_START_TEST=1 fits a tiny Mode benchmark
# (first TEST_N = 100 respondents, TEST_DRAWS = 25 draws, the recipe of
# dev/rerun_mmnl_benchmark.R) and runs steps 1-3 on it with MARGIN_BIC = Inf
# and MIN_GAIN = -Inf, so every fit, replacement and refit runs. It reads
# output/ only for the LCMNL BICs and writes only under
# KLUE_SECOND_START_SMOKE (default: klue_second_start_test in the system temp
# directory). A second run must estimate nothing (no START line):
#   KLUE_SECOND_START_TEST=1 KLUE_SECOND_START_SMOKE=<dir> KLUE_CORES=1 \
#     OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/run_mmnl_second_start.R
# Writes: output/mmnl_bench_std/robust/<ds>_M3_zero.rds and <ds>_M1_zero.rds;
#         kept copies in robust/ (<ds>_M3_<start>.rds,
#         <ds>_M4_before_second_start.rds, and the same for M1 and M2);
#         possibly replaced <ds>_M1.rds .. <ds>_M4.rds;
#         output/mmnl_second_start.csv (not in plan mode): the decision before
#         and after, then per step the zero-start LL (M3_zero_LL, M1_zero_LL),
#         whether the benchmark M3 / M1 is its zero start (M3_replaced,
#         M1_replaced), whether M4 / M2 was refitted from it (M4_refit,
#         M2_refit), their nested_ok, and whether step 3 ran (M1_step). These
#         describe the benchmark after the run, so a rerun still reports a
#         replacement made by an earlier run.
# =============================================================================

MARGIN_BIC <- 10
N_DRAWS    <- 3000L
MIN_GAIN   <- 0.01           # LL gain needed for a zero start to replace a benchmark fit
BENCH      <- "output/mmnl_bench_std"
OUT_CSV    <- "output/mmnl_second_start.csv"
PLAN_ONLY  <- identical(Sys.getenv("KLUE_SECOND_START_PLAN"), "1")
TEST       <- identical(Sys.getenv("KLUE_SECOND_START_TEST"), "1")
TEST_N     <- 100L           # TEST: the first respondents kept
TEST_DRAWS <- 25L

stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     paste0(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

SPECS <- list(   # as in dev/rerun_mmnl_benchmark.R
  Vittel      = list(script = "R/empirical_application.R", skip = "empirical.skip",
                     loader = "load_empirical_database",   dgp = "EMP_DGP",
                     lcmnl  = "output/vittel_respec_lcmnl_ext.csv"),
  Mode        = list(script = "R/empirical_mode_choice.R", skip = "mode.skip",
                     loader = "load_mode_choice_database", dgp = "MODE_DGP",
                     lcmnl  = "output/mode_lcmnl_results.csv"),
  SwissRoute  = list(script = "R/empirical_swiss_route.R", skip = "swiss.skip",
                     loader = "load_swiss_database",       dgp = "SWISS_DGP",
                     lcmnl  = "output/swiss_lcmnl_results.csv"),
  Electricity = list(script = "R/empirical_electricity.R", skip = "elec.skip",
                     loader = "load_electricity_database", dgp = "ELEC_DGP",
                     lcmnl  = "output/emp_cmax_extend_Electricity.csv"),
  Swissmetro  = list(script = "R/empirical_swissmetro.R",  skip = "sm.skip",
                     loader = "load_swissmetro_database",  dgp = "SM_DGP",
                     lcmnl  = "output/emp_cmax_extend_Swissmetro.csv")
)

if (TEST) {
  SMOKE <- Sys.getenv("KLUE_SECOND_START_SMOKE")
  if (!nzchar(SMOKE)) SMOKE <- file.path(dirname(tempdir()), "klue_second_start_test")
  SPECS      <- SPECS["Mode"]
  BENCH      <- file.path(SMOKE, "mmnl_bench_std")
  OUT_CSV    <- file.path(SMOKE, "mmnl_second_start.csv")
  N_DRAWS    <- TEST_DRAWS
  MARGIN_BIC <- Inf          # flag the test dataset whatever its margin
  MIN_GAIN   <- -Inf         # any converged zero start replaces, so every refit runs
  stamp("TEST: Mode, first ", TEST_N, " respondents, ", N_DRAWS, " draws, files under ", SMOKE)
}
ROBUST <- file.path(BENCH, "robust")

stage_file <- function(ds, m) file.path(BENCH, sprintf("%s_%s.rds", ds, m))
read_stage <- function(ds, m) { f <- stage_file(ds, m); if (file.exists(f)) readRDS(f) else NULL }
zero_file  <- function(ds, m) file.path(ROBUST, sprintf("%s_%s_zero.rds", ds, m))
conv_bic   <- function(o) if (!is.null(o) && isTRUE(o$converged)) o$BIC else NA_real_
conv_ll    <- function(o) if (!is.null(o) && isTRUE(o$converged)) o$LL else -Inf
# The benchmark `ind` (M3 or M1) is its zero start but `cor` (M4 or M2) was not
# refitted from it yet: an interrupted run.
pending    <- function(ds, ind, cor)
  isTRUE(read_stage(ds, ind)$second_start) && !isTRUE(read_stage(ds, cor)$second_start)
# Step 3 runs when the best MMNL has fixed constants (TEST: always).
fixed_best <- function(best) TEST || best %in% c("M1", "M2")

# Margin = best converged MMNL BIC - BIC-best LCMNL BIC (> 0 favours the LCMNL).
decide <- function(ds) {
  lc  <- read.csv(SPECS[[ds]]$lcmnl)
  lc  <- lc[lc$converged %in% TRUE, ]
  bic <- vapply(paste0("M", 1:4), function(m) conv_bic(read_stage(ds, m)), numeric(1))
  has_m3 <- file.exists(stage_file(ds, "M3"))
  best   <- if (all(is.na(bic))) NA_real_ else min(bic, na.rm = TRUE)
  margin <- best - min(lc$BIC)
  data.frame(dataset = ds, lcmnl_C = lc$C[which.min(lc$BIC)], lcmnl_BIC = round(min(lc$BIC), 1),
             best_mmnl = if (is.na(best)) NA else names(which.min(bic)),
             best_mmnl_BIC = round(best, 1), margin = round(margin, 1),
             benchmark_ready = has_m3,
             second_start = has_m3 && (is.na(margin) || margin > -MARGIN_BIC),
             stringsAsFactors = FALSE)
}

# What the step for `ind` (M3 or M1) and `cor` (M4 or M2) would estimate, from
# the files and with the tests of pair_step(). Prints one line and returns the
# number of fits listed.
pair_plan <- function(ds, ind, cor, label = "") {
  mi <- read_stage(ds, ind); mc <- read_stage(ds, cor); zf <- zero_file(ds, ind)
  z  <- if (file.exists(zf)) readRDS(zf)
  took    <- function(o) if (is.null(o$secs)) "" else sprintf(" [stage fit %.0f min]", o$secs / 60)
  is_zero <- isTRUE(mi$second_start)
  beats   <- !is_zero && !is.null(z) && conv_ll(z) > conv_ll(mi) + MIN_GAIN
  cor_due <- !isTRUE(mc$second_start)
  fits <- c(if (is.null(z)) sprintf("%s from start = \"zero\"%s", ind, took(mi)),
            if (cor_due && (is_zero || beats))
              sprintf("%s from the zero-start %s%s", cor, ind, took(mc)),
            if (cor_due && !is_zero && is.null(z))
              sprintf("then %s from it if it beats the benchmark %s (LL %.3f) by more than %g%s",
                      cor, ind, conv_ll(mi), MIN_GAIN, took(mc)))
  note <- if (is_zero) sprintf("the benchmark %s is its zero start%s", ind,
                               if (cor_due) "" else sprintf(" and %s was refitted from it", cor))
          else if (beats) sprintf("the %s zero start on file beats the benchmark %s", ind, ind)
          else if (!is.null(z))
            sprintf("the %s zero start on file (LL %.3f) does not beat the benchmark %s (LL %.3f)",
                    ind, conv_ll(z), ind, conv_ll(mi))
  stamp("plan ", ds, " ", ind, "/", cor, label, ": ",
        if (length(fits)) paste("estimate", paste(fits, collapse = "; ")) else "estimate nothing",
        if (length(note)) paste0(" (", note, ")"))
  invisible(length(fits))
}

# ---- Estimation (not reached in plan mode) -----------------------------------
setup_klue <- function() {
  if (exists("klue_mmnl", mode = "function")) return(invisible())
  n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else 1L
  options(klue.cores = n_cores, klue.mmnl.n_cores = 1L, mc.cores = n_cores)
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
  stopifnot(packageVersion("klue") >= "0.10.0")
  source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "mclust", "cluster"))  # stops on drift
}

# Random sets and fitting closure of a dataset; its adapter is sourced once.
DATA <- list()
load_ds <- function(ds) {
  if (!is.null(DATA[[ds]])) return(DATA[[ds]])
  s <- SPECS[[ds]]
  options(structure(list(TRUE), names = s$skip))   # load the adapter, skip its autorun
  suppressWarnings(suppressMessages(source(s$script, local = FALSE)))
  db  <- get(s$loader)()
  if (TEST) db <- db[db$ID %in% head(unique(db$ID), TEST_N), ]
  dgp <- get(s$dgp)
  fixed_r <- c(paste0("x", seq_len(dgp$n_generic)), "price")
  DATA[[ds]] <<- list(fixed_r = fixed_r, rand_r = c(fixed_r, paste0("asc", seq_len(dgp$n_asc))),
                      mm = function(...) klue_mmnl(db, n_draws = N_DRAWS, dgp = dgp, ...))
  DATA[[ds]]
}

fit <- function(label, fn) {
  stamp("START ", label)
  t <- system.time(obj <- fn())[3]
  obj$secs <- unname(t)
  stamp(sprintf("DONE  %s  LL=%.3f k=%d BIC=%.1f conv=%s nested_ok=%s (%.1f min)",
                label, obj$LL, obj$k, obj$BIC, obj$converged,
                if (is.null(obj$nested_ok)) "-" else obj$nested_ok, t / 60))
  obj
}

# Second start for `ind` (M3 or M1) and `cor` (M4 or M2, which nests it and is
# warm-started from it), with random parameters `rnd`: steps 1-2 of the header.
# Returns the step's columns of output/mmnl_second_start.csv.
pair_step <- function(ds, ind, cor, rnd, mm) {
  mi <- read_stage(ds, ind)
  stopifnot(identical(mi$random, rnd))
  zf <- zero_file(ds, ind)
  if (file.exists(zf)) { stamp("SKIP ", ds, " ", ind, " zero start (exists)"); z <- readRDS(zf) }
  else { z <- fit(paste(ds, ind, "zero start"), function() mm(random = rnd, start = "zero"))
         saveRDS(z, zf) }

  if (!isTRUE(mi$second_start) && conv_ll(z) > conv_ll(mi) + MIN_GAIN) {
    prev <- mi$selected_from %||% "pooled_mnl"
    keep <- file.path(ROBUST, sprintf("%s_%s_%s.rds", ds, ind, prev))
    if (!file.exists(keep)) file.copy(stage_file(ds, ind), keep)
    z$selected_from <- "zero"; z$second_start <- TRUE
    saveRDS(z, stage_file(ds, ind))
    stamp(sprintf("%s benchmark %s replaced by the zero start (LL %.3f > %.3f); previous kept as %s",
                  ds, ind, z$LL, mi$LL, keep))
    mi <- z
  } else if (isTRUE(mi$second_start)) {
    stamp(ds, " benchmark ", ind, " is the zero start already (replaced by an earlier run)")
  } else stamp(ds, " benchmark ", ind, " unchanged: the zero start does not reach a higher LL")

  # `cor` must start from the benchmark `ind`; refit it when `ind` came from the second start.
  mc <- read_stage(ds, cor)
  if (isTRUE(mi$second_start) && isTRUE(mi$converged) && !isTRUE(mc$second_start)) {
    if (!is.null(mc)) file.copy(stage_file(ds, cor),
                                file.path(ROBUST, sprintf("%s_%s_before_second_start.rds", ds, cor)))
    mc <- fit(sprintf("%s %s from the second-start %s", ds, cor, ind), function()
      mm(random = rnd, correlation = TRUE, warm_start = mi))
    mc$second_start <- TRUE
    saveRDS(mc, stage_file(ds, cor))
  } else if (isTRUE(mc$second_start)) stamp(ds, " ", cor, " already refitted from the second-start ", ind)
  r <- list(round(conv_ll(z), 3), isTRUE(mi$second_start), isTRUE(mc$second_start),
            if (is.null(mc$nested_ok)) NA else mc$nested_ok)
  names(r) <- c(paste0(ind, c("_zero_LL", "_replaced")), paste0(cor, c("_refit", "_nested_ok")))
  r
}

RESULT_ROW <- data.frame(dataset = NA_character_, M3_zero_LL = NA_real_, M3_replaced = NA,
                         M4_refit = NA, M4_nested_ok = NA, M1_step = FALSE, M1_zero_LL = NA_real_,
                         M1_replaced = NA, M2_refit = NA, M2_nested_ok = NA,
                         stringsAsFactors = FALSE)

second_start <- function(ds, flagged) {
  d   <- load_ds(ds)
  row <- RESULT_ROW; row$dataset <- ds
  if (flagged || pending(ds, "M3", "M4")) {
    r <- pair_step(ds, "M3", "M4", d$rand_r, d$mm); row[names(r)] <- r
  }
  # Step 3 is decided on the benchmark as steps 1-2 left it.
  best <- decide(ds)$best_mmnl
  row$M1_step <- (flagged && fixed_best(best)) || pending(ds, "M1", "M2")
  if (row$M1_step) {
    stamp(ds, " M1/M2 step (best MMNL after M3/M4: ", best, ")")
    r <- pair_step(ds, "M1", "M2", d$fixed_r, d$mm); row[names(r)] <- r
  }
  row
}

# TEST: M1-M4 at the test size, as dev/rerun_mmnl_benchmark.R fits them.
test_bench <- function(ds) {
  d <- load_ds(ds)
  dir.create(BENCH, recursive = TRUE, showWarnings = FALSE)
  for (p in list(c("M1", "M2"), c("M3", "M4"))) {
    rnd <- if (p[1] == "M1") d$fixed_r else d$rand_r
    if (!file.exists(stage_file(ds, p[1])))
      saveRDS(fit(paste(ds, p[1], "test benchmark"), function() d$mm(random = rnd)),
              stage_file(ds, p[1]))
    mi <- read_stage(ds, p[1])
    if (isTRUE(mi$converged) && !file.exists(stage_file(ds, p[2])))
      saveRDS(fit(paste(ds, p[2], "test benchmark"), function()
        d$mm(random = rnd, correlation = TRUE, warm_start = mi)), stage_file(ds, p[2]))
  }
}

# ---- Plan --------------------------------------------------------------------
if (TEST && !PLAN_ONLY) { setup_klue(); test_bench("Mode") }

plan <- do.call(rbind, lapply(names(SPECS), decide))
stamp(sprintf("rule: second start where best MMNL BIC - LCMNL BIC > -%g", MARGIN_BIC))
print(plan, row.names = FALSE)
for (ds in plan$dataset[!plan$benchmark_ready])
  stamp("WARNING ", ds, ": no benchmark M3 yet; run dev/rerun_mmnl_benchmark.R first")

# Datasets to process: the flagged ones, and any with a refit that an
# interrupted run left pending (the better fit may have taken the dataset out
# of the flagged set).
p34  <- vapply(plan$dataset, pending, logical(1), ind = "M3", cor = "M4")
p12  <- vapply(plan$dataset, pending, logical(1), ind = "M1", cor = "M2")
todo <- which(plan$second_start | p34 | p12)
for (i in todo) {
  ds <- plan$dataset[i]; fl <- plan$second_start[i]
  n34 <- if (fl || p34[i]) pair_plan(ds, "M3", "M4") else 0L
  if ((fl && fixed_best(plan$best_mmnl[i])) || p12[i]) {
    pair_plan(ds, "M1", "M2", sprintf(" (best MMNL %s%s)", plan$best_mmnl[i],
                                      if (n34 > 0) ", decided again after M3/M4" else ""))
  } else if (fl) stamp("plan ", ds, " M1/M2: nothing (best MMNL ", plan$best_mmnl[i],
                       " has random constants)")
}
if (!length(todo)) stamp("plan: no dataset flagged, nothing to estimate")
if (PLAN_ONLY) { stamp("plan mode: nothing estimated or written"); quit(save = "no") }

# ---- Run ---------------------------------------------------------------------
setup_klue()
dir.create(ROBUST, recursive = TRUE, showWarnings = FALSE)
res <- lapply(todo, function(i) second_start(plan$dataset[i], plan$second_start[i]))
after <- do.call(rbind, lapply(plan$dataset, decide))[, c("dataset", "best_mmnl", "best_mmnl_BIC", "margin")]
names(after)[-1] <- paste0(names(after)[-1], "_after")
out <- merge(plan, after, by = "dataset", sort = FALSE)
out <- merge(out, if (length(res)) do.call(rbind, res) else RESULT_ROW[0, ],
             by = "dataset", all.x = TRUE, sort = FALSE)
out <- out[match(names(SPECS), out$dataset), ]
write.csv(out, OUT_CSV, row.names = FALSE)
print(out, row.names = FALSE)
stamp("DONE")
