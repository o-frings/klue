#!/usr/bin/env Rscript
# =============================================================================
# dev/check_ladder_carryover.R
#
# Why. Each cell of the two package stress ladders keeps, per arm, the LL of
# the run that first fitted that arm. From 29 September to 5 October 2026 only
# the Apollo arms were rerun with the current harness (apollo_searchStart at
# its defaults and apollo_lcEM on both ladders, apollo_clust on the blocked
# one). The other arms keep LLs written from June to
# September by older klue versions and an older harness: klue_ml, klue_em,
# klue_em_rp, gmnl_lc and apollo_ss_published (Apollo's reduced search). This
# script refits those arms with the current code and the pinned packages and
# compares each LL with the stored one.
#
# Backs. The LLs of these arms in the blocked ladder (tab:stress_blocked) and
# the randomised ladder (tab:stress), and the gaps and miss counts built on
# them; the appendix statement that the reduced-search results come from an
# earlier harness and that the current code with Apollo 0.3.5 reproduces them
# on the cell checked (blocked design, moderate rung, seed 2; stored LL
# -2374.2962306), here extended to seed 1 of every rung of both ladders; the
# availability statement that the analyses ran in R 4.3.1 with klue 0.10.0,
# Apollo 0.3.5 and gmnl 1.1-3.2 with mlogit 1.0-3.1; and the cover letter's
# statement that the analyses run from seeded scripts.
#
# Cells and arms, fixed in advance:
#   blocked ladder (dev/stress_replicate_blocked.R; moderate, hard, very_hard
#   x seeds 1-10 = 30 cells on one shared blocked D-efficient design)
#     klue_ml, klue_em, klue_em_rp, gmnl_lc  every cell
#     apollo_ss_published                    moderate|seed2 (the cell the paper
#                                            cites), then seed 1 of each rung
#   randomised ladder (dev/stress_replicate.R; moderate, hard, very_hard,
#   hard_corr x seeds 1-5 = 20 cells)
#     klue_ml, klue_em, gmnl_lc              every cell; this ladder has no
#                                            klue_em_rp arm, so nothing to compare
#     apollo_ss_published                    seed 1 of each rung
# Each arm is the ladder's own METHODS entry, which calls the current runner in
# dev/compare_packages.R: run_klue ("ml", "em"), run_klue_em_rp (partition seed
# = the cell seed), run_gmnl, and run_apollo_ss_published_blocked or
# run_apollo_ss_published_randomised (SS_PUBLISHED_BLOCKED or
# SS_PUBLISHED_RANDOMISED: generic start, class-major order, delta_1 fixed).
# The runners loop over their cells at top level and cannot be sourced, so this
# script parses each runner and evaluates its top-level DGP, DESIGN, SEEDS,
# RUNGS, N_PER_CLASS and METHODS assignments in file order.
#
# Data. Regenerated as the runners' loops generate them: cell seed
# as.integer(1000*K + 100*(kap*100) + seed), written exactly as there
# (moderate 11500 + seed, hard and hard_corr 10000 + seed, very_hard 8000 +
# seed), 60 respondents per class. Blocked: klue_simulate(design = DESIGN),
# DESIGN = klue_design(n_cards = 48, n_blocks = 4), 12 cards per respondent.
# Randomised: T_tasks = 12, attr_corr = 0.6 on hard_corr. klue_ml reproducing
# its stored LL in a cell is the check that the data are the ladder's (column
# data_ok).
#
# Seeds and the RNG state. klue_design calls set.seed(20240601) and
# klue_simulate set.seed(cell seed) before drawing, so design and data do not
# depend on earlier draws. No carried-over arm reads the random-number stream
# it inherits:
#   klue_ml, klue_em     kmeans, Mclust and pam run under set.seed(123) (klue's
#                        .with_seed, which since 0.9.2 restores the caller's
#                        stream afterwards; before, set.seed(123) at the same
#                        points); the three hclust starts and BFGS or EM use no
#                        RNG
#   klue_em_rp           partition p (1-6) is drawn under
#                        set.seed(cell_seed*100000 + 50000 + p)
#                        (get_random_partition_starts), then EM
#   gmnl_lc              gmnl 1.1-3.2 starts the latent-class model at the
#                        pooled MNL (maxLik Newton-Raphson) shifted by
#                        seq(-0.02, 0.02, length.out = Q), equal shares, then
#                        BFGS; its LC likelihood has no draws
#   apollo_ss_published  apollo_searchStart calls set.seed(17) (Apollo's default
#                        seed 13, plus 4) before it draws its candidates; before
#                        that only apollo_makeLogLike draws, a cache key for
#                        iteration counting; apollo_estimate's draw (its test
#                        that each parameter moves the LL) comes after the
#                        reseed and only decides whether to stop
# So the order of the arms within a cell, which differed between the runs that
# wrote the stored LLs, does not matter. As a guard, every fit here starts from
# the RNG state the cell's klue_simulate call leaves, which is the state the
# runner's first arm of the cell started from, whatever the order of the fits
# or a resume. The test mode refits each of its arms from an unrelated RNG
# state and reports the LL difference (expected 0).
#
# Rule. reproduced = |LL_now - LL_stored| <= 1e-3. A klue arm is fitted first
# with options(klue.screen = TRUE), klue's default since 0.9.1. Where it
# misses, it is refitted on that cell with options(klue.screen = FALSE), the
# pre-0.9.1 path that README_REPRODUCE.md ("Package history") prescribes for
# results produced before 0.9.1, such as the clustering arms of both ladders
# (since 2026-10-06 the ladder runners set it). Both rows are
# kept (column klue_screen); the summary counts a (cell, arm) as reproduced
# when either row is. An arm the runner fits with klue.screen = FALSE (its
# SCREEN_OFF, since 2026-10-06 klue_ml) is also refitted that way in every
# cell, so the runner's own setting is checked everywhere. gmnl_lc and
# apollo_ss_published do not read the option (klue_screen NA).
#
# Run from the repository root, one lane per ladder, one core each (Apollo
# nCores = 1; R's reference BLAS is single-threaded):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 KLUE_CORES=1 \
#     nohup nice -n 15 Rscript dev/check_ladder_carryover.R blocked \
#     > output/ladder_carryover_blocked.log 2>&1 &
#   (the same with randomised, logging to output/ladder_carryover_randomised.log)
#   Rscript dev/check_ladder_carryover.R summary   # no fitting; safe while lanes run
# A lane fits klue_ml in every cell, then klue_em, klue_em_rp and gmnl_lc, and
# the reduced Apollo searches last; it ends with the summary.
# Cost (measured 2026-10-06, one core per lane): blocked 1 h 40 min (four
# reduced searches of 7.4-14 min), randomised 1 h 6 min.
# Resumable: every finished (cell, arm, klue.screen setting) is written to the
# lane's checkpoint (to a temporary file, then renamed) and skipped by a rerun;
# a fit whose runner catches an error (LL NA, with the runner's note) is saved
# as not reproduced; an error that escapes the runner, such as the version-pin
# stop, is logged, not saved, so a rerun retries it.
# The checkpoint records the versions in use; a rerun under other versions
# stops (delete the checkpoint to redo). One process per lane (lock file).
#
# Reads output/stress_replicate_blocked.rds and output/stress_replicate.rds
# (copied first; never written) and the two runners (parsed, not run).
# Writes:
#   output/ladder_carryover_check_<lane>.rds  checkpoint per lane: one record
#     per (cell, arm, klue.screen setting), the cells planned per arm, versions
#   output/ladder_carryover_check.csv  one row per record: ladder, rung, seed,
#     cell_seed, method, LL_stored, LL_now, diff, reproduced, klue_screen,
#     secs_now, converged, note (the runner's), data_ok, R_version,
#     klue_version, apollo_version, gmnl_version, mlogit_version (as loaded),
#     date
#   output/.ladder_carryover_<lane>.lock  while a lane runs
#
# Test (about 3 min for blocked, most of it the design; under 1 min for
# randomised):
# CARRY_TEST=1 fits moderate|seed1 of the lane with klue_ml (also with
# klue.screen = FALSE, to exercise that path) and gmnl_lc, refits each from an
# unrelated RNG state, and writes only under CARRY_TEST_DIR (default:
# ladder_carryover_test next to R's tempdir()), starting fresh each time:
#   CARRY_TEST=1 CARRY_TEST_DIR=<dir> Rscript dev/check_ladder_carryover.R blocked
#   CARRY_TEST=1 CARRY_TEST_DIR=<dir> Rscript dev/check_ladder_carryover.R randomised
#   CARRY_TEST=1 CARRY_TEST_DIR=<dir> Rscript dev/check_ladder_carryover.R summary
# Stored LLs of that cell: blocked klue_ml -2351.7176329, gmnl_lc -2351.7176428;
# randomised klue_ml -2767.0703019, gmnl_lc -2767.0708177.
#
# Result (run 2026-10-06, output/ladder_carryover_check.csv): all 188 planned
# (cell, arm) checks reproduce. At klue's default, klue_ml reproduces the
# blocked ladder to within 4e-5, but on the randomised ladder it ends elsewhere
# in 5 of 20 cells (hard seed 5 -3.62, very_hard seeds 1 and 4 -1.52 and
# -2.29, hard_corr seed 1 +2.71, hard_corr seed 2 -0.0013); with
# klue.screen = FALSE, the runners' setting for klue_ml, all 50 klue_ml cells
# of both ladders match exactly (max |diff| 0). klue_em (both ladders, within
# 2e-6), klue_em_rp, gmnl_lc and the reduced search reproduce at the default,
# the last three exactly. Four of the eight reduced-search cells sit at the
# cell's best LL, where any harness that reaches the optimum agrees; the other
# four (blocked moderate seed 2 at -61.8 LL; randomised moderate, very_hard and
# hard_corr seed 1) are local optima and do tell the harnesses apart.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
LANE <- if (length(args) == 1L) args else ""
if (!LANE %in% c("blocked", "randomised", "summary"))
  stop("usage: Rscript dev/check_ladder_carryover.R blocked|randomised|summary", call. = FALSE)
if (!file.exists("dev/compare_packages.R")) stop("run from the repository root", call. = FALSE)
TEST <- identical(Sys.getenv("CARRY_TEST"), "1")

OUT_DIR <- "output"
if (TEST) {
  OUT_DIR <- Sys.getenv("CARRY_TEST_DIR")
  if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(dirname(tempdir()), "ladder_carryover_test")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (TEST && normalizePath(OUT_DIR) == normalizePath("output", mustWork = FALSE))
  stop("CARRY_TEST_DIR points at output/; the test writes elsewhere", call. = FALSE)
CSV     <- file.path(OUT_DIR, "ladder_carryover_check.csv")
ckpt_of <- function(lane) file.path(OUT_DIR, sprintf("ladder_carryover_check_%s.rds", lane))

# The two ladders: runner, stored results, and the cells of the reduced Apollo
# search (fixed in advance; the paper's cell first).
LADDERS <- list(
  blocked    = list(runner  = "dev/stress_replicate_blocked.R",
                    file    = "output/stress_replicate_blocked.rds",
                    reduced = c("moderate|seed2", "moderate|seed1", "hard|seed1",
                                "very_hard|seed1")),
  randomised = list(runner  = "dev/stress_replicate.R",
                    file    = "output/stress_replicate.rds",
                    reduced = c("moderate|seed1", "hard|seed1", "very_hard|seed1",
                                "hard_corr|seed1")))
ARMS       <- c("klue_ml", "klue_em", "klue_em_rp", "gmnl_lc",
                "apollo_ss_published")                 # fitted in this order
KLUE_ARMS  <- c("klue_ml", "klue_em", "klue_em_rp")   # the arms that read options(klue.screen)
RUNG_ORDER <- c("moderate", "hard", "very_hard", "hard_corr")
TOL        <- 1e-3                                    # reproduced: |LL_now - LL_stored| <= TOL
TEST_CELL  <- "moderate|seed1"
TEST_ARMS  <- c("klue_ml", "gmnl_lc")
TEST_RNG_SEED <- 20261006L                            # the unrelated RNG state of the test refits
VERSION_COLS  <- c("R_version", "klue_version", "apollo_version", "gmnl_version",
                   "mlogit_version")

stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     sprintf(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a)) b else a
read_ladder <- function(file) {   # copy first, in case a ladder job is rewriting it
  for (i in 1:5) {
    tmp <- tempfile(fileext = ".rds")
    ok <- file.copy(file, tmp, overwrite = TRUE)
    x <- if (ok) tryCatch(readRDS(tmp), error = function(e) NULL)
    unlink(tmp)
    if (!is.null(x)) return(x)
    Sys.sleep(5)
  }
  stop("could not read ", file, call. = FALSE)
}
# Write to a temporary name in the same folder, then rename: readers (the
# summary, the other lane) never see a half-written file. The name starts with
# a dot, so a leftover is not published.
write_atomic <- function(path, writer) {
  tmp <- tempfile(paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  writer(tmp)
  if (!file.rename(tmp, path)) stop("could not write ", path, call. = FALSE)
}

# ---- Summary: merge the lanes' checkpoints into the CSV, print per arm ------
summarise <- function() {
  ck <- list()
  for (ln in names(LADDERS)) {
    f <- ckpt_of(ln)
    if (file.exists(f)) ck[[ln]] <- tryCatch(readRDS(f), error = function(e) {
      stamp("NOTE could not read %s: %s", f, conditionMessage(e)); NULL })
  }
  recs <- unname(unlist(lapply(ck, `[[`, "records"), recursive = FALSE))
  if (!length(recs)) { stamp("no records under %s yet", OUT_DIR); return(invisible(NULL)) }
  tab <- do.call(rbind, lapply(recs, as.data.frame, stringsAsFactors = FALSE))
  cell <- paste(tab$ladder, tab$rung, tab$seed)
  ml <- tab$method == "klue_ml"
  ok <- if (any(ml)) tapply(tab$reproduced[ml], cell[ml], any) else logical(0)
  tab$data_ok <- unname(ok[cell])          # NA until klue_ml has run in the cell
  tab <- tab[order(match(tab$ladder, names(LADDERS)), match(tab$rung, RUNG_ORDER), tab$seed,
                   match(tab$method, ARMS), !(tab$klue_screen %in% TRUE)), ]
  cols <- c("ladder", "rung", "seed", "cell_seed", "method", "LL_stored", "LL_now", "diff",
            "reproduced", "klue_screen", "secs_now", "converged", "note", "data_ok",
            VERSION_COLS, "date", if (TEST) "rng_check_diff")
  out <- tab[, cols]
  out$LL_stored <- round(out$LL_stored, 7); out$LL_now <- round(out$LL_now, 7)
  out$diff <- signif(out$diff, 3); out$secs_now <- round(out$secs_now, 1)
  write_atomic(CSV, function(f) write.csv(out, f, row.names = FALSE))

  # The row that decides a (cell, arm): the default fit, or the klue.screen =
  # FALSE refit where the default missed.
  key <- paste(tab$ladder, tab$rung, tab$seed, tab$method)
  dec <- do.call(rbind, lapply(split(tab, factor(key, unique(key))), function(d) {
    d0 <- d[!(d$klue_screen %in% FALSE), , drop = FALSE]
    d1 <- d[d$klue_screen %in% FALSE, , drop = FALSE]
    if (!nrow(d0)) d0 <- d1
    pick <- if (!isTRUE(d0$reproduced[1]) && nrow(d1)) d1[1, ] else d0[1, ]
    pick$reproduced_default <- isTRUE(d0$reproduced[1])
    pick$diff_default <- d0$diff[1]
    pick
  }))
  mx <- function(v) if (any(is.finite(v))) max(abs(v[is.finite(v)])) else NA_real_
  per_arm <- do.call(rbind, lapply(names(ck), function(ln) {
    plan <- ck[[ln]]$plan %||% integer(0)
    do.call(rbind, lapply(intersect(ARMS, union(names(plan), dec$method[dec$ladder == ln])),
                          function(m) {
      d <- dec[dec$ladder == ln & dec$method == m, , drop = FALSE]
      rep <- d$reproduced %in% TRUE
      data.frame(ladder = ln, method = m, planned = unname(plan[m]) %||% NA_integer_,
                 checked = nrow(d), reproduced = sum(rep),
                 at_default = sum(d$reproduced_default),
                 via_screen_false = sum(rep & !d$reproduced_default),
                 not_reproduced = sum(!rep),
                 max_abs_diff = signif(mx(d$diff), 3),
                 max_abs_diff_default = signif(mx(d$diff_default), 3),
                 not_reproduced_in = paste(paste0(d$rung, "|seed", d$seed)[!rep], collapse = " "),
                 stringsAsFactors = FALSE)
    }))
  }))
  old <- options(width = 220); on.exit(options(old), add = TRUE)
  cat("\n===== carried-over arms refitted with the current code (reproduced: |diff| <= 1e-3) =====\n")
  cat("(deciding row per cell: the default fit, or the klue.screen = FALSE refit where the default missed;\n",
      " max_abs_diff is over the deciding rows, max_abs_diff_default over the default fits)\n", sep = "")
  print(per_arm, row.names = FALSE)
  sf <- tab[tab$method %in% KLUE_ARMS & tab$klue_screen %in% FALSE, ]
  for (ln in unique(sf$ladder)) for (m in unique(sf$method[sf$ladder == ln])) {
    d <- sf[sf$ladder == ln & sf$method == m, ]
    stamp("%s %s with klue.screen = FALSE: reproduced in %d of the %d cells fitted that way (max |diff| %.1e)",
          ln, m, sum(d$reproduced %in% TRUE), nrow(d), mx(d$diff))
  }
  for (ln in unique(dec$ladder)) {
    d <- dec[dec$ladder == ln & dec$method == "klue_ml", ]
    stamp("data check, %s: klue_ml reproduces its stored LL in %d of %d cells checked", ln,
          sum(d$reproduced %in% TRUE), nrow(d))
  }
  p <- dec[dec$ladder == "blocked" & dec$rung == "moderate" & dec$seed == 2 &
             dec$method == "apollo_ss_published", ]
  if (nrow(p))
    stamp("the paper's cell (blocked moderate|seed2, apollo_ss_published): LL %.7f, stored %.7f, diff %+.1e -> %s",
          p$LL_now, p$LL_stored, p$diff, if (isTRUE(p$reproduced)) "reproduced" else "NOT reproduced")
  miss <- dec[!(dec$reproduced %in% TRUE), ]
  if (nrow(miss)) {
    cat("\nNot reproduced (deciding row):\n")
    print(miss[, c("ladder", "rung", "seed", "method", "klue_screen", "LL_stored", "LL_now",
                   "diff", "converged", "note", "data_ok")], row.names = FALSE)
  }
  n_plan <- sum(unlist(lapply(ck, function(x) x$plan)), na.rm = TRUE)
  stamp("wrote %s (%d rows; %d of %d planned (cell, arm) checks done%s)", CSV, nrow(out),
        nrow(dec), n_plan, if (nrow(dec) < n_plan) ", partial" else "")
  invisible(out)
}

if (LANE == "summary") { summarise(); quit(save = "no", status = 0) }

# ---- One lane -----------------------------------------------------------------
LAD  <- LADDERS[[LANE]]
CKPT <- ckpt_of(LANE)
LOCK <- file.path(OUT_DIR, sprintf(".ladder_carryover_%s.lock", LANE))
if (file.exists(LOCK)) {
  pid <- suppressWarnings(as.integer(readLines(LOCK, n = 1, warn = FALSE)))
  if (!is.na(pid) && pid != Sys.getpid() && isTRUE(tools::pskill(pid, 0L))) {
    stamp("the %s lane is already running (pid %d); exiting", LANE, pid)
    quit(save = "no", status = 0)
  }
}
writeLines(as.character(Sys.getpid()), LOCK)

# The runner's own top-level definitions, evaluated in file order: DGP, DESIGN
# (blocked only; klue_design, about 2.5 min), SEEDS, RUNGS, N_PER_CLASS, METHODS,
# SCREEN_OFF (the arms the runner fits with options(klue.screen = FALSE)).
runner_defs <- function(file) {
  want <- c("DGP", "DESIGN", "SEEDS", "RUNGS", "N_PER_CLASS", "METHODS", "SCREEN_OFF")
  env <- new.env(parent = globalenv())
  for (e in parse(file, keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% want) eval(e, env)
  miss <- setdiff(want, c(ls(env), if (LANE == "randomised") "DESIGN"))
  if (length(miss)) stop(file, " has no top-level ", paste(miss, collapse = ", "), call. = FALSE)
  env
}

# A cell's data, as the runner's loop generates them, and the RNG state the
# generation leaves (kept for the lane).
DATA <- new.env()
cell_data <- function(key) {
  if (is.null(DATA[[key]])) {
    rn <- sub("\\|seed[0-9]+$", "", key); sd_i <- as.integer(sub("^.*\\|seed", "", key))
    r <- RUNGS[[rn]]
    seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + sd_i)   # the runners' formula
    d <- if (LANE == "blocked")          # as in dev/stress_replicate_blocked.R
      klue_simulate(N_per_class = RD$N_PER_CLASS, true_K = r$K,
                    separation = r$kap, heterogeneity = r$sig,
                    seed = seed_i, dgp = RD$DGP, design = RD$DESIGN)
    else                                 # as in dev/stress_replicate.R
      klue_simulate(N_per_class = RD$N_PER_CLASS, T_tasks = 12, true_K = r$K,
                    separation = r$kap, heterogeneity = r$sig,
                    seed = seed_i, dgp = RD$DGP, attr_corr = r$corr)
    assign(key, list(key = key, rung = rn, seed = sd_i, K = r$K, seed_i = seed_i,
                     db = d$database, rng = get(".Random.seed", envir = globalenv())),
           envir = DATA)
  }
  DATA[[key]]
}

# One call of the ladder's METHODS entry, as the runner calls it. rng: the
# global RNG state to start from (NULL keeps the current one).
run_arm <- function(cd, arm, rng = cd$rng) {
  f <- RD$METHODS[[arm]]
  if (!is.null(rng)) assign(".Random.seed", rng, envir = globalenv())
  if (length(formals(f)) >= 3L) f(cd$db, cd$K, cd$seed_i) else f(cd$db, cd$K)
}

# One arm on one cell under one klue.screen setting (NA: the arm does not read
# it); the record goes to the checkpoint.
rec_id <- function(key, arm, screen)
  paste(key, arm, if (is.na(screen)) "-" else paste0("screen=", screen), sep = "|")
fit_arm <- function(key, arm, screen) {
  id <- rec_id(key, arm, screen)
  if (!is.null(CK$records[[id]])) return(CK$records[[id]])
  cd <- cell_data(key)
  if (!is.na(screen)) { op <- options(klue.screen = screen); on.exit(options(op), add = TRUE) }
  stamp("RUN  %s %s%s", key, arm, if (is.na(screen)) "" else sprintf(" (klue.screen = %s)", screen))
  fit <- tryCatch(run_arm(cd, arm), error = function(e) {
    stamp("ERROR %s %s: %s (not saved; a rerun retries it)", key, arm, conditionMessage(e))
    NULL })
  if (is.null(fit)) return(NULL)
  LL_stored <- unname(lad[[key]]$LL[[arm]])
  LL_now <- as.numeric(fit$LL)
  rng_diff <- NA_real_
  if (TEST) {   # the same fit from an unrelated RNG state: an arm that reads the stream moves
    set.seed(TEST_RNG_SEED); invisible(runif(1000))
    fit2 <- tryCatch(run_arm(cd, arm, rng = NULL), error = function(e) NULL)
    rng_diff <- if (is.null(fit2)) NA_real_ else as.numeric(fit2$LL) - LL_now
    stamp("  RNG check %s %s: refit from an unrelated RNG state, LL difference %s", key, arm,
          format(rng_diff))
  }
  rec <- c(list(ladder = LANE, rung = cd$rung, seed = cd$seed, cell_seed = cd$seed_i,
                method = arm, LL_stored = LL_stored, LL_now = LL_now,
                diff = LL_now - LL_stored,
                reproduced = isTRUE(abs(LL_now - LL_stored) <= TOL), klue_screen = screen,
                secs_now = as.numeric(fit$seconds), converged = isTRUE(fit$converged),
                note = as.character(fit$note %||% "")),
           VERSIONS, list(date = format(Sys.time(), "%Y-%m-%d %H:%M"), rng_check_diff = rng_diff))
  CK$records[[id]] <<- rec
  write_atomic(CKPT, function(f) saveRDS(CK, f))
  stamp("  %s %s%s: LL %.7f  stored %.7f  diff %+.1e  %s  (%.1f s)%s", key, arm,
        if (is.na(screen)) "" else sprintf(" [screen %s]", screen), LL_now, LL_stored,
        rec$diff, if (rec$reproduced) "reproduced" else "NOT reproduced", rec$secs_now,
        if (nzchar(rec$note)) paste0("  note: ", rec$note) else "")
  rec
}

tryCatch({
  sys.source("dev/compare_packages.R", envir = globalenv())   # functions only (guarded); checks the pins
  options(klue.cores = 1L, klue.mmnl.n_cores = 1L, mc.cores = 1L)   # Apollo: nCores = 1 in .apollo_lc_setup
  lad <- read_ladder(LAD$file)
  stamp("%s lane%s: evaluating the settings and arms of %s%s", LANE, if (TEST) " (TEST)" else "",
        LAD$runner, if (LANE == "blocked") ", which builds the shared design" else "")
  RD <- runner_defs(LAD$runner)
  RUNGS <- setNames(RD$RUNGS, vapply(RD$RUNGS, function(r) r$label, ""))
  arms <- intersect(if (TEST) TEST_ARMS else ARMS, names(RD$METHODS))
  if (length(gone <- setdiff(if (TEST) TEST_ARMS else ARMS, arms)))
    stamp("NOTE %s has no arm %s: nothing to compare", LAD$runner, paste(gone, collapse = ", "))
  if ("gmnl_lc" %in% arms) check_pinned_versions(c("gmnl", "mlogit"))   # stop now, not at the first gmnl fit
  VERSIONS <- list(R_version = paste(R.version$major, R.version$minor, sep = "."),
                   klue_version = .version_in_use("klue"),
                   apollo_version = .version_in_use("apollo"),
                   gmnl_version = .version_in_use("gmnl"),
                   mlogit_version = .version_in_use("mlogit"))

  cells <- unlist(lapply(names(RUNGS), function(rn) paste0(rn, "|seed", RD$SEEDS)))
  if (length(gone <- setdiff(cells, names(lad))))
    stamp("NOTE not in %s, skipped: %s", LAD$file, paste(gone, collapse = " "))
  cells <- intersect(cells, names(lad))
  plan <- setNames(lapply(arms, function(a) {
    ks <- if (TEST) intersect(TEST_CELL, cells)
          else if (a == "apollo_ss_published") intersect(LAD$reduced, cells) else cells
    Filter(function(k) is.finite(lad[[k]]$LL[a] %||% NA_real_), ks)   # a stored LL to compare with
  }), arms)

  CK <- if (!TEST && file.exists(CKPT)) readRDS(CKPT) else list(records = list())   # a test starts fresh
  if (length(CK$records) && !identical(CK$versions, VERSIONS))
    stop(sprintf("%s holds records made under other versions (%s); delete it to redo",
                 CKPT, paste(unlist(CK$versions), collapse = ", ")), call. = FALSE)
  CK$versions <- VERSIONS
  CK$plan <- vapply(plan, length, integer(1))
  write_atomic(CKPT, function(f) saveRDS(CK, f))
  stamp("R %s, klue %s, apollo %s, gmnl %s, mlogit %s (mlogit from %s); BLAS %s", VERSIONS$R_version,
        VERSIONS$klue_version, VERSIONS$apollo_version, VERSIONS$gmnl_version,
        VERSIONS$mlogit_version,
        if (isNamespaceLoaded("mlogit")) dirname(getNamespaceInfo("mlogit", "path")) else "-",
        basename(extSoftVersion()[["BLAS"]]))
  stamp("plan: %s; %d record(s) already in %s", paste(sprintf("%s %d cells", arms, CK$plan),
        collapse = ", "), length(CK$records), CKPT)

  for (arm in arms) {
    stamp("arm %s: %d cell(s)", arm, length(plan[[arm]]))
    for (key in plan[[arm]]) {
      rec <- fit_arm(key, arm, if (arm %in% KLUE_ARMS) TRUE else NA)
      if (arm %in% KLUE_ARMS && !is.null(rec) &&
          (TEST || !rec$reproduced || arm %in% RD$SCREEN_OFF))
        fit_arm(key, arm, FALSE)
    }
  }
  stamp("%s lane finished", LANE)
  summarise()
}, finally = unlink(LOCK))
