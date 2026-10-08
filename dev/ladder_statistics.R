#!/usr/bin/env Rscript
# =============================================================================
# dev/ladder_statistics.R
#
# Why. The two package stress ladders (dev/stress_replicate_blocked.R,
# dev/stress_replicate.R) end their logs with the mean gap to the best
# log-likelihood and the number of seeds more than 0.5 LL below it, per rung
# and arm; dev/run_apollo_isolation_blocked.R prints the misses and median
# times of its arms. The manuscript quotes further statistics from the same
# result files: miss counts at other thresholds, paired tests, miss sizes,
# time ranges, and how often Apollo's EM routine stopped at its iteration cap.
# Until 2026-10-06 these came from scratch scripts and from grep on the logs.
# This script recomputes all of them from the stored files; it estimates
# nothing. From the refits of dev/check_ladder_carryover.R it also scores
# clustering+ML at klue's default on the randomised ladder (the note to Table
# tab:stress).
#
# Statements it backs (version12.tex):
#   Blocked ladder (Table tab:stress_blocked and Figure fig:time, the H1b,
#   Speed, ML-or-EM and Scope paragraphs, the discussion):
#   - the miss counts on the hard and very-hard rungs at 0.1, 1, 2 and 5 LL
#     below the best (table note), and how many of the 30 seeds each arm
#     reaches;
#   - the exact McNemar test of each arm against clustering+ML on those 20
#     seeds, with the seeds missed by only one of the two: two only by
#     clustering+ML and none only by Apollo's default search; Apollo's
#     estimator from the clustering starts misses only seeds that
#     clustering+ML also misses;
#   - the sizes of the misses (clustering+ML 1.0 to 2.0 LL, Apollo's default
#     search 2.1 LL, Apollo's reduced search 62 LL on the moderate rung) and
#     their BIC equivalents of 2 x LL ("about 4", "as much as 124");
#   - that the best of clustering+ML, clustering+EM and RP+EM misses one seed;
#   - the 95% Wilson interval of a miss rate of 3 in 20;
#   - the Time per fit column (the range of the per-rung mean times), the range
#     of the two EM arms together, RP+EM's time relative to clustering+ML (also
#     in the cover letter), and Apollo's reduced search timed on seed 1 of each
#     rung;
#   - how often apollo_lcEM reported its iteration cap and skipped its closing
#     maximum-likelihood step, and how many of its misses that covers.
#   Randomised ladder (Appendix Table tab:stress and the paragraphs on it):
#   - the exact McNemar tests against clustering+ML over the 20 seeds, with 4
#     seeds missed only by clustering+ML and 2 only by Apollo's default search;
#   - the sizes of the misses of clustering+ML and of Apollo's default search;
#   - the cells where Apollo's default search alone reached the best, and by
#     how much;
#   - the Time per fit column, and the shortest, median and longest fit of
#     Apollo's default search, Apollo's reduced search and clustering+ML;
#   - the apollo_lcEM cap count;
#   - the table note on klue 0.10.0's default multistart: clustering+ML
#     refitted at that default ends elsewhere in 5 of the 20 cells and would
#     miss 5 seeds instead of 4, with the exact McNemar test against Apollo's
#     default search (5 seeds against 2, p = 0.45).
#   Isolation runs (Appendix Table tab:apollo_iso, Apollo appendix): the paired
#   exact tests A3 against A1, A4 against A2 and A4 against clustering+ML, and
#   how many times faster the default search ran with analytic gradients (A5)
#   than as shipped, cell by cell.
#
# Definitions, as in the runners. In a ladder cell, best = the highest LL any
# arm reached and gap = LL - best; an arm misses the cell when gap < -0.5 (the
# runners' "stuck"), and misses it at threshold t when gap < -t. The shortfall
# of a miss is -gap; every arm fits the same model, so a shortfall of d LL is a
# BIC difference of 2d. Exact McNemar test: two-sided exact binomial test
# (p = 0.5) on the cells missed by only one of the two arms. Wilson interval:
# wilson_ci() of R/statistical_tests.R. An arm is alone at the best when every
# other arm misses the cell; its margin is its LL minus the next best arm's.
# The isolation runs are scored as in their CSV: against the best LL of any run
# in that file or on the ladder (its miss column). apollo_lcEM, from the cell's
# section of the ladder log ("RUN <cell>" to "done <cell>"): the EM block runs
# from "Initialising EM algorithm" to "EM algorithm stopped: ..."; capped =
# Apollo printed "maximum number of iterations reached!"; last_EM_iteration =
# the last "Starting iteration"; stopped_by_criterion = the last improvement is
# below apollo_lcEM's default EMstoppingCriterion of 1e-5; continued_classical
# = Apollo printed "Continuing with classical estimation...". Apollo 0.3.5
# counts an iteration before it tests the cap of 100, so it reports the cap
# after 99 or more iterations, also when the criterion was met at iteration 99,
# and then skips the closing ML step. Each block is checked against that rule.
# Clustering+ML at klue's default (randomised ladder): the ladder's klue_ml ran
# klue's full-precision multistart, the default before 0.9.1
# (options(klue.screen = FALSE)). dev/check_ladder_carryover.R refitted it at
# the current default, which screens the six starts at a loose tolerance and
# refines only the best (its rows with klue_screen TRUE). That LL (LL_now)
# replaces the stored klue_ml LL as the arm klue_ml_screen_default, and each
# cell's best is recomputed over the arms. The arm moved in a cell when
# |LL_now - LL_stored| > 1e-3, which the carryover check calls not reproduced.
#
# Inputs (read only):
#   output/stress_replicate_blocked.rds   blocked ladder: 30 cells x 8 arms
#   output/stress_replicate.rds           randomised ladder: 20 cells x 6 arms
#   output/stress_replicate_blocked_std.log, output/stress_replicate_std.log
#                                         the two ladder runs' console logs
#   output/apollo_isolation_blocked.csv   isolation runs: 10 cells x 12 arms
#   output/apollo_smartstart_diag.csv     Apollo's reduced search timed on seed
#                                         1 of each blocked rung (bounded_secs)
#   output/ladder_carryover_check.csv     its randomised klue_ml rows with
#                                         klue_screen TRUE: klue_ml refitted at
#                                         klue's default
#   R/statistical_tests.R                 its wilson_ci() only
# On the blocked ladder the stored seconds of Apollo's reduced search
# (apollo_ss_published) are stale: dev/rerun_apollo_blocked.R rewrote that
# arm's LLs but not its seconds. They are written with a note saying so, next
# to the seed-1 times, which the paper reports. dev/publish_filter.zsh
# publishes both logs. Without a log the script reads that ladder's
# apollo_lcEM rows from output/apollo_lcem_cap.csv, skips the log checks and
# leaves output/apollo_lcem_cap.csv as it is.
#
# Checks, made before any statistic is written (the script stops if one
# fails): best and gap recomputed from the LLs equal the stored ones; the gaps
# each log printed per cell, and its closing table of mean gaps and stuck
# counts, agree with the rds to the 2 dp printed; the diagnostic's seed-1 LLs
# are the ladder's; the isolation CSV's ladder rows carry the rds's LLs, its
# data check (chk, klue refitted on the regenerated data) reproduces the
# ladder's clustering+ML LL, its gaps are taken from the best of its runs and
# the ladder's, and its miss column equals gap < -0.5; every apollo_lcEM block
# follows the rule above; the carryover CSV has one klue_ml row with
# klue_screen TRUE per randomised cell, all from one klue version and with
# data_ok (the cell's data are the ladder's), its LL_stored is the rds's
# klue_ml LL to the 7 dp printed, and its reproduced column is
# |LL_now - LL_stored| <= 1e-3.
#
# Writes:
#   output/ladder_statistics.csv  long format, one row per statistic: ladder
#     (blocked, randomised, or isolation for the isolation runs on ten blocked
#     cells), set (a rung; hard+very_hard; all = every cell of the ladder or of
#     the isolation runs; or one cell such as hard|seed5), method (the arm as
#     named in the rds or the isolation CSV; best_of(klue_ml+klue_em+klue_em_rp)
#     is the best of the three klue arms, klue_em+klue_em_rp the two EM arms
#     together, klue_ml_screen_default clustering+ML at klue's default),
#     comparator (the other arm of a paired test, ratio or difference, or the
#     next best arm), threshold (the LL below the best that counts as a miss),
#     statistic, value, source (input file or files) and note. On the
#     randomised ladder, klue_ml_screen_default has moved_cells and, in each
#     cell where it moved, change_LL (LL_now - LL_stored), both with comparator
#     klue_ml; then its misses at 0.5 LL and its exact McNemar test against
#     apollo_searchStart at the recomputed bests. apollo_searchStart's misses
#     and test at those bests have comparator klue_ml_screen_default
#   output/apollo_lcem_cap.csv  one row per ladder cell: ladder, rung, seed,
#     capped, last_EM_iteration, stopped_by_criterion, continued_classical and
#     missed (apollo_lcEM more than 0.5 LL below the best), then cell,
#     last_improvement, gap_LL and log_line (the line of the "EM algorithm
#     stopped" message)
# and prints a summary to the console.
#
# Seeds: deterministic; no RNG. Base R only.
# Cost: about a second on one core; reads about 60 KB of rds and CSV files and
# 2.5 MB of logs.
# Run from the repository root:
#   nice -n 15 Rscript dev/ladder_statistics.R
# =============================================================================

if (!file.exists("dev/ladder_statistics.R")) stop("run from the repository root", call. = FALSE)

# ---- Inputs, outputs and definitions ----------------------------------------

LADDERS <- list(
  blocked    = list(rds = "output/stress_replicate_blocked.rds",
                    log = "output/stress_replicate_blocked_std.log"),
  randomised = list(rds = "output/stress_replicate.rds",
                    log = "output/stress_replicate_std.log"))
F_ISO     <- "output/apollo_isolation_blocked.csv"
F_DIAG    <- "output/apollo_smartstart_diag.csv"
F_CARRY   <- "output/ladder_carryover_check.csv"
F_WILSON  <- "R/statistical_tests.R"
OUT_STATS <- "output/ladder_statistics.csv"
OUT_CAP   <- "output/apollo_lcem_cap.csv"
for (f in c(F_ISO, F_DIAG, F_CARRY, F_WILSON, vapply(LADDERS, `[[`, "", "rds")))
  if (!file.exists(f)) stop("input not found: ", f, call. = FALSE)

MISS       <- 0.5                    # a miss: more than 0.5 LL below the best
THRESHOLDS <- c(0.1, 0.5, 1, 2, 5)   # the thresholds of the blocked table's note
PRINT_TOL  <- 0.005 + 1e-9           # the logs print gaps, the diagnostic CSV LLs, to 2 dp
CARRY_TOL  <- 5e-8 + 1e-9            # the carryover CSV prints LLs to 7 dp
MOVE_TOL   <- 1e-3                   # moved: |LL_now - LL_stored| > 1e-3, as in dev/check_ladder_carryover.R
EM_MAX_IT  <- 100L                   # EMmaxIterations of run_apollo_lcEM() (dev/compare_packages.R)
EM_CRIT    <- 1e-5                   # apollo_lcEM's default EMstoppingCriterion (Apollo 0.3.5)

RUNG_ORDER <- c("moderate", "hard", "very_hard", "hard_corr")
SETS <- list(`hard+very_hard` = c("hard", "very_hard"), all = RUNG_ORDER)  # cell sets beyond one rung
SET_TITLE <- c(`hard+very_hard` = "Hard and very-hard rungs", all = "All rungs")
# The ladder arms in the order of the paper's tables, with the paper's names.
LABEL <- c(klue_ml = "Clustering+ML", klue_em = "Clustering+EM", klue_em_rp = "RP+EM",
           gmnl_lc = "gmnl", apollo_searchStart = "Apollo's default search",
           apollo_ss_published = "Apollo's reduced search", apollo_lcEM = "Apollo's EM routine",
           apollo_clust = "Apollo's estimator from the clustering starts")
REF     <- "klue_ml"                               # the comparator of the paired tests
BEST3   <- c("klue_ml", "klue_em", "klue_em_rp")   # the three klue arms
B3_NAME <- "best_of(klue_ml+klue_em+klue_em_rp)"
EM2     <- c("klue_em", "klue_em_rp")              # the two EM arms, whose times the text pools
SCREEN    <- "klue_ml_screen_default"              # klue_ml refitted at klue's default (randomised ladder)
SCREEN_VS <- "apollo_searchStart"                  # the arm the table note tests it against
STALE   <- paste("stale: dev/rerun_apollo_blocked.R rewrote this arm's LLs but not its",
                 "seconds; the paper reports seed1_secs")
SEED1   <- paste("bounded_secs of dev/run_apollo_smartstart.R, the reduced search on seed 1",
                 "of the rung: the times the paper reports")
ISO_LABEL <- c(
  A0  = "one apollo_estimate fit from the pooled-MNL start",
  A4  = "one apollo_estimate fit from the generic start",
  A1  = "reduced budget within +-0.1 of the pooled-MNL start",
  A3  = "reduced budget from the pooled-MNL start, candidates on [-3, 3]",
  A2  = "Apollo's default search from the generic start",
  A5  = "Apollo's default search from the pooled-MNL start, analytic gradients",
  A6  = "Apollo's estimator from the clustering starts, one covariance matrix",
  chk = "clustering+ML refitted on the regenerated data (data check)")
ISO_TESTS <- list(c("A3", "A1"), c("A4", "A2"), c("A4", "klue_ml"))   # arm, comparator

# wilson_ci() of R/statistical_tests.R. Sourcing that file runs its analyses,
# so only the function's definition is evaluated.
local({
  def <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                  identical(e[[2]], as.name("wilson_ci")),
                as.list(parse(F_WILSON, keep.source = FALSE)))
  if (length(def) != 1L) stop(F_WILSON, ": expected one definition of wilson_ci", call. = FALSE)
  eval(def[[1]], envir = globalenv())
})

# ---- Helpers ----------------------------------------------------------------

ROWS <- list()
emit <- function(ladder, set, method, statistic, value, source,
                 comparator = NA_character_, threshold = NA_real_, note = "")
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    ladder = ladder, set = set, method = method, comparator = comparator,
    threshold = threshold, statistic = statistic, value = unname(as.numeric(value)),
    source = source, note = note, stringsAsFactors = FALSE)

# Exact McNemar test of two paired miss vectors: two-sided exact binomial test
# on the cells missed by only one of them.
mcnemar_exact <- function(a, b) {
  n_a <- sum(a & !b); n_b <- sum(!a & b)
  c(only_method = n_a, only_comparator = n_b, both = sum(a & b), neither = sum(!a & !b),
    p_exact = if (n_a + n_b > 0L) binom.test(n_a, n_a + n_b, 0.5)$p.value else 1)
}

f2   <- function(v) sub("^-(0\\.00)$", "\\1", sprintf("%.2f", v))   # mean gaps as the tables print them
pad  <- function(s, w) formatC(s, width = -w)                        # left-justified
lab  <- function(a) unname(ifelse(a %in% names(LABEL), LABEL[a], a))
cells_txt <- function(v) if (length(v)) paste(v, collapse = ", ") else "none"
secs_txt <- function(lo, hi) {                   # seconds, and minutes from two minutes up
  s <- sprintf("%.1f-%.1f s", lo, hi)
  if (hi >= 120) s <- sprintf("%s (%.1f-%.1f min)", s, lo / 60, hi / 60)
  s
}
fail <- function(...) stop(..., call. = FALSE)

# ---- A ladder: the rds, checked against the runner's definitions ------------

read_ladder <- function(lad) {
  f <- LADDERS[[lad]]$rds
  x <- readRDS(f)
  arms <- names(x[[1]]$LL)
  ok <- vapply(names(x), function(cl) {
    z <- x[[cl]]
    identical(names(z$LL), arms) && identical(names(z$gap), arms) &&
      identical(names(z$secs), arms) && identical(z$harness, "apollo_std") &&
      identical(cl, paste0(z$rung, "|seed", z$seed)) &&
      all(is.finite(z$LL)) && all(is.finite(z$secs))
  }, logical(1))
  if (!all(ok)) fail(f, ": cells ", paste(names(x)[!ok], collapse = ", "), " lack an arm, ",
                     "a finite LL or time, a matching name or the harness tag apollo_std")
  if (!all(arms %in% names(LABEL))) fail(f, ": arm without a label: ",
                                         paste(setdiff(arms, names(LABEL)), collapse = ", "))
  rung  <- vapply(x, `[[`, "", "rung")
  rungs <- intersect(RUNG_ORDER, rung)
  seeds <- sort(unique(vapply(x, function(z) as.integer(z$seed), 1L)))
  cells <- paste0(rep(rungs, each = length(seeds)), "|seed", seeds)
  if (!setequal(cells, names(x)) || length(cells) != length(x))
    fail(f, ": the cells are not a complete grid of rungs ", paste(RUNG_ORDER, collapse = ", "),
         " x seeds")
  x <- x[cells]
  arms <- intersect(names(LABEL), arms)
  cells_x_arms <- function(field) {
    m <- t(vapply(x, function(z) unname(z[[field]][arms]), numeric(length(arms))))
    dimnames(m) <- list(cells, arms)
    m
  }
  LL   <- cells_x_arms("LL")
  best <- apply(LL, 1, max)
  GAP  <- LL - best
  d    <- max(abs(best - vapply(x, `[[`, 0, "best")), abs(GAP - cells_x_arms("gap")))
  if (d > 1e-9) fail(f, ": the stored best or gap differ from the recomputed ones by ", d)
  list(lad = lad, file = f, cells = cells, arms = arms, rungs = rungs, seeds = seeds,
       rung = setNames(rep(rungs, each = length(seeds)), cells),
       seed = setNames(rep(seeds, length(rungs)), cells),
       LL = LL, GAP = GAP, SEC = cells_x_arms("secs"), d_stored = d)
}

# ---- The runner's log: per-cell gaps, closing table, apollo_lcEM blocks -----

RE_RUN  <- "^\\[[0-9:]+\\] RUN  (\\S+)  .*$"
RE_DONE <- "^\\[[0-9:]+\\] +done (\\S+) : gaps +(.*)$"

# The runner prints each cell's gaps when the cell is done, and at the end the
# mean gap and stuck count per rung and arm; both must agree with the rds.
check_log <- function(D, lg, f) {
  i_done <- grep(RE_DONE, lg)
  cl_done <- sub(RE_DONE, "\\1", lg[i_done])
  if (anyDuplicated(cl_done) || !setequal(cl_done, D$cells))
    fail(f, ": expected one 'done' line per cell of ", D$file)
  worst_cell <- 0
  for (k in seq_along(i_done)) {
    tok <- strsplit(sub(RE_DONE, "\\2", lg[i_done[k]]), " +")[[1]]
    g <- setNames(as.numeric(sub("^[^=]*=", "", tok)), sub("=.*$", "", tok))
    if (!setequal(names(g), D$arms)) fail(f, ", line ", i_done[k], ": arms differ from the rds")
    worst_cell <- max(worst_cell, abs(g[D$arms] - D$GAP[cl_done[k], D$arms]))
  }
  h <- grep("MEAN gap_to_best by rung x method", lg, fixed = TRUE)
  if (length(h) != 1L) fail(f, ": expected one closing table of mean gaps")
  arms_log <- trimws(strsplit(lg[h + 1L], "|", fixed = TRUE)[[1]])[-1]
  if (!setequal(arms_log, D$arms)) fail(f, ": the closing table's arms differ from the rds")
  worst_mean <- 0; bad <- 0L; seen <- character(0)
  for (ln in lg[h + 1L + seq_along(D$rungs)]) {
    fld <- trimws(strsplit(ln, "|", fixed = TRUE)[[1]])
    if (!fld[1] %in% D$rungs) fail(f, ": unexpected row of the closing table: ", ln)
    seen <- c(seen, fld[1])
    for (j in seq_along(arms_log)) {
      m <- regmatches(fld[j + 1L], regexec("^(-?[0-9.]+) \\(stuck ([0-9]+)/([0-9]+)\\)$",
                                           fld[j + 1L]))[[1]]
      if (length(m) != 4L) fail(f, ": cannot read the entry '", fld[j + 1L], "'")
      g <- D$GAP[D$rung == fld[1], arms_log[j]]
      worst_mean <- max(worst_mean, abs(as.numeric(m[2]) - mean(g)))
      bad <- bad + (as.integer(m[3]) != sum(g < -MISS)) + (as.integer(m[4]) != length(g))
    }
  }
  if (!setequal(seen, D$rungs)) fail(f, ": the closing table lacks a rung")
  if (worst_cell > PRINT_TOL || worst_mean > PRINT_TOL || bad > 0L)
    fail(f, " disagrees with ", D$file, ": per-cell gaps by up to ", signif(worst_cell, 3),
         ", mean gaps by up to ", signif(worst_mean, 3), ", ", bad, " stuck counts")
  c(cells = length(i_done), entries = length(seen) * length(arms_log),
    worst = max(worst_cell, worst_mean))
}

# One row per cell: apollo_lcEM's EM block in that cell's section of the log.
parse_lcem <- function(D, lg, f) {
  run_cl  <- ifelse(grepl(RE_RUN, lg), sub(RE_RUN, "\\1", lg), NA_character_)
  done_cl <- ifelse(grepl(RE_DONE, lg), sub(RE_DONE, "\\1", lg), NA_character_)
  rows <- lapply(D$cells, function(cl) {
    i <- which(run_cl %in% cl); j <- which(done_cl %in% cl)
    if (length(i) != 1L || length(j) != 1L || j < i || any(!is.na(run_cl[seq.int(i + 1L, j)])))
      fail(f, ", ", cl, ": expected one section from 'RUN' to 'done'")
    sec  <- i:j
    init <- sec[startsWith(lg[sec], "Initialising EM algorithm")]
    stp  <- sec[startsWith(lg[sec], "EM algorithm stopped: ")]
    if (length(init) != 1L || length(stp) != 1L || stp < init)
      fail(f, ", ", cl, ": expected one apollo_lcEM block")
    em  <- lg[init:stp]
    it  <- as.integer(sub("^Starting iteration: ", "", em[startsWith(em, "Starting iteration: ")]))
    imp <- as.numeric(sub("^- Improvement: ", "", em[startsWith(em, "- Improvement: ")]))
    if (!length(it) || !identical(it, seq_along(it)) || length(imp) != length(it))
      fail(f, ", ", cl, ": the EM iterations are not numbered 1, 2, ... with one improvement each")
    data.frame(ladder = D$lad, rung = D$rung[[cl]], seed = D$seed[[cl]],
               capped = grepl("maximum number of iterations reached", lg[stp], fixed = TRUE),
               last_EM_iteration = length(it),
               stopped_by_criterion = imp[length(imp)] < EM_CRIT,
               continued_classical = any(startsWith(lg[stp:j], "Continuing with classical estimation")),
               missed = D$GAP[cl, "apollo_lcEM"] < -MISS,
               cell = cl, last_improvement = imp[length(imp)],
               gap_LL = D$GAP[cl, "apollo_lcEM"], log_line = stp,
               criterion_message = grepl("smaller than convergence", lg[stp], fixed = TRUE),
               stringsAsFactors = FALSE)
  })
  cap <- do.call(rbind, rows)
  # Apollo 0.3.5 (apollo_lcEM): the loop ends when the improvement falls below
  # the criterion or after iteration EM_MAX_IT; the counter is then one past
  # the last iteration, and a counter >= EM_MAX_IT means the cap message and
  # no closing ML step.
  ok <- cap$capped == (cap$last_EM_iteration >= EM_MAX_IT - 1L) &
        cap$capped != cap$criterion_message &
        (cap$last_EM_iteration == EM_MAX_IT | cap$stopped_by_criterion) &
        (cap$capped | cap$stopped_by_criterion) &
        cap$capped != cap$continued_classical
  if (!all(ok)) fail(f, ": apollo_lcEM in ", paste(cap$cell[!ok], collapse = ", "),
                     " does not follow Apollo 0.3.5's stopping rule")
  cap$criterion_message <- NULL
  cap
}

# Without the log: the rows an earlier run wrote to OUT_CAP, checked against the rds.
stored_cap <- function(D) {
  if (!file.exists(OUT_CAP))
    fail(LADDERS[[D$lad]]$log, " and ", OUT_CAP, " not found: no source for the apollo_lcEM cap")
  cap <- read.csv(OUT_CAP, stringsAsFactors = FALSE)
  cap <- cap[cap$ladder == D$lad, ]
  cap <- cap[match(D$cells, cap$cell), ]
  if (anyNA(cap$cell) || !identical(cap$missed, unname(D$GAP[D$cells, "apollo_lcEM"] < -MISS)))
    fail(OUT_CAP, ": the ", D$lad, " rows do not match ", D$file)
  rownames(cap) <- NULL
  cap
}

# ---- One ladder: statistics and console summary ----------------------------

ladder_section <- function(lad) {
  D   <- read_ladder(lad)
  src <- D$file
  f_log <- LADDERS[[lad]]$log
  have_log <- file.exists(f_log)
  if (have_log) {
    lg  <- readLines(f_log, warn = FALSE)
    chk <- check_log(D, lg, f_log)
    cap <- parse_lcem(D, lg, f_log)
  } else cap <- stored_cap(D)
  cap_src <- if (have_log) f_log else OUT_CAP
  miss  <- D$GAP < -MISS                                     # cells x arms
  RM    <- sapply(D$arms, function(a) tapply(D$SEC[, a], D$rung, mean)[D$rungs])  # rungs x arms
  stale <- function(a) if (lad == "blocked" && a == "apollo_ss_published") STALE else ""
  w <- max(nchar(LABEL[D$arms]))

  cat(sprintf("\n==== %s ladder: %s ====\n", if (lad == "blocked") "Blocked" else "Randomised", src))
  cat(sprintf("%d cells (%s; seeds %d-%d), %d arms\n", length(D$cells), paste(D$rungs, collapse = ", "),
              min(D$seeds), max(D$seeds), length(D$arms)))
  cat(sprintf("Checks passed: best and gaps reproduce the stored ones (max |diff| %.1e)", D$d_stored))
  if (have_log) cat(sprintf(";\n  the log's %d per-cell gap lines and its closing table (%d entries) agree with\n  the rds to the 2 dp printed (max |diff| %.4f).\n",
                            chk[["cells"]], chk[["entries"]], chk[["worst"]]))
  else cat(sprintf(".\n  %s not found: log checks skipped; apollo_lcEM rows read from %s.\n", f_log, OUT_CAP))

  # -- per rung: the table's cells, threshold counts, time per fit, ratios
  for (rg in D$rungs) {
    i <- D$rung == rg
    emit(lad, rg, NA, "n_cells", sum(i), src)
    for (a in D$arms) {
      g <- D$GAP[i, a]
      emit(lad, rg, a, "mean_gap_LL", mean(g), src)
      for (t in THRESHOLDS) emit(lad, rg, a, "misses", sum(g < -t), src, threshold = t)
      emit(lad, rg, a, "mean_secs", RM[rg, a], src, note = stale(a))
      if (a != REF)
        emit(lad, rg, a, "mean_secs_ratio", RM[rg, a] / RM[rg, REF], src, comparator = REF,
             note = paste(c("the arm's mean seconds over the comparator's", if (nzchar(stale(a))) "stale"),
                          collapse = "; "))
    }
  }
  for (a in D$arms) {
    emit(lad, "all", a, "rung_mean_secs_min", min(RM[, a]), src, note = stale(a))
    emit(lad, "all", a, "rung_mean_secs_max", max(RM[, a]), src, note = stale(a))
    emit(lad, "all", a, "secs_min", min(D$SEC[, a]), src, note = stale(a))
    emit(lad, "all", a, "secs_median", median(D$SEC[, a]), src, note = stale(a))
    emit(lad, "all", a, "secs_max", max(D$SEC[, a]), src, note = stale(a))
  }
  em2 <- all(EM2 %in% D$arms)
  if (em2) {
    emit(lad, "all", paste(EM2, collapse = "+"), "rung_mean_secs_min", min(RM[, EM2]), src)
    emit(lad, "all", paste(EM2, collapse = "+"), "rung_mean_secs_max", max(RM[, EM2]), src)
  }
  s1 <- NULL
  if (lad == "blocked") {         # Apollo's reduced search timed on seed 1 of each rung
    dg <- read.csv(F_DIAG, stringsAsFactors = FALSE)
    dg <- dg[match(D$rungs, dg$rung), ]
    if (anyNA(dg$rung)) fail(F_DIAG, ": a rung of the blocked ladder is missing")
    # the diagnostic reran the arm on the ladder's seed-1 data: same LL (2 dp)
    d1 <- max(abs(dg$bounded_LL - D$LL[paste0(D$rungs, "|seed1"), "apollo_ss_published"]))
    if (d1 > PRINT_TOL) fail(F_DIAG, ": bounded_LL is not the ladder's seed-1 LL (|diff| ", d1, ")")
    s1 <- setNames(dg$bounded_secs, D$rungs)
    for (rg in D$rungs) emit(lad, rg, "apollo_ss_published", "seed1_secs", s1[[rg]], F_DIAG, note = SEED1)
    emit(lad, "all", "apollo_ss_published", "seed1_secs_min", min(s1), F_DIAG, note = SEED1)
    emit(lad, "all", "apollo_ss_published", "seed1_secs_max", max(s1), F_DIAG, note = SEED1)
  }

  cat("\nMean gap to the best LL by rung, misses (> 0.5 LL below) in parentheses; time per\nfit = range of the per-rung mean times\n")
  cat(sprintf("  %s %s  time per fit\n", pad("", w), paste(formatC(D$rungs, width = 11), collapse = " ")))
  for (a in D$arms) {
    tc <- vapply(D$rungs, function(rg) sprintf("%s (%d)", f2(mean(D$GAP[D$rung == rg, a])),
                                               sum(miss[D$rung == rg, a])), "")
    tt <- secs_txt(min(RM[, a]), max(RM[, a]))
    if (!is.null(s1) && a == "apollo_ss_published")
      tt <- sprintf("%s on seed 1 [stored: %s, stale]", secs_txt(min(s1), max(s1)), tt)
    cat(sprintf("  %s %s  %s\n", pad(LABEL[[a]], w), paste(formatC(tc, width = 11), collapse = " "), tt))
  }
  if (em2) cat(sprintf("  The two EM arms together: %s\n", secs_txt(min(RM[, EM2]), max(RM[, EM2]))))
  cat(sprintf("Per-rung mean time relative to clustering+ML (%s):\n", paste(D$rungs, collapse = " / ")))
  for (a in setdiff(D$arms, REF))
    cat(sprintf("  %s %s%s\n", pad(LABEL[[a]], w),
                paste(sprintf("%7.1f", RM[, a] / RM[, REF]), collapse = " "),
                if (nzchar(stale(a))) "  [stale]" else ""))

  # -- cell sets: threshold counts, reach, Wilson, miss sizes, McNemar tests
  for (sn in names(SETS)) {
    i <- D$rung %in% SETS[[sn]]; n <- sum(i)
    emit(lad, sn, NA, "n_cells", n, src)
    cat(sprintf("\n%s, %d cells: misses at %s LL below the best; exact McNemar test\nagainst clustering+ML at 0.5 LL (cells missed only by the arm / only by clustering+ML / by both)\n",
                SET_TITLE[[sn]], n, paste(THRESHOLDS, collapse = " / ")))
    cat(sprintf("  %s %s   %s\n", pad("", w), paste(formatC(THRESHOLDS, width = 4), collapse = " "),
                "only arm / only C+ML / both       p"))
    for (a in D$arms) {
      g <- D$GAP[i, a]; k <- sum(g < -MISS)
      for (t in THRESHOLDS) emit(lad, sn, a, "misses", sum(g < -t), src, threshold = t)
      emit(lad, sn, a, "reached", n - k, src, threshold = MISS)
      ci <- wilson_ci(k, n)
      emit(lad, sn, a, "miss_rate_wilson95_lower", ci[["lower"]], src, threshold = MISS)
      emit(lad, sn, a, "miss_rate_wilson95_upper", ci[["upper"]], src, threshold = MISS)
      if (k > 0L) {
        sh <- -g[g < -MISS]
        emit(lad, sn, a, "shortfall_min_LL", min(sh), src, threshold = MISS)
        emit(lad, sn, a, "shortfall_max_LL", max(sh), src, threshold = MISS)
        emit(lad, sn, a, "bic_equivalent_min", 2 * min(sh), src, threshold = MISS)
        emit(lad, sn, a, "bic_equivalent_max", 2 * max(sh), src, threshold = MISS)
      }
      tst <- ""
      if (a != REF) {
        tt <- mcnemar_exact(miss[i, a], miss[i, REF])
        for (st in names(tt))
          emit(lad, sn, a, paste0("mcnemar_", st), tt[[st]], src, comparator = REF, threshold = MISS)
        tst <- sprintf("%8d / %9d / %4d  %8.4f", tt[["only_method"]], tt[["only_comparator"]],
                       tt[["both"]], tt[["p_exact"]])
      }
      cat(sprintf("  %s %s   %s\n", pad(LABEL[[a]], w),
                  paste(formatC(vapply(THRESHOLDS, function(t) sum(g < -t), 0L), width = 4), collapse = " "),
                  tst))
    }
    if (all(BEST3 %in% D$arms)) {
      b3 <- apply(D$GAP[i, BEST3], 1, max)
      emit(lad, sn, B3_NAME, "misses", sum(b3 < -MISS), src, threshold = MISS)
      cat(sprintf("  Best of clustering+ML, clustering+EM and RP+EM: misses %d of %d%s\n", sum(b3 < -MISS), n,
                  if (any(b3 < -MISS)) paste0(" (", paste(sprintf("%s by %.2f LL", names(b3)[b3 < -MISS],
                                                                  -b3[b3 < -MISS]), collapse = ", "), ")") else ""))
    }
    k <- sum(miss[i, REF]); ci <- wilson_ci(k, n)
    cat(sprintf("  Clustering+ML misses %d of %d: Wilson 95%% interval %.1f%% to %.1f%%\n",
                k, n, 100 * ci[["lower"]], 100 * ci[["upper"]]))
    cat(sprintf("  Reached the best (within 0.5 LL), of %d: %s\n", n,
                paste(sprintf("%s %d", lab(D$arms), n - colSums(miss[i, , drop = FALSE])), collapse = ", ")))
  }

  # -- every miss, with its BIC equivalent
  cat("\nMisses (all cells): LL below the best [BIC equivalent, 2 x LL]\n")
  for (a in D$arms) {
    mc <- D$cells[miss[, a]]
    for (cl in mc) {
      emit(lad, cl, a, "shortfall_LL", -D$GAP[cl, a], src, threshold = MISS)
      emit(lad, cl, a, "bic_equivalent", -2 * D$GAP[cl, a], src, threshold = MISS)
    }
    cat(sprintf("  %s %s\n", pad(LABEL[[a]], w), if (length(mc))
      paste(sprintf("%s %.2f [%.1f]", mc, -D$GAP[mc, a], -2 * D$GAP[mc, a]), collapse = ", ") else "none"))
  }
  if (all(BEST3 %in% D$arms)) {
    b3 <- apply(D$GAP[, BEST3], 1, max)
    for (cl in D$cells[b3 < -MISS]) {
      emit(lad, cl, B3_NAME, "shortfall_LL", -b3[[cl]], src, threshold = MISS)
      emit(lad, cl, B3_NAME, "bic_equivalent", -2 * b3[[cl]], src, threshold = MISS)
    }
  }

  # -- cells where one arm alone is at the best (every other arm misses)
  sole <- data.frame(cell = character(0), arm = character(0), next_arm = character(0),
                     margin = numeric(0), stringsAsFactors = FALSE)
  for (cl in D$cells) {
    at <- D$arms[!miss[cl, ]]
    if (length(at) == 1L) {
      o <- sort(D$LL[cl, ], decreasing = TRUE)
      sole[nrow(sole) + 1L, ] <- list(cl, at, names(o)[2], o[[1]] - o[[2]])
      emit(lad, cl, at, "sole_best_margin_LL", o[[1]] - o[[2]], src, comparator = names(o)[2],
           threshold = MISS, note = "LL above the next best arm (comparator); every other arm misses")
    }
  }
  for (a in D$arms) {
    s <- sole[sole$arm == a, ]
    emit(lad, "all", a, "sole_best_cells", nrow(s), src, threshold = MISS)
    if (nrow(s)) {
      emit(lad, "all", a, "sole_best_margin_min_LL", min(s$margin), src, threshold = MISS)
      emit(lad, "all", a, "sole_best_margin_max_LL", max(s$margin), src, threshold = MISS)
    }
  }
  cat("\nCells where one arm alone is at the best (every other arm misses):\n")
  if (!nrow(sole)) cat("  none\n")
  for (a in unique(sole$arm)) {
    s <- sole[sole$arm == a, ]
    cat(sprintf("  %s, %d: %s%s\n", LABEL[[a]], nrow(s),
                paste(sprintf("%s (+%.2f LL over %s)", s$cell, s$margin, lab(s$next_arm)), collapse = ", "),
                if (nrow(s) > 1L) sprintf("; margins %.2f to %.2f LL", min(s$margin), max(s$margin)) else ""))
    if (a != REF) {
      mk <- D$cells[miss[, REF]]
      cat(sprintf("    these are all of clustering+ML's misses: %s\n", if (setequal(s$cell, mk)) "yes"
                  else paste("no; clustering+ML misses", paste(mk, collapse = ", "))))
    }
  }

  cat("\nSeconds per fit over all cells, min / median / max:\n")
  for (a in D$arms) {
    s <- D$SEC[, a]
    cat(sprintf("  %s %.1f / %.1f / %.1f s%s%s\n", pad(LABEL[[a]], w), min(s), median(s), max(s),
                if (max(s) >= 120) sprintf(" (%.1f / %.1f / %.1f min)", min(s) / 60, median(s) / 60,
                                           max(s) / 60) else "",
                if (nzchar(stale(a))) "  [stale]" else ""))
  }

  # -- apollo_lcEM: the iteration cap and the closing ML step
  for (sn in c(D$rungs, names(SETS))) {
    i <- cap$rung %in% (if (sn %in% names(SETS)) SETS[[sn]] else sn)
    emit(lad, sn, "apollo_lcEM", "lcem_capped_cells", sum(cap$capped[i]), cap_src)
    emit(lad, sn, "apollo_lcEM", "lcem_capped_criterion_met_cells",
         sum(cap$capped[i] & cap$stopped_by_criterion[i]), cap_src,
         note = "reported the cap although the last improvement was below 1e-5")
    emit(lad, sn, "apollo_lcEM", "lcem_continued_classical_cells", sum(cap$continued_classical[i]), cap_src)
    emit(lad, sn, "apollo_lcEM", "lcem_missed_capped_cells", sum(cap$missed[i] & cap$capped[i]),
         paste(cap_src, "+", src), threshold = MISS)
  }
  cat(sprintf("\napollo_lcEM (%s): capped in %d of %d cells: %s\n", cap_src, sum(cap$capped), nrow(cap),
              paste(cap$cell[cap$capped], collapse = ", ")))
  cat(sprintf("  continued with classical estimation in the other %d\n", sum(cap$continued_classical)))
  for (sn in names(SETS)) {
    i <- cap$rung %in% SETS[[sn]]
    cat(sprintf("  %s: %d misses, %d of them capped; misses not capped: %s\n", SET_TITLE[[sn]],
                sum(cap$missed[i]), sum(cap$missed[i] & cap$capped[i]),
                cells_txt(cap$cell[i & cap$missed & !cap$capped])))
  }
  sc <- cap[cap$capped & cap$stopped_by_criterion, ]
  cat(sprintf("  capped although the last improvement was below %g: %s\n", EM_CRIT,
              cells_txt(sprintf("%s (iteration %d, improvement %.3g, log line %d)", sc$cell,
                                sc$last_EM_iteration, sc$last_improvement, sc$log_line))))
  list(D = D, cap = cap, have_log = have_log)
}

# ---- Randomised ladder: clustering+ML at klue's default ---------------------

# The refit at klue's default (dev/check_ladder_carryover.R) replaces the stored
# klue_ml LL as the arm SCREEN; each cell's best is recomputed over the arms,
# and SCREEN and SCREEN_VS are scored and tested at those bests.
screen_section <- function(D) {
  cc <- read.csv(F_CARRY, stringsAsFactors = FALSE)
  cc <- cc[cc$ladder == D$lad & cc$method == REF & cc$klue_screen %in% TRUE, ]
  cc_cell <- paste0(cc$rung, "|seed", cc$seed)
  if (anyDuplicated(cc_cell) || !setequal(cc_cell, D$cells))
    fail(F_CARRY, ": expected one ", REF, " row with klue_screen TRUE per cell of ", D$file)
  cc  <- cc[match(D$cells, cc_cell), ]
  ver <- unique(cc$klue_version)
  if (length(ver) != 1L || !all(is.finite(cc$LL_now)) || !all(cc$data_ok %in% TRUE))
    fail(F_CARRY, ": the ", D$lad, " ", REF, " rows with klue_screen TRUE need one klue version, ",
         "a finite LL_now and data_ok TRUE")
  d_st <- max(abs(cc$LL_stored - D$LL[, REF]))
  if (d_st > CARRY_TOL) fail(F_CARRY, ": LL_stored is not the ", REF, " LL of ", D$file, " (|diff| ", d_st, ")")
  chg   <- setNames(cc$LL_now - cc$LL_stored, D$cells)
  moved <- abs(chg) > MOVE_TOL
  if (!identical(cc$reproduced, unname(!moved)))
    fail(F_CARRY, ": its reproduced column is not |LL_now - LL_stored| <= ", MOVE_TOL)
  LL <- D$LL
  LL[, REF] <- cc$LL_now
  colnames(LL)[colnames(LL) == REF] <- SCREEN
  miss <- LL - apply(LL, 1, max) < -MISS                     # cells x arms, at the recomputed bests

  srcb <- paste(D$file, "+", F_CARRY)
  def  <- sprintf("klue_ml refitted at klue %s's default (klue.screen = TRUE)", ver)
  nb   <- sprintf("each cell's best recomputed with %s, %s, in place of the stored klue_ml", SCREEN, def)
  for (cl in D$cells[moved])
    emit(D$lad, cl, SCREEN, "change_LL", chg[[cl]], F_CARRY, comparator = REF,
         note = paste("LL of", def, "minus the stored LL (comparator)"))
  emit(D$lad, "all", SCREEN, "moved_cells", sum(moved), F_CARRY, comparator = REF,
       note = sprintf("cells where %s ends more than %g LL from the stored LL (comparator)", def, MOVE_TOL))
  emit(D$lad, "all", SCREEN, "misses", sum(miss[, SCREEN]), srcb, threshold = MISS, note = nb)
  emit(D$lad, "all", SCREEN_VS, "misses", sum(miss[, SCREEN_VS]), srcb, comparator = SCREEN,
       threshold = MISS, note = nb)
  for (pr in list(c(SCREEN, SCREEN_VS), c(SCREEN_VS, SCREEN))) {
    tt <- mcnemar_exact(miss[, pr[1]], miss[, pr[2]])
    for (st in names(tt))
      emit(D$lad, "all", pr[1], paste0("mcnemar_", st), tt[[st]], srcb, comparator = pr[2],
           threshold = MISS, note = nb)
  }

  tt <- mcnemar_exact(miss[, SCREEN], miss[, SCREEN_VS])            # as the table note states it
  cat(sprintf("\n==== Randomised ladder at klue %s's default: %s ====\n", ver, F_CARRY))
  cat(sprintf("Checks passed: one %s row with klue_screen TRUE per cell, each with data_ok; LL_stored\n  is the rds's (max |diff| %.1e); reproduced = |LL_now - LL_stored| <= %g.\n",
              REF, d_st, MOVE_TOL))
  cat(sprintf("Clustering+ML refitted at the default ends more than %g LL from the stored LL in %d of %d\n  cells (LL_now - LL_stored): %s\n",
              MOVE_TOL, sum(moved), length(moved), cells_txt(sprintf("%s %+.4f", D$cells[moved], chg[moved]))))
  cat(sprintf("Misses at the recomputed bests, with those at the stored LLs in parentheses:\n  clustering+ML %d (%d), Apollo's default search %d (%d)\n",
              sum(miss[, SCREEN]), sum(D$GAP[, REF] < -MISS), sum(miss[, SCREEN_VS]), sum(D$GAP[, SCREEN_VS] < -MISS)))
  cat(sprintf("Exact McNemar test against Apollo's default search: p = %.4f\n  missed only by clustering+ML (%d): %s\n  only by Apollo's default search (%d): %s\n  by both: %d\n",
              tt[["p_exact"]], tt[["only_method"]], cells_txt(D$cells[miss[, SCREEN] & !miss[, SCREEN_VS]]),
              tt[["only_comparator"]], cells_txt(D$cells[miss[, SCREEN_VS] & !miss[, SCREEN]]), tt[["both"]]))
}

# ---- Isolation runs (ten blocked cells) -------------------------------------

iso_section <- function(D) {
  iso <- read.csv(F_ISO, stringsAsFactors = FALSE)
  iso$arm <- sub("^ladder:", "", iso$arm)
  cells <- unique(iso$cell); arms <- unique(iso$arm)
  if (nrow(iso) != length(cells) * length(arms) || anyDuplicated(iso[c("cell", "arm")]) ||
      !all(cells %in% D$cells) || !all(is.finite(iso$LL)) || !all(is.finite(iso$secs)) ||
      !all(c(unlist(ISO_TESTS), "apollo_searchStart", "A5", "chk") %in% arms))
    fail(F_ISO, ": expected one finite row per cell and arm, with the arms of the tests")
  if (!identical(iso$miss, iso$gap < -MISS)) fail(F_ISO, ": its miss column disagrees with its gap")
  # the ladder rows are the ladder's fits (LL rounded to 3 dp) ...
  lr <- iso$arm %in% D$arms
  d_ll <- max(abs(iso$LL[lr] - D$LL[cbind(iso$cell[lr], iso$arm[lr])]))
  if (d_ll > 5e-4 + 1e-9) fail(F_ISO, ": the ladder rows do not carry the LLs of ", D$file)
  # ... and the gaps are taken from the best of the isolation runs and the ladder
  best_iso <- pmax(tapply(iso$LL[!lr], iso$cell[!lr], max)[cells], apply(D$LL[cells, ], 1, max))
  d_best <- max(abs(iso$LL - iso$gap - best_iso[iso$cell]))   # three roundings to 3 dp at most
  if (d_best > 1.5e-3 + 1e-9) fail(F_ISO, ": its gaps are not taken from the best of its runs and the ladder")
  by_cell <- function(v, a) { d <- iso[iso$arm == a, ]; setNames(d[[v]], d$cell)[cells] }
  M <- sapply(arms, by_cell, v = "miss")
  S <- sapply(arms, by_cell, v = "secs")
  # the data check of the isolation runs: klue refitted on the regenerated data
  # (chk) reaches the ladder's clustering+ML LL in every cell
  d_chk <- max(abs(by_cell("LL", "chk") - by_cell("LL", "klue_ml")))
  if (d_chk > 1e-3) fail(F_ISO, ": chk does not reproduce the ladder's klue_ml LL (|diff| ", d_chk, ")")
  desc <- function(a) if (a %in% names(ISO_LABEL)) ISO_LABEL[[a]] else paste("ladder run:", LABEL[[a]])
  isrc <- F_ISO
  stale <- function(a) if (a == "apollo_ss_published") STALE else ""

  iso_seeds <- sort(unique(D$seed[cells]))
  cat(sprintf("\n==== Isolation runs: %s ====\n", F_ISO))
  cat(sprintf("%d blocked cells (%s; seeds %s), %d arms. A miss = more than 0.5 LL below\nthe best LL of any run in the file or on the ladder.\n",
              length(cells), paste(unique(D$rung[cells]), collapse = " and "),
              if (all(diff(iso_seeds) == 1L)) sprintf("%d-%d", min(iso_seeds), max(iso_seeds))
              else paste(iso_seeds, collapse = ", "),
              length(arms)))
  cat(sprintf("Checks passed: the ladder rows carry the rds's LLs (max |diff| %.1e); chk reproduces\n  klue_ml (max |diff| %.1e); the gaps are taken from the best of the runs and the ladder\n  (max |diff| %.1e); miss = gap < -0.5.\n",
              d_ll, d_chk, d_best))
  emit("isolation", "all", NA, "n_cells", length(cells), isrc)
  w <- max(nchar(arms))
  cat(sprintf("  %s %s  %s\n", pad("arm", w), pad("misses", 7), "time per run: median (min-max)"))
  for (a in arms) {
    s <- S[, a]
    emit("isolation", "all", a, "misses", sum(M[, a]), isrc, threshold = MISS, note = desc(a))
    emit("isolation", "all", a, "secs_min", min(s), isrc, note = stale(a))
    emit("isolation", "all", a, "secs_median", median(s), isrc, note = stale(a))
    emit("isolation", "all", a, "secs_max", max(s), isrc, note = stale(a))
    tm <- if (median(s) >= 120) sprintf("%.1f min (%.1f-%.1f)", median(s) / 60, min(s) / 60, max(s) / 60)
          else sprintf("%.1f s (%.1f-%.1f)", median(s), min(s), max(s))
    cat(sprintf("  %s %s  %s %s%s\n", pad(a, w), pad(sprintf("%d/%d", sum(M[, a]), length(cells)), 7),
                pad(tm, 22), desc(a), if (nzchar(stale(a))) " [time stale]" else ""))
  }
  cat("Paired exact tests at 0.5 LL (cells missed only by the arm / only by the comparator):\n")
  for (pr in ISO_TESTS) {
    tt <- mcnemar_exact(M[, pr[1]], M[, pr[2]])
    for (st in names(tt))
      emit("isolation", "all", pr[1], paste0("mcnemar_", st), tt[[st]], isrc,
           comparator = pr[2], threshold = MISS)
    cat(sprintf("  %s against %s: %d / %d, p = %.4f\n", pr[1], lab(pr[2]),
                tt[["only_method"]], tt[["only_comparator"]], tt[["p_exact"]]))
  }
  r <- S[, "apollo_searchStart"] / S[, "A5"]
  rn <- "the comparator's seconds over the arm's, on the same cell"
  for (cl in cells) emit("isolation", cl, "A5", "speedup", r[[cl]], isrc, comparator = "apollo_searchStart", note = rn)
  emit("isolation", "all", "A5", "speedup_min", min(r), isrc, comparator = "apollo_searchStart", note = rn)
  emit("isolation", "all", "A5", "speedup_max", max(r), isrc, comparator = "apollo_searchStart", note = rn)
  cat(sprintf("A5 against the ladder's default search: %.1f-%.1f min against %.1f-%.1f min per fit,\n  %.2f to %.2f times faster, cell by cell; A5 misses %s, the default search %s\n",
              min(S[, "A5"]) / 60, max(S[, "A5"]) / 60, min(S[, "apollo_searchStart"]) / 60,
              max(S[, "apollo_searchStart"]) / 60, min(r), max(r),
              cells_txt(cells[M[, "A5"]]), cells_txt(cells[M[, "apollo_searchStart"]])))
}

# ---- Main -------------------------------------------------------------------

cat("dev/ladder_statistics.R: stress-ladder and isolation statistics from output/ (no RNG, no estimation)\n")
res <- lapply(setNames(names(LADDERS), names(LADDERS)), ladder_section)
screen_section(res$randomised$D)
iso_section(res$blocked$D)

tab <- do.call(rbind, ROWS)
key <- do.call(paste, c(tab[c("ladder", "set", "method", "comparator", "threshold", "statistic")], sep = "\r"))
if (anyDuplicated(key)) fail("duplicate statistic: ", gsub("\r", " / ", key[anyDuplicated(key)]))
tab$value <- signif(tab$value, 7)
write.csv(tab, OUT_STATS, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)", OUT_STATS, nrow(tab)))
if (all(vapply(res, `[[`, TRUE, "have_log"))) {
  cap <- do.call(rbind, lapply(res, `[[`, "cap"))
  cap$last_improvement <- signif(cap$last_improvement, 7)
  cap$gap_LL <- signif(cap$gap_LL, 7)
  write.csv(cap, OUT_CAP, row.names = FALSE)
  cat(sprintf(" and %s (%d rows).\n", OUT_CAP, nrow(cap)))
} else cat(sprintf("; %s left as it is (a log is missing).\n", OUT_CAP))
