#!/usr/bin/env Rscript
# =============================================================================
# dev/check_empirical_lcmnl_maxima.R
#
# Why. The second-start rule (dev/run_mmnl_second_start.R) re-estimates the
# MMNL benchmark from a second start, because a local maximum there biases the
# LCMNL-versus-MMNL reading toward the LCMNL. A local maximum on the LCMNL side
# biases it toward the MMNL, and the stored LCMNL ladders were fitted under
# several klue versions. This script is the LCMNL side of the rule: an
# independent search of each dataset's BIC-best LCMNL and its neighbours.
#
# Design. Datasets, adapters, loaders and DGPs as in dev/rerun_mmnl_benchmark.R;
# stored ladders as in dev/run_mmnl_second_start.R. C* is the BIC-best C of the
# stored ladder; the checked set is {C*-1, C*, C*+1} within the ladder's range,
# or above it up to C = 12 (see Above the stored ladder; Swissmetro's C* = 12
# is the cap, so {11, 12}). At each C:
#   (i)  klue_lcmnl(db, C, dgp = dgp) at klue 0.10.0 defaults (six clustering
#        starts, screened, the winner polished): does the stored LL reproduce?
#   (ii) R random-partition starts (get_random_partition_starts(): respondents
#        assigned to classes at random, then class-wise MNL fits), each
#        followed by direct ML: estimate_lcmnl(), BFGS with analytic gradients,
#        klue's reltol 1e-10, up to RP_MAXIT = 2000 iterations (klue's default
#        cap is 500; these starts begin further from a maximum). R = 60 at C*,
#        20 at each neighbour.
# If the search moves the BIC-best C, the same rule is applied around the new
# BIC-best (60 starts there, 20 at its neighbours), until the BIC-best is a C
# already searched as such.
#
# Above the stored ladder. The paper's class-count search (step (iii) of the
# empirical procedure) goes on until BIC turns or C = 12. So when the BIC-best
# C is the highest C searched and below C_MAX = 12, C + 1 is searched as its
# neighbour even where the stored ladder has no fit: (i) and (ii) as above,
# with 20 random starts, or 60 if C + 1 becomes the BIC-best (then C + 2 is
# searched the same way). The search stops when BIC turns (the BIC-best C is
# below the highest C searched) or at C = 12. Of the five datasets only Swiss
# route gets there: its stored ladder stops at C = 6, and the search moved its
# BIC-best from 5 to 6. A C above the stored ladder has no stored fit, so its
# klue default fit is the reference instead: such rows have stored = FALSE,
# LL_stored and BIC_stored hold the klue fit, and every comparison with "the
# stored LL" is with it (dLL_klue = 0, klue_reproduces = NA). k = C * npc +
# C - 1 and BIC = -2 LL + k log(N), N respondents, as in the stored rows
# (checked against every stored row before such a cell is fitted). A stored LL
# always counts towards the best LL at its C; the klue fit at a C above the
# stored ladder counts only if admissible (see Degenerate fits).
#
# Recorded per (dataset, C): whether C is on the stored ladder (stored), the
# stored, klue-default and best random-start LL and the BIC each implies (same
# k and N), the number of random starts within TOL = 0.01 LL of the best one
# and of the stored LL, the largest gain over the stored LL, a flag when it
# exceeds TOL, and the best random start's smallest class share and largest
# |coefficient| (a gain at a vanishing class with exploding coefficients is a
# degenerate point). Per dataset: the LCMNL-versus-MMNL margin, best converged
# MMNL BIC (output/mmnl_bench_std/<dataset>_M1..M4.rds) minus BIC-best LCMNL
# BIC (> 0 favours the LCMNL, as in dev/run_mmnl_second_start.R), for the
# stored ladder (C_best_stored, lcmnl_BIC_stored, margin_stored) and for the
# ladder with the best LL found at each checked C, those above the stored
# ladder included (C_best_checked, margin_checked and the *_any columns).
#
# Degenerate fits. A fit whose largest class-specific |coefficient| (tastes and
# constants, not the share parameters) exceeds MAX_ABS_COEF = 500 is not
# admissible: the best LL found at each C, the BIC-best C and the margin use
# admissible fits only (columns *_any keep every fit). Among the fits that
# decide the table, the bound separates two groups with a wide gap. At Swiss
# route C = 6, eleven random starts end within 0.4 LL of -1406.81 with a 3%
# class and largest coefficients from 916 to 2,676: the likelihood is flat
# along that coefficient, which diverges. At Swiss route C = 7, above the
# stored ladder, klue's default fit is not admissible either (a 3% class,
# largest coefficient 2,116), nor are 10 of the 20 random starts. The largest
# coefficient of any fit the paper's tables use is 261 (klue's stored
# Swissmetro C = 12 fit). Any bound between 262 and 915 gives the same table,
# C = 7 included, except rp_n_degenerate: some starts below the best LL at
# their C have largest coefficients in that range (467, 626 and 714 at Swiss
# route C = 7).
#
# Seeds. The random starts of (dataset, C) use seed = SPECS[[dataset]]$seed + C
# (Mode 100, SwissRoute 200, Vittel 300, Electricity 400, Swissmetro 500).
# get_random_partition_starts() draws start s under seed * 1e5 + 5e4 + s, so a
# start does not depend on R, on the order of the cells or on a resume.
# klue_lcmnl() seeds its own clustering. Deterministic.
#
# Resumable. An rds checkpoint holds every cell (dataset, C): the klue-default
# fit, then the random starts (saved at least once a minute and at the end of
# the cell). A rerun continues where the last one stopped. Once every cell is
# done, a rerun refits nothing and only rebuilds the CSV, e.g. after the MMNL
# benchmark files change. The checkpoint records the klue version and RP_MAXIT;
# a run under other values stops (delete the checkpoint to redo). One instance
# at a time (lock file holding the pid).
#
# Run from the repository root, one core (LCMNL only, no Apollo):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nohup nice -n 15 Rscript dev/check_empirical_lcmnl_maxima.R \
#     >> output/empirical_lcmnl_maxima.log 2>&1 &
# Some datasets only: prefix DATASETS=Mode,SwissRoute (names as in SPECS).
# Test, about a minute: EMP_MAX_TEST=1 runs Mode at its C* with R = 2 and
# writes only to EMP_MAX_TEST_DIR (default <TMPDIR>/emp_max_test):
#   EMP_MAX_TEST=1 EMP_MAX_TEST_DIR=/path/to/dir OMP_NUM_THREADS=1 \
#     OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 15 \
#     Rscript dev/check_empirical_lcmnl_maxima.R
# Cost (measured, one core): 33 min for the five datasets on 2026-10-05, 29 of
# them Vittel, Electricity and Swissmetro; the Swiss route search above the
# stored ladder (C = 7) added 41 s on 2026-10-06. A rerun on a complete
# checkpoint takes seconds.
# Writes: output/empirical_lcmnl_maxima.csv (one row per checked dataset and C,
#           those above the stored ladder included, with stored = FALSE),
#         output/empirical_lcmnl_maxima.rds (checkpoint: per-start LLs, best fits),
#         output/.empirical_lcmnl_maxima.lock (while running)
# =============================================================================

TEST <- identical(Sys.getenv("EMP_MAX_TEST"), "1")
OUT  <- "output"
if (TEST) {
  OUT <- Sys.getenv("EMP_MAX_TEST_DIR")
  if (!nzchar(OUT)) OUT <- file.path(dirname(tempdir()), "emp_max_test")
}
CSV   <- file.path(OUT, "empirical_lcmnl_maxima.csv")
CKPT  <- file.path(OUT, "empirical_lcmnl_maxima.rds")
LOCK  <- file.path(OUT, ".empirical_lcmnl_maxima.lock")
BENCH <- "output/mmnl_bench_std"         # MMNL benchmark fits, read only
R_BEST    <- if (TEST) 2L else 60L       # random starts at the BIC-best C
R_NB      <- 20L                         # random starts at each neighbour
RP_MAXIT  <- 2000L                       # BFGS iteration cap per random start
TOL       <- 0.01                        # LL tolerance: same maximum / above it
MAX_ABS_COEF <- 500                      # admissible fits only (header, "Degenerate fits")
C_MAX     <- 12L                         # top C of the class-count search (header)
SAVE_SECS <- 60                          # checkpoint at least this often in a cell

stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     paste0(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# Adapters and loaders as in dev/rerun_mmnl_benchmark.R, ladders as in
# dev/run_mmnl_second_start.R; cheapest datasets first.
SPECS <- list(
  Mode        = list(script = "R/empirical_mode_choice.R", skip = "mode.skip",
                     loader = "load_mode_choice_database", dgp = "MODE_DGP",
                     lcmnl  = "output/mode_lcmnl_results.csv",          seed = 100L),
  SwissRoute  = list(script = "R/empirical_swiss_route.R", skip = "swiss.skip",
                     loader = "load_swiss_database",       dgp = "SWISS_DGP",
                     lcmnl  = "output/swiss_lcmnl_results.csv",         seed = 200L),
  Vittel      = list(script = "R/empirical_application.R", skip = "empirical.skip",
                     loader = "load_empirical_database",   dgp = "EMP_DGP",
                     lcmnl  = "output/vittel_respec_lcmnl_ext.csv",     seed = 300L),
  Electricity = list(script = "R/empirical_electricity.R", skip = "elec.skip",
                     loader = "load_electricity_database", dgp = "ELEC_DGP",
                     lcmnl  = "output/emp_cmax_extend_Electricity.csv", seed = 400L),
  Swissmetro  = list(script = "R/empirical_swissmetro.R",  skip = "sm.skip",
                     loader = "load_swissmetro_database",  dgp = "SM_DGP",
                     lcmnl  = "output/emp_cmax_extend_Swissmetro.csv",  seed = 500L)
)

DATASETS <- if (TEST) "Mode" else {
  d <- trimws(strsplit(Sys.getenv("DATASETS"), ",")[[1]])
  if (any(nzchar(d))) d[nzchar(d)] else names(SPECS)
}
bad <- setdiff(DATASETS, names(SPECS))
if (length(bad)) stop("unknown dataset(s): ", paste(bad, collapse = ", "),
                      "; choose from ", paste(names(SPECS), collapse = ", "))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (file.exists(LOCK)) {
  pid <- suppressWarnings(as.integer(readLines(LOCK, n = 1, warn = FALSE)))
  if (!is.na(pid) && pid != Sys.getpid() && isTRUE(tools::pskill(pid, 0L))) {
    stamp("another check_empirical_lcmnl_maxima.R is running (pid ", pid, "); exiting")
    quit(save = "no", status = 0)
  }
}
writeLines(as.character(Sys.getpid()), LOCK)

# ---- Helpers ----------------------------------------------------------------

read_ladder <- function(ds) {
  lad <- read.csv(SPECS[[ds]]$lcmnl, stringsAsFactors = FALSE)
  lad <- lad[lad$converged %in% TRUE, c("C", "LL", "k", "BIC")]
  lad$stored <- TRUE
  lad[order(lad$C), ]
}

# The reference row of a C above the stored ladder: its klue default fit, with
# k and BIC as in the stored rows (N respondents).
ext_row <- function(C, LL, k, N)
  data.frame(C = C, LL = LL, k = k, BIC = -2 * LL + k * log(N), stored = FALSE)

# The stored ladder plus the reference row of each checkpointed C above it.
ref_ladder <- function(ds, lad = read_ladder(ds)) {
  for (cell in Filter(function(x) identical(x$dataset, ds) && x$C > max(lad$C), CK$cells))
    lad <- rbind(lad, ext_row(cell$C, cell$klue$LL, cell$k, cell$N))
  lad[order(lad$C), ]
}

best_mmnl <- function(ds) {             # min BIC over the converged M1..M4
  bic <- vapply(paste0("M", 1:4), function(m) {
    f <- file.path(BENCH, sprintf("%s_%s.rds", ds, m))
    if (!file.exists(f)) return(NA_real_)
    o <- tryCatch(readRDS(f), error = function(e) NULL)   # the MMNL lane may be writing it
    if (is.null(o)) { stamp("WARNING could not read ", f, "; margin without it"); return(NA_real_) }
    if (isTRUE(o$converged)) o$BIC else NA_real_
  }, numeric(1))
  if (all(is.na(bic))) list(model = NA_character_, BIC = NA_real_)
  else list(model = names(which.min(bic)), BIC = min(bic, na.rm = TRUE))
}

save_ck <- function() {                 # write-then-rename: never a half-written file
  tmp <- paste0(CKPT, ".tmp")
  saveRDS(CK, tmp)
  if (!file.rename(tmp, CKPT)) stop("could not write ", CKPT)
}
put <- function(key, cell) { CK$cells[[key]] <<- cell; save_ck() }

CUR <- NULL                             # the dataset in memory (one at a time)
data_for <- function(ds) {
  if (!identical(CUR$ds, ds)) {
    s <- SPECS[[ds]]
    options(structure(list(TRUE), names = s$skip))   # load the adapter, skip its autorun
    suppressWarnings(suppressMessages(source(s$script, local = FALSE)))
    db  <- get(s$loader)()
    dgp <- get(s$dgp)
    CUR <<- list(ds = ds, db = db, dgp = dgp, N = length(unique(db$ID)))
    stamp(sprintf("%s: N=%d, %d rows, %d generic + price, %d constants",
                  ds, CUR$N, nrow(db), dgp$n_generic, dgp$n_asc))
  }
  CUR
}

# One output row for a cell; `row` is its reference row (ref_ladder()): the
# stored ladder row at the cell's C, or above the ladder its klue default fit.
cell_row <- function(cell, row) {
  fin  <- function(x) if (length(x) == 1L && is.finite(x)) x else NA_real_
  fmax <- function(x) if (any(is.finite(x))) max(x[is.finite(x)]) else NA_real_
  st   <- isTRUE(row$stored)            # FALSE above the stored ladder
  bic  <- function(ll)                  # same k and N as the reference row
    if (st) row$BIC - 2 * (ll - row$LL) else -2 * ll + row$k * log(cell$N)
  rp   <- cell[["rp"]]
  ll   <- if (is.null(rp)) numeric(0) else rp$LL[is.finite(rp$LL)]
  LLk  <- fin(cell$klue$LL)
  LLr  <- if (length(ll)) max(ll) else NA_real_
  best <- cell$rp_best
  gain <- c(LLk, LLr) - row$LL
  gain <- if (all(is.na(gain))) NA_real_ else max(gain, na.rm = TRUE)
  LLb_any <- fmax(c(row$LL, LLk, LLr))
  # admissible fits only: klue's (tastes and constants are the first C * npc
  # entries of $par) and the random starts with a bounded largest coefficient;
  # a stored LL always counts, a klue reference fit above the ladder only if admissible
  npc  <- (row$k - cell$C + 1) / cell$C
  kbig <- if (is.null(cell$klue$par)) NA_real_ else max(abs(cell$klue$par[seq_len(cell$C * npc)]))
  adm  <- if (is.null(rp)) rp else rp[is.finite(rp$LL) & rp$max_abs_par <= MAX_ABS_COEF, ]
  ia   <- if (NROW(adm)) which.max(adm$LL) else integer(0)
  LLa  <- if (length(ia)) adm$LL[ia] else NA_real_
  LLb  <- fmax(c(if (st) row$LL, if (isTRUE(kbig <= MAX_ABS_COEF)) LLk, LLa))
  data.frame(
    dataset = cell$dataset, C = cell$C, k = row$k, N = cell$N %||% NA_integer_,
    stored = st, LL_stored = row$LL, BIC_stored = row$BIC,
    LL_klue = LLk, BIC_klue = bic(LLk), dLL_klue = LLk - row$LL,
    klue_reproduces = if (st) abs(LLk - row$LL) <= TOL else NA,
    klue_converged = cell$klue$converged %||% NA,
    klue_best_start = cell$klue$best_start %||% NA_character_,
    klue_singular = cell$klue$singular %||% NA, klue_secs = cell$klue$secs %||% NA_real_,
    R = NROW(rp), R_target = cell$R_target %||% NA_integer_, rp_seed = cell$seed,
    rp_n_converged = if (is.null(rp)) 0L else sum(rp$converged),
    rp_n_failed = if (is.null(rp)) 0L else sum(!is.finite(rp$LL)),
    LL_rp_best = LLr, BIC_rp_best = bic(LLr), dLL_rp_best = LLr - row$LL,
    rp_n_at_best = sum(ll >= LLr - TOL),
    rp_n_reach_stored = sum(ll >= row$LL - TOL),
    rp_n_above_stored = sum(ll > row$LL + TOL),
    rp_best_start = if (is.null(best)) NA_integer_ else best$start,
    rp_best_min_share = if (is.null(best)) NA_real_ else min(best$class_probs),
    rp_best_max_abs_par = if (is.null(best)) NA_real_ else best$max_abs_par,
    rp_best_admissible = if (is.null(best)) NA else best$max_abs_par <= MAX_ABS_COEF,
    rp_n_degenerate = if (is.null(rp)) 0L else sum(is.finite(rp$LL) & rp$max_abs_par > MAX_ABS_COEF),
    LL_rp_adm = LLa, BIC_rp_adm = bic(LLa), dLL_rp_adm = LLa - row$LL,
    rp_adm_start = if (length(ia)) adm$start[ia] else NA_integer_,
    rp_adm_min_share = if (length(ia)) adm$min_share[ia] else NA_real_,
    rp_adm_max_abs_par = if (length(ia)) adm$max_abs_par[ia] else NA_real_,
    klue_max_abs_par = kbig,
    rp_secs = if (is.null(rp)) 0 else sum(rp$secs) + (cell$rp_gen_secs %||% 0),
    max_gain_LL = gain, above_stored = !is.na(gain) && gain > TOL,
    LL_best = LLb, BIC_best = bic(LLb), LL_best_any = LLb_any, BIC_best_any = bic(LLb_any),
    klue_version = KLUE_VERSION,
    stringsAsFactors = FALSE)
}

# The reference ladder (ref_ladder()) with the best admissible LL found so far
# at each checked C (any = TRUE: the best LL of any fit, degenerate ones included).
checked_ladder <- function(ds, lad, any = FALSE) {
  lad$LL_best <- lad$LL; lad$BIC_best <- lad$BIC
  sfx <- if (any) "_any" else ""
  for (cell in Filter(function(x) identical(x$dataset, ds), CK$cells)) {
    i <- match(cell$C, lad$C)
    r <- cell_row(cell, lad[i, ])
    lad$LL_best[i] <- r[[paste0("LL_best", sfx)]]; lad$BIC_best[i] <- r[[paste0("BIC_best", sfx)]]
  }
  lad
}

# ---- One cell: klue default, then random-partition starts + direct ML ------

run_cell <- function(ds, C, R, lad) {
  if (C < 2L) R <- 0L                   # C = 1 is the MNL: no partition to draw
  key  <- sprintf("%s_C%d", ds, C)
  cell <- CK$cells[[key]] %||% list(dataset = ds, C = C, seed = SPECS[[ds]]$seed + C,
                                   rp_gen_secs = 0)
  cell$R_target <- max(cell$R_target %||% 0L, R)
  n_done <- NROW(cell[["rp"]])
  if (!is.null(cell$klue) && n_done >= R) {
    stamp(sprintf("SKIP  %s (checkpointed: klue default, %d random starts)", key, n_done))
    return(invisible(NULL))
  }
  d   <- data_for(ds)
  row <- lad[lad$C == C, ]              # none above the stored ladder
  k   <- C * d$dgp$npc + C - 1L
  if (NROW(row) && (row$k != k || abs(-2 * row$LL + k * log(d$N) - row$BIC) > 1e-3))
    stop(sprintf("%s C=%d: the stored ladder row (k = %d, BIC = %.3f) does not match these data (k = %d, N = %d)",
                 ds, C, row$k, row$BIC, k, d$N))
  kl  <- lad$C * d$dgp$npc + lad$C - 1L  # above the ladder, score C as the stored rows are
  if (!NROW(row) && (any(lad$k != kl) || any(abs(-2 * lad$LL + kl * log(d$N) - lad$BIC) > 1e-3)))
    stop(sprintf("%s: k = C * npc + C - 1 and BIC = -2 LL + k log(N) (N = %d) do not reproduce the stored ladder, so C = %d above it cannot be scored like it",
                 ds, d$N, C))
  cell$N <- d$N; cell$k <- k

  if (is.null(cell$klue)) {             # (i) klue 0.10.0 defaults
    t0 <- proc.time()[["elapsed"]]
    m  <- suppressMessages(klue_lcmnl(d$db, C, dgp = d$dgp, n_cores = 1L))
    ps <- vapply(m$method_results, function(r) r$LL, numeric(1))
    cell$klue <- list(LL = m$LL, converged = isTRUE(m$converged),
                      best_start = m$best_method %||% NA_character_,
                      singular = !is.null(m$vcov) && anyNA(m$vcov), per_start = ps,
                      par = m$par, class_probs = m$class_probs,
                      secs = proc.time()[["elapsed"]] - t0)
    put(key, cell)
    ref <- if (NROW(row)) sprintf("stored %.6f, %+.6f", row$LL, m$LL - row$LL) else
      sprintf("no stored fit, so the reference; smallest class share %.4f, largest |coefficient| %.1f",
              min(m$class_probs), max(abs(m$par[seq_len(C * d$dgp$npc)])))
    stamp(sprintf("%s klue default LL %.6f (%s) from %s, converged %s, singular %s (%.1f s)",
                  key, m$LL, ref, cell$klue$best_start,
                  cell$klue$converged, cell$klue$singular, cell$klue$secs))
    stamp("      screened LL per start: ",
          paste(sprintf("%s %.2f", names(ps), ps), collapse = ", "))
  }
  if (!NROW(row)) row <- ext_row(C, cell$klue$LL, k, d$N)   # above the stored ladder
  lab <- if (row$stored) "stored" else "klue default"

  if (n_done < R) {                     # (ii) random partitions + direct ML
    t0 <- proc.time()[["elapsed"]]
    starts <- get_random_partition_starts(d$db, C, n_starts = R, seed = cell$seed, dgp = d$dgp)
    cell$rp_gen_secs <- cell$rp_gen_secs + proc.time()[["elapsed"]] - t0
    saved <- proc.time()[["elapsed"]]
    for (i in seq.int(n_done + 1L, R)) {
      t1 <- proc.time()[["elapsed"]]
      f  <- tryCatch(estimate_lcmnl(d$db, C, start_betas = starts[[i]]$betas,
                                    start_shares = starts[[i]]$shares, dgp = d$dgp,
                                    maxit = RP_MAXIT, vcov = FALSE),
                     error = function(e) NULL)
      ok  <- !is.null(f) && is.finite(f$LL)
      big <- if (ok) max(abs(f$par[seq_len(C * d$dgp$npc)])) else NA_real_  # tastes + constants
      cell$rp <- rbind(cell[["rp"]], data.frame(
        start = i, LL = if (ok) f$LL else NA_real_, converged = ok && isTRUE(f$converged),
        min_share = if (ok) min(f$class_probs) else NA_real_, max_abs_par = big,
        secs = proc.time()[["elapsed"]] - t1))
      if (ok && f$LL > (cell$rp_best$LL %||% -Inf))
        cell$rp_best <- list(LL = f$LL, start = i, converged = isTRUE(f$converged),
                             par = f$par, class_probs = f$class_probs, max_abs_par = big)
      if (i == R || proc.time()[["elapsed"]] - saved >= SAVE_SECS) {
        put(key, cell); saved <- proc.time()[["elapsed"]]
      }
      if (i %% 10L == 0L && i < R)
        stamp(sprintf("      %s start %d/%d, best LL so far %.6f (%s %.6f)",
                      key, i, R, cell$rp_best$LL %||% NA_real_, lab, row$LL))
    }
  }

  r <- cell_row(cell, row)
  stamp(sprintf("%s %d random starts: best LL %.6f (%+.6f vs %s), %d within %.2f of it, %d reach the %s LL, %d above it, %d converged (%.0f s)",
                key, r$R, r$LL_rp_best, r$dLL_rp_best, lab, r$rp_n_at_best, TOL,
                r$rp_n_reach_stored, lab, r$rp_n_above_stored, r$rp_n_converged, r$rp_secs))
  if (isTRUE(r$above_stored))
    stamp(sprintf("FLAG  %s: best LL %.6f is %.4f above the %s %.6f (klue default %.6f, best random start %.6f with smallest class share %.4f and largest |coefficient| %.1f)",
                  key, r$LL_best_any, r$max_gain_LL, lab, r$LL_stored, r$LL_klue, r$LL_rp_best,
                  r$rp_best_min_share, r$rp_best_max_abs_par))
  invisible(NULL)
}

# ---- One dataset: C* and its neighbours, again if the BIC-best moves -------
# A neighbour can be any C of the stored ladder or, above it, any C up to C_MAX;
# it is above only when the BIC-best is the highest C searched.

check_dataset <- function(ds) {
  lad  <- read_ladder(ds)
  top  <- max(lad$C)
  nbr  <- c(lad$C[lad$C >= 2L], setdiff(seq_len(C_MAX), seq_len(top)))
  tg   <- integer(max(top, C_MAX))      # random starts wanted at each C (0 = not checked)
  seen <- integer(0)
  cs   <- lad$C[which.min(lad$BIC)]
  stamp(sprintf("%s: stored ladder C = %d..%d, BIC-best C* = %d (BIC %.2f)",
                ds, min(lad$C), max(lad$C), cs, min(lad$BIC)))
  repeat {
    seen <- c(seen, cs)
    tg[cs] <- max(tg[cs], R_BEST)
    if (!TEST) for (C in intersect(c(cs - 1L, cs + 1L), nbr))
      tg[C] <- max(tg[C], R_NB)
    if (!TEST && cs + 1L > top && cs + 1L <= C_MAX)
      stamp(sprintf("%s: the BIC-best C = %d is the highest C searched, so C = %d above the stored ladder is searched too (until BIC turns or C = %d)",
                    ds, cs, cs + 1L, C_MAX))
    for (C in c(cs, setdiff(which(tg > 0L), cs))) run_cell(ds, C, tg[C], lad)
    chk <- checked_ladder(ds, ref_ladder(ds, lad))
    cs  <- chk$C[which.min(chk$BIC_best)]
    if (TEST || cs %in% seen) break
    stamp(sprintf("%s: the search moved the BIC-best LCMNL to C = %d; searching it and its neighbours too",
                  ds, cs))
  }
}

# ---- Table: every checkpointed cell, margins from the current MMNL files ---

build_table <- function() {
  dss <- intersect(names(SPECS), vapply(CK$cells, function(x) x$dataset, ""))
  if (!length(dss)) return(NULL)
  rows <- list()
  for (ds in dss) {                     # *_stored: stored ladder; *_checked*: C above it included
    lad <- read_ladder(ds); ref <- ref_ladder(ds, lad)
    chk <- checked_ladder(ds, ref); mm <- best_mmnl(ds)
    cha <- checked_ladder(ds, ref, any = TRUE)
    per_ds <- data.frame(
      C_best_stored = lad$C[which.min(lad$BIC)], C_best_checked = chk$C[which.min(chk$BIC_best)],
      lcmnl_BIC_stored = min(lad$BIC), lcmnl_BIC_checked = min(chk$BIC_best, na.rm = TRUE),
      best_mmnl = mm$model, mmnl_BIC = mm$BIC,
      margin_stored = mm$BIC - min(lad$BIC), margin_checked = mm$BIC - min(chk$BIC_best, na.rm = TRUE),
      C_best_checked_any = cha$C[which.min(cha$BIC_best)],
      margin_checked_any = mm$BIC - min(cha$BIC_best, na.rm = TRUE),
      stringsAsFactors = FALSE)
    for (cell in Filter(function(x) identical(x$dataset, ds), CK$cells))
      rows[[length(rows) + 1L]] <- cbind(cell_row(cell, ref[ref$C == cell$C, ]), per_ds)
  }
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$dataset, names(SPECS)), tab$C), ]
  rownames(tab) <- NULL
  tab
}

write_table <- function(tab) {
  if (is.null(tab)) return(invisible(NULL))
  for (cc in grep("^(LL_|dLL_|max_gain)", names(tab), value = TRUE)) tab[[cc]] <- round(tab[[cc]], 6)
  for (cc in grep("BIC|^margin_", names(tab), value = TRUE)) tab[[cc]] <- round(tab[[cc]], 3)
  for (cc in c("rp_best_min_share", "rp_best_max_abs_par", "rp_adm_min_share",
              "rp_adm_max_abs_par", "klue_max_abs_par")) tab[[cc]] <- round(tab[[cc]], 4)
  for (cc in c("klue_secs", "rp_secs")) tab[[cc]] <- round(tab[[cc]], 1)
  write.csv(tab, CSV, row.names = FALSE)
  invisible(tab)
}

# ---- Main ------------------------------------------------------------------

tryCatch({
  options(klue.cores = 1L, klue.mmnl.n_cores = 1L, mc.cores = 1L)
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
  stopifnot(packageVersion("klue") >= "0.10.0")
  KLUE_VERSION <- as.character(packageVersion("klue"))

  SETUP <- list(klue = KLUE_VERSION, rp_maxit = RP_MAXIT)
  CK <- if (file.exists(CKPT)) readRDS(CKPT) else list(setup = SETUP, cells = list())
  if (!identical(CK$setup, SETUP))
    stop(sprintf("%s was written under klue %s with RP_MAXIT = %s; this run has klue %s and %d. Delete it to start again.",
                 CKPT, CK$setup$klue, CK$setup$rp_maxit, KLUE_VERSION, RP_MAXIT))

  stamp(sprintf("klue %s; datasets %s; random starts: %d at the BIC-best C%s, up to %d BFGS iterations each",
                KLUE_VERSION, paste(DATASETS, collapse = ", "), R_BEST,
                if (TEST) " only (TEST)" else
                  sprintf(", %d at each neighbour (above the stored ladder up to C = %d)", R_NB, C_MAX),
                RP_MAXIT))
  if (TEST) stamp("TEST: writing only to ", OUT)
  failed <- character(0)                # a failing dataset does not stop the others
  for (ds in DATASETS) {
    ok <- tryCatch({ check_dataset(ds); TRUE }, error = function(e) {
      stamp("ERROR ", ds, ": ", conditionMessage(e), "; skipping this dataset"); FALSE })
    if (!ok) failed <- c(failed, ds)
    write_table(build_table())          # rebuilt after every dataset
  }

  tab <- write_table(build_table())
  if (!is.null(tab)) {
    cat("\nPer dataset and C (LL; random starts within", TOL, "of their best / reaching the stored LL):\n")
    print(tab[, c("dataset", "C", "stored", "LL_stored", "LL_klue", "LL_rp_best", "R", "rp_n_at_best",
                  "rp_n_reach_stored", "max_gain_LL", "above_stored")], row.names = FALSE)
    cat("\nAdmissible fits (largest |coefficient| <=", MAX_ABS_COEF, "):\n")
    print(tab[, c("dataset", "C", "LL_rp_adm", "dLL_rp_adm", "rp_adm_min_share",
                  "rp_adm_max_abs_par", "rp_n_degenerate", "LL_best", "BIC_best")], row.names = FALSE)
    cat("\nPer dataset (margin = best MMNL BIC - BIC-best LCMNL BIC; > 0 favours the LCMNL):\n")
    print(unique(tab[, c("dataset", "C_best_stored", "C_best_checked", "best_mmnl", "mmnl_BIC",
                         "margin_stored", "margin_checked", "C_best_checked_any",
                         "margin_checked_any")]), row.names = FALSE)
  }
  if (length(failed)) stop("failed: ", paste(failed, collapse = ", "), " (see the ERROR lines)")
  stamp("DONE -> ", CSV)
}, finally = unlink(LOCK))
