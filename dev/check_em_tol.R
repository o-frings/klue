# dev/check_em_tol.R
#
# Why. On the randomised package ladder (dev/stress_replicate.R) clustering+EM
# (arm klue_em) ends below clustering+ML (arm klue_ml) in some cells, although
# both fit the same model from the same six clustering starts. EM stops when
# the LL changes by less than 1e-6 between iterations, or after 500 iterations
# (estimate_lcmnl_em's defaults, which the ladder uses). Is the gap a different
# optimum, or EM stopping early? This script refits EM from the six starts at
# that setting (loose) and at a tolerance of 1e-12 with up to 5000 iterations
# (tight), and compares the best of the six with klue_ml's LL. If the tight
# fit closes the gap, the gap was early stopping; if not, EM ends at a
# different, worse stationary point.
#
# Backs. version12.tex, Section 3.1, paragraph "Either estimator: ML or EM":
# in three of the five seeds of the randomised ladder's correlated rung EM
# ends at a worse stationary point than direct ML from the same clustering
# starts, by 2.3 to 5.7 LL, and running EM to a stopping tolerance of 1e-12
# changes none of these gaps by more than 1e-4 LL.
#
# Cells. Every seed (1-5) of the correlated rung hard_corr (K = 5,
# kap = 0.50, sig = 0.30, attr_corr = 0.6), and very_hard seed 4 (K = 5,
# kap = 0.30, sig = 0.35). Until 2026-10-07 the script checked only very_hard
# seed 4 and hard_corr seeds 1 and 4, and wrote no file.
#
# Data and seeds; the script is deterministic. RUNGS, N_PER_CLASS (60) and DGP
# are read from dev/stress_replicate.R (parsed, not run), and the data are
# generated as its loop does: klue_simulate(N_per_class = 60, T_tasks = 12,
# ...), attr_corr = 0.6 on hard_corr, data seed
# as.integer(1000*K + 100*(kap*100) + seed) (hard_corr 10000 + seed,
# very_hard seed 4 8004), which klue_simulate passes to set.seed. Of klue's six
# clustering starts, kmeans, Mclust and pam run under set.seed(123) (klue's
# .with_seed, which restores the caller's stream); the three hclust starts,
# BFGS and EM use no random numbers. options(klue.screen = FALSE), set first,
# is klue's pre-0.9.1 path, with which the ladder fits its klue_ml arm
# (SCREEN_OFF in the runner): tight cluster-wise MNL fits for the starts, and
# every ML start fitted at full precision. With it ml_best reproduces the
# cell's stored klue_ml LL in output/stress_replicate.rds, the check that the
# data are the ladder's (data_ok: |ml_best - ml_stored| <= 1e-3, the rule of
# dev/check_ladder_carryover.R). EM here starts from the same tight start
# fits, while the runner fits klue_em at klue's default, with looser start
# fits. Both give the stored klue_em LLs: em_loose here, and the default in
# dev/check_ladder_carryover.R.
#
# Run from the repository root (one core, about 25 min; 3 to 5.5 min a cell):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/check_em_tol.R
#
# Reads dev/stress_replicate.R (parsed) and output/stress_replicate.rds.
# Writes output/em_tol_check.csv, one row per cell: rung, seed, cell_seed;
# ml_best (klue_ml refitted), ml_stored, data_ok; em_stored (the ladder's
# klue_em LL); em_loose, em_tight (best of the six EM fits at each setting);
# gap_loose, gap_tight (EM minus ml_best); change_tight_minus_loose;
# max_start_change (largest |tight - loose| over the six starts); loose_at_cap
# (starts whose loose fit used all 500 iterations); max_iter_tight (most
# iterations of a tight fit; below 5000, every tight fit stopped on the
# tolerance); klue_version; date. The console shows each start's fits and the
# table.
#
# Result (run 2026-10-07, klue 0.10.0, 25 min): in all six cells ml_best and
# em_loose equal the stored klue_ml and klue_em LLs (to 1e-7). On the
# correlated rung EM ends below ML in seeds 1, 3 and 5 (gaps -5.72, -2.32,
# -4.19 LL) and ties it in seeds 2 and 4; very_hard seed 4: -2.43. The tight
# setting changes no gap by more than 2.3e-5 LL. Six loose fits stopped at the
# 500-iteration cap, none of them its cell's best start; the largest change at
# one start is 0.014 LL (hard_corr seed 3, gmm start). Every tight fit stopped
# on the tolerance (at most 3447 iterations).

options(klue.screen = FALSE)   # before any fit: the ladder's klue_ml path (see above)
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))

RUNNER   <- "dev/stress_replicate.R"
STORED   <- "output/stress_replicate.rds"
OUT      <- "output/em_tol_check.csv"
LOOSE    <- list(tol = 1e-6,  maxit = 500L)    # estimate_lcmnl_em's defaults, as in the ladder
TIGHT    <- list(tol = 1e-12, maxit = 5000L)   # tight: rule out early stopping
TOL_DATA <- 1e-3                               # data_ok: |ml_best - ml_stored| <= TOL_DATA
CELLS    <- data.frame(rung = c("very_hard", rep("hard_corr", 5)), seed = c(4L, 1:5))
if (!file.exists(RUNNER)) stop("run from the repository root", call. = FALSE)

# The runner's own RUNGS, N_PER_CLASS and DGP (parsed, not run).
rd <- new.env()
for (e in parse(RUNNER, keep.source = FALSE))
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
      as.character(e[[2]]) %in% c("RUNGS", "N_PER_CLASS", "DGP")) eval(e, rd)
RUNGS <- setNames(rd$RUNGS, vapply(rd$RUNGS, function(r) r$label, ""))
N_PER_CLASS <- rd$N_PER_CLASS; DGP <- rd$DGP
stored <- readRDS(STORED)

# EM from each start at one setting: the LL and the number of iterations.
em_fits <- function(db, K, starts, set)
  vapply(starts, function(st) {
    r <- if (!is.null(st)) tryCatch(
      estimate_lcmnl_em(db, K, start_betas = st$betas, start_shares = st$shares,
                        dgp = DGP, max_em_iter = set$maxit, tol = set$tol),
      error = function(e) NULL)
    if (is.null(r) || !is.finite(r$LL)) c(LL = NA_real_, iters = NA_real_)
    else c(LL = r$LL, iters = r$em_iters)
  }, c(LL = 0, iters = 0))

t_all <- proc.time()[["elapsed"]]
rows <- list()
for (i in seq_len(nrow(CELLS))) {
  rung <- CELLS$rung[i]; s <- CELLS$seed[i]; r <- RUNGS[[rung]]
  key <- paste0(rung, "|seed", s)
  t0 <- proc.time()[["elapsed"]]
  # the cell's data, as the runner's loop generates them
  seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + s)
  d <- klue_simulate(N_per_class = N_PER_CLASS, T_tasks = 12, true_K = r$K,
                     separation = r$kap, heterogeneity = r$sig,
                     seed = seed_i, dgp = DGP, attr_corr = r$corr)
  db <- d$database
  ml <- klue_lcmnl(db, r$K, dgp = DGP, estimator = "ml")$LL   # the ladder's klue_ml arm
  starts <- get_all_starts(db, r$K, dgp = DGP)                  # its six clustering starts
  lo <- em_fits(db, r$K, starts, LOOSE)
  ti <- em_fits(db, r$K, starts, TIGHT)
  cat(sprintf("\n%s (cell seed %d): %.1f min\n", key, seed_i,
              (proc.time()[["elapsed"]] - t0) / 60))
  for (nm in colnames(lo))
    cat(sprintf("  %-12s loose %.6f (%4.0f it)   tight %.6f (%4.0f it)   change %.1e\n",
                nm, lo["LL", nm], lo["iters", nm], ti["LL", nm], ti["iters", nm],
                ti["LL", nm] - lo["LL", nm]))
  em_loose <- max(lo["LL", ], na.rm = TRUE); em_tight <- max(ti["LL", ], na.rm = TRUE)
  ml_stored <- unname(stored[[key]]$LL[["klue_ml"]])
  rows[[key]] <- data.frame(
    rung = rung, seed = s, cell_seed = seed_i,
    ml_best = ml, ml_stored = ml_stored, data_ok = abs(ml - ml_stored) <= TOL_DATA,
    em_stored = unname(stored[[key]]$LL[["klue_em"]]),
    em_loose = em_loose, em_tight = em_tight,
    gap_loose = em_loose - ml, gap_tight = em_tight - ml,
    change_tight_minus_loose = em_tight - em_loose,
    max_start_change = max(abs(ti["LL", ] - lo["LL", ]), na.rm = TRUE),
    loose_at_cap = sum(lo["iters", ] >= LOOSE$maxit, na.rm = TRUE),
    max_iter_tight = max(ti["iters", ], na.rm = TRUE),
    klue_version = as.character(utils::packageVersion("klue")),
    date = format(Sys.Date()))
}

tab <- do.call(rbind, rows)
out <- tab
ll  <- c("ml_best", "ml_stored", "em_stored", "em_loose", "em_tight", "gap_loose", "gap_tight")
chg <- c("change_tight_minus_loose", "max_start_change")
out[ll] <- lapply(out[ll], round, 7); out[chg] <- lapply(out[chg], signif, 3)
write.csv(out, OUT, row.names = FALSE)

cat(sprintf("\n%-10s %4s | %10s %10s | %10s %10s | %8s %8s | %8s %9s\n", "rung", "seed",
            "ml_best", "ml_stored", "em_loose", "em_tight", "loose-ml", "tight-ml",
            "change", "max_start"))
for (i in seq_len(nrow(tab))) with(tab[i, ], cat(sprintf(
  "%-10s %4d | %10.2f %10.2f | %10.2f %10.2f | %+8.2f %+8.2f | %8.1e %9.1e\n", rung, seed,
  ml_best, ml_stored, em_loose, em_tight, gap_loose, gap_tight,
  change_tight_minus_loose, max_start_change)))
if (!all(tab$data_ok))
  cat("NOTE data check failed (ml_best is not the stored klue_ml LL):",
      paste(rownames(tab)[!tab$data_ok], collapse = " "), "\n")
cat("\nchange = em_tight - em_loose; max_start = largest |tight - loose| over the six starts.\n")
cat("Reading: if tight-ml ~ 0 but loose-ml < 0 => premature stopping (artifact).\n")
cat("         if tight-ml still < 0          => EM ends at a worse optimum from these starts.\n")
cat(sprintf("Wrote %s; %.1f min in all.\n", OUT, (proc.time()[["elapsed"]] - t_all) / 60))
