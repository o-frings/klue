#!/usr/bin/env Rscript
# =============================================================================
# dev/run_ladder_reference_optimum.R   (Stage 1)
#
# Why. The blocked stress ladder (dev/stress_replicate_blocked.R) fits each
# method at the true class count K* only and scores it against the best LL any
# method reached in that cell. That best is not an independent optimum, and a
# shortfall at K* does not by itself say whether the method would select a
# different number of classes. This script computes an independent reference
# optimum per cell at C = 1..6 and, from its BIC profile, whether each ladder
# method's shortfall at K* could change the BIC-selected class count.
#
# Cells (21): hard|seed1-10, very_hard|seed1-10, moderate|seed2. The data are
# regenerated as dev/stress_replicate_blocked.R generates them: the shared
# blocked design klue_design(n_cards = 48, n_blocks = 4) (its fixed internal
# seed), klue_dgp(n_generic = 4, n_alternatives = 3), 60 respondents per class,
# the rung's K*, separation and heterogeneity, and the data seed
# as.integer(1000*K + 100*(kap*100) + seed), written exactly as there.
#
# Per cell:
#   1. Check, for all cells before any other fit: klue_ml at C = K*
#      (klue_lcmnl(db, K*, estimator = "ml"), the call inside run_klue in
#      dev/compare_packages.R) must reproduce the ladder's klue_ml LL within
#      1e-3; the difference is stored. If any cell fails, the script stops
#      before fitting; REF_ALLOW_MISMATCH=1 fits those cells anyway, flagged.
#   2. Fits at C = 1..6 from three sources:
#        klue_ml     the six clustering starts + direct ML (klue's default
#                    screen and polish); at K* the check fit is reused
#        klue_em_rp  run_klue_em_rp (dev/compare_packages.R): six random
#                    partitions + EM, seed = the data seed, so at K* it is the
#                    ladder's klue_em_rp arm
#        rp_ml       50 random partitions + direct ML: respondents drawn
#                    uniformly into C classes (redrawn, up to 20 times, until
#                    each class has at least 3, as in klue's
#                    get_random_partition_starts), class-wise MNL starts
#                    (fit_cluster_mnls), then estimate_lcmnl to full precision
#                    (reltol 1e-10, maxit 500, continued once from where it
#                    stopped if it hit the cap); one start at C = 1 (MNL,
#                    concave log-likelihood)
#   3. Reference LL per C = the maximum over the three sources. BIC per C =
#      -2 LL + k log N, with k and N read off klue's fits: the script stops
#      unless they are k = 8C - 1 and N = number of respondents. C_ref =
#      argmin of the reference BIC; the runner-up margin is the BIC gap from
#      C_ref to the C with the second-lowest BIC.
#   4. Decision effect, analytic (no fits), for each ladder method m (the names
#      in the ladder's rds): g = refLL(K*) - LL_m(K*). The ladder fits K* only.
#      If m reaches the reference LL at every other C, its BIC profile is the
#      reference one with BIC(K*) raised by 2g; with d = min over C != K* of
#      BIC_ref(C) - BIC_ref(K*), m then selects K* iff 2g < d (C_sel). If
#      C_ref = K* (d > 0), the selection changes iff 2g > d, the margin to the
#      runner-up; 2g > d stays necessary for a change if m falls short at
#      other C too (lower LLs there only raise their BIC), so long as m does
#      not beat the reference there. If C_ref != K*, a shortfall at K*
#      (g > 0) cannot change the selection; only an LL above the reference by
#      more than -d/2 could.
#
# Seeds. Data and design as above; klue's clustering starts use its internal
# seed 123; klue_em_rp partition s (1..6) uses seed_i*100000 + 50000 + s
# (klue's generator); rp_ml partition s (1..50) at class count C uses
# set.seed(seed_i*1000 + 100*C + s), disjoint from those. Deterministic.
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   nohup nice -n 15 Rscript dev/run_ladder_reference_optimum.R \
#     > output/ladder_reference_optimum.log 2>&1 &
# Expected 8-20 min per cell, about 4-5 h in all and up to 7 h (extrapolated
# from logged klue timings; EM at C = 6 is the main unknown; the log prints the
# seconds per source and C). Resumable: the rds is saved after each cell and a
# rerun skips finished cells (a cell stored under other settings is refitted).
# Every run rebuilds the CSV from the rds and the current ladder file, so a
# rerun after the ladder is topped up refreshes the decision table without
# fitting.
#
# Reads output/stress_replicate_blocked.rds (copied first; never written).
# Writes (files of this script only):
#   output/ladder_reference_optimum.rds  named list by cell: per source and C
#     the LL, seconds and convergence; k per C; the rp_ml LL of every start;
#     parameter vectors of the two ML sources; reference LL, BIC and the
#     sources reaching it per C; the K* check
#   output/ladder_reference_optimum.csv  one row per (cell, ladder method): the
#     method's LL at K* (LL_K), g, the class count it would select (C_sel),
#     whether that differs from C_ref (could_change), bic_margin_K (d above),
#     and per-cell columns: C_ref, C_runner_up, bic_margin (runner-up minus
#     C_ref), the reference LL at K* and the sources within 1e-3 of it, the
#     rp_ml starts within 1e-3 of it, ladder best minus reference, the K*
#     check differences for klue_ml and klue_em_rp, each source's own BIC
#     choice over the fitted C, and the reference BIC at C = 1..6
#
# REF_TEST=1: hard|seed5 only, C in {K*-1, K*}, 3 rp_ml starts; reads the
# ladder file and writes only to REF_TEST_DIR (default: a folder in R's
# tempdir()), starting fresh each time; a few minutes.
# =============================================================================

TEST <- identical(Sys.getenv("REF_TEST"), "1")
ALLOW_MISMATCH <- identical(Sys.getenv("REF_ALLOW_MISMATCH"), "1")
if (!file.exists("dev/compare_packages.R")) stop("run from the repository root")

LADDER  <- "output/stress_replicate_blocked.rds"
OUT_DIR <- "output"
if (TEST) {
  OUT_DIR <- Sys.getenv("REF_TEST_DIR")
  if (!nzchar(OUT_DIR)) OUT_DIR <- file.path(tempdir(), "ladder_reference_optimum_test")
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (TEST && normalizePath(OUT_DIR) == normalizePath("output", mustWork = FALSE))
  stop("REF_TEST_DIR points at output/; the test writes elsewhere")
OUT_RDS <- file.path(OUT_DIR, "ladder_reference_optimum.rds")
OUT_CSV <- file.path(OUT_DIR, "ladder_reference_optimum.csv")

# Rungs as in dev/stress_replicate_blocked.R (K double, as there).
RUNGS <- list(moderate  = list(K = 4, kap = 0.75, sig = 0.25),
              hard      = list(K = 5, kap = 0.50, sig = 0.30),
              very_hard = list(K = 5, kap = 0.30, sig = 0.35))
CELLS <- if (TEST) "hard|seed5" else
  c(sprintf("hard|seed%d", 1:10), sprintf("very_hard|seed%d", 1:10), "moderate|seed2")
N_PER_CLASS <- 60L
C_MAX     <- 6L
N_RP      <- if (TEST) 3L else 50L   # rp_ml starts per C (one at C = 1)
RP_MAXIT  <- 500L                    # klue's MAX_ITER; a start that hits it is continued once
REPRO_TOL <- 1e-3                    # tolerance of the K* check
HIT_TOL   <- 1e-3                    # an LL within this of the reference reaches it
SOURCES   <- c("klue_ml", "klue_em_rp", "rp_ml")
cs_of       <- function(K) if (TEST) c(K - 1, K) else seq_len(C_MAX)
settings_of <- function(K) list(C = cs_of(K), n_rp = N_RP, rp_maxit = RP_MAXIT)

stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     sprintf(...))); flush.console() }
`%||%` <- function(a, b) if (is.null(a)) b else a
pick <- function(v, m) if (m %in% names(v)) unname(v[[m]]) else NA_real_
ladder_best <- function(lc) { v <- lc$LL[is.finite(lc$LL)]; if (length(v)) max(v) else NA_real_ }
key_parts <- function(key) list(rung = sub("\\|seed[0-9]+$", "", key),
                                seed = as.integer(sub("^.*\\|seed", "", key)))
read_ladder <- function() {   # copy first: the ladder script rewrites its file after each cell
  for (i in 1:5) {
    tmp <- tempfile(tmpdir = if (TEST) OUT_DIR else tempdir(), fileext = ".rds")
    ok <- file.copy(LADDER, tmp, overwrite = TRUE)
    x <- if (ok) tryCatch(readRDS(tmp), error = function(e) NULL)
    unlink(tmp)
    if (!is.null(x)) return(x)
    Sys.sleep(5)
  }
  stop("could not read ", LADDER)
}
save_rds <- function(x, path) {   # write then rename, so an interrupted save keeps the old file
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp)
  invisible(file.rename(tmp, path))
}

sys.source("dev/compare_packages.R", envir = globalenv())   # functions only (guarded)
options(mc.cores = 1L)
DGP <- klue_dgp(n_generic = 4, n_alternatives = 3)

# The ladder's data for one cell, from its seed formula and the shared design.
cell_data <- function(key) {
  kp <- key_parts(key); r <- RUNGS[[kp$rung]]
  seed_i <- as.integer(1000 * r$K + 100 * (r$kap * 100) + kp$seed)
  d <- klue_simulate(N_per_class = N_PER_CLASS, true_K = r$K, separation = r$kap,
                     heterogeneity = r$sig, seed = seed_i, dgp = DGP, design = DESIGN)
  c(kp, r, list(seed_i = seed_i, db = d$database))
}

# ---- the three sources ------------------------------------------------------
fit_klue_ml <- function(db, C) {
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(klue_lcmnl(db, C, dgp = DGP, estimator = "ml"),
                error = function(e) { message("klue_ml C=", C, " error: ", conditionMessage(e)); NULL })
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(f) || !is.finite(f$LL)) return(list(LL = NA_real_, converged = FALSE, secs = secs))
  list(LL = f$LL, k = f$k, BIC = f$BIC, converged = isTRUE(f$converged), par = f$par,
       secs = secs)
}

fit_em_rp <- function(db, C, seed_i) {
  f <- run_klue_em_rp(db, C, DGP, seed = seed_i)
  list(LL = as.numeric(f$LL), converged = isTRUE(f$converged), secs = as.numeric(f$seconds))
}

# One random-partition start, drawn as klue's get_random_partition_starts draws
# one, from a seed of this script.
rp_start <- function(db, C, seed) {
  N <- length(unique(db$ID))
  set.seed(seed)
  lab <- sample.int(C, N, replace = TRUE)
  tries <- 0L
  while (any(tabulate(lab, C) < 3L) && tries < 20L) {
    lab <- sample.int(C, N, replace = TRUE); tries <- tries + 1L
  }
  fit_cluster_mnls(lab, db, dgp = DGP)
}

fit_rp_ml <- function(db, C, seed_i, n_starts) {
  t0 <- proc.time()[["elapsed"]]
  if (C == 1) n_starts <- 1L
  lls <- rep(NA_real_, n_starts); conv <- rep(NA, n_starts); best <- NULL
  for (i in seq_len(n_starts)) {
    st <- rp_start(db, C, seed_i * 1000 + 100 * C + i)
    f <- tryCatch(estimate_lcmnl(db, C, start_betas = st$betas, start_shares = st$shares,
                                 dgp = DGP, maxit = RP_MAXIT, vcov = FALSE),
                  error = function(e) NULL)
    if (!is.null(f) && is.finite(f$LL) && !isTRUE(f$converged)) {   # stopped on the cap
      f2 <- tryCatch(estimate_lcmnl(db, C, dgp = DGP, start_par = f$par, maxit = RP_MAXIT,
                                    vcov = FALSE),
                     error = function(e) NULL)
      if (!is.null(f2) && is.finite(f2$LL) && f2$LL >= f$LL) f <- f2
    }
    if (is.null(f) || !is.finite(f$LL)) next
    lls[i] <- f$LL; conv[i] <- isTRUE(f$converged)
    if (is.null(best) || f$LL > best$LL) best <- f
  }
  secs <- proc.time()[["elapsed"]] - t0
  if (is.null(best)) return(list(LL = NA_real_, converged = FALSE, secs = secs,
                                 lls = lls, conv = conv))
  list(LL = best$LL, k = best$k, BIC = best$BIC, converged = isTRUE(best$converged),
       par = best$par, lls = lls, conv = conv, secs = secs)
}

# ---- selection and the decision table ---------------------------------------
# From a BIC profile named by C: C_ref, the runner-up and its margin, and
# d = min over C != K* of BIC(C) minus BIC(K*) with C_other its minimiser.
ref_selection <- function(bic, K) {
  b <- bic[is.finite(bic)]
  if (!length(b)) return(list(C_ref = NA_integer_, C_runner_up = NA_integer_,
                              bic_margin = NA_real_, d = NA_real_, C_other = NA_integer_))
  o <- order(b); Cs <- as.integer(names(b))[o]; bo <- unname(b[o])
  kk <- as.character(K); other <- b[names(b) != kk]
  list(C_ref       = Cs[1],
       C_runner_up = if (length(Cs) > 1L) Cs[2] else NA_integer_,
       bic_margin  = if (length(Cs) > 1L) bo[2] - bo[1] else NA_real_,
       d           = if (kk %in% names(b) && length(other)) min(other) - b[[kk]] else NA_real_,
       C_other     = if (length(other)) as.integer(names(other))[which.min(other)] else NA_integer_)
}

# Rows of the CSV for one cell: one per ladder method, with the cell's columns.
cell_rows <- function(z, lc) {
  kk <- as.character(z$K)
  sel <- ref_selection(z$ref_BIC, z$K)
  refK <- unname(z$ref_LL[kk])
  ll <- lc$LL %||% numeric(0)
  if (!length(ll)) ll <- c(none = NA_real_)   # cell missing from the ladder file
  g <- refK - unname(ll)
  C_sel <- ifelse(is.finite(g) & is.finite(sel$d),
                  ifelse(2 * g < sel$d, as.integer(z$K), sel$C_other), NA_integer_)
  own <- vapply(SOURCES, function(src) {   # each source's own BIC choice over the fitted C
    b <- -2 * z$LL[, src] + z$k * log(z$N)
    if (any(is.finite(b))) as.integer(z$C[which.min(b)]) else NA_integer_
  }, integer(1))
  lb <- ladder_best(lc); rp <- z$rp_lls[[kk]]
  chk_diff <- z$LL[kk, "klue_ml"] - pick(ll, "klue_ml")
  out <- data.frame(
    cell = z$cell, rung = z$rung, seed = z$seed, K = z$K, N = z$N,
    method = names(ll), LL_K = round(unname(ll), 4), g = round(g, 4),
    C_sel = C_sel, could_change = C_sel != sel$C_ref,
    C_ref = sel$C_ref, C_runner_up = sel$C_runner_up,
    bic_margin = round(sel$bic_margin, 3), bic_margin_K = round(sel$d, 3),
    LL_ref_K = round(refK, 4), ref_source_K = unname(z$ref_source[kk]),
    rp_hits_K = sum(rp >= refK - HIT_TOL, na.rm = TRUE), rp_starts = length(rp),
    ladder_best_K = round(lb, 4), ladder_best_minus_ref = round(lb - refK, 4),
    repro_diff = signif(chk_diff, 3), repro_ok = abs(chk_diff) <= REPRO_TOL,
    em_rp_diff = signif(z$LL[kk, "klue_em_rp"] - pick(ll, "klue_em_rp"), 3),
    C_bic_klue_ml = own[["klue_ml"]], C_bic_klue_em_rp = own[["klue_em_rp"]],
    C_bic_rp_ml = own[["rp_ml"]],
    row.names = NULL, stringsAsFactors = FALSE)
  for (C in seq_len(C_MAX))
    out[[paste0("BIC_ref_C", C)]] <- round(unname(z$ref_BIC[as.character(C)]), 3)
  out
}

write_table <- function(results, lad) {
  done <- intersect(CELLS, names(results))
  if (!length(done)) return(NULL)
  tab <- do.call(rbind, lapply(done, function(key) cell_rows(results[[key]], lad[[key]])))
  write.csv(tab, OUT_CSV, row.names = FALSE)
  tab
}

# ---- cells to fit -------------------------------------------------------------
results <- if (!TEST && file.exists(OUT_RDS)) readRDS(OUT_RDS) else list()
todo <- Filter(function(key) !identical(results[[key]]$settings,
                                        settings_of(RUNGS[[key_parts(key)$rung]]$K)), CELLS)
for (key in intersect(todo, names(results)))
  stamp("NOTE %s is stored with other settings; refitting it", key)
lad <- read_ladder()
stamp("%s: %d of %d cells to fit -> %s", if (TEST) "TEST" else "run", length(todo),
      length(CELLS), OUT_RDS)

# ---- step 1: the K* check, all cells first -------------------------------------
chk <- list()
if (length(todo)) {
  stamp("building the ladder's blocked design")
  DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)
  stamp("design built; checking klue_ml at K* against the ladder")
  for (key in todo) {
    cd <- cell_data(key)
    f <- fit_klue_ml(cd$db, cd$K)
    lad_ml <- pick(lad[[key]]$LL, "klue_ml")
    chk[[key]] <- list(fit = f, ladder_LL = lad_ml, diff = f$LL - lad_ml,
                       ok = isTRUE(abs(f$LL - lad_ml) <= REPRO_TOL))
    stamp("CHECK %-16s klue_ml at K*=%d  LL %.4f  ladder %.4f  diff %+.1e  %s", key,
          as.integer(cd$K), f$LL, lad_ml, f$LL - lad_ml, if (chk[[key]]$ok) "ok" else "FAIL")
  }
  bad <- names(chk)[!vapply(chk, function(x) x$ok, logical(1))]
  if (length(bad) && !ALLOW_MISMATCH)
    stop(sprintf(paste("klue_ml at K* does not reproduce the ladder's LL in %d cell(s): %s.",
                       "Fix the data regeneration first, or set REF_ALLOW_MISMATCH=1 to fit",
                       "these cells and flag them."), length(bad), paste(bad, collapse = ", ")),
         call. = FALSE)
  if (length(bad)) stamp("NOTE REF_ALLOW_MISMATCH=1: fitting %d failing cell(s), flagged", length(bad))
}

# ---- steps 2-3: fits at each C and the reference ---------------------------------
for (key in todo) {
  cd <- cell_data(key); db <- cd$db; K <- cd$K
  N_resp <- length(unique(db$ID))
  if (N_resp != N_PER_CLASS * K) stop("unexpected number of respondents in ", key)
  CS <- cs_of(K); cc <- as.character(CS)
  LLm <- secs <- matrix(NA_real_, length(CS), length(SOURCES), dimnames = list(cc, SOURCES))
  conv <- matrix(NA, length(CS), length(SOURCES), dimnames = list(cc, SOURCES))
  k_C <- setNames(rep(NA_real_, length(CS)), cc)
  rp_lls <- rp_conv <- par_ml <- par_rp <- setNames(vector("list", length(CS)), cc)
  stamp("RUN  %s  (K*=%d, N=%d, C = %s)", key, as.integer(K), N_resp, paste(CS, collapse = ","))
  for (i in seq_along(CS)) {
    C <- CS[i]
    fits <- list(klue_ml    = if (C == K) chk[[key]]$fit else fit_klue_ml(db, C),
                 klue_em_rp = fit_em_rp(db, C, cd$seed_i),
                 rp_ml      = fit_rp_ml(db, C, cd$seed_i, N_RP))
    for (src in SOURCES) {
      LLm[i, src] <- fits[[src]]$LL; secs[i, src] <- fits[[src]]$secs
      conv[i, src] <- fits[[src]]$converged
    }
    rp_lls[i] <- list(fits$rp_ml$lls); rp_conv[i] <- list(fits$rp_ml$conv)
    par_ml[i] <- list(fits$klue_ml$par); par_rp[i] <- list(fits$rp_ml$par)
    # k and N as klue's BIC uses them, read off an ML fit at this C.
    kf <- if (is.finite(fits$klue_ml$LL)) fits$klue_ml else fits$rp_ml
    if (!is.finite(kf$LL)) stop(sprintf("%s C=%d: neither ML source returned a fit", key, as.integer(C)))
    N_bic <- exp((kf$BIC + 2 * kf$LL) / kf$k)
    if (!isTRUE(kf$k == 8 * C - 1) || !isTRUE(abs(N_bic / N_resp - 1) < 1e-8))
      stop(sprintf("%s C=%d: klue's BIC uses k = %s and N = %s, not k = 8C - 1 and N = %d respondents",
                   key, as.integer(C), format(kf$k), format(N_bic), N_resp))
    k_C[i] <- kf$k
    ref <- max(LLm[i, is.finite(LLm[i, ])], -Inf)
    stamp("  C=%d  klue_ml %.3f (%.0f s)  klue_em_rp %.3f (%.0f s)  rp_ml %.3f (%.0f s; %d of %d starts at the reference, %d capped)  ref %.3f",
          as.integer(C), LLm[i, "klue_ml"], secs[i, "klue_ml"], LLm[i, "klue_em_rp"],
          secs[i, "klue_em_rp"], LLm[i, "rp_ml"], secs[i, "rp_ml"],
          sum(fits$rp_ml$lls >= ref - HIT_TOL, na.rm = TRUE), length(fits$rp_ml$lls),
          sum(!fits$rp_ml$conv, na.rm = TRUE), ref)
  }
  ref_LL  <- apply(LLm, 1, function(v) if (any(is.finite(v))) max(v[is.finite(v)]) else NA_real_)
  ref_src <- vapply(cc, function(r) {
    v <- LLm[r, ]; ok <- is.finite(v)
    if (any(ok)) paste(SOURCES[ok & v >= max(v[ok]) - HIT_TOL], collapse = "+") else ""
  }, character(1))
  ref_BIC <- -2 * ref_LL + k_C * log(N_resp)
  results[[key]] <- list(
    cell = key, rung = cd$rung, seed = cd$seed, K = K, kap = cd$kap, sig = cd$sig,
    seed_i = cd$seed_i, N = N_resp, settings = settings_of(K), C = CS, k = k_C,
    LL = LLm, secs = secs, converged = conv, rp_lls = rp_lls, rp_conv = rp_conv,
    par = list(klue_ml = par_ml, rp_ml = par_rp),
    ref_LL = ref_LL, ref_BIC = ref_BIC, ref_source = ref_src,
    check = chk[[key]][c("ladder_LL", "diff", "ok")],
    klue_version = as.character(utils::packageVersion("klue")),
    R_version = R.version.string, date = format(Sys.time(), "%Y-%m-%d %H:%M"))
  save_rds(results, OUT_RDS)
  tryCatch(write_table(results, lad),   # the fits are saved; a table error must not end the run
           error = function(e) stamp("NOTE CSV not written: %s", conditionMessage(e)))
  sel <- ref_selection(ref_BIC, K)
  stamp("done %s: C_ref = %s (runner-up %s, BIC margin %.2f); ladder best minus reference at K* %+.3f",
        key, sel$C_ref, sel$C_runner_up, sel$bic_margin,
        ladder_best(lad[[key]]) - ref_LL[[as.character(K)]])
}

# ---- step 4: decision table and summary ---------------------------------------------
lad <- read_ladder()   # current file: the ladder may have been topped up meanwhile
tab <- write_table(results, lad)
if (is.null(tab)) { stamp("no finished cells"); quit(save = "no", status = 0) }

cells <- tab[!duplicated(tab$cell), ]
cat("\n===== reference optimum per cell (BIC on respondents, k = 8C - 1) =====\n")
for (i in seq_len(nrow(cells))) {
  x <- cells[i, ]
  cat(sprintf(paste0("%-16s K*=%d  C_ref=%s  runner-up=%s  margin %7.2f  ladder best - ref %+8.3f",
                     "  check %+.1e  rp_ml at ref %s/%s  own choice: klue_ml %s, klue_em_rp %s, rp_ml %s\n"),
              x$cell, as.integer(x$K), x$C_ref, x$C_runner_up, x$bic_margin,
              x$ladder_best_minus_ref, x$repro_diff, x$rp_hits_K, x$rp_starts,
              x$C_bic_klue_ml, x$C_bic_klue_em_rp, x$C_bic_rp_ml))
}
cat("\n===== shortfall at K* and the selected class count, per ladder method =====\n")
cat("(cells whose K* check passed; C_sel assumes the method reaches the reference at every other C)\n")
ok <- tab[tab$repro_ok %in% TRUE, ]
for (m in unique(ok$method)) {
  d <- ok[ok$method == m, ]
  ch <- d$cell[d$could_change %in% TRUE]
  cat(sprintf("%-20s cells %2d  g > 0.5 in %2d  selection could change in %2d  max g %7.2f%s\n",
              m, nrow(d), sum(d$g > 0.5, na.rm = TRUE), length(ch),
              suppressWarnings(max(d$g, na.rm = TRUE)),
              if (length(ch)) paste0("  [", paste(ch, collapse = " "), "]") else ""))
}
above <- cells$cell[!is.na(cells$ladder_best_minus_ref) & cells$ladder_best_minus_ref > HIT_TOL]
if (length(above))
  stamp("NOTE a ladder arm is above the reference at K* in: %s", paste(above, collapse = ", "))
failed <- cells$cell[!(cells$repro_ok %in% TRUE)]
if (length(failed))
  stamp("NOTE K* check failed (excluded above): %s", paste(failed, collapse = ", "))
stamp("wrote %s and %s", OUT_RDS, OUT_CSV)
