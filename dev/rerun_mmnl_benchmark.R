#!/usr/bin/env Rscript
# =============================================================================
# MMNL benchmark re-run on the five public datasets (klue 0.10.0, standard
# Apollo estimation; first version klue 0.9.5, 2026-09-26).
#
# Why. The published benchmark fixed every constant while the LCMNL estimates
# them per class, so the LCMNL-versus-MMNL comparison partly measured random
# versus fixed constants. Two further defects: the correlated fits started from
# a mis-mapped point (log standard deviations on the linear Cholesky diagonal,
# constants reset to zero), which left the Swissmetro correlated MMNL 297 LL
# below the independent model it nests; and Vittel was estimated under a
# rank-deficient coding (see R/empirical_application.R). This script fits four
# benchmarks per dataset, Vittel under the identified coding:
#   M1  independent normals, lognormal price, fixed constants (the published one)
#   M2  M1 plus a full Cholesky covariance, warm-started from M1
#   M3  M1 with the constants random too
#   M4  M3 plus a full Cholesky covariance, warm-started from M3
# Every fit is one apollo_estimate call at 3000 MLHS draws, as in Apollo's
# example scripts (klue 0.10.0: sigma unconstrained, no coarse-draw stage, no
# fallbacks, convergence = Apollo's successfulEstimation). M1 and M3 start from
# the pooled MNL; M2 and M4 start where they coincide with M1 and M3 (Apollo's
# apollo_readBeta step, same draws), so they cannot end below them; `nested_ok`
# records the check. Vittel also extends its identified LCMNL ladder to
# C = 10 (vittel_respec_lcmnl.rds stops where BIC first rose, at C = 8, by 0.4).
#
# Stages are resumable: each fit is saved to
# output/mmnl_bench_std/<dataset>_<M>.rds and skipped if present (delete a file
# to redo it). output/mmnl_bench/ holds the superseded klue 0.9.5 fits (log-sd,
# two-stage) and is never read here. Deterministic: fixed MLHS draws (Apollo's
# default seed) and a seeded LCMNL multistart.
#
# Run (the machine's cap is 2 cores; Apollo stays single-core when detached,
# since a 2-core detached correlated fit died with a bus error on 2026-09-21).
# Two lanes in parallel, one core each:
#   KLUE_CORES=1 nohup nice -n 15 Rscript dev/rerun_mmnl_benchmark.R Vittel \
#     > output/mmnl_bench_vittel.log 2>&1 &
#   KLUE_CORES=1 nohup nice -n 15 Rscript dev/rerun_mmnl_benchmark.R \
#     Swissmetro Mode SwissRoute Electricity > output/mmnl_bench_others.log 2>&1 &
#   Rscript dev/rerun_mmnl_benchmark.R summary     # collect whatever exists
# Writes: output/mmnl_bench_std/*.rds, output/mmnl_benchmark_v12.csv,
#         output/vittel_respec_lcmnl_ext.csv
# =============================================================================

n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else 1L
options(klue.cores = n_cores, klue.mmnl.n_cores = 1L, mc.cores = n_cores)
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
stopifnot(packageVersion("klue") >= "0.10.0")
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "mclust", "cluster"))  # stops on drift

OUT      <- "output"
BENCH    <- file.path(OUT, "mmnl_bench_std")
N_DRAWS  <- 3000L
dir.create(BENCH, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     paste0(...))); flush.console() }

SPECS <- list(
  Swissmetro  = list(script = "R/empirical_swissmetro.R",  skip = "sm.skip",
                     loader = "load_swissmetro_database",  dgp = "SM_DGP"),
  Mode        = list(script = "R/empirical_mode_choice.R", skip = "mode.skip",
                     loader = "load_mode_choice_database", dgp = "MODE_DGP"),
  SwissRoute  = list(script = "R/empirical_swiss_route.R", skip = "swiss.skip",
                     loader = "load_swiss_database",       dgp = "SWISS_DGP"),
  Electricity = list(script = "R/empirical_electricity.R", skip = "elec.skip",
                     loader = "load_electricity_database", dgp = "ELEC_DGP"),
  Vittel      = list(script = "R/empirical_application.R", skip = "empirical.skip",
                     loader = "load_empirical_database",   dgp = "EMP_DGP")
)
MODELS <- c(M1 = "independent, fixed constants",   M2 = "correlated, fixed constants",
            M3 = "independent, random constants",  M4 = "correlated, random constants")

fit_stage <- function(name, model, fn) {
  f <- file.path(BENCH, sprintf("%s_%s.rds", name, model))
  if (file.exists(f)) { stamp("SKIP ", name, " ", model, " (exists)"); return(readRDS(f)) }
  stamp("START ", name, " ", model, " (", MODELS[[model]], ")")
  t <- system.time(obj <- fn())[3]
  obj$secs <- unname(t)
  saveRDS(obj, f)
  stamp(sprintf("DONE  %s %s  LL=%.2f k=%d BIC=%.1f conv=%s nested_ok=%s (%.1f min)",
                name, model, obj$LL, obj$k, obj$BIC, obj$converged,
                if (is.null(obj$nested_ok)) "-" else obj$nested_ok, t / 60))
  obj
}

run_dataset <- function(name) {
  s <- SPECS[[name]]
  options(structure(list(TRUE), names = s$skip))   # load the adapter, skip its autorun
  suppressWarnings(suppressMessages(source(s$script, local = FALSE)))
  db  <- get(s$loader)()
  dgp <- get(s$dgp)
  fixed_r <- c(paste0("x", seq_len(dgp$n_generic)), "price")
  rand_r  <- c(fixed_r, paste0("asc", seq_len(dgp$n_asc)))
  stamp(sprintf("%s: N=%d, %d generic + price, %d constants", name,
                length(unique(db$ID)), dgp$n_generic, dgp$n_asc))
  mm <- function(...) klue_mmnl(db, n_draws = N_DRAWS, dgp = dgp, ...)
  m1 <- fit_stage(name, "M1", function() mm(random = fixed_r))
  if (m1$converged) fit_stage(name, "M2", function()
    mm(random = fixed_r, correlation = TRUE, warm_start = m1))
  m3 <- fit_stage(name, "M3", function() mm(random = rand_r))
  if (m3$converged) fit_stage(name, "M4", function()
    mm(random = rand_r, correlation = TRUE, warm_start = m3))

  if (identical(name, "Vittel")) {
    # Extend the identified LCMNL ladder to C = 10 (seeded multistart per C).
    f <- file.path(OUT, "vittel_respec_lcmnl_ext.csv")
    if (!file.exists(f)) {
      base <- readRDS(file.path(OUT, "vittel_respec_lcmnl.rds"))$table
      rows <- list()
      for (C in setdiff(9:10, base$C)) {
        t0 <- Sys.time()
        m  <- klue_lcmnl(db, C, dgp = dgp)
        rows[[length(rows) + 1]] <- data.frame(
          C = C, converged = m$converged, LL = m$LL, k = m$k, BIC = m$BIC,
          AIC = m$AIC, ICL = m$ICL, ICL_BIC = m$ICL_BIC,
          min_class_prob = min(m$class_probs), singular = isTRUE(m$singular),
          best_method = if (is.null(m$best_method)) "pooled" else m$best_method,
          secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
          stringsAsFactors = FALSE)
        stamp(sprintf("  Vittel LCMNL C=%d  LL=%.2f k=%d BIC=%.1f", C, m$LL, m$k, m$BIC))
      }
      write.csv(rbind(base, do.call(rbind, rows)), f, row.names = FALSE)
    }
  }
}

summarise <- function() {
  # Stage files only: kept copies such as Swissmetro_M3_pooled_mnl.rds (written
  # by the second-start checks) are not benchmark rows.
  files <- list.files(BENCH, pattern = "^[A-Za-z]+_M[1-4]\\.rds$", full.names = TRUE)
  if (!length(files)) return(invisible(NULL))
  rows <- lapply(files, function(f) {
    o  <- readRDS(f)
    id <- sub("\\.rds$", "", basename(f))
    data.frame(dataset = sub("_M[1-4]$", "", id), model = sub("^.*_", "", id),
               description = MODELS[[sub("^.*_", "", id)]],
               converged = o$converged, LL = round(o$LL, 2), k = o$k,
               BIC = round(o$BIC, 1),
               nested_ok = if (is.null(o$nested_ok)) NA else o$nested_ok,
               apollo_code = if (is.null(o$apollo_status)) NA else o$apollo_status$code,
               apollo_message = if (is.null(o$apollo_status)) NA else o$apollo_status$message,
               vcov_available = !is.null(o$vcov) && !anyNA(o$vcov),
               minutes = round(o$secs / 60, 1), stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$dataset, names(SPECS)), tab$model), ]
  write.csv(tab, file.path(OUT, "mmnl_benchmark_v12.csv"), row.names = FALSE)
  print(tab, row.names = FALSE)
  invisible(tab)
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- names(SPECS)
for (a in setdiff(args, "summary")) {
  if (!a %in% names(SPECS)) stop("unknown dataset: ", a)
  run_dataset(a)
}
summarise()
stamp("DONE: ", paste(args, collapse = ", "))
