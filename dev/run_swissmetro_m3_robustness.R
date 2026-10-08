#!/usr/bin/env Rscript
# =============================================================================
# dev/run_swissmetro_m3_robustness.R
#
# Why. The Swissmetro M3 benchmark (independent normals, negative log-normal
# cost, random constants; car is the reference alternative), estimated once
# from the pooled MNL by dev/rerun_mmnl_benchmark.R under klue 0.10.0, stopped
# at LL -3583.65, 18.8 LL below the klue 0.9.5 fit of the same model with the
# same 3000 MLHS draws (-3564.90, output/mmnl_bench/Swissmetro_M3.rds): a lower
# local maximum. The Swissmetro verdict (continuous) does not depend on it, but
# the margin and the start of M4 do (decision 2026-10-01: check before M4).
#
# Fits (each one apollo_estimate at 3000 MLHS draws, klue 0.10.0, Apollo's
# memorySaver on):
#   from_old   warm start at the 0.9.5 optimum (its log standard deviations
#              converted to standard deviations: the same model and point)
#   zero       neutral start (start = "zero")
#   ref_train  train as the reference alternative (asc_map c(0, 1, 2)),
#              pooled-MNL start: a different model when the constants are
#              independent normals (Apollo manual, J-1 random constants)
#   ref_sm     Swissmetro as the reference (asc_map c(1, 0, 2))
# The benchmark (car reference) keeps the highest-LL fit among the stage file
# and from_old / zero: if one of those is higher, the stage file is kept as
# output/mmnl_bench_std/Swissmetro_M3_pooled_mnl.rds and replaced, so the
# queued M4 starts from the best M3. ref_train / ref_sm are sensitivity fits
# only (they measure how much M3 depends on the reference; M4 does not).
#
# Deterministic: Apollo's draws use its default seed; the pooled MNL start is
# deterministic. Resumable: a variant whose file exists is skipped.
# Run (one core, memorySaver):
#   KLUE_CORES=1 KLUE_MMNL_MEMORY_SAVER=TRUE nice -n 15 \
#     Rscript dev/run_swissmetro_m3_robustness.R
# Writes: output/mmnl_bench_std/robust/Swissmetro_M3_<variant>.rds,
#         output/swissmetro_m3_robustness.csv
# =============================================================================

n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else 1L
options(klue.cores = n_cores, klue.mmnl.n_cores = 1L, mc.cores = n_cores, sm.skip = TRUE)
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
stopifnot(packageVersion("klue") >= "0.10.0")
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "mclust", "cluster"))  # stops on drift
suppressWarnings(suppressMessages(source("R/empirical_swissmetro.R", local = FALSE)))

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
BENCH  <- "output/mmnl_bench_std"
ROBUST <- file.path(BENCH, "robust")
N_DRAWS <- 3000L
# KLUE_SM_ROBUST_TEST=1: a 30-draw dry run into a temporary directory that
# leaves the benchmark files and output/ untouched (for checking the script).
TEST <- identical(Sys.getenv("KLUE_SM_ROBUST_TEST"), "1")
if (TEST) { N_DRAWS <- 30L; ROBUST <- file.path(tempdir(), "sm_robust_test") }
RANDOM  <- c("x1", "x2", "price", "asc1", "asc2")
dir.create(ROBUST, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"),
                                     paste0(...))); flush.console() }

db    <- load_swissmetro_database()
# The pooled-MNL-start fit: the stage file, or its kept copy once replaced.
KEEP  <- file.path(BENCH, "Swissmetro_M3_pooled_mnl.rds")
stage <- readRDS(if (file.exists(KEEP)) KEEP else file.path(BENCH, "Swissmetro_M3.rds"))
old   <- readRDS("output/mmnl_bench/Swissmetro_M3.rds")
stopifnot(identical(stage$random, RANDOM), identical(old$random, RANDOM))

# The 0.9.5 optimum as a klue 0.10.0 warm start: sigma_* were log standard
# deviations, so the same point has sigma = exp(log sd).
old_ws <- old
sg <- grep("^sigma_", names(old_ws$par))
old_ws$par[sg] <- exp(old_ws$par[sg])
old_ws$sd_scale <- "linear"; old_ws$correlation <- FALSE; old_ws$price <- "lognormal"

VARIANTS <- list(
  from_old  = function() klue_mmnl(db, n_draws = N_DRAWS, random = RANDOM, dgp = SM_DGP,
                                   warm_start = old_ws),
  zero      = function() klue_mmnl(db, n_draws = N_DRAWS, random = RANDOM, dgp = SM_DGP,
                                   start = "zero"),
  ref_train = function() klue_mmnl(db, n_draws = N_DRAWS, random = RANDOM,
                                   dgp = klue_dgp(n_generic = 2, n_alternatives = 3,
                                                  asc_map = c(0, 1, 2))),
  ref_sm    = function() klue_mmnl(db, n_draws = N_DRAWS, random = RANDOM,
                                   dgp = klue_dgp(n_generic = 2, n_alternatives = 3,
                                                  asc_map = c(1, 0, 2)))
)

fits <- list()
for (v in names(VARIANTS)) {
  f <- file.path(ROBUST, sprintf("Swissmetro_M3_%s.rds", v))
  if (file.exists(f)) { stamp("SKIP ", v, " (exists)"); fits[[v]] <- readRDS(f); next }
  stamp("START ", v)
  t <- system.time(obj <- VARIANTS[[v]]())[3]
  obj$secs <- unname(t)
  saveRDS(obj, f)
  fits[[v]] <- obj
  stamp(sprintf("DONE  %s  LL=%.3f k=%d BIC=%.1f conv=%s code=%s iter=%s (%.1f min)",
                v, obj$LL, obj$k, obj$BIC, obj$converged,
                obj$apollo_status$code %||% NA, obj$apollo_status$nIter %||% NA, t / 60))
}

row <- function(name, o, reference, role)
  data.frame(fit = name, reference = reference, role = role,
             converged = isTRUE(o$converged), LL = if (isTRUE(o$converged)) o$LL else NA,
             k = o$k, BIC = if (isTRUE(o$converged)) o$BIC else NA,
             start_used = o$start_used %||% NA,
             apollo_code = o$apollo_status$code %||% NA,
             n_iter = o$apollo_status$nIter %||% NA,
             minutes = round((o$secs %||% NA) / 60, 1), stringsAsFactors = FALSE)
tab <- rbind(row("stage (pooled-MNL start)", stage, "car", "benchmark candidate"),
             row("klue 0.9.5 fit", old, "car", "superseded (log-sd); reference point"),
             row("from_old", fits$from_old, "car", "benchmark candidate"),
             row("zero", fits$zero, "car", "benchmark candidate"),
             row("ref_train", fits$ref_train, "train", "sensitivity"),
             row("ref_sm", fits$ref_sm, "Swissmetro", "sensitivity"))

# Keep the best car-reference fit as the benchmark M3 (so M4 starts from it).
cand <- list(stage = stage, from_old = fits$from_old, zero = fits$zero)
ll   <- vapply(cand, function(o) if (isTRUE(o$converged)) o$LL else -Inf, numeric(1))
best <- names(which.max(ll))
tab$selected <- tab$fit == c(stage = "stage (pooled-MNL start)", from_old = "from_old",
                             zero = "zero")[[best]]
if (TEST) {
  stamp("TEST run: best candidate would be '", best, "'; nothing replaced or written")
} else if (best != "stage") {
  if (!file.exists(KEEP)) file.copy(file.path(BENCH, "Swissmetro_M3.rds"), KEEP)
  sel <- cand[[best]]; sel$selected_from <- best
  saveRDS(sel, file.path(BENCH, "Swissmetro_M3.rds"))
  stamp(sprintf("benchmark M3 replaced by '%s' (LL %.3f > %.3f); original kept as %s",
                best, ll[[best]], ll[["stage"]], KEEP))
} else stamp("benchmark M3 unchanged: the pooled-MNL start gives the highest LL")
if (!TEST) write.csv(tab, "output/swissmetro_m3_robustness.csv", row.names = FALSE)
print(tab, row.names = FALSE)
stamp("DONE")
