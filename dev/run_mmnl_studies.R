#!/usr/bin/env Rscript
# =============================================================================
# Re-run the two randomised-design MMNL studies of the paper with klue >= 0.10.0
# (standard Apollo estimation: one apollo_estimate call at 3000 MLHS draws,
# sigma unconstrained, convergence = Apollo's successfulEstimation, correlated
# fits started from the independent estimates).
#
#   mmnl       klue_study_mmnl()       -> output/mmnl_results.csv        (tab:mmnl)
#   mmnl_corr  klue_study_mmnl_corr()  -> output/mmnl_correlated_results.csv (tab:mmnl_corr)
#
# Deterministic: per-condition seeds (studies/klue_studies.R, .cond_seed) and
# Apollo's fixed MLHS draws. Not resumable within a study: each study writes its
# CSV at the end.
# Cores: KLUE_CORES (default 1) for the LCMNL phase; Apollo runs on one core.
#
# Guards (added 2026-10-03, so that a second launch, e.g. by the rerun queue
# while a parallel run is under way, costs nothing):
#   - a study whose CSV was already written by this klue version is skipped
#     (marker output/<csv>.done, first line = klue version); KLUE_FORCE=1 redoes it;
#   - while another instance is running (lock output/.run_mmnl_studies.lock
#     holding a live pid), a new instance exits at once.
# Both checks run before klue is loaded, so a skipped launch is instant.
#
# Run: KLUE_CORES=1 nice -n 15 Rscript dev/run_mmnl_studies.R mmnl mmnl_corr
# =============================================================================

OUT  <- "output"   # = klue's OUTPUT_DIR (klue/R/00_config.R), asserted below
STUDY_CSV <- c(mmnl = "mmnl_results.csv", mmnl_corr = "mmnl_correlated_results.csv")
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- names(STUDY_CSV)
bad <- setdiff(args, names(STUDY_CSV))
if (length(bad)) stop("unknown study: ", bad[1])
stamp <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M"), paste0(...)))

KLUE_VER  <- unname(read.dcf("klue/DESCRIPTION", fields = "Version")[1, 1])
done_file <- function(a) file.path(OUT, paste0(STUDY_CSV[[a]], ".done"))
is_done   <- function(a) file.exists(done_file(a)) &&
  identical(readLines(done_file(a), n = 1, warn = FALSE), KLUE_VER)
if (!identical(Sys.getenv("KLUE_FORCE"), "1")) {
  skip <- args[vapply(args, is_done, logical(1))]
  for (a in skip) stamp("skip ", a, ": output/", STUDY_CSV[[a]], " already written by klue ",
                        KLUE_VER, " (KLUE_FORCE=1 to redo)")
  args <- setdiff(args, skip)
}
if (!length(args)) { stamp("nothing to do"); quit(save = "no", status = 0) }

LOCK <- file.path(OUT, ".run_mmnl_studies.lock")
if (file.exists(LOCK)) {
  pid <- suppressWarnings(as.integer(readLines(LOCK, n = 1, warn = FALSE)))
  if (!is.na(pid) && pid != Sys.getpid() && isTRUE(tools::pskill(pid, 0L))) {
    stamp("another run_mmnl_studies.R is running (pid ", pid, "); exiting")
    quit(save = "no", status = 0)
  }
}
writeLines(as.character(Sys.getpid()), LOCK)

n_cores <- if (nzchar(Sys.getenv("KLUE_CORES"))) as.integer(Sys.getenv("KLUE_CORES")) else 1L
options(klue.cores = n_cores, klue.mmnl.n_cores = 1L, mc.cores = n_cores)
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
stopifnot(packageVersion("klue") >= "0.10.0")
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "mclust", "cluster"))  # stops on drift
source("studies/klue_studies.R")
if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- OUT
stopifnot(identical(normalizePath(OUTPUT_DIR), normalizePath(OUT)))

STUDIES <- list(mmnl = klue_study_mmnl, mmnl_corr = klue_study_mmnl_corr)
tryCatch({
  for (a in args) {
    t0 <- Sys.time()
    res <- STUDIES[[a]](verbose = TRUE, n_cores = n_cores)
    write.csv(res, file.path(OUT, STUDY_CSV[[a]]), row.names = FALSE)
    writeLines(c(KLUE_VER, format(Sys.time(), "%Y-%m-%d %H:%M:%S")), done_file(a))
    stamp(sprintf("%s done in %.1f h -> output/%s", a,
                  as.numeric(difftime(Sys.time(), t0, units = "hours")), STUDY_CSV[[a]]))
  }
}, finally = unlink(LOCK))
