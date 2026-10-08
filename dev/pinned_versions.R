# dev/pinned_versions.R
# ============================================================================
# The package versions that produced the paper's results, and a check that the
# running library matches them. Sourced by the comparison harness
# (dev/compare_packages.R, hence by both stress ladders, the Apollo isolation
# arms, the benchmark and the reference-optimum check) and by the MMNL runners.
#
# A mismatch stops the run. CRAN now serves newer apollo, gmnl, mlogit and
# mclust releases, and apollo_searchStart switched its optimiser to BGW in
# Apollo 0.3.8, so a current library does not reproduce the tables. For
# exploratory runs only, KLUE_ALLOW_VERSION_DRIFT=1 turns the stop into a
# warning. A different R version only warns.
#
# Install the pinned versions (README_REPRODUCE.md, "Package versions"):
#   remotes::install_version("apollo", "0.3.5")
#   remotes::install_version("gmnl",   "1.1-3.2")
#   remotes::install_version("mlogit", "1.0-3.1", lib = "~/R/oldmlogit_lib",
#                            dependencies = FALSE)   # side library, gmnl only
# klue itself is loaded from this repository's klue/ tree (0.10.0).
# dev/record_versions.R writes the full list to output/package_versions.csv
# and renv.lock. Deterministic; no RNG.
# ============================================================================

PINNED_R <- "4.3.1"
PINNED_VERSIONS <- c(
  apollo  = "0.3.5",    # klue_mmnl() and every Apollo arm
  gmnl    = "1.1-3.2",  # gmnl_lc arm of both stress ladders
  mlogit  = "1.0-3.1",  # gmnl's data preparation only, from the side library
  idefix  = "1.1.0",    # klue_design(): the blocked D-efficient design
  mclust  = "6.1.1",    # klue's Gaussian-mixture clustering start
  cluster = "2.1.4",    # klue's PAM clustering start
  logistf = "1.26.1"    # the Firth regressions (dev/run_firth_h3.R)
)

# Version actually in use: the loaded namespace if there is one (so a side
# library put first on .libPaths() is what counts), else the first library
# copy on the path.
.version_in_use <- function(pkg) {
  if (isNamespaceLoaded(pkg)) return(as.character(getNamespaceVersion(pkg)))
  tryCatch(as.character(utils::packageVersion(pkg)),
           error = function(e) NA_character_)
}

check_pinned_versions <- function(pkgs) {
  r_have <- paste(R.version$major, R.version$minor, sep = ".")
  if (r_have != PINNED_R)
    warning(sprintf("R %s; the paper's results were produced with R %s.",
                    r_have, PINNED_R), call. = FALSE)
  bad <- character(0)
  for (p in pkgs) {
    want <- PINNED_VERSIONS[[p]]
    have <- .version_in_use(p)
    if (is.na(have) || package_version(have) != package_version(want))
      bad <- c(bad, sprintf("%s %s (pinned %s)", p,
                            if (is.na(have)) "not installed" else have, want))
  }
  if (length(bad)) {
    msg <- paste0("package versions differ from the pinned ones: ",
                  paste(bad, collapse = "; "),
                  ". See dev/pinned_versions.R; KLUE_ALLOW_VERSION_DRIFT=1 runs anyway.")
    if (identical(Sys.getenv("KLUE_ALLOW_VERSION_DRIFT"), "1"))
      warning(msg, call. = FALSE) else stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}
