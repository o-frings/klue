# dev/record_versions.R
# ============================================================================
# Records the software versions behind the paper's results. Writes
#   output/package_versions.csv  one row per package the analyses use: version
#                                in use, library it comes from, where it is
#                                used, the pinned version (dev/pinned_versions.R)
#                                and whether they match
#   output/session_info.txt      sessionInfo() with those packages loaded
#   renv.lock                    the full dependency closure (renv 1.x format),
#                                at the repository root
# The gmnl arm needs mlogit 1.0-3.1 from a side library (dev/compare_packages.R);
# the CSV lists that copy next to the main-library mlogit, which only supplies
# the Electricity data. renv.lock holds the main-library copy, because a
# lockfile has one version per package. klue is loaded from this repository's
# klue/ tree, so it is listed in the CSV but not in the lockfile.
#
# Apart from klue and the side-library mlogit 1.0-3.1 (the same version,
# reinstalled on 2026-10-06; the earlier copy was in /tmp/oldmlogit_lib), the R
# library behind the results has not changed since 2025-07-10 (installation
# dates in each package's DESCRIPTION), and R 4.3.1 is the only R on the
# machine, so one record covers every run in the paper.
#
# Deterministic; no RNG; no estimation. About 10 s.
# Run from the repo root:  Rscript dev/record_versions.R
# ============================================================================

if (!file.exists("dev/record_versions.R")) stop("run from the repository root")
OLD_MLOGIT_LIB <- path.expand(Sys.getenv("KLUE_OLD_MLOGIT_LIB", "~/R/oldmlogit_lib"))
source("dev/pinned_versions.R")

ROLES <- c(
  klue        = "every analysis (loaded from klue/ with pkgload::load_all)",
  apollo      = "klue_mmnl() and every Apollo arm; the mode-choice and Swiss route data",
  gmnl        = "gmnl_lc arm of both stress ladders",
  mlogit      = "Electricity data (main library); gmnl's data preparation (side library)",
  idefix      = "klue_design(): the blocked D-efficient design",
  mclust      = "klue's Gaussian-mixture clustering start",
  cluster     = "klue's PAM clustering start",
  logistf     = "Firth regressions of enumeration accuracy (dev/run_firth_h3.R)",
  maxLik      = "inside Apollo: BFGS in apollo_searchStart and apollo_lcEM",
  bgw         = "inside Apollo: the BGW optimiser of apollo_estimate",
  randtoolbox = "inside Apollo: MLHS and Halton draws",
  numDeriv    = "inside Apollo; dev/check_lcmnl_standard_errors.R",
  pkgload     = "loads klue from the source tree in every script",
  testthat    = "klue's unit tests only"
)

desc_field <- function(pkg, field, lib = NULL) {
  d <- tryCatch(utils::packageDescription(pkg, lib.loc = lib, fields = field),
                error = function(e) NA)
  if (is.null(d) || length(d) == 0) NA_character_ else as.character(d)
}

rows <- list()
add <- function(pkg, version, library, role) {
  pin <- if (pkg %in% names(PINNED_VERSIONS)) PINNED_VERSIONS[[pkg]] else NA_character_
  rows[[length(rows) + 1L]] <<- data.frame(
    package = pkg, version = version, library = library, used_for = role,
    pinned = pin,
    matches_pin = if (is.na(pin) || is.na(version)) NA else
      package_version(version) == package_version(pin),
    stringsAsFactors = FALSE)
}

add("R", paste(R.version$major, R.version$minor, sep = "."), R.home(),
    paste("all analyses;", R.version$platform))
rows[[1]]$pinned <- PINNED_R; rows[[1]]$matches_pin <- rows[[1]]$version == PINNED_R
add("klue", read.dcf("klue/DESCRIPTION", fields = "Version")[1, 1],
    "repository klue/ tree", ROLES[["klue"]])
main_lib <- .libPaths()[1]
for (p in setdiff(names(ROLES), "klue")) {
  v <- desc_field(p, "Version", main_lib)
  role <- if (p == "mlogit") "Electricity data" else ROLES[[p]]
  if (p == "mlogit") {
    # the main-library copy is not the one gmnl uses; its pin is the side copy
    rows[[length(rows) + 1L]] <- data.frame(
      package = p, version = v, library = main_lib, used_for = role,
      pinned = NA_character_, matches_pin = NA, stringsAsFactors = FALSE)
  } else add(p, v, main_lib, role)
}
add("mlogit", desc_field("mlogit", "Version", OLD_MLOGIT_LIB), OLD_MLOGIT_LIB,
    "gmnl's data preparation (gmnl 1.1-3.2 fails on mlogit >= 1.1)")
tab <- do.call(rbind, rows)
# Installation date = modification time of the library copy's DESCRIPTION.
tab$installed <- vapply(seq_len(nrow(tab)), function(i) {
  if (tab$package[i] %in% c("R", "klue")) return(NA_character_)
  d <- file.path(tab$library[i], tab$package[i], "DESCRIPTION")
  if (file.exists(d)) format(as.Date(file.info(d)$mtime)) else NA_character_
}, character(1))
# write the library paths with ~ for the home directory
tab$library <- sub(path.expand("~"), "~", tab$library, fixed = TRUE)
write.csv(tab, "output/package_versions.csv", row.names = FALSE)
print(tab[, c("package", "version", "pinned", "matches_pin")], row.names = FALSE)

# Session info with the analysis packages loaded the way the scripts load them
# (side-library mlogit first, so gmnl sees it).
.libPaths(c(OLD_MLOGIT_LIB, .libPaths()))
suppressWarnings(suppressMessages({
  pkgload::load_all("klue", quiet = TRUE)
  for (p in c("apollo", "gmnl", "idefix", "mclust", "cluster", "logistf"))
    requireNamespace(p, quietly = TRUE)
}))
writeLines(c(paste("Recorded", format(Sys.Date()), "by dev/record_versions.R"), "",
             utils::capture.output(print(utils::sessionInfo()))),
           "output/session_info.txt")
.libPaths(setdiff(.libPaths(), OLD_MLOGIT_LIB))

# Lockfile of the dependency closure, from the main library. renv writes it
# without activating a project; RENV_PATHS_ROOT goes to a temporary folder so
# nothing is cached in the user's renv root.
Sys.setenv(RENV_PATHS_ROOT = file.path(tempdir(), "renv_root"))
renv::lockfile_create(
  packages = c("apollo", "gmnl", "mlogit", "idefix", "mclust", "cluster",
               "logistf", "numDeriv", "pkgload", "testthat"),
  libpaths = main_lib, project = getwd(), prompt = FALSE) |>
  renv::lockfile_write(file = "renv.lock")
lock <- renv::lockfile_read("renv.lock")
cat(sprintf("\nrenv.lock: R %s, %d packages\n", lock$R$Version, length(lock$Packages)))
bad <- tab[!is.na(tab$matches_pin) & !tab$matches_pin, ]
if (nrow(bad)) {
  print(bad)
  stop("versions in use differ from dev/pinned_versions.R")
}
cat("all pinned versions match; wrote output/package_versions.csv, output/session_info.txt, renv.lock\n")
