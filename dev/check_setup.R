# =============================================================================
# Readiness check for the parked heavy reruns. Confirms (a) package
# dependencies, (b) klue is loadable and the key estimators exist, (c) the raw
# data files are present, and (d) each empirical loader builds its database.
# Estimates NOTHING -- no heavy compute. A green run means the parked scripts
# will load and reach their estimation step; convergence on large data / high C
# is only verified by an actual run.
#   Run:  Rscript dev/check_setup.R
# =============================================================================

ok_pkg <- function(p) isTRUE(requireNamespace(p, quietly = TRUE))

cat("== dependencies ==\n")
deps <- c("pkgload", "apollo", "idefix", "mclust", "cluster", "mlogit",
          "logistf", "gmnl")
for (p in deps) cat(sprintf("  %-9s %s\n", p, ok_pkg(p)))

cat("== klue + key functions ==\n")
tryCatch(suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE))),
         error = function(e) cat("  load_all error:", conditionMessage(e), "\n"))
for (f in c("klue_mmnl", "klue_mmnl_corr",
            "klue_lcmnl", "klue_design"))
  cat(sprintf("  %-26s %s\n", f, exists(f)))

cat("== data files ==\n")
cat(sprintf("  Vittel data.csv : %s\n", file.exists("test_data_files/data.csv")))
sm <- c("test_data_files/swissmetro/swissmetro.dat", "test_data_files/swissmetro.dat")
cat(sprintf("  Swissmetro .dat : %s\n", any(file.exists(sm))))

cat("== empirical loaders (build database only; no estimation) ==\n")
SPECS <- list(
  Mode        = c("R/empirical_mode_choice.R", "mode.skip",      "load_mode_choice_database"),
  SwissRoute  = c("R/empirical_swiss_route.R", "swiss.skip",     "load_swiss_database"),
  Electricity = c("R/empirical_electricity.R", "elec.skip",      "load_electricity_database"),
  Swissmetro  = c("R/empirical_swissmetro.R",  "sm.skip",        "load_swissmetro_database"),
  Vittel      = c("R/empirical_application.R", "empirical.skip", "load_empirical_database")
)
for (nm in names(SPECS)) {
  s <- SPECS[[nm]]
  res <- tryCatch({
    options(structure(list(TRUE), names = s[2]))   # disable that script's autorun
    source(s[1], local = FALSE)
    db <- get(s[3])()
    sprintf("OK   rows=%d  cols=%d  respondents=%d",
            nrow(db), ncol(db), length(unique(db$ID)))
  }, error = function(e) paste("ERR ", conditionMessage(e)))
  cat(sprintf("  %-12s %s\n", nm, res))
}
