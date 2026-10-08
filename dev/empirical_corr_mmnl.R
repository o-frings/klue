# dev/empirical_corr_mmnl.R  (WS3)
# SUPERSEDED for version12 by dev/rerun_mmnl_benchmark.R (models M2, M4), which
# starts each correlated fit exactly at its independent solution. Kept because
# it produced the version11 correlated-MMNL numbers; under klue >= 0.9.5 it no
# longer reproduces them (0.9.5 fixed the correlated starting values).
# Correlated MMNL on the five public datasets, alongside the independent-normal
# MMNL already reported in the paper. Answers the reviewer objection that the
# "genuinely discrete" verdicts (esp. Swissmetro) rest on a too-thin continuous
# benchmark: a full Cholesky-correlated MMNL is the richer continuous model.
#
# Reuses each empirical script's data loader + DGP (autorun disabled via its
# skip option). Run from repo root, prefixing each run with
# OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 15
# (Apollo's nCores is KLUE_CORES, default 2):
#   Rscript dev/empirical_corr_mmnl.R            # full: 3000 MLHS draws
#   Rscript dev/empirical_corr_mmnl.R 200 Mode   # smoke: 200 draws, Mode only
# Writes: output/empirical_corr_mmnl_<draws>draws.csv, one row per dataset,
#   saved after each (output/empirical_corr_mmnl_3000draws.csv for the full
#   run). Superseded, so not rerun for version12.
# Seeds: Apollo makes the MLHS draws under its default apollo_control$seed
#   (13); nothing else draws random numbers (the independent MMNL starts from
#   the pooled MNL, the correlated one from the independent fit).

# From the repo, load the dev source (latest code); else the installed package.
if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else {
  suppressWarnings(suppressMessages(library(klue)))
}
if (!dir.exists("output")) dir.create("output")

# Thermal cap: keep Apollo's MMNL parallelism low so the laptop stays cool/usable.
# Override with env var KLUE_CORES.
options(klue.mmnl.n_cores = max(1L, as.integer(Sys.getenv("KLUE_CORES", "2"))))

SPECS <- list(
  Mode        = list(script = "R/empirical_mode_choice.R", skip = "mode.skip",
                     loader = "load_mode_choice_database", dgp = "MODE_DGP"),
  SwissRoute  = list(script = "R/empirical_swiss_route.R", skip = "swiss.skip",
                     loader = "load_swiss_database",       dgp = "SWISS_DGP"),
  Electricity = list(script = "R/empirical_electricity.R", skip = "elec.skip",
                     loader = "load_electricity_database", dgp = "ELEC_DGP"),
  Swissmetro  = list(script = "R/empirical_swissmetro.R",  skip = "sm.skip",
                     loader = "load_swissmetro_database",  dgp = "SM_DGP"),
  Vittel      = list(script = "R/empirical_application.R", skip = "empirical.skip",
                     loader = "load_empirical_database",   dgp = "EMP_DGP")
)

run_one <- function(name, n_draws = 3000L) {
  s <- SPECS[[name]]
  options(structure(list(TRUE), names = s$skip))   # disable that script's autorun
  source(s$script, local = FALSE)                  # defines loader + DGP globally
  db  <- get(s$loader)()
  dgp <- get(s$dgp)
  cat(sprintf("\n--- %s: indep + correlated MMNL (%d draws) ---\n", name, n_draws))
  mi <- klue_mmnl(db, n_draws = n_draws, dgp = dgp)
  mc <- klue_mmnl_corr(db, n_draws = n_draws, dgp = dgp,
                       warm_start = if (isTRUE(mi$converged)) mi)
  data.frame(
    dataset    = name,
    k_indep    = mi$k, BIC_indep = round(mi$BIC, 1),
    k_corr     = mc$k, BIC_corr  = round(mc$BIC, 1),
    dBIC_corr  = round(mc$BIC - mi$BIC, 1),   # negative => correlated MMNL fits better
    conv       = mi$converged && mc$converged,
    stringsAsFactors = FALSE)
}

# Autorun only when invoked directly (Rscript), not when sourced for its
# functions (e.g. by dev/run_revision_batch.R).
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_draws <- if (length(args) >= 1) as.integer(args[1]) else 3000L
  which    <- if (length(args) >= 2) args[-1] else names(SPECS)

  out <- file.path("output", sprintf("empirical_corr_mmnl_%ddraws.csv", n_draws))
  # Incremental + resumable: process datasets one at a time, append to the CSV
  # after each, and skip any already present (so a hang on one dataset cannot
  # lose the others, and a re-run continues where it stopped).
  done <- if (file.exists(out)) read.csv(out, stringsAsFactors = FALSE) else NULL
  res <- done
  for (nm in which) {
    if (!is.null(done) && nm %in% done$dataset) {
      cat(sprintf("SKIP %s (already in %s)\n", nm, out)); next
    }
    row <- run_one(nm, n_draws = n_draws)
    res <- rbind(res, row)
    write.csv(res, out, row.names = FALSE)   # save after every dataset
    cat(sprintf(">> %s: BIC_indep=%.1f BIC_corr=%.1f dBIC_corr=%+.1f conv=%s [saved]\n",
                nm, row$BIC_indep, row$BIC_corr, row$dBIC_corr, row$conv))
  }
  cat("\n==== CORRELATED vs INDEPENDENT MMNL (empirical) ====\n")
  print(res, row.names = FALSE)
  cat(sprintf("\nwritten: %s\n", out))
}
