#!/usr/bin/env Rscript
# =============================================================================
# Regime classification of the five public-dataset applications.
#
# Runs klue::klue(run_mmnl_corr = TRUE) on each dataset
# and classifies the result into one of three regimes:
#   (i)   well-separated discrete: BIC and ICL agree
#   (ii)  overlapping discrete:    BIC keeps falling, ICL turns,
#                                  correlated MMNL does NOT close the LCMNL gap
#   (iii) approximately continuous: BIC keeps falling, ICL turns,
#                                   correlated MMNL closes the LCMNL gap
#
# Outputs:
#   output/regime_classification.rds   -- list of full workflow results per dataset
#   output/regime_classification.csv   -- one row per dataset, the summary line
#   output/regime_classification.done  -- marker file (touched on completion)
#
# Not yet run to completion: the May run stopped inside the correlated MMNL,
# so no result file is shipped (README_REPRODUCE.md).
#
# Run from the repository root, one core (needs test_data_files/data.csv and
# test_data_files/swissmetro/swissmetro.dat; README_REPRODUCE.md, "Data"):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/run_regime_classification.R
# Seeds: none set here. klue() fits klue_lcmnl() at C = 1-6, whose k-means,
# Gaussian-mixture and PAM clustering starts run under seed 123, then the
# independent MMNL and the correlated MMNL started from it, both on 500 MLHS
# draws that Apollo makes under its default apollo_control$seed, 13.
# =============================================================================

# Set skip flags so sourcing each empirical script only loads its helpers,
# not its auto-run block. Each existing empirical script checks its own flag.
options(LCMNL_NO_AUTORUN = TRUE,
        empirical.skip   = TRUE,
        mode.skip        = TRUE,
        swiss.skip       = TRUE,
        elec.skip        = TRUE,
        sm.skip          = TRUE)

if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else suppressMessages(library(klue))

cat("=== REGIME CLASSIFICATION runner ===\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Source each empirical script -- defines its load_*_database() and *_DGP
suppressMessages({
  source("R/empirical_application.R")    # Vittel
  source("R/empirical_mode_choice.R")    # Apollo mode
  source("R/empirical_swiss_route.R")    # Apollo Swiss route
  source("R/empirical_electricity.R")    # mlogit Electricity
  source("R/empirical_swissmetro.R")     # Swissmetro
})

dir.create("output", showWarnings = FALSE)

# Per-dataset specification: a loader + C_cands range + MMNL draw counts
#   (heavier MMNL for larger N*T to keep the apollo integer-overflow at bay
#    we cap n_draws so that N*T*n_draws <= 6.5e6)
datasets <- list(
  vittel = list(
    label  = "Vittel water DCE",
    loader = load_empirical_database,
    Cmax   = 6
  ),
  mode_choice = list(
    label  = "Apollo mode choice",
    loader = load_mode_choice_database,
    Cmax   = 6
  ),
  swiss_route = list(
    label  = "Apollo Swiss route",
    loader = load_swiss_database,
    Cmax   = 6
  ),
  electricity = list(
    label  = "Electricity",
    loader = load_electricity_database,
    Cmax   = 6
  ),
  swissmetro = list(
    label  = "Swissmetro",
    loader = load_swissmetro_database,
    Cmax   = 6
  )
)

# Safe n_draws given N, T (apollo integer-overflow guardrail at ~6.5e6)
safe_n_draws <- function(N, T_tasks, want = 500L) {
  total <- as.numeric(N) * as.numeric(T_tasks) * as.numeric(want)
  if (total > 6.5e6) {
    new_draws <- floor(6.5e6 / (as.numeric(N) * as.numeric(T_tasks)))
    new_draws <- max(100L, as.integer(new_draws))
    cat(sprintf("  N=%d T=%d: n_draws %d -> %d (apollo integer-overflow guardrail)\n",
                N, T_tasks, want, new_draws))
    return(new_draws)
  }
  as.integer(want)
}

# Regime classification logic
classify_regime <- function(summary_df, mmnl_indep_bic, mmnl_corr_bic,
                            Cmax, tol_bic = 5, tol_icl_lag = 1) {
  bic_argmin <- summary_df$C[which.min(summary_df$BIC)]
  icl_argmin <- summary_df$C[which.min(summary_df$ICL)]

  bic_hit_ceiling <- (bic_argmin == Cmax)
  icl_below_bic   <- (icl_argmin < bic_argmin - tol_icl_lag)

  if (!bic_hit_ceiling && !icl_below_bic) {
    return(list(regime = "(i) well-separated discrete",
                bic_argmin = bic_argmin,
                icl_argmin = icl_argmin,
                rationale = "BIC and ICL agree; classes cleanly separable"))
  }

  lcmnl_best_bic <- min(summary_df$BIC, na.rm = TRUE)
  if (is.null(mmnl_corr_bic) || is.na(mmnl_corr_bic)) {
    return(list(regime = "(ii) or (iii) indeterminate",
                bic_argmin = bic_argmin,
                icl_argmin = icl_argmin,
                rationale = "Correlated MMNL did not converge"))
  }
  gap <- mmnl_corr_bic - lcmnl_best_bic
  if (gap <= tol_bic) {
    return(list(regime = "(iii) approximately continuous",
                bic_argmin = bic_argmin,
                icl_argmin = icl_argmin,
                rationale = sprintf("corr. MMNL within %.1f BIC of LCMNL-best (gap %.1f)",
                                    tol_bic, gap)))
  }
  list(regime = "(ii) overlapping discrete",
       bic_argmin = bic_argmin,
       icl_argmin = icl_argmin,
       rationale = sprintf("corr. MMNL falls short of LCMNL-best by %.1f BIC", gap))
}

# Run each dataset
all_results <- list()
summary_rows <- list()

for (key in names(datasets)) {
  spec  <- datasets[[key]]
  label <- spec$label

  cat(sprintf("\n=================== %s ===================\n", label))
  t0 <- Sys.time()

  database <- tryCatch(spec$loader(), error = function(e) {
    cat("  loader failed:", conditionMessage(e), "\n"); NULL
  })
  if (is.null(database)) next

  N <- length(unique(database$ID))
  T_tasks <- length(unique(database$TASK))
  cat(sprintf("  N=%d, T=%d\n", N, T_tasks))

  n_draws <- safe_n_draws(N, T_tasks, want = 500L)

  res <- tryCatch(
    klue(
      database      = database,
      C_cands       = 1:spec$Cmax,
      run_mmnl      = TRUE,
      run_mmnl_corr = TRUE,
      mmnl_opts     = list(n_draws = n_draws),
      write_csv     = FALSE,
      verbose       = TRUE
    ),
    error = function(e) {
      cat("  workflow failed:", conditionMessage(e), "\n"); NULL
    }
  )
  if (is.null(res)) next

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("\n  Workflow elapsed: %.1f minutes\n", dt / 60))

  mmnl_bic       <- if (isTRUE(res$mmnl$converged))      res$mmnl$BIC      else NA_real_
  mmnl_corr_bic  <- if (isTRUE(res$mmnl_corr$converged)) res$mmnl_corr$BIC else NA_real_

  reg <- classify_regime(res$summary, mmnl_bic, mmnl_corr_bic, spec$Cmax)

  cat(sprintf("\n  LCMNL-best C=%d  BIC=%.1f  |  ICL-best C=%d  ICL=%.1f\n",
              reg$bic_argmin, min(res$summary$BIC),
              reg$icl_argmin, min(res$summary$ICL)))
  cat(sprintf("  MMNL (indep) BIC=%.1f  |  MMNL (corr) BIC=%.1f\n",
              mmnl_bic, mmnl_corr_bic))
  cat(sprintf("  REGIME: %s  (%s)\n", reg$regime, reg$rationale))

  all_results[[key]] <- list(
    label = label, N = N, T = T_tasks, n_draws = n_draws,
    summary_df = res$summary,
    mmnl_indep_bic = mmnl_bic,
    mmnl_corr_bic  = mmnl_corr_bic,
    bic_argmin = reg$bic_argmin,
    icl_argmin = reg$icl_argmin,
    regime = reg$regime,
    rationale = reg$rationale,
    full_result = res,
    elapsed_minutes = dt / 60
  )

  summary_rows[[key]] <- data.frame(
    dataset            = label,
    N                  = N,
    T                  = T_tasks,
    n_draws            = n_draws,
    lcmnl_best_C       = reg$bic_argmin,
    lcmnl_best_BIC     = min(res$summary$BIC),
    icl_best_C         = reg$icl_argmin,
    icl_best_ICL       = min(res$summary$ICL),
    mmnl_indep_BIC     = mmnl_bic,
    mmnl_corr_BIC      = mmnl_corr_bic,
    regime             = reg$regime,
    rationale          = reg$rationale,
    elapsed_minutes    = round(dt / 60, 1),
    stringsAsFactors   = FALSE
  )

  # Incremental save: write after each dataset in case of failure later
  saveRDS(all_results, "output/regime_classification.rds")
  write.csv(do.call(rbind, summary_rows),
            "output/regime_classification.csv", row.names = FALSE)
}

cat("\n========================== SUMMARY ==========================\n")
print(do.call(rbind, summary_rows), row.names = FALSE)

cat("\nEnd:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
file.create("output/regime_classification.done")
cat("\nWritten: output/regime_classification.rds and .csv\n")
