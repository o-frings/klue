# =============================================================================
# FIFTH EMPIRICAL APPLICATION: Swissmetro (Bierlaire et al. 2001)
#
# Stated-preference mode-choice data widely used as a benchmark for discrete
# choice modelling. Sfeir et al. (2021) report LCMNL optimal at K = 5 on this
# dataset using a different specification and sample treatment.
#
# Standard Biogeme sample: commute and business trips (PURPOSE 1 and 3) with
# a recorded choice, 752 respondents x 9 tasks. Unavailable alternatives (e.g.
# no car) are masked out of the choice set (av_* columns) rather than the
# tasks dropped.
#
# Specification (2 generic + cost + 2 ASCs per class):
#   asc_alt1  Train intercept             (car = reference, alt 3)
#   asc_alt2  Swissmetro intercept
#   x1        time     (minutes / 100)
#   x2        headway  (minutes / 100; 0 for train and car)
#   price     cost     (CHF / 100)
#
# Run from the repository root, one core (the data file goes in
# test_data_files/swissmetro/swissmetro.dat; README_REPRODUCE.md, "Data"):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/empirical_swissmetro.R
# With options(sm.skip = TRUE), sourcing it only defines the loader and the
# helpers (dev/rerun_mmnl_benchmark.R does so).
# Writes: output/swissmetro_lcmnl_results.csv (C = 1-6) and
#         output/swissmetro_model_comparison.csv (MNL, BIC-best LCMNL and the
#         independent MMNL; the MMNL row is superseded by
#         dev/rerun_mmnl_benchmark.R and was not rerun).
# Seeds: none set here. klue_lcmnl() runs its k-means, Gaussian-mixture and
# PAM clustering starts under seed 123 and restores the caller's random
# stream; klue_mmnl() leaves Apollo's apollo_control$seed at its default, 13,
# under which Apollo makes the 3000 MLHS draws.
# =============================================================================

options(LCMNL_NO_AUTORUN = TRUE)
# Canonical engine is the klue package. From the repo, load the dev source
# (latest code); otherwise use the installed package.
if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else {
  suppressWarnings(suppressMessages(library(klue)))
}
if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- "output"

SM_DGP <- klue_dgp(n_generic = 2, n_alternatives = 3)

# ---- 1. Build database -----------------------------------------------------

load_swissmetro_database <- function() {
  d <- read.table("test_data_files/swissmetro/swissmetro.dat", header = TRUE)
  # Standard Swissmetro sample (Bierlaire's canonical Biogeme treatment): keep
  # commute (PURPOSE 1) and business (PURPOSE 3) trips with a recorded choice.
  # Availability is NOT filtered here: it is passed to the estimator via av_*
  # columns and masked out of the choice set, so respondents without a car
  # (CAR_AV = 0) are retained rather than dropped.
  d <- d[d$CHOICE != 0 & d$PURPOSE %in% c(1, 3), ]
  d$resp_id <- as.integer(factor(d$ID))
  d <- d[order(d$resp_id), ]
  d$TASK <- ave(d$resp_id, d$resp_id, FUN = seq_along)
  # Balanced panel: keep respondents with the full (modal) task count.
  tpid <- table(d$resp_id)
  keep <- as.integer(names(tpid)[tpid == max(tpid)])
  d <- d[d$resp_id %in% keep, ]
  d$resp_id <- as.integer(factor(d$resp_id))
  d$TASK <- ave(d$resp_id, d$resp_id, FUN = seq_along)

  N <- max(d$resp_id); T_tasks <- max(d$TASK)
  cat(sprintf("Loaded N=%d, T=%d, total rows=%d\n", N, T_tasks, nrow(d)))

  db <- data.frame(ID = d$resp_id, TASK = d$TASK)

  # Alt 1 = train, Alt 2 = SM, Alt 3 = car
  db$x1_1 <- d$TRAIN_TT / 100
  db$x1_2 <- d$SM_TT    / 100
  db$x1_3 <- d$CAR_TT   / 100
  db$x2_1 <- 0  # train has no headway in this coding
  db$x2_2 <- d$SM_HE    / 100
  db$x2_3 <- 0  # car has no headway
  # Standard Swissmetro treatment (Bierlaire's canonical Biogeme spec): season-
  # ticket (GA) holders face zero marginal rail/Swissmetro cost.
  db$price_1 <- d$TRAIN_CO * (d$GA == 0) / 100
  db$price_2 <- d$SM_CO    * (d$GA == 0) / 100
  db$price_3 <- d$CAR_CO   / 100
  db$CHOICE <- d$CHOICE
  # Availability (1 = available): masked out of each task's choice set.
  db$av_1 <- as.integer(d$TRAIN_AV != 0)
  db$av_2 <- as.integer(d$SM_AV != 0)
  db$av_3 <- as.integer(d$CAR_AV != 0)
  db
}

# ---- 2. Pipeline -----------------------------------------------------------

run_sm <- function(database, C_cands = 1:6, dgp = SM_DGP, verbose = TRUE) {
  results <- list()
  for (cc in C_cands) {
    if (verbose) cat(sprintf("\n=== C = %d ===\n", cc))
    t0 <- Sys.time()
    m <- klue_lcmnl(database, cc, dgp = dgp)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose) {
      cat(sprintf("  conv=%s LL=%.2f BIC=%.2f AIC=%.2f ICL=%.2f ICL-BIC=%.2f method=%s time=%.1fs\n",
                  m$converged, m$LL, m$BIC, m$AIC, m$ICL,
                  ifelse(is.na(m$ICL_BIC), 0, m$ICL_BIC),
                  ifelse(is.null(m$best_method) || is.na(m$best_method),
                         "MNL", m$best_method), dt))
    }
    results[[as.character(cc)]] <- m
  }
  results
}

summarise_sm <- function(results) {
  df <- do.call(rbind, lapply(names(results), function(nm) {
    m <- results[[nm]]
    data.frame(C = as.integer(nm), converged = m$converged,
               LL = m$LL, k = m$k, BIC = m$BIC, AIC = m$AIC,
               ICL = m$ICL,
               ICL_BIC = ifelse(is.null(m$ICL_BIC), NA_real_, m$ICL_BIC),
               best_method = ifelse(is.null(m$best_method) ||
                                    is.na(m$best_method),
                                    "MNL", m$best_method),
               stringsAsFactors = FALSE)
  }))
  df$dBIC <- df$BIC - min(df$BIC, na.rm = TRUE)
  df$dAIC <- df$AIC - min(df$AIC, na.rm = TRUE)
  df$dICL <- df$ICL - min(df$ICL, na.rm = TRUE)
  df
}

# ---- 3. Auto-run -----------------------------------------------------------

if (!isTRUE(getOption("sm.skip", FALSE))) {
  database <- load_swissmetro_database()
  results <- run_sm(database)
  summary_df <- summarise_sm(results)

  cat("\n==== SWISSMETRO LCMNL ====\n"); print(summary_df, row.names = FALSE)
  best_C <- summary_df$C[which.min(summary_df$BIC)]
  cat(sprintf("\nBIC-best K* = %d\n", best_C))
  write.csv(summary_df, file.path(OUTPUT_DIR, "swissmetro_lcmnl_results.csv"),
            row.names = FALSE)

  best <- results[[as.character(best_C)]]

  cat("\n=== MMNL ===\n")
  t0 <- Sys.time()
  mmnl <- klue_mmnl(database, dgp = SM_DGP)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  conv=%s LL=%.2f BIC=%.2f AIC=%.2f k=%d time=%.1fs\n",
              mmnl$converged, mmnl$LL, mmnl$BIC, mmnl$AIC, mmnl$k, dt))

  cmp <- data.frame(
    model = c("MNL (C=1)",
              sprintf("LCMNL (C=%d)", best_C),
              "MMNL (indep. normals)"),
    LL  = c(results[["1"]]$LL, best$LL, mmnl$LL),
    k   = c(results[["1"]]$k, best$k, mmnl$k),
    BIC = c(results[["1"]]$BIC, best$BIC, mmnl$BIC),
    AIC = c(results[["1"]]$AIC, best$AIC, mmnl$AIC),
    stringsAsFactors = FALSE
  )
  cmp$dBIC <- cmp$BIC - min(cmp$BIC)
  cat("\n==== SWISSMETRO COMPARISON ====\n"); print(cmp, row.names = FALSE)
  cat(sprintf("BIC-preferred: %s\n", cmp$model[which.min(cmp$BIC)]))
  write.csv(cmp, file.path(OUTPUT_DIR, "swissmetro_model_comparison.csv"),
            row.names = FALSE)
}
