# =============================================================================
# THIRD EMPIRICAL APPLICATION: apollo_swissRouteChoiceData
#
# Swiss commuting route choice (Hess & Palma 2019). Balanced panel of
# 388 respondents x 9 tasks x 2 alternatives.
#
# Specification (3 generic + cost + 1 ASC per class):
#   asc_alt1  intercept on route 1 (route 2 = reference)
#   x1        travel time   (minutes / 60)
#   x2        headway       (minutes / 60)
#   x3        changes
#   price     travel cost   (CHF / 10)
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/empirical_swiss_route.R
# With options(swiss.skip = TRUE), sourcing it only defines the loader and the
# helpers (dev/rerun_mmnl_benchmark.R does so).
# Writes: output/swiss_lcmnl_results.csv (C = 1-6) and
#         output/swiss_model_comparison.csv (MNL, BIC-best LCMNL and the
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

SWISS_DGP <- klue_dgp(n_generic = 3, n_alternatives = 2)

# ---- 1. Build database -----------------------------------------------------

load_swiss_database <- function() {
  data(apollo_swissRouteChoiceData, package = "apollo")
  d <- apollo_swissRouteChoiceData

  d$resp_id <- as.integer(factor(d$ID))
  N <- max(d$resp_id)
  T_tasks <- max(table(d$resp_id))
  cat(sprintf("Loaded N=%d, T=%d, total rows=%d\n",
              N, T_tasks, nrow(d)))

  d$TASK <- ave(d$resp_id, d$resp_id, FUN = seq_along)
  db <- data.frame(ID = d$resp_id, TASK = d$TASK)

  # Scale: time in minutes, cost in CHF
  for (j in 1:2) {
    db[[sprintf("x1_%d", j)]] <- d[[paste0("tt", j)]] / 60
    db[[sprintf("x2_%d", j)]] <- d[[paste0("hw", j)]] / 60
    db[[sprintf("x3_%d", j)]] <- d[[paste0("ch", j)]]
    db[[sprintf("price_%d", j)]] <- d[[paste0("tc", j)]] / 10
  }
  db$CHOICE <- d$choice
  db
}

# ---- 2. Pipeline -----------------------------------------------------------

run_swiss <- function(database, C_cands = 1:6, dgp = SWISS_DGP, verbose = TRUE) {
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

summarise_swiss <- function(results) {
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

if (!isTRUE(getOption("swiss.skip", FALSE))) {
  database <- load_swiss_database()
  results <- run_swiss(database)
  summary_df <- summarise_swiss(results)

  cat("\n==== SWISS ROUTE LCMNL ====\n")
  print(summary_df, row.names = FALSE)
  best_C <- summary_df$C[which.min(summary_df$BIC)]
  cat(sprintf("\nBIC-best K* = %d\n", best_C))
  write.csv(summary_df, file.path(OUTPUT_DIR, "swiss_lcmnl_results.csv"),
            row.names = FALSE)

  best <- results[[as.character(best_C)]]

  cat("\n=== MMNL ===\n")
  t0 <- Sys.time()
  mmnl <- klue_mmnl(database, dgp = SWISS_DGP)
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
  cat("\n==== SWISS COMPARISON ====\n"); print(cmp, row.names = FALSE)
  cat(sprintf("BIC-preferred: %s\n", cmp$model[which.min(cmp$BIC)]))
  write.csv(cmp, file.path(OUTPUT_DIR, "swiss_model_comparison.csv"),
            row.names = FALSE)
}
