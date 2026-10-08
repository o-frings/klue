# =============================================================================
# FOURTH EMPIRICAL APPLICATION: mlogit::Electricity (Train 1998)
#
# Classic LCMNL benchmark: 361 respondents choosing electricity suppliers.
# Wide-format DCE with 4 alternatives x 6 attributes. Filtered to balanced
# panel (respondents with 12 tasks), giving the largest fully-balanced subset:
# 348 respondents x 12 tasks.
#
# Specification (5 generic + price + 3 ASCs per class):
#   asc_alt1, asc_alt2, asc_alt3 (alt 4 = reference)
#   x1 = cl   (contract length, years)
#   x2 = loc  (1 if local utility)
#   x3 = wk   (1 if well-known supplier)
#   x4 = tod  (1 if time-of-day pricing)
#   x5 = seas (1 if seasonal pricing)
#   price = pf  (price in cents/kWh, scaled by 10)
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/empirical_electricity.R
# With options(elec.skip = TRUE), sourcing it only defines the loader and the
# helpers (dev/rerun_mmnl_benchmark.R does so).
# Writes: output/electricity_lcmnl_results.csv (C = 1-6) and
#         output/electricity_model_comparison.csv (MNL, BIC-best LCMNL and the
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

ELEC_DGP <- klue_dgp(n_generic = 5, n_alternatives = 4)

# ---- 1. Build database -----------------------------------------------------

load_electricity_database <- function() {
  data(Electricity, package = "mlogit")
  d <- Electricity

  tasks_per_id <- table(d$id)
  keep_ids <- as.integer(names(tasks_per_id)[tasks_per_id == max(tasks_per_id)])
  d <- d[d$id %in% keep_ids, ]
  d$resp_id <- as.integer(factor(d$id))
  N <- max(d$resp_id)
  T_tasks <- max(table(d$resp_id))
  cat(sprintf("Balanced subset: N=%d, T=%d, total rows=%d\n",
              N, T_tasks, nrow(d)))

  d$TASK <- ave(d$resp_id, d$resp_id, FUN = seq_along)
  db <- data.frame(ID = d$resp_id, TASK = d$TASK)

  for (j in 1:4) {
    db[[sprintf("x1_%d", j)]] <- d[[paste0("cl",  j)]]
    db[[sprintf("x2_%d", j)]] <- d[[paste0("loc", j)]]
    db[[sprintf("x3_%d", j)]] <- d[[paste0("wk",  j)]]
    db[[sprintf("x4_%d", j)]] <- d[[paste0("tod", j)]]
    db[[sprintf("x5_%d", j)]] <- d[[paste0("seas", j)]]
    db[[sprintf("price_%d", j)]] <- d[[paste0("pf",  j)]] / 10
  }
  db$CHOICE <- d$choice
  db
}

# ---- 2. Pipeline -----------------------------------------------------------

run_elec <- function(database, C_cands = 1:6, dgp = ELEC_DGP, verbose = TRUE) {
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

summarise_elec <- function(results) {
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

if (!isTRUE(getOption("elec.skip", FALSE))) {
  database <- load_electricity_database()
  results <- run_elec(database)
  summary_df <- summarise_elec(results)

  cat("\n==== ELECTRICITY LCMNL ====\n"); print(summary_df, row.names = FALSE)
  best_C <- summary_df$C[which.min(summary_df$BIC)]
  cat(sprintf("\nBIC-best K* = %d\n", best_C))
  write.csv(summary_df, file.path(OUTPUT_DIR, "electricity_lcmnl_results.csv"),
            row.names = FALSE)

  best <- results[[as.character(best_C)]]

  cat("\n=== MMNL ===\n")
  t0 <- Sys.time()
  mmnl <- klue_mmnl(database, dgp = ELEC_DGP)
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
  cat("\n==== ELECTRICITY COMPARISON ====\n"); print(cmp, row.names = FALSE)
  cat(sprintf("BIC-preferred: %s\n", cmp$model[which.min(cmp$BIC)]))
  write.csv(cmp, file.path(OUTPUT_DIR, "electricity_model_comparison.csv"),
            row.names = FALSE)
}
