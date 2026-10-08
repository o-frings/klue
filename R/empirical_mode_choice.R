# =============================================================================
# SECOND EMPIRICAL APPLICATION: apollo_modeChoiceData (Hess & Palma 2019)
#
# Applies the hybrid workflow to Apollo's bundled mode-choice dataset: all SP
# tasks, with unavailable alternatives masked out of each task's choice set
# (av_* columns) rather than the tasks dropped. Balanced panel of 500
# respondents x 14 tasks.
#
# Specification (3 generic + cost + 3 ASCs per class):
#   asc_alt1  car intercept            (rail = reference)
#   asc_alt2  bus intercept
#   asc_alt3  air intercept
#   x1        time   (minutes / 60)
#   x2        access (minutes / 60; 0 for car)
#   x3        service (level 0-3; 0 for car and bus)
#   price     cost   (GBP / 10)
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/empirical_mode_choice.R
# With options(mode.skip = TRUE), sourcing it only defines the loader and the
# helpers (dev/rerun_mmnl_benchmark.R does so).
# Writes: output/mode_lcmnl_results.csv (C = 1-6),
#         output/mode_class_betas.csv (the BIC-best fit) and
#         output/mode_model_comparison.csv (MNL, BIC-best LCMNL and the
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

MODE_DGP <- klue_dgp(n_generic = 3, n_alternatives = 4)

# ---- 1. Build wide-format database ----------------------------------------

load_mode_choice_database <- function() {
  data(apollo_modeChoiceData, package = "apollo")
  d <- apollo_modeChoiceData

  # Keep all SP tasks; availability (av_*) is masked out of the choice set by
  # the estimator rather than dropping tasks where an alternative is missing.
  d <- d[d$SP == 1, ]

  d$resp_id <- as.integer(factor(d$ID))
  N <- max(d$resp_id)
  d <- d[order(d$resp_id, d$SP_task), ]

  tasks_per_id <- table(d$resp_id)
  full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == max(tasks_per_id)])
  d <- d[d$resp_id %in% full_ids, ]
  d$resp_id <- as.integer(factor(d$resp_id))
  N <- max(d$resp_id)
  T_tasks <- max(table(d$resp_id))

  cat(sprintf("Balanced subset: N=%d, T=%d, total rows=%d\n",
              N, T_tasks, nrow(d)))

  db <- data.frame(ID = d$resp_id, TASK = d$SP_task)

  scale_time <- 60   # minutes -> hours
  scale_cost <- 10   # GBP / 10
  scale_acc  <- 60   # minutes -> hours
  alts <- c("car", "bus", "air", "rail")

  for (j in seq_along(alts)) {
    a <- alts[j]
    db[[sprintf("x1_%d", j)]] <- d[[paste0("time_", a)]] / scale_time
    db[[sprintf("x2_%d", j)]] <- if (a == "car") 0 else
                                  d[[paste0("access_", a)]] / scale_acc
    db[[sprintf("x3_%d", j)]] <- if (a %in% c("car", "bus")) 0 else
                                  d[[paste0("service_", a)]]
    db[[sprintf("price_%d", j)]] <- d[[paste0("cost_", a)]] / scale_cost
  }

  db$CHOICE <- d$choice
  # Availability (1 = available): masked out of each task's choice set.
  db$av_1 <- as.integer(d$av_car  != 0)
  db$av_2 <- as.integer(d$av_bus  != 0)
  db$av_3 <- as.integer(d$av_air  != 0)
  db$av_4 <- as.integer(d$av_rail != 0)
  db
}

# ---- 2. Multistart LCMNL --------------------------------------------------

run_mode_lcmnl <- function(database, C_cands = 1:6, dgp = MODE_DGP,
                           verbose = TRUE) {
  results <- list()
  for (cc in C_cands) {
    if (verbose) cat(sprintf("\n=== Estimating C = %d ===\n", cc))
    t0 <- Sys.time()
    m <- klue_lcmnl(database, cc, dgp = dgp)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose) {
      cat(sprintf("  converged=%s LL=%.2f BIC=%.2f AIC=%.2f ICL=%.2f ICL-BIC=%.2f method=%s time=%.1fs\n",
                  m$converged, m$LL, m$BIC, m$AIC, m$ICL,
                  ifelse(is.na(m$ICL_BIC), 0, m$ICL_BIC),
                  ifelse(is.null(m$best_method) || is.na(m$best_method),
                         "MNL", m$best_method), dt))
    }
    results[[as.character(cc)]] <- m
  }
  results
}

# ---- 3. MMNL benchmark ----------------------------------------------------

run_mode_mmnl <- function(database, dgp = MODE_DGP, verbose = TRUE) {
  if (verbose) cat("\n=== Estimating MMNL (independent normals) ===\n")
  t0 <- Sys.time()
  m <- klue_mmnl(database, dgp = dgp)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) {
    cat(sprintf("  converged=%s LL=%.2f BIC=%.2f AIC=%.2f k=%d time=%.1fs\n",
                m$converged, m$LL, m$BIC, m$AIC, m$k, dt))
  }
  m
}

# ---- 4. Tidy summary ------------------------------------------------------

summarise_mode <- function(results) {
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

# ---- 5. Auto-run when sourced ---------------------------------------------

if (!isTRUE(getOption("mode.skip", FALSE))) {
  database <- load_mode_choice_database()

  results <- run_mode_lcmnl(database)
  summary_df <- summarise_mode(results)

  cat("\n==== MODE-CHOICE LCMNL RESULTS ====\n")
  print(summary_df, row.names = FALSE)

  best_C <- summary_df$C[which.min(summary_df$BIC)]
  cat(sprintf("\nBIC-best K* = %d\n", best_C))

  write.csv(summary_df, file.path(OUTPUT_DIR, "mode_lcmnl_results.csv"),
            row.names = FALSE)

  best <- results[[as.character(best_C)]]
  attr_labels <- c("time_per60", "access_per60", "service", "cost_per10")
  betas_df <- as.data.frame(best$betas)
  names(betas_df) <- attr_labels
  betas_df$class <- seq_len(nrow(betas_df))
  betas_df$share <- best$class_probs
  betas_df <- betas_df[, c("class", "share", attr_labels)]
  write.csv(betas_df, file.path(OUTPUT_DIR, "mode_class_betas.csv"),
            row.names = FALSE)
  cat("\nClass-specific coefficients (BIC-best model):\n")
  print(betas_df, row.names = FALSE)

  mmnl <- run_mode_mmnl(database)
  cmp <- data.frame(
    model = c("MNL (C=1)",
              sprintf("LCMNL (C=%d, BIC-best)", best_C),
              "MMNL (independent normals)"),
    LL  = c(results[["1"]]$LL, best$LL, mmnl$LL),
    k   = c(results[["1"]]$k,  best$k,  mmnl$k),
    BIC = c(results[["1"]]$BIC, best$BIC, mmnl$BIC),
    AIC = c(results[["1"]]$AIC, best$AIC, mmnl$AIC),
    stringsAsFactors = FALSE
  )
  cmp$dBIC <- cmp$BIC - min(cmp$BIC)

  cat("\n==== MODE-CHOICE MODEL COMPARISON ====\n")
  print(cmp, row.names = FALSE)
  cat(sprintf("\nBIC-preferred: %s\n", cmp$model[which.min(cmp$BIC)]))

  write.csv(cmp, file.path(OUTPUT_DIR, "mode_model_comparison.csv"),
            row.names = FALSE)
}
