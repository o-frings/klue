# =============================================================================
# EMPIRICAL APPLICATION: Vittel water-quality DCE (930 respondents)
#
# Applies the hybrid workflow of the simulation study (klue_lcmnl over
# C = 1-6, then the klue_mmnl benchmark) to a real discrete choice experiment
# without modifying the estimation engine.
#
# Data: 930 respondents, 8 choice sets, 3 alternatives
#   alt 1 = status quo (ASCsq = 1); alts 2 and 3 = policy alternatives
#
# Specification (7 generic dummies + price + 2 ASCs per class):
#   asc_alt1 : status-quo intercept       (alt 3 reference)
#   asc_alt2 : alt-2 intercept            (alt 3 reference)
#   x1  Ban_pesticides
#   x2  RenaturationYes
#   x3  Local_Negative_Effect_Ferti       (fertilizer ref = No_Negative)
#   x4  Mixedfertilizers                  (fertilizer ref = No_Negative)
#   x5  Hedgesconservation                (hedges     ref = Hedgesremoval)
#   x6  Hedgesplantation                  (hedges     ref = Hedgesremoval)
#   x7  Forest_Management_For_Water       (forest     ref = For_Biodiv)
#   price = Waterbill (per month, scaled to per 10 for numerical stability)
#
# IDENTIFICATION (found 2026-09-21). The forest attribute has exactly two levels
# and every policy alternative carries one of them, so the two indicators sum to
# 1 on alternatives 1-2 and to 0 on the status quo, which carries asc_alt1:
# Forest_For_Water + Forest_For_Biodiv + ASC_sq is constant across alternatives.
# Entering both levels leaves the model rank-deficient by one parameter per
# class. The log-likelihood and the class assignment are unaffected (the
# direction is flat), but those three coefficients are then identified only up
# to a constant and k counts a parameter the design cannot identify, which
# over-penalises BIC/AIC/ICL by C*log(N). The forest attribute is now coded
# against a reference level, like the fertilizer and hedges attributes;
# forest_ref = "none" restores the pre-2026-09-21 coding, which attains an
# identical log-likelihood at every C.
#
# Run from the repository root, one core (the data file goes in
# test_data_files/data.csv; README_REPRODUCE.md, "Data"):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript R/empirical_application.R
# With options(empirical.skip = TRUE), sourcing it only defines the loader and
# the helpers (dev/rerun_vittel_respec.R and dev/rerun_mmnl_benchmark.R do so).
# Writes: output/empirical_results.csv (C = 1-6),
#         output/empirical_class_betas.csv (the BIC-best fit) and
#         output/empirical_model_comparison.csv (MNL, BIC-best LCMNL and the
#         independent MMNL; the MMNL row is superseded by
#         dev/rerun_mmnl_benchmark.R and was not rerun).
# The copies in output/ were written on 2026-07-15 (committed 2026-07-24)
# under the pre-2026-09-21 coding; dev/rerun_vittel_respec.R writes the
# identified-coding results to output/vittel_respec_*.
# Seeds: none set here. klue_lcmnl() runs its k-means, Gaussian-mixture and
# PAM clustering starts under seed 123 and restores the caller's random
# stream; klue_mmnl() leaves Apollo's apollo_control$seed at its default, 13,
# under which Apollo makes the 3000 MLHS draws.
# =============================================================================

options(LCMNL_NO_AUTORUN = TRUE)  # load functions only; skip the simulation auto-run
# Canonical engine is the klue package. From the repo, load the dev source
# (latest code); otherwise use the installed package.
if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else {
  suppressWarnings(suppressMessages(library(klue)))
}
if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- "output"

EMP_DGP           <- klue_dgp(n_generic = 7, n_alternatives = 3)  # identified coding
EMP_DGP_COLLINEAR <- klue_dgp(n_generic = 8, n_alternatives = 3)  # forest_ref = "none"

# ---- 1. Build wide-format database from long-format CSV --------------------

# forest_ref: reference level for the two-level forest attribute. "biodiv"
# (default) drops the Biodiversity indicator, so x7 contrasts Water against
# Biodiversity and the model is identified; "none" enters both, the rank-
# deficient pre-2026-09-21 coding (see the identification note in the header).
load_empirical_database <- function(path = "test_data_files/data.csv",
                                    forest_ref = c("biodiv", "none")) {
  forest_ref <- match.arg(forest_ref)
  raw <- read.csv(path, stringsAsFactors = FALSE)

  policy_cols <- c("Ban_pesticides", "RenaturationYes",
                   "Local_Negative_Effect_Ferti", "Mixedfertilizers",
                   "Hedgesconservation", "Hedgesplantation",
                   "Forest_Management_For_Water", "Forest_Managment_For_Biodiv")
  if (forest_ref == "biodiv")
    policy_cols <- setdiff(policy_cols, "Forest_Managment_For_Biodiv")

  price_col <- "Waterbill"

  for (col in c(policy_cols, price_col)) {
    raw[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
    raw[[col]][is.na(raw[[col]])] <- 0
  }
  raw$Waterbill <- raw$Waterbill / 10  # numerical stability

  raw$resp_id <- as.integer(factor(raw$ID))
  raw <- raw[order(raw$resp_id, raw$CS, raw$alt), ]

  N <- max(raw$resp_id); T_tasks <- max(raw$CS); J <- max(raw$alt)
  stopifnot(nrow(raw) == N * T_tasks * J)

  db <- data.frame(ID = rep(seq_len(N), each = T_tasks),
                   TASK = rep(seq_len(T_tasks), times = N))

  long_idx <- function(j) which(raw$alt == j)
  for (j in 1:J) {
    rows <- long_idx(j)
    for (a in seq_along(policy_cols)) {
      db[[sprintf("x%d_%d", a, j)]] <- raw[rows, policy_cols[a]]
    }
    db[[sprintf("price_%d", j)]] <- raw[rows, price_col]
  }

  chosen <- raw[raw$choice == 1, ]
  stopifnot(nrow(chosen) == N * T_tasks)
  chosen <- chosen[order(chosen$resp_id, chosen$CS), ]
  db$CHOICE <- chosen$alt

  attr(db, "attr_labels") <- c(policy_cols, "Waterbill_per10")
  attr(db, "n_generic")   <- length(policy_cols)
  attr(db, "forest_ref")  <- forest_ref
  db
}

# ---- 2. Run multistart LCMNL for C = 1..6 ----------------------------------

run_empirical_lcmnl <- function(database, C_cands = 1:6, dgp = EMP_DGP,
                                verbose = TRUE) {
  results <- list()
  for (cc in C_cands) {
    if (verbose) cat(sprintf("\n=== Estimating C = %d ===\n", cc))
    t0 <- Sys.time()
    m <- klue_lcmnl(database, cc, dgp = dgp)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose) {
      cat(sprintf("  converged=%s  LL=%.2f  BIC=%.2f  AIC=%.2f  ICL=%.2f  ICL-BIC=%.2f  method=%s  time=%.1fs\n",
                  m$converged, m$LL, m$BIC, m$AIC, m$ICL,
                  ifelse(is.na(m$ICL_BIC), 0, m$ICL_BIC),
                  ifelse(is.null(m$best_method) || is.na(m$best_method),
                         "MNL", m$best_method), dt))
    }
    results[[as.character(cc)]] <- m
  }
  results
}

# ---- 2b. MMNL comparison ---------------------------------------------------
# Standard MMNL with independent normal random coefficients (lognormal price).
# klue_mmnl() at its defaults (one Apollo estimation, 3000 MLHS draws).

run_empirical_mmnl <- function(database, dgp = EMP_DGP, verbose = TRUE) {
  if (verbose) cat("\n=== Estimating MMNL (independent normals) ===\n")
  t0 <- Sys.time()
  m <- klue_mmnl(database, dgp = dgp)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (verbose) {
    cat(sprintf("  converged=%s  LL=%.2f  BIC=%.2f  AIC=%.2f  k=%d  time=%.1fs\n",
                m$converged, m$LL, m$BIC, m$AIC, m$k, dt))
  }
  m
}

# ---- 3. Summarise into a tidy data.frame ----------------------------------

summarise_empirical <- function(results) {
  df <- do.call(rbind, lapply(names(results), function(nm) {
    m <- results[[nm]]
    data.frame(C = as.integer(nm),
               converged = m$converged,
               LL  = m$LL,
               k   = m$k,
               BIC = m$BIC,
               AIC = m$AIC,
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

# ---- 4. Run when sourced ---------------------------------------------------

if (isTRUE(getOption("empirical.run", TRUE)) &&
    !isTRUE(getOption("empirical.skip", FALSE))) {
  database <- load_empirical_database()
  cat(sprintf("\nLoaded %d respondents, %d choice sets each, %d alternatives.\n",
              length(unique(database$ID)),
              max(database$TASK),
              EMP_DGP$n_alternatives))

  results <- run_empirical_lcmnl(database)
  summary_df <- summarise_empirical(results)

  cat("\n==== EMPIRICAL RESULTS ====\n")
  print(summary_df, row.names = FALSE)

  best_C <- summary_df$C[which.min(summary_df$BIC)]
  cat(sprintf("\nBIC-best K* = %d\n", best_C))

  write.csv(summary_df, file.path(OUTPUT_DIR, "empirical_results.csv"),
            row.names = FALSE)

  best <- results[[as.character(best_C)]]
  # labels follow the coding: the identified coding (forest_ref = "biodiv")
  # has no Biodiversity indicator
  attr_labels <- c("Ban_pesticides", "RenaturationYes",
                   "Local_Neg_Effect_Ferti", "Mixedfertilizers",
                   "Hedgesconservation", "Hedgesplantation",
                   "Forest_For_Water",
                   if (identical(attr(database, "forest_ref"), "none")) "Forest_For_Biodiv",
                   "Waterbill_per10")
  stopifnot(length(attr_labels) == ncol(best$betas))
  betas_df <- as.data.frame(best$betas)
  names(betas_df) <- attr_labels
  betas_df$class <- seq_len(nrow(betas_df))
  betas_df$share <- best$class_probs
  betas_df <- betas_df[, c("class", "share", attr_labels)]
  write.csv(betas_df, file.path(OUTPUT_DIR, "empirical_class_betas.csv"),
            row.names = FALSE)

  cat("\nClass-specific coefficients (BIC-best model):\n")
  print(betas_df, row.names = FALSE)

  # ---- MMNL comparison (tests H4 in empirical data) -----------------------
  mmnl <- run_empirical_mmnl(database)

  cmp <- data.frame(
    model = c("MNL (C=1)",
              sprintf("LCMNL (C=%d, BIC-best)", best_C),
              "MMNL (independent normals)"),
    LL  = c(results[["1"]]$LL,    best$LL,    mmnl$LL),
    k   = c(results[["1"]]$k,     best$k,     mmnl$k),
    BIC = c(results[["1"]]$BIC,   best$BIC,   mmnl$BIC),
    AIC = c(results[["1"]]$AIC,   best$AIC,   mmnl$AIC),
    stringsAsFactors = FALSE
  )
  cmp$dBIC <- cmp$BIC - min(cmp$BIC)

  cat("\n==== MODEL COMPARISON (H4) ====\n")
  print(cmp, row.names = FALSE)
  cat(sprintf("\nBIC-preferred model: %s\n",
              cmp$model[which.min(cmp$BIC)]))

  write.csv(cmp, file.path(OUTPUT_DIR, "empirical_model_comparison.csv"),
            row.names = FALSE)
}
