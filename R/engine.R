# =============================================================================
# LCMNL SIMULATION STUDY - R/APOLLO IMPLEMENTATION
# =============================================================================
#
# Translates the Julia simulation (full-sumulation-study.jl) to R using the
# Apollo choice modelling package for all model estimation.
#
# Key differences from Julia version:
#   - Class-specific ASCs for alt1/alt2 (alt3 = reference, no ASC)
#   - Apollo-based estimation
#   - R random number generator (different sequences from Julia)
#   - Dynamic segment deviations via cosine rotation (no confound with K)
#   - Multi-start estimation from 6 clustering methods (K-means, GMM,
#     Hierarchical Ward/Complete/Average, PAM)
#
# Author: Oliver Frings
# Date: February 2025
# =============================================================================

library(apollo)
library(mclust)    # GMM clustering
library(cluster)   # PAM, silhouette
library(parallel)  # mclapply for LCMNL parallelisation

options(warn = 1, digits = 4, scipen = 999)

# Number of cores for LCMNL condition-level parallelism (not Apollo).
# Apollo uses its own nCores=10 for within-model parallelism.
# Cap at 4 to avoid excessive fork overhead on macOS.
N_CORES_LCMNL <- min(4L, max(1L, detectCores() - 2L))

# =============================================================================
# CONSTANTS AND DGP DESIGN
# =============================================================================

# DGP configuration: encapsulates problem dimensions so all functions can
# work with variable numbers of attributes and alternatives.
klue_dgp <- function(n_generic = 4, n_alternatives = 3) {
  n_beta <- n_generic + 1L         # generic attributes + price
  n_asc  <- n_alternatives - 1L    # ASCs (last alt = reference)
  npc    <- n_beta + n_asc         # parameters per class
  beta_bar <- c(rep(0.5, n_generic), -1.5)
  attr_names <- c(paste0("x", 1:n_generic), "price")
  list(
    n_generic      = n_generic,
    n_alternatives = n_alternatives,
    n_beta         = n_beta,
    n_asc          = n_asc,
    npc            = npc,
    beta_bar       = beta_bar,
    attr_names     = attr_names,
    price_idx      = n_beta
  )
}

DGP_DEFAULT <- klue_dgp(4, 3)

# Backward-compatible global (used where dgp config not yet threaded through)
BETA_BAR <- DGP_DEFAULT$beta_bar

# Generate class deviation matrix using equally-rotated cosines.
# For K_classes classes, places deviations at equal angular spacing on all
# n_beta attributes (including price). This avoids confounding cluster type
# with number of classes.
# Properties: equal norm sqrt(n_beta/2) for all classes, symmetric spacing.
generate_segment_deviations <- function(K_classes, n_beta = 5L) {
  dev <- matrix(0, nrow = K_classes, ncol = n_beta)
  for (c_idx in 1:K_classes) {
    theta <- 2 * pi * (c_idx - 1) / K_classes
    for (k in 1:n_beta) {
      dev[c_idx, k] <- cos(theta + 2 * pi * (k - 1) / n_beta)
    }
  }
  dev
}

# Configuration
MAX_ITER      <- 500L
OUTPUT_DIR    <- "output"
# Note: package version does not create OUTPUT_DIR at load time. The workflow
# creates whichever output_dir the user requests, lazily.

# =============================================================================
# MMNL DEFAULTS
# -----------------------------------------------------------------------------
# All values below are overridable per-call (klue_mmnl, klue_mmnl_corr,
# klue via `mmnl_opts`) and via global options
# (`getOption("klue.mmnl.<name>")`). Call `klue_mmnl_defaults()` for the
# active settings as a named list.
# -----------------------------------------------------------------------------
N_DRAWS_MMNL              <- 3000L     # interNDraws for the main stage
N_DRAWS_MMNL_STAGE1       <- 200L      # interNDraws for the warm-start stage
DRAWS_TYPE_MMNL           <- "mlhs"    # one of "mlhs","halton","pmc","sobol"
ESTIMATION_ROUTINE_MMNL   <- "bgw"     # one of "bgw","bfgs","nr"
# Box constraints on the log-normal price parameters. The price coefficient
# is b_price = -exp(mu_price + exp(sigma_price) * draws). Without bounds the
# optimiser can drift mu_price into a region where b_price is enormous and
# the likelihood becomes non-finite. Pass NULL to disable.
MU_PRICE_BOUNDS_MMNL      <- c(-5, 3)
SIGMA_PRICE_BOUNDS_MMNL   <- c(-3, 1)

# Default number of cores for Apollo's MMNL cluster. Picks physical cores
# minus one (clamped to >= 1). Override via `n_cores` arg or the
# `klue.mmnl.n_cores` option.
.klue_default_mmnl_cores <- function() {
  nc <- tryCatch(parallel::detectCores(logical = FALSE),
                 error = function(e) NA_integer_)
  if (is.na(nc) || !is.finite(nc) || nc < 1L) nc <- 1L
  max(1L, as.integer(nc) - 1L)
}

#' MMNL default settings
#'
#' Returns the active MMNL defaults as a named list. Values can be overridden
#' per-call (arguments to \code{klue_mmnl} / \code{klue_mmnl_corr} /
#' \code{klue}) or globally via \code{options()} entries
#' \code{klue.mmnl.n_draws}, \code{klue.mmnl.n_draws_stage1},
#' \code{klue.mmnl.draws_type}, \code{klue.mmnl.estimation_routine},
#' \code{klue.mmnl.n_cores}, \code{klue.mmnl.mu_price_bounds},
#' \code{klue.mmnl.sigma_price_bounds}, \code{klue.mmnl.quiet}.
#'
#' @return A named list of the currently active defaults.
#' @export
klue_mmnl_defaults <- function() {
  list(
    n_draws            = getOption("klue.mmnl.n_draws",            N_DRAWS_MMNL),
    n_draws_stage1     = getOption("klue.mmnl.n_draws_stage1",     N_DRAWS_MMNL_STAGE1),
    draws_type         = getOption("klue.mmnl.draws_type",         DRAWS_TYPE_MMNL),
    estimation_routine = getOption("klue.mmnl.estimation_routine", ESTIMATION_ROUTINE_MMNL),
    n_cores            = getOption("klue.mmnl.n_cores",            .klue_default_mmnl_cores()),
    mu_price_bounds    = getOption("klue.mmnl.mu_price_bounds",    MU_PRICE_BOUNDS_MMNL),
    sigma_price_bounds = getOption("klue.mmnl.sigma_price_bounds", SIGMA_PRICE_BOUNDS_MMNL),
    quiet              = getOption("klue.mmnl.quiet",              TRUE)
  )
}

# =============================================================================
# SECTION 1: DATA GENERATION
# =============================================================================

klue_simulate <- function(N_per_class = 150, T_tasks = 20, true_K = 2,
                          separation = 1.0, heterogeneity = 0.25,
                          seed = 12345, class_proportions = NULL,
                          dgp = DGP_DEFAULT, sep_profile = NULL) {
  set.seed(seed)
  J <- dgp$n_alternatives
  n_beta <- dgp$n_beta
  n_generic <- dgp$n_generic
  beta_bar <- dgp$beta_bar

  segment_dev <- generate_segment_deviations(true_K, n_beta)
  true_betas <- matrix(0, nrow = true_K, ncol = n_beta)
  sep_weights <- if (is.null(sep_profile)) rep(1, n_beta) else sep_profile
  for (cc in 1:true_K) {
    true_betas[cc, ] <- beta_bar + separation * sep_weights * segment_dev[cc, ]
  }

  if (is.null(class_proportions)) {
    N <- N_per_class * true_K
    true_class <- rep(1:true_K, each = N_per_class)
  } else {
    N <- N_per_class * true_K
    sizes <- round(N * class_proportions)
    sizes[length(sizes)] <- N - sum(sizes[-length(sizes)])
    true_class <- unlist(lapply(1:true_K, function(cc) rep(cc, sizes[cc])))
    N <- length(true_class)
  }

  sigma_vec <- c(rep(heterogeneity, n_generic), heterogeneity * 0.3)
  individual_betas <- matrix(0, nrow = N, ncol = n_beta)
  for (n in 1:N) {
    individual_betas[n, ] <- rnorm(n_beta, true_betas[true_class[n], ], sigma_vec)
    individual_betas[n, n_beta] <- min(individual_betas[n, n_beta], -0.1)
  }

  n_obs <- N * T_tasks
  database <- data.frame(
    ID   = rep(1:N, each = T_tasks),
    TASK = rep(1:T_tasks, times = N)
  )
  for (j in 1:J) {
    for (a in 1:n_generic) {
      database[[paste0("x", a, "_", j)]] <- runif(n_obs, -1, 1)
    }
    database[[paste0("price_", j)]] <- runif(n_obs, 0.1, 0.9)
  }

  # Vectorised choice simulation via matrix multiply
  beta_rows <- individual_betas[database$ID, ]  # N*T x n_beta
  V_mat <- matrix(0, nrow = n_obs, ncol = J)
  for (j in 1:J) {
    Xj <- matrix(0, nrow = n_obs, ncol = n_beta)
    for (a in 1:n_generic) Xj[, a] <- database[[paste0("x", a, "_", j)]]
    Xj[, n_beta] <- database[[paste0("price_", j)]]
    V_mat[, j] <- rowSums(beta_rows * Xj)
  }
  U_mat <- V_mat - log(-log(matrix(runif(n_obs * J), nrow = n_obs, ncol = J)))
  database$CHOICE <- max.col(U_mat)

  list(database = database, true_betas = true_betas, true_class = true_class,
       individual_betas = individual_betas, N = N, T = T_tasks, K = true_K,
       dgp = dgp)
}

klue_simulate_cov <- function(N_per_class = 150, T_tasks = 20,
                                          true_K = 2, separation = 1.0,
                                          heterogeneity = 0.25, seed = 12345,
                                          covariate_strength = 1.0,
                                          dgp = DGP_DEFAULT) {
  set.seed(seed)
  N <- N_per_class * true_K; J <- dgp$n_alternatives; n_beta <- dgp$n_beta
  n_generic <- dgp$n_generic; beta_bar <- dgp$beta_bar
  Z1 <- rnorm(N); Z2 <- as.numeric(runif(N) > 0.5)

  true_class <- integer(N)
  for (n in 1:N) {
    if (true_K == 2) {
      p <- 1 / (1 + exp(-covariate_strength * (0.5 * Z1[n] + 0.5 * Z2[n])))
      true_class[n] <- ifelse(runif(1) < p, 2, 1)
    } else {
      lp <- covariate_strength * (0.5 * Z1[n] + 0.3 * Z2[n])
      true_class[n] <- max(1, min(true_K, ceiling((true_K + 1) * pnorm(lp))))
    }
  }

  segment_dev <- generate_segment_deviations(true_K, n_beta)
  true_betas <- matrix(0, nrow = true_K, ncol = n_beta)
  for (cc in 1:true_K) true_betas[cc, ] <- beta_bar + separation * segment_dev[cc, ]

  sigma_vec <- c(rep(heterogeneity, n_generic), heterogeneity * 0.3)
  individual_betas <- matrix(0, nrow = N, ncol = n_beta)
  for (n in 1:N) {
    individual_betas[n, ] <- rnorm(n_beta, true_betas[true_class[n], ], sigma_vec)
    individual_betas[n, n_beta] <- min(individual_betas[n, n_beta], -0.1)
  }

  n_obs <- N * T_tasks
  database <- data.frame(ID = rep(1:N, each = T_tasks),
                         TASK = rep(1:T_tasks, times = N))
  for (j in 1:J) {
    for (a in 1:n_generic) database[[paste0("x", a, "_", j)]] <- runif(n_obs, -1, 1)
    database[[paste0("price_", j)]] <- runif(n_obs, 0.1, 0.9)
  }

  # Vectorised choice simulation via matrix multiply
  beta_rows <- individual_betas[database$ID, ]
  V_mat <- matrix(0, nrow = n_obs, ncol = J)
  for (j in 1:J) {
    Xj <- matrix(0, nrow = n_obs, ncol = n_beta)
    for (a in 1:n_generic) Xj[, a] <- database[[paste0("x", a, "_", j)]]
    Xj[, n_beta] <- database[[paste0("price_", j)]]
    V_mat[, j] <- rowSums(beta_rows * Xj)
  }
  U_mat <- V_mat - log(-log(matrix(runif(n_obs * J), nrow = n_obs, ncol = J)))
  database$CHOICE <- max.col(U_mat)
  database$Z1 <- Z1[database$ID]
  database$Z2 <- Z2[database$ID]

  list(database = database, true_betas = true_betas, true_class = true_class,
       individual_betas = individual_betas, N = N, T = T_tasks, K = true_K,
       Z1 = Z1, Z2 = Z2, dgp = dgp)
}

klue_simulate_deff <- function(N_per_class = 150, T_tasks = 20,
                                     true_K = 2, separation = 1.0,
                                     heterogeneity = 0.25, seed = 12345,
                                     dgp = DGP_DEFAULT) {
  set.seed(seed)
  N <- N_per_class * true_K; J <- dgp$n_alternatives; n_beta <- dgp$n_beta
  n_generic <- dgp$n_generic; beta_bar <- dgp$beta_bar

  segment_dev <- generate_segment_deviations(true_K, n_beta)
  true_betas <- matrix(0, nrow = true_K, ncol = n_beta)
  for (cc in 1:true_K) true_betas[cc, ] <- beta_bar + separation * segment_dev[cc, ]
  true_class <- rep(1:true_K, each = N_per_class)

  sigma_vec <- c(rep(heterogeneity, n_generic), heterogeneity * 0.3)
  individual_betas <- matrix(0, nrow = N, ncol = n_beta)
  for (n in 1:N) {
    individual_betas[n, ] <- rnorm(n_beta, true_betas[true_class[n], ], sigma_vec)
    individual_betas[n, n_beta] <- min(individual_betas[n, n_beta], -0.1)
  }

  lvl_c <- c(-0.8, -0.4, 0.0, 0.4, 0.8)
  lvl_p <- c(0.2, 0.35, 0.5, 0.65, 0.8)

  n_obs <- N * T_tasks
  database <- data.frame(ID = rep(1:N, each = T_tasks),
                         TASK = rep(1:T_tasks, times = N))
  for (j in 1:J) {
    for (a in 1:n_generic) database[[paste0("x", a, "_", j)]] <- numeric(n_obs)
    database[[paste0("price_", j)]] <- numeric(n_obs)
  }

  row_idx <- 0
  for (n in 1:N) {
    for (t in 1:T_tasks) {
      row_idx <- row_idx + 1
      for (j in 1:J) {
        for (k in 1:n_generic) {
          idx <- ((n + t + j + k) %% 5) + 1
          database[[paste0("x", k, "_", j)]][row_idx] <- lvl_c[idx] + 0.1 * (2 * runif(1) - 1)
        }
        pidx <- ((n + t + j) %% 5) + 1
        database[[paste0("price_", j)]][row_idx] <- lvl_p[pidx] + 0.05 * (2 * runif(1) - 1)
      }
    }
  }

  # Vectorised choice simulation via matrix multiply
  beta_rows <- individual_betas[database$ID, ]
  V_mat <- matrix(0, nrow = n_obs, ncol = J)
  for (j in 1:J) {
    Xj <- matrix(0, nrow = n_obs, ncol = n_beta)
    for (a in 1:n_generic) Xj[, a] <- database[[paste0("x", a, "_", j)]]
    Xj[, n_beta] <- database[[paste0("price_", j)]]
    V_mat[, j] <- rowSums(beta_rows * Xj)
  }
  U_mat <- V_mat - log(-log(matrix(runif(n_obs * J), nrow = n_obs, ncol = J)))
  database$CHOICE <- max.col(U_mat)

  list(database = database, true_betas = true_betas, true_class = true_class,
       individual_betas = individual_betas, N = N, T = T_tasks, K = true_K,
       dgp = dgp)
}

# =============================================================================
# SECTION 2: CLUSTERING FOR STARTING VALUES
# =============================================================================

# Helper: get column names for attribute a across all alternatives
attr_col_names <- function(attr_name, J) {
  paste0(attr_name, "_", 1:J)
}

compute_rp_features <- function(database, dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  n_obs <- nrow(database)
  J <- dgp$n_alternatives
  n_beta <- dgp$n_beta
  n_generic <- dgp$n_generic
  CH <- database$CHOICE
  ri <- 1:n_obs

  # Build attribute matrices (n_obs x J) per attribute, compute chosen-vs-unchosen diffs
  diffs <- matrix(0, nrow = n_obs, ncol = n_beta)
  for (a in 1:n_beta) {
    aname <- if (a <= n_generic) paste0("x", a) else "price"
    Xa <- matrix(0, nrow = n_obs, ncol = J)
    for (j in 1:J) Xa[, j] <- database[[paste0(aname, "_", j)]]
    rsa <- rowSums(Xa)
    cha <- Xa[cbind(ri, CH)]
    diffs[, a] <- cha - (rsa - cha) / (J - 1)
  }

  # Average per individual (balanced panel: each ID has same number of rows)
  T_per_n <- as.integer(n_obs / N)
  features <- matrix(0, nrow = N, ncol = n_beta)
  for (k in 1:n_beta) {
    features[, k] <- colSums(matrix(diffs[, k], T_per_n, N)) / T_per_n
  }
  features
}

# --- Helper: standardise RP features and compute cluster means ---
standardise_features <- function(features) {
  mu <- colMeans(features); sigma <- apply(features, 2, sd)
  sigma[sigma == 0] <- 1
  scaled <- scale(features, center = mu, scale = sigma)
  scaled[!is.finite(scaled)] <- 0
  list(scaled = scaled, mu = mu, sigma = sigma)
}

# Fit separate MNL (C=1) per cluster to obtain coefficient-level starting values.
# Instead of heuristically scaling cluster centroids, this fits a proper MNL model
# on each cluster's observations, producing starting values directly in the
# coefficient space.
fit_cluster_mnls <- function(labels, database, dgp = DGP_DEFAULT) {
  C <- max(labels)
  all_ids <- unique(database$ID)
  N <- length(all_ids)
  n_beta <- dgp$n_beta

  betas <- matrix(0, nrow = C, ncol = n_beta)
  shares <- as.numeric(table(factor(labels, levels = 1:C))) / N

  for (cc in 1:C) {
    cluster_ids <- all_ids[labels == cc]
    if (length(cluster_ids) < 3) {
      betas[cc, ] <- c(rep(0.5, dgp$n_generic), -0.5)
      next
    }
    db_sub <- database[database$ID %in% cluster_ids, , drop = FALSE]
    db_sub$ID <- match(db_sub$ID, cluster_ids)

    mnl_fit <- tryCatch(
      estimate_lcmnl(db_sub, C = 1,
                     start_betas = matrix(0, 1, n_beta), dgp = dgp),
      error = function(e) NULL
    )
    if (!is.null(mnl_fit) && mnl_fit$converged) {
      betas[cc, ] <- mnl_fit$betas[1, ]
    } else {
      betas[cc, ] <- c(rep(0.5, dgp$n_generic), -0.5)
    }
  }

  list(betas = betas, shares = shares)
}

# --- 1. K-Means ---
get_kmeans_starts <- function(database, C,
                              features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))

  sf <- standardise_features(features)
  set.seed(123)
  km <- kmeans(sf$scaled, centers = C, nstart = 25, iter.max = 100)
  fit_cluster_mnls(km$cluster, database, dgp)
}

# --- 2. GMM (Gaussian Mixture Model via mclust) ---
get_gmm_starts <- function(database, C,
                           features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))

  sf <- standardise_features(features)
  set.seed(123)
  gmm <- mclust::Mclust(sf$scaled, G = C, verbose = FALSE)
  if (is.null(gmm)) return(NULL)

  fit_cluster_mnls(gmm$classification, database, dgp)
}

# --- 3. Hierarchical Clustering (Ward's D2) ---
get_hc_ward_starts <- function(database, C,
                               features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))
  sf <- standardise_features(features)
  hc <- hclust(dist(sf$scaled), method = "ward.D2")
  labels <- cutree(hc, k = C)
  fit_cluster_mnls(labels, database, dgp)
}

# --- 4. Hierarchical Clustering (Complete linkage) ---
get_hc_complete_starts <- function(database, C,
                                   features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))
  sf <- standardise_features(features)
  hc <- hclust(dist(sf$scaled), method = "complete")
  labels <- cutree(hc, k = C)
  fit_cluster_mnls(labels, database, dgp)
}

# --- 5. Hierarchical Clustering (Average linkage) ---
get_hc_average_starts <- function(database, C,
                                  features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))
  sf <- standardise_features(features)
  hc <- hclust(dist(sf$scaled), method = "average")
  labels <- cutree(hc, k = C)
  fit_cluster_mnls(labels, database, dgp)
}

# --- 6. PAM (K-Medoids) ---
get_pam_starts <- function(database, C,
                           features = NULL, dgp = DGP_DEFAULT) {
  if (is.null(features)) features <- compute_rp_features(database, dgp)
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))
  sf <- standardise_features(features)
  set.seed(123)
  pam_res <- cluster::pam(sf$scaled, k = C)
  fit_cluster_mnls(pam_res$clustering, database, dgp)
}

# --- Multi-start wrapper: all 6 methods (computes RP features once) ---
get_all_starts <- function(database, C, dgp = DGP_DEFAULT) {
  features <- compute_rp_features(database, dgp)
  methods <- list(
    kmeans      = get_kmeans_starts,
    gmm         = get_gmm_starts,
    hc_ward     = get_hc_ward_starts,
    hc_complete = get_hc_complete_starts,
    hc_average  = get_hc_average_starts,
    pam         = get_pam_starts
  )
  starts_list <- list()
  for (nm in names(methods)) {
    starts_list[[nm]] <- tryCatch(
      methods[[nm]](database, C, features = features, dgp = dgp),
      error = function(e) NULL
    )
  }
  starts_list
}

# Helper: build design matrices X[[j]] (n_obs x n_beta) for each alternative j.
# Used by estimate_lcmnl, compute_lc_posteriors, and other estimation functions.
build_design_matrices <- function(database, dgp = DGP_DEFAULT) {
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  n_obs <- nrow(database)
  X <- vector("list", J)
  for (j in 1:J) {
    Xj <- matrix(0, nrow = n_obs, ncol = n_beta)
    for (a in 1:n_generic) Xj[, a] <- database[[paste0("x", a, "_", j)]]
    Xj[, n_beta] <- database[[paste0("price_", j)]]
    X[[j]] <- Xj
  }
  X
}

# =============================================================================
# SECTION 3: MNL/LCMNL ESTIMATION (DIRECT MLE) AND APOLLO HELPERS
# Note: MNL is the special case of LCMNL with C = 1 (single class).
# =============================================================================

cleanup_apollo <- function() {
  to_remove <- c(
    grep("^apollo_", ls(envir = .GlobalEnv), value = TRUE),
    grep("^asc_alt", ls(envir = .GlobalEnv), value = TRUE),
    grep("^b_x[0-9]", ls(envir = .GlobalEnv), value = TRUE),
    grep("^b_price", ls(envir = .GlobalEnv), value = TRUE),
    grep("^delta_", ls(envir = .GlobalEnv), value = TRUE),
    grep("^mu_", ls(envir = .GlobalEnv), value = TRUE),
    grep("^sigma_", ls(envir = .GlobalEnv), value = TRUE)
  )
  # Protect our own functions and constants
  protected <- c("cleanup_apollo", "compute_rp_features",
                  "compute_lc_posteriors",
                  "standardise_features",
                  "fit_cluster_mnls", "attr_col_names", "build_design_matrices",
                  "klue_dgp", "DGP_DEFAULT",
                  "get_kmeans_starts", "get_gmm_starts",
                  "get_hc_ward_starts", "get_hc_complete_starts",
                  "get_hc_average_starts", "get_pam_starts", "get_all_starts",
                  "generate_segment_deviations",
                  "make_apollo_lcPars", "make_apollo_probabilities_lc",
                  "klue_simulate", "klue_simulate_cov",
                  "klue_simulate_deff",
                  "estimate_lcmnl", "klue_lcmnl", "klue_mmnl",
                  "compute_onehot_features", "get_all_starts_onehot",
                  "estimate_lcmnl_multistart_onehot", "run_initialisation_ablation",
                  "compute_ari", "compute_recovery",
                  "run_main_simulation", "summarise_main_results",
                  "run_mmnl_comparison", "run_convergence_ablation",
                  "run_sample_sensitivity", "run_correlated_mmnl_robustness",
                  "klue_mmnl_corr", ".run_apollo_mmnl_corr",
                  "run_unbalanced_analysis",
                  "run_design_comparison", "run_concomitant_analysis",
                  "run_unconditional_recovery", "run_clustering_comparison",
                  "run_full_study")
  to_remove <- setdiff(unique(to_remove), protected)
  if (length(to_remove) > 0) rm(list = to_remove, envir = .GlobalEnv)
}

# Helper: build apollo_lcPars with NO for-loops (Apollo checkIndices rejects them)
# Uses apollo_classAlloc for class probabilities — this produces proper gradient
# structures so apollo_lc can compute analytic gradients (much faster than
# numerical gradients from manual softmax + manual mixing).
make_apollo_lcPars <- function(C, dgp = DGP_DEFAULT) {
  n_generic <- dgp$n_generic; J <- dgp$n_alternatives
  # ASC names: asc_alt1 .. asc_alt{J-1}, beta names: b_x1..b_x{n_generic}, b_price
  asc_names <- paste0("asc_alt", 1:(J - 1))
  beta_names <- c(paste0("b_x", 1:n_generic), "b_price")
  param_names <- c(asc_names, beta_names)

  lines <- c("function(apollo_beta, apollo_inputs) {", "  lcpars <- list()")
  for (pname in param_names) {
    items <- paste0(pname, "_", 1:C)
    lines <- c(lines, sprintf('  lcpars[["%s"]] <- list(%s)', pname, paste(items, collapse = ", ")))
  }
  class_entries <- paste(sprintf('class_%d = %d', 1:C, 1:C), collapse = ", ")
  util_entries  <- paste(sprintf('class_%d = delta_%d', 1:C, 1:C), collapse = ", ")
  lines <- c(lines,
    sprintf('  classAlloc_settings <- list(classes = c(%s), utilities = list(%s))',
            class_entries, util_entries),
    '  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)',
    '  return(lcpars)',
    '}'
  )
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

# Helper: build apollo_probabilities for LC model with NO for-loops
# Uses apollo_lc for proper gradient propagation (enables analytic gradients).
make_apollo_probabilities_lc <- function(C, dgp = DGP_DEFAULT) {
  J <- dgp$n_alternatives; n_generic <- dgp$n_generic
  # Build alternatives and availability
  alt_entries <- paste(sprintf('alt%d = %d', 1:J, 1:J), collapse = ", ")
  avail_entries <- paste(sprintf('alt%d = 1', 1:J), collapse = ", ")

  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    sprintf('  mnl_settings <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE)',
            alt_entries, avail_entries)
  )
  for (s in 1:C) {
    lines <- c(lines, '  V <- list()')
    for (j in 1:J) {
      # Build utility: ASC (if not reference alt) + sum of beta * x terms
      terms <- c()
      if (j < J) terms <- c(terms, sprintf('asc_alt%d[[%d]]', j, s))
      for (a in 1:n_generic) terms <- c(terms, sprintf('b_x%d[[%d]] * x%d_%d', a, s, a, j))
      terms <- c(terms, sprintf('b_price[[%d]] * price_%d', s, j))
      lines <- c(lines, sprintf('  V[["alt%d"]] <- %s', j, paste(terms, collapse = " + ")))
    }
    lines <- c(lines,
      '  mnl_settings$utilities <- V',
      sprintf('  P[["class_%d"]] <- apollo_mnl(mnl_settings, functionality)', s),
      sprintf('  P[["class_%d"]] <- apollo_panelProd(P[["class_%d"]], apollo_inputs, functionality)', s, s)
    )
  }
  inclass_entries <- paste(sprintf('class_%d = P[["class_%d"]]', 1:C, 1:C), collapse = ", ")
  lines <- c(lines,
    sprintf('  lcSettings <- list(inClassProb = list(%s), classProb = pi_values)',
            inclass_entries),
    '  P[["model"]] <- apollo_lc(lcSettings, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)',
    '  return(P)',
    '}'
  )
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

# Compute individual-level posterior class probabilities via Bayes' rule.
# Uses log-space arithmetic to avoid underflow from panel products over T tasks.
# Uses design matrices for dimension-agnostic computation.
compute_lc_posteriors <- function(database, C, est, dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  if (C == 1) return(matrix(1, nrow = N, ncol = 1))

  T_total <- nrow(database)
  T_per_n <- as.integer(T_total / N)
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  X <- build_design_matrices(database, dgp)

  # Class probabilities from deltas
  deltas <- numeric(C)
  for (ci in 1:C) {
    pn <- paste0("delta_", ci)
    deltas[ci] <- if (pn %in% names(est)) est[pn] else 0
  }
  exp_d <- exp(deltas - max(deltas))
  log_pi <- log(exp_d / sum(exp_d))

  log_lik <- matrix(0, nrow = N, ncol = C)

  for (ci in 1:C) {
    # Extract betas for this class
    betas_c <- numeric(n_beta)
    for (a in 1:n_generic) betas_c[a] <- est[paste0("b_x", a, "_", ci)]
    betas_c[n_beta] <- est[paste0("b_price_", ci)]

    # ASCs for this class
    ascs <- numeric(J)
    for (j in 1:(J - 1)) ascs[j] <- est[paste0("asc_alt", j, "_", ci)]

    # Compute utilities: V[,j] = X[[j]] %*% betas_c + asc_j
    V <- matrix(0, nrow = T_total, ncol = J)
    for (j in 1:J) V[, j] <- X[[j]] %*% betas_c + ascs[j]

    V_max <- do.call(pmax, lapply(1:J, function(j) V[, j]))
    lse <- V_max + log(rowSums(exp(V - V_max)))
    ch <- database$CHOICE
    choice_V <- V[cbind(1:T_total, ch)]
    log_prob <- choice_V - lse

    log_lik[, ci] <- colSums(matrix(log_prob, nrow = T_per_n, ncol = N))
  }

  log_joint <- sweep(log_lik, 2, log_pi, "+")
  log_max <- apply(log_joint, 1, max)
  log_denom <- log_max + log(rowSums(exp(log_joint - log_max)))
  exp(log_joint - log_denom)
}

# Direct MLE for MNL/LCMNL using optim() with BFGS and analytic gradients.
# When C = 1 this is a standard MNL; when C >= 2 it is a latent class MNL
# (discrete mixture of MNL components).  Bypasses Apollo entirely since
# LCMNL involves only discrete (not continuous) mixing.
# Uses matrix operations for dimension-agnostic computation.
# Parameter layout per class: [beta_1..beta_{n_beta}, asc_1..asc_{J-1}]
estimate_lcmnl <- function(database, C, start_betas = NULL, start_shares = NULL,
                           dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  T_total <- nrow(database)
  T_per_n <- as.integer(T_total / N)
  n_beta <- dgp$n_beta; n_asc <- dgp$n_asc; npc <- dgp$npc; J <- dgp$n_alternatives
  n_generic <- dgp$n_generic

  if (is.null(start_betas)) {
    starts <- get_kmeans_starts(database, C, dgp = dgp)
    start_betas  <- starts$betas
    start_shares <- starts$shares
  }
  if (is.null(start_shares)) start_shares <- rep(1 / C, C)

  fail_result <- list(
    converged = FALSE, C = C, LL = -Inf, BIC = Inf, AIC = Inf, ICL = Inf,
    ICL_BIC = NA_real_,
    k = 0, betas = matrix(0, C, n_beta), class_probs = rep(1 / C, C),
    posteriors = matrix(1 / C, N, C)
  )

  # Pre-build design matrices: X[[j]] is n_obs x n_beta
  X <- build_design_matrices(database, dgp)

  # Pre-compute choice indicators: ch_ind[,j] = 1 if CHOICE==j
  ch <- database$CHOICE
  ch_ind <- matrix(0, nrow = T_total, ncol = J)
  for (j in 1:J) ch_ind[, j] <- as.numeric(ch == j)

  # Pre-compute chosen-alternative design matrix: Xc[t,a] = X_{choice_t}[t,a]
  Xc <- matrix(0, nrow = T_total, ncol = n_beta)
  for (j in 1:J) Xc <- Xc + ch_ind[, j] * X[[j]]

  if (C == 1) {
    par0 <- c(start_betas[1, ], rep(0, n_asc))
    n_free <- npc

    neg_ll <- function(par) {
      betas <- par[1:n_beta]
      ascs <- c(par[(n_beta + 1):npc], 0)  # last alt = reference
      V <- matrix(0, T_total, J)
      for (j in 1:J) V[, j] <- X[[j]] %*% betas + ascs[j]
      Vm <- do.call(pmax, lapply(1:J, function(j) V[, j]))
      choice_V <- rowSums(ch_ind * V)
      -sum(choice_V - Vm - log(rowSums(exp(V - Vm))))
    }

    grad_ll <- function(par) {
      betas <- par[1:n_beta]
      ascs <- c(par[(n_beta + 1):npc], 0)
      V <- matrix(0, T_total, J)
      for (j in 1:J) V[, j] <- X[[j]] %*% betas + ascs[j]
      Vm <- do.call(pmax, lapply(1:J, function(j) V[, j]))
      eV <- exp(V - Vm); probs <- eV / rowSums(eV)
      # Beta gradient: sum_t (Xc - sum_j p_j * X_j)
      EX <- matrix(0, T_total, n_beta)
      for (j in 1:J) EX <- EX + probs[, j] * X[[j]]
      g_beta <- -colSums(Xc - EX)
      # ASC gradient: sum_t (ch_ind_j - p_j) for j = 1..(J-1)
      g_asc <- numeric(n_asc)
      for (j in 1:n_asc) g_asc[j] <- -sum(ch_ind[, j] - probs[, j])
      c(g_beta, g_asc)
    }
  } else {
    par0 <- numeric(C * npc + C - 1L)
    for (ci in 1:C) {
      off <- (ci - 1L) * npc
      par0[off + 1:n_beta] <- start_betas[ci, ]
    }
    for (ci in 1:(C - 1L)) {
      par0[C * npc + ci] <- log(max(start_shares[ci], 0.01) /
                                 max(start_shares[C], 0.01))
    }
    n_free <- C * npc + C - 1L

    neg_ll <- function(par) {
      deltas <- c(par[(C * npc + 1L):(C * npc + C - 1L)], 0)
      dm <- max(deltas)
      log_pi <- deltas - dm - log(sum(exp(deltas - dm)))

      log_panel <- matrix(0, N, C)
      for (ci in 1:C) {
        off <- (ci - 1L) * npc
        betas <- par[off + 1:n_beta]
        ascs <- c(par[off + (n_beta + 1):npc], 0)
        V <- matrix(0, T_total, J)
        for (j in 1:J) V[, j] <- X[[j]] %*% betas + ascs[j]
        Vm <- do.call(pmax, lapply(1:J, function(j) V[, j]))
        choice_V <- rowSums(ch_ind * V)
        log_panel[, ci] <- colSums(matrix(
          choice_V - Vm - log(rowSums(exp(V - Vm))), T_per_n, N))
      }
      log_joint <- sweep(log_panel, 2, log_pi, "+")
      lm <- apply(log_joint, 1, max)
      -sum(lm + log(rowSums(exp(log_joint - lm))))
    }

    grad_ll <- function(par) {
      deltas <- c(par[(C * npc + 1L):(C * npc + C - 1L)], 0)
      dm <- max(deltas)
      ed <- exp(deltas - dm); pi_c <- ed / sum(ed)

      log_panel <- matrix(0, N, C)
      pscores <- vector("list", C)

      for (ci in 1:C) {
        off <- (ci - 1L) * npc
        betas <- par[off + 1:n_beta]
        ascs <- c(par[off + (n_beta + 1):npc], 0)
        V <- matrix(0, T_total, J)
        for (j in 1:J) V[, j] <- X[[j]] %*% betas + ascs[j]
        Vm <- do.call(pmax, lapply(1:J, function(j) V[, j]))
        eV <- exp(V - Vm); probs <- eV / rowSums(eV)
        choice_V <- rowSums(ch_ind * V)

        log_panel[, ci] <- colSums(matrix(
          choice_V - Vm - log(rowSums(eV)), T_per_n, N))

        # MNL scores: chosen_x - E[x] for betas, chosen_ind - prob for ASCs
        EX <- matrix(0, T_total, n_beta)
        for (j in 1:J) EX <- EX + probs[, j] * X[[j]]
        score_beta <- Xc - EX  # T_total x n_beta
        score_asc <- matrix(0, T_total, n_asc)
        for (j in 1:n_asc) score_asc[, j] <- ch_ind[, j] - probs[, j]

        # Panel-sum scores: T_per_n x N -> N x npc
        ps <- matrix(0, N, npc)
        for (k in 1:n_beta) ps[, k] <- colSums(matrix(score_beta[, k], T_per_n, N))
        for (k in 1:n_asc) ps[, n_beta + k] <- colSums(matrix(score_asc[, k], T_per_n, N))
        pscores[[ci]] <- ps
      }

      log_joint <- sweep(log_panel, 2, log(pi_c), "+")
      lm <- apply(log_joint, 1, max)
      w <- exp(log_joint - (lm + log(rowSums(exp(log_joint - lm)))))

      g <- numeric(C * npc + C - 1L)
      for (ci in 1:C) {
        g[(ci - 1L) * npc + 1:npc] <- -colSums(w[, ci] * pscores[[ci]])
      }
      for (ci in 1:(C - 1L)) {
        g[C * npc + ci] <- -(sum(w[, ci]) - N * pi_c[ci])
      }
      g
    }
  }

  result <- tryCatch(
    suppressWarnings(optim(par0, neg_ll, gr = grad_ll, method = "BFGS",
                           control = list(maxit = MAX_ITER, reltol = 1e-10))),
    error = function(e) {
      label <- if (C == 1L) "MNL" else paste0("LCMNL C=", C)
      message("[", label, "] optim error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(result) || result$convergence != 0) return(fail_result)

  p <- result$par
  LL <- -result$value

  if (C == 1) {
    betas_mat   <- matrix(p[1:n_beta], nrow = 1)
    class_probs <- 1
    posteriors  <- matrix(1, nrow = N, ncol = 1)
    BIC <- -2 * LL + n_free * log(N)
    AIC <- -2 * LL + 2 * n_free
    ICL <- BIC
    ICL_BIC <- 0
  } else {
    est <- c()
    betas_mat <- matrix(0, C, n_beta)
    for (ci in 1:C) {
      off <- (ci - 1L) * npc
      betas_mat[ci, ] <- p[off + 1:n_beta]
      for (a in 1:n_generic) est[paste0("b_x", a, "_", ci)] <- p[off + a]
      est[paste0("b_price_", ci)] <- p[off + n_beta]
      for (j in 1:n_asc) est[paste0("asc_alt", j, "_", ci)] <- p[off + n_beta + j]
    }
    for (ci in 1:(C - 1L)) est[paste0("delta_", ci)] <- p[C * npc + ci]
    est[paste0("delta_", C)] <- 0

    deltas <- c(p[(C * npc + 1L):(C * npc + C - 1L)], 0)
    exp_d <- exp(deltas - max(deltas))
    class_probs <- as.numeric(exp_d / sum(exp_d))

    BIC <- -2 * LL + n_free * log(N)
    AIC <- -2 * LL + 2 * n_free

    posteriors <- compute_lc_posteriors(database, C, est, dgp)
    posteriors <- pmax(posteriors, 1e-100)
    H   <- -sum(posteriors * log(posteriors))
    ICL <- BIC + 2 * H
    ICL_BIC <- 2 * H
  }

  list(converged = TRUE, C = C,
       model_type = if (C == 1L) "MNL" else "LCMNL",
       LL = LL, BIC = BIC, AIC = AIC, ICL = ICL, ICL_BIC = ICL_BIC,
       k = n_free, betas = betas_mat, class_probs = class_probs,
       posteriors = posteriors)
}

# --- Multi-start MNL/LCMNL: for C=1 (MNL) a single run suffices;
# --- for C>=2 (LCMNL) run from all 6 clustering methods, keep best LL ---
klue_lcmnl <- function(database, C, dgp = DGP_DEFAULT) {
  all_starts <- get_all_starts(database, C, dgp = dgp)
  N <- length(unique(database$ID))

  best <- list(
    converged = FALSE, C = C, LL = -Inf, BIC = Inf, AIC = Inf, ICL = Inf,
    ICL_BIC = NA_real_,
    k = 0, betas = matrix(0, C, dgp$n_beta), class_probs = rep(1 / C, C),
    posteriors = matrix(1 / C, N, C), best_method = NA_character_
  )
  method_results <- list()

  for (nm in names(all_starts)) {
    if (is.null(all_starts[[nm]])) next
    res <- tryCatch(
      estimate_lcmnl(database, C,
                     start_betas = all_starts[[nm]]$betas,
                     start_shares = all_starts[[nm]]$shares,
                     dgp = dgp),
      error = function(e) NULL
    )
    if (is.null(res)) next
    method_results[[nm]] <- res
    if (res$converged && res$LL > best$LL) {
      best <- res
      best$best_method <- nm
    }
  }
  best$method_results <- method_results
  best
}

# -----------------------------------------------------------------------------
# One-hot initialisation variant (Level 2 ablation: compare RP contrasts to
# one-hot encoded choice indicators as the clustering feature). The structure
# parallels the main RP-contrast pipeline; only the feature matrix differs.
# -----------------------------------------------------------------------------

# N x (T*J) matrix; entry (i, (t-1)*J + j) is 1 if respondent i chose alt j
# at task t, else 0. Assumes a balanced panel.
compute_onehot_features <- function(database, dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  n_obs <- nrow(database)
  J <- dgp$n_alternatives
  T_per_n <- as.integer(n_obs / N)
  stopifnot(N * T_per_n == n_obs)

  ids <- unique(database$ID)
  features <- matrix(0, nrow = N, ncol = T_per_n * J)
  for (i in seq_along(ids)) {
    rows <- which(database$ID == ids[i])
    chosen <- database$CHOICE[rows]
    col_idx <- ((seq_len(T_per_n) - 1L) * J) + chosen
    features[i, col_idx] <- 1
  }
  features
}

get_all_starts_onehot <- function(database, C, dgp = DGP_DEFAULT) {
  features <- compute_onehot_features(database, dgp)
  methods <- list(
    kmeans      = get_kmeans_starts,
    gmm         = get_gmm_starts,
    hc_ward     = get_hc_ward_starts,
    hc_complete = get_hc_complete_starts,
    hc_average  = get_hc_average_starts,
    pam         = get_pam_starts
  )
  starts_list <- list()
  for (nm in names(methods)) {
    starts_list[[nm]] <- tryCatch(
      methods[[nm]](database, C, features = features, dgp = dgp),
      error = function(e) NULL
    )
  }
  starts_list
}

estimate_lcmnl_multistart_onehot <- function(database, C, dgp = DGP_DEFAULT) {
  all_starts <- get_all_starts_onehot(database, C, dgp = dgp)
  N <- length(unique(database$ID))

  best <- list(
    converged = FALSE, C = C, LL = -Inf, BIC = Inf, AIC = Inf, ICL = Inf,
    ICL_BIC = NA_real_,
    k = 0, betas = matrix(0, C, dgp$n_beta), class_probs = rep(1 / C, C),
    posteriors = matrix(1 / C, N, C), best_method = NA_character_
  )
  for (nm in names(all_starts)) {
    if (is.null(all_starts[[nm]])) next
    res <- tryCatch(
      estimate_lcmnl(database, C,
                     start_betas = all_starts[[nm]]$betas,
                     start_shares = all_starts[[nm]]$shares,
                     dgp = dgp),
      error = function(e) NULL
    )
    if (is.null(res)) next
    if (res$converged && res$LL > best$LL) {
      best <- res
      best$best_method <- nm
    }
  }
  best
}

# =============================================================================
# SECTION 4: MMNL ESTIMATION VIA APOLLO
# =============================================================================

# Helper: dynamically build apollo_randCoeff function for MMNL (independent)
.make_apollo_randCoeff <- function(dgp = DGP_DEFAULT) {
  n_generic <- dgp$n_generic
  lines <- c("function(apollo_beta, apollo_inputs) {", "  randcoeff <- list()")
  for (a in 1:n_generic) {
    lines <- c(lines, sprintf('  randcoeff[["b_x%d"]] <- mu_x%d + exp(sigma_x%d) * draws_x%d', a, a, a, a))
  }
  lines <- c(lines, '  randcoeff[["b_price"]] <- -exp(mu_price + exp(sigma_price) * draws_price)')
  lines <- c(lines, '  return(randcoeff)', '}')
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

# Helper: dynamically build apollo_probabilities function for MMNL
.make_apollo_prob_mmnl <- function(dgp = DGP_DEFAULT) {
  J <- dgp$n_alternatives; n_generic <- dgp$n_generic
  alt_entries <- paste(sprintf('alt%d = %d', 1:J, 1:J), collapse = ", ")
  avail_entries <- paste(sprintf('alt%d = 1', 1:J), collapse = ", ")

  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    '  V <- list()'
  )
  for (j in 1:J) {
    terms <- c()
    if (j < J) terms <- c(terms, sprintf('asc_alt%d', j))
    for (a in 1:n_generic) terms <- c(terms, sprintf('b_x%d * x%d_%d', a, a, j))
    terms <- c(terms, sprintf('b_price * price_%d', j))
    lines <- c(lines, sprintf('  V[["alt%d"]] <- %s', j, paste(terms, collapse = " + ")))
  }
  lines <- c(lines,
    sprintf('  mnl_settings <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE, utilities = V)',
            alt_entries, avail_entries),
    '  P[["model"]] <- apollo_mnl(mnl_settings, functionality)',
    '  P <- apollo_panelProd(P, apollo_inputs, functionality)',
    '  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)',
    '  return(P)',
    '}'
  )
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

# Helper: run one Apollo MMNL estimation with given draws and starting values
.run_apollo_mmnl <- function(database, n_draws, start_beta,
                             dgp                = DGP_DEFAULT,
                             n_cores            = NULL,
                             draws_type         = DRAWS_TYPE_MMNL,
                             estimation_routine = ESTIMATION_ROUTINE_MMNL,
                             bounds             = NULL) {
  cleanup_apollo()
  n_generic <- dgp$n_generic

  if (is.null(n_cores)) {
    n_cores <- getOption("klue.mmnl.n_cores", .klue_default_mmnl_cores())
  }
  n_cores <- max(1L, as.integer(n_cores))

  apollo_control <<- list(
    modelName       = paste0("MMNL_sim_", as.integer(Sys.time()) %% 100000,
                             "_", sample.int(10000, 1)),
    modelDescr      = "MMNL simulation",
    indivID         = "ID",
    nCores          = n_cores,
    mixing          = TRUE,
    outputDirectory = tempdir()
  )

  apollo_beta  <<- start_beta
  apollo_fixed <<- c()

  draw_names <- c(paste0("draws_x", 1:n_generic), "draws_price")
  apollo_draws <<- list(
    interDrawsType = draws_type,
    interNDraws    = as.integer(n_draws),
    interUnifDraws = c(),
    interNormDraws = draw_names,
    intraDrawsType = draws_type,
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )

  apollo_randCoeff <<- .make_apollo_randCoeff(dgp)
  apollo_probabilities <<- .make_apollo_prob_mmnl(dgp)

  if (exists("apollo_lcPars", envir = .GlobalEnv)) rm("apollo_lcPars", envir = .GlobalEnv)

  # NOTE: apollo_probabilities is set *before* apollo_validateInputs because
  # validateInputs runs an internal pre-processing check that inspects the
  # current apollo_probabilities in globalenv. If it isn't set yet, Apollo
  # emits "WARNING: The pre-processing of 'apollo_probabilities' failed in
  # initial testing." and the subsequent apollo_estimate call fails for
  # certain (J, n_generic) combinations.
  apollo_inputs <<- tryCatch(
    apollo_validateInputs(
      apollo_beta    = apollo_beta,
      apollo_fixed   = apollo_fixed,
      database       = database,
      apollo_control = apollo_control
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl: apollo_validateInputs error: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(apollo_inputs)) return(NULL)

  est_settings <- list(
    estimationRoutine = estimation_routine,
    writeIter         = FALSE,
    silent            = TRUE
  )
  if (!is.null(bounds)) {
    est_settings$bounds <- bounds
  }

  tryCatch(
    apollo_estimate(
      apollo_beta          = apollo_beta,
      apollo_fixed         = apollo_fixed,
      apollo_probabilities = apollo_probabilities,
      apollo_inputs        = apollo_inputs,
      estimate_settings    = est_settings
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl: apollo_estimate error: ",
              conditionMessage(e))
      NULL
    }
  )
}

klue_mmnl <- function(database,
                          n_draws            = NULL,
                          n_draws_stage1     = NULL,
                          draws_type         = NULL,
                          estimation_routine = NULL,
                          n_cores            = NULL,
                          quiet              = NULL,
                          mu_price_bounds    = NULL,
                          sigma_price_bounds = NULL,
                          dgp                = DGP_DEFAULT) {
  d <- klue_mmnl_defaults()
  if (is.null(n_draws))            n_draws            <- d$n_draws
  if (is.null(n_draws_stage1))     n_draws_stage1     <- d$n_draws_stage1
  if (is.null(draws_type))         draws_type         <- d$draws_type
  if (is.null(estimation_routine)) estimation_routine <- d$estimation_routine
  if (is.null(n_cores))            n_cores            <- d$n_cores
  if (is.null(quiet))              quiet              <- d$quiet
  # Bounds: pass NA explicitly to disable; NULL means "use defaults".
  if (is.null(mu_price_bounds))    mu_price_bounds    <- d$mu_price_bounds
  if (is.null(sigma_price_bounds)) sigma_price_bounds <- d$sigma_price_bounds
  if (length(mu_price_bounds)    == 1 && is.na(mu_price_bounds))    mu_price_bounds    <- NULL
  if (length(sigma_price_bounds) == 1 && is.na(sigma_price_bounds)) sigma_price_bounds <- NULL

  cleanup_apollo()

  # If `quiet`, redirect Apollo's stdout to a tempfile so we can attach the
  # last 40 lines to the result when something fails. The file is NOT deleted
  # on exit — its path is returned with the fail result so users can inspect.
  log_file <- tempfile(pattern = "klue_mmnl_", fileext = ".log")
  if (isTRUE(quiet)) {
    old_sink_out <- sink.number()
    old_sink_msg <- sink.number(type = "message")
    msg_con <- file(log_file, open = "at")
    sink(log_file)
    sink(msg_con, type = "message")
    on.exit({
      while (sink.number() > old_sink_out) sink()
      while (sink.number(type = "message") > old_sink_msg) sink(type = "message")
      try(close(msg_con), silent = TRUE)
      cleanup_apollo()
    }, add = TRUE)
  } else {
    on.exit(cleanup_apollo(), add = TRUE)
  }

  N <- length(unique(database$ID))
  n_beta <- dgp$n_beta; n_generic <- dgp$n_generic; J <- dgp$n_alternatives

  read_log_tail <- function() {
    if (!isTRUE(quiet) || !file.exists(log_file)) return(NULL)
    out <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) NULL)
    if (length(out) == 0) NULL else utils::tail(out, 40L)
  }

  fail_with <- function(reason) {
    list(converged = FALSE, LL = -Inf, BIC = Inf, AIC = Inf,
         k = 0, mu = rep(0, n_beta), sigma = rep(0, n_beta),
         reason          = reason,
         apollo_log_tail = read_log_tail(),
         apollo_log_path = if (isTRUE(quiet) && file.exists(log_file)) log_file else NULL,
         settings        = list(n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
                                draws_type = draws_type,
                                estimation_routine = estimation_routine,
                                n_cores = n_cores,
                                mu_price_bounds = mu_price_bounds,
                                sigma_price_bounds = sigma_price_bounds))
  }

  # Build a bounds list keyed by name for whichever start vector is in play.
  make_bounds <- function(beta_vec) {
    if (is.null(mu_price_bounds) && is.null(sigma_price_bounds)) return(NULL)
    lower <- rep(-Inf, length(beta_vec))
    upper <- rep( Inf, length(beta_vec))
    names(lower) <- names(upper) <- names(beta_vec)
    if (!is.null(mu_price_bounds) && "mu_price" %in% names(beta_vec)) {
      lower["mu_price"] <- mu_price_bounds[1]
      upper["mu_price"] <- mu_price_bounds[2]
    }
    if (!is.null(sigma_price_bounds) && "sigma_price" %in% names(beta_vec)) {
      lower["sigma_price"] <- sigma_price_bounds[1]
      upper["sigma_price"] <- sigma_price_bounds[2]
    }
    list(lower = lower, upper = upper)
  }

  # MNL-informed starting values (MNL = LCMNL with C=1)
  mnl_starts <- matrix(0, nrow = 1, ncol = n_beta)
  mnl_fit <- tryCatch(
    estimate_lcmnl(database, C = 1, start_betas = mnl_starts, dgp = dgp),
    error = function(e) NULL
  )

  # Build starting beta vector dynamically
  beta0 <- c()
  for (j in 1:(J - 1)) beta0[paste0("asc_alt", j)] <- 0
  if (!is.null(mnl_fit) && mnl_fit$converged) {
    mnl_b <- mnl_fit$betas[1, ]
    for (a in 1:n_generic) beta0[paste0("mu_x", a)] <- mnl_b[a]
    beta0["mu_price"] <- if (mnl_b[n_beta] < 0) log(-mnl_b[n_beta]) else 0.0
  } else {
    for (a in 1:n_generic) beta0[paste0("mu_x", a)] <- 0.5
    beta0["mu_price"] <- 0.0
  }
  for (a in 1:n_generic) beta0[paste0("sigma_x", a)] <- log(0.5)
  beta0["sigma_price"] <- log(0.3)

  # Stage 1: cheap warm-start. If it returns non-finite estimates or pins to a
  # bound (a sign the optimiser failed to converge sensibly) we fall back to
  # the MNL-derived beta0 for stage 2.
  stage1 <- .run_apollo_mmnl(database, n_draws_stage1, beta0,
                             dgp                = dgp,
                             n_cores            = n_cores,
                             draws_type         = draws_type,
                             estimation_routine = estimation_routine,
                             bounds             = make_bounds(beta0))

  beta1 <- beta0
  if (!is.null(stage1) && !is.null(stage1$estimate)) {
    est_s1 <- stage1$estimate
    ok <- all(is.finite(est_s1))
    if (ok) {
      bnd <- make_bounds(beta0)
      if (!is.null(bnd)) {
        rng <- bnd$upper - bnd$lower
        rng[!is.finite(rng)] <- 0
        slack <- 1e-3 * rng
        ok <- all(est_s1 >= bnd$lower + slack & est_s1 <= bnd$upper - slack,
                  na.rm = TRUE)
      }
    }
    if (ok) beta1 <- est_s1
  }

  # Stage 2: main estimation with full draws.
  model <- .run_apollo_mmnl(database, n_draws, beta1,
                            dgp                = dgp,
                            n_cores            = n_cores,
                            draws_type         = draws_type,
                            estimation_routine = estimation_routine,
                            bounds             = make_bounds(beta1))
  if (is.null(model))                                return(fail_with("stage2_apollo_estimate_failed"))
  if (is.null(model$estimate) ||
      !all(is.finite(model$estimate)))               return(fail_with("stage2_non_finite_estimates"))
  if (is.null(model$LLout) || !is.finite(model$LLout[1]))
                                                     return(fail_with("stage2_non_finite_LL"))

  est    <- model$estimate
  LL     <- model$LLout[1]
  n_free <- length(est)
  BIC    <- -2 * LL + n_free * log(N)
  AIC    <- -2 * LL + 2 * n_free

  mu_names    <- c(paste0("mu_x", 1:n_generic), "mu_price")
  sigma_names <- c(paste0("sigma_x", 1:n_generic), "sigma_price")
  list(converged = TRUE, LL = LL, BIC = BIC, AIC = AIC, k = n_free,
       mu = est[mu_names], sigma = exp(est[sigma_names]),
       reason = "ok", apollo_log_tail = NULL, apollo_log_path = NULL,
       settings = list(n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
                       draws_type = draws_type,
                       estimation_routine = estimation_routine,
                       n_cores = n_cores,
                       mu_price_bounds = mu_price_bounds,
                       sigma_price_bounds = sigma_price_bounds))
}

# =============================================================================
# SECTION 5: METRICS
# =============================================================================

compute_ari <- function(true_labels, pred_labels) {
  cont <- table(true_labels, pred_labels)
  a <- rowSums(cont); b <- colSums(cont); n <- sum(cont)
  sum_nij2 <- sum(cont * (cont - 1)) / 2
  sum_a2   <- sum(a * (a - 1)) / 2
  sum_b2   <- sum(b * (b - 1)) / 2
  expected <- sum_a2 * sum_b2 / (n * (n - 1) / 2)
  max_idx  <- (sum_a2 + sum_b2) / 2
  if (max_idx == expected) return(1)
  (sum_nij2 - expected) / (max_idx - expected)
}

compute_recovery <- function(true_betas, est_betas) {
  K <- nrow(true_betas); Ce <- nrow(est_betas)
  if (K != Ce) return(list(rmse = NA, bias = NA))
  if (K == 1) {
    diffs <- true_betas - est_betas
    return(list(rmse = sqrt(mean(diffs^2)), bias = mean(diffs)))
  }

  # Optimal permutation via exhaustive search (feasible for K <= 8, i.e. 40320 perms)
  all_perms <- function(n) {
    if (n == 1) return(list(1L))
    smaller <- all_perms(n - 1)
    result <- vector("list", factorial(n))
    idx <- 0
    for (p in smaller) {
      for (pos in 1:n) {
        idx <- idx + 1
        result[[idx]] <- append(p, n, after = pos - 1)
      }
    }
    result
  }

  perms <- all_perms(K)
  best_cost <- Inf
  best_perm <- seq_len(K)
  for (p in perms) {
    cost <- sum((true_betas - est_betas[p, , drop = FALSE])^2)
    if (cost < best_cost) {
      best_cost <- cost
      best_perm <- p
    }
  }
  diffs <- true_betas - est_betas[best_perm, , drop = FALSE]
  list(rmse = sqrt(mean(diffs^2)), bias = mean(diffs))
}

# =============================================================================
# SECTION 6: MAIN SIMULATION
# =============================================================================

run_main_simulation <- function(true_K_values = c(1, 2, 3, 4, 5),
                                kappa_values  = c(0.5, 0.75, 1.0, 1.25, 1.5),
                                sigma_values  = c(0.1, 0.15, 0.2, 0.25),
                                n_reps = 5, C_cands = 1:6,
                                dgp = DGP_DEFAULT, sep_profile = NULL,
                                verbose = TRUE) {
  # K=1: kappa irrelevant (no segments), fix at 0 to avoid redundant runs
  if (1L %in% true_K_values) {
    conds_k1 <- expand.grid(true_K = 1, kappa = 0,
                            sigma = sigma_values, rep = 1:n_reps)
    conds_kn <- expand.grid(true_K = setdiff(true_K_values, 1),
                            kappa = kappa_values,
                            sigma = sigma_values, rep = 1:n_reps)
    conditions <- rbind(conds_k1, conds_kn)
  } else {
    conditions <- expand.grid(true_K = true_K_values, kappa = kappa_values,
                              sigma = sigma_values, rep = 1:n_reps)
  }
  nc <- nrow(conditions)
  if (verbose) {
    cat(sprintf("==== MAIN SIMULATION (%d conditions, %d cores, %d attrs) ====\n",
                nc, N_CORES_LCMNL, dgp$n_beta))
  }

  .run_one_main <- function(i) {
    tK  <- conditions$true_K[i]; kap <- conditions$kappa[i]
    sig <- conditions$sigma[i];  rp  <- conditions$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)

    npc <- if (tK == 1L) 300L else 150L
    data <- klue_simulate(N_per_class = npc, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp, sep_profile = sep_profile)

    bics <- rep(Inf, length(C_cands))
    aics <- rep(Inf, length(C_cands))
    icls <- rep(Inf, length(C_cands))
    models <- list()
    for (j in seq_along(C_cands)) {
      m <- klue_lcmnl(data$database, C_cands[j], dgp = dgp)
      models[[j]] <- m
      if (m$converged) { bics[j] <- m$BIC; aics[j] <- m$AIC; icls[j] <- m$ICL }
    }

    sb <- C_cands[which.min(bics)]; sa <- C_cands[which.min(aics)]
    si <- C_cands[which.min(icls)]
    bic_best_idx <- which.min(bics)
    bic_method <- models[[bic_best_idx]]$best_method
    icl_bic_val <- models[[bic_best_idx]]$ICL_BIC

    ari_val <- NA
    tK_idx <- which(C_cands == tK)
    if (length(tK_idx) == 1 && models[[tK_idx]]$converged && tK > 1) {
      pred <- apply(models[[tK_idx]]$posteriors, 1, which.max)
      ari_val <- compute_ari(data$true_class, pred)
    }

    rmse_val <- NA; bias_val <- NA
    if (sb == tK && length(tK_idx) == 1 && models[[tK_idx]]$converged && tK > 1) {
      rec <- compute_recovery(data$true_betas, models[[tK_idx]]$betas)
      rmse_val <- rec$rmse; bias_val <- rec$bias
    }

    list(true_K = tK, kappa = kap, sigma = sig, rep = rp,
         selected_bic = sb, selected_aic = sa, selected_icl = si,
         bic_correct = as.integer(sb == tK),
         aic_correct = as.integer(sa == tK),
         icl_correct = as.integer(si == tK),
         best_method = ifelse(is.null(bic_method) || is.na(bic_method),
                              NA_character_, bic_method),
         ari = ari_val, rmse = rmse_val, bias = bias_val,
         icl_bic = icl_bic_val)
  }

  batch_size <- N_CORES_LCMNL * 3L
  res_list <- vector("list", nc)
  done <- 0L

  for (b_start in seq(1L, nc, by = batch_size)) {
    b_end <- min(b_start + batch_size - 1L, nc)
    idx <- b_start:b_end
    batch <- mclapply(idx, .run_one_main, mc.cores = N_CORES_LCMNL)
    res_list[idx] <- batch
    done <- b_end
    if (verbose) {
      for (i in idx) {
        r <- res_list[[i]]
        status <- ifelse(r$bic_correct == 1L, "OK", "X")
        meth <- ifelse(is.na(r$best_method), "?", r$best_method)
        c_label <- if (r$selected_bic == 1L) "C=1(MNL)" else paste0("C=", r$selected_bic)
        cat(sprintf("  [%3d/%d] K=%d kappa=%.2f sigma=%.2f rep=%d ... BIC->%s %s [%s]\n",
                    i, nc, r$true_K, r$kappa, r$sigma, r$rep,
                    c_label, status, meth))
      }
    }
  }

  results <- data.frame(
    condition    = 1:nc,
    true_K       = sapply(res_list, `[[`, "true_K"),
    kappa        = sapply(res_list, `[[`, "kappa"),
    sigma        = sapply(res_list, `[[`, "sigma"),
    rep          = sapply(res_list, `[[`, "rep"),
    selected_bic = sapply(res_list, `[[`, "selected_bic"),
    selected_aic = sapply(res_list, `[[`, "selected_aic"),
    selected_icl = sapply(res_list, `[[`, "selected_icl"),
    bic_correct  = sapply(res_list, `[[`, "bic_correct"),
    aic_correct  = sapply(res_list, `[[`, "aic_correct"),
    icl_correct  = sapply(res_list, `[[`, "icl_correct"),
    best_method  = sapply(res_list, `[[`, "best_method"),
    ari          = sapply(res_list, `[[`, "ari"),
    rmse         = sapply(res_list, `[[`, "rmse"),
    bias         = sapply(res_list, `[[`, "bias"),
    icl_bic      = sapply(res_list, `[[`, "icl_bic"),
    stringsAsFactors = FALSE
  )
  results
}

summarise_main_results <- function(df) {
  cat("\n==== RESULTS ====\n")
  cat(sprintf("  BIC: %.1f%%  AIC: %.1f%%  ICL: %.1f%%\n",
              100 * mean(df$bic_correct), 100 * mean(df$aic_correct),
              100 * mean(df$icl_correct)))
  cat("\n--- BY K ---\n")
  for (k in sort(unique(df$true_K)))
    cat(sprintf("  K=%d: %.1f%%\n", k, 100 * mean(df$bic_correct[df$true_K == k])))
  cat("\n--- BY kappa ---\n")
  for (kap in sort(unique(df$kappa)))
    cat(sprintf("  kappa=%.2f: %.1f%%\n", kap, 100 * mean(df$bic_correct[df$kappa == kap])))
  v <- !is.na(df$rmse)
  if (any(v))
    cat(sprintf("\n  RMSE: %.4f  Bias: %.4f  ARI: %.3f\n",
                mean(df$rmse[v]), mean(df$bias[v]), mean(df$ari[!is.na(df$ari)])))
}

# =============================================================================
# SECTION 7: MMNL COMPARISON
# =============================================================================

run_mmnl_comparison <- function(n_cond = 80, n_draws = N_DRAWS_MMNL,
                                C_cands = 1:5, verbose = TRUE,
                                dgp = DGP_DEFAULT) {
  # K=1: pure continuous heterogeneity (no segments); kappa irrelevant, fix at 0.
  # Include sigma=0.35 to test higher heterogeneity where MMNL advantage is clearest.
  conds_k1 <- expand.grid(true_K = 1, kappa = 0,
                           sigma = c(0.15, 0.25, 0.35), rep = 1:3)
  conds_kn <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                           sigma = c(0.15, 0.25), rep = 1:3)
  conds <- rbind(conds_k1, conds_kn)
  conds <- conds[1:min(n_cond, nrow(conds)), ]
  nc <- nrow(conds)
  if (verbose) cat(sprintf("\n==== MNL vs LCMNL vs MMNL (%d conditions, incl. K=1) ====\n", nc))

  # Phase 1: Run LCMNL in parallel (no Apollo global state)
  if (verbose) cat("  Phase 1: LCMNL estimation (parallel)...\n")
  .run_lcmnl_part <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)
    npc <- if (tK == 1L) 300L else 150L
    data <- klue_simulate(N_per_class = npc, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    mnl_bic <- Inf  # C=1 result, i.e. standard MNL
    best_lc_bic <- Inf; best_lc_C <- NA_integer_; best_lc_method <- NA_character_
    for (Cc in C_cands) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged) {
        if (Cc == 1L) mnl_bic <- m$BIC  # MNL = LCMNL with C=1
        if (Cc >= 2L && m$BIC < best_lc_bic) {
          best_lc_bic <- m$BIC; best_lc_C <- Cc
          best_lc_method <- m$best_method
        }
      }
    }
    list(tK = tK, kap = kap, sig = sig, seed = seed, npc = npc,
         mnl_bic = mnl_bic, lc_bic = best_lc_bic,
         lc_C = best_lc_C, lc_method = best_lc_method)
  }
  lc_results <- mclapply(1:nc, .run_lcmnl_part, mc.cores = N_CORES_LCMNL)

  # Phase 2: Run MMNL sequentially (Apollo global state — must be serial)
  if (verbose) cat("  Phase 2: MMNL estimation (sequential, Apollo)...\n")
  rK <- integer(nc); rkap <- numeric(nc); rsig <- numeric(nc)
  rmnl <- rep(NA, nc); rlc <- rep(NA, nc); rmm <- rep(NA, nc)
  rlc_C <- rep(NA_integer_, nc)
  rlc_method <- rep(NA_character_, nc); valid <- logical(nc)

  for (i in 1:nc) {
    lr <- lc_results[[i]]
    rK[i] <- lr$tK; rkap[i] <- lr$kap; rsig[i] <- lr$sig

    if (verbose) cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... ", i, nc, lr$tK, lr$kap, lr$sig))

    # Regenerate data with same seed for MMNL
    data <- klue_simulate(N_per_class = lr$npc, T_tasks = 20, true_K = lr$tK,
                          separation = lr$kap, heterogeneity = lr$sig, seed = lr$seed,
                          dgp = dgp)
    mm <- klue_mmnl(data$database, n_draws = n_draws, dgp = dgp)

    if (is.finite(lr$mnl_bic) && mm$converged) {
      rmnl[i] <- lr$mnl_bic; rlc[i] <- lr$lc_bic; rmm[i] <- mm$BIC
      rlc_C[i] <- lr$lc_C; rlc_method[i] <- lr$lc_method
      valid[i] <- TRUE
      bics <- c(MNL = unname(lr$mnl_bic), LCMNL = unname(lr$lc_bic), MMNL = unname(mm$BIC))
      winner <- names(which.min(bics))
      meth <- ifelse(is.na(lr$lc_method), "-", lr$lc_method)
      if (!is.finite(lr$lc_bic)) {
        lc_str <- "LCMNL(-)"
      } else {
        lc_str <- if (lr$lc_C == 1L) "MNL(C=1)" else sprintf("LCMNL(C=%d)", lr$lc_C)
      }
      if (verbose) cat(sprintf("%s -> %s [%s]\n", lc_str, winner, meth))
    } else {
      if (verbose) cat("FAILED\n")
    }
  }

  idx <- which(valid)
  bic_winner <- mapply(function(mnl, lc, mm) {
    bics <- c(MNL = mnl, LCMNL = lc, MMNL = mm)
    names(which.min(bics))
  }, rmnl[idx], rlc[idx], rmm[idx])
  df <- data.frame(true_K = rK[idx], kappa = rkap[idx], sigma = rsig[idx],
                   mnl_BIC = rmnl[idx],
                   lcmnl_BIC = rlc[idx], lcmnl_C = rlc_C[idx],
                   lcmnl_method = rlc_method[idx],
                   mmnl_BIC = rmm[idx],
                   bic_prefers = bic_winner,
                   stringsAsFactors = FALSE)
  if (verbose && nrow(df) > 0) {
    n_tot <- nrow(df)
    n_mnl  <- sum(df$bic_prefers == "MNL")
    n_lc   <- sum(df$bic_prefers == "LCMNL")
    n_mm   <- sum(df$bic_prefers == "MMNL")
    cat(sprintf("  Overall: MNL %d/%d (%.0f%%)  LCMNL %d/%d (%.0f%%)  MMNL %d/%d (%.0f%%)\n",
                n_mnl, n_tot, 100*n_mnl/n_tot,
                n_lc,  n_tot, 100*n_lc/n_tot,
                n_mm,  n_tot, 100*n_mm/n_tot))
    k1 <- df[df$true_K == 1, ]
    if (nrow(k1) > 0) {
      cat(sprintf("  K=1 (continuous DGP): MNL %d  LCMNL %d  MMNL %d  (of %d)\n",
                  sum(k1$bic_prefers == "MNL"),
                  sum(k1$bic_prefers == "LCMNL"),
                  sum(k1$bic_prefers == "MMNL"), nrow(k1)))
    }
    kn <- df[df$true_K > 1, ]
    if (nrow(kn) > 0) {
      cat(sprintf("  K>1 (discrete DGP):   MNL %d  LCMNL %d  MMNL %d  (of %d)\n",
                  sum(kn$bic_prefers == "MNL"),
                  sum(kn$bic_prefers == "LCMNL"),
                  sum(kn$bic_prefers == "MMNL"), nrow(kn)))
    }
  }
  df
}

# =============================================================================
# SECTION 8: SUPPLEMENTARY ANALYSES
# =============================================================================

run_convergence_ablation <- function(n_random = 50, n_cond = 40, verbose = TRUE,
                                     dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 1: CONVERGENCE ====\n")
  conds <- expand.grid(true_K = c(3, 4, 5), kappa = c(0.50, 0.75, 1.00),
                       sigma = c(0.15, 0.25), rep = 1:2)
  conds <- conds[1:min(n_cond, nrow(conds)), ]

  .run_one_conv <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)

    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)

    t1 <- system.time(cres <- klue_lcmnl(data$database, tK, dgp = dgp))[3]

    # Uninformed random starts: N(0,2) for attributes, -exp(N(0,1)) for price
    rLLs <- rep(-Inf, n_random)
    t2 <- system.time({
      for (r in 1:n_random) {
        set.seed(seed * 1000 + r)
        rb <- matrix(rnorm(tK * dgp$n_generic, mean = 0, sd = 2), nrow = tK, ncol = dgp$n_generic)
        rb <- cbind(rb, -exp(rnorm(tK, mean = 0, sd = 1)))
        rr <- estimate_lcmnl(data$database, tK, start_betas = rb, dgp = dgp)
        if (rr$converged) rLLs[r] <- rr$LL
      }
    })[3]

    clust_LL <- if (cres$converged) cres$LL else -Inf
    n_conv <- sum(rLLs > -Inf)
    n_global <- sum(abs(rLLs - clust_LL) < 0.1)
    n_local <- sum(rLLs > -Inf & rLLs < clust_LL - 0.1)
    match_at <- NA_integer_
    for (r in 1:n_random) {
      if (rLLs[r] >= clust_LL - 0.1) { match_at <- r; break }
    }

    if (cres$converged) {
      list(valid = TRUE, K = tK, kappa = kap,
           cluster_best = is.na(match_at) || abs(clust_LL - max(rLLs)) < 0.1,
           pct_global = n_global / n_random,
           pct_local = n_local / n_random,
           match_at = match_at,
           time_ratio = t2 / max(t1, 0.01))
    } else {
      list(valid = FALSE)
    }
  }

  res_list <- mclapply(1:nrow(conds), .run_one_conv, mc.cores = N_CORES_LCMNL)

  valid <- sapply(res_list, `[[`, "valid")
  vres <- res_list[valid]
  cluster_best <- sapply(vres, `[[`, "cluster_best")
  pct_global   <- sapply(vres, `[[`, "pct_global")
  pct_local    <- sapply(vres, `[[`, "pct_local")
  match_at     <- sapply(vres, `[[`, "match_at")
  time_ratio   <- sapply(vres, `[[`, "time_ratio")
  Ks           <- sapply(vres, `[[`, "K")

  df <- data.frame(K = Ks, cluster_best = cluster_best,
                   pct_global = pct_global, pct_local = pct_local,
                   match_at = match_at, time_ratio = time_ratio)
  if (verbose && nrow(df) > 0) {
    matched <- !is.na(df$match_at)
    cat(sprintf("  Clust optimal: %.1f%%  Global: %.1f%%  Local: %.1f%%\n",
                100 * mean(df$cluster_best), 100 * mean(df$pct_global),
                100 * mean(df$pct_local)))
    cat(sprintf("  Match: %d/%d  Median starts: %s\n",
                sum(matched), nrow(df),
                ifelse(any(matched), as.character(median(df$match_at[matched])), "N/A")))
  }
  df
}

# Three-arm initialisation ablation: RP contrasts vs one-hot choice indicators
# vs uninformed random starts, on the H1 convergence design. The "global"
# log-likelihood for a condition is the max across all three arms (within tol).
run_initialisation_ablation <- function(n_random = 50, n_cond = 40,
                                        verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== INITIALISATION ABLATION (RP / one-hot / random) ====\n")
  conds <- expand.grid(true_K = c(3, 4, 5), kappa = c(0.50, 0.75, 1.00),
                       sigma = c(0.15, 0.25), rep = 1:2)
  conds <- conds[1:min(n_cond, nrow(conds)), ]

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)

    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)

    # Arm A: RP-contrast clustering (main workflow)
    rp_res <- tryCatch(
      klue_lcmnl(data$database, tK, dgp = dgp),
      error = function(e) NULL
    )
    rp_LL <- if (!is.null(rp_res) && isTRUE(rp_res$converged)) rp_res$LL else -Inf

    # Arm B: One-hot choice indicators clustering
    oh_res <- tryCatch(
      estimate_lcmnl_multistart_onehot(data$database, tK, dgp = dgp),
      error = function(e) NULL
    )
    oh_LL <- if (!is.null(oh_res) && isTRUE(oh_res$converged)) oh_res$LL else -Inf

    # Arm C: Uninformed random starts
    rLLs <- rep(-Inf, n_random)
    for (r in 1:n_random) {
      set.seed(seed * 1000 + r)
      rb <- matrix(rnorm(tK * dgp$n_generic, mean = 0, sd = 2),
                   nrow = tK, ncol = dgp$n_generic)
      rb <- cbind(rb, -exp(rnorm(tK, mean = 0, sd = 1)))
      rr <- tryCatch(
        estimate_lcmnl(data$database, tK, start_betas = rb, dgp = dgp),
        error = function(e) NULL
      )
      if (!is.null(rr) && isTRUE(rr$converged)) rLLs[r] <- rr$LL
    }

    # Global LL = best across the three arms; tolerance for "at global"
    best_LL <- max(c(rp_LL, oh_LL, rLLs, -Inf))
    tol <- 0.1
    any_converged <- (rp_LL > -Inf) || (oh_LL > -Inf) || any(rLLs > -Inf)

    list(valid = any_converged,
         K = tK, kappa = kap, sigma = sig,
         rp_at_global = (rp_LL >= best_LL - tol),
         oh_at_global = (oh_LL >= best_LL - tol),
         random_pct_global = mean(rLLs >= best_LL - tol),
         rp_LL = rp_LL, oh_LL = oh_LL, best_LL = best_LL,
         gap_rp_oh = rp_LL - oh_LL)
  }

  res_list <- mclapply(1:nrow(conds), .run_one, mc.cores = N_CORES_LCMNL)

  valid <- sapply(res_list, function(x) isTRUE(x$valid))
  vres <- res_list[valid]
  df <- data.frame(
    K                 = sapply(vres, `[[`, "K"),
    kappa             = sapply(vres, `[[`, "kappa"),
    sigma             = sapply(vres, `[[`, "sigma"),
    rp_at_global      = sapply(vres, `[[`, "rp_at_global"),
    oh_at_global      = sapply(vres, `[[`, "oh_at_global"),
    random_pct_global = sapply(vres, `[[`, "random_pct_global"),
    rp_LL             = sapply(vres, `[[`, "rp_LL"),
    oh_LL             = sapply(vres, `[[`, "oh_LL"),
    best_LL           = sapply(vres, `[[`, "best_LL"),
    gap_rp_oh         = sapply(vres, `[[`, "gap_rp_oh")
  )

  if (verbose && nrow(df) > 0) {
    cat(sprintf("  RP contrasts reach global:    %.1f%% (%d/%d conditions)\n",
                100 * mean(df$rp_at_global),
                sum(df$rp_at_global), nrow(df)))
    cat(sprintf("  One-hot indicators at global: %.1f%% (%d/%d)\n",
                100 * mean(df$oh_at_global),
                sum(df$oh_at_global), nrow(df)))
    cat(sprintf("  Random starts at global:      %.1f%% per start\n",
                100 * mean(df$random_pct_global)))
    by_K <- aggregate(cbind(rp_at_global, oh_at_global, random_pct_global)
                      ~ K, data = df, FUN = mean)
    cat("\n  By true K:\n")
    print(by_K, row.names = FALSE)
  }
  df
}

run_unbalanced_analysis <- function(verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 3: UNBALANCED ====\n")
  configs <- list(
    list(name = "balanced", props = c(1/3, 1/3, 1/3)),
    list(name = "mild",     props = c(0.5, 0.3, 0.2)),
    list(name = "moderate", props = c(0.6, 0.25, 0.15)),
    list(name = "severe",   props = c(0.7, 0.2, 0.1))
  )
  conds <- expand.grid(kappa = c(0.75, 1.0, 1.25), rep = 1:5)

  # Flatten all config x condition combos for parallel dispatch
  all_jobs <- expand.grid(ci = seq_along(configs), cond_i = 1:nrow(conds))
  .run_one_unbal <- function(row) {
    ci <- all_jobs$ci[row]; cfg <- configs[[ci]]
    kap <- conds$kappa[all_jobs$cond_i[row]]
    rp  <- conds$rep[all_jobs$cond_i[row]]
    seed <- as.integer(3000 + 1000 * ci + 100 * (kap * 100) + rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 3,
                          separation = kap, heterogeneity = 0.2, seed = seed,
                          class_proportions = cfg$props, dgp = dgp)
    bics <- rep(Inf, 5)
    for (Cc in 1:5) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged) bics[Cc] <- m$BIC
    }
    list(ci = ci, correct = as.integer(which.min(bics) == 3))
  }

  res_list <- mclapply(1:nrow(all_jobs), .run_one_unbal, mc.cores = N_CORES_LCMNL)

  config_names <- character(0); accuracies <- numeric(0)
  for (ci in seq_along(configs)) {
    hits <- sapply(res_list[sapply(res_list, `[[`, "ci") == ci], `[[`, "correct")
    acc <- 100 * mean(hits)
    config_names <- c(config_names, configs[[ci]]$name)
    accuracies   <- c(accuracies, acc)
    if (verbose) cat(sprintf("  %s: %.1f%%\n", configs[[ci]]$name, acc))
  }
  data.frame(config = config_names, accuracy = accuracies)
}

run_design_comparison <- function(verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 4: D-EFFICIENT ====\n")
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0), rep = 1:5)

  .run_one_design <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]; rp <- conds$rep[i]
    seed <- as.integer(4000 + 100 * tK + 10 * (kap * 100) + rp)
    dr <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                        separation = kap, heterogeneity = 0.2, seed = seed,
                        dgp = dgp)
    dd <- klue_simulate_deff(N_per_class = 150, T_tasks = 20, true_K = tK,
                                   separation = kap, heterogeneity = 0.2, seed = seed,
                                   dgp = dgp)
    br <- rep(Inf, 5); bd <- rep(Inf, 5)
    for (Cc in 1:5) {
      m <- klue_lcmnl(dr$database, Cc, dgp = dgp)
      if (m$converged) br[Cc] <- m$BIC
      m <- klue_lcmnl(dd$database, Cc, dgp = dgp)
      if (m$converged) bd[Cc] <- m$BIC
    }
    list(r_ok = as.integer(which.min(br) == tK),
         d_ok = as.integer(which.min(bd) == tK))
  }

  res_list <- mclapply(1:nrow(conds), .run_one_design, mc.cores = N_CORES_LCMNL)
  rok <- sum(sapply(res_list, `[[`, "r_ok"))
  dok <- sum(sapply(res_list, `[[`, "d_ok"))

  if (verbose) cat(sprintf("  Random: %.1f%%  D-eff: %.1f%%\n",
                            100 * rok / nrow(conds), 100 * dok / nrow(conds)))
  data.frame(random = 100 * rok / nrow(conds), deff = 100 * dok / nrow(conds))
}

run_concomitant_analysis <- function(verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 5: CONCOMITANT ====\n")
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0, 1.25),
                       cs = c(0.5, 1.0, 1.5), rep = 1:3)

  .run_one_conc <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    cs <- conds$cs[i]; rp <- conds$rep[i]
    seed <- as.integer(5000 + 100 * tK + 10 * (kap * 100) + rp + cs * 1000)
    data <- klue_simulate_cov(
      N_per_class = 150, T_tasks = 20, true_K = tK, separation = kap,
      heterogeneity = 0.2, seed = seed, covariate_strength = cs,
      dgp = dgp
    )
    bics <- rep(Inf, 5)
    for (Cc in 1:5) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged) bics[Cc] <- m$BIC
    }
    correct <- as.integer(which.min(bics) == tK)
    mt <- klue_lcmnl(data$database, tK, dgp = dgp)
    ari_val <- NA
    if (mt$converged && tK > 1) {
      pred <- apply(mt$posteriors, 1, which.max)
      ari_val <- compute_ari(data$true_class, pred)
    }
    list(correct = correct, ari = ari_val)
  }

  res_list <- mclapply(1:nrow(conds), .run_one_conc, mc.cores = N_CORES_LCMNL)
  ok   <- sum(sapply(res_list, `[[`, "correct"))
  aris <- na.omit(sapply(res_list, `[[`, "ari"))

  if (verbose) cat(sprintf("  Accuracy: %.1f%%  ARI: %.3f\n",
                            100 * ok / nrow(conds), mean(aris)))
  data.frame(accuracy = 100 * ok / nrow(conds), mean_ari = mean(aris))
}

run_unconditional_recovery <- function(n_cond = 80, verbose = TRUE,
                                       dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 6: RECOVERY ====\n")
  conds <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                       sigma = c(0.15, 0.25), rep = 1:3)
  conds <- conds[1:min(n_cond, nrow(conds)), ]

  .run_one_recov <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    m <- klue_lcmnl(data$database, tK, dgp = dgp)
    if (m$converged && tK > 1) {
      r <- compute_recovery(data$true_betas, m$betas)
      pred <- apply(m$posteriors, 1, which.max)
      list(valid = TRUE, rmse = r$rmse, bias = r$bias,
           ari = compute_ari(data$true_class, pred))
    } else {
      list(valid = FALSE)
    }
  }

  res_list <- mclapply(1:nrow(conds), .run_one_recov, mc.cores = N_CORES_LCMNL)
  valid <- sapply(res_list, `[[`, "valid")
  rmses  <- sapply(res_list[valid], `[[`, "rmse")
  biases <- sapply(res_list[valid], `[[`, "bias")
  aris   <- sapply(res_list[valid], `[[`, "ari")

  if (verbose) cat(sprintf("  RMSE: %.4f  Bias: %.4f  ARI: %.3f\n",
                            mean(rmses), mean(biases), mean(aris)))
  data.frame(rmse = mean(rmses), bias = mean(biases), ari = mean(aris))
}

# --- SUPP 7: Clustering Method Comparison ---
run_clustering_comparison <- function(verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 7: CLUSTERING METHODS ====\n")

  method_names <- c("kmeans", "gmm", "hc_ward", "hc_complete", "hc_average", "pam")
  conds <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                       sigma = c(0.15, 0.25), rep = 1:3)
  nc <- nrow(conds)

  .run_one_clust <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(6000 + 1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)

    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)

    all_starts <- get_all_starts(data$database, tK, dgp = dgp)

    m_correct <- setNames(integer(length(method_names)), method_names)
    m_conv    <- setNames(integer(length(method_names)), method_names)
    m_ll      <- setNames(rep(NA_real_, length(method_names)), method_names)

    for (mi in seq_along(method_names)) {
      nm <- method_names[mi]
      if (is.null(all_starts[[nm]])) next

      res <- tryCatch(
        estimate_lcmnl(data$database, tK,
                       start_betas = all_starts[[nm]]$betas,
                       start_shares = all_starts[[nm]]$shares,
                       dgp = dgp),
        error = function(e) NULL
      )
      if (!is.null(res) && res$converged) {
        m_conv[mi] <- 1L
        m_ll[mi] <- res$LL

        bics <- rep(Inf, 5)
        for (Cc in 1:5) {
          if (Cc == tK) {
            bics[Cc] <- res$BIC
          } else {
            s_cc <- tryCatch(
              get(paste0("get_", nm, "_starts"))(data$database, Cc, dgp = dgp),
              error = function(e) NULL
            )
            if (is.null(s_cc)) next
            m_cc <- tryCatch(
              estimate_lcmnl(data$database, Cc,
                             start_betas = s_cc$betas,
                             start_shares = s_cc$shares,
                             dgp = dgp),
              error = function(e) NULL
            )
            if (!is.null(m_cc) && m_cc$converged) bics[Cc] <- m_cc$BIC
          }
        }
        if (which.min(bics) == tK) m_correct[mi] <- 1L
      }
    }

    best_idx <- which.max(m_ll)
    best_m <- if (length(best_idx) > 0) method_names[best_idx] else NA_character_
    list(true_K = tK, kappa = kap, sigma = sig,
         m_correct = m_correct, m_conv = m_conv, m_ll = m_ll,
         best_method = best_m)
  }

  res_list <- mclapply(1:nc, .run_one_clust, mc.cores = N_CORES_LCMNL)

  # Reassemble matrices
  method_correct <- t(sapply(res_list, `[[`, "m_correct"))
  method_converged <- t(sapply(res_list, `[[`, "m_conv"))
  cond_K   <- sapply(res_list, `[[`, "true_K")
  cond_kap <- sapply(res_list, `[[`, "kappa")
  cond_sig <- sapply(res_list, `[[`, "sigma")
  best_methods <- sapply(res_list, `[[`, "best_method")

  if (verbose) {
    for (i in 1:nc)
      cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... best=%s\n",
                  i, nc, cond_K[i], cond_kap[i], cond_sig[i], best_methods[i]))
    cat("\n  Per-method BIC accuracy (%):\n")
    for (nm in method_names) {
      valid <- method_converged[, nm] == 1
      if (any(valid)) {
        cat(sprintf("    %-15s: %.1f%% (converged %d/%d)\n", nm,
                    100 * mean(method_correct[valid, nm]),
                    sum(valid), nc))
      }
    }
    cat(sprintf("\n  Best method distribution: %s\n",
                paste(names(table(best_methods)), table(best_methods),
                      sep = "=", collapse = ", ")))
  }

  data.frame(
    true_K = cond_K, kappa = cond_kap, sigma = cond_sig,
    best_method = best_methods,
    method_correct, method_converged
  )
}

# --- SUPP 8: Sample Size and Panel Length Sensitivity ---
run_sample_sensitivity <- function(verbose = TRUE, dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 8: SAMPLE SIZE / PANEL LENGTH ====\n")

  T_values   <- c(8L, 12L, 20L)
  Npc_values <- c(50L, 100L, 150L)
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0),
                       sigma = 0.20, rep = 1:3)

  # Flatten all T x Npc x cond combos for parallel dispatch
  all_jobs <- expand.grid(Tv_idx = seq_along(T_values),
                          Npc_idx = seq_along(Npc_values),
                          cond_i = 1:nrow(conds))

  .run_one_sens <- function(row) {
    Tv  <- T_values[all_jobs$Tv_idx[row]]
    Npc <- Npc_values[all_jobs$Npc_idx[row]]
    i   <- all_jobs$cond_i[row]
    tK <- conds$true_K[i]; kap <- conds$kappa[i]; rp <- conds$rep[i]
    seed <- as.integer(7000 + 100 * tK + 10 * (kap * 100) + rp + Tv * 100 + Npc)
    data <- klue_simulate(N_per_class = Npc, T_tasks = Tv, true_K = tK,
                          separation = kap, heterogeneity = 0.2, seed = seed,
                          dgp = dgp)
    bics <- rep(Inf, 5)
    for (Cc in 1:5) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged) bics[Cc] <- m$BIC
    }
    list(Tv = Tv, Npc = Npc, correct = as.integer(which.min(bics) == tK))
  }

  res_list <- mclapply(1:nrow(all_jobs), .run_one_sens, mc.cores = N_CORES_LCMNL)

  T_col <- integer(0); Npc_col <- integer(0); accuracies <- numeric(0)
  for (Tv in T_values) {
    for (Npc in Npc_values) {
      hits <- sapply(res_list[sapply(res_list, function(x) x$Tv == Tv & x$Npc == Npc)],
                     `[[`, "correct")
      acc <- 100 * mean(hits)
      T_col <- c(T_col, Tv); Npc_col <- c(Npc_col, Npc)
      accuracies <- c(accuracies, acc)
      if (verbose) cat(sprintf("  T=%2d  N/class=%3d  (N_total=%d-%d): %.1f%%\n",
                               Tv, Npc, Npc * 2, Npc * 3, acc))
    }
  }
  data.frame(T_tasks = T_col, N_per_class = Npc_col, accuracy = accuracies)
}

# --- SUPP 9: Correlated MMNL Robustness ---
# Cholesky parameterization: b = mu + L*draws (full correlation structure)

# Helper: dynamically build apollo_randCoeff for Cholesky-correlated MMNL
.make_apollo_randCoeff_corr <- function(dgp = DGP_DEFAULT) {
  n_generic <- dgp$n_generic; n_beta <- dgp$n_beta
  # Attribute names for draws: draws_x1..draws_x{n_generic}, draws_price
  attr_short <- c(paste0("x", 1:n_generic), "price")
  draw_names <- paste0("draws_", attr_short)
  # Cholesky element naming: for generic attrs, s_x{a}_{col}; for price, s_pr_{col}
  chol_prefix <- c(paste0("s_x", 1:n_generic), "s_pr")

  lines <- c("function(apollo_beta, apollo_inputs) {", "  randcoeff <- list()")
  for (a in 1:n_beta) {
    # b_a = mu_a + sum_{col=1}^{a} L_{a,col} * draws_{col}
    chol_terms <- paste0(chol_prefix[a], "_", 1:a, " * ", draw_names[1:a])
    inner <- paste(chol_terms, collapse = " + ")
    if (a < n_beta) {
      # Generic attribute: normal distribution
      lines <- c(lines, sprintf('  randcoeff[["b_%s"]] <- mu_%s + %s',
                                 attr_short[a], attr_short[a], inner))
    } else {
      # Price: negative log-normal
      lines <- c(lines, sprintf('  randcoeff[["b_price"]] <- -exp(mu_price + %s)', inner))
    }
  }
  lines <- c(lines, '  return(randcoeff)', '}')
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- asNamespace("apollo")
  fn
}

.run_apollo_mmnl_corr <- function(database, n_draws, start_beta,
                                  dgp                = DGP_DEFAULT,
                                  n_cores            = NULL,
                                  draws_type         = DRAWS_TYPE_MMNL,
                                  estimation_routine = ESTIMATION_ROUTINE_MMNL,
                                  bounds             = NULL) {
  cleanup_apollo()
  n_generic <- dgp$n_generic

  if (is.null(n_cores)) {
    n_cores <- getOption("klue.mmnl.n_cores", .klue_default_mmnl_cores())
  }
  n_cores <- max(1L, as.integer(n_cores))

  apollo_control <<- list(
    modelName       = paste0("MMNL_corr_", as.integer(Sys.time()) %% 100000,
                             "_", sample.int(10000, 1)),
    modelDescr      = "Correlated MMNL simulation",
    indivID         = "ID",
    nCores          = n_cores,
    mixing          = TRUE,
    outputDirectory = tempdir()
  )

  apollo_beta  <<- start_beta
  apollo_fixed <<- c()

  draw_names <- c(paste0("draws_x", 1:n_generic), "draws_price")
  apollo_draws <<- list(
    interDrawsType = draws_type,
    interNDraws    = as.integer(n_draws),
    interUnifDraws = c(),
    interNormDraws = draw_names,
    intraDrawsType = draws_type,
    intraNDraws    = 0,
    intraUnifDraws = c(),
    intraNormDraws = c()
  )

  apollo_randCoeff <<- .make_apollo_randCoeff_corr(dgp)
  apollo_probabilities <<- .make_apollo_prob_mmnl(dgp)

  if (exists("apollo_lcPars", envir = .GlobalEnv)) rm("apollo_lcPars", envir = .GlobalEnv)

  # See note in .run_apollo_mmnl: apollo_probabilities must be set before
  # apollo_validateInputs.
  apollo_inputs <<- tryCatch(
    apollo_validateInputs(
      apollo_beta    = apollo_beta,
      apollo_fixed   = apollo_fixed,
      database       = database,
      apollo_control = apollo_control
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl_corr: apollo_validateInputs error: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(apollo_inputs)) return(NULL)

  est_settings <- list(
    estimationRoutine = estimation_routine,
    writeIter         = FALSE,
    silent            = TRUE
  )
  if (!is.null(bounds)) {
    est_settings$bounds <- bounds
  }

  tryCatch(
    apollo_estimate(
      apollo_beta          = apollo_beta,
      apollo_fixed         = apollo_fixed,
      apollo_probabilities = apollo_probabilities,
      apollo_inputs        = apollo_inputs,
      estimate_settings    = est_settings
    ),
    error = function(e) {
      message("klue:.run_apollo_mmnl_corr: apollo_estimate error: ",
              conditionMessage(e))
      NULL
    }
  )
}

klue_mmnl_corr <- function(database,
                               n_draws            = NULL,
                               n_draws_stage1     = NULL,
                               draws_type         = NULL,
                               estimation_routine = NULL,
                               n_cores            = NULL,
                               quiet              = NULL,
                               dgp                = DGP_DEFAULT) {
  d <- klue_mmnl_defaults()
  if (is.null(n_draws))            n_draws            <- d$n_draws
  if (is.null(n_draws_stage1))     n_draws_stage1     <- d$n_draws_stage1
  if (is.null(draws_type))         draws_type         <- d$draws_type
  if (is.null(estimation_routine)) estimation_routine <- d$estimation_routine
  if (is.null(n_cores))            n_cores            <- d$n_cores
  if (is.null(quiet))              quiet              <- d$quiet

  cleanup_apollo()

  log_file <- tempfile(pattern = "klue_mmnl_corr_", fileext = ".log")
  if (isTRUE(quiet)) {
    old_sink_out <- sink.number()
    old_sink_msg <- sink.number(type = "message")
    msg_con <- file(log_file, open = "at")
    sink(log_file)
    sink(msg_con, type = "message")
    on.exit({
      while (sink.number() > old_sink_out) sink()
      while (sink.number(type = "message") > old_sink_msg) sink(type = "message")
      try(close(msg_con), silent = TRUE)
      cleanup_apollo()
    }, add = TRUE)
  } else {
    on.exit(cleanup_apollo(), add = TRUE)
  }

  N <- length(unique(database$ID))
  n_beta <- dgp$n_beta; n_generic <- dgp$n_generic; J <- dgp$n_alternatives

  read_log_tail <- function() {
    if (!isTRUE(quiet) || !file.exists(log_file)) return(NULL)
    out <- tryCatch(readLines(log_file, warn = FALSE), error = function(e) NULL)
    if (length(out) == 0) NULL else utils::tail(out, 40L)
  }
  fail_with <- function(reason) {
    list(converged = FALSE, LL = -Inf, BIC = Inf, AIC = Inf, k = 0,
         reason          = reason,
         apollo_log_tail = read_log_tail(),
         apollo_log_path = if (isTRUE(quiet) && file.exists(log_file)) log_file else NULL,
         settings        = list(n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
                                draws_type = draws_type,
                                estimation_routine = estimation_routine,
                                n_cores = n_cores))
  }

  indep_fit <- klue_mmnl(database,
                             n_draws            = n_draws_stage1,
                             n_draws_stage1     = 100L,
                             draws_type         = draws_type,
                             estimation_routine = estimation_routine,
                             n_cores            = n_cores,
                             quiet              = quiet,
                             dgp                = dgp)

  attr_short  <- c(paste0("x", 1:n_generic), "price")
  chol_prefix <- c(paste0("s_x", 1:n_generic), "s_pr")

  if (indep_fit$converged) {
    mu_starts  <- unname(indep_fit$mu)
    sig_starts <- log(unname(indep_fit$sigma))
  } else {
    mu_starts  <- c(rep(0.5, n_generic), 0.0)
    sig_starts <- rep(log(0.5), n_beta)
  }

  # Build starting beta: ASCs + means + Cholesky lower-triangular elements
  beta0 <- c()
  for (j in 1:(J - 1)) beta0[paste0("asc_alt", j)] <- 0
  for (a in 1:n_beta) beta0[paste0("mu_", attr_short[a])] <- mu_starts[a]
  for (row in 1:n_beta) {
    for (col in 1:row) {
      name <- paste0(chol_prefix[row], "_", col)
      beta0[name] <- if (row == col) sig_starts[row] else 0
    }
  }

  stage1 <- .run_apollo_mmnl_corr(database, n_draws_stage1, beta0,
                                  dgp                = dgp,
                                  n_cores            = n_cores,
                                  draws_type         = draws_type,
                                  estimation_routine = estimation_routine)
  beta1 <- beta0
  if (!is.null(stage1) && !is.null(stage1$estimate) &&
      all(is.finite(stage1$estimate))) {
    beta1 <- stage1$estimate
  }

  model <- .run_apollo_mmnl_corr(database, n_draws, beta1,
                                 dgp                = dgp,
                                 n_cores            = n_cores,
                                 draws_type         = draws_type,
                                 estimation_routine = estimation_routine)
  if (is.null(model))                                 return(fail_with("stage2_apollo_estimate_failed"))
  if (is.null(model$estimate) ||
      !all(is.finite(model$estimate)))                return(fail_with("stage2_non_finite_estimates"))
  if (is.null(model$LLout) || !is.finite(model$LLout[1]))
                                                      return(fail_with("stage2_non_finite_LL"))

  est    <- model$estimate
  LL     <- model$LLout[1]
  n_free <- length(est)
  BIC    <- -2 * LL + n_free * log(N)
  AIC    <- -2 * LL + 2 * n_free

  list(converged = TRUE, LL = LL, BIC = BIC, AIC = AIC, k = n_free,
       reason = "ok", apollo_log_tail = NULL, apollo_log_path = NULL,
       settings = list(n_draws = n_draws, n_draws_stage1 = n_draws_stage1,
                       draws_type = draws_type,
                       estimation_routine = estimation_routine,
                       n_cores = n_cores))
}

run_correlated_mmnl_robustness <- function(n_draws = N_DRAWS_MMNL, verbose = TRUE,
                                            dgp = DGP_DEFAULT) {
  if (verbose) cat("\n==== SUPP 9: CORRELATED MMNL ROBUSTNESS ====\n")

  # Subset: K=1 (continuous DGP) + K=2,3 (discrete DGP), moderate conditions
  conds_k1 <- expand.grid(true_K = 1, kappa = 0,
                           sigma = c(0.15, 0.25, 0.35), rep = 1:2)
  conds_kn <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0),
                           sigma = c(0.15, 0.25), rep = 1:2)
  conds <- rbind(conds_k1, conds_kn)
  nc <- nrow(conds)
  if (verbose) cat(sprintf("  %d conditions (K=1: %d, K>=2: %d)\n",
                           nc, nrow(conds_k1), nrow(conds_kn)))

  # Phase 1: LCMNL in parallel (no Apollo)
  if (verbose) cat("  Phase 1: LCMNL estimation (parallel)...\n")
  .run_lcmnl_corr_part <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- as.integer(1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)
    npc <- if (tK == 1L) 300L else 150L
    data <- klue_simulate(N_per_class = npc, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    best_lc_bic <- Inf; best_lc_C <- NA_integer_
    for (Cc in 1:5) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged && m$BIC < best_lc_bic) {
        best_lc_bic <- m$BIC; best_lc_C <- Cc
      }
    }
    list(tK = tK, kap = kap, sig = sig, seed = seed, npc = npc,
         lc_bic = best_lc_bic, lc_C = best_lc_C)
  }
  lc_results <- mclapply(1:nc, .run_lcmnl_corr_part, mc.cores = N_CORES_LCMNL)

  # Phase 2: MMNL (independent + correlated) sequentially via Apollo
  if (verbose) cat("  Phase 2: MMNL estimation (sequential, Apollo)...\n")
  rK <- integer(nc); rkap <- numeric(nc); rsig <- numeric(nc)
  r_lc <- rep(NA, nc); r_indep <- rep(NA, nc); r_corr <- rep(NA, nc)
  r_lc_C <- rep(NA_integer_, nc); valid <- logical(nc)

  for (i in 1:nc) {
    lr <- lc_results[[i]]
    rK[i] <- lr$tK; rkap[i] <- lr$kap; rsig[i] <- lr$sig

    if (verbose) cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... ", i, nc, lr$tK, lr$kap, lr$sig))

    data <- klue_simulate(N_per_class = lr$npc, T_tasks = 20, true_K = lr$tK,
                          separation = lr$kap, heterogeneity = lr$sig, seed = lr$seed,
                          dgp = dgp)

    mm_indep <- klue_mmnl(data$database, n_draws = n_draws, dgp = dgp)
    mm_corr  <- klue_mmnl_corr(data$database, n_draws = n_draws, dgp = dgp)

    if (is.finite(lr$lc_bic) && mm_indep$converged && mm_corr$converged) {
      r_lc[i] <- lr$lc_bic; r_lc_C[i] <- lr$lc_C
      r_indep[i] <- mm_indep$BIC; r_corr[i] <- mm_corr$BIC
      valid[i] <- TRUE

      bics <- c(LCMNL = lr$lc_bic, MMNL_indep = mm_indep$BIC, MMNL_corr = mm_corr$BIC)
      winner <- names(which.min(bics))
      lc_lab <- if (lr$lc_C == 1L) "MNL(C=1)" else sprintf("LCMNL(C=%d)", lr$lc_C)
      if (verbose) cat(sprintf("%s vs indep vs corr -> %s\n", lc_lab, winner))
    } else {
      if (verbose) cat("FAILED\n")
    }
  }

  idx <- which(valid)
  bic_winner <- mapply(function(lc, ind, cor) {
    bics <- c(LCMNL = lc, MMNL_indep = ind, MMNL_corr = cor)
    names(which.min(bics))
  }, r_lc[idx], r_indep[idx], r_corr[idx])

  df <- data.frame(true_K = rK[idx], kappa = rkap[idx], sigma = rsig[idx],
                   lcmnl_BIC = r_lc[idx], lcmnl_C = r_lc_C[idx],
                   mmnl_indep_BIC = r_indep[idx], mmnl_corr_BIC = r_corr[idx],
                   bic_prefers = bic_winner, stringsAsFactors = FALSE)

  if (verbose && nrow(df) > 0) {
    n_tot <- nrow(df)
    cat(sprintf("  LCMNL: %d/%d  MMNL_indep: %d/%d  MMNL_corr: %d/%d\n",
                sum(df$bic_prefers == "LCMNL"), n_tot,
                sum(df$bic_prefers == "MMNL_indep"), n_tot,
                sum(df$bic_prefers == "MMNL_corr"), n_tot))
    kn <- df[df$true_K > 1, ]
    if (nrow(kn) > 0)
      cat(sprintf("  K>1 (discrete DGP): LCMNL %d/%d  MMNL_corr %d/%d\n",
                  sum(kn$bic_prefers == "LCMNL"), nrow(kn),
                  sum(kn$bic_prefers == "MMNL_corr"), nrow(kn)))
  }
  df
}

# =============================================================================
# SECTION 9: MASTER EXECUTION
# =============================================================================

run_full_study <- function(run_main = TRUE, run_mmnl = TRUE, run_supp = TRUE,
                           verbose = TRUE, dgp = DGP_DEFAULT) {
  results <- list()
  t_start <- Sys.time()

  if (run_main) {
    results[["main"]] <- run_main_simulation(verbose = verbose, dgp = dgp)
    summarise_main_results(results[["main"]])
    write.csv(results[["main"]], file.path(OUTPUT_DIR, "main_results.csv"),
              row.names = FALSE)
  }

  if (run_mmnl) {
    results[["mmnl"]] <- run_mmnl_comparison(verbose = verbose, dgp = dgp)
    write.csv(results[["mmnl"]], file.path(OUTPUT_DIR, "mmnl_results.csv"),
              row.names = FALSE)
  }

  if (run_supp) {
    results[["convergence"]] <- run_convergence_ablation(verbose = verbose, dgp = dgp)
    write.csv(results[["convergence"]],
              file.path(OUTPUT_DIR, "convergence_results.csv"), row.names = FALSE)

    results[["unbalanced"]] <- run_unbalanced_analysis(verbose = verbose, dgp = dgp)
    write.csv(results[["unbalanced"]],
              file.path(OUTPUT_DIR, "unbalanced_results.csv"), row.names = FALSE)

    results[["design"]] <- run_design_comparison(verbose = verbose, dgp = dgp)
    write.csv(results[["design"]],
              file.path(OUTPUT_DIR, "design_results.csv"), row.names = FALSE)

    results[["concomitant"]] <- run_concomitant_analysis(verbose = verbose, dgp = dgp)
    write.csv(results[["concomitant"]],
              file.path(OUTPUT_DIR, "concomitant_results.csv"), row.names = FALSE)

    results[["recovery"]] <- run_unconditional_recovery(verbose = verbose, dgp = dgp)
    write.csv(results[["recovery"]],
              file.path(OUTPUT_DIR, "recovery_results.csv"), row.names = FALSE)

    results[["clustering"]] <- run_clustering_comparison(verbose = verbose, dgp = dgp)
    write.csv(results[["clustering"]],
              file.path(OUTPUT_DIR, "clustering_comparison_results.csv"), row.names = FALSE)

    results[["sample_sens"]] <- run_sample_sensitivity(verbose = verbose, dgp = dgp)
    write.csv(results[["sample_sens"]],
              file.path(OUTPUT_DIR, "sample_sensitivity_results.csv"), row.names = FALSE)

    results[["mmnl_corr"]] <- run_correlated_mmnl_robustness(verbose = verbose, dgp = dgp)
    write.csv(results[["mmnl_corr"]],
              file.path(OUTPUT_DIR, "mmnl_correlated_results.csv"), row.names = FALSE)
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  if (verbose) cat(sprintf("\n==== DONE: %.1f min ====\n", elapsed))

  return(results)
}

# Package version: simulation auto-run block intentionally removed.
# Users can call run_full_study() / run_main_simulation() directly if needed.

# =============================================================================
# Backward-compatibility aliases.
# Old names map to the new canonical klue_* names; both are exported so
# existing call sites continue to work unchanged. Slated for removal in a
# future major release.
# =============================================================================

make_dgp_config           <- klue_dgp
estimate_lcmnl_multistart <- klue_lcmnl
estimate_mmnl             <- klue_mmnl
estimate_mmnl_corr        <- klue_mmnl_corr

# Data-generation aliases
generate_data                 <- klue_simulate
generate_data_with_covariates <- klue_simulate_cov
generate_data_defficient      <- klue_simulate_deff
