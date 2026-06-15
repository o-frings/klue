# Starting values: respondent features, the six clustering methods behind one
# registry, cluster-wise MNL fits, and the two non-clustering baselines
# (pooled-MNL perturbation, random partition) used by the estimator study.

# Per-respondent revealed-preference signatures: mean chosen-vs-unchosen
# attribute contrasts, averaged over the (balanced) panel.
compute_rp_features <- function(database, dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  n_obs <- nrow(database)
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  CH <- database$CHOICE
  ri <- 1:n_obs
  diffs <- matrix(0, nrow = n_obs, ncol = n_beta)
  for (a in 1:n_beta) {
    aname <- if (a <= n_generic) paste0("x", a) else "price"
    Xa <- matrix(0, nrow = n_obs, ncol = J)
    for (j in 1:J) Xa[, j] <- database[[paste0(aname, "_", j)]]
    rsa <- rowSums(Xa)
    cha <- Xa[cbind(ri, CH)]
    diffs[, a] <- cha - (rsa - cha) / (J - 1)
  }
  T_per_n <- as.integer(n_obs / N)
  features <- matrix(0, nrow = N, ncol = n_beta)
  for (k in 1:n_beta) {
    features[, k] <- colSums(matrix(diffs[, k], T_per_n, N)) / T_per_n
  }
  features
}

# One-hot choice indicators: N x (T*J), 1 where respondent i chose alt j at
# task t. The ablation alternative to RP contrasts.
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
    features[i, ((seq_len(T_per_n) - 1L) * J) + chosen] <- 1
  }
  features
}

standardise_features <- function(features) {
  mu <- colMeans(features); sigma <- apply(features, 2, sd)
  sigma[sigma == 0] <- 1
  scaled <- scale(features, center = mu, scale = sigma)
  scaled[!is.finite(scaled)] <- 0
  list(scaled = scaled, mu = mu, sigma = sigma)
}

# Fit a separate MNL (C=1) per cluster: starting values directly in
# coefficient space rather than heuristically scaled centroids.
fit_cluster_mnls <- function(labels, database, dgp = DGP_DEFAULT) {
  C <- max(labels)
  all_ids <- unique(database$ID)
  n_beta <- dgp$n_beta
  fallback <- c(rep(0.5, dgp$n_generic), -0.5)

  betas <- matrix(0, nrow = C, ncol = n_beta)
  shares <- as.numeric(table(factor(labels, levels = 1:C))) / length(all_ids)

  for (cc in 1:C) {
    cluster_ids <- all_ids[labels == cc]
    if (length(cluster_ids) < 3) { betas[cc, ] <- fallback; next }
    db_sub <- database[database$ID %in% cluster_ids, , drop = FALSE]
    db_sub$ID <- match(db_sub$ID, cluster_ids)
    mnl_fit <- tryCatch(
      estimate_lcmnl(db_sub, C = 1,
                     start_betas = matrix(0, 1, n_beta), dgp = dgp),
      error = function(e) NULL
    )
    betas[cc, ] <- if (!is.null(mnl_fit) && mnl_fit$converged)
      mnl_fit$betas[1, ] else fallback
  }
  list(betas = betas, shares = shares)
}

KLUE_CLUSTER_METHODS <- c("kmeans", "gmm", "hc_ward", "hc_complete",
                          "hc_average", "pam")

# Cluster labels for one method. The stochastic methods are seeded exactly as
# in 0.6.x (set.seed(123) immediately before the call).
.cluster_labels <- function(method, scaled, C) {
  switch(method,
    kmeans = {
      set.seed(123)
      kmeans(scaled, centers = C, nstart = 25, iter.max = 100)$cluster
    },
    gmm = {
      set.seed(123)
      g <- mclust::Mclust(scaled, G = C, verbose = FALSE)
      if (is.null(g)) NULL else g$classification
    },
    hc_ward     = cutree(hclust(dist(scaled), method = "ward.D2"), k = C),
    hc_complete = cutree(hclust(dist(scaled), method = "complete"), k = C),
    hc_average  = cutree(hclust(dist(scaled), method = "average"), k = C),
    pam = {
      set.seed(123)
      cluster::pam(scaled, k = C)$clustering
    },
    stop("Unknown clustering method: ", method)
  )
}

#' Clustering-based starting values for one method
#'
#' Clusters respondents by their choice features, then fits a class-wise MNL to
#' produce starting coefficients and class shares for the latent-class
#' estimator.
#'
#' @param database choice data frame with respondent ID, CHOICE, and attribute
#'   columns, as consumed by the estimator.
#' @param C number of latent classes (number of clusters to form). With C = 1
#'   the data are treated as a single class and no clustering is performed.
#' @param method clustering method, one of "kmeans", "gmm", "hc_ward",
#'   "hc_complete", "hc_average", "pam".
#' @param features optional pre-computed N-by-p feature matrix. When NULL, the
#'   features are computed from `database` according to `feature_type`.
#' @param feature_type "rp" (revealed-preference contrasts, the default) or
#'   "onehot" (choice indicators; the ablation arm).
#' @param dgp data-generating-process specification supplying dimensions such
#'   as the number of alternatives, coefficients, and generic coefficients.
#' @return A list with `betas` (a C-by-n_beta matrix of class-wise starting
#'   coefficients) and `shares` (a length-C vector of class shares), or NULL if
#'   the chosen clustering method failed to return labels.
#' @export
klue_starts <- function(database, C, method = "kmeans", features = NULL,
                        feature_type = c("rp", "onehot"), dgp = DGP_DEFAULT) {
  feature_type <- match.arg(feature_type)
  if (is.null(features)) {
    features <- if (feature_type == "rp") compute_rp_features(database, dgp)
                else compute_onehot_features(database, dgp)
  }
  if (C == 1) return(fit_cluster_mnls(rep(1L, nrow(features)), database, dgp))
  sf <- standardise_features(features)
  labels <- .cluster_labels(method, sf$scaled, C)
  if (is.null(labels)) return(NULL)
  fit_cluster_mnls(labels, database, dgp)
}

# All six methods, computing the feature matrix once. NULL entries mark
# methods that errored (the multistart skips them).
get_all_starts <- function(database, C, dgp = DGP_DEFAULT,
                           feature_type = c("rp", "onehot")) {
  feature_type <- match.arg(feature_type)
  features <- if (feature_type == "rp") compute_rp_features(database, dgp)
              else compute_onehot_features(database, dgp)
  sapply(KLUE_CLUSTER_METHODS, function(m) {
    tryCatch(klue_starts(database, C, m, features = features, dgp = dgp),
             error = function(e) NULL)
  }, simplify = FALSE)
}

# MNL-perturbation starts: pooled MNL + jitter. The informed-but-not-
# clustering baseline (Apollo/gmnl/lclogit warm-start practice). Generic
# coefficients get additive Gaussian jitter; price gets multiplicative
# log-normal jitter so it stays negative. Deterministic given `seed`.
get_mnl_perturbation_starts <- function(database, C, n_starts = 50L,
                                        jitter_sd = 0.5, seed = 1L,
                                        dgp = DGP_DEFAULT) {
  n_generic <- dgp$n_generic; n_beta <- dgp$n_beta
  pooled <- estimate_lcmnl(database, 1L, dgp = dgp)
  beta_hat <- if (pooled$converged) as.numeric(pooled$betas[1, ]) else dgp$beta_bar
  b_gen   <- beta_hat[1:n_generic]
  b_price <- beta_hat[n_beta]
  if (!is.finite(b_price) || b_price >= 0) b_price <- -abs(dgp$beta_bar[n_beta])

  lapply(seq_len(n_starts), function(s) {
    set.seed(seed * 100000L + s)
    bs <- matrix(0, C, n_beta)
    for (ci in 1:C) {
      bs[ci, 1:n_generic] <- b_gen + rnorm(n_generic, 0, jitter_sd)
      bs[ci, n_beta] <- b_price * exp(rnorm(1, 0, jitter_sd))
    }
    list(betas = bs, shares = rep(1 / C, C))
  })
}

# Random-partition starts: what Stata lclogit/Latent GOLD do by default --
# random assignment of respondents to classes, then class-wise MNL fits.
# Deterministic given `seed`; RNG stream distinct from the perturbation arm.
get_random_partition_starts <- function(database, C, n_starts = 50L,
                                        seed = 1L, dgp = DGP_DEFAULT) {
  N <- length(unique(database$ID))
  lapply(seq_len(n_starts), function(s) {
    set.seed(seed * 100000L + 50000L + s)
    labels <- sample.int(C, N, replace = TRUE)
    tries <- 0L
    while (any(tabulate(labels, C) < 3L) && tries < 20L) {
      labels <- sample.int(C, N, replace = TRUE); tries <- tries + 1L
    }
    st <- fit_cluster_mnls(labels, database, dgp = dgp)
    list(betas = st$betas, shares = st$shares)
  })
}
