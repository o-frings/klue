# Data generation: blocked D-efficient design + one simulator covering the
# plain, concomitant-covariate, and legacy fake-D-efficient DGPs.
#
# RNG discipline: every code path consumes the random number stream in exactly
# the same order as klue 0.6.x, so seeded datasets are bit-identical.

# Class deviation matrix using equally-rotated cosines: deviations at equal
# angular spacing on all n_beta attributes (including price), equal norm
# sqrt(n_beta/2). Avoids confounding cluster type with number of classes.
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

# Gaussian-copula attribute draws in [-1, 1]: Z ~ MVN(0, R) -> U = pnorm(Z) ->
# X = 2U - 1. Marginals stay symmetric (mean 0); only cross-attribute
# dependence changes. R must have unit diagonal.
.draw_correlated_attrs <- function(n_obs, n_generic, R) {
  L <- chol(R)
  Z <- matrix(rnorm(n_obs * n_generic), n_obs, n_generic) %*% L
  2 * pnorm(Z) - 1
}

# attr_corr may be NULL (independent), a scalar rho (equicorrelation), or a
# full correlation matrix.
.attr_corr_matrix <- function(attr_corr, n_generic) {
  if (is.null(attr_corr)) return(NULL)
  if (length(attr_corr) == 1L) {
    R <- matrix(attr_corr, n_generic, n_generic); diag(R) <- 1
    return(R)
  }
  R <- as.matrix(attr_corr)
  stopifnot(nrow(R) == n_generic, ncol(R) == n_generic)
  R
}

#' Blocked D-efficient design from a population-mean prior
#'
#' One locally-D-efficient design is constructed for the whole study (the
#' analyst designs around population means, NOT the unknown latent classes),
#' then split into blocks. Each respondent answers one block, so all
#' respondents in a block face identical cards. This is the realistic-DCE
#' baseline. Requires the idefix package.
#'
#' @param n_cards Total number of choice cards (sets) in the design. Must be a
#'   multiple of \code{n_blocks}.
#' @param n_blocks Number of blocks the cards are split into. Each respondent
#'   answers one block.
#' @param dgp Data-generating-process specification supplying
#'   \code{n_alternatives}, \code{n_beta}, \code{n_generic}, and
#'   \code{beta_bar}.
#' @param priors Numeric prior coefficient vector used to construct the
#'   locally-D-efficient design. Defaults to \code{dgp$beta_bar} when
#'   \code{NULL}.
#' @param n_lvls Number of attribute levels per attribute.
#' @param n_start Number of random starting designs for the Modfed search.
#' @param seed Integer seed for reproducible design construction and blocking.
#' @return A list with elements \code{cards} (an \code{n_cards} by
#'   \code{n_alternatives} by \code{n_beta} array of attribute values),
#'   \code{blocks} (a list mapping each block to its card indices), \code{T}
#'   (cards per block), \code{n_cards}, \code{n_blocks}, \code{dgp}, and
#'   \code{Derror} (the Bayesian D-error of the chosen design).
#' @export
klue_design <- function(n_cards = 48L, n_blocks = 4L,
                        dgp = DGP_DEFAULT, priors = NULL,
                        n_lvls = 4L, n_start = 8L, seed = 20240601L) {
  stopifnot(n_cards %% n_blocks == 0L)
  if (!requireNamespace("idefix", quietly = TRUE))
    stop("Package 'idefix' is required for the blocked D-efficient design (klue_design).")
  if (is.null(priors)) priors <- dgp$beta_bar
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  set.seed(seed)
  lvl_g  <- seq(-1, 1, length.out = n_lvls)
  lvl_p  <- seq(0.1, 0.9, length.out = n_lvls)
  c.lvls <- c(rep(list(lvl_g), n_generic), list(lvl_p))
  cs <- idefix::Profiles(lvls = rep(n_lvls, n_beta),
                         coding = rep("C", n_beta), c.lvls = c.lvls)
  D  <- idefix::Modfed(cand.set = cs, n.sets = n_cards, n.alts = J,
                       alt.cte = rep(0, J),
                       par.draws = matrix(priors, nrow = 1), n.start = n_start)
  des <- D$BestDesign$design          # (n_cards*J) x n_beta, rows stacked set.alt
  cards <- array(0, dim = c(n_cards, J, n_beta))
  for (s in 1:n_cards) for (j in 1:J) cards[s, j, ] <- des[(s - 1) * J + j, ]
  perm   <- sample(n_cards)
  blocks <- split(perm, rep(1:n_blocks, each = n_cards %/% n_blocks))
  list(cards = cards, blocks = blocks, T = n_cards %/% n_blocks,
       n_cards = n_cards, n_blocks = n_blocks, dgp = dgp,
       Derror = D$BestDesign$DB.error)
}

# Individual betas around the class means; price truncated at -0.1.
.draw_individual_betas <- function(true_class, true_betas, heterogeneity, dgp) {
  N <- length(true_class); n_beta <- dgp$n_beta
  sigma_vec <- c(rep(heterogeneity, dgp$n_generic), heterogeneity * 0.3)
  ib <- matrix(0, nrow = N, ncol = n_beta)
  for (n in 1:N) {
    ib[n, ] <- rnorm(n_beta, true_betas[true_class[n], ], sigma_vec)
    ib[n, n_beta] <- min(ib[n, n_beta], -0.1)
  }
  ib
}

# Attribute columns: blocked-design path (respondents share their block's
# cards; consumes one sample(N) call) or random path (fresh draws per row,
# optionally copula-correlated). Returns the database skeleton and T_tasks
# (which the design overrides with its block size).
.build_attr_database <- function(N, T_tasks, dgp, design = NULL, attr_corr = NULL) {
  J <- dgp$n_alternatives; n_generic <- dgp$n_generic; n_beta <- dgp$n_beta
  if (!is.null(design)) {
    T_tasks <- design$T
    n_obs   <- N * T_tasks
    resp_block <- rep(1:design$n_blocks, length.out = N)[sample(N)]  # balanced
    database <- data.frame(ID = rep(1:N, each = T_tasks),
                           TASK = rep(1:T_tasks, times = N))
    for (j in 1:J) {
      for (a in 1:n_generic) database[[paste0("x", a, "_", j)]] <- numeric(n_obs)
      database[[paste0("price_", j)]] <- numeric(n_obs)
    }
    for (n in 1:N) {
      cardset <- design$blocks[[resp_block[n]]]
      rows    <- ((n - 1) * T_tasks + 1):(n * T_tasks)
      for (j in 1:J) {
        for (a in 1:n_generic)
          database[[paste0("x", a, "_", j)]][rows] <- design$cards[cardset, j, a]
        database[[paste0("price_", j)]][rows] <- design$cards[cardset, j, n_beta]
      }
    }
  } else {
    n_obs <- N * T_tasks
    database <- data.frame(ID = rep(1:N, each = T_tasks),
                           TASK = rep(1:T_tasks, times = N))
    R <- .attr_corr_matrix(attr_corr, n_generic)
    for (j in 1:J) {
      if (is.null(R)) {
        for (a in 1:n_generic) database[[paste0("x", a, "_", j)]] <- runif(n_obs, -1, 1)
      } else {
        Xc <- .draw_correlated_attrs(n_obs, n_generic, R)
        for (a in 1:n_generic) database[[paste0("x", a, "_", j)]] <- Xc[, a]
      }
      database[[paste0("price_", j)]] <- runif(n_obs, 0.1, 0.9)
    }
  }
  list(database = database, T_tasks = T_tasks)
}

# Type-1-EV choice simulation via matrix multiply; adds CHOICE in place.
.simulate_choices <- function(database, individual_betas, dgp) {
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  n_obs <- nrow(database)
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
  database
}

#' Simulate panel choice data from the Frings (2026) DGP
#'
#' One entry point for the plain DGP and the concomitant-covariate DGP
#' (\code{covariates = TRUE}; class membership driven by Z1 ~ N(0,1) and a
#' binary Z2). With \code{design} from \code{klue_design()}, respondents share
#' their block's cards; otherwise attributes are drawn fresh per row,
#' optionally correlated via \code{attr_corr}.
#'
#' @param N_per_class Number of respondents per latent class.
#' @param T_tasks Number of choice tasks per respondent. Overridden by the
#'   block size when \code{design} is supplied.
#' @param true_K Number of latent classes.
#' @param separation Scalar controlling the distance of class means from the
#'   grand mean \code{beta_bar}.
#' @param heterogeneity Within-class standard deviation of the individual
#'   coefficients.
#' @param seed Integer seed for reproducible simulation.
#' @param class_proportions Optional numeric vector of class shares summing to
#'   one. When \code{NULL}, classes are equally sized. Ignored when
#'   \code{covariates = TRUE}.
#' @param dgp Data-generating-process specification supplying \code{n_beta},
#'   \code{beta_bar}, \code{n_alternatives}, and \code{n_generic}.
#' @param sep_profile Optional numeric vector of per-attribute separation
#'   weights. When \code{NULL}, all attributes are weighted equally. Ignored
#'   when \code{covariates = TRUE}.
#' @param attr_corr Optional attribute correlation: \code{NULL} (independent),
#'   a scalar equicorrelation, or a full correlation matrix. Applies only to
#'   the random (non-design) path and is ignored when \code{covariates = TRUE}.
#' @param design Optional blocked design from \code{klue_design()}. When
#'   supplied, respondents share their block's cards.
#' @param seg_scale Optional per-class scale (recycled to length
#'   \code{true_K}) placing segments at unequal distances from the grand mean.
#'   Ignored when \code{covariates = TRUE}.
#' @param covariates Logical; when \code{TRUE}, class membership is driven by
#'   covariates Z1 ~ N(0, 1) and a binary Z2 rather than fixed proportions.
#' @param covariate_strength Scalar scaling the covariate effect on class
#'   membership when \code{covariates = TRUE}.
#' @return A list with elements \code{database} (the simulated panel with a
#'   \code{CHOICE} column), \code{true_betas}, \code{true_class},
#'   \code{individual_betas}, \code{N}, \code{T}, \code{K}, and \code{dgp};
#'   when \code{covariates = TRUE} it also contains \code{Z1} and \code{Z2}
#'   (and adds these columns to \code{database}).
#' @export
klue_simulate <- function(N_per_class = 150, T_tasks = 20, true_K = 2,
                          separation = 1.0, heterogeneity = 0.25,
                          seed = 12345, class_proportions = NULL,
                          dgp = DGP_DEFAULT, sep_profile = NULL,
                          attr_corr = NULL, design = NULL, seg_scale = NULL,
                          covariates = FALSE, covariate_strength = 1.0) {
  set.seed(seed)
  n_beta <- dgp$n_beta; beta_bar <- dgp$beta_bar
  N <- N_per_class * true_K
  Z1 <- NULL; Z2 <- NULL

  if (covariates) {
    # Covariate DGP ignores sep_profile / seg_scale / class_proportions /
    # attr_corr (as in 0.6.x); class membership comes from the covariates.
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
    attr_corr <- NULL
  } else {
    segment_dev <- generate_segment_deviations(true_K, n_beta)
    true_betas <- matrix(0, nrow = true_K, ncol = n_beta)
    sep_weights <- if (is.null(sep_profile)) rep(1, n_beta) else sep_profile
    # seg_scale != 1 places segments at unequal distances from the grand mean
    # (asymmetric geometry), breaking the equal-norm cosine structure.
    sc <- if (is.null(seg_scale)) rep(1, true_K) else rep_len(seg_scale, true_K)
    for (cc in 1:true_K) {
      true_betas[cc, ] <- beta_bar + separation * sc[cc] * sep_weights * segment_dev[cc, ]
    }
    if (is.null(class_proportions)) {
      true_class <- rep(1:true_K, each = N_per_class)
    } else {
      sizes <- round(N * class_proportions)
      sizes[length(sizes)] <- N - sum(sizes[-length(sizes)])
      true_class <- unlist(lapply(1:true_K, function(cc) rep(cc, sizes[cc])))
      N <- length(true_class)
    }
  }

  individual_betas <- .draw_individual_betas(true_class, true_betas,
                                             heterogeneity, dgp)
  ad <- .build_attr_database(N, T_tasks, dgp, design = design,
                             attr_corr = attr_corr)
  database <- .simulate_choices(ad$database, individual_betas, dgp)

  out <- list(database = database, true_betas = true_betas,
              true_class = true_class, individual_betas = individual_betas,
              N = N, T = ad$T_tasks, K = true_K, dgp = dgp)
  if (covariates) {
    out$database$Z1 <- Z1[database$ID]
    out$database$Z2 <- Z2[database$ID]
    out$Z1 <- Z1; out$Z2 <- Z2
  }
  out
}

#' @rdname klue_simulate
#' @export
klue_simulate_cov <- function(N_per_class = 150, T_tasks = 20, true_K = 2,
                              separation = 1.0, heterogeneity = 0.25,
                              seed = 12345, covariate_strength = 1.0,
                              dgp = DGP_DEFAULT, design = NULL) {
  klue_simulate(N_per_class = N_per_class, T_tasks = T_tasks, true_K = true_K,
                separation = separation, heterogeneity = heterogeneity,
                seed = seed, dgp = dgp, design = design,
                covariates = TRUE, covariate_strength = covariate_strength)
}

#' Legacy grid-with-jitter design (deprecated)
#'
#' Kept only so klue_study_design() and old scripts reproduce 0.6.x output.
#' This is NOT a D-optimal design (fixed level grid + jitter); use
#' \code{klue_design()} + the \code{design} argument of \code{klue_simulate}.
#'
#' @param N_per_class Number of respondents per latent class.
#' @param T_tasks Number of choice tasks per respondent.
#' @param true_K Number of latent classes.
#' @param separation Scalar controlling the distance of class means from the
#'   grand mean \code{beta_bar}.
#' @param heterogeneity Within-class standard deviation of the individual
#'   coefficients.
#' @param seed Integer seed for reproducible simulation.
#' @param dgp Data-generating-process specification supplying
#'   \code{n_alternatives}, \code{n_beta}, \code{n_generic}, and
#'   \code{beta_bar}.
#' @return A list with elements \code{database} (the simulated panel with a
#'   \code{CHOICE} column), \code{true_betas}, \code{true_class},
#'   \code{individual_betas}, \code{N}, \code{T}, \code{K}, and \code{dgp}.
#' @export
klue_simulate_deff <- function(N_per_class = 150, T_tasks = 20, true_K = 2,
                               separation = 1.0, heterogeneity = 0.25,
                               seed = 12345, dgp = DGP_DEFAULT) {
  set.seed(seed)
  N <- N_per_class * true_K; J <- dgp$n_alternatives; n_beta <- dgp$n_beta
  n_generic <- dgp$n_generic

  segment_dev <- generate_segment_deviations(true_K, n_beta)
  true_betas <- matrix(0, nrow = true_K, ncol = n_beta)
  for (cc in 1:true_K) true_betas[cc, ] <- dgp$beta_bar + separation * segment_dev[cc, ]
  true_class <- rep(1:true_K, each = N_per_class)

  individual_betas <- .draw_individual_betas(true_class, true_betas,
                                             heterogeneity, dgp)

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
  database <- .simulate_choices(database, individual_betas, dgp)

  list(database = database, true_betas = true_betas, true_class = true_class,
       individual_betas = individual_betas, N = N, T = T_tasks, K = true_K,
       dgp = dgp)
}
