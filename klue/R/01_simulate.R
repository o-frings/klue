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
#' @param priors Prior coefficients: a numeric vector gives a locally
#'   D-efficient design at that point; a matrix (rows = prior draws, columns =
#'   coefficients) gives a Bayesian D-efficient design with the D-error
#'   averaged over the draws. Defaults to \code{dgp$beta_bar} when \code{NULL}.
#' @param n_lvls Number of attribute levels per attribute.
#' @param n_start Number of random starting designs for the Modfed search.
#' @param seed Integer seed for reproducible design construction and blocking.
#' @return A list with elements \code{cards} (an \code{n_cards} by
#'   \code{n_alternatives} by \code{n_beta} array of attribute values),
#'   \code{blocks} (a list mapping each block to its card indices), \code{T}
#'   (cards per block), \code{n_cards}, \code{n_blocks}, \code{dgp}, and
#'   \code{Derror} (the D-error of the chosen design: Bayesian when
#'   \code{priors} is a matrix of draws, local otherwise).
#' @export
klue_design <- function(n_cards = 48L, n_blocks = 4L,
                        dgp = DGP_DEFAULT, priors = NULL,
                        n_lvls = 4L, n_start = 8L, seed = 20240601L) {
  stopifnot(n_cards %% n_blocks == 0L)
  if (!requireNamespace("idefix", quietly = TRUE))
    stop("Package 'idefix' is required for the blocked D-efficient design (klue_design).")
  if (is.null(priors)) priors <- dgp$beta_bar
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta; n_generic <- dgp$n_generic
  pd <- if (is.matrix(priors)) priors else matrix(priors, nrow = 1)
  if (ncol(pd) != n_beta)
    stop("klue_design: priors must have ", n_beta,
         " columns (n_beta); got ", ncol(pd))
  set.seed(seed)
  lvl_g  <- seq(-1, 1, length.out = n_lvls)
  lvl_p  <- seq(0.1, 0.9, length.out = n_lvls)
  c.lvls <- c(rep(list(lvl_g), n_generic), list(lvl_p))
  cs <- idefix::Profiles(lvls = rep(n_lvls, n_beta),
                         coding = rep("C", n_beta), c.lvls = c.lvls)
  D  <- idefix::Modfed(cand.set = cs, n.sets = n_cards, n.alts = J,
                       alt.cte = rep(0, J),
                       par.draws = pd, n.start = n_start,
                       parallel = FALSE)  # idefix's default spawns detectCores() - 1 workers;
                                          # the design is identical either way
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

#' Simulate choices for a canonical database from given coefficients
#'
#' Draws type-1 extreme-value choices for every (respondent, task) row of a
#' canonical wide database from the utilities \eqn{V_j = x_j'\beta}. There is
#' no separate ASC block: alternative-specific constants and framing shifts are
#' carried as indicator attributes of the database.
#'
#' @param database Canonical wide database (e.g. from
#'   \code{\link{klue_database_long}}); an existing \code{CHOICE} column is
#'   replaced.
#' @param betas Numeric vector of length \code{dgp$n_beta} (generic attributes
#'   in order, then price), or an N x \code{n_beta} matrix of respondent-specific
#'   coefficients whose row k belongs to \code{ID == k}.
#' @param dgp Data-generating-process specification giving
#'   \code{n_alternatives}, \code{n_generic} and \code{n_beta}
#'   (\code{\link{klue_dgp}}).
#' @return The database with a simulated \code{CHOICE} column.
#' @export
klue_simulate_choices <- function(database, betas, dgp) {
  N <- max(database$ID); k <- if (is.matrix(betas)) ncol(betas) else length(betas)
  ib <- if (is.matrix(betas)) betas else matrix(betas, nrow = N, ncol = k, byrow = TRUE)
  if (nrow(ib) != N || !k %in% c(dgp$n_beta, dgp$npc))
    stop(sprintf("betas must have %d columns (taste) or %d (taste and constants), one row per respondent (%d)",
                 dgp$n_beta, dgp$npc, N))
  if (k == dgp$n_beta) return(.simulate_choices(database, ib, dgp))
  V <- .row_utilities(.design_full(database, dgp), ib[database$ID, , drop = FALSE])
  U <- V - log(-log(matrix(runif(length(V)), nrow = nrow(V))))
  database$CHOICE <- max.col(U)
  database
}

#' Draw respondent coefficients from a fitted model
#'
#' Returns an \code{N} x \code{npc} matrix of coefficients (taste
#' coefficients, then constants, named by label) for
#' \code{\link{klue_simulate_choices}}: the MNL vector repeated, or for a
#' latent class fit each row the coefficients of a class drawn with the class
#' shares (common coefficients share their one value). Uses the caller's RNG;
#' the drawn classes are in \code{attr(, "class_draw")}.
#' @param fit An MNL or constant-share LCMNL klue fit.
#' @param N Number of respondents.
#' @export
klue_draw_betas <- function(fit, N) {
  if (fit$model_type == "MNL")
    return(matrix(fit$par, N, length(fit$par), byrow = TRUE, dimnames = list(NULL, names(fit$par))))
  L <- .fit_layout(fit); nm <- names(fit$par)
  lab <- sub("_class1$", "", nm[grepl("_class1$", nm)])
  B <- t(vapply(seq_len(L$C), function(c) fit$par[paste0(lab, "_class", c)], numeric(length(lab))))
  cls <- sample.int(L$C, N, replace = TRUE, prob = .klue_shares(fit$par, L$delta))
  out <- B[cls, , drop = FALSE]; dimnames(out) <- list(NULL, lab); attr(out, "class_draw") <- cls
  out
}

#' Expected choice probabilities from a fitted model
#'
#' Returns the \code{n_obs} x \code{J} matrix of choice probabilities of every
#' task row of \code{database}: the MNL probabilities, or for a latent class
#' fit their mixture over the classes with the class shares.
#' @param fit An MNL or constant-share LCMNL klue fit.
#' @param database Canonical database (\code{CHOICE} is not used).
#' @param dgp The dgp of the fit.
#' @export
klue_predict <- function(fit, database, dgp) {
  Xf <- .design_full(database, dgp); n <- nrow(database); npc <- dgp$npc
  cls <- if (fit$model_type == "MNL") 1L else seq_len(.fit_layout(fit)$C)
  pr <- lapply(cls, function(c) .row_probs(Xf, matrix(fit$par[(c - 1L) * npc + seq_len(npc)], n, npc, byrow = TRUE)))
  if (length(cls) == 1L) return(pr[[1]])
  pi <- .klue_shares(fit$par, .fit_layout(fit)$delta)
  Reduce(`+`, Map(`*`, pr, pi))
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


