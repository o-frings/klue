# MNL/LCMNL estimation by direct maximum likelihood (BFGS, analytic
# gradients) and by EM. Both maximise the same LCMNL log-likelihood; MNL is
# the C = 1 special case. One shared MNL kernel serves the C=1 path, each
# class inside the C>=2 likelihood, the EM M-step, and fit_cluster_mnls.
# Parameter layout per class: [beta_1..beta_{n_beta}, asc_1..asc_{J-1}],
# class-share deltas appended (last class = reference).

# Design matrices X[[j]] (n_obs x n_beta) per alternative.
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

# Everything the likelihood needs, computed once per estimation call.
# Data are sorted by respondent in blocks of T_per_n rows (the format all
# klue generators and database builders produce).
.lcmnl_context <- function(database, dgp) {
  N <- length(unique(database$ID))
  T_total <- nrow(database)
  J <- dgp$n_alternatives; n_beta <- dgp$n_beta
  X <- build_design_matrices(database, dgp)
  ch <- database$CHOICE
  ch_ind <- matrix(0, nrow = T_total, ncol = J)
  for (j in 1:J) ch_ind[, j] <- as.numeric(ch == j)
  Xc <- matrix(0, nrow = T_total, ncol = n_beta)         # chosen-alt attributes
  for (j in 1:J) Xc <- Xc + ch_ind[, j] * X[[j]]
  list(X = X, ch_ind = ch_ind, Xc = Xc, N = N, T_total = T_total,
       T_per_n = as.integer(T_total / N), J = J, n_beta = n_beta,
       n_asc = dgp$n_asc, npc = dgp$npc)
}

# Per-task log-likelihood (and choice probabilities) for one class parameter
# vector par = c(betas, ascs[1:n_asc]); reference-alt ASC = 0.
.mnl_eval <- function(ctx, par, probs = FALSE) {
  betas <- par[1:ctx$n_beta]
  ascs  <- c(par[(ctx$n_beta + 1):ctx$npc], 0)
  V <- matrix(0, ctx$T_total, ctx$J)
  for (j in 1:ctx$J) V[, j] <- ctx$X[[j]] %*% betas + ascs[j]
  Vm <- do.call(pmax, lapply(1:ctx$J, function(j) V[, j]))
  eV <- exp(V - Vm)
  out <- list(tll = rowSums(ctx$ch_ind * V) - Vm - log(rowSums(eV)))
  if (probs) out$probs <- eV / rowSums(eV)
  out
}

# MNL score per task: chosen_x - E[x] for betas, chosen_ind - prob for ASCs.
.mnl_score <- function(ctx, probs) {
  EX <- matrix(0, ctx$T_total, ctx$n_beta)
  for (j in 1:ctx$J) EX <- EX + probs[, j] * ctx$X[[j]]
  list(beta = ctx$Xc - EX,
       asc  = ctx$ch_ind[, 1:ctx$n_asc, drop = FALSE] -
              probs[, 1:ctx$n_asc, drop = FALSE])
}

# Sum a per-task column over each respondent's T_per_n rows -> length N.
.panel_sum <- function(x, ctx) colSums(matrix(x, ctx$T_per_n, ctx$N))

# Single-class (optionally row-weighted) MNL fit via BFGS.
.fit_mnl <- function(ctx, par0, rw = NULL) {
  neg_ll <- function(par) {
    tll <- .mnl_eval(ctx, par)$tll
    if (is.null(rw)) -sum(tll) else -sum(rw * tll)
  }
  grad_ll <- function(par) {
    e <- .mnl_eval(ctx, par, probs = TRUE)
    s <- .mnl_score(ctx, e$probs)
    if (is.null(rw)) c(-colSums(s$beta), -colSums(s$asc))
    else c(-colSums(rw * s$beta), -colSums(rw * s$asc))
  }
  suppressWarnings(optim(par0, neg_ll, gr = grad_ll, method = "BFGS",
                         control = list(maxit = MAX_ITER, reltol = 1e-10)))
}

# Per-respondent panel log-likelihood matrix (N x C) for stacked parameters.
.panel_loglik <- function(ctx, par, C) {
  log_panel <- matrix(0, ctx$N, C)
  for (ci in 1:C) {
    tll <- .mnl_eval(ctx, par[(ci - 1L) * ctx$npc + 1:ctx$npc])$tll
    log_panel[, ci] <- .panel_sum(tll, ctx)
  }
  log_panel
}

# Row-wise softmax of log_panel + log_pi: posterior class probabilities and
# the sample log-likelihood.
.posterior_weights <- function(log_panel, log_pi) {
  log_joint <- sweep(log_panel, 2, log_pi, "+")
  lm <- apply(log_joint, 1, max)
  denom <- lm + log(rowSums(exp(log_joint - lm)))
  list(w = exp(log_joint - denom), LL = sum(denom))
}

.lcmnl_fail <- function(C, N, n_beta, extra = NULL) {
  c(list(converged = FALSE, C = C, LL = -Inf, BIC = Inf, AIC = Inf, ICL = Inf,
         ICL_BIC = NA_real_, k = 0, betas = matrix(0, C, n_beta),
         class_probs = rep(1 / C, C), posteriors = matrix(1 / C, N, C)),
    extra)
}

#' Direct-MLE LCMNL (MNL when C = 1)
#'
#' BFGS on the full parameter vector with analytic gradients. Bypasses Apollo
#' entirely: LCMNL involves only discrete mixing.
#' @param database Data frame in long format, sorted by respondent in blocks of
#'   \code{T_per_n} rows, with columns \code{ID}, \code{CHOICE}, and the
#'   attribute columns named by the DGP (\code{x*_j}, \code{price_j}).
#' @param C Integer number of latent classes. \code{C = 1} fits a plain MNL.
#' @param start_betas Optional \code{C} x \code{n_beta} matrix of starting taste
#'   coefficients. If \code{NULL}, k-means starting values are computed via
#'   \code{klue_starts}.
#' @param start_shares Optional length-\code{C} vector of starting class shares.
#'   If \code{NULL}, equal shares (\code{1 / C}) are used.
#' @param dgp Data-generating-process specification list giving the design
#'   dimensions (\code{n_alternatives}, \code{n_beta}, \code{n_generic},
#'   \code{n_asc}, \code{npc}). Defaults to \code{DGP_DEFAULT}.
#' @return A list with the fit. \code{converged} (logical) flags a successful
#'   optimisation; \code{C} is the number of classes; \code{LL} the maximised
#'   log-likelihood; \code{BIC}, \code{AIC}, and \code{ICL} the corresponding
#'   information criteria (\code{ICL_BIC} is the entropy penalty
#'   \code{ICL - BIC}); \code{k} the number of free parameters; \code{betas} a
#'   \code{C} x \code{n_beta} matrix of estimated taste coefficients;
#'   \code{class_probs} the length-\code{C} class shares; and \code{posteriors}
#'   the \code{N} x \code{C} matrix of posterior class-membership probabilities.
#'   \code{model_type} is \code{"MNL"} when \code{C = 1} and \code{"LCMNL"}
#'   otherwise. A failed fit returns the same fields with \code{converged} set
#'   to \code{FALSE}.
#' @export
estimate_lcmnl <- function(database, C, start_betas = NULL, start_shares = NULL,
                           dgp = DGP_DEFAULT) {
  ctx <- .lcmnl_context(database, dgp)
  N <- ctx$N; npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc

  if (is.null(start_betas)) {
    starts <- klue_starts(database, C, "kmeans", dgp = dgp)
    start_betas  <- starts$betas
    start_shares <- starts$shares
  }
  if (is.null(start_shares)) start_shares <- rep(1 / C, C)

  fail_result <- .lcmnl_fail(C, N, n_beta)

  if (C == 1) {
    par0 <- c(start_betas[1, ], rep(0, n_asc))
    n_free <- npc
    result <- tryCatch(.fit_mnl(ctx, par0), error = function(e) {
      message("[MNL] optim error: ", conditionMessage(e)); NULL
    })
    if (is.null(result) || result$convergence != 0) return(fail_result)
    LL <- -result$value
    BIC <- -2 * LL + n_free * log(N)
    return(list(converged = TRUE, C = 1L, model_type = "MNL",
                LL = LL, BIC = BIC, AIC = -2 * LL + 2 * n_free,
                ICL = BIC, ICL_BIC = 0, k = n_free,
                betas = matrix(result$par[1:n_beta], nrow = 1),
                class_probs = 1, posteriors = matrix(1, nrow = N, ncol = 1)))
  }

  par0 <- numeric(C * npc + C - 1L)
  for (ci in 1:C) par0[(ci - 1L) * npc + 1:n_beta] <- start_betas[ci, ]
  for (ci in 1:(C - 1L)) {
    par0[C * npc + ci] <- log(max(start_shares[ci], 0.01) /
                              max(start_shares[C], 0.01))
  }
  n_free <- C * npc + C - 1L

  log_pi_of <- function(deltas) {
    dm <- max(deltas)
    deltas - dm - log(sum(exp(deltas - dm)))
  }
  neg_ll <- function(par) {
    log_pi <- log_pi_of(c(par[(C * npc + 1L):(C * npc + C - 1L)], 0))
    log_panel <- .panel_loglik(ctx, par, C)
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
      e <- .mnl_eval(ctx, par[(ci - 1L) * npc + 1:npc], probs = TRUE)
      log_panel[, ci] <- .panel_sum(e$tll, ctx)
      s <- .mnl_score(ctx, e$probs)
      ps <- matrix(0, N, npc)
      for (k in 1:n_beta) ps[, k] <- .panel_sum(s$beta[, k], ctx)
      for (k in 1:n_asc) ps[, n_beta + k] <- .panel_sum(s$asc[, k], ctx)
      pscores[[ci]] <- ps
    }
    w <- .posterior_weights(log_panel, log(pi_c))$w
    g <- numeric(n_free)
    for (ci in 1:C) g[(ci - 1L) * npc + 1:npc] <- -colSums(w[, ci] * pscores[[ci]])
    for (ci in 1:(C - 1L)) g[C * npc + ci] <- -(sum(w[, ci]) - N * pi_c[ci])
    g
  }

  result <- tryCatch(
    suppressWarnings(optim(par0, neg_ll, gr = grad_ll, method = "BFGS",
                           control = list(maxit = MAX_ITER, reltol = 1e-10))),
    error = function(e) {
      message("[LCMNL C=", C, "] optim error: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(result) || result$convergence != 0) return(fail_result)

  p <- result$par
  LL <- -result$value
  betas_mat <- t(vapply(1:C, function(ci) p[(ci - 1L) * npc + 1:n_beta],
                        numeric(n_beta)))
  deltas <- c(p[(C * npc + 1L):(C * npc + C - 1L)], 0)
  exp_d <- exp(deltas - max(deltas))
  class_probs <- as.numeric(exp_d / sum(exp_d))

  # Posteriors at the optimum (Bayes' rule in log space, same formula as the
  # 0.6.x compute_lc_posteriors helper).
  posteriors <- .posterior_weights(.panel_loglik(ctx, p, C), log(class_probs))$w
  posteriors_c <- pmax(posteriors, 1e-100)
  H <- -sum(posteriors_c * log(posteriors_c))
  BIC <- -2 * LL + n_free * log(N)

  list(converged = TRUE, C = C, model_type = "LCMNL",
       LL = LL, BIC = BIC, AIC = -2 * LL + 2 * n_free,
       ICL = BIC + 2 * H, ICL_BIC = 2 * H,
       k = n_free, betas = betas_mat, class_probs = class_probs,
       posteriors = posteriors_c)
}

#' EM estimator for the same LCMNL likelihood
#'
#' E-step: posterior class weights; M-step: shares + per-class weighted MNL
#' (warm-started at the class's current parameters). From identical starts,
#' EM and direct ML should reach the same optimum. C = 1 delegates to
#' \code{estimate_lcmnl} so the two estimators agree exactly there.
#' @param database Data frame in long format, sorted by respondent in blocks of
#'   \code{T_per_n} rows, with columns \code{ID}, \code{CHOICE}, and the
#'   attribute columns named by the DGP (\code{x*_j}, \code{price_j}).
#' @param C Integer number of latent classes. \code{C = 1} delegates to
#'   \code{estimate_lcmnl}.
#' @param start_betas Optional \code{C} x \code{n_beta} matrix of starting taste
#'   coefficients. If \code{NULL}, k-means starting values are computed via
#'   \code{klue_starts}.
#' @param start_shares Optional length-\code{C} vector of starting class shares.
#'   If \code{NULL}, equal shares (\code{1 / C}) are used.
#' @param dgp Data-generating-process specification list giving the design
#'   dimensions (\code{n_alternatives}, \code{n_beta}, \code{n_generic},
#'   \code{n_asc}, \code{npc}). Defaults to \code{DGP_DEFAULT}.
#' @param max_em_iter Maximum number of EM iterations. Default 500.
#' @param tol Convergence tolerance on the change in log-likelihood between
#'   successive EM iterations. Default 1e-6.
#' @param verbose If \code{TRUE}, print the log-likelihood at each iteration.
#'   Default \code{FALSE}.
#' @return A list with the fit. \code{converged} (logical) flags a successful
#'   run; \code{C} is the number of classes; \code{LL} the maximised
#'   log-likelihood; \code{BIC}, \code{AIC}, and \code{ICL} the corresponding
#'   information criteria (\code{ICL_BIC} is the entropy penalty
#'   \code{ICL - BIC}); \code{k} the number of free parameters; \code{betas} a
#'   \code{C} x \code{n_beta} matrix of estimated taste coefficients;
#'   \code{class_probs} the length-\code{C} class shares; \code{posteriors} the
#'   \code{N} x \code{C} matrix of posterior class-membership probabilities;
#'   \code{em_iters} the number of EM iterations run; and \code{estimator} the
#'   string \code{"em"}. \code{model_type} is \code{"MNL"} when \code{C = 1} and
#'   \code{"LCMNL"} otherwise. A failed fit returns the same fields with
#'   \code{converged} set to \code{FALSE}.
#' @export
estimate_lcmnl_em <- function(database, C, start_betas = NULL,
                              start_shares = NULL, dgp = DGP_DEFAULT,
                              max_em_iter = 500L, tol = 1e-6, verbose = FALSE) {
  if (is.null(start_betas)) {
    starts <- klue_starts(database, C, "kmeans", dgp = dgp)
    start_betas  <- starts$betas
    start_shares <- starts$shares
  }
  if (is.null(start_shares)) start_shares <- rep(1 / C, C)

  if (C == 1L) {
    res <- estimate_lcmnl(database, C, start_betas, start_shares, dgp = dgp)
    res$em_iters <- 0L; res$estimator <- "em"
    return(res)
  }

  ctx <- .lcmnl_context(database, dgp)
  N <- ctx$N; npc <- ctx$npc; n_beta <- ctx$n_beta
  fail_result <- .lcmnl_fail(C, N, n_beta,
                             extra = list(em_iters = 0L, estimator = "em"))
  row_resp <- rep(seq_len(N), each = ctx$T_per_n)

  par_c <- matrix(0, C, npc)
  for (ci in 1:C) par_c[ci, 1:n_beta] <- start_betas[ci, ]   # ASCs start at 0
  pi_c <- pmax(start_shares, 1e-6); pi_c <- pi_c / sum(pi_c)

  LL_prev <- -Inf; w <- matrix(1 / C, N, C); em_iters <- 0L; LL <- -Inf
  for (it in seq_len(max_em_iter)) {
    em_iters <- it
    # E-step
    log_panel <- matrix(0, N, C)
    for (ci in 1:C) log_panel[, ci] <- .panel_sum(.mnl_eval(ctx, par_c[ci, ])$tll, ctx)
    pw <- .posterior_weights(log_panel, log(pi_c))
    w <- pw$w; LL <- pw$LL
    if (!is.finite(LL)) return(fail_result)
    if (verbose) cat(sprintf("  [EM %3d] LL = %.4f\n", it, LL))
    if (abs(LL - LL_prev) < tol) break
    LL_prev <- LL
    # M-step
    pi_c <- pmax(colMeans(w), 1e-8); pi_c <- pi_c / sum(pi_c)
    for (ci in 1:C) {
      rw <- w[row_resp, ci]
      if (sum(rw) < 1e-6) next   # collapsed class: keep current parameters
      opt <- tryCatch(.fit_mnl(ctx, par_c[ci, ], rw = rw),
                      error = function(e) NULL)
      if (!is.null(opt)) par_c[ci, ] <- opt$par
    }
  }

  n_free <- C * npc + (C - 1L)
  BIC <- -2 * LL + n_free * log(N)
  post <- pmax(w, 1e-100)
  H <- -sum(post * log(post))

  list(converged = TRUE, C = C, model_type = "LCMNL",
       LL = LL, BIC = BIC, AIC = -2 * LL + 2 * n_free,
       ICL = BIC + 2 * H, ICL_BIC = 2 * H,
       k = n_free, betas = par_c[, 1:n_beta, drop = FALSE],
       class_probs = as.numeric(pi_c), posteriors = w,
       em_iters = em_iters, estimator = "em")
}

#' Multi-start LCMNL: best of the six clustering initialisations
#'
#' For C = 1 (MNL) a single run suffices; for C >= 2 the model is estimated
#' from all six clustering starts and the best converged log-likelihood wins.
#' @param database Data frame in long format, sorted by respondent in blocks of
#'   \code{T_per_n} rows, with columns \code{ID}, \code{CHOICE}, and the
#'   attribute columns named by the DGP (\code{x*_j}, \code{price_j}).
#' @param C Integer number of latent classes. \code{C = 1} fits a plain MNL.
#' @param dgp Data-generating-process specification list giving the design
#'   dimensions (\code{n_alternatives}, \code{n_beta}, \code{n_generic},
#'   \code{n_asc}, \code{npc}). Defaults to \code{DGP_DEFAULT}.
#' @param estimator "ml" (direct BFGS, default) or "em".
#' @param feature_type "rp" (default) or "onehot" clustering features.
#' @param n_cores number of cores for the six per-start fits. Default 1
#'   (sequential). With \code{n_cores > 1} the starts are fit concurrently via
#'   \code{parallel::mclapply} (fork-based; not available on Windows). The
#'   best-of-six selection is order-deterministic and so is independent of
#'   \code{n_cores}. Leave at 1 inside the study drivers, which already
#'   parallelise across conditions -- nesting would oversubscribe the cores.
#' @return The best-fitting per-start result, a list with the same fields as
#'   the chosen estimator (\code{estimate_lcmnl} for \code{"ml"},
#'   \code{estimate_lcmnl_em} for \code{"em"}): \code{converged}, \code{C},
#'   \code{LL}, \code{BIC}, \code{AIC}, \code{ICL} (and \code{ICL_BIC}),
#'   \code{k}, \code{betas}, \code{class_probs}, \code{posteriors}, and for the
#'   EM estimator \code{em_iters} and \code{estimator}. Two extra fields are
#'   added: \code{best_method}, the name of the winning clustering start, and
#'   \code{method_results}, a named list of the individual fits from every
#'   start. If no start converges, the failure result is returned with
#'   \code{converged = FALSE} and \code{best_method = NA}.
#' @export
klue_lcmnl <- function(database, C, dgp = DGP_DEFAULT,
                       estimator = c("ml", "em"),
                       feature_type = c("rp", "onehot"),
                       n_cores = 1L) {
  estimator <- match.arg(estimator)
  feature_type <- match.arg(feature_type)
  fit_one <- if (estimator == "em") estimate_lcmnl_em else estimate_lcmnl
  all_starts <- get_all_starts(database, C, dgp = dgp,
                               feature_type = feature_type)
  N <- length(unique(database$ID))

  # The six per-start fits are independent and ~90% of the runtime; optionally
  # fit them concurrently. Order is preserved so the best-of-six tie-break is
  # identical to the sequential path.
  nm_ok <- names(all_starts)[!vapply(all_starts, is.null, logical(1))]
  fit_start <- function(nm) tryCatch(
    fit_one(database, C, start_betas = all_starts[[nm]]$betas,
            start_shares = all_starts[[nm]]$shares, dgp = dgp),
    error = function(e) NULL)
  fits <- if (n_cores > 1L && length(nm_ok) > 1L) {
    parallel::mclapply(nm_ok, fit_start, mc.cores = min(as.integer(n_cores),
                                                        length(nm_ok)))
  } else {
    lapply(nm_ok, fit_start)
  }
  names(fits) <- nm_ok

  best <- .lcmnl_fail(C, N, dgp$n_beta,
                      extra = list(best_method = NA_character_))
  method_results <- list()
  for (nm in nm_ok) {                 # sequential scan -> deterministic winner
    res <- fits[[nm]]
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
