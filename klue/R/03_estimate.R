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

# Per-alternative design matrices with the constants: n_obs x npc, the taste
# columns then one indicator column per ASC group (the per-class parameter layout).
.design_full <- function(database, dgp) {
  X <- build_design_matrices(database, dgp); n <- nrow(database)
  lapply(seq_len(dgp$n_alternatives), function(j) {
    A <- matrix(0, n, dgp$n_asc); g <- dgp$asc_map[j]
    if (g > 0L) A[, g] <- 1
    cbind(X[[j]], A)
  })
}

# Utilities (n_obs x J) and choice probabilities for row-wise coefficients B (n_obs x npc).
.row_utilities <- function(Xf, B) matrix(vapply(Xf, function(Xj) rowSums(Xj * B), numeric(nrow(B))), nrow(B))
.row_probs <- function(Xf, B) {
  V <- .row_utilities(Xf, B); V <- V - apply(V, 1, max); E <- exp(V); E / rowSums(E)
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
  # Optional per-task availability: wide columns av_1..av_J (1 = available,
  # 0 = not). Absent, or present but all ones -> full availability, and the
  # estimator takes the unchanged (default) code path. When some alternative is
  # unavailable it is masked out of that task's choice set.
  A <- NULL; has_avail <- FALSE
  av_cols <- paste0("av_", 1:J)
  if (all(av_cols %in% names(database))) {
    A <- vapply(av_cols, function(cc) as.numeric(database[[cc]]), numeric(T_total))
    if (any(A == 0)) {
      has_avail <- TRUE
      if (any(rowSums(A * ch_ind) == 0))
        stop("chosen alternative marked unavailable (av_j = 0) in some tasks")
    } else A <- NULL
  }
  list(X = X, ch_ind = ch_ind, Xc = Xc, N = N, T_total = T_total,
       T_per_n = as.integer(T_total / N), J = J, n_beta = n_beta,
       n_asc = dgp$n_asc, npc = dgp$npc, asc_map = dgp$asc_map,
       A = A, has_avail = has_avail)
}

# Per-task log-likelihood (and choice probabilities) for one class parameter
# vector par = c(betas, asc[1:n_asc]); ASCs are shared within groups per
# ctx$asc_map (group 0 = reference alternative, ASC fixed at 0).
.mnl_eval <- function(ctx, par, probs = FALSE) {
  betas <- par[1:ctx$n_beta]
  asc_par <- if (ctx$n_asc > 0) par[(ctx$n_beta + 1):ctx$npc] else numeric(0)
  m <- ctx$asc_map
  asc_by_alt <- numeric(ctx$J)
  asc_by_alt[m > 0] <- asc_par[m[m > 0]]              # map group ASC onto its alternatives
  V <- matrix(0, ctx$T_total, ctx$J)
  for (j in 1:ctx$J) V[, j] <- ctx$X[[j]] %*% betas + asc_by_alt[j]
  if (isTRUE(ctx$has_avail)) {
    # Mask unavailable alternatives for the normaliser (Vc), but keep the
    # unmasked V for the chosen-utility term to avoid 0 * -Inf = NaN.
    Vc <- V; Vc[ctx$A == 0] <- -Inf
    Vm <- do.call(pmax, lapply(1:ctx$J, function(j) Vc[, j]))
    eV <- exp(Vc - Vm)                       # unavailable alternatives -> 0
    out <- list(tll = rowSums(ctx$ch_ind * V) - Vm - log(rowSums(eV)))
  } else {
    Vm <- do.call(pmax, lapply(1:ctx$J, function(j) V[, j]))
    eV <- exp(V - Vm)
    out <- list(tll = rowSums(ctx$ch_ind * V) - Vm - log(rowSums(eV)))
  }
  if (probs) out$probs <- eV / rowSums(eV)
  out
}

# MNL score per task: chosen_x - E[x] for betas; for ASC group g, the sum of
# (chosen_ind - prob) over the alternatives sharing that group.
.mnl_score <- function(ctx, probs) {
  EX <- matrix(0, ctx$T_total, ctx$n_beta)
  for (j in 1:ctx$J) EX <- EX + probs[, j] * ctx$X[[j]]
  resid <- ctx$ch_ind - probs                          # T_total x J
  asc <- matrix(0, ctx$T_total, ctx$n_asc)
  for (g in seq_len(ctx$n_asc))
    asc[, g] <- rowSums(resid[, ctx$asc_map == g, drop = FALSE])
  list(beta = ctx$Xc - EX, asc = asc)
}

# Sum a per-task column over each respondent's T_per_n rows -> length N.
.panel_sum <- function(x, ctx) colSums(matrix(x, ctx$T_per_n, ctx$N))

# Names of the full parameter vector: the labels (C = 1), or <label>_class<c>
# for every class followed by the share parameters delta1..delta{C-1}.
.full_pnames <- function(labels, C)
  if (C == 1L) labels else c(paste0(rep(labels, C), "_class", rep(seq_len(C), each = length(labels))),
                             .lab("delta", C - 1L))

# Coefficients common to all classes: the full parameter vector is M %*% theta,
# where theta holds one entry per common coefficient (named by its label) and
# one per other parameter (named as in the full vector). `extra` names the
# parameters after the C taste blocks (share deltas or membership slopes).
# NULL when C = 1 or nothing is common, so the default path is unchanged.
.common_map <- function(dgp, C, extra) {
  cm <- dgp$common
  if (C < 2L || !length(cm)) return(NULL)
  l <- .par_labels(dgp); np <- length(l)
  key <- c(unlist(lapply(seq_len(C), function(ci)
           ifelse(seq_len(np) %in% cm, l, paste0(l, "_class", ci)))), extra)
  th <- unique(key)
  list(M = outer(key, th, `==`) * 1, theta_names = th, common = l[cm])
}

# Optimise `obj` over theta (par = M theta) or, without a map, over par itself.
.optim_mapped <- function(obj, par0, cmap, maxit, reltol) {
  if (is.null(cmap))
    return(suppressWarnings(optim(par0, obj$neg_ll, gr = obj$grad_ll, method = "BFGS",
                                  control = list(maxit = maxit, reltol = reltol))))
  M <- cmap$M
  th0 <- drop(solve(crossprod(M), crossprod(M, par0)))   # a common entry starts at the mean of its class starts
  res <- suppressWarnings(optim(th0, function(t) obj$neg_ll(drop(M %*% t)),
                                gr = function(t) drop(crossprod(M, obj$grad_ll(drop(M %*% t)))),
                                method = "BFGS", control = list(maxit = maxit, reltol = reltol)))
  res$par <- drop(M %*% res$par)
  res
}

# Covariance with common coefficients: observed information and sandwich in
# theta, mapped back to the full vector (M V M'), so functions of the full
# parameter vector keep working. G holds the per-respondent full-vector scores.
.mapped_vcov <- function(par, cmap, obj, G, pnames) {
  M <- cmap$M; nf <- length(pnames)
  na <- matrix(NA_real_, nf, nf, dimnames = list(pnames, pnames))
  th <- drop(solve(crossprod(M), crossprod(M, par)))
  H <- tryCatch(stats::optimHess(th, function(t) obj$neg_ll(drop(M %*% t)),
                                 function(t) drop(crossprod(M, obj$grad_ll(drop(M %*% t))))),
                error = function(e) NULL)
  aliased <- .aliased_params(H, cmap$theta_names)
  Vt <- if (is.null(H) || length(aliased)) NULL
        else tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
  if (is.null(Vt)) return(list(vcov = na, robust_vcov = na, pnames = pnames, aliased = aliased))
  Gt <- G %*% M
  full <- function(A) { B <- M %*% A %*% t(M); dimnames(B) <- list(pnames, pnames); B }
  list(vcov = full(Vt), robust_vcov = full(Vt %*% crossprod(Gt) %*% Vt),
       pnames = pnames, aliased = aliased)
}

# Single-class (optionally row-weighted) MNL fit via BFGS. `hessian = TRUE`
# additionally returns the observed-information Hessian at the optimum.
.fit_mnl <- function(ctx, par0, rw = NULL, maxit = MAX_ITER, reltol = 1e-10,
                     hessian = FALSE) {
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
                         control = list(maxit = maxit, reltol = reltol),
                         hessian = hessian))
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

# Factory for the C >= 2 LCMNL negative log-likelihood and its analytic
# gradient, shared by estimate_lcmnl (optimisation) and .lcmnl_vcov (observed
# information) so the likelihood has a single definition. par layout:
# [class1(npc), .., classC(npc), delta_1..delta_{C-1}] (last class reference).
.lcmnl_objective <- function(ctx, C) {
  npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc; N <- ctx$N
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
      for (k in seq_len(n_asc)) ps[, n_beta + k] <- .panel_sum(s$asc[, k], ctx)
      pscores[[ci]] <- ps
    }
    w <- .posterior_weights(log_panel, log(pi_c))$w
    g <- numeric(n_free)
    for (ci in 1:C) g[(ci - 1L) * npc + 1:npc] <- -colSums(w[, ci] * pscores[[ci]])
    for (ci in 1:(C - 1L)) g[C * npc + ci] <- -(sum(w[, ci]) - N * pi_c[ci])
    g
  }
  list(neg_ll = neg_ll, grad_ll = grad_ll, n_free = n_free)
}

# Covariance of a converged C >= 2 LCMNL fit at parameter vector `par`.
# vcov = inverse observed information (numerical Hessian of the negative
# log-likelihood formed from the analytic gradient); robust_vcov = sandwich
# clustered on respondent, from the analytic per-respondent scores of the
# LCMNL log-likelihood. Returns NA-filled matrices if the information is
# singular (e.g. a rank-deficient design), matching the C = 1 path.
# One message for every singular-covariance path, naming the parameters that are
# not separately identified so the caller can respecify rather than guess.
.singular_msg <- function(what, C, aliased) {
  head <- if (identical(what, "MNL")) "[MNL]" else sprintf("[%s C=%d]", what, C)
  det <- if (length(aliased) && length(aliased) <= 12)
           sprintf(" Not separately identified: %s.", paste(aliased, collapse = ", "))
         else ""
  sprintf(paste0("%s information matrix is singular: covariance and standard errors set to NA.%s",
                 " Check for collinear attributes (e.g. dummy-coded levels that sum to a constant",
                 " alongside an alternative-specific constant). Estimates, log-likelihood and",
                 " class assignment are unaffected; the reported k counts every parameter."),
          head, det)
}

# Parameters that lie on a null direction of an observed-information matrix.
# A design with an exact linear dependency -- dummy-coded levels that sum to a
# constant alongside an ASC, say -- leaves an eigenvalue at ~0. solve() still
# "inverts" such a matrix (LAPACK only errors on exact computational
# singularity), returning huge negative variances and NaN standard errors while
# nothing signals a problem, so rank is tested here instead. Returns the names
# of the parameters loading on a null or negative-curvature direction, or
# character(0) when the information matrix is positive definite.
.aliased_params <- function(H, pnames, tol = sqrt(.Machine$double.eps)) {
  if (is.null(H) || anyNA(H)) return(pnames)
  H <- (H + t(H)) / 2
  ev <- tryCatch(eigen(H, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev)) return(pnames)
  bad <- which(ev$values <= tol * max(abs(ev$values)))
  if (!length(bad)) return(character(0))
  load <- sqrt(rowSums(ev$vectors[, bad, drop = FALSE]^2))
  al <- pnames[load > 0.1]
  if (length(al)) al else pnames
}

.lcmnl_vcov <- function(ctx, par, C, labels = c(.lab("b", ctx$n_beta), .lab("asc", ctx$n_asc)),
                        cmap = NULL) {
  obj <- .lcmnl_objective(ctx, C); n_free <- obj$n_free
  pnames <- .full_pnames(labels, C)
  if (!is.null(cmap)) return(.mapped_vcov(par, cmap, obj, .lcmnl_scores(ctx, par, C), pnames))
  na <- matrix(NA_real_, n_free, n_free, dimnames = list(pnames, pnames))
  H <- tryCatch(stats::optimHess(par, obj$neg_ll, obj$grad_ll),
                error = function(e) NULL)
  aliased <- .aliased_params(H, pnames)
  vcov <- if (is.null(H) || length(aliased)) na
          else tryCatch({ V <- solve((H + t(H)) / 2)
                          dimnames(V) <- list(pnames, pnames); V },
                        error = function(e) na)
  G <- .lcmnl_scores(ctx, par, C)
  robust_vcov <- tryCatch({ RV <- vcov %*% crossprod(G) %*% vcov
                            dimnames(RV) <- list(pnames, pnames); RV },
                          error = function(e) na)
  list(vcov = vcov, robust_vcov = robust_vcov, pnames = pnames, aliased = aliased)
}

# Per-respondent scores G (N x n_free) of the LCMNL log-likelihood. For class
# c the taste block is w_nc * (panel-summed MNL score); for share delta_j the
# entry is w_nj - pi_j. colSums(G) equals the total score (~0 at the optimum).
.lcmnl_scores <- function(ctx, par, C) {
  npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc; N <- ctx$N
  n_free <- C * npc + C - 1L
  deltas <- c(par[(C * npc + 1L):(C * npc + C - 1L)], 0)
  dm <- max(deltas); ed <- exp(deltas - dm); pi_c <- ed / sum(ed)
  log_panel <- matrix(0, N, C); pscores <- vector("list", C)
  for (ci in 1:C) {
    e <- .mnl_eval(ctx, par[(ci - 1L) * npc + 1:npc], probs = TRUE)
    log_panel[, ci] <- .panel_sum(e$tll, ctx)
    s <- .mnl_score(ctx, e$probs)
    ps <- matrix(0, N, npc)
    for (k in 1:n_beta) ps[, k] <- .panel_sum(s$beta[, k], ctx)
    for (k in seq_len(n_asc)) ps[, n_beta + k] <- .panel_sum(s$asc[, k], ctx)
    pscores[[ci]] <- ps
  }
  w <- .posterior_weights(log_panel, log(pi_c))$w
  G <- matrix(0, N, n_free)
  for (ci in 1:C) G[, (ci - 1L) * npc + 1:npc] <- w[, ci] * pscores[[ci]]
  for (j in seq_len(C - 1L)) G[, C * npc + j] <- w[, j] - pi_c[j]
  G
}

# --- Concomitant (covariate-driven) class membership -------------------------
# Generalises the constant-share model: class-membership probabilities become
# pi_{nc} = softmax_c(t(Z_n) gamma_c), Z_n a respondent-level design (intercept +
# covariates), gamma_C = 0 (reference). p = 1 (intercept only) reproduces the
# constant-share model exactly. Taste-parameter scores are unchanged; the
# membership score for class c is (w_nc - pi_nc) * Z_n.

# Respondent-level membership design Z (N x p), intercept first, then the named
# covariate columns (taken from each respondent's first row; they must be
# constant within respondent).
.lcmnl_membership_design <- function(database, ctx, membership) {
  idx <- (seq_len(ctx$N) - 1L) * ctx$T_per_n + 1L
  Z <- matrix(1, ctx$N, 1L, dimnames = list(NULL, "(Intercept)"))
  if (length(membership)) {
    covs <- vapply(membership, function(v) as.numeric(database[[v]][idx]),
                   numeric(ctx$N))
    covs <- matrix(covs, nrow = ctx$N, dimnames = list(NULL, membership))
    Z <- cbind(Z, covs)
  }
  Z
}

# Factory: concomitant LCMNL negative log-likelihood, its gradient, and the
# per-respondent score matrix (single source, so gradient and covariance agree).
# par layout: [class1(npc), .., classC(npc), gamma_1(p), .., gamma_{C-1}(p)].
.lcmnl_cov_objective <- function(ctx, C, Z) {
  npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc; N <- ctx$N
  p_mem <- ncol(Z); n_taste <- C * npc; n_free <- n_taste + (C - 1L) * p_mem
  state <- function(par) {
    log_panel <- matrix(0, N, C); pscores <- vector("list", C)
    for (ci in 1:C) {
      e <- .mnl_eval(ctx, par[(ci - 1L) * npc + 1:npc], probs = TRUE)
      log_panel[, ci] <- .panel_sum(e$tll, ctx)
      s <- .mnl_score(ctx, e$probs)
      ps <- matrix(0, N, npc)
      for (k in 1:n_beta) ps[, k] <- .panel_sum(s$beta[, k], ctx)
      for (k in seq_len(n_asc)) ps[, n_beta + k] <- .panel_sum(s$asc[, k], ctx)
      pscores[[ci]] <- ps
    }
    eta <- matrix(0, N, C)
    if (C > 1L) eta[, 1:(C - 1L)] <- Z %*% matrix(par[(n_taste + 1L):n_free],
                                                  nrow = p_mem)
    em <- apply(eta, 1, max); ex <- exp(eta - em); pimat <- ex / rowSums(ex)
    list(log_panel = log_panel, pscores = pscores, pimat = pimat)
  }
  posterior <- function(st) {
    lj <- st$log_panel + log(st$pimat)
    lm <- apply(lj, 1, max)
    denom <- lm + log(rowSums(exp(lj - lm)))
    list(w = exp(lj - denom), LL = sum(denom))
  }
  scores <- function(par) {
    st <- state(par); w <- posterior(st)$w
    G <- matrix(0, N, n_free)
    for (ci in 1:C) G[, (ci - 1L) * npc + 1:npc] <- w[, ci] * st$pscores[[ci]]
    if (C > 1L) for (ci in 1:(C - 1L))
      G[, n_taste + (ci - 1L) * p_mem + 1:p_mem] <- (w[, ci] - st$pimat[, ci]) * Z
    G
  }
  neg_ll <- function(par) -posterior(state(par))$LL
  list(neg_ll = neg_ll, grad_ll = function(par) -colSums(scores(par)),
       scores = scores, state = state, posterior = posterior,
       n_free = n_free, n_taste = n_taste, p_mem = p_mem)
}

# Parameter names for the concomitant model: taste blocks, then membership slopes.
.lcmnl_cov_mem_names <- function(C, Z)
  if (C > 1L) unlist(lapply(1:(C - 1L), function(cc) paste0(colnames(Z), "_class", cc))) else character(0)
.lcmnl_cov_pnames <- function(labels, C, Z)
  c(paste0(rep(labels, C), "_class", rep(seq_len(C), each = length(labels))), .lcmnl_cov_mem_names(C, Z))

# Covariance for a converged concomitant fit (same construction as .lcmnl_vcov).
.lcmnl_cov_vcov <- function(ctx, C, Z, par, labels = c(.lab("b", ctx$n_beta), .lab("asc", ctx$n_asc)),
                            cmap = NULL) {
  obj <- .lcmnl_cov_objective(ctx, C, Z); n_free <- obj$n_free
  pn <- .lcmnl_cov_pnames(labels, C, Z)
  if (!is.null(cmap)) return(.mapped_vcov(par, cmap, obj, obj$scores(par), pn))
  na <- matrix(NA_real_, n_free, n_free, dimnames = list(pn, pn))
  H <- tryCatch(stats::optimHess(par, obj$neg_ll, obj$grad_ll),
                error = function(e) NULL)
  aliased <- .aliased_params(H, pn)
  vcov <- if (is.null(H) || length(aliased)) na
          else tryCatch({ V <- solve((H + t(H)) / 2)
                          dimnames(V) <- list(pn, pn); V },
                        error = function(e) na)
  G <- obj$scores(par)
  robust_vcov <- tryCatch({ RV <- vcov %*% crossprod(G) %*% vcov
                            dimnames(RV) <- list(pn, pn); RV },
                          error = function(e) na)
  list(vcov = vcov, robust_vcov = robust_vcov, pnames = pn, aliased = aliased)
}

.lcmnl_fail <- function(C, N, n_beta, extra = NULL) {
  c(list(converged = FALSE, C = C, LL = -Inf, BIC = Inf, AIC = Inf, ICL = Inf,
         ICL_BIC = NA_real_, k = 0, betas = matrix(0, C, n_beta),
         class_probs = rep(1 / C, C), posteriors = matrix(1 / C, N, C),
         failed = TRUE, entropy_norm = NA_real_, max_abs_par = NA_real_),
    extra)
}

# Status fields shared by every fit (0.10.0): `failed` when the optimiser did not
# converge or, once a covariance is attached, when it is singular or has NA
# entries; the largest absolute estimate; and the normalised classification
# entropy 1 - H / (N ln C) (NA at C = 1).
.fit_status <- function(fit, N, H = NA_real_) {
  has_v <- !is.null(fit$vcov)
  fit$failed <- !isTRUE(fit$converged) ||
    (has_v && (isTRUE(fit$singular) || is.null(fit$robust_vcov) || anyNA(fit$robust_vcov)))
  fit$max_abs_par <- if (length(fit$par)) max(abs(fit$par)) else NA_real_
  fit$entropy_norm <- if (fit$C >= 2L) 1 - H / (N * log(fit$C)) else NA_real_
  fit
}

#' Direct-MLE LCMNL (MNL when C = 1)
#'
#' BFGS on the full parameter vector with analytic gradients. Bypasses Apollo
#' entirely: LCMNL involves only discrete mixing.
#' @param database Data frame in long format, sorted by respondent in blocks of
#'   \code{T_per_n} rows, with columns \code{ID}, \code{CHOICE}, and the
#'   attribute columns named by the DGP (\code{x*_j}, \code{price_j}). Optional
#'   availability columns \code{av_1..av_J} (1 = available, 0 = not) mask
#'   unavailable alternatives out of each task's choice set; if absent or all
#'   ones, every alternative is available (the default). The chosen alternative
#'   must be available in every task.
#' @param C Integer number of latent classes. \code{C = 1} fits a plain MNL.
#' @param start_betas Optional \code{C} x \code{n_beta} matrix of starting taste
#'   coefficients. If \code{NULL}, k-means starting values are computed via
#'   \code{klue_starts}.
#' @param start_shares Optional length-\code{C} vector of starting class shares.
#'   If \code{NULL}, equal shares (\code{1 / C}) are used.
#' @param dgp Data-generating-process specification list giving the design
#'   dimensions (\code{n_alternatives}, \code{n_beta}, \code{n_generic},
#'   \code{n_asc}, \code{npc}). Defaults to \code{DGP_DEFAULT}.
#' @param maxit,reltol BFGS iteration cap and relative tolerance. The defaults
#'   fit to full precision; the multistart screen passes looser values.
#' @param start_par Optional full starting parameter vector (betas, ASCs, and
#'   share deltas as one vector), used to warm-start from a previous fit.
#'   Overrides \code{start_betas}/\code{start_shares}.
#' @param vcov If \code{TRUE} (default), compute the parameter covariance at the
#'   optimum. Set \code{FALSE} to skip the Hessian on throw-away fits (the
#'   multistart uses this on its screening fits, attaching the covariance only
#'   to the winner).
#' @return A list with the fit. \code{converged} (logical) flags a successful
#'   optimisation; \code{C} is the number of classes; \code{LL} the maximised
#'   log-likelihood; \code{BIC}, \code{AIC}, and \code{ICL} the corresponding
#'   information criteria (\code{ICL_BIC} is the entropy penalty
#'   \code{ICL - BIC}); \code{k} the number of free parameters; \code{betas} a
#'   \code{C} x \code{n_beta} matrix of estimated taste coefficients;
#'   \code{class_probs} the length-\code{C} class shares; \code{posteriors}
#'   the \code{N} x \code{C} matrix of posterior class-membership
#'   probabilities; and \code{par} the full parameter vector (reusable as
#'   \code{start_par}). When \code{vcov = TRUE}, \code{vcov} (inverse observed
#'   information) and \code{robust_vcov} (sandwich estimator clustered on
#'   respondent) give the covariance of \code{par} for any \code{C}; if the
#'   observed information is rank-deficient both are set to \code{NA},
#'   \code{singular} is \code{TRUE} and \code{aliased} names the parameters
#'   that are not separately identified (estimates, \code{LL} and class
#'   assignment are unaffected, but \code{k} still counts every parameter, so
#'   the information criteria are penalised for parameters the design cannot
#'   identify). For
#'   \code{C = 1} the parameters are labelled \code{b1..b{n_beta}} (taste
#'   coefficients, price last) and \code{asc1..asc{n_asc}}; for \code{C >= 2}
#'   each class \code{c} contributes \code{b1_class{c}..b{n_beta}_class{c}} and
#'   \code{asc1_class{c}..asc{n_asc}_class{c}}, followed by the \code{C - 1}
#'   share parameters \code{delta1..delta{C-1}} (last class is the reference).
#'   \code{model_type} is \code{"MNL"} when
#'   \code{C = 1} and \code{"LCMNL"} otherwise. A fit that hits \code{maxit} returns its
#'   (finite) values with \code{converged = FALSE}; an optimiser error returns
#'   the failure skeleton with \code{LL = -Inf}.
#' @export
estimate_lcmnl <- function(database, C, start_betas = NULL, start_shares = NULL,
                           dgp = DGP_DEFAULT, maxit = MAX_ITER,
                           reltol = 1e-10, start_par = NULL, vcov = TRUE) {
  ctx <- .lcmnl_context(database, dgp)
  N <- ctx$N; npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc

  if (is.null(start_par)) {
    if (is.null(start_betas)) {
      starts <- klue_starts(database, C, "kmeans", dgp = dgp)
      start_betas  <- starts$betas
      start_shares <- starts$shares
    }
    if (is.null(start_shares)) start_shares <- rep(1 / C, C)
  }

  fail_result <- .lcmnl_fail(C, N, n_beta)

  if (C == 1) {
    par0 <- if (!is.null(start_par)) start_par
            else c(start_betas[1, ], rep(0, n_asc))
    n_free <- npc
    result <- tryCatch(.fit_mnl(ctx, par0, maxit = maxit, reltol = reltol,
                                hessian = vcov),
                       error = function(e) {
      message("[MNL] optim error: ", conditionMessage(e)); NULL
    })
    if (is.null(result) || !result$convergence %in% c(0L, 1L))
      return(fail_result)
    LL <- -result$value
    BIC <- -2 * LL + n_free * log(N)
    # Covariance of the full parameter vector par = c(betas[1:n_beta],
    # ascs[1:n_asc]). vcov is the inverse observed information; robust_vcov is
    # the sandwich estimator clustered on respondent, built from the analytic
    # per-task scores (.mnl_score) summed within respondent (.panel_sum).
    pnames <- .par_labels(dgp)
    V <- RV <- NULL
    aliased <- character(0)
    if (vcov) {
      aliased <- .aliased_params(result$hessian, pnames)
      V <- if (length(aliased)) matrix(NA_real_, npc, npc)
           else tryCatch(solve(result$hessian),
                         error = function(e) matrix(NA_real_, npc, npc))
      ev <- .mnl_eval(ctx, result$par, probs = TRUE)
      sc <- .mnl_score(ctx, ev$probs)
      G  <- apply(cbind(sc$beta, sc$asc), 2, .panel_sum, ctx = ctx)  # N x npc
      RV <- V %*% crossprod(G) %*% V
      dimnames(V) <- dimnames(RV) <- list(pnames, pnames)
    }
    if (vcov && anyNA(V)) message(.singular_msg("MNL", 1L, aliased))
    return(.fit_status(list(converged = result$convergence == 0L, singular = vcov && anyNA(V),
                aliased = aliased, C = 1L,
                model_type = "MNL",
                LL = LL, BIC = BIC, AIC = -2 * LL + 2 * n_free,
                ICL = BIC, ICL_BIC = 0, k = n_free,
                betas = matrix(result$par[1:n_beta], nrow = 1),
                class_probs = 1, posteriors = matrix(1, nrow = N, ncol = 1),
                par = setNames(result$par, pnames),
                vcov = V, robust_vcov = RV,
                common = character(0), price = pnames[n_beta]), N))
  }

  if (!is.null(start_par)) {
    par0 <- start_par
  } else {
    par0 <- numeric(C * npc + C - 1L)
    for (ci in 1:C) par0[(ci - 1L) * npc + 1:n_beta] <- start_betas[ci, ]
    for (ci in 1:(C - 1L)) {
      par0[C * npc + ci] <- log(max(start_shares[ci], 0.01) /
                                max(start_shares[C], 0.01))
    }
  }
  n_free <- C * npc + C - 1L
  labels <- .par_labels(dgp)
  cmap <- .common_map(dgp, C, .lab("delta", C - 1L))
  if (!is.null(cmap)) n_free <- ncol(cmap$M)

  obj <- .lcmnl_objective(ctx, C)

  result <- tryCatch(
    .optim_mapped(obj, par0, cmap, maxit, reltol),
    error = function(e) {
      message("[LCMNL C=", C, "] optim error: ", conditionMessage(e)); NULL
    }
  )
  if (is.null(result) || !result$convergence %in% c(0L, 1L))
    return(fail_result)

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

  # Covariance of the full parameter vector (per-class betas + ASCs and the
  # C-1 share deltas). Skipped during the multistart screen (vcov = FALSE);
  # klue_lcmnl attaches it once, to the polished winner.
  cv <- if (vcov) .lcmnl_vcov(ctx, p, C, labels, cmap) else NULL
  if (!is.null(cv)) names(p) <- cv$pnames

  if (!is.null(cv) && anyNA(cv$vcov)) message(.singular_msg("LCMNL", C, cv$aliased))
  .fit_status(list(converged = result$convergence == 0L, singular = !is.null(cv) && anyNA(cv$vcov),
       aliased = cv$aliased, C = C, model_type = "LCMNL",
       LL = LL, BIC = BIC, AIC = -2 * LL + 2 * n_free,
       ICL = BIC + 2 * H, ICL_BIC = 2 * H,
       k = n_free, betas = betas_mat, class_probs = class_probs,
       posteriors = posteriors_c, par = p,
       vcov = if (is.null(cv)) NULL else cv$vcov,
       robust_vcov = if (is.null(cv)) NULL else cv$robust_vcov,
       common = if (is.null(cmap)) character(0) else cmap$common,
       price = paste0(labels[n_beta], "_class", seq_len(C))), N, H)
}

#' Concomitant LCMNL: covariate-driven class membership
#'
#' Direct-MLE LCMNL in which the class-membership probabilities depend on
#' respondent-level covariates, \code{pi_{nc} = softmax_c(t(Z_n) gamma_c)} with the
#' last class as reference (\code{gamma_C = 0}). The class-specific taste model
#' is identical to \code{\link{estimate_lcmnl}}; only the share model changes.
#' With no covariates this reduces exactly to the constant-share LCMNL.
#'
#' @param database Long-format choice data (as for \code{estimate_lcmnl}), plus
#'   the membership covariate columns, which must be constant within respondent.
#' @param C Integer number of latent classes.
#' @param membership Character vector of covariate column names entering the
#'   class-membership model (an intercept is always included). Default none
#'   (recovers the constant-share model).
#' @param start_betas,start_shares Optional taste and share starting values; if
#'   \code{NULL}, k-means starts are computed. Membership slopes start at zero.
#' @param dgp,maxit,reltol,start_par,vcov As in \code{\link{estimate_lcmnl}}.
#' @return A list like \code{estimate_lcmnl} with, additionally, \code{gamma}
#'   (a \code{p} x \code{(C-1)} matrix of membership coefficients, rows named by
#'   \code{Z}'s columns) and \code{shares} the \code{N} x \code{C} matrix of
#'   respondent-level class probabilities (\code{class_probs} is their mean).
#'   \code{vcov}/\code{robust_vcov} cover taste params then the membership
#'   parameters \code{{cov}_class{c}}.
#' @seealso \code{\link{estimate_lcmnl}}
#' @export
estimate_lcmnl_cov <- function(database, C, membership = character(0),
                               start_betas = NULL, start_shares = NULL,
                               dgp = DGP_DEFAULT, maxit = MAX_ITER,
                               reltol = 1e-10, start_par = NULL, vcov = TRUE) {
  ctx <- .lcmnl_context(database, dgp)
  N <- ctx$N; npc <- ctx$npc; n_beta <- ctx$n_beta; n_asc <- ctx$n_asc
  Z <- .lcmnl_membership_design(database, ctx, membership)
  p_mem <- ncol(Z); n_taste <- C * npc; n_free <- n_taste + (C - 1L) * p_mem
  fail_result <- .lcmnl_fail(C, N, n_beta)

  if (is.null(start_par)) {
    if (is.null(start_betas)) {
      starts <- klue_starts(database, C, "kmeans", dgp = dgp)
      start_betas  <- starts$betas
      start_shares <- starts$shares
    }
    if (is.null(start_shares)) start_shares <- rep(1 / C, C)
    par0 <- numeric(n_free)
    for (ci in 1:C) par0[(ci - 1L) * npc + 1:n_beta] <- start_betas[ci, ]
    if (C > 1L) for (ci in 1:(C - 1L))           # intercept from starting shares
      par0[n_taste + (ci - 1L) * p_mem + 1L] <-
        log(max(start_shares[ci], 0.01) / max(start_shares[C], 0.01))
  } else {
    par0 <- start_par
  }

  labels <- .par_labels(dgp)
  cmap <- .common_map(dgp, C, .lcmnl_cov_mem_names(C, Z))
  k_free <- if (is.null(cmap)) n_free else ncol(cmap$M)
  obj <- .lcmnl_cov_objective(ctx, C, Z)
  result <- tryCatch(
    .optim_mapped(obj, par0, cmap, maxit, reltol),
    error = function(e) {
      message("[LCMNL-cov C=", C, "] optim error: ", conditionMessage(e)); NULL
    })
  if (is.null(result) || !result$convergence %in% c(0L, 1L)) return(fail_result)

  p <- result$par; LL <- -result$value
  betas_mat <- t(vapply(1:C, function(ci) p[(ci - 1L) * npc + 1:n_beta],
                        numeric(n_beta)))
  gamma_mat <- if (C > 1L) matrix(p[(n_taste + 1L):n_free], nrow = p_mem,
                                  dimnames = list(colnames(Z), paste0("class", 1:(C - 1L))))
               else matrix(nrow = p_mem, ncol = 0)
  st <- obj$state(p); pimat <- st$pimat
  post <- obj$posterior(st)$w; post_c <- pmax(post, 1e-100)
  H <- -sum(post_c * log(post_c))
  BIC <- -2 * LL + k_free * log(N)

  cv <- if (vcov) .lcmnl_cov_vcov(ctx, C, Z, p, labels, cmap) else NULL
  if (!is.null(cv)) names(p) <- cv$pnames

  if (!is.null(cv) && anyNA(cv$vcov)) message(.singular_msg("LCMNL-cov", C, cv$aliased))
  .fit_status(list(converged = result$convergence == 0L, singular = !is.null(cv) && anyNA(cv$vcov),
       aliased = cv$aliased, C = C, model_type = "LCMNL-cov",
       LL = LL, BIC = BIC, AIC = -2 * LL + 2 * k_free,
       ICL = BIC + 2 * H, ICL_BIC = 2 * H, k = k_free,
       betas = betas_mat, gamma = gamma_mat,
       class_probs = colMeans(pimat), shares = pimat,
       posteriors = post_c, par = p,
       vcov = if (is.null(cv)) NULL else cv$vcov,
       robust_vcov = if (is.null(cv)) NULL else cv$robust_vcov,
       common = if (is.null(cmap)) character(0) else cmap$common,
       price = paste0(labels[n_beta], "_class", seq_len(C))), N, H)
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
  if (length(dgp$common))
    stop("estimate_lcmnl_em does not support common coefficients; use estimate_lcmnl")
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
      if (sum(rw) < 1e-6) {
        message("[EM C=", C, " iter ", it, "] class ", ci,
                " collapsed; keeping its current parameters")
        next
      }
      opt <- tryCatch(.fit_mnl(ctx, par_c[ci, ], rw = rw),
                      error = function(e) NULL)
      if (is.null(opt)) {
        message("[EM C=", C, " iter ", it, "] M-step fit for class ", ci,
                " errored; keeping its current parameters")
      } else par_c[ci, ] <- opt$par
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
       par = c(as.vector(t(par_c)), log(pi_c[-C] / pi_c[C])),
       em_iters = em_iters, estimator = "em")
}

#' Multi-start LCMNL: best of the six clustering initialisations
#'
#' For C = 1 (MNL) a single run suffices; for C >= 2 the model is estimated
#' from all six clustering starts and the best log-likelihood wins. With the
#' \code{"ml"} estimator the starts are screened at a loose tolerance
#' (\code{maxit = 200}, \code{reltol = 1e-6}; fits hitting the cap still
#' count) and only the winner is polished to full precision, warm-started
#' from its own solution. Set \code{screen = FALSE} (or session-wide
#' \code{options(klue.screen = FALSE)}) to fit every start at full precision
#' with no polish step -- the pre-0.9.1 code path, kept for exact
#' reproduction of results produced with it. Each start logs its LL and
#' elapsed time via \code{message()}.
#' @param database Data frame in long format, sorted by respondent in blocks of
#'   \code{T_per_n} rows, with columns \code{ID}, \code{CHOICE}, and the
#'   attribute columns named by the DGP (\code{x*_j}, \code{price_j}). Optional
#'   availability columns \code{av_1..av_J} (1 = available, 0 = not) mask
#'   unavailable alternatives out of each task's choice set; if absent or all
#'   ones, every alternative is available (the default). The chosen alternative
#'   must be available in every task.
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
#' @param screen logical; the 0.9.1 fast-screening behaviour (default). With
#'   the \code{"ml"} estimator the six starts are screened at the loose
#'   tolerance and only the winner is polished; the cluster-wise MNL start
#'   fits also run loose (see \code{\link{klue_starts}}). \code{FALSE}
#'   restores the pre-0.9.1 path end to end: tight cluster-MNL start fits,
#'   every start estimated at full precision, no polish. The default reads
#'   \code{getOption("klue.screen", TRUE)}, so
#'   \code{options(klue.screen = FALSE)} restores the old path in every
#'   caller (\code{klue()}, the study drivers, and the start generators)
#'   without changing call sites.
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
#' @param clean logical; fit every start at full precision with its covariance
#'   and return the highest-LL start whose fit has not \code{failed} (no
#'   polish step). If none qualifies, the highest-LL start is returned with
#'   \code{failed = TRUE}. The fit then also reports \code{n_same_optimum},
#'   the number of starts within 1e-4 of its LL. Default \code{FALSE} keeps
#'   the 0.10.0 path. Needs the \code{"ml"} estimator.
#' @param maxit,reltol BFGS settings of the full-precision fits.
#' @param cluster_cap maximum respondents clustered for the starts (see
#'   \code{\link{get_all_starts}}).
#' @export
klue_lcmnl <- function(database, C, dgp = DGP_DEFAULT,
                       estimator = c("ml", "em"),
                       feature_type = c("rp", "onehot"),
                       n_cores = 1L,
                       screen = getOption("klue.screen", TRUE),
                       clean = FALSE, maxit = MAX_ITER, reltol = 1e-10,
                       cluster_cap = 2000L) {
  estimator <- match.arg(estimator)
  feature_type <- match.arg(feature_type)
  if (estimator == "em" && isTRUE(clean)) stop("clean needs the ml estimator")
  fit_one <- if (estimator == "em") estimate_lcmnl_em else estimate_lcmnl
  if (C == 1L) {                # all six clustering starts coincide at C = 1
    res <- if (estimator == "ml") fit_one(database, 1L, dgp = dgp, maxit = maxit, reltol = reltol)
           else fit_one(database, 1L, dgp = dgp)
    res$best_method <- "pooled"
    res$method_results <- list(pooled = res)
    return(res)
  }
  all_starts <- get_all_starts(database, C, dgp = dgp,
                               feature_type = feature_type, screen = screen,
                               cluster_cap = cluster_cap)
  .lcmnl_multistart(database, C, dgp, all_starts, fit_one, estimator,
                    feature_type, n_cores, screen, isTRUE(clean), maxit, reltol)
}

# A fit that estimates and, when alpha_cost is set, whose price coefficient is
# significantly negative in every class (one-sided, robust SE): the selection
# rule's admissibility check.
.lc_ok <- function(fit, alpha_cost = NULL) {
  if (isTRUE(fit$failed) || !is.finite(fit$LL) || is.null(fit$robust_vcov)) return(FALSE)
  if (is.null(alpha_cost)) return(TRUE)
  pr <- fit$price; z <- fit$par[pr] / sqrt(diag(fit$robust_vcov)[pr])
  all(is.finite(z) & z < stats::qnorm(alpha_cost))
}

# The six-start fit of klue_lcmnl for one dgp (C >= 2).
.lcmnl_multistart <- function(database, C, dgp, all_starts, fit_one, estimator,
                              feature_type, n_cores, screen, clean, maxit, reltol) {
  N <- length(unique(database$ID))

  # The six per-start fits are independent and ~90% of the runtime; optionally
  # fit them concurrently. Order is preserved so the best-of-six tie-break is
  # identical to the sequential path. The ml screen runs loose: only the
  # winner gets the tight tolerance, so maxit-hit fits stay in the race.
  # screen = FALSE fits each start tight and skips the polish (pre-0.9.1 path).
  # Covariance is skipped on every candidate fit (vcov = FALSE) and attached
  # once below, to the winner, so the Hessian is formed only for the fit kept.
  # clean = TRUE fits every start tight with its covariance and keeps the best
  # start that has not failed.
  screen_ctl <- if (estimator == "ml") {
                  if (clean) list(maxit = maxit, reltol = reltol, vcov = TRUE)
                  else if (isTRUE(screen)) list(maxit = 200L, reltol = 1e-6, vcov = FALSE)
                  else list(maxit = maxit, reltol = reltol, vcov = FALSE)
                } else list()
  nm_ok <- names(all_starts)[!vapply(all_starts, is.null, logical(1))]
  fit_start <- function(nm) {
    t0 <- proc.time()[["elapsed"]]
    res <- tryCatch(
      do.call(fit_one, c(list(database, C, start_betas = all_starts[[nm]]$betas,
                              start_shares = all_starts[[nm]]$shares,
                              dgp = dgp), screen_ctl)),
      error = function(e) NULL)
    message(sprintf("  [klue C=%d %s/%s] LL = %s (%.1f min)", C, feature_type,
                    nm, if (is.null(res)) "error" else sprintf("%.1f", res$LL),
                    (proc.time()[["elapsed"]] - t0) / 60))
    res
  }
  fits <- if (n_cores > 1L && length(nm_ok) > 1L) {
    invisible(gc(FALSE))   # shrink the heap the forked children inherit (COW)
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
    if (!is.list(res)) {   # NULL = fit errored (logged in the worker); other
      if (!is.null(res))   # non-lists are mclapply try-errors from dead workers
        warning("klue_lcmnl: start '", nm, "' lost to a worker failure: ",
                paste(as.character(res), collapse = " "), call. = FALSE)
      next
    }
    if (!is.finite(res$LL)) next
    method_results[[nm]] <- res
    if (!clean && res$LL > best$LL) {
      best <- res
      best$best_method <- nm
    }
  }
  if (clean && length(method_results)) {
    LLs <- vapply(method_results, `[[`, 0, "LL")
    ok <- !vapply(method_results, function(r) isTRUE(r$failed), NA)
    pick <- names(method_results)[if (any(ok)) which(ok)[which.max(LLs[ok])] else which.max(LLs)]
    best <- method_results[[pick]]
    best$best_method <- pick
    best$n_same_optimum <- sum(abs(LLs - best$LL) < 1e-4)
  }
  if (!clean && estimator == "ml" && !is.na(best$best_method) && is.finite(best$LL)) {
    if (isTRUE(screen)) {
      # Polish the loose winner to full precision (this fit also carries vcov).
      polished <- tryCatch(
        estimate_lcmnl(database, C, dgp = dgp, start_par = best$par,
                       maxit = maxit, reltol = reltol),
        error = function(e) NULL)
      if (!is.null(polished) && polished$converged) {
        polished$best_method <- best$best_method
        best <- polished
      }
    }
    if (is.null(best$vcov)) {
      # No-screen path (winner already tight), or polish failed: attach the
      # covariance at the winner's parameters without re-optimising, and
      # refresh the status fields that depend on it.
      cv <- tryCatch(.lcmnl_vcov(.lcmnl_context(database, dgp), best$par, C,
                                 .par_labels(dgp), .common_map(dgp, C, .lab("delta", C - 1L))),
                     error = function(e) NULL)
      if (!is.null(cv)) {
        names(best$par)    <- cv$pnames
        best$vcov          <- cv$vcov
        best$robust_vcov   <- cv$robust_vcov
        best$singular      <- anyNA(cv$vcov)
        best$aliased       <- cv$aliased
        best <- .fit_status(best, N, best$ICL_BIC / 2)
      }
    }
  }
  best$method_results <- method_results
  best
}
