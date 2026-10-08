# Post-estimation inference on a klue fit (estimate_lcmnl / estimate_lcmnl_cov).
# These operate on the fit's parameter vector and covariance (robust by default),
# so downstream code never hand-rolls Wald / Delta / Fieller / TOST on klue
# output. Parameters are referenced by the names on `fit$par` (and the dimnames
# of the covariance): e.g. "b1".."b{n_beta}", "asc1".., or for C >= 2 the
# per-class "b{k}_class{c}" / concomitant "{cov}_class{c}" labels (the labels of
# klue_dgp(labels = ) when set; klue_par() resolves a label). Nonlinear functions of
# the parameter vector go through the delta method (klue_nlcom) with a central-
# difference Jacobian; klue_mean_fn and klue_wtp_fn build the population means of
# a latent class fit as such functions.

# Pull the estimate vector and a usable covariance (robust unless robust=FALSE).
.klue_estV <- function(fit, robust = TRUE) {
  est <- fit$par
  if (is.null(est)) stop("fit has no $par (is this a klue fit?)")
  if (is.null(names(est))) stop("fit$par is unnamed; refit with vcov = TRUE for labelled parameters")
  V <- if (isTRUE(robust) && !is.null(fit$robust_vcov)) fit$robust_vcov else fit$vcov
  if (is.null(V) || anyNA(V)) stop("fit has no usable covariance (refit with vcov = TRUE)")
  list(est = est, V = V)
}

# Expand a weight spec into a full-length weight vector over all parameters.
# `w` is a named numeric vector (parameter -> weight) or a character vector of
# parameter names (unit weights).
.klue_wvec <- function(w, nm) {
  if (is.character(w)) w <- stats::setNames(rep(1, length(w)), w)
  if (is.null(names(w))) stop("weights must be named by parameter")
  bad <- setdiff(names(w), nm)
  if (length(bad)) stop("unknown parameter(s): ", paste(bad, collapse = ", "))
  full <- stats::setNames(numeric(length(nm)), nm); full[names(w)] <- as.numeric(w); full
}

#' Linear combination of coefficients from a klue fit
#'
#' @param fit A klue fit (from \code{\link{estimate_lcmnl}} /
#'   \code{\link{estimate_lcmnl_cov}}) carrying \code{par} and a covariance.
#' @param weights Named numeric vector \code{c(param = weight, ...)} (or a
#'   character vector of parameter names for unit weights) defining \eqn{c'\beta}.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @return Named vector: \code{estimate}, \code{se}, \code{z}, two-sided \code{p}.
#' @export
klue_lincom <- function(fit, weights, robust = TRUE) {
  cv <- .klue_estV(fit, robust); w <- .klue_wvec(weights, names(cv$est))
  val <- sum(w * cv$est); se <- sqrt(as.numeric(t(w) %*% cv$V %*% w))
  c(estimate = val, se = se, z = val / se, p = 2 * stats::pnorm(-abs(val / se)))
}

#' Joint Wald test that a set of coefficients equals a value
#'
#' @param fit A klue fit.
#' @param terms Character vector of parameter names to test jointly.
#' @param value Null value(s) (recycled). Default \code{0}.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @param fn Optional function of the parameter vector returning the vector
#'   tested against \code{value} (delta method). The statistic then uses a
#'   generalised inverse of its covariance, and \code{df} is that
#'   covariance's rank (eigenvalues above \code{tol} times the largest).
#'   \code{terms} is ignored.
#' @param tol Relative eigenvalue tolerance of the rank (with \code{fn}).
#' @return Named vector: \code{chisq}, \code{df}, \code{p}.
#' @export
klue_wald <- function(fit, terms = NULL, value = 0, robust = TRUE, fn = NULL, tol = 1e-8) {
  if (!is.null(fn)) {
    cv <- .klue_estV(fit, robust); g <- fn(cv$est) - value; J <- .klue_jacobian(fn, cv$est)
    S <- J %*% cv$V %*% t(J); e <- eigen((S + t(S)) / 2, symmetric = TRUE)
    keep <- e$values > tol * max(e$values)
    W <- sum(drop(crossprod(e$vectors[, keep, drop = FALSE], g))^2 / e$values[keep]); df <- sum(keep)
    return(c(chisq = W, df = df, p = stats::pchisq(W, df, lower.tail = FALSE)))
  }
  cv <- .klue_estV(fit, robust); i <- match(terms, names(cv$est))
  if (anyNA(i)) stop("unknown parameter(s): ", paste(terms[is.na(i)], collapse = ", "))
  d <- cv$est[i] - value
  W <- as.numeric(t(d) %*% solve(cv$V[i, i, drop = FALSE]) %*% d)
  c(chisq = W, df = length(i), p = stats::pchisq(W, length(i), lower.tail = FALSE))
}

#' Delta-method willingness-to-pay (ratio of linear combinations)
#'
#' Returns \eqn{-\,\mathrm{scale}\times (a'\beta)/(b'\beta)} with a Delta-method
#' standard error and 95\% interval. For a plain WTP pass \code{num = c(attr = 1)},
#' \code{den = c(cost = 1)}.
#' @param fit A klue fit.
#' @param num,den Named weight vectors for the numerator and denominator.
#' @param scale Multiplier (e.g. currency scaling). Default \code{1}.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @return Named vector: \code{wtp}, \code{se}, \code{lo}, \code{hi}.
#' @export
klue_wtp <- function(fit, num, den, scale = 1, robust = TRUE) {
  cv <- .klue_estV(fit, robust)
  a <- .klue_wvec(num, names(cv$est)); b <- .klue_wvec(den, names(cv$est))
  A <- sum(a * cv$est); B <- sum(b * cv$est); wtp <- -scale * A / B
  grad <- scale * (-a / B + b * A / B^2)                 # d(-A/B)/dbeta
  use <- which(grad != 0)
  se <- sqrt(as.numeric(t(grad[use]) %*% cv$V[use, use] %*% grad[use]))
  c(wtp = wtp, se = se, lo = wtp - stats::qnorm(.975) * se, hi = wtp + stats::qnorm(.975) * se)
}

#' Fieller confidence interval for a ratio of linear combinations
#'
#' 95\%-by-default CI for \eqn{(a'\beta)/(b'\beta)}. Returns \code{NA} bounds when
#' the denominator is not significantly different from zero (unbounded CI).
#' @param fit A klue fit.
#' @param num,den Named weight vectors for the numerator and denominator, or
#'   functions of the parameter vector (a generalised Fieller interval from
#'   their delta-method covariance), e.g. from \code{\link{klue_mean_fn}}.
#' @param level Confidence level. Default \code{0.95}.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @return Named vector: \code{ratio}, \code{lo}, \code{hi}.
#' @export
klue_fieller <- function(fit, num, den, level = 0.95, robust = TRUE) {
  cv <- .klue_estV(fit, robust)
  if (is.function(num) || is.function(den)) {
    lin <- function(w) { a <- .klue_wvec(w, names(cv$est)); function(p) sum(a * p) }
    fa <- if (is.function(num)) num else lin(num); fb <- if (is.function(den)) den else lin(den)
    A <- fa(cv$est); B <- fb(cv$est)
    J <- .klue_jacobian(function(p) c(fa(p), fb(p)), cv$est); S <- J %*% cv$V %*% t(J)
    vAA <- S[1, 1]; vBB <- S[2, 2]; vAB <- S[1, 2]
  } else {
    a <- .klue_wvec(num, names(cv$est)); b <- .klue_wvec(den, names(cv$est))
    A <- sum(a * cv$est); B <- sum(b * cv$est)
    vAA <- as.numeric(t(a) %*% cv$V %*% a); vBB <- as.numeric(t(b) %*% cv$V %*% b)
    vAB <- as.numeric(t(a) %*% cv$V %*% b)
  }
  z2 <- stats::qnorm(1 - (1 - level) / 2)^2
  qa <- B^2 - z2 * vBB; qb <- -2 * (A * B - z2 * vAB); qc <- A^2 - z2 * vAA
  if (qa <= 0) return(c(ratio = A / B, lo = NA, hi = NA))   # denominator ~ 0: unbounded
  disc <- qb^2 - 4 * qa * qc
  if (disc < 0) return(c(ratio = A / B, lo = NA, hi = NA))
  r <- sort(c((-qb - sqrt(disc)) / (2 * qa), (-qb + sqrt(disc)) / (2 * qa)))
  c(ratio = A / B, lo = r[1], hi = r[2])
}

#' Two-one-sided-tests (TOST) equivalence test for a contrast
#'
#' Tests whether \eqn{c'\beta} lies within \eqn{\pm}\code{margin}. A small
#' \code{p_tost} rejects non-equivalence (i.e. supports equivalence).
#' @param fit A klue fit.
#' @param weights Named weight vector defining the contrast \eqn{c'\beta}, or a
#'   function of the parameter vector (delta method, \code{\link{klue_nlcom}}).
#' @param margin Positive equivalence margin.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @return Named vector: \code{estimate}, \code{se}, \code{margin}, \code{p_tost}.
#' @export
klue_tost <- function(fit, weights, margin, robust = TRUE) {
  lc <- if (is.function(weights)) klue_nlcom(fit, weights, robust) else klue_lincom(fit, weights, robust)
  est <- unname(lc["estimate"]); se <- unname(lc["se"])
  p_lo <- stats::pnorm((est + margin) / se, lower.tail = FALSE)   # H0: c'b <= -margin
  p_hi <- stats::pnorm((est - margin) / se)                       # H0: c'b >=  margin
  c(estimate = est, se = se, margin = margin, p_tost = max(p_lo, p_hi))
}

#' Parameter name of a labelled coefficient in a klue fit
#'
#' Resolves a coefficient label (see the \code{labels} of \code{\link{klue_dgp}})
#' to the name it carries in \code{fit$par}: the label itself in an MNL, and
#' \code{<label>_class<c>} in a latent class fit. A coefficient common to all
#' classes resolves to its class-1 entry when \code{class} is not given (every
#' class holds the same value). Use the result as weight names in
#' \code{\link{klue_lincom}}, \code{\link{klue_wtp}} and the other tests.
#' @param fit A klue fit.
#' @param label Coefficient label.
#' @param class Optional class number(s) for a latent class fit.
#' @return Character vector of parameter names.
#' @export
klue_par <- function(fit, label, class = NULL) {
  nm <- names(fit$par)
  if (is.null(class) && label %in% nm) return(label)
  if (is.null(class)) {
    if (!label %in% fit$common) stop("give the class of class-specific coefficient '", label, "'")
    class <- 1L
  }
  out <- paste0(label, "_class", class)
  if (!all(out %in% nm)) stop("unknown parameter(s): ", paste(setdiff(out, nm), collapse = ", "))
  out
}

# Central-difference Jacobian (rows = outputs of fn) with step 1e-6 * max(|x_i|, 1).
.klue_jacobian <- function(fn, x) {
  m <- length(fn(x)); J <- matrix(0, m, length(x))
  for (i in seq_along(x)) {
    h <- 1e-6 * max(abs(x[i]), 1); e <- replace(numeric(length(x)), i, h)
    J[, i] <- (fn(x + e) - fn(x - e)) / (2 * h)
  }
  J
}

#' Delta-method inference on a function of a klue fit's parameters
#'
#' @param fit A klue fit carrying \code{par} and a covariance.
#' @param fn Function of the full parameter vector returning one number, e.g.
#'   from \code{\link{klue_wtp_fn}} or \code{\link{klue_mean_fn}}.
#' @param robust Use the cluster-robust covariance (default \code{TRUE}).
#' @param level Confidence level of \code{lo} and \code{hi}. Default \code{0.95}.
#' @return Named vector: \code{estimate}, \code{se}, \code{z}, two-sided
#'   \code{p}, \code{lo}, \code{hi}. The gradient is a central difference with
#'   step \code{1e-6 * max(|par_i|, 1)}.
#' @export
klue_nlcom <- function(fit, fn, robust = TRUE, level = 0.95) {
  cv <- .klue_estV(fit, robust); est <- fn(cv$est)
  if (length(est) != 1L) stop("fn must return one number; use klue_wald(fn = ) for a vector")
  J <- .klue_jacobian(fn, cv$est); se <- sqrt(drop(J %*% cv$V %*% t(J)))
  q <- stats::qnorm(1 - (1 - level) / 2)
  c(estimate = est, se = se, z = est / se, p = 2 * stats::pnorm(-abs(est / se)),
    lo = est - q * se, hi = est + q * se)
}

# Class layout of an MNL or constant-share LCMNL fit from its parameter names.
.fit_layout <- function(fit) {
  nm <- names(fit$par)
  if (is.null(nm)) stop("fit$par is unnamed; refit with vcov = TRUE")
  if (!fit$model_type %in% c("MNL", "LCMNL"))
    stop("population means are implemented for MNL and constant-share LCMNL fits")
  list(nm = nm, C = as.integer(fit$C), delta = grep("^delta[0-9]+$", nm))
}
.klue_idx <- function(L, labels, class) {
  key <- if (L$C == 1L) labels else paste0(labels, "_class", class)
  i <- match(key, L$nm)
  if (anyNA(i)) stop("unknown parameter(s): ", paste(key[is.na(i)], collapse = ", "))
  i
}
.klue_shares <- function(p, delta) { d <- c(unname(p[delta]), 0); e <- exp(d - max(d)); e / sum(e) }

#' Population mean of a coefficient as a function of the parameters
#'
#' Returns \code{function(par)} giving the share-weighted mean of a labelled
#' coefficient over the latent classes, \eqn{\sum_c \pi_c \beta_c} (the
#' coefficient itself in an MNL). Pass it to \code{\link{klue_nlcom}}, or as
#' numerator and denominator to \code{\link{klue_fieller}} for a ratio of means.
#' @param fit An MNL or constant-share LCMNL klue fit.
#' @param label Coefficient label.
#' @export
klue_mean_fn <- function(fit, label) {
  L <- .fit_layout(fit); i <- .klue_idx(L, label, seq_len(L$C))
  if (L$C == 1L) function(p) p[[i]] else function(p) sum(.klue_shares(p, L$delta) * p[i])
}

#' Willingness to pay as a function of the parameters
#'
#' Returns \code{function(par)} giving \eqn{-\mathrm{scale}\, a'\beta / \beta_{den}}:
#' in an MNL the ratio, in an LCMNL its share-weighted mean over the classes
#' (the population mean), or one class's ratio with \code{class}. Common
#' coefficients enter every class with the same value.
#' @param fit An MNL or constant-share LCMNL klue fit.
#' @param num Named weights over coefficient labels (constants included).
#' @param den Label of the denominator; default the fit's price coefficient.
#' @param scale Multiplier (e.g. currency units). Default \code{1}.
#' @param class Optional class number for that class's ratio.
#' @export
klue_wtp_fn <- function(fit, num, den = NULL, scale = 1, class = NULL) {
  L <- .fit_layout(fit)
  if (is.null(den)) den <- sub("_class1$", "", fit$price[1])
  cls <- if (L$C == 1L) 1L else seq_len(L$C)
  ia <- lapply(cls, function(c) .klue_idx(L, names(num), c)); id <- .klue_idx(L, den, cls)
  a <- as.numeric(num)
  ratio <- function(p) vapply(seq_along(cls), function(k) -scale * sum(a * p[ia[[k]]]) / p[id[k]], 0)
  if (L$C == 1L) return(function(p) ratio(p))
  if (!is.null(class)) return(function(p) ratio(p)[class])
  function(p) sum(.klue_shares(p, L$delta) * ratio(p))
}
