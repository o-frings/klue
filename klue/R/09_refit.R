# Re-estimating a fitted model on new data (0.10.0): subsets, one arm of an
# experiment, added or dropped terms, a membership model, simulated replications.

#' Re-estimate a fitted MNL or latent class model
#'
#' Fits the same model class to \code{database}, starting from the fit's
#' estimates matched by label: an MNL, or a constant-share LCMNL at
#' \code{fit$C} (with the common coefficients of \code{dgp}), or with
#' \code{membership} the concomitant model at \code{fit$C}. Labels new to
#' \code{dgp} start at 0, labels it drops are ignored, and \code{start}
#' overrides any start value. A common coefficient starts at the mean of its
#' class values.
#' @param fit An MNL or constant-share LCMNL klue fit.
#' @param database Canonical database.
#' @param dgp The dgp of the model to fit (\code{\link{klue_dgp}} with labels).
#' @param start Optional named start values: full parameter names, or labels
#'   (set in every class).
#' @param multistart Latent class only: also run
#'   \code{\link{klue_lcmnl}(clean = TRUE)} and keep it if it has not failed
#'   and has the higher log-likelihood.
#' @param membership Optional respondent-level covariate columns: fit the
#'   concomitant model (\code{\link{estimate_lcmnl_cov}}); the intercepts
#'   start at the fit's share parameters and the slopes at 0.
#' @param maxit,reltol BFGS settings.
#' @param n_cores Cores of the multistart.
#' @return The fit, with \code{failed}.
#' @export
klue_refit <- function(fit, database, dgp, start = NULL, multistart = FALSE, membership = NULL,
                       maxit = MAX_ITER, reltol = 1e-10, n_cores = 1L) {
  if (!fit$model_type %in% c("MNL", "LCMNL")) stop("klue_refit takes an MNL or constant-share LCMNL fit")
  C <- as.integer(fit$C)
  one <- function(d) {
    s0 <- .refit_start(fit, d, C, start)
    if (!length(membership))
      return(estimate_lcmnl(database, C, dgp = d, start_par = s0, maxit = maxit, reltol = reltol))
    npt <- C * d$npc
    mem <- unlist(lapply(seq_len(C - 1L), function(j) c(s0[npt + j], numeric(length(membership)))))
    estimate_lcmnl_cov(database, C, membership, dgp = d, start_par = c(s0[seq_len(npt)], mem),
                       maxit = maxit, reltol = reltol)
  }
  f <- one(dgp)
  if (isTRUE(multistart) && C >= 2L && !length(membership)) {
    m <- klue_lcmnl(database, C, dgp = dgp, clean = TRUE, maxit = maxit, reltol = reltol, n_cores = n_cores)
    if (.lc_ok(m) && (!.lc_ok(f) || m$LL > f$LL + 1e-8)) f <- m
  }
  f
}

# Start vector in the layout of `dgp` with C classes, from a fit's estimates by label.
.refit_start <- function(fit, dgp, C, start = NULL) {
  l <- .par_labels(dgp); full <- .full_pnames(l, C); old <- fit$par
  s <- stats::setNames(numeric(length(full)), full)
  if (fit$C == 1L) {
    for (lab in intersect(l, names(old)))
      s[if (C == 1L) lab else paste0(lab, "_class", seq_len(C))] <- old[[lab]]
  } else {
    hit <- intersect(full, names(old)); s[hit] <- old[hit]
  }
  for (k in names(start)) {
    tgt <- if (k %in% full) k else if (C > 1L && k %in% l) paste0(k, "_class", seq_len(C))
           else stop("unknown start value: ", k)
    s[tgt] <- start[[k]]
  }
  unname(s)
}
