# Recovery metrics: adjusted Rand index and label-permutation-invariant
# RMSE/bias of class coefficients.

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
  K <- nrow(true_betas)
  if (K != nrow(est_betas)) return(list(rmse = NA, bias = NA))
  if (K == 1) {
    diffs <- true_betas - est_betas
    return(list(rmse = sqrt(mean(diffs^2)), bias = mean(diffs)))
  }
  # Optimal label permutation by exhaustive search: the K^K enumeration
  # below explodes past K = 8 (9^9 = 387M rows), so refuse loudly.
  if (K > 8) stop("compute_recovery: exhaustive permutation matching supports ",
                  "K <= 8; got K = ", K)
  perms <- as.matrix(expand.grid(rep(list(1:K), K)))
  perms <- perms[apply(perms, 1, function(p) length(unique(p)) == K), , drop = FALSE]
  costs <- apply(perms, 1, function(p) sum((true_betas - est_betas[p, , drop = FALSE])^2))
  best_perm <- perms[which.min(costs), ]
  diffs <- true_betas - est_betas[best_perm, , drop = FALSE]
  list(rmse = sqrt(mean(diffs^2)), bias = mean(diffs))
}

#' D-error of a choice design
#'
#' MNL D-error \eqn{\det(I(\beta))^{-1/K}} of the task rows of a canonical
#' database (one row per choice task; \code{ID} and \code{CHOICE} are not
#' used), with the constants among the \code{K = npc} parameters (Ngene's
#' \code{;con}). The information sums over the rows, so the value matches
#' Ngene's D-error for one respondent answering every row.
#' @param database Canonical database of the design's tasks.
#' @param dgp Its dgp (\code{\link{klue_dgp}}).
#' @param par Optional prior vector (taste coefficients, then constants);
#'   default zeros.
#' @param draws Optional matrix of prior draws (one per row); the D-error is
#'   then their mean (Bayesian D-error).
#' @param tol Relative eigenvalue tolerance of the rank.
#' @return A list: \code{derror}, \code{rank} and \code{full_rank} of the
#'   information at the first prior, and \code{K}.
#' @export
klue_derror <- function(database, dgp, par = NULL, draws = NULL, tol = 1e-10) {
  Xf <- .design_full(database, dgp); n <- nrow(database); K <- dgp$npc
  info <- function(b) {
    P <- .row_probs(Xf, matrix(b, n, K, byrow = TRUE))
    M <- Reduce(`+`, Map(function(Xj, j) P[, j] * Xj, Xf, seq_along(Xf)))
    Reduce(`+`, lapply(seq_along(Xf), function(j) { Z <- Xf[[j]] - M; crossprod(Z * P[, j], Z) }))
  }
  B <- if (!is.null(draws)) as.matrix(draws) else matrix(if (is.null(par)) 0 else par, 1, K)
  ev <- eigen(info(B[1, ]), symmetric = TRUE, only.values = TRUE)$values
  rank <- sum(ev > tol * max(ev))
  d <- apply(B, 1, function(b) { dt <- determinant(info(b)); if (dt$sign <= 0) Inf else exp(-dt$modulus[[1]] / K) })
  list(derror = mean(d), rank = rank, full_rank = rank == K, K = K)
}
