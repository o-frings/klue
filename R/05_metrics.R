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
  # Optimal label permutation by exhaustive search (fine for K <= 8)
  perms <- as.matrix(expand.grid(rep(list(1:K), K)))
  perms <- perms[apply(perms, 1, function(p) length(unique(p)) == K), , drop = FALSE]
  costs <- apply(perms, 1, function(p) sum((true_betas - est_betas[p, , drop = FALSE])^2))
  best_perm <- perms[which.min(costs), ]
  diffs <- true_betas - est_betas[best_perm, , drop = FALSE]
  list(rmse = sqrt(mean(diffs^2)), bias = mean(diffs))
}
