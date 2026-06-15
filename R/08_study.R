# Monte Carlo study drivers reproducing Frings (2026). All drivers share the
# same skeleton: a condition grid with a deterministic per-condition seed,
# a run-one-condition worker dispatched over .klue_cores() forks, and a
# results frame. Per-condition seeds and grids are identical to klue 0.6.x.

# ---- shared plumbing --------------------------------------------------------

.study_mclapply <- function(n, run_one, n_cores = .klue_cores()) {
  parallel::mclapply(seq_len(n), run_one, mc.cores = n_cores)
}

# The standard seed formula used by the H1-style studies.
.cond_seed <- function(tK, kap, sig, rp, base = 0L) {
  as.integer(base + 1000 * tK + 100 * (kap * 100) + 10 * (sig * 100) + rp)
}

# The H1 ablation grid (convergence / initialisation / estimator studies).
.h1_conditions <- function(n_cond = NULL, cond_idx = NULL) {
  conds <- expand.grid(true_K = c(3, 4, 5), kappa = c(0.50, 0.75, 1.00),
                       sigma = c(0.15, 0.25), rep = 1:2)
  if (!is.null(cond_idx)) conds[cond_idx, , drop = FALSE]
  else conds[1:min(n_cond, nrow(conds)), ]
}

# BIC for klue_lcmnl over a range of C; Inf where not converged.
.bic_scan <- function(database, C_range, dgp) {
  vapply(C_range, function(Cc) {
    m <- klue_lcmnl(database, Cc, dgp = dgp)
    if (m$converged) m$BIC else Inf
  }, numeric(1))
}

# Diffuse random start spec r for condition seed `seed`: N(0,2) on generic
# attributes, -exp(N(0,1)) on price.
.random_start <- function(seed, r, tK, dgp) {
  set.seed(seed * 1000L + r)
  rb <- matrix(rnorm(tK * dgp$n_generic, mean = 0, sd = 2),
               nrow = tK, ncol = dgp$n_generic)
  list(betas = cbind(rb, -exp(rnorm(tK, mean = 0, sd = 1))),
       shares = rep(1 / tK, tK))
}

.col <- function(res_list, name) sapply(res_list, `[[`, name)

# ---- main simulation --------------------------------------------------------

#' Main Monte Carlo: class-count recovery across the condition grid
#'
#' Runs the primary arm: simulate data over a grid of true class counts,
#' segment separation, and within-class heterogeneity, then estimate LCMNL for
#' each candidate class count and record whether BIC, AIC, and ICL select the
#' true class count, plus parameter and clustering recovery.
#'
#' @param true_K_values true class counts to simulate.
#' @param kappa_values segment-separation values (kappa).
#' @param sigma_values within-class heterogeneity values (sigma).
#' @param n_reps replications per condition.
#' @param C_cands candidate class counts fitted in each condition.
#' @param dgp data-generating-process specification.
#' @param sep_profile optional separation profile passed to the simulator.
#' @param blocked if \code{TRUE}, use one shared blocked D-efficient design with
#'   T per respondent equal to the block size; if \code{FALSE}, use the
#'   per-respondent random-attribute path with T equal to \code{T_random}.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param seg_scale optional segment-membership scaling passed to the simulator.
#' @param T_random number of tasks per respondent in the random-attribute path.
#' @param verbose if \code{TRUE}, print per-condition progress.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per condition: condition index, true_K,
#'   kappa, sigma, rep, the BIC/AIC/ICL selected class counts and correctness
#'   flags, the best clustering method, and ARI, RMSE, bias, and ICL_BIC values.
#' @export
klue_study_main <- function(true_K_values = c(1, 2, 3, 4, 5),
                            kappa_values  = c(0.5, 0.75, 1.0, 1.25, 1.5),
                            sigma_values  = c(0.1, 0.15, 0.2, 0.25),
                            n_reps = 5, C_cands = 1:6,
                            dgp = DGP_DEFAULT, sep_profile = NULL,
                            blocked = TRUE, n_cards = 48L, n_blocks = 4L,
                            seg_scale = NULL, T_random = 20L,
                            verbose = TRUE, n_cores = .klue_cores()) {
  # Realistic-DCE baseline (blocked = TRUE, default): ONE blocked D-efficient
  # design shared across the study; T per respondent = block size. blocked =
  # FALSE reproduces the per-respondent random-attribute path (T = T_random).
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  # K=1: kappa irrelevant (no segments), fixed at 0 to avoid redundant runs.
  if (1L %in% true_K_values) {
    conditions <- rbind(
      expand.grid(true_K = 1, kappa = 0, sigma = sigma_values, rep = 1:n_reps),
      expand.grid(true_K = setdiff(true_K_values, 1), kappa = kappa_values,
                  sigma = sigma_values, rep = 1:n_reps))
  } else {
    conditions <- expand.grid(true_K = true_K_values, kappa = kappa_values,
                              sigma = sigma_values, rep = 1:n_reps)
  }
  nc <- nrow(conditions)
  if (verbose) {
    cat(sprintf("==== MAIN SIMULATION (%d conditions, %d cores, %d attrs) ====\n",
                nc, n_cores, dgp$n_beta))
  }

  .run_one <- function(i) {
    tK  <- conditions$true_K[i]; kap <- conditions$kappa[i]
    sig <- conditions$sigma[i];  rp  <- conditions$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    npc <- if (tK == 1L) 300L else 150L
    data <- klue_simulate(N_per_class = npc, T_tasks = T_random, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp, sep_profile = sep_profile,
                          design = design, seg_scale = seg_scale)

    bics <- aics <- icls <- rep(Inf, length(C_cands))
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

  # Batched dispatch so progress prints as conditions finish.
  batch_size <- n_cores * 3L
  res_list <- vector("list", nc)
  for (b_start in seq(1L, nc, by = batch_size)) {
    idx <- b_start:min(b_start + batch_size - 1L, nc)
    res_list[idx] <- parallel::mclapply(idx, .run_one, mc.cores = n_cores)
    if (verbose) {
      for (i in idx) {
        r <- res_list[[i]]
        c_label <- if (r$selected_bic == 1L) "C=1(MNL)" else paste0("C=", r$selected_bic)
        cat(sprintf("  [%3d/%d] K=%d kappa=%.2f sigma=%.2f rep=%d ... BIC->%s %s [%s]\n",
                    i, nc, r$true_K, r$kappa, r$sigma, r$rep, c_label,
                    ifelse(r$bic_correct == 1L, "OK", "X"),
                    ifelse(is.na(r$best_method), "?", r$best_method)))
      }
    }
  }

  data.frame(
    condition    = 1:nc,
    true_K       = .col(res_list, "true_K"),
    kappa        = .col(res_list, "kappa"),
    sigma        = .col(res_list, "sigma"),
    rep          = .col(res_list, "rep"),
    selected_bic = .col(res_list, "selected_bic"),
    selected_aic = .col(res_list, "selected_aic"),
    selected_icl = .col(res_list, "selected_icl"),
    bic_correct  = .col(res_list, "bic_correct"),
    aic_correct  = .col(res_list, "aic_correct"),
    icl_correct  = .col(res_list, "icl_correct"),
    best_method  = .col(res_list, "best_method"),
    ari          = .col(res_list, "ari"),
    rmse         = .col(res_list, "rmse"),
    bias         = .col(res_list, "bias"),
    icl_bic      = .col(res_list, "icl_bic"),
    stringsAsFactors = FALSE
  )
}

#' Print a summary of main Monte Carlo results
#'
#' Prints overall and per-condition BIC selection accuracy (by class count and
#' by separation) along with mean parameter recovery from a
#' \code{klue_study_main} results frame.
#'
#' @param df a results data.frame returned by \code{klue_study_main}.
#' @return Invisibly returns \code{NULL}; called for the printed summary.
#' @export
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

# ---- MNL vs LCMNL vs MMNL ----------------------------------------------------

#' Monte Carlo: MNL vs LCMNL vs MMNL model selection
#'
#' Compares fixed-coefficient MNL, latent-class MNL, and mixed (random
#' coefficient) MNL on each condition, recording which model BIC prefers. Covers
#' both a continuous-heterogeneity DGP (K=1) and discrete-segment DGPs (K>1).
#'
#' @param n_cond number of conditions to run from the condition grid.
#' @param n_draws number of draws for MMNL estimation.
#' @param C_cands candidate class counts fitted for the LCMNL arm.
#' @param verbose if \code{TRUE}, print per-condition progress.
#' @param dgp data-generating-process specification.
#' @param n_cores cores for the parallel LCMNL phase.
#' @return A data.frame with one row per converged condition: true_K, kappa,
#'   sigma, the MNL/LCMNL/MMNL BIC values, the selected LCMNL class count and
#'   method, and which model BIC prefers.
#' @export
klue_study_mmnl <- function(n_cond = 80, n_draws = N_DRAWS_MMNL,
                            C_cands = 1:5, verbose = TRUE, dgp = DGP_DEFAULT,
                            n_cores = .klue_cores()) {
  # K=1: pure continuous heterogeneity; sigma=0.35 included where the MMNL
  # advantage is clearest.
  conds <- rbind(
    expand.grid(true_K = 1, kappa = 0, sigma = c(0.15, 0.25, 0.35), rep = 1:3),
    expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                sigma = c(0.15, 0.25), rep = 1:3))
  conds <- conds[1:min(n_cond, nrow(conds)), ]
  nc <- nrow(conds)
  if (verbose) cat(sprintf("\n==== MNL vs LCMNL vs MMNL (%d conditions, incl. K=1) ====\n", nc))

  # Phase 1: LCMNL in parallel (no Apollo global state).
  if (verbose) cat("  Phase 1: LCMNL estimation (parallel)...\n")
  .run_lcmnl_part <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    npc <- if (tK == 1L) 300L else 150L
    data <- klue_simulate(N_per_class = npc, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    mnl_bic <- Inf
    best_lc_bic <- Inf; best_lc_C <- NA_integer_; best_lc_method <- NA_character_
    for (Cc in C_cands) {
      m <- klue_lcmnl(data$database, Cc, dgp = dgp)
      if (m$converged) {
        if (Cc == 1L) mnl_bic <- m$BIC
        if (Cc >= 2L && m$BIC < best_lc_bic) {
          best_lc_bic <- m$BIC; best_lc_C <- Cc; best_lc_method <- m$best_method
        }
      }
    }
    list(tK = tK, kap = kap, sig = sig, seed = seed, npc = npc,
         mnl_bic = mnl_bic, lc_bic = best_lc_bic,
         lc_C = best_lc_C, lc_method = best_lc_method)
  }
  lc_results <- .study_mclapply(nc, .run_lcmnl_part, n_cores)

  # Phase 2: MMNL sequentially (Apollo global state -- must be serial).
  if (verbose) cat("  Phase 2: MMNL estimation (sequential, Apollo)...\n")
  rows <- vector("list", nc)
  for (i in 1:nc) {
    lr <- lc_results[[i]]
    if (verbose) cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... ",
                             i, nc, lr$tK, lr$kap, lr$sig))
    data <- klue_simulate(N_per_class = lr$npc, T_tasks = 20, true_K = lr$tK,
                          separation = lr$kap, heterogeneity = lr$sig,
                          seed = lr$seed, dgp = dgp)
    mm <- klue_mmnl(data$database, n_draws = n_draws, dgp = dgp)
    if (is.finite(lr$mnl_bic) && mm$converged) {
      bics <- c(MNL = unname(lr$mnl_bic), LCMNL = unname(lr$lc_bic),
                MMNL = unname(mm$BIC))
      winner <- names(which.min(bics))
      rows[[i]] <- data.frame(true_K = lr$tK, kappa = lr$kap, sigma = lr$sig,
                              mnl_BIC = lr$mnl_bic, lcmnl_BIC = lr$lc_bic,
                              lcmnl_C = lr$lc_C, lcmnl_method = lr$lc_method,
                              mmnl_BIC = mm$BIC, bic_prefers = winner,
                              stringsAsFactors = FALSE)
      if (verbose) {
        lc_str <- if (!is.finite(lr$lc_bic)) "LCMNL(-)"
                  else if (lr$lc_C == 1L) "MNL(C=1)"
                  else sprintf("LCMNL(C=%d)", lr$lc_C)
        cat(sprintf("%s -> %s [%s]\n", lc_str, winner,
                    ifelse(is.na(lr$lc_method), "-", lr$lc_method)))
      }
    } else if (verbose) cat("FAILED\n")
  }
  df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  if (verbose && !is.null(df) && nrow(df) > 0) {
    n_tot <- nrow(df)
    counts <- table(factor(df$bic_prefers, levels = c("MNL", "LCMNL", "MMNL")))
    cat(sprintf("  Overall: MNL %d/%d (%.0f%%)  LCMNL %d/%d (%.0f%%)  MMNL %d/%d (%.0f%%)\n",
                counts["MNL"], n_tot, 100 * counts["MNL"] / n_tot,
                counts["LCMNL"], n_tot, 100 * counts["LCMNL"] / n_tot,
                counts["MMNL"], n_tot, 100 * counts["MMNL"] / n_tot))
    for (lbl in c("K=1 (continuous DGP):", "K>1 (discrete DGP):  ")) {
      sub <- if (startsWith(lbl, "K=1")) df[df$true_K == 1, ] else df[df$true_K > 1, ]
      if (nrow(sub) > 0) {
        cat(sprintf("  %s MNL %d  LCMNL %d  MMNL %d  (of %d)\n", lbl,
                    sum(sub$bic_prefers == "MNL"), sum(sub$bic_prefers == "LCMNL"),
                    sum(sub$bic_prefers == "MMNL"), nrow(sub)))
      }
    }
  }
  df
}

# ---- H1 ablations: convergence / initialisation / estimator -----------------

#' Monte Carlo: clustering start vs random restarts (convergence)
#'
#' Tests whether the clustering-based start reaches the same log-likelihood
#' optimum as a budget of diffuse random restarts, and how many random starts
#' are needed to match it.
#'
#' @param n_random number of diffuse random restarts per condition.
#' @param n_cond number of conditions to run from the H1 grid.
#' @param verbose if \code{TRUE}, print a summary.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a shared blocked D-efficient design;
#'   otherwise use the per-respondent random-attribute path.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param seg_scale optional segment-membership scaling passed to the simulator.
#' @param T_random number of tasks per respondent in the random-attribute path.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per converged condition: true class count
#'   K, whether the clustering start is best, the fraction of random starts at
#'   the global and at local optima, the first random start that matches the
#'   clustering log-likelihood, and the clustering-to-random time ratio.
#' @export
klue_study_convergence <- function(n_random = 50, n_cond = 40, verbose = TRUE,
                                   dgp = DGP_DEFAULT,
                                   blocked = FALSE, n_cards = 48L, n_blocks = 4L,
                                   seg_scale = NULL, T_random = 20L,
                                   n_cores = .klue_cores()) {
  if (verbose) cat(sprintf("\n==== SUPP 1: CONVERGENCE%s ====\n",
                           if (blocked) " (blocked design)" else ""))
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  conds <- .h1_conditions(n_cond)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = T_random, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp, design = design, seg_scale = seg_scale)

    t1 <- system.time(cres <- klue_lcmnl(data$database, tK, dgp = dgp))[3]
    rLLs <- rep(-Inf, n_random)
    t2 <- system.time({
      for (r in 1:n_random) {
        st <- .random_start(seed, r, tK, dgp)
        rr <- estimate_lcmnl(data$database, tK, start_betas = st$betas, dgp = dgp)
        if (rr$converged) rLLs[r] <- rr$LL
      }
    })[3]

    if (!cres$converged) return(list(valid = FALSE))
    clust_LL <- cres$LL
    match_at <- NA_integer_
    for (r in 1:n_random) {
      if (rLLs[r] >= clust_LL - 0.1) { match_at <- r; break }
    }
    list(valid = TRUE, K = tK, kappa = kap,
         cluster_best = is.na(match_at) || abs(clust_LL - max(rLLs)) < 0.1,
         pct_global = sum(abs(rLLs - clust_LL) < 0.1) / n_random,
         pct_local = sum(rLLs > -Inf & rLLs < clust_LL - 0.1) / n_random,
         match_at = match_at,
         time_ratio = t2 / max(t1, 0.01))
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  vres <- res_list[.col(res_list, "valid")]
  df <- data.frame(K = .col(vres, "K"), cluster_best = .col(vres, "cluster_best"),
                   pct_global = .col(vres, "pct_global"),
                   pct_local = .col(vres, "pct_local"),
                   match_at = .col(vres, "match_at"),
                   time_ratio = .col(vres, "time_ratio"))
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

#' Monte Carlo: initialisation feature ablation
#'
#' Three-arm ablation comparing clustering starts built from RP contrasts, from
#' one-hot indicators, and uninformed random starts, recording how often each
#' arm reaches the global log-likelihood (best across all arms).
#'
#' @param n_random number of uninformed random starts per condition.
#' @param n_cond number of conditions to run from the H1 grid.
#' @param attr_corr optional attribute correlation passed to the simulator;
#'   \code{NULL} for uncorrelated attributes.
#' @param verbose if \code{TRUE}, print a summary.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a shared blocked D-efficient design;
#'   otherwise use the per-respondent random-attribute path.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param seg_scale optional segment-membership scaling passed to the simulator.
#' @param T_random number of tasks per respondent in the random-attribute path.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per valid condition: true class count K,
#'   kappa, sigma, whether the RP-contrast and one-hot starts reach the global
#'   optimum, the per-start fraction of random starts at the global optimum, the
#'   RP, one-hot, and best log-likelihoods, and the RP-minus-one-hot gap.
#' @export
klue_study_initialisation <- function(n_random = 50, n_cond = 40,
                                      attr_corr = NULL,
                                      verbose = TRUE, dgp = DGP_DEFAULT,
                                      blocked = FALSE, n_cards = 48L, n_blocks = 4L,
                                      seg_scale = NULL, T_random = 20L,
                                      n_cores = .klue_cores()) {
  # Three-arm feature ablation: RP contrasts vs one-hot indicators vs
  # uninformed random starts. "Global" LL = best across all arms (tol 0.1).
  if (verbose) cat(sprintf("\n==== INITIALISATION ABLATION (RP / one-hot / random)%s%s ====\n",
                           if (is.null(attr_corr)) "" else sprintf(" (attr_corr=%.2f)", attr_corr),
                           if (blocked) " (blocked design)" else ""))
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  conds <- .h1_conditions(n_cond)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = T_random, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp, attr_corr = attr_corr,
                          design = design, seg_scale = seg_scale)

    arm_LL <- function(res) if (!is.null(res) && isTRUE(res$converged)) res$LL else -Inf
    rp_LL <- arm_LL(tryCatch(klue_lcmnl(data$database, tK, dgp = dgp),
                             error = function(e) NULL))
    oh_LL <- arm_LL(tryCatch(klue_lcmnl(data$database, tK, dgp = dgp,
                                        feature_type = "onehot"),
                             error = function(e) NULL))
    rLLs <- rep(-Inf, n_random)
    for (r in 1:n_random) {
      st <- .random_start(seed, r, tK, dgp)
      rr <- tryCatch(estimate_lcmnl(data$database, tK, start_betas = st$betas,
                                    dgp = dgp), error = function(e) NULL)
      if (!is.null(rr) && isTRUE(rr$converged)) rLLs[r] <- rr$LL
    }

    best_LL <- max(c(rp_LL, oh_LL, rLLs, -Inf))
    tol <- 0.1
    list(valid = (rp_LL > -Inf) || (oh_LL > -Inf) || any(rLLs > -Inf),
         K = tK, kappa = kap, sigma = sig,
         rp_at_global = (rp_LL >= best_LL - tol),
         oh_at_global = (oh_LL >= best_LL - tol),
         random_pct_global = mean(rLLs >= best_LL - tol),
         rp_LL = rp_LL, oh_LL = oh_LL, best_LL = best_LL,
         gap_rp_oh = rp_LL - oh_LL)
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  vres <- res_list[sapply(res_list, function(x) isTRUE(x$valid))]
  df <- data.frame(
    K                 = .col(vres, "K"),
    kappa             = .col(vres, "kappa"),
    sigma             = .col(vres, "sigma"),
    rp_at_global      = .col(vres, "rp_at_global"),
    oh_at_global      = .col(vres, "oh_at_global"),
    random_pct_global = .col(vres, "random_pct_global"),
    rp_LL             = .col(vres, "rp_LL"),
    oh_LL             = .col(vres, "oh_LL"),
    best_LL           = .col(vres, "best_LL"),
    gap_rp_oh         = .col(vres, "gap_rp_oh")
  )
  if (verbose && nrow(df) > 0) {
    cat(sprintf("  RP contrasts reach global:    %.1f%% (%d/%d conditions)\n",
                100 * mean(df$rp_at_global), sum(df$rp_at_global), nrow(df)))
    cat(sprintf("  One-hot indicators at global: %.1f%% (%d/%d)\n",
                100 * mean(df$oh_at_global), sum(df$oh_at_global), nrow(df)))
    cat(sprintf("  Random starts at global:      %.1f%% per start\n",
                100 * mean(df$random_pct_global)))
    by_K <- aggregate(cbind(rp_at_global, oh_at_global, random_pct_global)
                      ~ K, data = df, FUN = mean)
    cat("\n  By true K:\n")
    print(by_K, row.names = FALSE)
  }
  df
}

#' Monte Carlo: estimator x start-strategy benchmark
#'
#' Crosses two estimators (maximum likelihood and EM) with four start
#' strategies (clustering, MNL perturbation, random partition, diffuse random),
#' recording per-start at-global rate, best-of-budget success, and the number of
#' starts needed to reach the global optimum (the best log-likelihood across all
#' cells of a condition).
#'
#' @param n_random number of diffuse random starts per condition.
#' @param n_perturb number of perturbation and random-partition starts per
#'   condition.
#' @param n_cond number of conditions to run from the H1 grid.
#' @param jitter_sd standard deviation of the MNL perturbation jitter.
#' @param attr_corr optional attribute correlation passed to the simulator;
#'   \code{NULL} for uncorrelated attributes.
#' @param verbose if \code{TRUE}, print a summary.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a shared blocked D-efficient design;
#'   otherwise use the per-respondent random-attribute path.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param seg_scale optional segment-membership scaling passed to the simulator.
#' @param T_random number of tasks per respondent in the random-attribute path.
#' @param cond_idx optional vector of condition indices for chunked or
#'   checkpointed runs; \code{NULL} uses the first \code{n_cond} conditions.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per estimator x strategy x condition cell:
#'   true class count K, kappa, sigma, rep, estimator, strategy, number of
#'   starts, per-start at-global rate, whether the global optimum was reached,
#'   starts-to-global, and the best and global log-likelihoods.
#' @export
klue_study_estimator <- function(n_random = 50L, n_perturb = 50L, n_cond = 36L,
                                 jitter_sd = 0.5, attr_corr = NULL,
                                 verbose = TRUE, dgp = DGP_DEFAULT,
                                 blocked = FALSE, n_cards = 48L, n_blocks = 4L,
                                 seg_scale = NULL, T_random = 20L,
                                 cond_idx = NULL, n_cores = .klue_cores()) {
  # {ML, EM} x {clustering, perturbation, random partition, diffuse random}.
  # Per condition: per-start at-global rate, best-of-budget success, and
  # starts-to-global; the global optimum is the best LL across all cells.
  # cond_idx selects specific conditions for chunked/checkpointed runs.
  if (verbose) cat(sprintf("\n==== ESTIMATOR x START-STRATEGY BENCHMARK%s%s ====\n",
                           if (is.null(attr_corr)) "" else sprintf(" (attr_corr=%.2f)", attr_corr),
                           if (blocked) " (blocked design)" else ""))
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  conds <- .h1_conditions(n_cond, cond_idx)
  tol <- 0.1
  est_fns <- list(ml = estimate_lcmnl, em = estimate_lcmnl_em)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = T_random, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp, attr_corr = attr_corr,
                          design = design, seg_scale = seg_scale)
    db <- data$database

    # Pre-draw start specs ONCE so both estimators see identical starts.
    perturb_starts <- get_mnl_perturbation_starts(db, tK, n_starts = n_perturb,
                                                  jitter_sd = jitter_sd,
                                                  seed = seed, dgp = dgp)
    partition_starts <- get_random_partition_starts(db, tK, n_starts = n_perturb,
                                                    seed = seed, dgp = dgp)
    random_starts <- lapply(seq_len(n_random), function(r)
      .random_start(seed, r, tK, dgp))

    arm_LLs <- list()
    for (est in names(est_fns)) {
      fit_one <- est_fns[[est]]
      cl <- tryCatch(klue_lcmnl(db, tK, dgp = dgp, estimator = est),
                     error = function(e) NULL)
      cl_LLs <- if (!is.null(cl) && length(cl$method_results)) {
        vapply(cl$method_results,
               function(r) if (isTRUE(r$converged)) r$LL else -Inf, numeric(1))
      } else -Inf
      fit_starts <- function(starts) vapply(starts, function(s) {
        r <- tryCatch(fit_one(db, tK, start_betas = s$betas,
                              start_shares = s$shares, dgp = dgp),
                      error = function(e) NULL)
        if (!is.null(r) && isTRUE(r$converged)) r$LL else -Inf
      }, numeric(1))
      arm_LLs[[est]] <- list(clustering       = unname(cl_LLs),
                             perturbation     = fit_starts(perturb_starts),
                             random_partition = fit_starts(partition_starts),
                             random           = fit_starts(random_starts))
    }

    global_LL <- max(unlist(arm_LLs), -Inf)
    if (!is.finite(global_LL)) return(NULL)
    rows <- list()
    for (est in names(arm_LLs)) {
      for (strat in names(arm_LLs[[est]])) {
        lls <- arm_LLs[[est]][[strat]]
        at_global <- lls >= global_LL - tol
        first_hit <- which(at_global)
        rows[[length(rows) + 1L]] <- data.frame(
          K = tK, kappa = kap, sigma = sig, rep = rp,
          estimator = est, strategy = strat,
          n_starts = length(lls),
          at_global_rate = mean(at_global),
          reached_global = any(at_global),
          starts_to_global = if (length(first_hit)) first_hit[1] else NA_integer_,
          best_LL = max(lls, -Inf), global_LL = global_LL,
          stringsAsFactors = FALSE)
      }
    }
    do.call(rbind, rows)
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  df <- do.call(rbind, res_list[!vapply(res_list, is.null, logical(1))])

  if (verbose && !is.null(df) && nrow(df) > 0) {
    cat(sprintf("  %d conditions x 2 estimators x 3 strategies\n",
                length(unique(paste(df$K, df$kappa, df$sigma, df$rep)))))
    agg <- aggregate(cbind(at_global_rate, reached_global) ~ estimator + strategy,
                     data = df, FUN = mean)
    agg <- agg[order(agg$estimator, -agg$at_global_rate), ]
    cat("\n  Per-start at-global rate and best-of-budget success:\n")
    for (k in seq_len(nrow(agg))) {
      cat(sprintf("    %-3s x %-12s  per-start %5.1f%%   reached-global %5.1f%%\n",
                  agg$estimator[k], agg$strategy[k],
                  100 * agg$at_global_rate[k], 100 * agg$reached_global[k]))
    }
  }
  df
}

# ---- robustness arms ----------------------------------------------------------

#' Monte Carlo: robustness to unbalanced class proportions
#'
#' Tests BIC class-count recovery for a three-class DGP as class proportions
#' move from balanced to severely unbalanced.
#'
#' @param verbose if \code{TRUE}, print per-configuration accuracy.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a shared blocked D-efficient design;
#'   otherwise use the per-respondent random-attribute path.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per imbalance configuration: the
#'   configuration name and the BIC class-count selection accuracy (percent).
#' @export
klue_study_unbalanced <- function(verbose = TRUE, dgp = DGP_DEFAULT,
                                  blocked = FALSE, n_cards = 48L, n_blocks = 4L,
                                  n_cores = .klue_cores()) {
  if (verbose) cat(sprintf("\n==== SUPP 3: UNBALANCED%s ====\n",
                           if (blocked) " (blocked design)" else ""))
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  configs <- list(
    list(name = "balanced", props = c(1/3, 1/3, 1/3)),
    list(name = "mild",     props = c(0.5, 0.3, 0.2)),
    list(name = "moderate", props = c(0.6, 0.25, 0.15)),
    list(name = "severe",   props = c(0.7, 0.2, 0.1))
  )
  conds <- expand.grid(kappa = c(0.75, 1.0, 1.25), rep = 1:5)
  all_jobs <- expand.grid(ci = seq_along(configs), cond_i = 1:nrow(conds))

  .run_one <- function(row) {
    ci <- all_jobs$ci[row]; cfg <- configs[[ci]]
    kap <- conds$kappa[all_jobs$cond_i[row]]
    rp  <- conds$rep[all_jobs$cond_i[row]]
    seed <- as.integer(3000 + 1000 * ci + 100 * (kap * 100) + rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 3,
                          separation = kap, heterogeneity = 0.2, seed = seed,
                          class_proportions = cfg$props, dgp = dgp, design = design)
    bics <- .bic_scan(data$database, 1:5, dgp)
    list(ci = ci, correct = as.integer(which.min(bics) == 3))
  }

  res_list <- .study_mclapply(nrow(all_jobs), .run_one, n_cores)
  accs <- vapply(seq_along(configs), function(ci) {
    100 * mean(.col(res_list[.col(res_list, "ci") == ci], "correct"))
  }, numeric(1))
  if (verbose) {
    for (ci in seq_along(configs))
      cat(sprintf("  %s: %.1f%%\n", configs[[ci]]$name, accs[ci]))
  }
  data.frame(config = vapply(configs, `[[`, character(1), "name"),
             accuracy = accs)
}

#' Monte Carlo: random vs D-efficient experimental design
#'
#' Compares BIC class-count recovery on data generated with a random-attribute
#' design against data generated with a D-efficient design, holding the DGP
#' otherwise fixed.
#'
#' @param verbose if \code{TRUE}, print recovery rates for each design.
#' @param dgp data-generating-process specification.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row: the BIC class-count recovery accuracy
#'   (percent) for the random design and for the D-efficient design.
#' @export
klue_study_design <- function(verbose = TRUE, dgp = DGP_DEFAULT,
                              n_cores = .klue_cores()) {
  if (verbose) cat("\n==== SUPP 4: D-EFFICIENT ====\n")
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0), rep = 1:5)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]; rp <- conds$rep[i]
    seed <- as.integer(4000 + 100 * tK + 10 * (kap * 100) + rp)
    dr <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                        separation = kap, heterogeneity = 0.2, seed = seed,
                        dgp = dgp)
    dd <- klue_simulate_deff(N_per_class = 150, T_tasks = 20, true_K = tK,
                             separation = kap, heterogeneity = 0.2, seed = seed,
                             dgp = dgp)
    list(r_ok = as.integer(which.min(.bic_scan(dr$database, 1:5, dgp)) == tK),
         d_ok = as.integer(which.min(.bic_scan(dd$database, 1:5, dgp)) == tK))
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  rok <- sum(.col(res_list, "r_ok")); dok <- sum(.col(res_list, "d_ok"))
  if (verbose) cat(sprintf("  Random: %.1f%%  D-eff: %.1f%%\n",
                           100 * rok / nrow(conds), 100 * dok / nrow(conds)))
  data.frame(random = 100 * rok / nrow(conds), deff = 100 * dok / nrow(conds))
}

#' Monte Carlo: concomitant covariates on class membership
#'
#' Tests BIC class-count recovery and clustering recovery (ARI) when class
#' membership depends on a concomitant covariate of varying strength.
#'
#' @param verbose if \code{TRUE}, print accuracy and ARI.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a shared blocked D-efficient design;
#'   otherwise use the per-respondent random-attribute path.
#' @param n_cards number of design cards for the blocked design.
#' @param n_blocks number of blocks for the blocked design.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row: the BIC class-count selection accuracy
#'   (percent) and the mean ARI across conditions.
#' @export
klue_study_concomitant <- function(verbose = TRUE, dgp = DGP_DEFAULT,
                                   blocked = FALSE, n_cards = 48L, n_blocks = 4L,
                                   n_cores = .klue_cores()) {
  if (verbose) cat(sprintf("\n==== SUPP 5: CONCOMITANT%s ====\n",
                           if (blocked) " (blocked design)" else ""))
  design <- if (blocked)
    klue_design(n_cards = n_cards, n_blocks = n_blocks, dgp = dgp) else NULL
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0, 1.25),
                       cs = c(0.5, 1.0, 1.5), rep = 1:3)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    cs <- conds$cs[i]; rp <- conds$rep[i]
    seed <- as.integer(5000 + 100 * tK + 10 * (kap * 100) + rp + cs * 1000)
    data <- klue_simulate_cov(
      N_per_class = 150, T_tasks = 20, true_K = tK, separation = kap,
      heterogeneity = 0.2, seed = seed, covariate_strength = cs,
      dgp = dgp, design = design
    )
    bics <- .bic_scan(data$database, 1:5, dgp)
    mt <- klue_lcmnl(data$database, tK, dgp = dgp)
    ari_val <- NA
    if (mt$converged && tK > 1) {
      pred <- apply(mt$posteriors, 1, which.max)
      ari_val <- compute_ari(data$true_class, pred)
    }
    list(correct = as.integer(which.min(bics) == tK), ari = ari_val)
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  ok   <- sum(.col(res_list, "correct"))
  aris <- na.omit(.col(res_list, "ari"))
  if (verbose) cat(sprintf("  Accuracy: %.1f%%  ARI: %.3f\n",
                           100 * ok / nrow(conds), mean(aris)))
  data.frame(accuracy = 100 * ok / nrow(conds), mean_ari = mean(aris))
}

#' Monte Carlo: parameter and clustering recovery at the true class count
#'
#' Fits LCMNL at the true class count and reports how well the class-specific
#' coefficients and the class assignments are recovered.
#'
#' @param n_cond number of conditions to run from the condition grid.
#' @param verbose if \code{TRUE}, print mean recovery metrics.
#' @param dgp data-generating-process specification.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row: mean coefficient RMSE, mean bias, and mean
#'   ARI across valid conditions.
#' @export
klue_study_recovery <- function(n_cond = 80, verbose = TRUE, dgp = DGP_DEFAULT,
                                n_cores = .klue_cores()) {
  if (verbose) cat("\n==== SUPP 6: RECOVERY ====\n")
  conds <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                       sigma = c(0.15, 0.25), rep = 1:3)
  conds <- conds[1:min(n_cond, nrow(conds)), ]

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    m <- klue_lcmnl(data$database, tK, dgp = dgp)
    if (!m$converged || tK <= 1) return(list(valid = FALSE))
    r <- compute_recovery(data$true_betas, m$betas)
    pred <- apply(m$posteriors, 1, which.max)
    list(valid = TRUE, rmse = r$rmse, bias = r$bias,
         ari = compute_ari(data$true_class, pred))
  }

  res_list <- .study_mclapply(nrow(conds), .run_one, n_cores)
  vres <- res_list[.col(res_list, "valid")]
  out <- data.frame(rmse = mean(.col(vres, "rmse")),
                    bias = mean(.col(vres, "bias")),
                    ari = mean(.col(vres, "ari")))
  if (verbose) cat(sprintf("  RMSE: %.4f  Bias: %.4f  ARI: %.3f\n",
                           out$rmse, out$bias, out$ari))
  out
}

#' Monte Carlo: comparison of clustering start methods
#'
#' Compares the clustering methods in \code{KLUE_CLUSTER_METHODS} as sources of
#' starting values, recording per-method convergence, BIC class-count accuracy,
#' and which method attains the highest log-likelihood per condition.
#'
#' @param verbose if \code{TRUE}, print per-condition and per-method summaries.
#' @param dgp data-generating-process specification.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per condition: true_K, kappa, sigma, the
#'   best method by log-likelihood, and per-method BIC-correct and converged
#'   indicator columns.
#' @export
klue_study_clustering <- function(verbose = TRUE, dgp = DGP_DEFAULT,
                                  n_cores = .klue_cores()) {
  if (verbose) cat("\n==== SUPP 7: CLUSTERING METHODS ====\n")
  method_names <- KLUE_CLUSTER_METHODS
  conds <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                       sigma = c(0.15, 0.25), rep = 1:3)
  nc <- nrow(conds)

  .run_one <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp, base = 6000L)
    data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = tK,
                          separation = kap, heterogeneity = sig, seed = seed,
                          dgp = dgp)
    all_starts <- get_all_starts(data$database, tK, dgp = dgp)

    m_correct <- setNames(integer(length(method_names)), method_names)
    m_conv    <- setNames(integer(length(method_names)), method_names)
    m_ll      <- setNames(rep(NA_real_, length(method_names)), method_names)

    for (nm in method_names) {
      if (is.null(all_starts[[nm]])) next
      res <- tryCatch(
        estimate_lcmnl(data$database, tK,
                       start_betas = all_starts[[nm]]$betas,
                       start_shares = all_starts[[nm]]$shares, dgp = dgp),
        error = function(e) NULL
      )
      if (is.null(res) || !res$converged) next
      m_conv[nm] <- 1L
      m_ll[nm] <- res$LL
      # BIC scan across C with THIS method's starts only.
      bics <- rep(Inf, 5)
      for (Cc in 1:5) {
        if (Cc == tK) { bics[Cc] <- res$BIC; next }
        s_cc <- tryCatch(klue_starts(data$database, Cc, method = nm, dgp = dgp),
                         error = function(e) NULL)
        if (is.null(s_cc)) next
        m_cc <- tryCatch(
          estimate_lcmnl(data$database, Cc, start_betas = s_cc$betas,
                         start_shares = s_cc$shares, dgp = dgp),
          error = function(e) NULL
        )
        if (!is.null(m_cc) && m_cc$converged) bics[Cc] <- m_cc$BIC
      }
      if (which.min(bics) == tK) m_correct[nm] <- 1L
    }

    best_idx <- which.max(m_ll)
    list(true_K = tK, kappa = kap, sigma = sig,
         m_correct = m_correct, m_conv = m_conv, m_ll = m_ll,
         best_method = if (length(best_idx) > 0) method_names[best_idx] else NA_character_)
  }

  res_list <- .study_mclapply(nc, .run_one, n_cores)
  method_correct   <- t(sapply(res_list, `[[`, "m_correct"))
  method_converged <- t(sapply(res_list, `[[`, "m_conv"))
  best_methods     <- .col(res_list, "best_method")
  cond_K <- .col(res_list, "true_K"); cond_kap <- .col(res_list, "kappa")
  cond_sig <- .col(res_list, "sigma")

  if (verbose) {
    for (i in 1:nc)
      cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... best=%s\n",
                  i, nc, cond_K[i], cond_kap[i], cond_sig[i], best_methods[i]))
    cat("\n  Per-method BIC accuracy (%):\n")
    for (nm in method_names) {
      valid <- method_converged[, nm] == 1
      if (any(valid)) {
        cat(sprintf("    %-15s: %.1f%% (converged %d/%d)\n", nm,
                    100 * mean(method_correct[valid, nm]), sum(valid), nc))
      }
    }
    cat(sprintf("\n  Best method distribution: %s\n",
                paste(names(table(best_methods)), table(best_methods),
                      sep = "=", collapse = ", ")))
  }

  data.frame(true_K = cond_K, kappa = cond_kap, sigma = cond_sig,
             best_method = best_methods, method_correct, method_converged)
}

#' Monte Carlo: sensitivity to sample size and panel length
#'
#' Tests BIC class-count recovery across a grid of respondents per class and
#' tasks per respondent. Under blocking, the number of blocks is held fixed and
#' cards per block vary with panel length so panel length is not confounded with
#' design diversity.
#'
#' @param verbose if \code{TRUE}, print accuracy for each sample-size cell.
#' @param dgp data-generating-process specification.
#' @param blocked if \code{TRUE}, use a per-panel-length blocked D-efficient
#'   design; otherwise use the per-respondent random-attribute path.
#' @param n_cores cores for condition-level parallelism.
#' @return A data.frame with one row per cell: tasks per respondent (T_tasks),
#'   respondents per class (N_per_class), and BIC class-count accuracy (percent).
#' @export
klue_study_sample <- function(verbose = TRUE, dgp = DGP_DEFAULT, blocked = FALSE,
                              n_cores = .klue_cores()) {
  if (verbose) cat(sprintf("\n==== SUPP 8: SAMPLE SIZE / PANEL LENGTH%s ====\n",
                           if (blocked) " (blocked design)" else ""))
  T_values   <- c(8L, 12L, 20L)
  Npc_values <- c(50L, 100L, 150L)
  # Under blocking, hold the NUMBER OF BLOCKS fixed at 4 and vary cards-per-
  # block (= T), so panel length is not confounded with design diversity.
  blk_spec <- list(`8` = c(32L, 4L), `12` = c(48L, 4L), `20` = c(80L, 4L))
  designs <- if (blocked) setNames(lapply(T_values, function(Tt) {
               p <- blk_spec[[as.character(Tt)]]
               klue_design(n_cards = p[1], n_blocks = p[2], dgp = dgp)
             }), as.character(T_values)) else NULL
  conds <- expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0),
                       sigma = 0.20, rep = 1:3)
  all_jobs <- expand.grid(Tv_idx = seq_along(T_values),
                          Npc_idx = seq_along(Npc_values),
                          cond_i = 1:nrow(conds))

  .run_one <- function(row) {
    Tv  <- T_values[all_jobs$Tv_idx[row]]
    Npc <- Npc_values[all_jobs$Npc_idx[row]]
    i   <- all_jobs$cond_i[row]
    tK <- conds$true_K[i]; kap <- conds$kappa[i]; rp <- conds$rep[i]
    seed <- as.integer(7000 + 100 * tK + 10 * (kap * 100) + rp + Tv * 100 + Npc)
    dsg  <- if (blocked) designs[[as.character(Tv)]] else NULL
    data <- klue_simulate(N_per_class = Npc, T_tasks = Tv, true_K = tK,
                          separation = kap, heterogeneity = 0.2, seed = seed,
                          dgp = dgp, design = dsg)
    bics <- .bic_scan(data$database, 1:5, dgp)
    list(Tv = Tv, Npc = Npc, correct = as.integer(which.min(bics) == tK))
  }

  res_list <- .study_mclapply(nrow(all_jobs), .run_one, n_cores)
  grid <- expand.grid(Npc = Npc_values, Tv = T_values)[, 2:1]
  grid <- data.frame(Tv = rep(T_values, each = length(Npc_values)),
                     Npc = rep(Npc_values, times = length(T_values)))
  accs <- vapply(seq_len(nrow(grid)), function(g) {
    sel <- sapply(res_list, function(x) x$Tv == grid$Tv[g] & x$Npc == grid$Npc[g])
    100 * mean(.col(res_list[sel], "correct"))
  }, numeric(1))
  if (verbose) {
    for (g in seq_len(nrow(grid)))
      cat(sprintf("  T=%2d  N/class=%3d  (N_total=%d-%d): %.1f%%\n",
                  grid$Tv[g], grid$Npc[g], grid$Npc[g] * 2, grid$Npc[g] * 3, accs[g]))
  }
  data.frame(T_tasks = grid$Tv, N_per_class = grid$Npc, accuracy = accs)
}

#' Monte Carlo: correlated MMNL robustness
#'
#' Compares LCMNL against MMNL with independent random coefficients and MMNL
#' with correlated random coefficients, recording which model BIC prefers.
#'
#' @param n_draws number of draws for MMNL estimation.
#' @param verbose if \code{TRUE}, print per-condition progress and a summary.
#' @param dgp data-generating-process specification.
#' @param n_cores cores for the parallel LCMNL phase.
#' @return A data.frame with one row per converged condition: true_K, kappa,
#'   sigma, the LCMNL BIC and selected class count, the independent and
#'   correlated MMNL BIC values, and which model BIC prefers.
#' @export
klue_study_mmnl_corr <- function(n_draws = N_DRAWS_MMNL, verbose = TRUE,
                                 dgp = DGP_DEFAULT, n_cores = .klue_cores()) {
  if (verbose) cat("\n==== SUPP 9: CORRELATED MMNL ROBUSTNESS ====\n")
  conds <- rbind(
    expand.grid(true_K = 1, kappa = 0, sigma = c(0.15, 0.25, 0.35), rep = 1:2),
    expand.grid(true_K = c(2, 3), kappa = c(0.75, 1.0),
                sigma = c(0.15, 0.25), rep = 1:2))
  nc <- nrow(conds)
  if (verbose) cat(sprintf("  %d conditions (K=1: %d, K>=2: %d)\n",
                           nc, sum(conds$true_K == 1), sum(conds$true_K > 1)))

  if (verbose) cat("  Phase 1: LCMNL estimation (parallel)...\n")
  .run_lcmnl_part <- function(i) {
    tK <- conds$true_K[i]; kap <- conds$kappa[i]
    sig <- conds$sigma[i]; rp <- conds$rep[i]
    seed <- .cond_seed(tK, kap, sig, rp)
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
  lc_results <- .study_mclapply(nc, .run_lcmnl_part, n_cores)

  if (verbose) cat("  Phase 2: MMNL estimation (sequential, Apollo)...\n")
  rows <- vector("list", nc)
  for (i in 1:nc) {
    lr <- lc_results[[i]]
    if (verbose) cat(sprintf("  [%2d/%d] K=%d kappa=%.2f sigma=%.2f ... ",
                             i, nc, lr$tK, lr$kap, lr$sig))
    data <- klue_simulate(N_per_class = lr$npc, T_tasks = 20, true_K = lr$tK,
                          separation = lr$kap, heterogeneity = lr$sig,
                          seed = lr$seed, dgp = dgp)
    mm_indep <- klue_mmnl(data$database, n_draws = n_draws, dgp = dgp)
    mm_corr  <- klue_mmnl_corr(data$database, n_draws = n_draws, dgp = dgp)
    if (is.finite(lr$lc_bic) && mm_indep$converged && mm_corr$converged) {
      bics <- c(LCMNL = lr$lc_bic, MMNL_indep = mm_indep$BIC,
                MMNL_corr = mm_corr$BIC)
      winner <- names(which.min(bics))
      rows[[i]] <- data.frame(true_K = lr$tK, kappa = lr$kap, sigma = lr$sig,
                              lcmnl_BIC = lr$lc_bic, lcmnl_C = lr$lc_C,
                              mmnl_indep_BIC = mm_indep$BIC,
                              mmnl_corr_BIC = mm_corr$BIC,
                              bic_prefers = winner, stringsAsFactors = FALSE)
      if (verbose) {
        lc_lab <- if (lr$lc_C == 1L) "MNL(C=1)" else sprintf("LCMNL(C=%d)", lr$lc_C)
        cat(sprintf("%s vs indep vs corr -> %s\n", lc_lab, winner))
      }
    } else if (verbose) cat("FAILED\n")
  }
  df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])

  if (verbose && !is.null(df) && nrow(df) > 0) {
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

# ---- master driver ------------------------------------------------------------

#' Run the full Frings (2026) Monte Carlo study
#'
#' Master driver that runs the main arm, the MMNL comparison, and the
#' supplementary robustness arms, writing each result to a CSV in
#' \code{OUTPUT_DIR} and returning all results in one list.
#'
#' @param run_main if \code{TRUE}, run the main class-count recovery arm.
#' @param run_mmnl if \code{TRUE}, run the MNL vs LCMNL vs MMNL comparison.
#' @param run_supp if \code{TRUE}, run the supplementary robustness arms.
#' @param verbose if \code{TRUE}, print progress and per-arm summaries.
#' @param dgp data-generating-process specification forwarded to every driver.
#' @param n_cores cores for condition-level parallelism, forwarded to every
#'   driver. \code{n_cores = 1} runs fully sequentially (lowest memory,
#'   reproducible timing, friendly to a shared machine); the default
#'   \code{.klue_cores()} parallelises across conditions. Output is identical
#'   either way (per-condition seeds). Set it to match your time budget.
#' @return A named list of the result data.frames for each arm that was run
#'   (\code{main}, \code{mmnl}, and the supplementary arms).
#' @export
klue_study <- function(run_main = TRUE, run_mmnl = TRUE, run_supp = TRUE,
                       verbose = TRUE, dgp = DGP_DEFAULT,
                       n_cores = .klue_cores()) {
  results <- list()
  t_start <- Sys.time()
  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE,
                                          showWarnings = FALSE)
  save_csv <- function(name, df) {
    write.csv(df, file.path(OUTPUT_DIR, name), row.names = FALSE)
  }

  if (run_main) {
    results[["main"]] <- klue_study_main(verbose = verbose, dgp = dgp,
                                         n_cores = n_cores)
    summarise_main_results(results[["main"]])
    save_csv("main_results.csv", results[["main"]])
  }
  if (run_mmnl) {
    results[["mmnl"]] <- klue_study_mmnl(verbose = verbose, dgp = dgp,
                                         n_cores = n_cores)
    save_csv("mmnl_results.csv", results[["mmnl"]])
  }
  if (run_supp) {
    supp <- list(
      convergence = list(fn = klue_study_convergence, csv = "convergence_results.csv"),
      unbalanced  = list(fn = klue_study_unbalanced,  csv = "unbalanced_results.csv"),
      design      = list(fn = klue_study_design,      csv = "design_results.csv"),
      concomitant = list(fn = klue_study_concomitant, csv = "concomitant_results.csv"),
      recovery    = list(fn = klue_study_recovery,    csv = "recovery_results.csv"),
      clustering  = list(fn = klue_study_clustering,  csv = "clustering_comparison_results.csv"),
      sample_sens = list(fn = klue_study_sample,      csv = "sample_sensitivity_results.csv"),
      mmnl_corr   = list(fn = klue_study_mmnl_corr,   csv = "mmnl_correlated_results.csv")
    )
    for (nm in names(supp)) {
      results[[nm]] <- supp[[nm]]$fn(verbose = verbose, dgp = dgp,
                                     n_cores = n_cores)
      save_csv(supp[[nm]]$csv, results[[nm]])
    }
  }

  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  if (verbose) cat(sprintf("\n==== DONE: %.1f min ====\n", elapsed))
  results
}
