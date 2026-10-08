# =============================================================================
# STATISTICAL TESTS AND CONFIDENCE INTERVALS FOR MANUSCRIPT
# =============================================================================
#
# Reads output/{main,mmnl,mmnl_correlated,recovery,convergence,unbalanced,
# design,concomitant,sample_sensitivity}_results.csv and computes:
#   1. Wilson confidence intervals for all proportions
#   2. McNemar test for BIC vs ICL (paired comparison)
#   3. Logistic regression of BIC success on design factors
#   4. Exact tests for MMNL comparison
#   5. Bootstrap CIs for means (RMSE, bias, ARI)
#   6. Wilcoxon test for convergence comparison
#
# Run from the repository root (about a second):  Rscript R/statistical_tests.R
# Prints to the console and writes no file.
# Seeds: boot_ci() calls set.seed(42) before each bootstrap (R = 10,000
# resamples), so each bootstrap CI reproduces on its own; the other tests use
# no RNG.
#
# Author: Oliver Frings
# =============================================================================

OUTPUT_DIR <- "output"

# -----------------------------------------------------------------------------
# Helper: Wilson score confidence interval for a proportion
# -----------------------------------------------------------------------------
wilson_ci <- function(x, n, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  p_hat <- x / n
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
  c(lower = max(0, center - margin),
    estimate = p_hat,
    upper = min(1, center + margin))
}

# Format as percentage with CI
fmt_pct_ci <- function(x, n, digits = 1) {
  ci <- wilson_ci(x, n)
  sprintf("%.1f\\%% [%.1f, %.1f]", ci["estimate"] * 100,
          ci["lower"] * 100, ci["upper"] * 100)
}

# Bootstrap CI for a mean
boot_ci <- function(vals, R = 10000, alpha = 0.05, seed = 42) {
  set.seed(seed)
  n <- length(vals)
  boot_means <- replicate(R, mean(sample(vals, n, replace = TRUE)))
  q <- quantile(boot_means, c(alpha / 2, 1 - alpha / 2))
  c(lower = unname(q[1]), estimate = mean(vals), upper = unname(q[2]))
}

fmt_mean_ci <- function(vals, digits = 3) {
  ci <- boot_ci(vals)
  sprintf("%.3f [%.3f, %.3f]", ci["estimate"], ci["lower"], ci["upper"])
}

cat("=" , rep("=", 78), "\n", sep = "")
cat("  STATISTICAL ANALYSIS FOR MANUSCRIPT\n")
cat("=", rep("=", 78), "\n\n", sep = "")

# =============================================================================
# 1. MAIN RESULTS: CIs for overall accuracy (Table 1)
# =============================================================================
cat("--- TABLE 1: Overall criterion comparison ---\n\n")

main <- read.csv(file.path(OUTPUT_DIR, "main_results.csv"))
n_total <- nrow(main)

for (crit in c("bic", "aic", "icl")) {
  col <- paste0(crit, "_correct")
  x <- sum(main[[col]])
  ci <- wilson_ci(x, n_total)
  cat(sprintf("  %s accuracy: %d/%d = %.1f%% [%.1f, %.1f]\n",
              toupper(crit), x, n_total, ci["estimate"] * 100,
              ci["lower"] * 100, ci["upper"] * 100))
}

# McNemar test: BIC vs ICL
cat("\n  McNemar test (BIC vs ICL):\n")
tab_bi <- table(BIC = main$bic_correct, ICL = main$icl_correct)
cat("    Contingency table:\n")
print(tab_bi)
mcn <- mcnemar.test(tab_bi)
cat(sprintf("    chi-sq = %.2f, p = %.2e\n", mcn$statistic, mcn$p.value))

# McNemar test: BIC vs AIC
cat("\n  McNemar test (BIC vs AIC):\n")
tab_ba <- table(BIC = main$bic_correct, AIC = main$aic_correct)
cat("    Contingency table:\n")
print(tab_ba)
mcn2 <- mcnemar.test(tab_ba)
cat(sprintf("    chi-sq = %.2f, p = %.2e\n", mcn2$statistic, mcn2$p.value))

# =============================================================================
# 2. BIC accuracy by K (Table 2)
# =============================================================================
cat("\n--- TABLE 2: BIC accuracy by true K ---\n\n")

for (k in sort(unique(main$true_K))) {
  sub <- main[main$true_K == k, ]
  x <- sum(sub$bic_correct)
  n <- nrow(sub)
  ci <- wilson_ci(x, n)
  cat(sprintf("  K=%d: %d/%d = %.1f%% [%.1f, %.1f]\n",
              k, x, n, ci["estimate"] * 100,
              ci["lower"] * 100, ci["upper"] * 100))
}

# =============================================================================
# 3. BIC accuracy by kappa (Table 3)
# =============================================================================
cat("\n--- TABLE 3: BIC accuracy by separation (kappa) ---\n\n")

# K>=2 only (K=1 has kappa=0)
main_multi <- main[main$true_K >= 2, ]
for (kap in sort(unique(main_multi$kappa))) {
  sub <- main_multi[main_multi$kappa == kap, ]
  x <- sum(sub$bic_correct)
  n <- nrow(sub)
  ci <- wilson_ci(x, n)
  cat(sprintf("  kappa=%.2f: %d/%d = %.1f%% [%.1f, %.1f]\n",
              kap, x, n, ci["estimate"] * 100,
              ci["lower"] * 100, ci["upper"] * 100))
}

# =============================================================================
# 4. BIC accuracy by K x kappa (Table 4)
# =============================================================================
cat("\n--- TABLE 4: BIC accuracy by K x kappa ---\n\n")

for (k in c(2, 3, 4, 5)) {
  for (kap in sort(unique(main_multi$kappa))) {
    sub <- main_multi[main_multi$true_K == k & main_multi$kappa == kap, ]
    if (nrow(sub) == 0) next
    x <- sum(sub$bic_correct)
    n <- nrow(sub)
    ci <- wilson_ci(x, n)
    cat(sprintf("  K=%d, kappa=%.2f: %d/%d = %.1f%% [%.1f, %.1f]\n",
                k, kap, x, n, ci["estimate"] * 100,
                ci["lower"] * 100, ci["upper"] * 100))
  }
}

# =============================================================================
# 5. BIC accuracy by sigma (Table 5)
# =============================================================================
cat("\n--- TABLE 5: BIC accuracy by sigma ---\n\n")

for (sig in sort(unique(main$sigma))) {
  sub <- main[main$sigma == sig, ]
  x <- sum(sub$bic_correct)
  n <- nrow(sub)
  ci <- wilson_ci(x, n)
  cat(sprintf("  sigma=%.2f: %d/%d = %.1f%% [%.1f, %.1f]\n",
              sig, x, n, ci["estimate"] * 100,
              ci["lower"] * 100, ci["upper"] * 100))
}

# =============================================================================
# 6. LOGISTIC REGRESSION: bic_correct ~ true_K + kappa + sigma
# =============================================================================
cat("\n--- LOGISTIC REGRESSION: BIC success on design factors ---\n\n")

# For K>=2 only (K=1 always succeeds)
main_multi$true_K_fac <- factor(main_multi$true_K)
fit <- glm(bic_correct ~ true_K_fac + kappa + sigma,
           data = main_multi, family = binomial)
cat("  Logistic regression (K>=2 only):\n")
print(summary(fit))

# Odds ratios with CIs
cat("\n  Odds ratios (95% CI):\n")
or <- exp(cbind(OR = coef(fit), confint.default(fit)))
print(round(or, 3))

# =============================================================================
# 7. MMNL COMPARISON (Table 7)
# =============================================================================
cat("\n--- TABLE 7: MMNL comparison ---\n\n")

mmnl <- read.csv(file.path(OUTPUT_DIR, "mmnl_results.csv"))

# Discrete DGP
disc <- mmnl[mmnl$true_K >= 2, ]
n_disc <- nrow(disc)
x_lcmnl_disc <- sum(disc$bic_prefers == "LCMNL")
ci_disc <- wilson_ci(x_lcmnl_disc, n_disc)
cat(sprintf("  Discrete DGP (K>=2): LCMNL preferred %d/%d = %.1f%% [%.1f, %.1f]\n",
            x_lcmnl_disc, n_disc, ci_disc["estimate"] * 100,
            ci_disc["lower"] * 100, ci_disc["upper"] * 100))

# Binomial test for discrete DGP
bt <- binom.test(x_lcmnl_disc, n_disc, p = 1/3,
                 alternative = "greater")
cat(sprintf("  Binomial test (H0: p=1/3): p = %.2e\n", bt$p.value))

# Continuous DGP
cont <- mmnl[mmnl$true_K == 1, ]
n_cont <- nrow(cont)
x_lcmnl_cont <- sum(cont$bic_prefers == "LCMNL")
x_mmnl_cont  <- sum(cont$bic_prefers == "MMNL")
x_mnl_cont   <- sum(cont$bic_prefers == "MNL")
cat(sprintf("  Continuous DGP (K=1): MNL=%d, LCMNL=%d, MMNL=%d (n=%d)\n",
            x_mnl_cont, x_lcmnl_cont, x_mmnl_cont, n_cont))

# Multinomial exact test (chi-squared goodness of fit)
cat("  Chi-squared test (H0: uniform selection):\n")
observed <- c(MNL = x_mnl_cont, LCMNL = x_lcmnl_cont, MMNL = x_mmnl_cont)
ct <- chisq.test(observed)
cat(sprintf("    chi-sq = %.2f, df = %d, p = %.4f\n",
            ct$statistic, ct$parameter, ct$p.value))

# Fisher exact test on full 2x3 table
cat("\n  Fisher exact test (DGP x model preference):\n")
tab_mmnl <- table(
  DGP = ifelse(mmnl$true_K >= 2, "discrete", "continuous"),
  Preferred = mmnl$bic_prefers
)
print(tab_mmnl)
ft <- fisher.test(tab_mmnl)
cat(sprintf("  Fisher p = %.2e\n", ft$p.value))

# =============================================================================
# 8. PARAMETER RECOVERY (Table 8)
# =============================================================================
cat("\n--- TABLE 8: Parameter recovery ---\n\n")

recov <- read.csv(file.path(OUTPUT_DIR, "recovery_results.csv"))
cat(sprintf("  RMSE = %.4f\n", recov$rmse))
cat(sprintf("  Bias = %.4f\n", recov$bias))

# Bootstrap CIs from individual conditions where bic_correct & K>=2
recov_conds <- main[main$bic_correct == 1 & main$true_K >= 2, ]
rmse_vals <- recov_conds$rmse[!is.na(recov_conds$rmse)]
bias_vals <- recov_conds$bias[!is.na(recov_conds$bias)]

if (length(rmse_vals) > 0) {
  ci_rmse <- boot_ci(rmse_vals)
  cat(sprintf("  RMSE: %.4f [%.4f, %.4f] (n=%d conditions)\n",
              ci_rmse["estimate"], ci_rmse["lower"], ci_rmse["upper"],
              length(rmse_vals)))
  ci_bias <- boot_ci(bias_vals)
  cat(sprintf("  Bias: %.4f [%.4f, %.4f] (n=%d conditions)\n",
              ci_bias["estimate"], ci_bias["lower"], ci_bias["upper"],
              length(bias_vals)))
}

# =============================================================================
# 9. ARI BY K (Table 9)
# =============================================================================
cat("\n--- TABLE 9: ARI by true K ---\n\n")

ari_conds <- main[!is.na(main$ari) & main$true_K >= 2, ]

for (k in sort(unique(ari_conds$true_K))) {
  vals <- ari_conds$ari[ari_conds$true_K == k]
  if (length(vals) > 1) {
    ci <- boot_ci(vals)
    cat(sprintf("  K=%d: mean=%.3f [%.3f, %.3f] (n=%d)\n",
                k, ci["estimate"], ci["lower"], ci["upper"], length(vals)))
  }
}

all_ari <- ari_conds$ari
ci_all <- boot_ci(all_ari)
cat(sprintf("  Overall: mean=%.3f [%.3f, %.3f] (n=%d)\n",
            ci_all["estimate"], ci_all["lower"], ci_all["upper"],
            length(all_ari)))

# =============================================================================
# 10. CONVERGENCE (Table 10)
# =============================================================================
cat("\n--- TABLE 10: Convergence comparison ---\n\n")

conv <- read.csv(file.path(OUTPUT_DIR, "convergence_results.csv"))

# Overall clustering success rate
n_conv <- nrow(conv)
x_cluster_best <- sum(conv$cluster_best)
ci_clust <- wilson_ci(x_cluster_best, n_conv)
cat(sprintf("  Clustering at global: %d/%d = %.1f%% [%.1f, %.1f]\n",
            x_cluster_best, n_conv,
            ci_clust["estimate"] * 100,
            ci_clust["lower"] * 100,
            ci_clust["upper"] * 100))

# Mean pct_global of random starts
ci_pctg <- boot_ci(conv$pct_global)
cat(sprintf("  Random starts at global: mean=%.1f%% [%.1f, %.1f]\n",
            ci_pctg["estimate"] * 100,
            ci_pctg["lower"] * 100,
            ci_pctg["upper"] * 100))

# By K
for (k in c(3, 4, 5)) {
  sub <- conv[conv$K == k, ]
  n_k <- nrow(sub)
  x_k <- sum(sub$cluster_best)
  ci_k <- wilson_ci(x_k, n_k)
  cat(sprintf("  K=%d: clustering at global %d/%d = %.1f%% [%.1f, %.1f]\n",
              k, x_k, n_k, ci_k["estimate"] * 100,
              ci_k["lower"] * 100, ci_k["upper"] * 100))

  ci_rand_k <- boot_ci(sub$pct_global)
  cat(sprintf("       random at global: mean=%.1f%% [%.1f, %.1f]\n",
              ci_rand_k["estimate"] * 100,
              ci_rand_k["lower"] * 100,
              ci_rand_k["upper"] * 100))
}

# Wilcoxon signed-rank test: cluster_best (1/0) vs pct_global
# Each condition: did clustering reach global? (binary) vs what fraction of
# random starts did?
# Better: compare via paired test on whether clustering beats random starts
cat("\n  Wilcoxon signed-rank test (pct_global vs 50% baseline):\n")
wt <- wilcox.test(conv$pct_global, mu = 0.5, alternative = "less")
cat(sprintf("    V = %.0f, p = %.2e\n", wt$statistic, wt$p.value))

# =============================================================================
# 11. SUPPLEMENTARY: Unbalanced (Table 11)
# =============================================================================
cat("\n--- TABLE 11: Unbalanced class proportions ---\n\n")

unbal <- read.csv(file.path(OUTPUT_DIR, "unbalanced_results.csv"))
for (i in 1:nrow(unbal)) {
  # These are aggregate accuracy percentages from 15 conditions each
  cat(sprintf("  %s: %.1f%%\n", unbal$config[i], unbal$accuracy[i]))
}
cat("  Note: Each based on 15 conditions; Wilson CIs not computable from aggregate.\n")
cat("  Need per-condition data for CIs.\n")

# =============================================================================
# 12. SUPPLEMENTARY: Design (Table 12)
# =============================================================================
cat("\n--- TABLE 12: Design type ---\n\n")

design <- read.csv(file.path(OUTPUT_DIR, "design_results.csv"))
cat(sprintf("  Random: %.1f%%, D-efficient: %.1f%%\n",
            design$random, design$deff))

# =============================================================================
# 13. SUPPLEMENTARY: Concomitant (Table 13)
# =============================================================================
cat("\n--- TABLE 13: Concomitant variables ---\n\n")

conc <- read.csv(file.path(OUTPUT_DIR, "concomitant_results.csv"))
cat(sprintf("  Accuracy: %.1f%%\n", conc$accuracy))
cat(sprintf("  Mean ARI: %.3f\n", conc$mean_ari))

# =============================================================================
# 14. SUPPLEMENTARY: Sample sensitivity (Table 14)
# =============================================================================
cat("\n--- TABLE 14: Sample sensitivity ---\n\n")

samp <- read.csv(file.path(OUTPUT_DIR, "sample_sensitivity_results.csv"))
for (i in 1:nrow(samp)) {
  cat(sprintf("  T=%d, N=%d: %.1f%%\n",
              samp$T_tasks[i], samp$N_per_class[i], samp$accuracy[i]))
}

# =============================================================================
# 15. CORRELATED MMNL (Table 15)
# =============================================================================
cat("\n--- TABLE 15: Correlated MMNL ---\n\n")

corr <- read.csv(file.path(OUTPUT_DIR, "mmnl_correlated_results.csv"))

# By DGP
disc_c <- corr[corr$true_K >= 2, ]
cont_c <- corr[corr$true_K == 1, ]

x_lc_dc <- sum(disc_c$bic_prefers == "LCMNL")
n_dc <- nrow(disc_c)
ci_lc_dc <- wilson_ci(x_lc_dc, n_dc)
cat(sprintf("  Discrete DGP: LCMNL preferred %d/%d = %.1f%% [%.1f, %.1f]\n",
            x_lc_dc, n_dc, ci_lc_dc["estimate"] * 100,
            ci_lc_dc["lower"] * 100, ci_lc_dc["upper"] * 100))

cat(sprintf("  Continuous DGP (n=%d):\n", nrow(cont_c)))
cat(sprintf("    LCMNL: %d, MMNL_indep: %d, MMNL_corr: %d\n",
            sum(cont_c$bic_prefers == "LCMNL"),
            sum(cont_c$bic_prefers == "MMNL_indep"),
            sum(cont_c$bic_prefers == "MMNL_corr")))

# Binomial test: discrete DGP
bt2 <- binom.test(x_lc_dc, n_dc, p = 1/3, alternative = "greater")
cat(sprintf("  Binomial test (discrete, H0: p=1/3): p = %.2e\n", bt2$p.value))

# =============================================================================
# SUMMARY OF KEY STATISTICS FOR MANUSCRIPT
# =============================================================================
cat("\n")
cat("=", rep("=", 78), "\n", sep = "")
cat("  KEY STATISTICS READY FOR MANUSCRIPT\n")
cat("=", rep("=", 78), "\n\n")

# Overall BIC
x_bic <- sum(main$bic_correct)
ci_bic <- wilson_ci(x_bic, n_total)
cat(sprintf("BIC overall: %.1f%% [%.1f, %.1f], n=%d\n",
            ci_bic["estimate"] * 100, ci_bic["lower"] * 100,
            ci_bic["upper"] * 100, n_total))

# BIC for K<=3
k123 <- main[main$true_K <= 3, ]
x_k123 <- sum(k123$bic_correct)
n_k123 <- nrow(k123)
ci_k123 <- wilson_ci(x_k123, n_k123)
cat(sprintf("BIC for K<=3: %.1f%% [%.1f, %.1f], n=%d\n",
            ci_k123["estimate"] * 100, ci_k123["lower"] * 100,
            ci_k123["upper"] * 100, n_k123))

cat("\nDone.\n")
