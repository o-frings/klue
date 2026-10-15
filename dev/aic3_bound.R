# =============================================================================
# BOUNDS ON AIC3 SELECTION ACCURACY FROM THE STORED AIC AND BIC CHOICES
#
# Andrews & Currim (2003, JMR 40(2)) found AIC3 = -2LL + 3k the most accurate
# segment-retention criterion for finite-mixture logit models on simulated
# scanner panels. The main simulation studies stored only each criterion's
# selected C, not the per-C log-likelihoods, so AIC3 cannot be recomputed
# without refitting. It can be bounded exactly.
#
# Argument: within a condition every criterion ranks the same fitted models,
# and k increases with C. For penalties lambda1 < lambda2 the minimiser of
# -2LL_C + lambda k_C is weakly decreasing in lambda (add the two optimality
# inequalities: (lambda2 - lambda1)(k2 - k1) <= 0). AIC's penalty is 2,
# AIC3's 3 and BIC's log(N) >= log(300) = 5.70, so
#   C_BIC <= C_AIC3 <= C_AIC   in every condition.
# Hence AIC3 can be correct only where C_BIC <= K* <= C_AIC (upper bound),
# is certainly correct where C_BIC = C_AIC = K* (lower bound), and
# over-selects wherever BIC does.
#
# Inputs: output/main_blocked.csv (1,680 blocked conditions) and
#   output/main_results.csv (420 randomised conditions).
# Writes: output/aic3_bound.csv. Deterministic; no random numbers.
# Run from the repository root: Rscript dev/aic3_bound.R
# =============================================================================

if (!exists("OUTPUT_DIR")) OUTPUT_DIR <- "output"

bound_one <- function(file, design) {
  d <- read.csv(file.path(OUTPUT_DIR, file))
  stopifnot(all(d$selected_bic <= d$selected_aic))   # monotonicity holds
  K <- d$true_K; b <- d$selected_bic; a <- d$selected_aic
  data.frame(
    design = design, n = nrow(d),
    bic_correct = mean(b == K), aic_correct = mean(a == K),
    aic3_correct_lower = mean(b == K & a == K),
    aic3_correct_upper = mean(b <= K & K <= a),
    bic_over = mean(b > K),                 # AIC3 over-selects at least here
    k1_bic_over = mean(b[K == 1] > 1),      # continuous K* = 1 conditions
    k1_n = sum(K == 1))
}

out <- rbind(bound_one("main_blocked.csv", "blocked"),
             bound_one("main_results.csv", "randomised"))
print(out, row.names = FALSE, digits = 4)
write.csv(out, file.path(OUTPUT_DIR, "aic3_bound.csv"), row.names = FALSE)
