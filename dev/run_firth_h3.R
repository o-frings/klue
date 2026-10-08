# =============================================================================
# Gate 2 + Gate 3 (H3): Firth-penalised logistic regression of correct BIC
# enumeration on the design factors.
#
# Replaces the ordinary-ML glm() in R/statistical_tests.R, whose Wald SEs
# are inflated by quasi-complete separation (100% accuracy in several cells).
# Firth penalisation (Heinze & Schemper 2002) gives finite estimates and
# profile-likelihood (PL) CIs that are valid under separation.
#
# RHS matches version10.tex (eq. near line 260):
#   logit Pr(correct) = a + sum_{j=3}^{5} gamma_j * 1[K*=j] + delta*kappa + eta*sigma
# Estimated on the multi-class conditions (K* >= 2; K*=1 always succeeds).
#
# Light job -- no compute window needed. Reads output/main_results.csv and
# output/main_blocked.csv. Run from the repository root:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/run_firth_h3.R
# Writes: output/firth_h3_orthogonal.csv  (randomized / high-diversity regime)
#         output/firth_h3_blocked.csv     (blocked D-efficient regime)
# Deterministic; no RNG (glm and logistf, with profile-likelihood CIs).
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("logistf", quietly = TRUE))
    stop("Package 'logistf' is required: install.packages('logistf')")
  library(logistf)
})

OUTPUT_DIR <- "output"
FORM <- bic_correct ~ true_K_fac + kappa + sigma

fit_h3 <- function(csv, label) {
  d <- read.csv(file.path(OUTPUT_DIR, csv))
  d <- d[d$true_K >= 2, ]            # H3 is about the multi-class conditions
  d$true_K_fac <- factor(d$true_K)

  # (1) Ordinary ML logit -- the current model (Wald CIs, separation-inflated)
  g    <- glm(FORM, data = d, family = binomial)
  gsum <- summary(g)$coefficients
  gci  <- confint.default(g)        # Wald
  trm  <- rownames(gsum)

  # (2) Firth-penalised logit with profile-likelihood CIs.
  # Raise maxit: under separation the PL limits for some terms need many iters.
  f <- logistf(FORM, data = d,
               control   = logistf.control(maxit = 10000),
               plcontrol = logistpl.control(maxit = 10000))

  out <- data.frame(
    regime          = label,
    term            = trm,
    coef_ml         = gsum[trm, "Estimate"],
    se_ml           = gsum[trm, "Std. Error"],
    ci_lower_ml     = gci[trm, 1],
    ci_upper_ml     = gci[trm, 2],
    coef_firth      = as.numeric(coef(f))[match(trm, names(coef(f)))],
    se_firth        = sqrt(diag(f$var))[match(trm, names(coef(f)))],
    ci_lower_firth  = f$ci.lower[match(trm, names(coef(f)))],
    ci_upper_firth  = f$ci.upper[match(trm, names(coef(f)))],
    p_firth         = f$prob[match(trm, names(coef(f)))],
    row.names = NULL
  )
  out$or_firth       <- exp(out$coef_firth)
  out$or_lower_firth <- exp(out$ci_lower_firth)
  out$or_upper_firth <- exp(out$ci_upper_firth)

  cat(sprintf("\n===== H3 logistic regression: %s  (N=%d, K* >= 2) =====\n",
              label, nrow(d)))
  disp_cols <- c("coef_ml","se_ml","coef_firth","se_firth",
                 "ci_lower_firth","ci_upper_firth","p_firth")
  print(cbind(term = out$term, round(out[, disp_cols], 4)), row.names = FALSE)
  sep <- which(out$se_ml > 10)
  if (length(sep))
    cat(sprintf("  separation symptom: ML Wald SE > 10 for %s (Firth finite)\n",
                paste(out$term[sep], collapse = ", ")))
  out
}

res_o <- fit_h3("main_results.csv", "randomized")   # high-diversity regime (was "orthogonal")
res_b <- fit_h3("main_blocked.csv", "blocked")

write.csv(res_o, file.path(OUTPUT_DIR, "firth_h3_orthogonal.csv"), row.names = FALSE)
write.csv(res_b, file.path(OUTPUT_DIR, "firth_h3_blocked.csv"),    row.names = FALSE)
cat("\nWrote output/firth_h3_orthogonal.csv and output/firth_h3_blocked.csv\n")

# ---- Gate 2 bar: does the sigma null result survive under Firth? ----
cat("\n--- sigma (within-class heterogeneity) under Firth ---\n")
for (r in list(res_o, res_b)) {
  s <- r[r$term == "sigma", ]
  cat(sprintf("[%-10s] coef=%+.3f  PL 95%% CI [%+.3f, %+.3f]  p=%.3f  -> %s\n",
              s$regime, s$coef_firth, s$ci_lower_firth, s$ci_upper_firth, s$p_firth,
              ifelse(s$p_firth > 0.05, "null (claim holds)", "SIGNIFICANT (revise)")))
}
