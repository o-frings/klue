# =============================================================================
# Gate 3 (H5): Parameter recovery (RMSE, bias) and classification quality (ARI),
# aggregated per design regime.
#
# Reproduces the randomized-regime tab:recovery / tab:ari numbers (verification)
# and produces the blocked-regime panels required by Path A of the revision.
#
# Recovery (RMSE, bias): multi-class, BIC-correct conditions (the paper's N=357
#   for the randomized regime). ARI: all multi-class conditions (the paper's
#   N=400). Bootstrap 95% CIs (R=10,000); boot_ci() calls set.seed(42) before
#   each bootstrap.
#
# Light job -- no compute window needed. Reads output/main_results.csv and
# output/main_blocked.csv.
#   Run (repository root):  Rscript dev/aggregate_recovery.R
# Writes: output/recovery_blocked.csv, output/ari_blocked.csv
# =============================================================================

OUTPUT_DIR <- "output"
BOOT_R <- 10000L
SEED   <- 42L

boot_ci <- function(vals, R = BOOT_R, alpha = 0.05, seed = SEED) {
  vals <- vals[is.finite(vals)]
  set.seed(seed)
  bm <- replicate(R, mean(sample(vals, length(vals), replace = TRUE)))
  q  <- quantile(bm, c(alpha / 2, 1 - alpha / 2))
  c(estimate = mean(vals), lower = unname(q[1]), upper = unname(q[2]),
    n = length(vals))
}

summarise_regime <- function(csv, label) {
  d <- read.csv(file.path(OUTPUT_DIR, csv))

  # Recovery: multi-class & BIC-correct (where parameter recovery is defined)
  rec <- d[d$true_K >= 2 & d$bic_correct == 1 &
             is.finite(d$rmse) & is.finite(d$bias), ]
  rmse_ci <- boot_ci(rec$rmse)
  bias_ci <- boot_ci(rec$bias)

  # ARI: all multi-class conditions
  ad <- d[d$true_K >= 2 & is.finite(d$ari), ]
  ari_all <- boot_ci(ad$ari)
  ks <- sort(unique(ad$true_K))
  ari_tab <- do.call(rbind, lapply(ks, function(k) {
    ci <- boot_ci(ad$ari[ad$true_K == k])
    data.frame(group = paste0("K", k), ari = ci["estimate"],
               ari_lo = ci["lower"], ari_hi = ci["upper"], n = ci["n"])
  }))
  ari_tab <- rbind(ari_tab,
                   data.frame(group = "all", ari = ari_all["estimate"],
                              ari_lo = ari_all["lower"], ari_hi = ari_all["upper"],
                              n = ari_all["n"]))
  ari_tab$regime <- label
  rownames(ari_tab) <- NULL

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("recovery N=%d  RMSE=%.4f [%.4f, %.4f]  bias=%.4f [%.4f, %.4f]\n",
              nrow(rec), rmse_ci["estimate"], rmse_ci["lower"], rmse_ci["upper"],
              bias_ci["estimate"], bias_ci["lower"], bias_ci["upper"]))
  print(cbind(group = ari_tab$group,
              round(ari_tab[, c("ari", "ari_lo", "ari_hi", "n")], 3)),
        row.names = FALSE)

  rec_row <- data.frame(regime = label, N = nrow(rec),
                        rmse = rmse_ci["estimate"], rmse_lo = rmse_ci["lower"],
                        rmse_hi = rmse_ci["upper"], bias = bias_ci["estimate"],
                        bias_lo = bias_ci["lower"], bias_hi = bias_ci["upper"],
                        row.names = NULL)
  list(recovery = rec_row, ari = ari_tab)
}

o <- summarise_regime("main_results.csv", "randomized")  # verify vs paper
b <- summarise_regime("main_blocked.csv",  "blocked")     # new Path A panels

write.csv(b$recovery, file.path(OUTPUT_DIR, "recovery_blocked.csv"), row.names = FALSE)
write.csv(b$ari[, c("regime", "group", "ari", "ari_lo", "ari_hi", "n")],
          file.path(OUTPUT_DIR, "ari_blocked.csv"), row.names = FALSE)
cat("\nWrote output/recovery_blocked.csv and output/ari_blocked.csv\n")
