#!/usr/bin/env Rscript
# =============================================================================
# Swissmetro class profiles at C = 4 (ICL minimum) and C = 5 (runner-up; the
# class count of Sfeir et al. 2022) -- the substantive output behind the
# empirical section (2026-09-26).
#
# For each class: share, coefficients with respondent-clustered (robust) SEs,
# and the value of travel time VTT = b_time / b_cost in CHF per hour (time and
# cost both enter per 100 units, so the ratio is CHF per minute, times 60),
# with a delta-method SE. Flags classes whose time or cost coefficient has the
# wrong sign. Seeded multistart (klue_lcmnl), canonical sample as in
# R/empirical_swissmetro.R.
#
# Run: KLUE_CORES=1 nice -n 15 Rscript dev/swissmetro_class_profiles.R
# Writes: output/swissmetro_class_profiles.csv
# =============================================================================
options(klue.cores = 1L, sm.skip = TRUE)
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(source("R/empirical_swissmetro.R")))

db <- load_swissmetro_database()
rows <- list()
for (C in 4:5) {
  m <- klue_lcmnl(db, C, dgp = SM_DGP)
  se <- sqrt(diag(m$robust_vcov)); est <- m$par
  for (cl in seq_len(C)) {
    nm <- function(k) sprintf("b%d_class%d", k, cl)
    vtt <- klue_wtp(m, num = setNames(-1, nm(1)), den = setNames(1, nm(3)), scale = 60)
    rows[[length(rows) + 1]] <- data.frame(
      C = C, class = cl, share = m$class_probs[cl],
      b_time = est[nm(1)], se_time = se[nm(1)],
      b_headway = est[nm(2)], se_headway = se[nm(2)],
      b_cost = est[nm(3)], se_cost = se[nm(3)],
      asc_train = est[sprintf("asc1_class%d", cl)], asc_sm = est[sprintf("asc2_class%d", cl)],
      vtt_chf_h = vtt[["wtp"]], vtt_se = vtt[["se"]],
      wrong_sign = est[nm(1)] > 0 || est[nm(3)] > 0,
      LL = m$LL, BIC = m$BIC, row.names = NULL)
  }
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/swissmetro_class_profiles.csv", row.names = FALSE)
print(within(tab, { share <- round(share, 3) }), digits = 3, row.names = FALSE)
