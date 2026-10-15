# =============================================================================
# CROSS-CHECK AGAINST HESS & BUNCH (arXiv:2506.02722v2, 2026), SWISS VTT DATA
#
# Hess & Bunch fit a two-class LCMNL to apollo_swissRouteChoiceData (388
# respondents x 9 tasks, J = 2) with four linear level-of-service attributes
# (travel time, travel cost, headway, interchanges), no constant and
# constant-only class shares. Across eight estimation strategies they find
# eight distinct optima; the best is LL = -1,552.53 (their Tables 3 and 5).
# This script fits the same model with klue from the six clustering starts
# and reports where each start ends, so the paper can say whether the
# clustering starts reach their best solution.
#
# Specification: klue_dgp(n_generic = 3, n_alternatives = 2, asc_map = 0:0),
# i.e. no constant; attribute scaling as in R/empirical_swiss_route.R
# (time and headway / 60, cost / 10), which leaves the likelihood unchanged.
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/check_hess_bunch_swiss.R
# Writes: output/hess_bunch_swiss_check.csv (one row per clustering start,
#   full-precision fit, plus the default screened multistart).
# Seeds: none set here. klue_lcmnl() runs its k-means, Gaussian-mixture and
# PAM clustering starts under seed 123 and restores the caller's stream.
# =============================================================================

options(swiss.skip = TRUE)
source("R/empirical_swiss_route.R")   # loads klue, defines load_swiss_database()

HB_BEST_LL <- -1552.53                # Hess & Bunch's best of eight optima
database <- load_swiss_database()
dgp_noasc <- klue_dgp(n_generic = 3, n_alternatives = 2, asc_map = c(0L, 0L))

# Every start at full precision, so each start's own endpoint is visible.
full <- klue_lcmnl(database, 2L, dgp = dgp_noasc, screen = FALSE)
rows <- lapply(names(full$method_results), function(nm) {
  m <- full$method_results[[nm]]
  data.frame(run = "full_precision", start = nm,
             converged = isTRUE(m$converged), LL = m$LL, k = m$k,
             gap_to_hess_bunch = m$LL - HB_BEST_LL)
})
# The default path (screened at a loose tolerance, winner polished).
dflt <- klue_lcmnl(database, 2L, dgp = dgp_noasc)
rows[[length(rows) + 1L]] <- data.frame(
  run = "default_screened", start = dflt$best_method,
  converged = isTRUE(dflt$converged), LL = dflt$LL, k = dflt$k,
  gap_to_hess_bunch = dflt$LL - HB_BEST_LL)

out <- do.call(rbind, rows)
print(out, row.names = FALSE, digits = 8)
cat(sprintf("\nClass shares (default fit): %s\n",
            paste(round(dflt$class_probs, 4), collapse = ", ")))
write.csv(out, file.path(OUTPUT_DIR, "hess_bunch_swiss_check.csv"),
          row.names = FALSE)
