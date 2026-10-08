#!/usr/bin/env Rscript
# =============================================================================
# H4 beyond the extremes, blocked design (2026-09-26).
#
# The published H4 test pits a single normal (the MMNL's own form) against
# well-separated segments (kappa >= 0.75). This arm adds the two harder cases a
# referee will ask about, on the same shared 48-card blocked D-efficient design
# (T = 12, N = 300):
#   skew   K* = 1, continuous but right-skewed tastes: each generic coefficient is
#          beta_bar + sigma * u with u a standardised lognormal (log-sd 0.75,
#          skewness ~2.9, mean 0, sd 1); price as in the main DGP
#          (normal, sd 0.3*sigma, truncated at -0.1). sigma in {0.25, 0.40}.
#          Correct reading: continuous (MMNL or MNL), not segments.
#   overlap K* = 2 with overlapping segments: kappa = 0.5, sigma in {0.30, 0.40}.
#          Correct reading: discrete (LCMNL).
# Per replication: LCMNL C = 1..5 (C = 1 is the MNL) and the independent MMNL
# (normals + lognormal price, one Apollo estimation at 3,000 MLHS); winner = lowest BIC.
# 10 replications per cell, 40 conditions.
#
# Resumable (skips cells already in the CSV).
# Seeds: data seed
#   as.integer(3e4 + 1000*(cell == "overlap") + 100*(sigma*100) + rep),
# kept in the CSV's seed column; the shared design
# klue_design(n_cards = 48, n_blocks = 4) under set.seed(20240601). A skew
# replication takes its block assignment from klue_simulate under that seed,
# then draws its skewed tastes and its choices after set.seed(seed + 1).
# klue's clustering starts use seed 123; Apollo makes the MLHS draws under its
# default apollo_control$seed (13). The formula gives skew sigma = 0.40 and
# overlap sigma = 0.30 the same seeds (34001-34010); their data still differ.
# Run from the repository root, one core (dev/rerun_standard_apollo.sh adds
# Apollo's memorySaver):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   nohup nice -n 15 Rscript dev/run_h4_misspec_blocked.R \
#     > output/h4_misspec_blocked.log 2>&1 &
# Writes: output/h4_misspec_blocked.csv
# =============================================================================
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "idefix", "mclust", "cluster"))  # stops on drift
options(klue.cores = 1L, klue.mmnl.n_cores = 1L)
OUT <- "output"; CSV <- file.path(OUT, "h4_misspec_blocked.csv")
stamp <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%m-%d %H:%M:%S"), m)); flush.console() }

DGP    <- klue_dgp(n_generic = 4, n_alternatives = 3)
DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)   # the main blocked design
LOGSD  <- 0.75
std_lognormal <- function(n) {
  m <- exp(LOGSD^2 / 2); s <- sqrt((exp(LOGSD^2) - 1) * exp(LOGSD^2))
  (exp(LOGSD * rnorm(n)) - m) / s
}

conds <- rbind(
  expand.grid(cell = "skew",    true_K = 1L, kappa = 0,   sigma = c(0.25, 0.40), rep = 1:10,
              stringsAsFactors = FALSE),
  expand.grid(cell = "overlap", true_K = 2L, kappa = 0.5, sigma = c(0.30, 0.40), rep = 1:10,
              stringsAsFactors = FALSE))
simulate_cell <- function(cell, sig, seed) {
  if (cell == "overlap")
    return(klue_simulate(N_per_class = 150, true_K = 2, separation = 0.5,
                         heterogeneity = sig, seed = seed, dgp = DGP,
                         design = DESIGN)$database)
  # skew: attribute layout from the blocked design, then skewed tastes
  db <- klue_simulate(N_per_class = 300, true_K = 1, separation = 0,
                      heterogeneity = 0, seed = seed, dgp = DGP,
                      design = DESIGN)$database
  set.seed(seed + 1L)
  N  <- max(db$ID); nb <- DGP$n_beta; ng <- DGP$n_generic
  ib <- matrix(DGP$beta_bar, nrow = N, ncol = nb, byrow = TRUE)
  for (a in seq_len(ng)) ib[, a] <- ib[, a] + sig * std_lognormal(N)
  ib[, nb] <- pmin(rnorm(N, DGP$beta_bar[nb], 0.3 * sig), -0.1)
  klue_simulate_choices(db, ib, DGP)
}

# Main loop only when run as a script, so the functions above can be sourced
# for testing.
if (sys.nframe() == 0L) {
rows <- if (file.exists(CSV)) read.csv(CSV, stringsAsFactors = FALSE) else NULL
done <- if (is.null(rows)) character(0) else paste(rows$cell, rows$sigma, rows$rep)
stamp(sprintf("H4 misspecification arm (blocked): %d conditions", nrow(conds)))
for (i in seq_len(nrow(conds))) {
  cd <- conds[i, ]
  key <- paste(cd$cell, cd$sigma, cd$rep)
  if (key %in% done) { stamp(paste("SKIP", key)); next }
  seed <- as.integer(3e4 + 1000 * (cd$cell == "overlap") + 100 * (cd$sigma * 100) + cd$rep)
  db <- simulate_cell(cd$cell, cd$sigma, seed)
  bic <- rep(Inf, 5)
  for (C in 1:5) {
    m <- tryCatch(klue_lcmnl(db, C, dgp = DGP), error = function(e) NULL)
    if (!is.null(m) && isTRUE(m$converged)) bic[C] <- m$BIC
  }
  lc_C <- (2:5)[which.min(bic[2:5])]
  mm <- tryCatch(klue_mmnl(db, n_cores = 1L, dgp = DGP), error = function(e) NULL)
  mm_bic <- if (!is.null(mm) && isTRUE(mm$converged)) mm$BIC else Inf
  bics <- c(MNL = bic[1], LCMNL = min(bic[2:5]), MMNL = mm_bic)
  row <- data.frame(cell = cd$cell, true_K = cd$true_K, kappa = cd$kappa,
                    sigma = cd$sigma, rep = cd$rep, seed = seed,
                    mnl_bic = bic[1], lc_bic = min(bic[2:5]), lc_C = lc_C,
                    mmnl_bic = mm_bic, winner = names(which.min(bics)),
                    lcmnl_overselects = min(bic[2:5]) < bic[1],
                    stringsAsFactors = FALSE)
  rows <- rbind(rows, row)
  write.csv(rows, CSV, row.names = FALSE)
  stamp(sprintf("DONE %s -> %s (MNL=%.0f LCMNL[C%d]=%.0f MMNL=%.0f)", key, row$winner,
                bic[1], lc_C, row$lc_bic, mm_bic))
}
stamp("ALL DONE")
print(with(rows, table(paste(cell, sigma), winner)))
}
