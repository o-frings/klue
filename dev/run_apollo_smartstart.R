# =============================================================================
# Gate 1: does apollo_searchStart with smartStart = TRUE (NOT Apollo's default;
# the 0.3.5 default is smartStart = FALSE) converge on the blocked, low-diversity
# design once attributes are rescaled?
#
# Historical diagnostic for the non-default arms published in version12. The
# smartStart failure is a numerical overflow in Apollo 0.3.5: its candidate
# weights exp(0.5 * eigenvalue) are Inf when the log-likelihood Hessian at the
# start has eigenvalues above about 1400, so sample() gets NaN probabilities.
# The stress ladders now run apollo_searchStart at its defaults.
# This script, on the blocked stress rungs, runs:
#   (a) smartStart = TRUE on RAW attributes      -> reproduce the failure;
#   (b) smartStart = TRUE on RESCALED attributes -> the proposed fix;
#   (c) smartStart = FALSE bounded [-3,3]         -> Apollo's reduced search
#                                                    (apollo_ss_published);
# and compares each to the clustering+ML global LL. Utility is bilinear, so
# standardising attributes leaves the optimum LL unchanged (only the coefficient
# scale changes), making (b) directly comparable to the global.
#
# Needs apollo (installed); does NOT need gmnl/mlogit. Moderate compute.
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/run_apollo_smartstart.R
# Writes: output/apollo_smartstart_diag.csv (its bounded_secs column holds the
#         times of Apollo's reduced search that the paper reports)
# Seeds: the shared design klue_design(n_cards = 48, n_blocks = 4) runs under
# set.seed(20240601); each rung uses seed 1, data seed
# as.integer(1000*K + 100*(kap*100) + 1), so the data are the seed-1 cells of
# dev/stress_replicate_blocked.R; klue's clustering starts use seed 123;
# apollo_searchStart draws its candidates after set.seed(17) (Apollo's default
# seed 13, plus 4) in all three Apollo runs. LLs reproduce; seconds do not.
# =============================================================================

if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else suppressWarnings(suppressMessages(library(klue)))
sys.source("dev/compare_packages.R", envir = globalenv())   # method runners (guarded)
options(mc.cores = 1L, klue.mmnl.n_cores = 1L)

DGP    <- klue_dgp(n_generic = 4, n_alternatives = 3)
DESIGN <- klue_design(n_cards = 48L, n_blocks = 4L, dgp = DGP)

# Standardise each attribute (pooled across alternatives) to ~unit scale so the
# searchStart Hessian is well-conditioned. The optimum LL is invariant.
rescale_db <- function(db, dgp) {
  J <- dgp$n_alternatives; ng <- dgp$n_generic
  std <- function(cols) {
    v <- unlist(db[, cols]); m <- mean(v); s <- stats::sd(v)
    if (s > 0) for (cc in cols) db[[cc]] <<- (db[[cc]] - m) / s
  }
  for (a in 1:ng) std(paste0("x", a, "_", 1:J))
  std(paste0("price_", 1:J))
  db
}

RUNGS <- list(
  list(label = "moderate",  K = 4, kap = 0.75, sig = 0.25),
  list(label = "hard",      K = 5, kap = 0.50, sig = 0.30),
  list(label = "very_hard", K = 5, kap = 0.30, sig = 0.35)
)

rows <- list()
for (r in RUNGS) {
  seed <- as.integer(1000 * r$K + 100 * (r$kap * 100) + 1)
  d  <- klue_simulate(N_per_class = 60, T_tasks = 12, true_K = r$K,
                      separation = r$kap, heterogeneity = r$sig,
                      seed = seed, design = DESIGN)
  db  <- d$database
  dbs <- rescale_db(db, DGP)

  glob <- run_klue(db, r$K, DGP, "ml")$LL                       # clustering+ML global
  raw  <- run_apollo_ss_published_randomised(db,  r$K, DGP) # smartStart=TRUE, raw attrs
  scl  <- run_apollo_ss_published_randomised(dbs, r$K, DGP) # smartStart=TRUE, rescaled
  bnd  <- run_apollo_ss_published_blocked(db,  r$K, DGP)  # smartStart=FALSE bounded (paper)

  rows[[length(rows) + 1L]] <- data.frame(
    rung = r$label, K = r$K, global_LL = round(glob, 2),
    ss_raw_LL = round(raw$LL, 2), ss_raw_note = raw$note,
    ss_scaled_LL = round(scl$LL, 2), ss_scaled_note = scl$note,
    ss_scaled_secs = round(scl$seconds),
    bounded_LL = round(bnd$LL, 2), bounded_secs = round(bnd$seconds),
    stringsAsFactors = FALSE)
  cat(sprintf("[%-9s K=%d] global=%.1f | ss_raw=%s | ss_scaled=%s (%.0fs) | bounded=%.1f\n",
              r$label, r$K, glob,
              ifelse(is.finite(raw$LL), sprintf("%.1f", raw$LL), paste0("NA[", raw$note, "]")),
              ifelse(is.finite(scl$LL), sprintf("%.1f", scl$LL), paste0("NA[", scl$note, "]")),
              scl$seconds, bnd$LL))
  flush.console()
}
out <- do.call(rbind, rows)
write.csv(out, "output/apollo_smartstart_diag.csv", row.names = FALSE)
cat("\nwritten: output/apollo_smartstart_diag.csv\n"); print(out, row.names = FALSE)
