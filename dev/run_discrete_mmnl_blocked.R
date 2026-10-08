# dev/run_discrete_mmnl_blocked.R
# Gap B: LCMNL-vs-MMNL on DISCRETE DGPs (K* in {2,3,4}) under the blocked
# D-efficient baseline -- the realistic-design counterpart of the orthogonal
# discrete result (tab:mmnl). Confirms LCMNL is preferred when heterogeneity is
# genuinely discrete, under the realistic design. (The continuous K*=1 case is
# already covered by run_k1_mmnl_blocked.R / tab:mmnl_blocked.)
#
# Low-intensity: single core, single-thread BLAS, nice; MMNL via Apollo (n_cores=1
# to avoid the detached-cluster deadlock). Resumable: skips (K,kappa,sigma,rep)
# cells already in output/discrete_mmnl_blocked.csv.
# Launch from the repository root:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   nohup nice -n 19 Rscript dev/run_discrete_mmnl_blocked.R > /tmp/discrete_mmnl_blocked.log 2>&1 &
# dev/rerun_standard_apollo.sh runs it with nice -n 15 and Apollo's memorySaver,
# logging to output/discrete_mmnl_blocked_std.log.
# Writes: output/discrete_mmnl_blocked.csv (one row per cell, saved after each).
# Seeds: data seed
#   as.integer(2e4 + 1000*K + 100*(kappa*100) + 10*(sigma*100) + rep);
# the shared design klue_design(n_cards = 48, n_blocks = 4) under
# set.seed(20240601); klue's clustering starts use seed 123; Apollo makes the
# MLHS draws under its default apollo_control$seed (13).
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "idefix", "mclust", "cluster"))  # stops on drift
OUT <- "output"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
CSV <- file.path(OUT, "discrete_mmnl_blocked.csv")
options(klue.cores = 1L)   # LCMNL side single-core (the old N_CORES_LCMNL namespace poke was a silent no-op)
stamp <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }

design <- klue_design(n_cards = 48L, n_blocks = 4L)
conds <- expand.grid(true_K = c(2, 3, 4), kappa = c(0.75, 1.0, 1.25),
                     sigma = c(0.15, 0.25), rep = 1:2)   # 36 discrete cells
done <- if (file.exists(CSV)) { d0 <- read.csv(CSV); paste(d0$true_K,d0$kappa,d0$sigma,d0$rep) } else character(0)
rows <- if (file.exists(CSV)) read.csv(CSV) else NULL
stamp(sprintf("Gap B: discrete LCMNL-vs-MMNL (blocked), %d cells, 1 core", nrow(conds)))

for (i in seq_len(nrow(conds))) {
  tK <- conds$true_K[i]; kap <- conds$kappa[i]; sig <- conds$sigma[i]; rp <- conds$rep[i]
  key <- paste(tK, kap, sig, rp)
  if (key %in% done) { stamp(sprintf("SKIP K=%d kap=%.2f sig=%.2f rep=%d", tK,kap,sig,rp)); next }
  seed <- as.integer(2e4 + 1000*tK + 100*(kap*100) + 10*(sig*100) + rp)
  stamp(sprintf("START K=%d kap=%.2f sig=%.2f rep=%d", tK, kap, sig, rp))
  data <- klue_simulate(N_per_class = 150, true_K = tK, separation = kap,
                        heterogeneity = sig, seed = seed, design = design)
  bic <- rep(Inf, 5)
  for (C in 1:5) { m <- tryCatch(klue_lcmnl(data$database, C), error=function(e) NULL)
                   if (!is.null(m) && isTRUE(m$converged)) bic[C] <- m$BIC }
  mnl_bic <- bic[1]; lc_C <- (2:5)[which.min(bic[2:5])]; lc_bic <- min(bic[2:5])
  mm <- tryCatch(klue_mmnl(data$database, n_cores = 1L), error=function(e) NULL)
  mm_bic <- if (!is.null(mm) && isTRUE(mm$converged)) mm$BIC else Inf
  bics <- c(MNL = mnl_bic, LCMNL = lc_bic, MMNL = mm_bic)
  winner <- names(which.min(bics))
  row <- data.frame(true_K=tK, kappa=kap, sigma=sig, rep=rp, mnl_bic=mnl_bic,
                    lc_bic=lc_bic, lc_C=lc_C, mmnl_bic=mm_bic, winner=winner,
                    lcmnl_preferred = (lc_bic < mnl_bic && lc_bic < mm_bic), stringsAsFactors=FALSE)
  rows <- rbind(rows, row); write.csv(rows, CSV, row.names = FALSE)
  stamp(sprintf("DONE  K=%d kap=%.2f sig=%.2f rep=%d -> winner=%s (LCMNL[C%d]=%.0f MMNL=%.0f MNL=%.0f)",
                tK,kap,sig,rp, winner, lc_C, lc_bic, mm_bic, mnl_bic))
}
stamp("ALL DONE")
rows$winner <- sub("\\..*","",rows$winner)
cat("\n=== Discrete DGPs, blocked design: model preference ===\n"); print(table(rows$winner))
cat(sprintf("LCMNL preferred in %d/%d discrete cells\n", sum(rows$lcmnl_preferred), nrow(rows)))
