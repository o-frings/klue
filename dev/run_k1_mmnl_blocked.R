# =============================================================================
# B: K*=1 (continuous) MMNL confirmation under the BLOCKED baseline.
# Goal: show that although BIC over-selects WITHIN the LCMNL family on continuous
# data under the blocked design (31% pick C>1), the workflow's triad catches it --
# MNL or MMNL attains a lower BIC than the over-selected LCMNL, so no spurious
# discrete segments are reported.
#
# LOW-INTENSITY / COOL run (per request): single-threaded BLAS + Apollo nCores=1 +
# LCMNL single-core + lowest priority + a cooling pause between conditions. Only
# one core is ever busy; the run is deliberately slow. Needs NO internet.
# Resumable: skips (sigma,rep) cells already in output/k1_mmnl_blocked_full.csv.
# Writes: output/k1_mmnl_blocked_full.csv (one row per cell, saved after each).
# Seeds: data seed as.integer(1e4 + 100*(sigma*100) + rep); the shared design
# klue_design(n_cards = 48, n_blocks = 4) under set.seed(20240601); klue's
# clustering starts use seed 123; Apollo makes the MLHS draws under its default
# apollo_control$seed (13).
#
# Launch from the repository root (cool, detached):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#   KLUE_SLEEP=15 nohup nice -n 19 Rscript dev/run_k1_mmnl_blocked.R \
#     > /tmp/k1_mmnl_blocked.log 2>&1 &
# dev/rerun_standard_apollo.sh runs it with nice -n 15, the default 120 s pause
# and Apollo's memorySaver, logging to output/k1_mmnl_blocked_std.log.
# =============================================================================
suppressMessages(suppressWarnings(pkgload::load_all("klue", quiet = TRUE)))
source("dev/pinned_versions.R"); check_pinned_versions(c("apollo", "idefix", "mclust", "cluster"))  # stops on drift
OUT <- "output"; if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
CSV <- file.path(OUT, "k1_mmnl_blocked_full.csv")  # full K*=1 grid; keeps the old 9-cell file
options(klue.cores = 1L)   # LCMNL side single-core (the old N_CORES_LCMNL namespace poke was a silent no-op)
sleep_s <- as.numeric(Sys.getenv("KLUE_SLEEP", "120"))                    # cooling pause
stamp <- function(m) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), m)); flush.console() }

# One blocked D-efficient design, shared (matches the main blocked baseline).
design <- klue_design(n_cards = 48L, n_blocks = 4L)

conds <- expand.grid(sigma = c(0.10, 0.15, 0.20, 0.25), rep = 1:20)  # all K*=1 blocked conditions (4 sigma x 20 reps = 80), matching the main blocked study
done <- if (file.exists(CSV)) { d0 <- read.csv(CSV); paste(d0$sigma, d0$rep) } else character(0)
rows <- if (file.exists(CSV)) read.csv(CSV) else NULL

stamp(sprintf("B: K*=1 MMNL confirmation (blocked), %d cells, 1 core, BLAS=1, sleep=%gs between cells",
              nrow(conds), sleep_s))

for (i in seq_len(nrow(conds))) {
  sig <- conds$sigma[i]; rp <- conds$rep[i]
  key <- paste(sig, rp)
  if (key %in% done) { stamp(sprintf("SKIP sigma=%.2f rep=%d (done)", sig, rp)); next }
  seed <- as.integer(1e4 + 100 * (sig * 100) + rp)
  stamp(sprintf("START sigma=%.2f rep=%d", sig, rp))

  data <- klue_simulate(N_per_class = 300, T_tasks = 12, true_K = 1,
                        separation = 0, heterogeneity = sig, seed = seed, design = design)

  # LCMNL family C=1..5 (C=1 == MNL); record MNL BIC and best discrete (C>=2) BIC.
  bic <- rep(Inf, 5)
  for (C in 1:5) { m <- tryCatch(klue_lcmnl(data$database, C), error = function(e) NULL)
                   if (!is.null(m) && isTRUE(m$converged)) bic[C] <- m$BIC }
  mnl_bic <- bic[1]
  lc_C    <- (2:5)[which.min(bic[2:5])]; lc_bic <- min(bic[2:5])

  # MMNL: independent normals + lognormal price, one Apollo estimation at 3000 MLHS, 1 core.
  mm <- tryCatch(klue_mmnl(data$database, n_cores = 1L), error = function(e) NULL)
  mm_bic <- if (!is.null(mm) && isTRUE(mm$converged)) mm$BIC else Inf

  bics   <- c(MNL = mnl_bic, LCMNL = lc_bic, MMNL = mm_bic)
  winner <- names(which.min(bics))
  row <- data.frame(sigma = sig, rep = rp, mnl_bic = mnl_bic, lc_bic = lc_bic,
                    lc_C = lc_C, mmnl_bic = mm_bic, winner = winner,
                    lcmnl_overselects = lc_bic < mnl_bic, stringsAsFactors = FALSE)
  rows <- rbind(rows, row)
  write.csv(rows, CSV, row.names = FALSE)   # incremental save
  stamp(sprintf("DONE sigma=%.2f rep=%d -> winner=%s (MNL=%.0f LCMNL[C%d]=%.0f MMNL=%.0f)",
                sig, rp, winner, mnl_bic, lc_C, lc_bic, mm_bic))

  if (i < nrow(conds)) { stamp(sprintf("cooling %gs ...", sleep_s)); Sys.sleep(sleep_s) }
}

stamp("ALL DONE")
cat("\n=== K*=1 (continuous) model preference under blocked design ===\n")
print(table(rows$winner))
cat(sprintf("LCMNL over-selects (best C>=2 beats MNL) in %d/%d cells; ",
            sum(rows$lcmnl_overselects), nrow(rows)))
cat(sprintf("LCMNL is the OVERALL winner in %d/%d (expect ~0: triad catches over-selection)\n",
            sum(rows$winner == "LCMNL"), nrow(rows)))
