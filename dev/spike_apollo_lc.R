# dev/spike_apollo_lc.R
# Spike: drive Apollo's latent-class EM (apollo_lcEM), written as Apollo's
# EM_LC_no_covariates.r example (explicit-logit class allocation, the class loop
# for(s in 1:length(pi_values)); dev/compare_packages.R, run_apollo_lcEM), on one
# simulated dataset and check its log-likelihood against klue's direct ML.
# Single core; tiny data. (Until 2026-09-29 this spike used the unrolled
# apollo_estimate spec, which apollo_lcEM rejects; the failure was the spec, not
# apollo_lcEM.)

suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))
sys.source("dev/compare_packages.R", envir = globalenv())   # Apollo LC builders; autorun guarded

dgp <- klue_dgp(n_generic = 4, n_alternatives = 3)
C   <- 3
d   <- klue_simulate(N_per_class = 60, T_tasks = 8, true_K = C,
                     separation = 1.0, heterogeneity = 0.2, seed = 42, dgp = dgp)
db  <- d$database

ref <- run_klue(db, C, dgp, "ml")
em  <- run_apollo_lcEM(db, C, dgp)
cat(sprintf("[klue ML]      LL = %.4f\n[apollo_lcEM]  LL = %.4f  converged = %s  (%.1f s)\n",
            ref$LL, em$LL, em$converged, em$seconds))
cat(sprintf("gap (apollo_lcEM - klue) = %+.4f\n", em$LL - ref$LL))
