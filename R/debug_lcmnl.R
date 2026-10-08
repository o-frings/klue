# Diagnostic: run ONE 2-class LCMNL estimation through Apollo with full output,
# started from klue's k-means starting values. Standalone: loads klue and the
# Apollo LC spec builders (make_apollo_lcPars / make_apollo_probabilities_lc)
# from dev/compare_packages.R. Run from the repo root:  Rscript R/debug_lcmnl.R
suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
suppressWarnings(suppressMessages(library(apollo)))
sys.source("dev/compare_packages.R", envir = globalenv())   # builders only; its autorun is guarded

cat("=== DIAGNOSTIC: Testing C=2 LCMNL estimation ===\n\n")

# 1. Generate test data
cat("--- Generating data (K=2, kappa=1.0, sigma=0.2) ---\n")
data <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 2,
                      separation = 1.0, heterogeneity = 0.2, seed = 12345)
cat("  N =", length(unique(data$database$ID)), "\n")
cat("  True betas:\n"); print(data$true_betas)

# 2. Get clustering starts
cat("\n--- Computing k-means starts for C=2 ---\n")
starts <- klue_starts(data$database, C = 2, method = "kmeans")
cat("  Start betas:\n"); print(starts$betas)
cat("  Start shares:", starts$shares, "\n")

# 3. Set up Apollo (NO sink, NO tryCatch)
cleanup_apollo()
C <- 2
N <- length(unique(data$database$ID))

apollo_beta_vec <- c()
for (ci in 1:C) {
  apollo_beta_vec[paste0("asc_alt1_", ci)] <- 0
  apollo_beta_vec[paste0("asc_alt2_", ci)] <- 0
  apollo_beta_vec[paste0("b_x1_", ci)]     <- starts$betas[ci, 1]
  apollo_beta_vec[paste0("b_x2_", ci)]     <- starts$betas[ci, 2]
  apollo_beta_vec[paste0("b_x3_", ci)]     <- starts$betas[ci, 3]
  apollo_beta_vec[paste0("b_x4_", ci)]     <- starts$betas[ci, 4]
  apollo_beta_vec[paste0("b_price_", ci)]  <- starts$betas[ci, 5]
  apollo_beta_vec[paste0("asc_alt3_", ci)] <- 0   # reference constant, fixed (as in Apollo's EM example)
}
for (ci in 1:C) {
  apollo_beta_vec[paste0("delta_", ci)] <- if (ci < C) {
    log(max(starts$shares[ci], 0.01) / max(starts$shares[C], 0.01))
  } else 0
}
apollo_fixed_vec <- c(paste0("asc_alt3_", 1:C), paste0("delta_", C))

cat("\n--- apollo_beta ---\n"); print(apollo_beta_vec)
cat("--- apollo_fixed ---\n"); print(apollo_fixed_vec)

apollo_control <<- list(
  modelName       = "DEBUG_LC2",
  modelDescr      = "2-class LCMNL debug",
  indivID         = "ID",
  nCores          = 1,
  outputDirectory = tempdir()
)
apollo_beta  <<- apollo_beta_vec
apollo_fixed <<- apollo_fixed_vec

cat("\n--- Defining apollo_lcPars (unrolled, apollo_classAlloc) ---\n")
apollo_lcPars <<- make_apollo_lcPars(C, DGP_DEFAULT)
cat("  Environment:", environmentName(environment(apollo_lcPars)), "\n")
cat("  Generated function body:\n")
print(body(apollo_lcPars))
cat("  Done.\n")

cat("\n--- Calling apollo_validateInputs ---\n")
apollo_inputs <<- apollo_validateInputs(
  apollo_beta    = apollo_beta,
  apollo_fixed   = apollo_fixed,
  database       = data$database,
  apollo_control = apollo_control
)
cat("  Validation OK!\n")

cat("\n--- Defining apollo_probabilities (unrolled, apollo_lc) ---\n")
apollo_probabilities <<- make_apollo_probabilities_lc(C, DGP_DEFAULT)
cat("  Environment:", environmentName(environment(apollo_probabilities)), "\n")
cat("  Generated function body:\n")
print(body(apollo_probabilities))
cat("  Done.\n")

cat("\n--- Calling apollo_estimate ---\n")
model <- apollo_estimate(
  apollo_beta          = apollo_beta,
  apollo_fixed         = apollo_fixed,
  apollo_probabilities = apollo_probabilities,
  apollo_inputs        = apollo_inputs,
  estimate_settings    = list(maxIterations = 500, printLevel = 0)
)
cat("\n--- Estimation complete! ---\n")
cat("  LL:", model$LLout[1], "\n")
cat("  Estimates:\n"); print(model$estimate)

cleanup_apollo()
cat("\n=== DIAGNOSTIC DONE ===\n")
