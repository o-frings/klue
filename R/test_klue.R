# =============================================================================
# Validation: run the packaged workflow klue() on all 5 empirical datasets and
# confirm C=1..2 LLs match the reference numbers produced by the original
# empirical_*.R adapters. Tolerance: 0.01 for C=1 (deterministic MNL), 0.20 for
# C=2 (LCMNL multistart - small numerical jitter possible across runs).
# Run from the repo root:  Rscript R/test_klue.R   (a few minutes, single core)
# =============================================================================

if (requireNamespace("pkgload", quietly = TRUE) && dir.exists("klue")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else {
  suppressWarnings(suppressMessages(library(klue)))
}

PASS_LL_TOL_C1 <- 0.01
PASS_LL_TOL_C2 <- 0.20

check <- function(label, summary_df, ref_c1, ref_c2) {
  got_c1 <- summary_df$LL[summary_df$C == 1]
  got_c2 <- summary_df$LL[summary_df$C == 2]
  ok_c1 <- abs(got_c1 - ref_c1) < PASS_LL_TOL_C1
  ok_c2 <- abs(got_c2 - ref_c2) < PASS_LL_TOL_C2
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  C=1: got %.4f, ref %.4f, diff %.6f  %s\n",
              got_c1, ref_c1, got_c1 - ref_c1, if (ok_c1) "PASS" else "FAIL"))
  cat(sprintf("  C=2: got %.4f, ref %.4f, diff %.6f  %s\n",
              got_c2, ref_c2, got_c2 - ref_c2, if (ok_c2) "PASS" else "FAIL"))
  ok_c1 && ok_c2
}

all_pass <- TRUE

# -----------------------------------------------------------------------------
# 1. VITTEL (long format, CSV, indicator choice)
# -----------------------------------------------------------------------------

vittel_attrs <- c("Ban_pesticides", "RenaturationYes",
                  "Local_Negative_Effect_Ferti", "Mixedfertilizers",
                  "Hedgesconservation", "Hedgesplantation",
                  "Forest_Management_For_Water", "Forest_Managment_For_Biodiv")
res_vittel <- klue(
  data           = "test_data_files/data.csv",
  format         = "long",
  id_col         = "ID", task_col = "CS", alt_col = "alt",
  choice_col     = "choice",
  attribute_cols = vittel_attrs,
  price_col      = "Waterbill",
  choice_format  = "indicator",
  price_scaling  = 10,
  C_cands        = 1:2, run_mmnl = FALSE,
  output_prefix  = "validation_vittel", verbose = TRUE
)
all_pass <- check("Vittel", res_vittel$summary, -7234.92460881401, -6328.56968284772) && all_pass

# -----------------------------------------------------------------------------
# 2. MODE CHOICE (wide format, Apollo data, availability)
# -----------------------------------------------------------------------------

# Mode choice uses the NEW canonical names (attribute_cols, price_col, avail_col).
data(apollo_modeChoiceData, package = "apollo")
d_mode <- apollo_modeChoiceData[apollo_modeChoiceData$SP == 1, ]
res_mode <- klue(
  data           = d_mode,
  format         = "wide",
  id_col         = "ID", task_col = "SP_task", choice_col = "choice",
  attribute_cols = list(
    time    = c("time_car", "time_bus", "time_air", "time_rail"),
    access  = c(NA,         "access_bus", "access_air", "access_rail"),
    service = c(NA, NA,                  "service_air", "service_rail")
  ),
  price_col      = c("cost_car", "cost_bus", "cost_air", "cost_rail"),
  avail_col      = c("av_car",   "av_bus",   "av_air",   "av_rail"),
  scalings       = list(time = 60, access = 60, price = 10),
  C_cands        = 1:2, run_mmnl = FALSE,
  output_prefix  = "validation_mode", verbose = TRUE
)
all_pass <- check("Mode choice", res_mode$summary, -3096.56887135365, -2966.02262237717) && all_pass

# -----------------------------------------------------------------------------
# 3. SWISS ROUTE (wide format, Apollo data, no availability)
# -----------------------------------------------------------------------------

# Swiss route uses the NEW canonical names.
data(apollo_swissRouteChoiceData, package = "apollo")
res_swiss <- klue(
  data           = apollo_swissRouteChoiceData,
  format         = "wide",
  id_col         = "ID", task_col = NULL,           # synthesised
  choice_col     = "choice",
  attribute_cols = list(
    tt = c("tt1", "tt2"),
    hw = c("hw1", "hw2"),
    ch = c("ch1", "ch2")
  ),
  price_col      = c("tc1", "tc2"),
  scalings       = list(tt = 60, hw = 60, price = 10),
  C_cands        = 1:2, run_mmnl = FALSE,
  output_prefix  = "validation_swiss", verbose = TRUE
)
all_pass <- check("Swiss route", res_swiss$summary, -1665.61994629559, -1551.58940227507) && all_pass

# -----------------------------------------------------------------------------
# 4. ELECTRICITY (wide format, mlogit data, balanced-panel filter)
# -----------------------------------------------------------------------------

# Electricity uses the NEW canonical names.
data(Electricity, package = "mlogit")
res_elec <- klue(
  data           = Electricity,
  format         = "wide",
  id_col         = "id", task_col = NULL,
  choice_col     = "choice",
  attribute_cols = list(
    cl   = c("cl1",  "cl2",  "cl3",  "cl4"),
    loc  = c("loc1", "loc2", "loc3", "loc4"),
    wk   = c("wk1",  "wk2",  "wk3",  "wk4"),
    tod  = c("tod1", "tod2", "tod3", "tod4"),
    seas = c("seas1","seas2","seas3","seas4")
  ),
  price_col      = c("pf1", "pf2", "pf3", "pf4"),
  scalings       = list(price = 10),
  C_cands        = 1:2, run_mmnl = FALSE,
  output_prefix  = "validation_electricity", verbose = TRUE
)
all_pass <- check("Electricity", res_elec$summary, -4799.51007020693, -4361.72429162753) && all_pass

# -----------------------------------------------------------------------------
# 5. SWISSMETRO (wide format, .dat file, availability)
# -----------------------------------------------------------------------------

# Swissmetro uses the OLD aliases (attributes/price/availability) to verify
# backward compatibility. Should produce identical numbers to the canonical
# names.
d_sm <- read.table("test_data_files/swissmetro/swissmetro.dat", header = TRUE,
                   stringsAsFactors = FALSE)
d_sm <- d_sm[d_sm$CHOICE != 0, ]
res_sm <- klue(
  data         = d_sm,
  format       = "wide",
  id_col       = "ID", task_col = NULL,
  choice_col   = "CHOICE",
  attributes   = list(                                    # OLD name (alias)
    time    = c("TRAIN_TT", "SM_TT", "CAR_TT"),
    headway = c(NA,         "SM_HE", NA)
  ),
  price        = c("TRAIN_CO", "SM_CO", "CAR_CO"),        # OLD name (alias)
  availability = c("TRAIN_AV", "SM_AV", "CAR_AV"),        # OLD name (alias)
  scalings     = list(time = 100, headway = 100, price = 100),
  C_cands      = 1:2, run_mmnl = FALSE,
  output_prefix = "validation_swissmetro", verbose = TRUE
)
all_pass <- check("Swissmetro", res_sm$summary, -7551.62319970711, -6307.88687513661) && all_pass

# -----------------------------------------------------------------------------
# Verdict
# -----------------------------------------------------------------------------

cat("\n=========================================================\n")
if (all_pass) {
  cat("ALL 5 ADAPTERS VALIDATE: LLs match reference within tolerance.\n")
} else {
  cat("ONE OR MORE ADAPTERS FAILED VALIDATION.\n")
  quit(status = 1)
}
cat("=========================================================\n")
