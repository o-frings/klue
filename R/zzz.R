# Silence R CMD check NOTEs for variable names that Apollo expects to find
# in the caller's environment (apollo_* objects), plus a few helpers used
# inside generated character-string utility code that the static analyser
# cannot see.

utils::globalVariables(c(
  "apollo_beta", "apollo_fixed", "apollo_control", "apollo_inputs",
  "apollo_draws", "apollo_randCoeff", "apollo_probabilities",
  "apollo_lcPars",
  "DGP_DEFAULT", "N_DRAWS_MMNL"
))
