# Silence R CMD check NOTEs: apollo_probabilities is assigned into the
# caller's global environment (Apollo requires it there), and the generated
# model functions reference variables the static analyser cannot see.

utils::globalVariables(c(
  "apollo_probabilities",
  "rp_at_global", "oh_at_global", "random_pct_global",
  "at_global_rate", "reached_global", "K"
))
