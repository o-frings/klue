# Sanity test: load the installed package and run a minimal workflow call.

library(klue)

# Use the smallest dataset (swiss route, 2 alts, 348 respondents x 9 tasks)
data(apollo_swissRouteChoiceData, package = "apollo")

res <- klue(
  data       = apollo_swissRouteChoiceData,
  format     = "wide",
  id_col     = "ID", task_col = NULL, choice_col = "choice",
  attributes = list(
    tt = c("tt1", "tt2"),
    hw = c("hw1", "hw2"),
    ch = c("ch1", "ch2")
  ),
  price      = c("tc1", "tc2"),
  scalings   = list(tt = 60, hw = 60, price = 10),
  C_cands    = 1:2, run_mmnl = FALSE,
  output_prefix = "pkgtest", write_csv = FALSE, verbose = TRUE
)

ref_c1 <- -1665.61994629559
ref_c2 <- -1551.58940227507
ok_c1 <- abs(res$summary$LL[res$summary$C == 1] - ref_c1) < 0.01
ok_c2 <- abs(res$summary$LL[res$summary$C == 2] - ref_c2) < 0.20

if (ok_c1 && ok_c2) {
  cat("\n[PASS] Installed klue package reproduces reference numbers.\n")
} else {
  cat("\n[FAIL] Installed package produced unexpected numbers.\n")
  quit(status = 1)
}
