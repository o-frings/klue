# klue rebuild plan — fewer lines, faster core

Date: 2026-06-12. Target: rebuild klue (currently ~3,900 R lines:
`engine.R` 3,064 + `workflow.R` 797) at roughly **1,300–1,500 R lines plus one
~150-line C++ kernel**, with identical exported behaviour and a 3–10× faster
LCMNL path. The package stays a normal R package.

## Where the lines currently go

| Block | Lines (approx.) | Diagnosis |
|---|---|---|
| `klue_study_*` replication drivers (Sections 6–9) | ~1,100 | 12 functions, all the same shape: `expand.grid` → seed → simulate → fit loop → `mclapply` → aggregate. Paper-replication code living inside the user package. |
| MMNL / Apollo plumbing | ~700 | `klue_mmnl` and `klue_mmnl_corr` (+ their `.run_apollo_*` and randCoeff builders) are ~90% identical: same sink/log machinery, same two-stage warm start, same fail handling. |
| Data simulation | ~440 | `klue_simulate`, `klue_simulate_cov`, `klue_simulate_deff` share ~80% of their bodies (segment deviations, individual betas, design-vs-random attribute path, choice simulation). |
| Database builders (`workflow.R`) | ~430 | Long and wide builders duplicate availability/balanced-panel/scaling/label logic. |
| Clustering starts | ~290 | Six `get_*_starts` functions differing only in one clustering call, plus a near-duplicate one-hot variant and a duplicate multistart wrapper. |
| LCMNL estimators (BFGS + EM) | ~360 | Both estimators re-implement the same MNL log-likelihood/gradient kernel three times (C=1 branch, C≥2 branch, EM weighted M-step). |
| `klue()` driver, metrics, cleanup, aliases, config | ~580 | 23 backward-compat aliases; `cleanup_apollo` maintains a 30-name "protected" list by hand; hand-written NAMESPACE and .Rd files. |

## Phase 0 — Lock current behaviour (do this before touching code)

Build a golden-output regression harness (`dev/klue_golden.R`):

- `klue_demo()` and `klue_demo(full = TRUE)`: record `summary`, `best_C`,
  `class_betas`, MMNL `LL/BIC`.
- The five empirical validation runs already scripted in `R/test_klue.R`.
- Three `klue_simulate` seeds (K = 1, 2, 3) → fit → record LL/BIC/betas/ARI.
- Save to `dev/golden/klue_0.6.3.rds`. Every later phase must reproduce these
  within 1e-6 (estimation) and bit-exactly (data generation), seeded per
  README_REPRODUCE conventions.

## Phase 1 — Move the replication harness out of the user-facing core

The 12 `klue_study_*` drivers are paper code, not workflow code. Parent-repo
scripts (`R/simulation_study.R`, `R/run_init_ablation.R`, `dev/run_*_blocked.R`,
`dev/run_revision_batch.R`) call them, so they cannot simply vanish. Plan:

1. Write one generic driver in `R/study.R` (~150 lines):
   `klue_study_run(conditions, simulate_args_fn, fit_fn, summarise_fn, n_cores)`
   handling the shared grid/seed/mclapply/aggregate/verbose plumbing.
2. Re-express each study as a ~15–30 line config on top of it
   (~250 lines total for all 12). Keep the exported `klue_study_*` names so the
   paper scripts keep working; numerical output must match the golden runs
   (the per-condition seed formulas stay byte-identical).
3. Net: ~1,100 → ~400 lines.

Backward-compat aliases: generate the 23 aliases in a 5-line loop with a
`.Deprecated()` warning, schedule removal for 1.0. Sweep parent-repo scripts to
canonical names (they are listed in NAMESPACE, easy grep).

## Phase 2 — Deduplicate the core

1. **One simulator.** `klue_simulate(..., covariates = FALSE, design = NULL)`
   absorbs `klue_simulate_cov`. Delete `klue_simulate_deff`: it is the fake
   D-efficient generator (level grid + jitter, no optimisation) superseded by
   `klue_design()` + `design=`; keep a deprecated alias that warns and calls
   the design path. CAUTION: merging must preserve each function's exact RNG
   call order under its old entry point, or the simulated datasets (and every
   paper number) change — verify against Phase 0 goldens per seed.
   ~440 → ~180 lines.
2. **One starts function.** `klue_starts(database, C, method = c("kmeans",
   "gmm","hc_ward","hc_complete","hc_average","pam"), features = c("rp",
   "onehot"))` with a small named-list registry mapping method → label vector.
   Replaces 8 functions + the duplicated one-hot multistart wrapper.
   ~400 → ~120 lines.
3. **One MNL kernel.** A single internal `mnl_ll_grad(par, X, ch_ind, weights)`
   used by: pooled MNL (C=1), each class inside the C≥2 likelihood, the EM
   weighted M-step, and `fit_cluster_mnls`. ~360 → ~200 lines, and it becomes
   the natural seam for the C++ kernel in Phase 3.
4. **One MMNL entry.** `klue_mmnl(database, correlation = FALSE, ...)`; the
   correlated variant differs only in the randCoeff builder and the start
   vector. Shared: sink/log capture, `fail_with`, bounds, two-stage warm
   start, `.run_apollo_mmnl`. Keep `klue_mmnl_corr` as a one-line wrapper.
   ~700 → ~300 lines.
5. **Apollo state hygiene.** Replace `cleanup_apollo`'s hand-maintained
   protected list: record exactly the names the package assigns into
   `.GlobalEnv` (a small character vector constant) and remove only those.
   Drops ~50 fragile lines and a maintenance trap.
6. **Database builders.** Factor shared availability-filter, balanced-panel,
   scaling, and label/attr stamping helpers out of long/wide.
   ~430 → ~300 lines.

## Phase 3 — C++ kernel for the hot path (RcppArmadillo)

**Recommendation: Rcpp/RcppArmadillo, not Julia.** The hot path is the LCMNL
log-likelihood + analytic gradient, called O(100) optim iterations × 6 starts ×
|C_cands| × (for studies) hundreds of conditions. The R version is already
vectorised; the remaining cost is repeated large-matrix allocation per optim
call. One `src/lcmnl.cpp` (~150 lines) exposing:

- `lcmnl_nll(par, X_cube, choice, T_per_n, C)` → negative LL
- `lcmnl_grad(...)` → gradient (same layout as the current R `grad_ll`)
- a row-weighted variant for the EM M-step

Both estimators and `fit_cluster_mnls` then call the same compiled kernel.
Expected 3–10× on LCMNL estimation; the full `klue_study_main` run drops
proportionally. MMNL stays in Apollo (already compiled internally).

Why not Julia: it cannot ship inside an R package — JuliaConnectoR/JuliaCall
require users to install a separate ~1 GB runtime, break `R CMD check`
portability and r-universe binary builds, and complicate cross-language seed
reproducibility. RcppArmadillo is compiled by r-universe/CRAN infrastructure;
end users installing binaries need no toolchain. The `Julia/` folder in the
parent repo can stay as research code.

Guard-rails:

- Keep the pure-R kernel as the reference implementation behind
  `engine = c("cpp", "r")` (option `klue.engine`), with a testthat test
  asserting R and C++ agree to 1e-8 on fixed inputs.
- The README's "reproduce bit-exactly" validation claim becomes
  "agree to numerical tolerance" unless the paper runs are pinned to 0.6.3
  (see Risks).

## Phase 4 — Documentation via roxygen2

NAMESPACE and the 7 .Rd files are hand-written ("Generated by hand. Keep in
sync…"). Convert exported functions to roxygen blocks and generate
NAMESPACE/man with `devtools::document()`. One source of truth, the hand-sync
burden disappears, and the docs get rebuilt automatically when the API
changes. DESCRIPTION already declares `RoxygenNote: 7.3.0`.

## Phase 5 — Verification

1. Golden regression harness from Phase 0 (both engines).
2. Existing testthat suite (`test-em.R`, `test-starts.R`, `test-corr-dgp.R`)
   plus new engine-equivalence and starts-registry tests.
3. `R CMD check --as-cran` clean.
4. Timing benchmark before/after: `klue_demo(full = TRUE)` and one
   `klue_study_main` slice (e.g. 8 conditions); record in the plan's PR.

## Risks and decision points

1. **Paper reproducibility (the big one).** v10 paper numbers were produced
   with klue 0.6.x. Any RNG-order change in the merged simulator, or
   float-ordering change from C++, perturbs results. Two safe options:
   (a) pin paper reproduction to the archived `klue_0.6.3.tar.gz` (already in
   the repo) and release the rewrite as 1.0; or (b) make the rewrite
   reproduce 0.6.3 outputs exactly with `engine = "r"` and treat C++ as an
   opt-in accelerator. Recommend (a) for the paper, (b)'s test discipline for
   the package.
2. **Parent-repo call sites** use deprecated alias names; one grep-and-rename
   sweep is part of Phase 1.
3. **Compiler dependency**: source installs now need a toolchain (Rtools/Xcode
   CLT). r-universe serves binaries, so most users are unaffected.

## Target line counts

| Component | Now | After |
|---|---|---|
| Study drivers | ~1,100 | ~400 |
| MMNL/Apollo | ~700 | ~300 |
| Simulators + design | ~440 | ~180 |
| Database builders + `klue()` driver | ~710 | ~550 |
| Starts | ~400 | ~120 |
| LCMNL estimators | ~360 | ~200 (R ref) |
| Misc (metrics, aliases, cleanup, config) | ~250 | ~120 |
| **R total** | **~3,960** | **~1,870 → ~1,500 after engine="r" path slims** |
| C++ | 0 | ~150 |

Sequencing: 0 → 1 → 2 are pure-R and each independently shippable (golden
tests green after every phase). 3 and 4 are independent of each other and can
follow in either order.
