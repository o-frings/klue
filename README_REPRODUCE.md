# Reproduction materials

Code and result files for *Classes or Continuum? A Pre-Registrable
Specification Rule for Latent Class Logit Models* (Frings).

Everything in the paper is produced by the **`klue`** R package (the estimation
engine), the simulation-study drivers in **`studies/klue_studies.R`**, and the
orchestration scripts in **`dev/`**. Legacy code (the pre-package `simulation_study.R`,
the 0.6.x engine and the version-comparison harness) is kept in `archive/code_backups/` of the
working repository and is not published.

> **AI assistance.** Several analysis and orchestration scripts in `dev/` (for
> example the package stress-ladder runners and the EM-variant comparison) were
> developed with AI assistance (Claude, Anthropic). All results reported in the
> paper were reviewed and verified by the author, who takes full responsibility
> for them.

## Contents

| Path | What it is |
|------|------------|
| `klue/` | The R package (0.10.0): LCMNL/MMNL estimators, `klue()` workflow, inference helpers, `klue_simulate`/`klue_simulate_choices`/`klue_design`. Has its own `README.md`, `DESCRIPTION`, `NAMESPACE`, `man/`, `tests/`. |
| `studies/klue_studies.R` | The `klue_study_*` drivers that run every simulation experiment; `source()` it after loading the package. |
| `dev/` | Orchestration scripts that run the paper's experiments under the randomised and blocked designs, and the package head-to-head. |
| `R/empirical_*.R` | The five empirical applications (Rhin-Meuse, file names `vittel`; Apollo mode, Apollo Swiss route, Electricity, Swissmetro). Each loads `klue` via `pkgload`. |
| `output/` | Result files (`.csv`/`.rds`) reported in the paper. |

## Requirements

- R 4.3.1 and the package versions below; `renv.lock` (written by
  `dev/record_versions.R`) lists the full dependency closure, and
  `output/package_versions.csv` the packages the analyses call. `klue` imports
  `mclust` and `cluster` and suggests `apollo`, `idefix`, `mlogit` and
  `testthat`; every script loads it from this repository's `klue/` tree with
  `pkgload::load_all("klue")` (installing the package works too:
  `R CMD INSTALL klue`, or `remotes::install_github("o-frings/klue@v0.10.0",
  subdir = "klue")`). `logistf` runs the Firth regressions
  (`dev/run_firth_h3.R`) and `numDeriv` the standard-error check.
- **`klue_mmnl` (the MMNL benchmark) uses `apollo` internally** as its simulated-ML
  (MLHS) engine. The LCMNL side is pure `klue` (`optim` + analytic gradients).
- **The package head-to-head additionally runs the competitors' own routines**
  (`dev/compare_packages.R`): Apollo's `apollo_searchStart` and `gmnl`'s latent-class
  estimator. `gmnl` 1.1-3.2 requires a pre-`dfidx` `mlogit` (< 1.1). The harness
  looks for it in the side library `~/R/oldmlogit_lib` (`KLUE_OLD_MLOGIT_LIB`
  points elsewhere) and puts that library first on the path, so only `gmnl`
  sees the old `mlogit`; the `gmnl` arm stops if it is missing.

### Package versions behind the results

| Package | Version | Used for |
|---|---|---|
| R | 4.3.1 | everything |
| `klue` | 0.10.0 | this repository's `klue/` tree |
| `apollo` | 0.3.5 | `klue_mmnl()`, every Apollo arm, the mode-choice and Swiss route data |
| `gmnl` | 1.1-3.2 | the `gmnl_lc` arm of both stress ladders |
| `mlogit` | 1.0-3.1, side library | `gmnl`'s data preparation |
| `mlogit` | 1.1-1 | the Electricity data |
| `idefix` | 1.1.0 | `klue_design()`, the blocked D-efficient design |
| `mclust` | 6.1.1 | `klue`'s Gaussian-mixture clustering start |
| `cluster` | 2.1.4 (ships with R 4.3.1) | `klue`'s PAM clustering start |
| `logistf` | 1.26.1 | the Firth regressions |

No package the analyses load from the main library changed while the results
were produced (the installation dates are in `output/package_versions.csv`;
R 4.3.1 was the only R installed). The side-library copy of `mlogit` 1.0-3.1
was reinstalled from CRAN on 2026-10-06, the date in that file (the earlier
copy lived in `/tmp`); `dev/check_ladder_carryover.R` then refitted every
`gmnl_lc` cell with it and reproduced the stored log-likelihoods exactly. CRAN now serves newer releases of `apollo` (0.3.9), `gmnl` (1.1-4),
`mlogit` (2.0-0) and `mclust` (6.1.3), and `apollo_searchStart` has used BGW
since Apollo 0.3.8, so `install.packages()` does not reproduce the tables.
Install the pinned versions,

```r
remotes::install_version("apollo",  "0.3.5")
remotes::install_version("gmnl",    "1.1-3.2")
remotes::install_version("idefix",  "1.1.0")
remotes::install_version("mclust",  "6.1.1")
remotes::install_version("logistf", "1.26.1")
dir.create("~/R/oldmlogit_lib", recursive = TRUE)  # install.packages() stops if lib does not exist
remotes::install_version("mlogit",  "1.0-3.1", lib = "~/R/oldmlogit_lib",
                         dependencies = FALSE)   # side library, for gmnl only
```

or restore the main library with `renv::restore(lockfile = "renv.lock")`; the
lockfile holds one `mlogit` (1.1-1, for the Electricity data), so the
side-library lines above are needed either way.
`dev/pinned_versions.R` holds the pins. The comparison harness
(`dev/compare_packages.R`, hence both stress ladders, the isolation arms, the
benchmark and the reference-optimum check) and the MMNL runners stop when an
installed version differs; `KLUE_ALLOW_VERSION_DRIFT=1` runs anyway, for
exploratory work only.

**Seeds.** Every script that draws random numbers sets its seeds and states
them in its header. `klue` seeds its own stochastic steps (the k-means, Gaussian
mixture and PAM clusterings with seed 123, `klue_design()` with 20240601,
`klue_simulate()` with its `seed` argument) and restores the caller's random
stream, so the `R/empirical_*.R` adapters need no `set.seed()`. Apollo makes
its draws under `apollo_control$seed` (default 13), and `apollo_searchStart`
draws its candidates under that seed plus 4. Wall-clock times will differ from
machine to machine; log-likelihoods and estimates reproduce.

## Data

No restricted-access data are redistributed here. Obtain the datasets from their
sources:

- Rhin-Meuse water-quality DCE: the survey data are deposited at `doi:10.57745/MXDUB2`
  (Amiri, Abildtrup, Garcia and Montagné-Huck 2024) and analysed in Amiri, Abildtrup and
  Garcia (2024), *Revue d'économie politique* 134(5): 739–770. Scripts, arguments and
  result files keep the case's earlier name, `vittel`/`Vittel`. The analysis file
  `test_data_files/data.csv` (loaded by `R/empirical_application.R`; long format, 930
  respondents × 8 tasks × 3 alternatives = 22,320 rows) is rebuilt from three inputs by
  `Rscript R/prepare_rhinmeuse_data.R --check`, which also compares the result with an
  existing `data.csv`:
  - the deposit's CSV, placed in `test_data_files/MXDUB2/`;
  - `test_data_files/rhinmeuse_design.csv`, the choice-card design (3 versions × 8 cards
    × 3 alternatives) from the survey's authors;
  - `test_data_files/rhinmeuse_kept.csv`, the 930 respondents kept by the authors'
    cleaning (the deposit's Stata files) and the design version each saw.
  TODO: the two small tables go into this repository with the authors' permission, or
  they are available from them on request.
- Apollo mode-choice and Swiss route-choice: the `apollo` R package.
- Electricity: the `mlogit` R package.
- Swissmetro: distributed with BIOGEME; place it at
  `test_data_files/swissmetro/swissmetro.dat` (`R/empirical_swissmetro.R`).

The Apollo and `mlogit` datasets load from the installed packages.

## Manuscript table → script map

**Simulation, randomised (idealised) design** — the paper's "randomised" regime;
some file and function names retain the historical `orthogonal` label.

| Table | Produced by | Result file |
|-------|-------------|-------------|
| Convergence (H1), `tab:convergence` | `klue_study_convergence()` | `output/convergence_results.csv` |
| Estimator × start (H1b), `tab:est_start` | `klue_study_estimator()` | `output/estimator_orthogonal.rds` |
| Feature ablation (H1a), RP contrasts vs one-hot | `R/run_init_ablation.R` (`klue_study_initialisation()`) | `output/init_ablation.rds`, `output/init_ablation_summary.txt` |
| Package stress ladder, `tab:stress` | `dev/stress_replicate.R` | `output/stress_replicate.rds` |
| IC performance / K×κ, `tab:overall`, `tab:K_kappa` | `klue_study_main()` | `output/main_results.csv` |
| LCMNL vs MMNL (H4), `tab:mmnl` | `klue_study_mmnl()` | `output/mmnl_results.csv` |
| Recovery / ARI, `tab:recovery`, `tab:ari` | `klue_study_main()` | `output/main_results.csv` |
| Robustness, `tab:unbalanced`/`design`/`concomitant`/`sample_sens` | `klue_study_unbalanced/design/concomitant/sample()` | `output/*_results.csv` |
| Correlated MMNL, `tab:mmnl_corr` | `klue_study_mmnl_corr()` | `output/mmnl_correlated_results.csv` |

**Simulation, blocked D-efficient (realistic) design** — `blocked = TRUE`

| Table | Produced by | Result file |
|-------|-------------|-------------|
| Enumeration K×κ, `tab:K_kappa_blocked` | `dev/run_blocked_baseline.R` (`klue_study_main(blocked=TRUE, n_reps=20)`) | `output/main_blocked.csv` |
| Convergence/init (H1/H1a) | `dev/run_blocked_baseline.R` | `output/{convergence,init_ablation}_blocked.csv` |
| Estimator × start (H1b), `tab:est_start_blocked` | `dev/run_estimator_blocked.R` (`klue_study_estimator(blocked=TRUE)`) | `output/estimator_blocked.csv` |
| Package stress ladder, `tab:stress_blocked` | `dev/stress_replicate_blocked.R` (klue ML/EM, `lclogit`-style EM+random-partition, Apollo `searchStart` at defaults, `apollo_lcEM`, Apollo from the clustering starts, gmnl; Apollo's reduced search, the non-default arm version12 reported, is kept as `apollo_ss_published` and runs through the same runner; `dev/rerun_apollo_blocked.R`, which first wrote it, is superseded and stops) | `output/stress_replicate_blocked.rds` |
| H4 MMNL, `tab:mmnl_blocked` | `dev/run_k1_mmnl_blocked.R` (K*=1), `dev/run_discrete_mmnl_blocked.R` (K*>=2) | `output/k1_mmnl_blocked_full.csv`, `output/discrete_mmnl_blocked.csv` |
| Robustness, `tab:robust_blocked` | `dev/run_robustness_blocked.R` | `output/{unbalanced,concomitant,sample}_blocked.csv` |

**Empirical applications** (Appendix): `R/empirical_*.R` → `output/*_lcmnl_results.csv`,
`output/*_model_comparison.csv` (Vittel: `R/empirical_application.R` →
`output/empirical_results.csv`, `output/empirical_class_betas.csv`,
`output/empirical_model_comparison.csv`, in the pre-2026-09-21 coding; the
paper's Vittel numbers come from `dev/rerun_vittel_respec.R`, below); the MMNL benchmark (independent and correlated,
fixed and random constants) via `dev/rerun_mmnl_benchmark.R` (see "Standard-Apollo
reruns" below).
Regime classification of the five applications: `R/run_regime_classification.R` →
`output/regime_classification.{rds,csv}` — written but **not yet run to
completion** (the May run stopped inside the correlated MMNL; only the log
survives), so no result file is shipped.

**Vittel, identified specification.** The published Vittel numbers entered both
levels of the two-level forest attribute alongside the status-quo ASC, so
`Forest_For_Water + Forest_For_Biodiv + ASC_sq` was constant across
alternatives and the model was rank-deficient by one parameter per class: the
likelihood is unaffected but `k` counted a parameter the design cannot
identify, over-penalising BIC/AIC/ICL by `C·log(N)`, and the three affected
coefficients were identified only up to a constant.
`R/empirical_application.R` now codes the forest attribute against a reference
level (`forest_ref = "biodiv"`; pass `"none"` to reproduce the earlier coding,
which attains the same log-likelihood). `dev/rerun_vittel_respec.R` re-runs the
whole Vittel arm — LCMNL ladder extended until BIC turns, independent MMNL and
correlated MMNL — under the identified coding, writing `output/vittel_respec_*`
alongside the published files rather than over them. klue ≥ 0.9.4 detects this
class of defect: a rank-deficient fit sets `vcov`/`robust_vcov` to `NA`, flags
`singular`, and names the aliased parameters in `aliased`.

**Confidence intervals and tests** quoted in the text (Wilson CIs, McNemar, seeded
bootstrap CIs, Wilcoxon): `R/statistical_tests.R`, which reads the `output/*.csv`
files above.

**Figures**: `dev/make_figures.R` regenerates `fig_ic.pdf`, `fig_reliability.pdf`,
`fig_time.pdf` in the repository root from `output/main_blocked.csv`,
`output/estimator_blocked.rds`, `output/stress_replicate_blocked.rds` and, for the
times of Apollo's reduced search (measured on seed 1 of each rung),
`output/apollo_smartstart_diag.csv` (deterministic, no RNG); it prints the plotted
percentages and per-rung mean times. Each figure is
saved at the exact physical width it is included at in the manuscript
(`\linewidth` = 6.925 in; `fig_time` at 0.62 `\linewidth`), so all figure text
prints at a uniform 9 pt with no LaTeX rescaling.

## Revision (version10) reruns

Scripts added for the version10 revision. Each is self-contained, seeded where it
uses RNG, and writes to `output/` with a stable filename. Light jobs run
immediately; **heavy jobs are parked** — launch them in a dedicated compute
window (they write fresh filenames and do not overwrite the v9 artifacts).

| Gate / item | Script | Run | Output | Cost |
|---|---|---|---|---|
| Firth H3 (G2, G3) | `dev/run_firth_h3.R` | `Rscript dev/run_firth_h3.R` | `output/firth_h3_{orthogonal,blocked}.csv` | light |
| Blocked recovery/ARI (G3, Path A) | `dev/aggregate_recovery.R` | `Rscript dev/aggregate_recovery.R` | `output/{recovery,ari}_blocked.csv` | light |
| Correlated MMNL ×5 (G4) | `dev/empirical_corr_mmnl.R` | `Rscript dev/empirical_corr_mmnl.R` | `output/empirical_corr_mmnl_3000draws.csv` | **parked/heavy** |
| Extend C until BIC turns (G5) | `dev/run_emp_cmax_extend.R` | `Rscript dev/run_emp_cmax_extend.R` | `output/emp_cmax_extend_<dataset>.csv` | **parked/heavy** |

The correlated-MMNL script is incremental/resumable (one dataset per row, skips
already-done datasets); the extend-C script resumes from the last fitted `C`, so
an interrupted compute window loses nothing. Both are deterministic: MLHS draws
are fixed and klue seeds its clustering starts itself (seed 123).

Gate 1 diagnostic done: `dev/run_apollo_smartstart.R` →
`output/apollo_smartstart_diag.csv` (rescaling does not rescue Apollo's
non-default `smartStart=TRUE`; the ladders now run `apollo_searchStart` at its
defaults). Minor 7
(full-grid k1 MMNL) done → `output/k1_mmnl_blocked_full.csv`. Still to run (need a
compute window and the pinned-`mlogit`+`gmnl` side-library): the proper-`idefix`
`tab:design` rerun (Minor 9). The randomised package stress ladder (5 seeds per rung)
(`dev/stress_replicate.R`, Minor 6) has since run; see "Standard-Apollo reruns".

## Standard-Apollo reruns (klue 0.10.0)

Since klue 0.10.0 (2026-09-29) every Apollo call follows Apollo's own example
scripts (apollochoicemodelling.com, examples page):

- **MMNL** (`klue_mmnl()`, as in `MMNL_preference_space.r` and
  `MMNL_preference_space_correlated.r`): random coefficients
  `mu + sigma * draws` with `sigma` unconstrained, a negative log-normal price
  `-exp(mu + sigma * draws)`, the correlated model written out as a
  lower-triangular Cholesky and started from the independent estimates
  (`apollo_readBeta` step: same means, standard deviations on the diagonal,
  zero correlations), one `apollo_estimate()` call at 3,000 MLHS draws with
  Apollo's default settings (BGW), and convergence as Apollo reports it
  (`successfulEstimation`). klue's departures from the examples: starts from
  the pooled MNL (standard deviations at the examples' 0.01), Cholesky
  off-diagonals started at 0 rather than 0.01 (so the correlated start
  reproduces the independent log-likelihood), 3,000 MLHS draws, BIC on
  respondents, and, on data with availability columns, utilities measured from
  the chosen alternative's (Apollo 0.3.5's `apollo_mnl` otherwise returns
  `Inf * 0 = NaN` for chosen utilities below about -709; the probabilities are
  unchanged). Hand-written scripts in the examples' form give the same
  estimates to machine precision.
- **Latent class, Apollo arms of the stress ladders** (`dev/compare_packages.R`):
  model code as in `LC_no_covariates.r` (`apollo_classAlloc`, literal class
  loop) for `apollo_estimate`/`apollo_searchStart`, and as in
  `EM_LC_no_covariates.r` (explicit-logit allocation,
  `for(s in 1:length(pi_values))`) for `apollo_lcEM`. `apollo_searchStart`
  runs at its 0.3.5 defaults (100 candidates within +-0.1 of the start,
  `smartStart = FALSE`, 5 stages); the start is the pooled MNL with class `c`
  scaled by `1/c` (the examples start class b at half of class a), listed
  parameter by parameter with the last class as reference, as in the
  examples. The non-default arms published in version12 are kept as
  `apollo_ss_published` (their start, parameter order and reference class
  included), so their numbers can be reproduced. Hand-written scripts in the
  examples' form give bitwise-identical results.

klue up to 0.9.5 differed in ways that change results: `mu + exp(sigma) *
draws` for the independent MMNL (it can stop at a saddle point when a standard
deviation is zero; `output/mmnl_bench/Vittel_M1.rds` shows the signature,
`sigma_x3 = -47.5` and no covariance matrix), a 200-draw first stage, automatic
retries and fallbacks, a price box passed as `estimate_settings$bounds`
(Apollo has no such setting, so it never applied), `converged` = finite
estimates rather than Apollo's own flag, and `silent = TRUE` (needed only
because an integer draw count overflowed Apollo's core-count hint). Every MMNL number in version12 came
from that code and is refitted below.

Run from the repository root:

```
bash dev/rerun_standard_apollo.sh            # one single-core lane (default)
bash dev/rerun_standard_apollo.sh --two-lanes
```

Every step runs with Apollo's `memorySaver` (`KLUE_MMNL_MEMORY_SAVER=TRUE`),
which leaves the estimates bitwise identical and roughly halves peak memory:
the Vittel correlated MMNL at 3,000 draws took about 12 GB without it, and two
lanes in parallel overflowed the 32 GB machine. It first archives the
superseded outputs to `output/superseded_klue095/`, then runs, each step seeded
and resumable (lane A then lane B in the default single lane):

| Lane | Script | Output | Replaces |
|---|---|---|---|
| A | `dev/rerun_mmnl_benchmark.R Vittel` | `output/mmnl_bench_std/Vittel_M{1..4}.rds`, `output/vittel_respec_lcmnl_ext.csv` | `output/mmnl_bench/` (kept as the 0.9.5 record) |
| A | `dev/rerun_mmnl_benchmark.R Swissmetro Mode SwissRoute Electricity` | `output/mmnl_bench_std/*.rds` | tab:emp_summary MMNL columns and the appendix MMNL rows |
| A | `dev/run_mmnl_second_start.R` | `output/mmnl_bench_std/robust/<dataset>_M3_zero.rds` and `_M1_zero.rds`, possibly replaced `<dataset>_M1.rds` to `_M4.rds` (previous fits kept in `robust/`), `output/mmnl_second_start.csv` | new (see "Second start" below) |
| A | `dev/rerun_mmnl_benchmark.R summary` | `output/mmnl_benchmark_v12.csv` | the benchmark table, refreshed after the second start |
| A | `dev/run_k1_mmnl_blocked.R` | `output/k1_mmnl_blocked_full.csv` | tab:mmnl_blocked, K* = 1 row |
| A | `dev/run_discrete_mmnl_blocked.R` | `output/discrete_mmnl_blocked.csv` | tab:mmnl_blocked, K* >= 2 row |
| A | `dev/run_h4_misspec_blocked.R` | `output/h4_misspec_blocked.csv` | new |
| B | `dev/stress_replicate_blocked.R` | `output/stress_replicate_blocked.rds` (topped up: `apollo_searchStart` at defaults, `apollo_lcEM`, `apollo_clust`) | tab:stress_blocked, fig:reliability(b), fig:time |
| B | `dev/stress_replicate.R` | `output/stress_replicate.rds` (topped up: `apollo_searchStart` at defaults, `apollo_lcEM`) | tab:stress |
| B | `dev/run_mmnl_studies.R mmnl mmnl_corr` | `output/mmnl_results.csv`, `output/mmnl_correlated_results.csv` (each with a `.done` marker) | tab:mmnl, tab:mmnl_corr |

**Swissmetro M3 check (2026-10-01).** The Swissmetro M3 benchmark, estimated
once from the pooled MNL, stopped 18.8 LL below the klue 0.9.5 fit of the same
model with the same draws. `dev/run_swissmetro_m3_robustness.R` refits it from
that optimum and from a neutral start, keeps the highest-LL fit as the
benchmark (the original stays as `output/mmnl_bench_std/Swissmetro_M3_pooled_mnl.rds`)
so that M4 starts from it, and adds the train- and Swissmetro-reference
variants as a sensitivity (`output/mmnl_bench_std/robust/`,
`output/swissmetro_m3_robustness.csv`). It runs before Swissmetro M4.

**Apollo isolation arms (2026-10-04).** On the blocked ladder, Apollo's default
search starts from a vector this project builds (klue's pooled MNL, class c
scaled by 1/c), while the reported arm used reduced settings, a [-3, 3] window
and a generic start. `dev/run_apollo_isolation_blocked.R` separates these on
ten cells fixed in advance (hard and very-hard seeds 1-5, the ladder's data):
one fit with no search from either start (A0, A4), the reduced budget inside
the default box (A1), the reported window from the MNL-based start (A3), and
Apollo's defaults from the generic start (A2). A klue check arm (`chk`)
confirms that the regenerated data match the ladder's. How each outcome is to
be read is written in the script header. Run per group of arms, e.g.
`ARMS=chk,A0,A4 Rscript dev/run_apollo_isolation_blocked.R`; each group writes
`output/apollo_isolation_blocked_<ARMS>.rds`, and
`Rscript dev/run_apollo_isolation_blocked.R summary` merges them into
`output/apollo_isolation_blocked.csv`.

**Apollo isolation arms A5 and A6 (2026-10-05).** Two more arms of
`dev/run_apollo_isolation_blocked.R` show where Apollo's time goes, on the same
ten cells. Apollo 0.3.5's `apollo_searchStart` skips the model pre-processing
that `apollo_estimate` applies, so its search runs on numerical gradients. A5
is the ladder's default search with that pre-processing
(`apollo_modifyUserDefFunc`) applied first: the same candidates, searched on
analytic gradients. A6 is the `apollo_clust` arm with one covariance matrix
instead of one per fit, computed at the best of the six fits as klue_ml does;
it times the clustering starts, the six fits and the covariance separately.
The script header fixes how A5 is to be read; A6 measures cost only.
`ARMS=A6,A5 Rscript dev/run_apollo_isolation_blocked.R` (one core, about 1 h
45 min) writes `output/apollo_isolation_blocked_A6A5.rds`; the `summary` mode
adds both arms to `output/apollo_isolation_blocked.csv`. `ISO_TEST=1` runs a
small version on one cell into a scratch directory in about 4 minutes.

**Reference optimum per ladder cell, Stage 1 (2026-10-05).** The blocked
ladder fits every method at the true class count K* only and scores it
against the best LL any method reached in the cell.
`dev/run_ladder_reference_optimum.R` computes an independent reference on the
ladder's regenerated data for the hard and very-hard cells (seeds 1-10) and
moderate seed 2: at each C = 1..6, the highest LL of klue's clustering starts
with direct ML, `run_klue_em_rp` (six random partitions + EM) and 50 random
partitions with direct ML. It first checks that klue_ml at K* reproduces the
ladder's LL in every cell. From the reference BIC (k = 8C - 1, N =
respondents) it records the selected C and the margin to the runner-up, and
for each ladder method the shortfall g at K* and whether 2g exceeds that
margin. One core, about 4-7 h, resumable per cell; writes
`output/ladder_reference_optimum.{rds,csv}`. `REF_TEST=1` fits one cell with
three random starts into `REF_TEST_DIR`.

**Well-separated benchmark at Apollo's defaults (2026-10-05).** version12
quoted one well-separated benchmark (LL -3021.2; klue 1.5 s, Apollo 173 s
with the reduced, non-default search); no log of that run survives.
`dev/run_benchmark_apollo_defaults.R` regenerates its dataset (the validation
gate in `dev/compare_packages.R`: `klue_simulate(N_per_class = 100, T_tasks =
12, true_K = 3, separation = 1.5, heterogeneity = 0.2, seed = 42)`, C = 3;
klue's LL must round to the published -3021.2) and fits it four ways: klue,
`apollo_searchStart` at Apollo 0.3.5's defaults, the same search on analytic
gradients, and the reduced configuration. One core, about 20-40 min,
resumable per arm; writes `output/benchmark_apollo_defaults.{csv,rds}`.
`BENCH_TEST=1` runs tiny searches into `BENCH_TEST_DIR`.

**MMNL studies on a second core (2026-10-03).** `dev/run_mmnl_studies.R` does
not depend on the other steps, so it may run beside the lane, on one core
(within the machine's core cap). It skips a study whose CSV was already
written by the current klue version (marker `output/<csv>.done`;
`KLUE_FORCE=1` redoes it), and a second instance exits while one is running
(lock `output/.run_mmnl_studies.lock`), so the lane's own launch of the step
costs nothing once the parallel run has started or finished.

**Second start (2026-10-01; symmetric since 2026-10-05).** A local maximum on
either side can decide a close reading: a lower MMNL maximum biases it toward
the LCMNL, a lower LCMNL maximum toward the MMNL. Wherever the best converged
MMNL is not at least 10 BIC points below the BIC-best LCMNL, both sides get a
second set of starts; the rule is the same for all five datasets.
MMNL side, `dev/run_mmnl_second_start.R`: M3 is refitted from a neutral start
(`start = "zero"`) and, if that reaches a higher log-likelihood, becomes the
benchmark M3 with M4 refitted from it; when the best MMNL has fixed constants
(M1 or M2), M1 is refitted the same way and M2 from it. Previous fits are kept
in `output/mmnl_bench_std/robust/`. `KLUE_SECOND_START_PLAN=1` lists the fits a
run would make without estimating anything; `KLUE_SECOND_START_TEST=1
KLUE_SECOND_START_SMOKE=<dir>` runs every step on a tiny Mode benchmark.
Afterwards `Rscript dev/rerun_mmnl_benchmark.R summary` refreshes
`output/mmnl_benchmark_v12.csv`.
LCMNL side, `dev/check_empirical_lcmnl_maxima.R`: for each dataset it refits
the BIC-best C* and its neighbours (C* - 1, C*, C* + 1; when the BIC-best C is
the highest C searched and below 12, also C + 1 above the stored ladder, which
added Swiss route C = 7, where BIC rises) with `klue_lcmnl()` defaults, to
check the stored log-likelihood
reproduces, and with random-partition starts (60 at C*, 20 at each
neighbour), each followed by direct ML; a gain above 0.01 is flagged and the
margin to the best MMNL recomputed. A fit whose largest class-specific
|coefficient| exceeds 500 counts as degenerate and is left out of the
BIC-best C and the margin (columns `*_any` keep it): at Swiss route C = 6,
eleven starts end within 0.4 LL of -1406.81 with a 3% class and coefficients
of 916 to 2,676, while no fit the tables use has a coefficient above 261. A
rerun after the checkpoint is complete refits nothing and only rebuilds the
CSV. It is run on all five datasets (since 2026-10-05), one core, about 35
minutes, resumable; writes
`output/empirical_lcmnl_maxima.{csv,rds}`. `EMP_MAX_TEST=1
EMP_MAX_TEST_DIR=<dir>` runs Mode at C* with two random starts.
The lane runs the empirical benchmark before the H4 steps (order changed
2026-10-01; the order does not change any number).

**Checks behind numbers in the text (2026-10-06).** These scripts compute
numbers the manuscript quotes that the runners above do not print, or check
claims about the software. Each is deterministic or seeded, runs from the
repository root on one core, and writes to `output/`.

| What it backs | Script | Run | Output | Cost |
|---|---|---|---|---|
| Ladder and isolation statistics the text quotes: miss counts at 0.1-5 LL, exact McNemar tests against clustering+ML, miss sizes and BIC equivalents, Wilson interval, per-rung times, cells where one arm alone reached the best, the `apollo_lcEM` iteration cap, and clustering+ML at `klue`'s screened default on the randomised ladder (cells that move, misses, exact McNemar test against Apollo's default search; from `output/ladder_carryover_check.csv`) | `dev/ladder_statistics.R` | `nice -n 15 Rscript dev/ladder_statistics.R` | `output/ladder_statistics.csv`, `output/apollo_lcem_cap.csv` | light (about 1 s; no estimation, no RNG) |
| `klue`'s C >= 2 standard errors: analytic gradient against numDeriv; classical SEs against an independently coded numerical Hessian and against Apollo 0.3.5 started at klue's solution; robust SEs against an independent sandwich and against Apollo (factor sqrt(N/(N-1))); SE error when price is rescaled | `dev/check_lcmnl_standard_errors.R` | `nice -n 15 Rscript dev/check_lcmnl_standard_errors.R` | `output/lcmnl_se_check.csv` | light (about 35 s) |
| Apollo 0.3.5 gradients: `apollo_searchStart` on models written as in Apollo's examples (MNL, mode-choice MNL, latent class) passes `grad = NULL` to maxLik, costs 2K + 1 evaluations per numerical gradient and prints "not able to compute analytical gradients"; analytic after `apollo_modifyUserDefFunc` and in `apollo_estimate`; `apollo_lcEM` numerical; `apollo_searchStart`'s default settings | `dev/check_apollo_searchstart_gradients.R` | `nice -n 15 Rscript dev/check_apollo_searchstart_gradients.R` | `output/apollo_searchstart_gradients.csv` | light (about 10 s) |
| Stress-ladder arms carried over from earlier code, refitted with the current code: klue_ml, klue_em and gmnl_lc on every cell, klue_em_rp on the blocked cells, Apollo's reduced search on blocked moderate seed 2 and seed 1 of every rung of both ladders; a klue arm that misses its stored LL by more than 1e-3 is refitted with `options(klue.screen = FALSE)` | `dev/check_ladder_carryover.R` | `KLUE_CORES=1 nohup nice -n 15 Rscript dev/check_ladder_carryover.R blocked > output/ladder_carryover_blocked.log 2>&1 &`, the same with `randomised`, then `Rscript dev/check_ladder_carryover.R summary` | `output/ladder_carryover_check.csv` (lane checkpoints `output/ladder_carryover_check_{blocked,randomised}.rds`) | one core per lane, about 1.5-2 h; resumable |
| EM from the clustering starts ending below direct ML on the randomised ladder's correlated rung (the ML-or-EM paragraph): EM refitted from the six starts at a stopping tolerance of 1e-6 (up to 500 iterations, the ladder's setting) and 1e-12 (up to 5000) on every hard_corr seed and very_hard seed 4, against klue_ml with `options(klue.screen = FALSE)`, which reproduces its stored LL (the data check) | `dev/check_em_tol.R` | `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 15 Rscript dev/check_em_tol.R` | `output/em_tol_check.csv` | one core, about 25 min |
| Package versions behind the results | `dev/record_versions.R` | `Rscript dev/record_versions.R` | `output/package_versions.csv`, `output/session_info.txt`, `renv.lock` | light (about 10 s) |

`output/ladder_statistics.csv` has one row per statistic (columns ladder,
set, method, comparator, threshold, statistic, value, source, note). The rows
with method `klue_ml_screen_default` put the refit of clustering+ML at
`klue`'s default (`output/ladder_carryover_check.csv`, klue_screen TRUE) in
place of the randomised ladder's stored fit and recompute each cell's best;
Apollo's default search scored against those bests has comparator
`klue_ml_screen_default`. On the
blocked ladder the stored times of Apollo's reduced search are stale and
marked so; the times the paper reports are the seed-1 runs in
`output/apollo_smartstart_diag.csv`. `output/apollo_lcem_cap.csv` records,
per ladder cell, what `apollo_lcEM` printed in the run log: whether it
reported the iteration cap, its last EM iteration and improvement, and whether
it continued with classical estimation. The two run logs it parses,
`output/stress_replicate_blocked_std.log` and `output/stress_replicate_std.log`,
are published too; without them the script reads `output/apollo_lcem_cap.csv`
back instead of rewriting it, and the 44 `apollo_lcEM` rows of
`output/ladder_statistics.csv` then name that file as their source (the values
are the same). `dev/check_lcmnl_standard_errors.R`
supersedes `dev/benchmark_lcmnl_se.R` for the paper's statement: there Apollo
converges from its own start, so the two optima differ slightly (1.1e-5 at
C = 3, against 3e-8 at the same point).

**Collecting the results.** `Rscript dev/collect_standard_reruns.R` (base R,
no estimation, safe beside a running lane) gathers every number the reruns
replace, next to the value published in version12, into
`output/standard_rerun_summary.{csv,md}`; unfinished pieces are marked
"pending", a discrete verdict reached while a benchmark is missing is marked
"provisional", and a new MMNL fit more than 0.5 LL below the superseded 0.9.5
fit of the same model is flagged. It is bookkeeping for the revision: it
reads the pre-revision `version12.tex`, the lane logs and
`output/superseded_klue095/`, none of which are published, so it runs only in
the working repository, and its outputs are not published either.

Superseded by `dev/rerun_mmnl_benchmark.R` and not rerun: the adapters'
independent-MMNL rows in `output/*_model_comparison.csv`,
`output/empirical_corr_mmnl_*draws.csv` (`dev/empirical_corr_mmnl.R`) and
`output/vittel_respec_mmnl*.rds` (`dev/rerun_vittel_respec.R`, whose MMNL stages
now carry a `_std` suffix).

## Publishing

| Script | What it does |
|--------|--------------|
| `dev/make_reproduction_archive.sh` | Builds `klue-reproduction-v<version>.zip` (code, result files, figures, `renv.lock`, `CITATION.cff`, the licence and the exact package tarball) for the GitHub release and the Zenodo deposit. The archive holds the package as `klue_<version>.tar.gz`: run `tar xzf klue_<version>.tar.gz` in its root first, which creates the `klue/` tree every script loads. |
| `dev/publish_filter.zsh` | The rule `dev/make_reproduction_archive.sh` uses for `output/`: which top-level files are published (not run logs other than the two stress-ladder logs, resume part-files, development cross-checks, revision bookkeeping, or `output/vittel_respec_lcmnl.rds`, whose fit objects carry each Vittel respondent's class posteriors; its table is published as `output/vittel_respec_lcmnl.csv`) and which subfolders (`mmnl_bench_std/`, `mmnl_bench/`). |

The public repository keeps the package in `klue/`, so both
`remotes::install_github("o-frings/klue", subdir = "klue")` and the r-universe
registry entry (`o-frings.r-universe.dev/packages.json`, which carries
`"subdir": "klue"`) resolve it there. The sync also copies `renv.lock`, the
figures and `klue/CITATION.cff` to the public root. Neither script
redistributes `test_data_files/`.

## Quick start

```r
pkgload::load_all("klue")
source("studies/klue_studies.R")
options(klue.screen = FALSE)   # the pre-0.9.1 multistart behind the simulation tables (see "Package history")
klue_study_main(blocked = TRUE, n_reps = 20)        # blocked enumeration baseline
klue_study_estimator(blocked = TRUE)                # H1b initialisation benchmark
# Package head-to-head (needs apollo + gmnl + mlogit 1.0-3.1 in
# ~/R/oldmlogit_lib, see Requirements; single-core):
#   Rscript dev/stress_replicate_blocked.R
```

Cores: the `dev/run_*` runners and `R/run_init_ablation.R` honour `KLUE_CORES=<n>`
(they set `options(klue.cores = n)`, the drivers' condition-level parallelism knob);
unset, the drivers use `.klue_cores()` = min(16, cores − 2). Output is identical at
any core count (per-condition seeds).

Notes on reproducibility: the stress ladders run `apollo_searchStart` at Apollo
0.3.5's defaults (`smartStart = FALSE` is the default). The optional
`smartStart = TRUE` stops on the blocked cards because its candidate weights
`exp(0.5 * eigenvalue)` overflow when the log-likelihood Hessian at the start has
eigenvalues above about 1,400 (`dev/run_apollo_smartstart.R`). Apollo's cluster is
kept to a single core (it can deadlock when detached at > 1 core).

## Package history

`klue/` is the single package (0.10.0, 2026-09-29). Version 0.9 rewrote the
0.6.x engine and was verified equivalent to it (data generation bit-exact;
LCMNL LL/BIC to 1e-5 with `options(klue.screen = FALSE)`; MMNL to 1e-2; the
last run matched 15 of 15 study components; the harness and its logs,
`output/compare_klue_*`, are not published), then trimmed to the estimation library: the
`klue_study_*` drivers moved to `studies/klue_studies.R`, the 0.6.x alias names
are gone, and every script here uses the `klue_*` names. Since 0.9.1 the
multistart screens the six starts at a loose tolerance and polishes the
winner; for table-exact reproduction of pre-0.9.1 runs set
`options(klue.screen = FALSE)`, which fits every start at full precision.

**Which results need it.** Every result file written before klue 0.9.1
(2026-07-07) came from the full-precision multistart: the randomised-design
studies (`output/main_results.csv`, `convergence_results.csv`,
`estimator_orthogonal.rds`, `init_ablation.rds`, `unbalanced_results.csv`,
`design_results.csv`, `concomitant_results.csv`), the blocked-design studies
(`output/main_blocked.csv`, `convergence_blocked.csv`, `estimator_blocked.*`,
`init_ablation_blocked.*`, `{unbalanced,concomitant,sample}_blocked.csv`, and
`recovery_blocked.csv`, `ari_blocked.csv` and `firth_h3_*.csv`, which are
computed from them), and the clustering arms of both stress ladders. Reproduce
these with `options(klue.screen = FALSE)`; their runners set it
(`dev/run_blocked_baseline.R`, `dev/run_estimator_blocked.R`,
`dev/run_robustness_blocked.R`, `R/run_init_ablation.R` and the two
stress-ladder runners), and the Quick start above shows it for the
randomised-design studies.
Everything written later (the empirical LCMNL ladders, run 2026-07-15 and
committed 2026-07-24, and the
October reruns: `tab:mmnl`, `tab:mmnl_corr`, `tab:mmnl_blocked`, the H4
misspecification runs, the empirical MMNL benchmark and the empirical LCMNL
search) used the default. `dev/check_ladder_carryover.R` shows what the switch
does on the stress ladders: with it, clustering+ML reproduces every stored
log-likelihood exactly; at the default it reproduces the blocked ladder to
within 4e-5, but on the randomised ladder it ends elsewhere in 5 of the 20
cells (3 lower by 1.5 to 3.6 LL, 1 higher by 2.7, 1 lower by 0.001), so that
it would miss 5 seeds instead of 4. Clustering+EM, RP+EM and gmnl reproduce
at the default. Version 0.10.0 makes every Apollo call follow
Apollo's example scripts ("Standard-Apollo reruns" above; `klue/NEWS.md`). The
old engine and the comparison harness are kept in `archive/code_backups/` of the
working repository, which is not published. The versions of every other package are under "Requirements".
