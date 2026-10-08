# klue

Your clue to K. A reusable workflow for specifying latent class multinomial
logit (LCMNL) models. Implements the hybrid ML / random-utility framework
from Frings (2026), *A Clustering-Initialised Specification Workflow for
Latent Class Choice Models*.

---

## Install

```r
install.packages("klue",
                 repos = c("https://o-frings.r-universe.dev",
                           "https://cloud.r-project.org"))
```

(Or from a tarball: `install.packages("klue_<version>.tar.gz", repos = NULL, type = "source")`.)

---

## 1. See it work — one line, no setup

```r
library(klue)
klue_demo()
```

Runs the full workflow on a bundled example dataset (~3 seconds) and prints
the BIC/AIC/ICL summary plus class-specific coefficients. Use this once
to see what the output looks like before plugging in your own data.

For the full version (C = 1..6 + MMNL benchmark with 500 draws, a few minutes): `klue_demo(full = TRUE)`.

---

## 2. Use it on your own data

The package needs to know which column in your data is which. Two common
data layouts are supported:

### Long format — one row per (respondent × task × alternative)

```r
library(klue)

res <- klue(
  data           = "my_data.csv",          # CSV path or data.frame
  id_col         = "respondent_id",        # who made the choice
  task_col       = "task",                 # which choice task
  alt_col        = "alternative",          # which alternative
  choice_col     = "chosen",               # 0/1 indicator of the picked row
  attribute_cols = c("attr1", "attr2"),    # generic attributes
  price_col      = "price"                 # the price/cost attribute
)
```

That's all you need. Defaults: estimates `C = 1..6`, runs MMNL, writes 3
CSVs to `output/`, returns a results list invisibly.

### Wide format — one row per (respondent × task), attributes by alternative

```r
res <- klue(
  data           = my_data,
  format         = "wide",
  id_col         = "ID",
  task_col       = "task",                 # or NULL to auto-number
  choice_col     = "choice",               # integer 1..J of chosen alt
  attribute_cols = list(
    time = c("time_alt1", "time_alt2", "time_alt3"),
    qual = c("qual_alt1", "qual_alt2", "qual_alt3")
  ),
  price_col      = c("cost_alt1", "cost_alt2", "cost_alt3"),
  avail_col      = c("av_alt1", "av_alt2", "av_alt3")   # optional
)
```

Use `NA` in an attribute slot to encode a structural zero on that
alternative (e.g. some attributes only apply to some alternatives).

Long and wide use the **same argument names** (`attribute_cols`,
`price_col`, `avail_col`); only the *type* differs (scalar in long,
length-J vector or named list in wide).

---

## 2b. Simulate data instead (Monte Carlo / methodology testing)

The data-generating process used in the Frings (2026) Monte Carlo study
is also exposed:

```r
sim <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 2,
                     separation = 1.0, heterogeneity = 0.25, seed = 42)

# Returned: sim$database (canonical wide format) + sim$true_betas, sim$true_class
res <- klue(database = sim$database, C_cands = 1:4)
res$best_C    # should recover 2 here
```

Variant: `klue_simulate(covariates = TRUE)` drives class membership by covariates.
Full reference: `?klue_simulate`.

The Monte Carlo study from Frings (2026) is not part of the package. Its
drivers — `klue_study()` (the 420-condition main simulation plus the
robustness analyses, ~12 hours) and the per-component `klue_study_main()`,
`klue_study_unbalanced()`, `klue_study_design()`, `klue_study_clustering()`,
etc. — live in `studies/klue_studies.R` of the reproduction repository:
load the package, then `source()` that file (see its `README_REPRODUCE.md`).

---

## 3. Inspect the results

```r
res$summary          # one row per C: LL, k, BIC, AIC, ICL, ΔBIC, best clustering method
res$class_betas      # class shares and coefficients for the BIC-best model
res$comparison       # MNL (C=1) vs LCMNL (BIC-best) vs MMNL
res$best_C           # the BIC-best number of classes
res$best_lcmnl       # the fit itself (posteriors, betas, etc.)
```

Three CSVs are written to `output/` by default (override with `output_dir`):

- `workflow_summary.csv`
- `workflow_class_betas.csv`
- `workflow_model_comparison.csv`

---

## Optional refinements

| Argument | What it does |
|---|---|
| `price_scaling = 10` | Divide price by 10 before estimation (numerical stability). |
| `scalings = list(time = 60, ...)` | Per-attribute scaling (e.g. seconds → minutes). |
| `avail_col` (long) / `availability` (wide) | Drop tasks where some alternative is unavailable. |
| `C_cands = 1:4` | Pick a different range of class counts. |
| `run_mmnl = FALSE` | Skip the independent-normals MMNL benchmark (faster). |
| `run_mmnl_corr = TRUE` | Additionally estimate a correlated-normals MMNL (full Cholesky covariance). Slower; tests whether heterogeneity is genuinely correlated across attributes. |
| `output_prefix = "myrun"` | Prefix for the output CSV filenames. |
| `output_dir = "results"` | Write CSVs somewhere other than `output/`. |
| `write_csv = FALSE` | Return results in memory only. |

Full reference: `?klue`. The package exports the `klue_*` names only; the paper's Monte Carlo drivers live in `../studies/klue_studies.R`.

---

## Availability filtering

The estimation engine assumes every alternative is available in every task.
If your data has availability columns, pass them via `avail_col` (long
format) or `availability` (wide format). The workflow filters to
fully-available tasks and reports the drop rate. Partial-availability
estimation is not yet supported.

## Validation

Five reference applications in Frings (2026) — Rhin-Meuse water-quality DCE,
Apollo mode/route choice, Electricity (Train 1998), Swissmetro — reproduce
bit-exactly through this workflow against hand-coded baselines. See
`R/test_klue.R` in the source tree for the test.

## Citing

If you use klue in published work, please cite **the methodology paper,
the software itself, and the upstream packages it builds on**:

- Frings, O. (2026). *A Clustering-Initialised Specification Workflow for
  Latent Class Choice Models*. Working paper.
- Frings, O. (2026). *klue: Hybrid Machine Learning and Random-Utility
  Workflow for Latent Class Multinomial Logit Model Specification*.
  R package. https://github.com/o-frings/klue
- Hess, S., & Palma, D. (2019). *Apollo: a flexible, powerful and
  customisable freeware package for choice model estimation and
  application*. Journal of Choice Modelling, 32, 100170.
  [doi:10.1016/j.jocm.2019.100170](https://doi.org/10.1016/j.jocm.2019.100170)
- Scrucca, L., Fop, M., Murphy, T.B., & Raftery, A.E. (2016). *mclust 5:
  clustering, classification and density estimation using Gaussian finite
  mixture models*. The R Journal, 8(1), 289-317.
  [doi:10.32614/RJ-2016-021](https://doi.org/10.32614/RJ-2016-021)
- Maechler, M., Rousseeuw, P., Struyf, A., Hubert, M., Hornik, K. (2024).
  *cluster: Cluster Analysis Basics and Extensions*. R package.

```r
citation("klue")    # returns all five BibTeX entries (auto-updates version)
```

## Acknowledgements

klue wraps the [Apollo](http://www.ApolloChoiceModelling.com/) choice-model
estimation engine (Hess & Palma 2019); two of the six starting-value
clusterings come from [mclust](https://mclust-org.github.io/mclust/)
(Scrucca et al. 2016, GMM) and [cluster](https://cran.r-project.org/package=cluster)
(Maechler et al., PAM). The hybrid workflow and the diagnostics layer are
klue's contribution; the underlying MLE machinery and the clustering
algorithms are not.

---

# Version history

## 0.9.3 — simplification (2026-09-17)

- The package is the estimation library only: `klue_database_*`, `klue_dgp`, `estimate_lcmnl`
  (+ `_cov`, `_em`), `klue_lcmnl`, `klue()`, `klue_mmnl` (+ `_corr`), the `klue_*` inference
  helpers, `klue_simulate`, `klue_simulate_choices`, `klue_design`, `klue_starts`, `klue_demo`. 22 exports.
- The paper's Monte Carlo drivers (`klue_study_*`, `klue_study`), the recovery metrics and the
  legacy fixed-design simulator moved to `../studies/klue_studies.R` (source it after loading
  the package). The covariate-DGP sugar wrapper was folded into `klue_simulate(covariates = TRUE)`.
- The 0.6.x compatibility aliases (`run_*`, `build_*`, `generate_*`, `estimate_*multistart*`,
  `make_dgp_config`) are removed; the scripts in `R/` and `dev/` were renamed accordingly.
- `klue_database_long()` orders respondents and tasks numerically (no collation dependence),
  takes `covariate_cols` and sets `attr(db, "ids")`; fits carry a `singular` flag with a message
  when the information matrix cannot be inverted.

This is the rebuilt klue: same exported API as the 0.6.x engine in
`../klue/`, at roughly 70% of the line count. User-facing documentation: see
`../klue/README.md` (everything there applies unchanged). Numerical behaviour
matched 0.6.x exactly up to 0.9.0; 0.9.1 changes the multistart mechanics
(below) — it targets the same optimum, but tie-breaks between equal-LL starts
can differ.

## 0.9.1 — large-N fixes (2026-07-07)

Motivated by a runaway run at N = 6,000: 55+ minutes stuck in the O(N^2)
clustering-starts phase with nothing written.

- `klue_starts()` and the multistart cluster a fixed-seed respondent
  subsample capped at `cluster_cap = 2000` (pam, Mclust, and hclust are
  O(N^2) or worse; starting values do not need every respondent). The LCMNL
  itself still fits the full sample.
- Two-stage multistart in `klue_lcmnl()`: all starts are screened at
  `maxit = 200, reltol = 1e-6`; only the best log-likelihood is polished at
  full precision, warm-started via the new `estimate_lcmnl(start_par =)`.
  Fits that hit the iteration cap keep their finite LL with
  `converged = FALSE` instead of being discarded.
- Fit objects return the full parameter vector as `fit$par` (betas, per-class
  ASCs, share deltas) — needed for parameter-level class criteria.
- Each start logs its LL and elapsed time via `message()`, so long runs are
  observable from a redirected console.

Net effect at N = 6,000 (17 attributes, 4 alternatives, C = 2): ~30 s per
multistart on 6 cores, previously unbounded.

## 0.9.2 — hardening (2026-07-07)

Principle applied throughout: **no silent fallbacks** — every degraded path
now warns, messages, or stops.

- apollo moved from Depends to Suggests: only `klue_mmnl()` needs it (hard
  `stop()` with install hint if absent); `library(klue)` no longer attaches
  apollo, so LCMNL-only pipelines load faster and without apollo's warnings.
- A dead `mclapply` worker (try-error in the result slot) now raises a
  warning naming the lost start instead of crashing the winner scan.
- Deterministic seeding (clustering, subsampling, baseline start generators)
  saves and restores `.Random.seed`, so klue no longer resets the caller's
  RNG stream.
- `C = 1` short-circuits to a single fit (the six clustering starts coincide
  there; previously six identical fits ran).
- `gc()` before forking workers shrinks the copy-on-write heap children
  inherit.
- `compute_recovery()` refuses K > 8 loudly instead of attempting a K^K
  permutation enumeration.
- Loud fallbacks: cluster-start fallback coefficients, EM collapsed-class /
  failed M-step skips, and all MMNL warm-start fallbacks (failed pooled MNL,
  non-negative price coefficient, unusable stage-1 estimates) now report
  what happened and what is used instead.
- `estimate_lcmnl_em()` returns `fit$par` like the ML path;
  `apollo_probabilities` cleanup is registered before the global assignment,
  so an apollo error cannot leak it into the global environment.

## What changed internally

- One shared MNL log-likelihood/gradient kernel (`R/03_estimate.R`) serves
  the pooled MNL, every class inside the C >= 2 likelihood, the EM M-step,
  and the cluster-wise start fits.
- Sequencing vs parallelisation is an explicit knob throughout:
  `klue()` / `klue_lcmnl(..., n_cores =)` fit the six clustering starts
  concurrently (~3x on a multicore box), and every `klue_study_*` driver plus
  the master `klue_study(..., n_cores =)` (in `../studies/klue_studies.R`)
  parallelise across conditions.
  `n_cores = 1` is fully sequential (lowest memory, shared-machine-friendly);
  higher uses more cores. Output is identical either way (per-condition seeds /
  order-independent multistart selection). Study drivers default to
  `.klue_cores()` = `min(16, cores - 2)`; the interactive entry points default
  to 1 so they never nest inside the drivers' parallelism.
- A C++ (RcppArmadillo) kernel was built and verified exact but is ~20%
  *slower* than this BLAS-bound R kernel, so it was parked in
  `../dev/cpp_experiment/` rather than shipped; see
  `../dev/cpp_kernel_findings.md` and `../dev/perf_levers.md`.
- One simulator (`klue_simulate`, with `covariates =` and `design =`)
  replaces three near-duplicates; every code path consumes the RNG stream in
  the same order as 0.6.x, so seeded datasets are bit-identical.
- The six clustering-start functions collapsed into a registry behind
  `klue_starts(method =, feature_type = c("rp", "onehot"))`.
- `klue_mmnl(correlation = TRUE/FALSE)` replaces the duplicated
  independent/correlated MMNL pair. Apollo components are passed to
  `apollo_validateInputs()` as explicit arguments; `apollo_beta`,
  `apollo_fixed` and `apollo_probabilities` go through the global environment
  during estimation, as in Apollo's examples, and the caller's own objects are
  restored afterwards (since 0.10.0; see NEWS.md).
- The 12 study drivers share grid/seed/mclapply/assembly plumbing
  (`R/08_study.R`); per-condition seeds and grids are unchanged.
- (0.9.3) The 0.6.x alias names are gone; use the `klue_*` names.

