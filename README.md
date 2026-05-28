# klue

Your clue to K. A reusable workflow for specifying latent class multinomial
logit (LCMNL) models. Implements the hybrid ML / random-utility framework
from Frings (2026), *A Hybrid Machine Learning and Random Utility Framework
for Latent Class Model Specification*.

---

## Install

```r
install.packages("klue",
                 repos = c("https://o-frings.r-universe.dev",
                           "https://cloud.r-project.org"))
```

(Or from a tarball: `install.packages("klue_0.1.0.tar.gz", repos = NULL, type = "source")`.)

---

## 1. See it work — one line, no setup

```r
library(klue)
klue_demo()
```

Runs the full workflow on a bundled example dataset (~3 seconds) and prints
the BIC/AIC/ICL summary plus class-specific coefficients. Use this once
to see what the output looks like before plugging in your own data.

For the full version (C = 1..6 + MMNL benchmark, ~30 sec): `klue_demo(full = TRUE)`.

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

Variants: `klue_simulate_cov()` (with concomitant covariates driving
class membership) and `klue_simulate_deff()` (D-efficient design).
Full reference: `?klue_simulate`.

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

Full reference: `?klue`. Old names (`run_lcmnl_workflow`, `build_database*`, `estimate_*`) are kept as silent aliases for backward compatibility.

---

## Availability filtering

The estimation engine assumes every alternative is available in every task.
If your data has availability columns, pass them via `avail_col` (long
format) or `availability` (wide format). The workflow filters to
fully-available tasks and reports the drop rate. Partial-availability
estimation is not yet supported.

## Validation

Five reference applications in Frings (2026) — Vittel water-quality DCE,
Apollo mode/route choice, Electricity (Train 1998), Swissmetro — reproduce
bit-exactly through this workflow against hand-coded baselines. See
`R/test_klue.R` in the source tree for the test.

## Citing

If you use klue in published work, please cite the klue paper **and** the
packages it builds on:

- Frings, O. (2026). *A Hybrid Machine Learning and Random Utility Framework
  for Latent Class Model Specification*. Working paper.
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
citation("klue")    # returns all four BibTeX entries
```

## Acknowledgements

klue wraps the [Apollo](http://www.ApolloChoiceModelling.com/) choice-model
estimation engine (Hess & Palma 2019); two of the six starting-value
clusterings come from [mclust](https://mclust-org.github.io/mclust/)
(Scrucca et al. 2016, GMM) and [cluster](https://cran.r-project.org/package=cluster)
(Maechler et al., PAM). The hybrid workflow and the diagnostics layer are
klue's contribution; the underlying MLE machinery and the clustering
algorithms are not.
