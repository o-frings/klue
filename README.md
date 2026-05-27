# klue

Your clue to K. A reusable workflow for specifying latent class multinomial
logit (LCMNL) models. Implements the hybrid ML / random-utility framework
from Frings (2026), *A Hybrid Machine Learning and Random Utility Framework
for Latent Class Model Specification*.

## What it does

Given any discrete-choice dataset and a column mapping, the package:

1. Initialises the LCMNL optimiser from six clusterings of respondents'
   revealed-preference signatures (kmeans, GMM, three flavours of
   hierarchical clustering, PAM).
2. Estimates LCMNL for a user-specified range of class counts (`C_cands`).
3. Estimates an MMNL benchmark (independent normals, log-normal on price).
4. Returns BIC / AIC / ICL / entropy summaries and class-specific
   coefficients for the BIC-best model.
5. Writes a tidy bundle of CSVs to disk.

## Install

Pick whichever fits the recipient's setup:

**From a built tarball (easiest — share `klue_0.1.0.tar.gz`):**
```r
install.packages("klue_0.1.0.tar.gz", repos = NULL, type = "source")
# or, with the remotes package:
remotes::install_local("klue_0.1.0.tar.gz")
```

**From a local source folder:**
```r
remotes::install_local("klue")        # path to the klue/ folder
# or
devtools::install("klue")
```

**From GitHub (once pushed):**
```r
remotes::install_github("oliverfrings/klue")
```

The package depends on `apollo`, `mclust`, and `cluster`. They install
automatically with the methods above.

## Quick start

### Long-format data

```r
library(klue)

res <- run_lcmnl_workflow(
  data           = "path/to/long_data.csv",
  format         = "long",
  id_col         = "respondent_id",
  task_col       = "task",
  alt_col        = "alternative",
  choice_col     = "chosen",            # 0/1 indicator
  attribute_cols = c("attr1", "attr2", "attr3"),
  price_col      = "price",
  price_scaling  = 10,
  C_cands        = 1:6
)

res$summary       # one row per C with LL / BIC / AIC / ICL
res$class_betas   # class-specific coefficients for BIC-best model
res$comparison    # MNL vs LCMNL vs MMNL
res$best_C        # BIC-best number of classes
```

### Wide-format data

```r
data(apollo_modeChoiceData, package = "apollo")
d <- apollo_modeChoiceData[apollo_modeChoiceData$SP == 1, ]

res <- run_lcmnl_workflow(
  data         = d,
  format       = "wide",
  id_col       = "ID", task_col = "SP_task", choice_col = "choice",
  attributes   = list(
    time    = c("time_car", "time_bus", "time_air", "time_rail"),
    access  = c(NA,         "access_bus", "access_air", "access_rail"),  # NA => structural 0
    service = c(NA, NA, "service_air", "service_rail")
  ),
  price        = c("cost_car", "cost_bus", "cost_air", "cost_rail"),
  availability = c("av_car",   "av_bus",   "av_air",   "av_rail"),
  scalings     = list(time = 60, access = 60, price = 10),
  C_cands      = 1:6
)
```

## Outputs

`run_lcmnl_workflow()` writes (when `write_csv = TRUE`, default):

- `<prefix>_summary.csv` — one row per C: LL, k, BIC, AIC, ICL, ΔBIC, ΔAIC, ΔICL, best clustering method.
- `<prefix>_class_betas.csv` — class shares and per-attribute coefficients for the BIC-best model.
- `<prefix>_model_comparison.csv` — MNL (C=1) vs LCMNL (BIC-best) vs MMNL.

## Availability filtering

The estimation engine assumes every alternative is available in every task
(balanced panel × balanced choice sets). If your data has availability
columns, pass them via `avail_col` (long format) or `availability` (wide
format). The workflow filters the panel to fully-available tasks and reports
the drop rate. Partial-availability estimation is not yet supported.

## Validation

All five reference applications in Frings (2026) — Vittel water-quality DCE,
Apollo mode/route choice, Electricity (Train 1998), Swissmetro — reproduce
bit-exactly through this workflow against their hand-coded original
adapters. See `R/test_klue.R` in the source tree for the test.

## Citing

If you use klue in published work, please cite the klue paper **and** the
packages it builds on:

- Frings, O. (2026). *A Hybrid Machine Learning and Random Utility Framework
  for Latent Class Model Specification*. Journal of Choice Modelling.
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
