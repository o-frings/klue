# klue

Your clue to K. An R package for specifying latent class multinomial logit
(LCMNL) models, plus the full reproduction materials for the paper it
implements:

> Frings, O. (2026). *A Clustering-Initialised Specification Workflow for
> Latent Class Choice Models*. Working paper.

## Layout

| Path | What it is |
|------|------------|
| `klue/` | The R package: clustering-initialised LCMNL estimation, MMNL benchmark, BIC/AIC/ICL and entropy diagnostics. See `klue/README.md` for the package tour and `klue/NEWS.md` for the changelog. |
| `README_REPRODUCE.md` | The reproduction guide: which script produces every table and figure in the paper, and where the result files live. |
| `studies/`, `dev/`, `R/`, `output/` | Simulation-study drivers, orchestration scripts, the five empirical applications, and the result files the paper reports. |

## Install the package

```r
install.packages("klue",
                 repos = c("https://o-frings.r-universe.dev",
                           "https://cloud.r-project.org"))
# or from this repository:
remotes::install_github("o-frings/klue", subdir = "klue")
```

## Quick start

```r
library(klue)
klue_demo()   # full workflow on a simulated example, ~3 s
```

## Citing

See `CITATION.cff`, or run `citation("klue")` for the full BibTeX bundle
(the methodology paper, the package, and the upstream packages it builds on).

## License

MIT (see `LICENSE`). No restricted-access data are redistributed here; the
Data section of `README_REPRODUCE.md` says where to obtain each dataset.
