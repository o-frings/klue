\name{klue-package}
\alias{klue-package}
\docType{package}
\title{klue: Hybrid Machine Learning and Random-Utility Workflow for Latent Class Multinomial Logit Model Specification}
\description{
  Implements the three-step workflow from Frings (2026, working paper)
  for specifying latent class multinomial logit (LCMNL) models:
  initialise the maximum-likelihood multi-start from six clusterings of
  respondents' revealed-preference signatures (k-means, Gaussian mixture,
  hierarchical clustering with Ward / complete / average linkage, and
  partitioning around medoids); estimate LCMNL across a user-specified
  range of class counts; and benchmark against a mixed multinomial logit
  (MMNL, independent and optionally correlated normals) with BIC, AIC,
  ICL, and a classification-entropy diagnostic.
}
\section{Public API by topic}{
  \describe{
    \item{Entry point}{
      \code{\link{klue}} -- main workflow runner.
      \code{\link{klue_demo}} -- zero-setup demo on bundled data.
    }
    \item{Data harmonising}{
      \code{\link{klue_database}}, \code{klue_database_long},
      \code{klue_database_wide} -- reshape long- or wide-format
      discrete-choice data to the canonical format the engine consumes;
      handle availability filtering and per-attribute scaling.
    }
    \item{Data simulation (Monte Carlo)}{
      \code{\link{klue_simulate}}, \code{klue_simulate_cov},
      \code{klue_simulate_deff} -- generate synthetic panel choice data
      from a known LCMNL DGP. The standard random-attribute design;
      a variant with concomitant covariates driving class membership;
      and a D-efficient design variant.
    }
    \item{Estimators (advanced)}{
      \code{\link{klue_lcmnl}} -- LCMNL multi-start.
      \code{\link{klue_mmnl}} -- MMNL with independent normals.
      \code{\link{klue_mmnl_corr}} -- MMNL with correlated normals
      (full Cholesky covariance).
    }
    \item{Configuration}{
      \code{\link{klue_dgp}} -- design / parameter config (number of
      attributes, alternatives, mean-beta vector).
      \code{\link{klue_mmnl_defaults}} -- active default settings for
      both MMNL estimators (draws, cores, routine, bounds).
    }
    \item{Paper replication}{
      \code{\link{klue_study}} -- runs all 11 component drivers below
      (~12 h end-to-end).
      \code{klue_study_main} -- 420-condition Monte Carlo.
      \code{klue_study_mmnl}, \code{klue_study_mmnl_corr} -- MMNL
      benchmark and correlated-MMNL robustness.
      \code{klue_study_convergence}, \code{klue_study_initialisation}
      -- starting-value ablations.
      \code{klue_study_unbalanced}, \code{klue_study_design},
      \code{klue_study_concomitant}, \code{klue_study_clustering},
      \code{klue_study_recovery}, \code{klue_study_sample}
      -- robustness analyses across nuisance design factors.
    }
  }
}
\section{Quick start}{
  \preformatted{
    library(klue)
    klue_demo()                              # see it work
    res <- klue(data = "my.csv", format = "long",
                id_col = "id", task_col = "task",
                alt_col = "alt", choice_col = "chosen",
                attribute_cols = c("a", "b"), price_col = "price")
    res$summary; res$best_C; res$class_betas
  }
}
\section{Backward compatibility}{
  Every \code{klue_*} function has a deprecated alias preserving the
  pre-v0.4.0 name (\code{run_lcmnl_workflow}, \code{build_database*},
  \code{estimate_*}, \code{make_dgp_config}, \code{generate_data*},
  \code{run_*}). Same function object, slated for removal in a future
  major release.
}
\section{Citation}{
  See \code{citation("klue")} for the BibTeX bundle: the klue paper plus
  \pkg{apollo} (Hess & Palma 2019; the underlying estimator),
  \pkg{mclust} (Scrucca et al. 2016; GMM clustering) and \pkg{cluster}
  (Maechler et al.; PAM).
}
\section{Validation}{
  The hybrid workflow has been validated bit-exactly against five
  public reference applications (Vittel water-quality discrete choice
  experiment, Apollo mode and route choice, Electricity, Swissmetro).
}
\author{
  Oliver Frings \email{oliver.frings@agroparistech.fr}
}
\references{
  Frings, O. (2026). \emph{A Hybrid Machine Learning and Random Utility
  Framework for Latent Class Model Specification}. Working paper.

  Hess, S., & Palma, D. (2019). Apollo: a flexible, powerful and
  customisable freeware package for choice model estimation and
  application. \emph{Journal of Choice Modelling}, 32, 100170.
}
\keyword{package}
