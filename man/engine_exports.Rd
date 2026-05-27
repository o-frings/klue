\name{engine_exports}
\alias{make_dgp_config}
\alias{estimate_lcmnl_multistart}
\alias{estimate_mmnl}
\title{Engine functions (advanced)}
\description{
  Low-level estimation engine used internally by
  \code{\link{run_lcmnl_workflow}}. Most users will not need to call these
  directly.

  \itemize{
    \item \code{make_dgp_config(n_generic, n_alternatives)} returns a
      design configuration (attribute names, parameter counts) consumed by
      the estimators.
    \item \code{estimate_lcmnl_multistart(database, C, dgp)} estimates an
      LCMNL with \code{C} classes from six clustering-based starting values,
      keeps the best by log-likelihood. Returns a fit object with
      \code{LL}, \code{BIC}, \code{AIC}, \code{ICL}, \code{ICL_BIC},
      \code{betas}, \code{class_probs}, \code{posteriors}, and
      \code{best_method}.
    \item \code{estimate_mmnl(database, dgp)} estimates a mixed
      multinomial logit with independent normal random coefficients
      (log-normal on price). Returns \code{LL}, \code{BIC}, \code{AIC},
      \code{mu}, \code{sigma}.
  }

  The expected \code{database} is canonical wide format: columns \code{ID},
  \code{TASK}, \code{x1_1..xN_J}, \code{price_1..price_J}, \code{CHOICE},
  with \code{CHOICE} an integer index in \code{1..J}.
}
\usage{
make_dgp_config(n_generic = 4, n_alternatives = 3)
estimate_lcmnl_multistart(database, C, dgp = DGP_DEFAULT)
estimate_mmnl(database, n_draws = N_DRAWS_MMNL, n_draws_stage1 = 200L,
              dgp = DGP_DEFAULT)
}
\arguments{
  \item{n_generic}{Number of generic attributes (excluding price).}
  \item{n_alternatives}{Number of alternatives \eqn{J}. Last alternative
    is the reference (no ASC).}
  \item{database}{Canonical wide data frame.}
  \item{C}{Number of latent classes (\code{C >= 1}; \code{C = 1} = MNL).}
  \item{n_draws, n_draws_stage1}{Number of simulation draws for MMNL.}
  \item{dgp}{Output of \code{make_dgp_config}.}
}
\seealso{\code{\link{run_lcmnl_workflow}}, \code{\link{build_database}}}
