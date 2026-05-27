\name{engine_exports}
\alias{make_dgp_config}
\alias{estimate_lcmnl_multistart}
\alias{estimate_mmnl}
\alias{estimate_mmnl_corr}
\alias{klue_mmnl_defaults}
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
    \item \code{estimate_mmnl(database, ...)} estimates a mixed
      multinomial logit with independent normal random coefficients
      (log-normal on price). Returns \code{LL}, \code{BIC}, \code{AIC},
      \code{mu}, \code{sigma}, plus \code{converged}, \code{reason},
      \code{apollo_log_tail}, \code{apollo_log_path}, and \code{settings}.
    \item \code{estimate_mmnl_corr(database, ...)} same as above with a
      fully Cholesky-correlated covariance structure.
    \item \code{klue_mmnl_defaults()} returns the active MMNL defaults as
      a named list. See \emph{MMNL settings} below.
  }

  The expected \code{database} is canonical wide format: columns \code{ID},
  \code{TASK}, \code{x1_1..xN_J}, \code{price_1..price_J}, \code{CHOICE},
  with \code{CHOICE} an integer index in \code{1..J}.
}
\section{MMNL settings}{
  All MMNL arguments default to \code{NULL}, in which case the value comes
  from \code{klue_mmnl_defaults()}, which in turn reads from \code{options()}
  (\code{klue.mmnl.<name>}) and falls back to the package built-ins:
  \tabular{lll}{
    \strong{argument} \tab \strong{built-in default} \tab \strong{option name} \cr
    \code{n_draws}            \tab \code{3000L}           \tab \code{klue.mmnl.n_draws} \cr
    \code{n_draws_stage1}     \tab \code{200L}            \tab \code{klue.mmnl.n_draws_stage1} \cr
    \code{draws_type}         \tab \code{"mlhs"}          \tab \code{klue.mmnl.draws_type} \cr
    \code{estimation_routine} \tab \code{"bgw"}           \tab \code{klue.mmnl.estimation_routine} \cr
    \code{n_cores}            \tab physical cores - 1    \tab \code{klue.mmnl.n_cores} \cr
    \code{quiet}              \tab \code{TRUE}            \tab \code{klue.mmnl.quiet} \cr
    \code{mu_price_bounds}    \tab \code{c(-5, 3)}        \tab \code{klue.mmnl.mu_price_bounds} \cr
    \code{sigma_price_bounds} \tab \code{c(-3, 1)}        \tab \code{klue.mmnl.sigma_price_bounds}
  }
  Pass \code{NA} to either bounds argument to disable the corresponding
  constraint. When \code{quiet = TRUE} Apollo's output is sunk to a tempfile
  whose path and last 40 lines are attached to the returned fit if the run
  does not converge.
}
\usage{
make_dgp_config(n_generic = 4, n_alternatives = 3)
estimate_lcmnl_multistart(database, C, dgp = DGP_DEFAULT)
estimate_mmnl(database,
              n_draws = NULL, n_draws_stage1 = NULL,
              draws_type = NULL, estimation_routine = NULL,
              n_cores = NULL, quiet = NULL,
              mu_price_bounds = NULL, sigma_price_bounds = NULL,
              dgp = DGP_DEFAULT)
estimate_mmnl_corr(database,
                   n_draws = NULL, n_draws_stage1 = NULL,
                   draws_type = NULL, estimation_routine = NULL,
                   n_cores = NULL, quiet = NULL,
                   dgp = DGP_DEFAULT)
klue_mmnl_defaults()
}
\arguments{
  \item{n_generic}{Number of generic attributes (excluding price).}
  \item{n_alternatives}{Number of alternatives \eqn{J}. Last alternative
    is the reference (no ASC).}
  \item{database}{Canonical wide data frame.}
  \item{C}{Number of latent classes (\code{C >= 1}; \code{C = 1} = MNL).}
  \item{n_draws, n_draws_stage1}{Simulation draws for the main and warm-start
    MMNL stages. \code{NULL} = use defaults.}
  \item{draws_type}{Draws-generator algorithm: one of \code{"mlhs"},
    \code{"halton"}, \code{"pmc"}, \code{"sobol"}.}
  \item{estimation_routine}{Optimiser passed to \code{apollo_estimate}: one
    of \code{"bgw"}, \code{"bfgs"}, \code{"nr"}.}
  \item{n_cores}{Number of cores for Apollo's likelihood evaluation.
    Default: physical cores minus one.}
  \item{quiet}{If \code{TRUE} (default), Apollo's stdout is captured to a
    tempfile and only surfaced (tail of 40 lines) when the run fails.}
  \item{mu_price_bounds, sigma_price_bounds}{Length-2 numeric vectors of
    lower/upper bounds for the log-normal price parameters. Pass \code{NA}
    to disable.}
  \item{dgp}{Output of \code{make_dgp_config}.}
}
\seealso{\code{\link{run_lcmnl_workflow}}, \code{\link{build_database}}}
