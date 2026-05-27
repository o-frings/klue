\name{klue_engine}
\alias{klue_dgp}
\alias{klue_lcmnl}
\alias{klue_mmnl}
\alias{klue_mmnl_corr}
\alias{klue_mmnl_defaults}
\alias{make_dgp_config}
\alias{estimate_lcmnl_multistart}
\alias{estimate_mmnl}
\alias{estimate_mmnl_corr}
\title{Engine functions (advanced)}
\description{
  Low-level estimation engine used internally by \code{\link{klue}}.
  Most users will not need to call these directly.

  \itemize{
    \item \code{klue_dgp(n_generic, n_alternatives)} returns a design
      configuration (attribute names, parameter counts) consumed by the
      estimators.
    \item \code{klue_lcmnl(database, C, dgp)} estimates an LCMNL with
      \code{C} classes from six clustering-based starting values and
      keeps the best by log-likelihood. Returns a fit object with
      \code{LL}, \code{BIC}, \code{AIC}, \code{ICL}, \code{ICL_BIC},
      \code{betas}, \code{class_probs}, \code{posteriors}, and
      \code{best_method}.
    \item \code{klue_mmnl(database, dgp, ...)} estimates a mixed
      multinomial logit with independent normal random coefficients
      (log-normal on price). Returns \code{LL}, \code{BIC}, \code{AIC},
      \code{mu}, \code{sigma} on success; \code{reason} +
      \code{apollo_log_tail} on failure.
    \item \code{klue_mmnl_corr(database, dgp, ...)} estimates a mixed
      multinomial logit with a full Cholesky-parameterised covariance
      matrix over random coefficients (log-normal on price). Same return
      shape as \code{klue_mmnl} but the number of free parameters scales
      as \eqn{n_{beta}(n_{beta}+1)/2}.
    \item \code{klue_mmnl_defaults()} returns the current default
      settings for both MMNL estimators (\code{n_draws},
      \code{n_draws_stage1}, \code{draws_type}, \code{estimation_routine},
      \code{n_cores}, \code{quiet}, \code{mu_price_bounds},
      \code{sigma_price_bounds}). Override individually via the
      \code{mmnl_opts} argument of \code{\link{klue}} or with
      \code{options(klue.mmnl.*)}.
  }

  The expected \code{database} is canonical wide format: columns \code{ID},
  \code{TASK}, \code{x1_1..xN_J}, \code{price_1..price_J}, \code{CHOICE},
  with \code{CHOICE} an integer index in \code{1..J}.

  \code{make_dgp_config}, \code{estimate_lcmnl_multistart},
  \code{estimate_mmnl}, and \code{estimate_mmnl_corr} are deprecated
  aliases for the above.
}
\usage{
klue_dgp(n_generic = 4, n_alternatives = 3)
klue_lcmnl(database, C, dgp = DGP_DEFAULT)
klue_mmnl(database, n_draws = NULL, n_draws_stage1 = NULL,
          draws_type = NULL, estimation_routine = NULL,
          n_cores = NULL, quiet = NULL, mu_price_bounds = NULL,
          sigma_price_bounds = NULL, dgp = DGP_DEFAULT)
klue_mmnl_corr(database, n_draws = NULL, n_draws_stage1 = NULL,
               draws_type = NULL, estimation_routine = NULL,
               n_cores = NULL, quiet = NULL, dgp = DGP_DEFAULT)
klue_mmnl_defaults()
}
\arguments{
  \item{n_generic}{Number of generic attributes (excluding price).}
  \item{n_alternatives}{Number of alternatives \eqn{J}. Last alternative
    is the reference (no ASC).}
  \item{database}{Canonical wide data frame.}
  \item{C}{Number of latent classes (\code{C >= 1}; \code{C = 1} = MNL).}
  \item{n_draws, n_draws_stage1, draws_type, estimation_routine, n_cores,
    quiet, mu_price_bounds, sigma_price_bounds}{Tuning knobs;
    \code{NULL} means "use the current \code{klue_mmnl_defaults()}".}
  \item{dgp}{Output of \code{klue_dgp}.}
}
\seealso{\code{\link{klue}}, \code{\link{klue_database}}}
