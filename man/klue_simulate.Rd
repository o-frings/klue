\name{klue_simulate}
\alias{klue_simulate}
\alias{klue_simulate_cov}
\alias{klue_simulate_deff}
\alias{generate_data}
\alias{generate_data_with_covariates}
\alias{generate_data_defficient}
\title{Simulate panel discrete-choice data}
\description{
  Generate synthetic panel choice data from a known latent-class DGP.
  The returned \code{database} element is already in canonical wide
  format and can be passed directly to \code{\link{klue}} for
  estimation. This is the data-generating process used in the
  \eqn{>400}-condition Monte Carlo study in Frings (2026, working paper).

  Three flavours:
  \itemize{
    \item \code{klue_simulate} -- the standard random-attribute design.
    \item \code{klue_simulate_cov} -- same design plus two concomitant
      covariates (\eqn{Z_1 \sim N(0,1)}, \eqn{Z_2 \sim Bernoulli(0.5)})
      driving class membership; useful for testing covariate-aware LCMNL.
    \item \code{klue_simulate_deff} -- D-efficient experimental design
      (lower observational noise per task than the random design).
  }

  Old names \code{generate_data}, \code{generate_data_with_covariates}
  and \code{generate_data_defficient} are deprecated aliases.
}
\usage{
klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 2,
              separation = 1.0, heterogeneity = 0.25,
              seed = 12345, class_proportions = NULL,
              dgp = DGP_DEFAULT, sep_profile = NULL)

klue_simulate_cov(N_per_class = 150, T_tasks = 20, true_K = 2,
                  separation = 1.0, heterogeneity = 0.25,
                  seed = 12345, covariate_strength = 1.0,
                  dgp = DGP_DEFAULT)

klue_simulate_deff(N_per_class = 150, T_tasks = 20, true_K = 2,
                   separation = 1.0, heterogeneity = 0.25,
                   seed = 12345, dgp = DGP_DEFAULT)
}
\arguments{
  \item{N_per_class}{Respondents per latent class. Total \eqn{N} =
    \code{N_per_class * true_K} (unless \code{class_proportions} is set).}
  \item{T_tasks}{Choice tasks per respondent.}
  \item{true_K}{Number of latent classes in the data-generating process.}
  \item{separation}{Inter-class distance multiplier (\eqn{\kappa} in the
    paper). Larger values make classes easier to recover.}
  \item{heterogeneity}{Within-class taste-variation s.d. (\eqn{\sigma}
    in the paper). Larger values blur the discrete-class structure.}
  \item{seed}{RNG seed for reproducibility.}
  \item{class_proportions}{Optional length-\code{true_K} vector of
    proportions summing to 1; default = balanced. (\code{klue_simulate}
    only.)}
  \item{covariate_strength}{Strength of the covariate-to-class mapping
    in \code{klue_simulate_cov}. \code{0} = covariates are noise; large
    values make class membership nearly deterministic in \eqn{(Z_1, Z_2)}.}
  \item{dgp}{Design configuration from \code{\link{klue_dgp}} (number of
    generic attributes, alternatives, mean-beta vector).}
  \item{sep_profile}{Optional per-attribute weights for the separation
    vector; \code{NULL} = uniform across attributes.}
}
\value{
  A list with elements:
  \item{database}{Canonical wide data frame ready for \code{\link{klue}}:
    \code{ID}, \code{TASK}, \code{x1_1..xN_J}, \code{price_1..price_J},
    \code{CHOICE}.}
  \item{true_betas}{\code{true_K} x \code{n_beta} matrix of class means.}
  \item{true_class}{Length-\eqn{N} integer vector of true class memberships.}
  \item{individual_betas}{\eqn{N} x \code{n_beta} matrix of person-level taste vectors.}
  \item{N, T, K, dgp}{Echoed dimensions and config.}

  \code{klue_simulate_cov} additionally returns the simulated covariates
  \code{Z1}, \code{Z2}.
}
\examples{
\dontrun{
# Generate 300 respondents from 2 classes, 20 tasks each, moderate separation
sim <- klue_simulate(N_per_class = 150, T_tasks = 20, true_K = 2,
                     separation = 1.0, heterogeneity = 0.25, seed = 42)

# Run the workflow on the simulated data
res <- klue(database = sim$database, C_cands = 1:4, run_mmnl = TRUE)
res$summary       # BIC should pick C = 2
res$best_C
}
}
\seealso{\code{\link{klue}}, \code{\link{klue_dgp}}}
