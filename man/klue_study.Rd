\name{klue_study}
\alias{klue_study}
\alias{klue_study_main}
\alias{klue_study_mmnl}
\alias{klue_study_convergence}
\alias{klue_study_initialisation}
\alias{klue_study_unbalanced}
\alias{klue_study_design}
\alias{klue_study_concomitant}
\alias{klue_study_recovery}
\alias{klue_study_clustering}
\alias{klue_study_sample}
\alias{klue_study_mmnl_corr}
\alias{run_full_study}
\alias{run_main_simulation}
\alias{run_mmnl_comparison}
\alias{run_convergence_ablation}
\alias{run_initialisation_ablation}
\alias{run_unbalanced_analysis}
\alias{run_design_comparison}
\alias{run_concomitant_analysis}
\alias{run_unconditional_recovery}
\alias{run_clustering_comparison}
\alias{run_sample_sensitivity}
\alias{run_correlated_mmnl_robustness}
\title{Replicate the Monte Carlo simulation study from Frings (2026)}
\description{
  The \code{klue_study*()} family reproduces the simulation study reported
  in Frings (2026, working paper). \code{klue_study()} is the umbrella
  driver that runs the 11 component analyses below in sequence; each
  component is also callable standalone for partial replication.

  All results are written to \code{getOption("klue.output_dir", "output")}
  as CSV files (one per analysis) and returned as a named list.

  \strong{Runtime warning}: \code{klue_study()} takes ~12 hours on a
  modern laptop. Individual components range from ~5 min
  (\code{klue_study_unbalanced}, \code{klue_study_design}) to ~2-3 hours
  (\code{klue_study_main}, \code{klue_study_mmnl_corr}).

  The \code{run_*} aliases (e.g. \code{run_full_study},
  \code{run_main_simulation}) are silent aliases retained for backward
  compatibility with the paper's original simulation scripts.
}
\usage{
klue_study(run_main = TRUE, run_mmnl = TRUE, run_supp = TRUE,
           verbose = TRUE, dgp = DGP_DEFAULT)

klue_study_main(true_K_values = c(1, 2, 3, 4, 5),
                kappa_values = c(0.5, 0.75, 1.0, 1.25, 1.5),
                sigma_values = c(0.1, 0.15, 0.2, 0.25),
                n_reps = 5, C_cands = 1:6,
                dgp = DGP_DEFAULT, sep_profile = NULL, verbose = TRUE)

klue_study_mmnl(n_cond = 80, n_draws = N_DRAWS_MMNL,
                dgp = DGP_DEFAULT, verbose = TRUE)

klue_study_convergence(n_random = 50, n_cond = 40,
                       verbose = TRUE, dgp = DGP_DEFAULT)

klue_study_initialisation(n_random = 50, n_cond = 40,
                          verbose = TRUE, dgp = DGP_DEFAULT)

klue_study_unbalanced(verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_design(verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_concomitant(verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_recovery(n_cond = 80, verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_clustering(verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_sample(verbose = TRUE, dgp = DGP_DEFAULT)
klue_study_mmnl_corr(n_draws = N_DRAWS_MMNL, verbose = TRUE,
                     dgp = DGP_DEFAULT)
}
\arguments{
  \item{run_main, run_mmnl, run_supp}{Logical switches for which blocks
    of \code{klue_study()} to execute. \code{run_supp} controls the
    nine supplementary robustness analyses.}
  \item{true_K_values, kappa_values, sigma_values, n_reps}{Main simulation
    factor grid. Default is the 420-condition design of the paper:
    5 true K's x 5 separations x 4 heterogeneities x 5 reps + the K=1
    cell. \code{n_reps} replications per cell.}
  \item{n_cond}{Number of conditions sampled in the MMNL-comparison and
    recovery analyses.}
  \item{n_draws}{Number of simulation draws for MMNL fits.}
  \item{n_random}{Number of random starting-value restarts in the
    convergence / initialisation ablations.}
  \item{C_cands}{Class counts to evaluate per condition.}
  \item{dgp}{Design configuration (see \code{\link{klue_dgp}}). Default
    is 4 generic attributes + price, 3 alternatives.}
  \item{sep_profile}{Optional per-attribute separation weights; \code{NULL}
    = uniform.}
  \item{verbose}{Print progress lines.}
}
\value{
  Each driver returns a \code{data.frame} (or list of data frames for
  \code{klue_study()}) with one row per condition / iteration.
  Columns vary by analysis but always include the factor levels that
  define the condition plus the BIC-selected K, log-likelihood, and
  any analysis-specific diagnostics.
}
\section{Files written}{
  Under \code{getOption("klue.output_dir", "output")}:
  \itemize{
    \item \code{main_results.csv}
    \item \code{mmnl_results.csv}
    \item \code{convergence_results.csv}
    \item \code{uninformed_convergence_results.csv}
    \item \code{unbalanced_results.csv}
    \item \code{design_results.csv}
    \item \code{concomitant_results.csv}
    \item \code{recovery_results.csv}
    \item \code{clustering_comparison_results.csv}
    \item \code{sample_sensitivity_results.csv}
    \item \code{mmnl_correlated_results.csv}
  }
}
\examples{
\dontrun{
# Quickest individual driver: ~5 min on a 4-class unbalanced design
res <- klue_study_unbalanced()

# Run a custom slice of the main simulation: 2-class only, 3 reps
res <- klue_study_main(true_K_values = 2, n_reps = 3, C_cands = 1:4)

# Full paper replication (~12 hours)
all_results <- klue_study()
}
}
\seealso{\code{\link{klue}}, \code{\link{klue_simulate}}, \code{\link{klue_dgp}}}
