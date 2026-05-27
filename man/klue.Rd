\name{klue}
\alias{klue}
\alias{run_lcmnl_workflow}
\title{Run the full LCMNL specification workflow on a discrete-choice dataset}
\description{
  Estimates LCMNL for a range of class counts, runs an MMNL benchmark
  (independent and optionally correlated normals), produces a
  BIC / AIC / ICL / entropy summary table and a class-specific coefficient
  table, and writes CSVs to disk.

  \code{run_lcmnl_workflow} is a deprecated alias for the same function.
}
\usage{
klue(database = NULL,
     data = NULL,
     format = c("auto", "long", "wide"),
     C_cands = 1:6,
     run_mmnl = TRUE,
     run_mmnl_corr = FALSE,
     mmnl_opts = list(),
     attr_labels = NULL,
     output_prefix = "workflow",
     output_dir = NULL,
     write_csv = TRUE,
     verbose = TRUE,
     \ldots)
}
\arguments{
  \item{database}{An already-built canonical wide data frame (columns
    \code{ID}, \code{TASK}, \code{x1_1..xN_J}, \code{price_1..price_J},
    \code{CHOICE}). If supplied, \code{data} and the column-mapping
    arguments are ignored.}
  \item{data}{A data frame or CSV path in long or wide format.
    \code{\link{klue_database}} is invoked internally; pass the
    appropriate column-mapping arguments via \code{\ldots}.}
  \item{format}{Forwarded to \code{\link{klue_database}}.}
  \item{C_cands}{Integer vector of class counts to estimate. Default
    \code{1:6}. \code{C = 1} is the MNL baseline.}
  \item{run_mmnl}{If \code{TRUE}, also estimates an MMNL benchmark with
    independent normal random coefficients (log-normal on price).}
  \item{run_mmnl_corr}{If \code{TRUE}, additionally estimates a correlated
    MMNL benchmark with a full Cholesky-parameterised covariance matrix
    over random coefficients (log-normal on price).}
  \item{mmnl_opts}{Named list of arguments forwarded to BOTH
    \code{\link{klue_mmnl}} and \code{\link{klue_mmnl_corr}}. See
    \code{\link{klue_mmnl_defaults}} for the active defaults.}
  \item{attr_labels}{Override display column names in the betas CSV.}
  \item{output_prefix}{Filename prefix for the written CSVs.}
  \item{output_dir}{Output directory. Default: \code{getOption("klue.output_dir", "output")}.}
  \item{write_csv}{Whether to write CSV files (\code{TRUE} by default).}
  \item{verbose}{Per-C progress lines and a final printout.}
  \item{\ldots}{Column-mapping arguments forwarded to \code{\link{klue_database}}.}
}
\value{
  Invisibly, a list:
  \item{database}{Canonical wide data frame used for estimation.}
  \item{dgp}{Output of \code{\link{klue_dgp}}.}
  \item{lcmnl}{Named list of per-C fits.}
  \item{mmnl}{Independent-normals MMNL fit (\code{NULL} if disabled or failed).}
  \item{mmnl_corr}{Correlated MMNL fit (\code{NULL} if \code{run_mmnl_corr = FALSE}
    or it failed).}
  \item{summary}{Data frame with one row per \code{C}: \code{LL}, \code{k},
    \code{BIC}, \code{AIC}, \code{ICL}, \code{ICL_BIC}, \code{dBIC},
    \code{dAIC}, \code{dICL}, \code{best_method}.}
  \item{class_betas}{Per-class coefficients for the BIC-best model.}
  \item{comparison}{MNL vs LCMNL-best vs MMNL (independent) vs MMNL
    (correlated), with whichever subset was estimated.}
  \item{best_C}{BIC-best number of classes.}
  \item{best_lcmnl}{The BIC-best fit.}
}
\section{Files written}{
  When \code{write_csv = TRUE}:
  \itemize{
    \item \code{<prefix>_summary.csv}
    \item \code{<prefix>_class_betas.csv}
    \item \code{<prefix>_model_comparison.csv} (if \code{run_mmnl = TRUE}
      or \code{run_mmnl_corr = TRUE})
  }
}
\examples{
\dontrun{
# Long-format CSV
res <- klue(
  data = "data.csv", format = "long",
  id_col = "id", task_col = "task",
  alt_col = "alt", choice_col = "chosen",
  attribute_cols = c("a", "b", "c"), price_col = "price",
  C_cands = 1:6
)
res$summary
res$best_C
}
}
\seealso{\code{\link{klue_database}}, \code{\link{klue_demo}}}
