\name{klue_demo}
\alias{klue_demo}
\title{Run a zero-setup demo of the klue workflow}
\description{
  Estimates the klue workflow on Apollo's bundled Swiss route-choice dataset
  (348 commuters, 9 binary tasks). Useful as a first call to see the output
  shape before wiring your own data.
}
\usage{
klue_demo(full = FALSE, verbose = TRUE)
}
\arguments{
  \item{full}{If \code{FALSE} (default), runs a fast slice: \code{C = 1:2},
    no MMNL benchmark (~3 sec). If \code{TRUE}, runs the full workflow:
    \code{C = 1:6} + MMNL (~30 sec).}
  \item{verbose}{Print progress lines and the final summary. Default \code{TRUE}.}
}
\value{
  Invisibly, the same results list returned by \code{\link{run_lcmnl_workflow}}.
}
\examples{
\dontrun{
library(klue)
klue_demo()              # fast slice
res <- klue_demo(full = TRUE)
res$summary
res$class_betas
res$best_C
}
}
\seealso{\code{\link{run_lcmnl_workflow}}}
