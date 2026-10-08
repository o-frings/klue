#' klue: Hybrid ML and Random-Utility Workflow for LCMNL Specification
#'
#' Implements the three-step workflow of Frings (2026) for specifying latent
#' class multinomial logit models: a maximum-likelihood multi-start initialised
#' from six clusterings of respondents' revealed-preference signatures; LCMNL
#' estimation across a range of class counts; and a mixed multinomial logit
#' benchmark, reported with BIC, AIC, ICL, and a classification-entropy
#' diagnostic. The 0.9 series is a rewrite of the 0.6.x engine with the same
#' exported API. The MMNL benchmark requires the suggested 'apollo' package;
#' everything else runs without it.
#'
#' Start with \code{\link{klue_demo}} for a zero-setup example, then
#' \code{\link{klue}} on your own long- or wide-format choice data. The
#' Monte Carlo study of the paper lives outside the package in
#' \code{studies/klue_studies.R} (see README_REPRODUCE.md).
#'
#' @docType package
#' @name klue-package
#' @aliases klue-package
#' @keywords internal
#' @importFrom mclust Mclust mclustBIC
#' @importFrom cluster pam
#' @importFrom stats kmeans hclust dist as.dist cutree optim rnorm runif qnorm
#'   pnorm sd median model.matrix ave aggregate na.omit setNames
#' @importFrom utils read.csv read.table write.csv head tail
#' @importFrom parallel mclapply detectCores
NULL
