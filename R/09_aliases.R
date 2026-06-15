# Backward-compatibility aliases. Old names map to the canonical klue_*
# names; both are exported so existing call sites keep working. Slated for
# removal in a future major release -- new code should use the klue_* names.

#' Deprecated backward-compatibility aliases
#'
#' Each function listed here forwards to its canonical \code{klue_} function.
#' These aliases are kept so older scripts that use the previous names keep
#' working.
#'
#' @name klue-deprecated
#' @keywords internal
#' @rdname klue-deprecated
#' @export
run_lcmnl_workflow  <- function(...) klue(...)
#' @rdname klue-deprecated
#' @export
build_database      <- function(...) klue_database(...)
#' @rdname klue-deprecated
#' @export
build_database_long <- function(...) klue_database_long(...)
#' @rdname klue-deprecated
#' @export
build_database_wide <- function(...) klue_database_wide(...)

#' @rdname klue-deprecated
#' @export
make_dgp_config           <- function(...) klue_dgp(...)
#' @rdname klue-deprecated
#' @export
estimate_lcmnl_multistart <- function(...) klue_lcmnl(...)
#' @rdname klue-deprecated
#' @export
estimate_mmnl             <- function(...) klue_mmnl(...)
#' @rdname klue-deprecated
#' @export
estimate_mmnl_corr        <- function(...) klue_mmnl_corr(...)
#' @rdname klue-deprecated
#' @export
estimate_lcmnl_multistart_onehot <- function(database, C, dgp = DGP_DEFAULT)
  klue_lcmnl(database, C, dgp = dgp, feature_type = "onehot")

#' @rdname klue-deprecated
#' @export
generate_data                 <- function(...) klue_simulate(...)
#' @rdname klue-deprecated
#' @export
generate_data_with_covariates <- function(...) klue_simulate_cov(...)
#' @rdname klue-deprecated
#' @export
generate_data_defficient      <- function(...) klue_simulate_deff(...)

#' @rdname klue-deprecated
#' @export
run_full_study                 <- function(...) klue_study(...)
#' @rdname klue-deprecated
#' @export
run_main_simulation            <- function(...) klue_study_main(...)
#' @rdname klue-deprecated
#' @export
run_mmnl_comparison            <- function(...) klue_study_mmnl(...)
#' @rdname klue-deprecated
#' @export
run_convergence_ablation       <- function(...) klue_study_convergence(...)
#' @rdname klue-deprecated
#' @export
run_initialisation_ablation    <- function(...) klue_study_initialisation(...)
#' @rdname klue-deprecated
#' @export
run_unbalanced_analysis        <- function(...) klue_study_unbalanced(...)
#' @rdname klue-deprecated
#' @export
run_design_comparison          <- function(...) klue_study_design(...)
#' @rdname klue-deprecated
#' @export
run_concomitant_analysis       <- function(...) klue_study_concomitant(...)
#' @rdname klue-deprecated
#' @export
run_unconditional_recovery     <- function(...) klue_study_recovery(...)
#' @rdname klue-deprecated
#' @export
run_clustering_comparison      <- function(...) klue_study_clustering(...)
#' @rdname klue-deprecated
#' @export
run_sample_sensitivity         <- function(...) klue_study_sample(...)
#' @rdname klue-deprecated
#' @export
run_correlated_mmnl_robustness <- function(...) klue_study_mmnl_corr(...)
