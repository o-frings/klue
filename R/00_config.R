# Constants, DGP configuration, and tuning defaults.

MAX_ITER   <- 500L
OUTPUT_DIR <- "output"

# DGP configuration: problem dimensions, so all functions work with variable
# numbers of attributes and alternatives.

#' Data-generating-process configuration
#'
#' Builds the problem-dimension specification shared by the estimators,
#' simulators, and study drivers, so they work with a variable number of
#' generic attributes and alternatives.
#'
#' @param n_generic number of generic attributes (the price attribute is added
#'   on top of these).
#' @param n_alternatives number of alternatives per choice task; the last
#'   alternative is the reference (no ASC).
#' @return A named list with \code{n_generic}, \code{n_alternatives},
#'   \code{n_beta} (generic + price), \code{n_asc} (alternatives - 1),
#'   \code{npc} (parameters per class), \code{beta_bar} (population-mean
#'   coefficients), \code{attr_names}, and \code{price_idx}.
#' @export
klue_dgp <- function(n_generic = 4, n_alternatives = 3) {
  n_beta <- n_generic + 1L         # generic attributes + price
  n_asc  <- n_alternatives - 1L    # ASCs (last alt = reference)
  list(
    n_generic      = n_generic,
    n_alternatives = n_alternatives,
    n_beta         = n_beta,
    n_asc          = n_asc,
    npc            = n_beta + n_asc,
    beta_bar       = c(rep(0.5, n_generic), -1.5),
    attr_names     = c(paste0("x", 1:n_generic), "price"),
    price_idx      = n_beta
  )
}

DGP_DEFAULT <- klue_dgp(4, 3)
BETA_BAR    <- DGP_DEFAULT$beta_bar   # backward-compatible global

# Cores for condition-level parallelism in the study drivers (not Apollo).
# Default leaves two cores free; the generous upper cap only guards against
# fork/memory blow-up on very-high-core servers. Benchmarks show throughput
# keeps rising to ~physical-core count, so the old hard cap of 4 left ~1.5x on
# the table for 8-12 core machines. Override with options(klue.cores = n) --
# e.g. set it to 1 when wrapping a driver in your own parallel loop.
.klue_cores <- function() {
  getOption("klue.cores",
            min(16L, max(1L, parallel::detectCores() - 2L)))
}

# ---- MMNL defaults ----------------------------------------------------------
# Overridable per-call and via options("klue.mmnl.<name>").
N_DRAWS_MMNL            <- 3000L
N_DRAWS_MMNL_STAGE1     <- 200L
DRAWS_TYPE_MMNL         <- "mlhs"
ESTIMATION_ROUTINE_MMNL <- "bgw"
# Box constraints on the log-normal price parameters: without bounds the
# optimiser can drift mu_price into a region where b_price is enormous and
# the likelihood becomes non-finite. Pass NULL to disable.
MU_PRICE_BOUNDS_MMNL    <- c(-5, 3)
SIGMA_PRICE_BOUNDS_MMNL <- c(-3, 1)

.klue_default_mmnl_cores <- function() {
  nc <- tryCatch(parallel::detectCores(logical = FALSE),
                 error = function(e) NA_integer_)
  if (is.na(nc) || !is.finite(nc) || nc < 1L) nc <- 1L
  max(1L, as.integer(nc) - 1L)
}

#' MMNL default settings
#'
#' Returns the active MMNL defaults as a named list. Values can be overridden
#' per-call (arguments to \code{klue_mmnl} / \code{klue}) or globally via
#' \code{options()} entries \code{klue.mmnl.<name>}.
#'
#' @return A named list of the currently active defaults.
#' @export
klue_mmnl_defaults <- function() {
  list(
    n_draws            = getOption("klue.mmnl.n_draws",            N_DRAWS_MMNL),
    n_draws_stage1     = getOption("klue.mmnl.n_draws_stage1",     N_DRAWS_MMNL_STAGE1),
    draws_type         = getOption("klue.mmnl.draws_type",         DRAWS_TYPE_MMNL),
    estimation_routine = getOption("klue.mmnl.estimation_routine", ESTIMATION_ROUTINE_MMNL),
    n_cores            = getOption("klue.mmnl.n_cores",            .klue_default_mmnl_cores()),
    mu_price_bounds    = getOption("klue.mmnl.mu_price_bounds",    MU_PRICE_BOUNDS_MMNL),
    sigma_price_bounds = getOption("klue.mmnl.sigma_price_bounds", SIGMA_PRICE_BOUNDS_MMNL),
    quiet              = getOption("klue.mmnl.quiet",              TRUE)
  )
}
