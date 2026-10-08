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
#' @param asc_map optional integer vector of length \code{n_alternatives}
#'   mapping each alternative to an ASC group (0 = reference / no ASC; group
#'   ids must be contiguous 1..K). Alternatives sharing a group share one ASC,
#'   e.g. \code{c(1, 1, 2, 0)} gives two offering tiles a common ASC, the third
#'   alternative its own, and the fourth as reference. Default \code{NULL}
#'   reproduces the per-alternative structure \code{c(1, .., n_alt-1, 0)}.
#' @param labels optional names of the \code{n_generic + 1} taste coefficients
#'   (price last). Fits then name their parameters \code{<label>} (MNL) or
#'   \code{<label>_class<c>} (latent class). Default \code{NULL} keeps
#'   \code{b1..b{n_beta}}.
#' @param asc_labels optional names of the \code{n_asc} constants. Default
#'   \code{NULL} keeps \code{asc1..asc{n_asc}}.
#' @param common optional coefficients shared by all latent classes, as labels
#'   or as positions in the per-class vector (taste coefficients, then
#'   constants). A latent class fit then estimates one value for each, and its
#'   \code{k} falls by \code{C - 1} per common coefficient. Default none.
#' @return A named list with \code{n_generic}, \code{n_alternatives},
#'   \code{n_beta} (generic + price), \code{n_asc} (number of ASC groups),
#'   \code{npc} (parameters per class), \code{asc_map}, \code{beta_bar}
#'   (population-mean coefficients), \code{attr_names}, \code{price_idx},
#'   \code{par_labels} (the per-class parameter names) and \code{common}
#'   (positions of the common coefficients).
#' @export
klue_dgp <- function(n_generic = 4, n_alternatives = 3, asc_map = NULL,
                     labels = NULL, asc_labels = NULL, common = NULL) {
  n_beta <- n_generic + 1L         # generic attributes + price
  if (is.null(asc_map)) asc_map <- c(seq_len(n_alternatives - 1L), 0L)
  asc_map <- as.integer(asc_map)
  if (length(asc_map) != n_alternatives)
    stop("asc_map must have length n_alternatives (", n_alternatives, ")")
  pos <- asc_map[asc_map > 0L]
  if (any(asc_map < 0L) || (length(pos) && !setequal(pos, seq_len(max(pos)))))
    stop("asc_map entries must be 0 (reference) or contiguous group ids 1..K")
  n_asc <- if (length(pos)) max(pos) else 0L   # number of ASC groups
  if (!is.null(labels) && length(labels) != n_beta)
    stop("labels must name the ", n_beta, " taste coefficients (price last)")
  if (!is.null(asc_labels) && length(asc_labels) != n_asc)
    stop("asc_labels must name the ", n_asc, " constants")
  par_labels <- c(if (is.null(labels)) .lab("b", n_beta) else as.character(labels),
                  if (is.null(asc_labels)) .lab("asc", n_asc) else as.character(asc_labels))
  if (anyDuplicated(par_labels)) stop("parameter labels must be unique")
  dgp <- list(
    n_generic      = n_generic,
    n_alternatives = n_alternatives,
    n_beta         = n_beta,
    n_asc          = n_asc,
    npc            = n_beta + n_asc,
    asc_map        = asc_map,
    beta_bar       = c(rep(0.5, n_generic), -1.5),
    attr_names     = c(paste0("x", 1:n_generic), "price"),
    price_idx      = n_beta,
    par_labels     = par_labels,
    common         = integer(0)
  )
  .with_common(dgp, common)
}

# Numbered names prefix1..prefixN; none when n = 0 (paste0 would return "prefix").
.lab <- function(prefix, n) if (n > 0L) paste0(prefix, seq_len(n)) else character(0)

# Per-class parameter labels of a dgp, including dgps built before 0.10.0.
.par_labels <- function(dgp)
  if (!is.null(dgp$par_labels)) dgp$par_labels else c(.lab("b", dgp$n_beta), .lab("asc", dgp$n_asc))

# The dgp with `common` (labels or per-class positions) as its common coefficients.
.with_common <- function(dgp, common) {
  if (is.null(common) || !length(common)) { dgp$common <- integer(0); return(dgp) }
  l <- .par_labels(dgp)
  idx <- if (is.character(common)) match(common, l) else as.integer(common)
  bad <- is.na(idx) | idx < 1L | idx > length(l)
  if (any(bad)) stop("unknown common coefficient(s): ", paste(common[bad], collapse = ", "))
  dgp$common <- sort(unique(idx)); dgp
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
# Overridable per-call and via options("klue.mmnl.<name>"). The estimation
# routine and the core count are Apollo's own defaults (BGW, 1 core); 3000 MLHS
# draws is the paper's protocol (MLHS rather than Halton, which Apollo's manual
# advises against beyond five random coefficients).
N_DRAWS_MMNL            <- 3000L
DRAWS_TYPE_MMNL         <- "mlhs"
ESTIMATION_ROUTINE_MMNL <- "bgw"
N_CORES_MMNL            <- 1L
# Apollo's apollo_control$memorySaver (default FALSE): TRUE computes the
# analytic gradient in chunks of about two respondents, which cuts its memory
# (otherwise about one observations-by-draws array per parameter) without
# changing the estimates. Also settable with the environment variable
# KLUE_MMNL_MEMORY_SAVER=TRUE.
MEMORY_SAVER_MMNL       <- FALSE

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
    draws_type         = getOption("klue.mmnl.draws_type",         DRAWS_TYPE_MMNL),
    estimation_routine = getOption("klue.mmnl.estimation_routine", ESTIMATION_ROUTINE_MMNL),
    n_cores            = getOption("klue.mmnl.n_cores",            N_CORES_MMNL),
    memory_saver       = getOption("klue.mmnl.memory_saver",
                                   isTRUE(as.logical(Sys.getenv("KLUE_MMNL_MEMORY_SAVER",
                                                                MEMORY_SAVER_MMNL)))),
    quiet              = getOption("klue.mmnl.quiet",              TRUE)
  )
}
