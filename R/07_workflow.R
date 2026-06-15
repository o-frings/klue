# klue(): the "plug in data -> all results" entry point, and klue_demo().

#' Run the full LCMNL specification workflow
#'
#' Three calling conventions: (1) `database =` canonical wide format;
#' (2) `data =` + long-format column mapping; (3) `data =` + wide-format
#' column mapping (both forwarded to \code{klue_database}). Loops over
#' `C_cands`, optionally estimates the MMNL benchmarks, prints and writes
#' summary CSVs, and returns the full results list invisibly.
#'
#' @param database canonical wide-format data frame with columns ID, TASK,
#'   CHOICE plus attribute columns. Supply this or `data`, not both.
#' @param data raw data frame passed to \code{klue_database} together with the
#'   column-mapping arguments in `...` when `database` is not supplied.
#' @param format one of "auto", "long", "wide": the layout of `data`, forwarded
#'   to \code{klue_database}. Default "auto".
#' @param C_cands integer vector of class counts to estimate. Default 1:6.
#' @param run_mmnl logical; if TRUE estimate the independent-normals MMNL
#'   benchmark. Default TRUE.
#' @param run_mmnl_corr logical; if TRUE estimate the correlated-normals MMNL
#'   benchmark. Default FALSE.
#' @param mmnl_opts named list of extra arguments forwarded to \code{klue_mmnl}
#'   for both MMNL flavours. Default an empty list.
#' @param attr_labels character vector of length `n_generic + 1` naming the
#'   attribute coefficients (the last entry is price). If NULL, taken from the
#'   database attribute or generated as x1, x2, ..., price.
#' @param output_prefix character prefix for the written CSV file names.
#'   Default "workflow".
#' @param output_dir directory for the written CSVs. If NULL, uses OUTPUT_DIR.
#' @param write_csv logical; if TRUE write the summary, class-coefficient and
#'   model-comparison CSVs. Default TRUE.
#' @param verbose logical; if TRUE print progress and result tables. Default
#'   TRUE.
#' @param n_cores cores for the per-start LCMNL fits within each C (forwarded
#'   to \code{klue_lcmnl}). Default 1 (sequential). Setting it to the number of
#'   physical cores fits the six clustering starts concurrently and is the main
#'   speed lever for a single dataset (~3x on a typical multicore machine).
#' @param ... column-mapping and scaling arguments forwarded to
#'   \code{klue_database} when `data` is supplied instead of `database`.
#' @return Invisibly, a list with elements `database`, `dgp`, `lcmnl` (the
#'   per-C fits), `mmnl`, `mmnl_corr`, `summary`, `class_betas`, `comparison`,
#'   `best_C` (the BIC-best class count) and `best_lcmnl` (the BIC-best fit).
#' @export
klue <- function(database = NULL,
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
                 n_cores = 1L,
                 ...) {
  format <- match.arg(format)

  if (is.null(database)) {
    if (is.null(data)) {
      stop("Provide either `database` (canonical wide format) or `data` plus ",
           "the column-mapping arguments (klue_database is called internally).")
    }
    database <- klue_database(data, format = format, verbose = verbose, ...)
  } else if (!all(c("ID", "TASK", "CHOICE") %in% names(database))) {
    stop("`database` must contain columns ID, TASK, CHOICE (canonical wide format).")
  }

  J  <- attr(database, "n_alternatives")
  Ng <- attr(database, "n_generic")
  if (is.null(J) || is.null(Ng)) {
    J  <- length(grep("^price_\\d+$", names(database)))
    Ng <- length(grep("^x\\d+_1$", names(database)))
    if (J < 2 || Ng < 1) {
      stop("Could not infer n_alternatives / n_generic from database columns.")
    }
  }
  dgp <- klue_dgp(n_generic = Ng, n_alternatives = J)

  if (is.null(attr_labels)) {
    attr_labels <- attr(database, "attr_labels")
    if (is.null(attr_labels)) attr_labels <- c(paste0("x", seq_len(Ng)), "price")
  }
  stopifnot(length(attr_labels) == Ng + 1)

  if (is.null(output_dir)) output_dir <- OUTPUT_DIR
  if (write_csv && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- LCMNL across C_cands -------------------------------------------------
  results <- list()
  for (cc in C_cands) {
    if (verbose) cat(sprintf("\n=== Estimating C = %d ===\n", cc))
    t0 <- Sys.time()
    m <- klue_lcmnl(database, cc, dgp = dgp, n_cores = n_cores)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (verbose) {
      cat(sprintf("  converged=%s  LL=%.2f  BIC=%.2f  AIC=%.2f  ICL=%.2f  ICL-BIC=%.2f  method=%s  time=%.1fs\n",
                  m$converged, m$LL, m$BIC, m$AIC, m$ICL,
                  ifelse(is.na(m$ICL_BIC), 0, m$ICL_BIC),
                  ifelse(is.null(m$best_method) || is.na(m$best_method),
                         "MNL", m$best_method), dt))
    }
    results[[as.character(cc)]] <- m
  }

  summary_df <- do.call(rbind, lapply(names(results), function(nm) {
    m <- results[[nm]]
    data.frame(
      C = as.integer(nm), converged = m$converged,
      LL = m$LL, k = m$k, BIC = m$BIC, AIC = m$AIC, ICL = m$ICL,
      ICL_BIC = ifelse(is.null(m$ICL_BIC), NA_real_, m$ICL_BIC),
      best_method = ifelse(is.null(m$best_method) || is.na(m$best_method),
                           "MNL", m$best_method),
      stringsAsFactors = FALSE
    )
  }))
  summary_df$dBIC <- summary_df$BIC - min(summary_df$BIC, na.rm = TRUE)
  summary_df$dAIC <- summary_df$AIC - min(summary_df$AIC, na.rm = TRUE)
  summary_df$dICL <- summary_df$ICL - min(summary_df$ICL, na.rm = TRUE)

  best_C <- summary_df$C[which.min(summary_df$BIC)]
  best   <- results[[as.character(best_C)]]

  betas_df <- as.data.frame(best$betas)
  names(betas_df) <- attr_labels
  betas_df$class  <- seq_len(nrow(betas_df))
  betas_df$share  <- best$class_probs
  betas_df <- betas_df[, c("class", "share", attr_labels)]

  # ---- MMNL benchmarks ------------------------------------------------------
  if (!is.list(mmnl_opts)) stop("`mmnl_opts` must be a named list.")
  .fit_mmnl_flavour <- function(label, correlation) {
    if (verbose) cat(sprintf("\n=== Estimating MMNL (%s) ===\n", label))
    t0 <- Sys.time()
    args <- c(list(database = database, correlation = correlation, dgp = dgp),
              mmnl_opts)
    fit <- tryCatch(do.call(klue_mmnl, args),
                    error = function(e) {
                      if (verbose) cat("  MMNL (", label, ") failed: ",
                                       conditionMessage(e), "\n", sep = "")
                      NULL
                    })
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(fit) && verbose) {
      cat(sprintf("  converged=%s  LL=%.2f  BIC=%.2f  AIC=%.2f  k=%d  time=%.1fs\n",
                  fit$converged, fit$LL, fit$BIC, fit$AIC, fit$k, dt))
      if (!isTRUE(fit$converged)) {
        cat(sprintf("  MMNL (%s) did not converge (reason: %s)\n", label,
                    if (is.null(fit$reason)) "unknown" else fit$reason))
        if (!is.null(fit$apollo_log_tail)) {
          cat("  Last lines of the Apollo log:\n")
          cat(paste0("    ", fit$apollo_log_tail), sep = "\n")
          cat("\n")
        }
        if (!is.null(fit$apollo_log_path)) {
          cat(sprintf("  Full Apollo log: %s\n", fit$apollo_log_path))
        }
      }
    }
    fit
  }

  mmnl_fit      <- if (run_mmnl) .fit_mmnl_flavour("independent normals", FALSE) else NULL
  mmnl_corr_fit <- if (run_mmnl_corr) .fit_mmnl_flavour("correlated normals", TRUE) else NULL

  comparison_df <- NULL
  if ((!is.null(mmnl_fit) || !is.null(mmnl_corr_fit)) && "1" %in% names(results)) {
    mnl <- results[["1"]]
    row_of <- function(model, m) data.frame(model = model, LL = m$LL, k = m$k,
                                            BIC = m$BIC, AIC = m$AIC,
                                            stringsAsFactors = FALSE)
    rows <- list(row_of("MNL (C=1)", mnl),
                 row_of(sprintf("LCMNL (C=%d, BIC-best)", best_C), best))
    if (!is.null(mmnl_fit))
      rows <- c(rows, list(row_of("MMNL (independent normals)", mmnl_fit)))
    if (!is.null(mmnl_corr_fit))
      rows <- c(rows, list(row_of("MMNL (correlated normals)", mmnl_corr_fit)))
    comparison_df <- do.call(rbind, rows)
    comparison_df$dBIC <- comparison_df$BIC - min(comparison_df$BIC)
  }

  # ---- Output ---------------------------------------------------------------
  if (write_csv) {
    out <- function(suffix) file.path(output_dir, paste0(output_prefix, suffix))
    write.csv(summary_df, out("_summary.csv"), row.names = FALSE)
    write.csv(betas_df, out("_class_betas.csv"), row.names = FALSE)
    if (!is.null(comparison_df)) {
      write.csv(comparison_df, out("_model_comparison.csv"), row.names = FALSE)
    }
  }

  if (verbose) {
    cat("\n==== LCMNL SUMMARY ====\n")
    print(summary_df, row.names = FALSE)
    cat(sprintf("\nBIC-best K* = %d\n", best_C))
    cat("\n==== CLASS COEFFICIENTS (BIC-best) ====\n")
    print(betas_df, row.names = FALSE)
    if (!is.null(comparison_df)) {
      cat("\n==== MODEL COMPARISON ====\n")
      print(comparison_df, row.names = FALSE)
      cat(sprintf("\nBIC-preferred model: %s\n",
                  comparison_df$model[which.min(comparison_df$BIC)]))
    }
    if (write_csv) {
      cat(sprintf("\nWrote %s_{summary,class_betas%s}.csv to %s/\n",
                  output_prefix,
                  if (!is.null(comparison_df)) ",model_comparison" else "",
                  output_dir))
    }
    if (!isTRUE(getOption("klue.suppress_citation", FALSE))) {
      pkg_ver <- tryCatch(as.character(utils::packageVersion("klue")),
                          error = function(e) "dev")
      cat("\n--- Please cite klue if you use it in published work ---\n")
      cat("  Frings (2026). A Hybrid Machine Learning and Random Utility\n")
      cat("    Framework for Latent Class Model Specification. Working paper.\n")
      cat(sprintf("  Frings (2026). klue: R package version %s.\n", pkg_ver))
      cat("    https://github.com/o-frings/klue\n")
      cat("\n  Plus the upstream packages: apollo (Hess & Palma 2019),\n")
      cat("  mclust (Scrucca et al. 2016), cluster (Maechler et al.).\n")
      cat("  citation(\"klue\")  # full BibTeX bundle\n")
    }
  }

  invisible(list(
    database    = database,
    dgp         = dgp,
    lcmnl       = results,
    mmnl        = mmnl_fit,
    mmnl_corr   = mmnl_corr_fit,
    summary     = summary_df,
    class_betas = betas_df,
    comparison  = comparison_df,
    best_C      = best_C,
    best_lcmnl  = best
  ))
}

#' Zero-setup demonstration on Apollo's Swiss route-choice data
#'
#' `full = FALSE` (default): C = 1..2, no MMNL (~3 sec). `full = TRUE`:
#' C = 1..6 + MMNL benchmark (~30 sec).
#'
#' @param full logical; if TRUE estimate C = 1..6 with the MMNL benchmark,
#'   otherwise C = 1..2 without MMNL. Default FALSE.
#' @param verbose logical; if TRUE print progress and result tables. Default
#'   TRUE.
#' @return Invisibly, the list returned by \code{klue} for the Swiss
#'   route-choice data.
#' @export
klue_demo <- function(full = FALSE, verbose = TRUE) {
  if (!requireNamespace("apollo", quietly = TRUE)) {
    stop("klue_demo() needs the apollo package installed.")
  }
  data_env <- new.env()
  utils::data("apollo_swissRouteChoiceData", package = "apollo", envir = data_env)
  d <- data_env$apollo_swissRouteChoiceData

  if (verbose) {
    cat("klue_demo: Swiss route choice (Apollo example data).\n")
    cat("  348 commuters x 9 binary route comparisons (after balanced-panel filter).\n")
    cat("  3 attributes (travel time, headway, # changes) + travel cost.\n")
    cat(sprintf("  Estimating LCMNL for C = 1..%d%s.\n\n",
                if (full) 6 else 2,
                if (full) " + MMNL benchmark" else ""))
  }

  klue(
    data           = d,
    format         = "wide",
    id_col         = "ID", task_col = NULL, choice_col = "choice",
    attribute_cols = list(
      travel_time = c("tt1", "tt2"),
      headway     = c("hw1", "hw2"),
      changes     = c("ch1", "ch2")
    ),
    price_col      = c("tc1", "tc2"),
    scalings       = list(travel_time = 60, headway = 60, price = 10),
    C_cands        = if (full) 1:6 else 1:2,
    run_mmnl       = full,
    write_csv      = FALSE,
    verbose        = verbose,
    output_prefix  = "demo"
  )
}
