# =============================================================================
# LCMNL_WORKFLOW.R
#
# Generic "plug in data -> all model outputs" entry point for the hybrid
# LCMNL specification workflow developed in `R/simulation_study.R`. Wraps the
# engine in two user-facing helpers:
#
#   build_database(...)         Long- or wide-format -> canonical wide format.
#   run_lcmnl_workflow(...)     Loops C, runs MMNL, summarises, writes CSVs.
#
# The engine (make_dgp_config, estimate_lcmnl_multistart, estimate_mmnl) lives
# in simulation_study.R and is reused unchanged.
#
# Canonical wide database format (what the engine consumes):
#   ID, TASK, x1_1..xN_J, price_1..price_J, CHOICE
# where N = n_generic (number of generic attributes), J = n_alternatives, and
# CHOICE is the integer index of the chosen alternative in 1..J. The reference
# alternative is the last one (no ASC).
#
# Availability note: the engine assumes every alternative is available in every
# task. Both build_database helpers accept availability columns and filter the
# panel to fully-available tasks (dropping tasks where any alt is unavailable).
# The drop rate is reported. Partial-availability estimation is not supported.
# =============================================================================



# =============================================================================
# Internal helpers
# =============================================================================

.coerce_numeric_zero_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

.filter_fully_available_long <- function(raw, task_idx_col, alt_col, avail_col,
                                         J, verbose = TRUE) {
  if (is.null(avail_col)) return(raw)
  if (!avail_col %in% names(raw)) {
    stop("Availability column `", avail_col, "` not found in data.")
  }
  avail <- suppressWarnings(as.numeric(raw[[avail_col]]))
  avail[is.na(avail)] <- 0
  raw[[avail_col]] <- avail

  total_avail_per_task <- tapply(raw[[avail_col]], raw[[task_idx_col]], sum)
  full_tasks <- as.integer(names(total_avail_per_task)[total_avail_per_task == J])
  n_total  <- length(unique(raw[[task_idx_col]]))
  n_full   <- length(full_tasks)
  if (verbose && n_full < n_total) {
    cat(sprintf("Availability filter: dropping %d of %d tasks (%.1f%%) where not all %d alternatives are available.\n",
                n_total - n_full, n_total, 100 * (n_total - n_full) / n_total, J))
  }
  if (n_full == 0) stop("No tasks remain after availability filter.")
  raw[raw[[task_idx_col]] %in% full_tasks, , drop = FALSE]
}

.filter_balanced_panel <- function(db, verbose = TRUE) {
  tasks_per_id <- table(db$ID)
  T_max <- max(tasks_per_id)
  full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == T_max])
  if (length(full_ids) < length(tasks_per_id)) {
    n_drop <- length(tasks_per_id) - length(full_ids)
    if (verbose) {
      cat(sprintf("Balanced-panel filter: dropping %d of %d respondents (%.1f%%) with fewer than %d tasks.\n",
                  n_drop, length(tasks_per_id),
                  100 * n_drop / length(tasks_per_id), T_max))
    }
    db <- db[db$ID %in% full_ids, , drop = FALSE]
    db$ID <- as.integer(factor(db$ID))
  }
  db
}


# =============================================================================
# build_database_long
# =============================================================================
#
# Long-format input: one row per (respondent, task, alternative).
#
# Args
#   data            data.frame or path to CSV.
#   id_col, task_col, alt_col, choice_col   column names.
#   attribute_cols  character vector of generic attribute column names. Order
#                   is preserved (becomes x1..xN). All must be present in data.
#   price_col       name of the price/cost column. Becomes the last attribute.
#   choice_format   "indicator" (default): choice_col is 0/1, exactly one
#                       chosen row per (id, task).
#                   "alt_index": choice_col is the chosen alt index, repeated
#                       on every row of the (id, task) group.
#   avail_col       optional column with 0/1 availability per row. Tasks where
#                   any alternative has avail==0 are dropped, with a warning.
#                   Default NULL (assume everything available).
#   scalings        optional named list of per-attribute scalings. Names match
#                   entries in `attribute_cols` and/or `price_col`. Each value
#                   is a scalar that the corresponding column is divided by.
#                   Default: NULL (no scaling). Example: list(time = 60, cost = 10).
#                   Backward compat: `price_scaling` is honoured as a shortcut
#                   for `scalings = list(<price_col> = price_scaling)`.
#   price_scaling   shortcut for scalings[[price_col]]. Default 1.
#   verbose         print progress and filter rates.
# =============================================================================

build_database_long <- function(data,
                                id_col, task_col, alt_col, choice_col,
                                attribute_cols, price_col,
                                choice_format = c("indicator", "alt_index"),
                                avail_col = NULL,
                                scalings = NULL,
                                price_scaling = 1,
                                verbose = TRUE) {
  choice_format <- match.arg(choice_format)
  stopifnot(length(attribute_cols) >= 1, length(price_col) == 1)

  raw <- if (is.character(data) && length(data) == 1) {
    read.csv(data, stringsAsFactors = FALSE)
  } else {
    as.data.frame(data)
  }

  required <- c(id_col, task_col, alt_col, choice_col,
                attribute_cols, price_col,
                if (!is.null(avail_col)) avail_col else character(0))
  miss <- setdiff(required, names(raw))
  if (length(miss) > 0) {
    stop("Columns missing from data: ", paste(miss, collapse = ", "))
  }

  if (is.null(scalings)) scalings <- list()
  if (price_scaling != 1 && is.null(scalings[[price_col]])) {
    scalings[[price_col]] <- price_scaling
  }

  for (col in c(attribute_cols, price_col)) {
    raw[[col]] <- .coerce_numeric_zero_na(raw[[col]])
    if (!is.null(scalings[[col]])) {
      stopifnot(is.numeric(scalings[[col]]), scalings[[col]] > 0)
      raw[[col]] <- raw[[col]] / scalings[[col]]
    }
  }

  raw$.resp <- as.integer(factor(raw[[id_col]]))
  raw$.alt  <- as.integer(as.factor(raw[[alt_col]]))
  raw$.task_idx <- as.integer(factor(paste(raw$.resp, raw[[task_col]], sep = "_")))
  raw <- raw[order(raw$.task_idx, raw$.alt), ]

  J <- max(raw$.alt)

  raw <- .filter_fully_available_long(raw,
                                      task_idx_col = ".task_idx",
                                      alt_col      = ".alt",
                                      avail_col    = avail_col,
                                      J            = J,
                                      verbose      = verbose)

  task_sizes <- table(raw$.task_idx)
  bad_tasks  <- as.integer(names(task_sizes)[task_sizes != J])
  if (length(bad_tasks) > 0) {
    if (verbose) {
      cat(sprintf("Dropping %d malformed tasks (not exactly %d alternatives present).\n",
                  length(bad_tasks), J))
    }
    raw <- raw[!raw$.task_idx %in% bad_tasks, , drop = FALSE]
  }
  if (nrow(raw) == 0) stop("No tasks remain after filtering.")

  raw$.task_idx <- as.integer(factor(raw$.task_idx))
  raw$.resp     <- as.integer(factor(raw$.resp))

  tasks_per_id <- tapply(raw$.task_idx, raw$.resp, function(x) length(unique(x)))
  if (length(unique(tasks_per_id)) != 1) {
    db_tmp <- data.frame(ID = raw$.resp,
                         task_count = tasks_per_id[as.character(raw$.resp)])
    raw <- .filter_balanced_panel_long(raw, verbose = verbose)
    raw$.task_idx <- as.integer(factor(raw$.task_idx))
    raw$.resp     <- as.integer(factor(raw$.resp))
  }

  N <- length(unique(raw$.resp))
  T_const <- as.integer(length(unique(raw$.task_idx)) / N)

  db <- data.frame(ID = rep(seq_len(N), each = T_const),
                   TASK = rep(seq_len(T_const), times = N))

  for (j in seq_len(J)) {
    rows_j <- raw[raw$.alt == j, , drop = FALSE]
    rows_j <- rows_j[order(rows_j$.task_idx), ]
    if (nrow(rows_j) != N * T_const) {
      stop(sprintf("Alternative %d appears %d times after filtering, expected %d.",
                   j, nrow(rows_j), N * T_const))
    }
    for (a in seq_along(attribute_cols)) {
      db[[sprintf("x%d_%d", a, j)]] <- rows_j[[attribute_cols[a]]]
    }
    db[[sprintf("price_%d", j)]] <- rows_j[[price_col]]
  }

  if (choice_format == "indicator") {
    chosen <- raw[raw[[choice_col]] == 1, , drop = FALSE]
    if (nrow(chosen) != N * T_const) {
      stop(sprintf("Choice indicator (`%s`): expected exactly one chosen row per (id, task), got %d rows out of %d task groups.",
                   choice_col, nrow(chosen), N * T_const))
    }
    chosen <- chosen[order(chosen$.task_idx), ]
    db$CHOICE <- chosen$.alt
  } else {
    task_choice <- raw[!duplicated(raw$.task_idx), ]
    task_choice <- task_choice[order(task_choice$.task_idx), ]
    db$CHOICE <- as.integer(task_choice[[choice_col]])
    if (any(is.na(db$CHOICE)) || any(db$CHOICE < 1) || any(db$CHOICE > J)) {
      stop(sprintf("Choice values must be integers in 1..%d.", J))
    }
  }

  scaled_labels <- function(col) {
    sc <- scalings[[col]]
    if (is.null(sc) || sc == 1) col else sprintf("%s_per%g", col, sc)
  }
  attr(db, "attr_labels")    <- c(vapply(attribute_cols, scaled_labels, character(1)),
                                  scaled_labels(price_col))
  attr(db, "n_alternatives") <- J
  attr(db, "n_generic")      <- length(attribute_cols)

  if (verbose) {
    cat(sprintf("Built canonical database: N=%d respondents, T=%d tasks, J=%d alternatives, %d generic attributes + price.\n",
                N, T_const, J, length(attribute_cols)))
  }
  db
}

# Variant of .filter_balanced_panel that operates on the long-format frame
.filter_balanced_panel_long <- function(raw, verbose = TRUE) {
  tasks_per_id <- tapply(raw$.task_idx, raw$.resp, function(x) length(unique(x)))
  T_max <- max(tasks_per_id)
  full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == T_max])
  n_drop <- length(tasks_per_id) - length(full_ids)
  if (verbose && n_drop > 0) {
    cat(sprintf("Balanced-panel filter: dropping %d of %d respondents (%.1f%%) with fewer than %d tasks (after availability filter).\n",
                n_drop, length(tasks_per_id),
                100 * n_drop / length(tasks_per_id), T_max))
  }
  raw[raw$.resp %in% full_ids, , drop = FALSE]
}


# =============================================================================
# build_database_wide
# =============================================================================
#
# Wide-format input: one row per (respondent, task), with attributes spread
# across alternative-suffixed columns.
#
# Args
#   data            data.frame or path to CSV.
#   id_col          column identifying the respondent.
#   task_col        column identifying the choice task within respondent. If
#                   missing or NULL, tasks are numbered 1..T per respondent.
#   choice_col      column with the chosen alternative index (1..J).
#   attributes      named list. Each entry is a length-J character vector
#                   giving the column name in `data` for that attribute on
#                   each alternative. Use NA in a slot to encode a structural
#                   zero (attribute does not apply to that alternative).
#                   Example for mode choice with 4 alts (car/bus/air/rail):
#                     list(
#                       time    = c("time_car","time_bus","time_air","time_rail"),
#                       access  = c(NA, "access_bus", "access_air", "access_rail"),
#                       service = c(NA, NA, "service_air", "service_rail")
#                     )
#   price           length-J character vector of price/cost column names.
#   availability    optional length-J character vector of availability columns.
#                   Tasks where any avail==0 are dropped, with a warning. NULL
#                   = assume all alternatives available everywhere.
#   scalings        optional named list of per-attribute scalings. Names match
#                   entries in `attributes` and/or the string "price". Each
#                   value is a scalar that the corresponding columns are
#                   divided by. Example: list(time = 60, price = 10).
#   verbose         print progress and filter rates.
# =============================================================================

build_database_wide <- function(data,
                                id_col, task_col = NULL, choice_col,
                                attributes, price,
                                availability = NULL,
                                scalings = NULL,
                                verbose = TRUE) {
  stopifnot(is.list(attributes), length(attributes) >= 1,
            !is.null(names(attributes)),
            all(nzchar(names(attributes))))

  J <- length(price)
  for (nm in names(attributes)) {
    if (length(attributes[[nm]]) != J) {
      stop(sprintf("attributes[['%s']] has length %d, expected %d (= length(price)).",
                   nm, length(attributes[[nm]]), J))
    }
  }
  if (!is.null(availability) && length(availability) != J) {
    stop(sprintf("availability has length %d, expected %d.",
                 length(availability), J))
  }

  raw <- if (is.character(data) && length(data) == 1) {
    if (grepl("\\.dat$", data, ignore.case = TRUE)) {
      read.table(data, header = TRUE, stringsAsFactors = FALSE)
    } else {
      read.csv(data, stringsAsFactors = FALSE)
    }
  } else {
    as.data.frame(data)
  }

  attr_cols_used <- unlist(lapply(attributes, function(v) v[!is.na(v)]),
                           use.names = FALSE)
  required <- unique(c(id_col, if (!is.null(task_col)) task_col else NULL,
                       choice_col,
                       attr_cols_used, price,
                       if (!is.null(availability)) availability else NULL))
  miss <- setdiff(required, names(raw))
  if (length(miss) > 0) {
    stop("Columns missing from data: ", paste(miss, collapse = ", "))
  }

  if (is.null(scalings)) scalings <- list()

  raw[[".resp"]] <- as.integer(factor(raw[[id_col]]))
  if (is.null(task_col)) {
    raw <- raw[order(raw$.resp), ]
    raw[[".task_raw"]] <- ave(raw$.resp, raw$.resp, FUN = seq_along)
  } else {
    raw[[".task_raw"]] <- raw[[task_col]]
  }
  raw[[".task_idx"]] <- as.integer(factor(paste(raw$.resp, raw$.task_raw, sep = "_")))

  if (!is.null(availability)) {
    for (j in seq_len(J)) {
      raw[[availability[j]]] <- .coerce_numeric_zero_na(raw[[availability[j]]])
    }
    avail_sum_per_row <- rowSums(as.matrix(raw[, availability, drop = FALSE]))
    keep <- avail_sum_per_row == J
    n_total <- nrow(raw); n_keep <- sum(keep)
    if (verbose && n_keep < n_total) {
      cat(sprintf("Availability filter: dropping %d of %d tasks (%.1f%%) where not all %d alternatives are available.\n",
                  n_total - n_keep, n_total,
                  100 * (n_total - n_keep) / n_total, J))
    }
    if (n_keep == 0) stop("No tasks remain after availability filter.")
    raw <- raw[keep, , drop = FALSE]
  }

  # Balanced panel
  raw[[".resp"]] <- as.integer(factor(raw[[".resp"]]))
  tasks_per_id <- table(raw[[".resp"]])
  T_max <- max(tasks_per_id)
  full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == T_max])
  if (length(full_ids) < length(tasks_per_id)) {
    n_drop <- length(tasks_per_id) - length(full_ids)
    if (verbose) {
      cat(sprintf("Balanced-panel filter: dropping %d of %d respondents (%.1f%%) with fewer than %d tasks.\n",
                  n_drop, length(tasks_per_id),
                  100 * n_drop / length(tasks_per_id), T_max))
    }
    raw <- raw[raw[[".resp"]] %in% full_ids, , drop = FALSE]
    raw[[".resp"]] <- as.integer(factor(raw[[".resp"]]))
  }

  raw <- raw[order(raw[[".resp"]], raw[[".task_idx"]]), ]
  N <- length(unique(raw[[".resp"]]))
  T_const <- as.integer(nrow(raw) / N)

  attr_names <- names(attributes)
  Ng <- length(attr_names)

  db <- data.frame(ID = rep(seq_len(N), each = T_const),
                   TASK = rep(seq_len(T_const), times = N))

  for (j in seq_len(J)) {
    for (a in seq_along(attr_names)) {
      col_nm <- attributes[[attr_names[a]]][j]
      vals <- if (is.na(col_nm)) {
        rep(0, nrow(raw))
      } else {
        .coerce_numeric_zero_na(raw[[col_nm]])
      }
      sc <- scalings[[attr_names[a]]]
      if (!is.null(sc) && sc != 1) vals <- vals / sc
      db[[sprintf("x%d_%d", a, j)]] <- vals
    }
    pvals <- .coerce_numeric_zero_na(raw[[price[j]]])
    psc <- scalings[["price"]]
    if (!is.null(psc) && psc != 1) pvals <- pvals / psc
    db[[sprintf("price_%d", j)]] <- pvals
  }

  ch <- as.integer(raw[[choice_col]])
  if (any(is.na(ch)) || any(ch < 1) || any(ch > J)) {
    stop(sprintf("`%s`: choice values must be integers in 1..%d.", choice_col, J))
  }
  db$CHOICE <- ch

  scaled_label <- function(name, default) {
    sc <- scalings[[name]]
    if (is.null(sc) || sc == 1) default else sprintf("%s_per%g", default, sc)
  }
  attr(db, "attr_labels") <- c(
    vapply(attr_names, function(nm) scaled_label(nm, nm), character(1)),
    scaled_label("price", "price")
  )
  attr(db, "n_alternatives") <- J
  attr(db, "n_generic")      <- Ng

  if (verbose) {
    cat(sprintf("Built canonical database: N=%d respondents, T=%d tasks, J=%d alternatives, %d generic attributes + price.\n",
                N, T_const, J, Ng))
  }
  db
}


# =============================================================================
# build_database  (dispatcher)
# =============================================================================
#
# Dispatches to build_database_long or build_database_wide based on `format`.
# If `format = "auto"`, infers: if `attributes` is a named list -> "wide";
# if `alt_col` is provided -> "long".
# =============================================================================

build_database <- function(data,
                          format = c("auto", "long", "wide"),
                          ...) {
  format <- match.arg(format)
  args <- list(...)
  if (format == "auto") {
    if (!is.null(args$attributes) && is.list(args$attributes)) {
      format <- "wide"
    } else if (!is.null(args$alt_col)) {
      format <- "long"
    } else {
      stop("Could not infer format. Pass format = 'long' or 'wide' explicitly, ",
           "or supply `alt_col` (long) or `attributes` as a named list (wide).")
    }
  }
  if (format == "long") {
    do.call(build_database_long, c(list(data = data), args))
  } else {
    do.call(build_database_wide, c(list(data = data), args))
  }
}


# =============================================================================
# run_lcmnl_workflow: "plug in data -> all results" entry point
# =============================================================================
#
# Three calling conventions:
#
# (1) Canonical database already built:
#         run_lcmnl_workflow(database = db, C_cands = 1:6, ...)
#
# (2) Long format + column mapping (forwards to build_database_long):
#         run_lcmnl_workflow(data = "x.csv", format = "long",
#                            id_col, task_col, alt_col, choice_col,
#                            attribute_cols, price_col, ...)
#
# (3) Wide format + column mapping (forwards to build_database_wide):
#         run_lcmnl_workflow(data = df, format = "wide",
#                            id_col, choice_col, attributes, price, ...)
#
# Additional args
#   C_cands         integer vector of class counts. Default 1:6.
#   run_mmnl        whether to estimate the MMNL benchmark. Default TRUE.
#   attr_labels     override display names in the betas CSV. Default: taken
#                   from attr(database, "attr_labels").
#   output_prefix   filename prefix for the CSVs. Default "workflow".
#   output_dir      output directory. Default OUTPUT_DIR ("output").
#   write_csv       write CSVs to disk. Default TRUE.
#   verbose         per-C progress lines and final printout. Default TRUE.
#
# Returns invisibly a list with: database, dgp, lcmnl, mmnl, summary,
# class_betas, comparison, best_C, best_lcmnl.
# =============================================================================

run_lcmnl_workflow <- function(database = NULL,
                               data = NULL,
                               format = c("auto", "long", "wide"),
                               C_cands = 1:6,
                               run_mmnl = TRUE,
                               attr_labels = NULL,
                               output_prefix = "workflow",
                               output_dir = NULL,
                               write_csv = TRUE,
                               verbose = TRUE,
                               ...) {
  format <- match.arg(format)

  if (is.null(database)) {
    if (is.null(data)) {
      stop("Provide either `database` (canonical wide format) or `data` plus ",
           "the column-mapping arguments (build_database is called internally).")
    }
    database <- build_database(data, format = format, verbose = verbose, ...)
  } else {
    if (!all(c("ID", "TASK", "CHOICE") %in% names(database))) {
      stop("`database` must contain columns ID, TASK, CHOICE (canonical wide format).")
    }
  }

  J  <- attr(database, "n_alternatives")
  Ng <- attr(database, "n_generic")
  if (is.null(J) || is.null(Ng)) {
    price_cols <- grep("^price_\\d+$", names(database), value = TRUE)
    J  <- length(price_cols)
    x1 <- grep("^x\\d+_1$", names(database), value = TRUE)
    Ng <- length(x1)
    if (J < 2 || Ng < 1) {
      stop("Could not infer n_alternatives / n_generic from database columns.")
    }
  }
  dgp <- make_dgp_config(n_generic = Ng, n_alternatives = J)

  if (is.null(attr_labels)) {
    attr_labels <- attr(database, "attr_labels")
    if (is.null(attr_labels)) attr_labels <- c(paste0("x", seq_len(Ng)), "price")
  }
  stopifnot(length(attr_labels) == Ng + 1)

  if (is.null(output_dir)) output_dir <- OUTPUT_DIR
  if (write_csv && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  results <- list()
  for (cc in C_cands) {
    if (verbose) cat(sprintf("\n=== Estimating C = %d ===\n", cc))
    t0 <- Sys.time()
    m <- estimate_lcmnl_multistart(database, cc, dgp = dgp)
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

  mmnl_fit <- NULL
  comparison_df <- NULL
  if (run_mmnl) {
    if (verbose) cat("\n=== Estimating MMNL (independent normals) ===\n")
    t0 <- Sys.time()
    mmnl_fit <- tryCatch(estimate_mmnl(database, dgp = dgp),
                         error = function(e) {
                           if (verbose) cat("  MMNL failed:", conditionMessage(e), "\n")
                           NULL
                         })
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(mmnl_fit) && verbose) {
      cat(sprintf("  converged=%s  LL=%.2f  BIC=%.2f  AIC=%.2f  k=%d  time=%.1fs\n",
                  mmnl_fit$converged, mmnl_fit$LL, mmnl_fit$BIC,
                  mmnl_fit$AIC, mmnl_fit$k, dt))
    }
    if (!is.null(mmnl_fit) && "1" %in% names(results)) {
      mnl <- results[["1"]]
      comparison_df <- data.frame(
        model = c("MNL (C=1)",
                  sprintf("LCMNL (C=%d, BIC-best)", best_C),
                  "MMNL (independent normals)"),
        LL  = c(mnl$LL,  best$LL,  mmnl_fit$LL),
        k   = c(mnl$k,   best$k,   mmnl_fit$k),
        BIC = c(mnl$BIC, best$BIC, mmnl_fit$BIC),
        AIC = c(mnl$AIC, best$AIC, mmnl_fit$AIC),
        stringsAsFactors = FALSE
      )
      comparison_df$dBIC <- comparison_df$BIC - min(comparison_df$BIC)
    }
  }

  if (write_csv) {
    write.csv(summary_df,
              file.path(output_dir, paste0(output_prefix, "_summary.csv")),
              row.names = FALSE)
    write.csv(betas_df,
              file.path(output_dir, paste0(output_prefix, "_class_betas.csv")),
              row.names = FALSE)
    if (!is.null(comparison_df)) {
      write.csv(comparison_df,
                file.path(output_dir, paste0(output_prefix, "_model_comparison.csv")),
                row.names = FALSE)
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
      cat("\n--- Please cite klue if you use it in published work ---\n")
      cat("  Frings, O. (2026). A Hybrid Machine Learning and Random Utility\n")
      cat("  Framework for Latent Class Model Specification.\n")
      cat("  Journal of Choice Modelling.\n")
      cat("  citation(\"klue\")  # full BibTeX entry\n")
    }
  }

  invisible(list(
    database    = database,
    dgp         = dgp,
    lcmnl       = results,
    mmnl        = mmnl_fit,
    summary     = summary_df,
    class_betas = betas_df,
    comparison  = comparison_df,
    best_C      = best_C,
    best_lcmnl  = best
  ))
}
