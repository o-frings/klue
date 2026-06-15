# Data harmonisers: long- or wide-format discrete-choice data -> the
# canonical wide format the engine consumes:
#   ID, TASK, x1_1..xN_J, price_1..price_J, CHOICE
# (N = n_generic, J = n_alternatives, CHOICE in 1..J, last alt = reference.)
# The engine assumes every alternative available in every task; tasks failing
# the availability filter and respondents breaking the balanced panel are
# dropped with a reported rate.

.coerce_numeric_zero_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
}

.read_raw <- function(data) {
  if (is.character(data) && length(data) == 1) {
    if (grepl("\\.dat$", data, ignore.case = TRUE)) {
      read.table(data, header = TRUE, stringsAsFactors = FALSE)
    } else {
      read.csv(data, stringsAsFactors = FALSE)
    }
  } else {
    as.data.frame(data)
  }
}

.require_cols <- function(raw, required) {
  miss <- setdiff(required, names(raw))
  if (length(miss) > 0)
    stop("Columns missing from data: ", paste(miss, collapse = ", "))
}

.report_drop <- function(verbose, what, n_drop, n_total, detail) {
  if (verbose && n_drop > 0) {
    cat(sprintf("%s: dropping %d of %d %s.\n",
                what, n_drop, n_total, detail))
  }
}

# Label helper: attribute name, suffixed by its scaling if != 1.
.scaled_label <- function(name, scalings, default = name) {
  sc <- scalings[[name]]
  if (is.null(sc) || sc == 1) default else sprintf("%s_per%g", default, sc)
}

.stamp_db <- function(db, attr_labels, J, Ng, verbose) {
  attr(db, "attr_labels")    <- attr_labels
  attr(db, "n_alternatives") <- J
  attr(db, "n_generic")      <- Ng
  if (verbose) {
    cat(sprintf("Built canonical database: N=%d respondents, T=%d tasks, J=%d alternatives, %d generic attributes + price.\n",
                length(unique(db$ID)), max(db$TASK), J, Ng))
  }
  db
}

#' Build the canonical database from long-format data
#'
#' One row per (respondent, task, alternative). `choice_format` is
#' "indicator" (0/1 on the chosen row) or "alt_index" (chosen alt repeated on
#' every row of the task). `scalings` is a named list of divisors keyed by
#' attribute/price column name; `price_scaling` is a shortcut for the price
#' column. Tasks where any alternative is unavailable (avail_col == 0) and
#' respondents with fewer tasks than the panel maximum are dropped.
#' @param data A data.frame, or a single file path to a `.dat`
#'   (whitespace-delimited, header) or `.csv` file to read.
#' @param id_col Name of the respondent identifier column.
#' @param task_col Name of the choice-task identifier column (unique within
#'   respondent).
#' @param alt_col Name of the alternative identifier column within a task.
#' @param choice_col Name of the choice column, interpreted per
#'   `choice_format`.
#' @param attribute_cols Character vector of generic attribute column names;
#'   each becomes the `x*_*` block in the order supplied.
#' @param price_col Length-one name of the price column.
#' @param choice_format Either "indicator" (0/1, one chosen row per task) or
#'   "alt_index" (the chosen alternative index repeated on every row of the
#'   task).
#' @param avail_col Optional name of an availability column (0/1); tasks where
#'   not all alternatives are available are dropped.
#' @param scalings Optional named list of positive divisors keyed by attribute
#'   or price column name; matching columns are divided before use.
#' @param price_scaling Positive divisor applied to the price column; shortcut
#'   for adding the price column to `scalings`.
#' @param verbose Logical; if `TRUE`, print drop counts and a build summary.
#' @return A data.frame in the canonical wide format with columns `ID`,
#'   `TASK`, the generic attribute blocks `x1_1..xN_J`, the price columns
#'   `price_1..price_J`, and `CHOICE` (the chosen alternative index in
#'   `1..J`). The last alternative is the reference. Attributes `attr_labels`
#'   (attribute and price labels, suffixed by scaling where applicable),
#'   `n_alternatives` (J), and `n_generic` (number of generic attributes) are
#'   attached.
#' @export
klue_database_long <- function(data,
                               id_col, task_col, alt_col, choice_col,
                               attribute_cols, price_col,
                               choice_format = c("indicator", "alt_index"),
                               avail_col = NULL,
                               scalings = NULL,
                               price_scaling = 1,
                               verbose = TRUE) {
  choice_format <- match.arg(choice_format)
  stopifnot(length(attribute_cols) >= 1, length(price_col) == 1)
  raw <- .read_raw(data)
  .require_cols(raw, c(id_col, task_col, alt_col, choice_col, attribute_cols,
                       price_col, avail_col))

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

  # Availability filter: keep tasks where all J alternatives are available.
  if (!is.null(avail_col)) {
    avail <- .coerce_numeric_zero_na(raw[[avail_col]])
    total_avail <- tapply(avail, raw$.task_idx, sum)
    full_tasks <- as.integer(names(total_avail)[total_avail == J])
    n_total <- length(unique(raw$.task_idx))
    .report_drop(verbose, "Availability filter", n_total - length(full_tasks),
                 n_total, sprintf("tasks where not all %d alternatives are available", J))
    if (length(full_tasks) == 0) stop("No tasks remain after availability filter.")
    raw <- raw[raw$.task_idx %in% full_tasks, , drop = FALSE]
  }

  # Malformed tasks: not exactly J rows.
  task_sizes <- table(raw$.task_idx)
  bad_tasks <- as.integer(names(task_sizes)[task_sizes != J])
  .report_drop(verbose, "Malformed-task filter", length(bad_tasks),
               length(task_sizes), sprintf("tasks without exactly %d alternatives", J))
  if (length(bad_tasks) > 0) raw <- raw[!raw$.task_idx %in% bad_tasks, , drop = FALSE]
  if (nrow(raw) == 0) stop("No tasks remain after filtering.")
  raw$.task_idx <- as.integer(factor(raw$.task_idx))
  raw$.resp     <- as.integer(factor(raw$.resp))

  # Balanced panel: keep respondents with the maximum task count.
  tasks_per_id <- tapply(raw$.task_idx, raw$.resp, function(x) length(unique(x)))
  if (length(unique(tasks_per_id)) != 1) {
    T_max <- max(tasks_per_id)
    full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == T_max])
    .report_drop(verbose, "Balanced-panel filter",
                 length(tasks_per_id) - length(full_ids), length(tasks_per_id),
                 sprintf("respondents with fewer than %d tasks", T_max))
    raw <- raw[raw$.resp %in% full_ids, , drop = FALSE]
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

  labels <- c(vapply(attribute_cols, function(cl) .scaled_label(cl, scalings),
                     character(1)),
              .scaled_label(price_col, scalings))
  .stamp_db(db, labels, J, length(attribute_cols), verbose)
}

#' Build the canonical database from wide-format data
#'
#' One row per (respondent, task). `attribute_cols` is a named list; each
#' entry is a length-J character vector of column names (NA = structural zero
#' on that alternative). `price_col` is a length-J vector; `avail_col`
#' optionally too. `scalings` is keyed by attribute name or "price".
#' @param data A data.frame, or a single file path to a `.dat`
#'   (whitespace-delimited, header) or `.csv` file to read.
#' @param id_col Name of the respondent identifier column.
#' @param task_col Optional name of the choice-task identifier column; if
#'   `NULL`, tasks are numbered sequentially within each respondent.
#' @param choice_col Name of the column holding the chosen alternative index
#'   in `1..J`.
#' @param attribute_cols Named list of length-J character vectors of column
#'   names, one vector per generic attribute; `NA` marks a structural zero for
#'   that alternative. The list element names become attribute labels.
#' @param price_col Length-J character vector of per-alternative price column
#'   names.
#' @param avail_col Optional length-J character vector of per-alternative
#'   availability column names (0/1); tasks where not all alternatives are
#'   available are dropped.
#' @param scalings Optional named list of positive divisors keyed by attribute
#'   name or "price"; matching columns are divided before use.
#' @param attributes Deprecated alias for `attribute_cols`.
#' @param price Deprecated alias for `price_col`.
#' @param availability Deprecated alias for `avail_col`.
#' @param verbose Logical; if `TRUE`, print drop counts and a build summary.
#' @return A data.frame in the canonical wide format with columns `ID`,
#'   `TASK`, the generic attribute blocks `x1_1..xN_J`, the price columns
#'   `price_1..price_J`, and `CHOICE` (the chosen alternative index in
#'   `1..J`). The last alternative is the reference. Attributes `attr_labels`
#'   (attribute and price labels, suffixed by scaling where applicable),
#'   `n_alternatives` (J), and `n_generic` (number of generic attributes) are
#'   attached.
#' @export
klue_database_wide <- function(data,
                               id_col, task_col = NULL, choice_col,
                               attribute_cols = NULL, price_col = NULL,
                               avail_col = NULL,
                               scalings = NULL,
                               # Deprecated aliases (kept for backward compat):
                               attributes = NULL, price = NULL,
                               availability = NULL,
                               verbose = TRUE) {
  if (is.null(attribute_cols) && !is.null(attributes))   attribute_cols <- attributes
  if (is.null(price_col)      && !is.null(price))        price_col      <- price
  if (is.null(avail_col)      && !is.null(availability)) avail_col      <- availability
  if (is.null(attribute_cols)) {
    stop("`attribute_cols` is required: a named list of length-J ",
         "character vectors of column names (one vector per attribute).")
  }
  if (is.null(price_col)) {
    stop("`price_col` is required: a length-J character vector of price ",
         "column names (one per alternative).")
  }
  stopifnot(is.list(attribute_cols), length(attribute_cols) >= 1,
            !is.null(names(attribute_cols)), all(nzchar(names(attribute_cols))))

  J <- length(price_col)
  for (nm in names(attribute_cols)) {
    if (length(attribute_cols[[nm]]) != J) {
      stop(sprintf("attribute_cols[['%s']] has length %d, expected %d (= length(price_col)).",
                   nm, length(attribute_cols[[nm]]), J))
    }
  }
  if (!is.null(avail_col) && length(avail_col) != J) {
    stop(sprintf("avail_col has length %d, expected %d.", length(avail_col), J))
  }

  raw <- .read_raw(data)
  attr_cols_used <- unlist(lapply(attribute_cols, function(v) v[!is.na(v)]),
                           use.names = FALSE)
  .require_cols(raw, unique(c(id_col, task_col, choice_col, attr_cols_used,
                              price_col, avail_col)))
  if (is.null(scalings)) scalings <- list()

  raw$.resp <- as.integer(factor(raw[[id_col]]))
  if (is.null(task_col)) {
    raw <- raw[order(raw$.resp), ]
    raw$.task_raw <- ave(raw$.resp, raw$.resp, FUN = seq_along)
  } else {
    raw$.task_raw <- raw[[task_col]]
  }
  raw$.task_idx <- as.integer(factor(paste(raw$.resp, raw$.task_raw, sep = "_")))

  if (!is.null(avail_col)) {
    for (j in seq_len(J)) raw[[avail_col[j]]] <- .coerce_numeric_zero_na(raw[[avail_col[j]]])
    keep <- rowSums(as.matrix(raw[, avail_col, drop = FALSE])) == J
    .report_drop(verbose, "Availability filter", sum(!keep), nrow(raw),
                 sprintf("tasks where not all %d alternatives are available", J))
    if (!any(keep)) stop("No tasks remain after availability filter.")
    raw <- raw[keep, , drop = FALSE]
  }

  # Balanced panel
  raw$.resp <- as.integer(factor(raw$.resp))
  tasks_per_id <- table(raw$.resp)
  T_max <- max(tasks_per_id)
  full_ids <- as.integer(names(tasks_per_id)[tasks_per_id == T_max])
  .report_drop(verbose, "Balanced-panel filter",
               length(tasks_per_id) - length(full_ids), length(tasks_per_id),
               sprintf("respondents with fewer than %d tasks", T_max))
  if (length(full_ids) < length(tasks_per_id)) {
    raw <- raw[raw$.resp %in% full_ids, , drop = FALSE]
    raw$.resp <- as.integer(factor(raw$.resp))
  }

  raw <- raw[order(raw$.resp, raw$.task_idx), ]
  N <- length(unique(raw$.resp))
  T_const <- as.integer(nrow(raw) / N)
  attr_names <- names(attribute_cols)

  db <- data.frame(ID = rep(seq_len(N), each = T_const),
                   TASK = rep(seq_len(T_const), times = N))
  for (j in seq_len(J)) {
    for (a in seq_along(attr_names)) {
      col_nm <- attribute_cols[[attr_names[a]]][j]
      vals <- if (is.na(col_nm)) rep(0, nrow(raw))
              else .coerce_numeric_zero_na(raw[[col_nm]])
      sc <- scalings[[attr_names[a]]]
      if (!is.null(sc) && sc != 1) vals <- vals / sc
      db[[sprintf("x%d_%d", a, j)]] <- vals
    }
    pvals <- .coerce_numeric_zero_na(raw[[price_col[j]]])
    psc <- scalings[["price"]]
    if (!is.null(psc) && psc != 1) pvals <- pvals / psc
    db[[sprintf("price_%d", j)]] <- pvals
  }

  ch <- as.integer(raw[[choice_col]])
  if (any(is.na(ch)) || any(ch < 1) || any(ch > J)) {
    stop(sprintf("`%s`: choice values must be integers in 1..%d.", choice_col, J))
  }
  db$CHOICE <- ch

  labels <- c(vapply(attr_names, function(nm) .scaled_label(nm, scalings),
                     character(1)),
              .scaled_label("price", scalings))
  .stamp_db(db, labels, J, length(attr_names), verbose)
}

#' Build the canonical database (format dispatcher)
#'
#' `format = "auto"` infers: named-list `attribute_cols` -> wide;
#' `alt_col` present -> long.
#' @param data A data.frame, or a single file path to a `.dat`
#'   (whitespace-delimited, header) or `.csv` file to read.
#' @param format One of "auto", "long", or "wide". With "auto", a named-list
#'   `attribute_cols` selects the wide builder and an `alt_col` argument
#'   selects the long builder.
#' @param ... Further arguments passed to [klue_database_long()] or
#'   [klue_database_wide()] according to the resolved format.
#' @return A data.frame in the canonical wide format with columns `ID`,
#'   `TASK`, the generic attribute blocks `x1_1..xN_J`, the price columns
#'   `price_1..price_J`, and `CHOICE` (the chosen alternative index in
#'   `1..J`). The last alternative is the reference. Attributes `attr_labels`,
#'   `n_alternatives`, and `n_generic` are attached.
#' @export
klue_database <- function(data, format = c("auto", "long", "wide"), ...) {
  format <- match.arg(format)
  args <- list(...)
  if (format == "auto") {
    attrs_arg <- args$attribute_cols
    if (is.null(attrs_arg)) attrs_arg <- args$attributes  # accept old name
    if (!is.null(attrs_arg) && is.list(attrs_arg)) {
      format <- "wide"
    } else if (!is.null(args$alt_col)) {
      format <- "long"
    } else {
      stop("Could not infer format. Pass format = 'long' or 'wide' explicitly, ",
           "or supply `alt_col` (long) or `attribute_cols` as a named list (wide).")
    }
  }
  fn <- if (format == "long") klue_database_long else klue_database_wide
  do.call(fn, c(list(data = data), args))
}
