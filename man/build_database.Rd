\name{build_database}
\alias{build_database}
\alias{build_database_long}
\alias{build_database_wide}
\title{Build a canonical wide database from long- or wide-format input}
\description{
  Reshape an arbitrary discrete-choice dataset into the canonical wide
  format expected by the estimation engine. Both helpers also handle
  availability filtering, per-attribute scaling, and balanced-panel
  enforcement.
}
\usage{
build_database(data, format = c("auto", "long", "wide"), \ldots)

build_database_long(data,
                    id_col, task_col, alt_col, choice_col,
                    attribute_cols, price_col,
                    choice_format = c("indicator", "alt_index"),
                    avail_col = NULL,
                    scalings = NULL,
                    price_scaling = 1,
                    verbose = TRUE)

build_database_wide(data,
                    id_col, task_col = NULL, choice_col,
                    attributes, price,
                    availability = NULL,
                    scalings = NULL,
                    verbose = TRUE)
}
\arguments{
  \item{data}{A data frame or path to a CSV (or \code{.dat}) file.}
  \item{format}{\code{"long"} (one row per (id, task, alternative)) or
    \code{"wide"} (one row per (id, task) with alternative-suffixed columns).
    \code{"auto"} infers from whether \code{attributes} (list) or
    \code{alt_col} is supplied.}
  \item{id_col, task_col, alt_col, choice_col}{Column names. For wide format,
    \code{task_col = NULL} synthesises 1..T per respondent. For wide format
    \code{choice_col} is the chosen alternative index (1..J).}
  \item{attribute_cols, price_col}{(long) Character vectors of attribute
    column names and a single price column.}
  \item{attributes, price}{(wide) Named list mapping each attribute to a
    length-J character vector of column names (\code{NA} encodes a
    structural zero on that alternative); and a length-J vector of price
    column names.}
  \item{choice_format}{(long) \code{"indicator"} = 0/1 chosen-row flag;
    \code{"alt_index"} = chosen-alt integer repeated on every row of the
    (id, task) group.}
  \item{avail_col, availability}{Optional availability columns. Tasks where
    any alternative has \code{avail == 0} are dropped, with a drop-rate
    message. The engine does not support partial availability.}
  \item{scalings}{Optional named list of per-attribute scalings; each
    column is divided by the corresponding scalar. Names match the entries
    of \code{attribute_cols} / \code{price_col} (long) or
    \code{names(attributes)} / the literal \code{"price"} (wide).}
  \item{price_scaling}{(long, shortcut) Equivalent to
    \code{scalings = list(<price_col> = price_scaling)}.}
  \item{verbose}{Print build progress, filter rates, and final dimensions.}
  \item{\ldots}{Passed to \code{build_database_long} or
    \code{build_database_wide}.}
}
\value{
  A data frame with columns \code{ID}, \code{TASK}, \code{x1_1..xN_J},
  \code{price_1..price_J}, \code{CHOICE}. Attributes set:
  \code{attr_labels} (display names for the betas table),
  \code{n_alternatives}, \code{n_generic}.
}
\examples{
\dontrun{
db <- build_database_long(
  data = "data.csv",
  id_col = "respondent_id", task_col = "task",
  alt_col = "alternative", choice_col = "chosen",
  attribute_cols = c("attr1", "attr2"),
  price_col = "price", price_scaling = 10
)
}
}
\seealso{\code{\link{run_lcmnl_workflow}}}
