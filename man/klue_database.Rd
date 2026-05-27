\name{klue_database}
\alias{klue_database}
\alias{klue_database_long}
\alias{klue_database_wide}
\alias{build_database}
\alias{build_database_long}
\alias{build_database_wide}
\title{Build a canonical wide database from long- or wide-format input}
\description{
  Reshape an arbitrary discrete-choice dataset into the canonical wide
  format expected by the estimation engine. Both helpers also handle
  availability filtering, per-attribute scaling, and balanced-panel
  enforcement.

  \code{build_database}, \code{build_database_long} and
  \code{build_database_wide} are deprecated aliases.

  Argument naming is homogenised across long and wide format: the same
  identifiers \code{attribute_cols}, \code{price_col}, \code{avail_col}
  are used in both, though their shape differs (scalar in long, vector
  of \eqn{J} column names in wide).
}
\usage{
klue_database(data, format = c("auto", "long", "wide"), \ldots)

klue_database_long(data,
                   id_col, task_col, alt_col, choice_col,
                   attribute_cols, price_col,
                   choice_format = c("indicator", "alt_index"),
                   avail_col = NULL,
                   scalings = NULL,
                   price_scaling = 1,
                   verbose = TRUE)

klue_database_wide(data,
                   id_col, task_col = NULL, choice_col,
                   attribute_cols = NULL, price_col = NULL,
                   avail_col = NULL,
                   scalings = NULL,
                   attributes = NULL, price = NULL,
                   availability = NULL,
                   verbose = TRUE)
}
\arguments{
  \item{data}{A data frame or path to a CSV (or \code{.dat}) file.}
  \item{format}{\code{"long"} (one row per (id, task, alternative)) or
    \code{"wide"} (one row per (id, task) with alternative-suffixed columns).
    \code{"auto"} infers from whether \code{attribute_cols} is a named list
    (\eqn{\Rightarrow} wide) or \code{alt_col} is supplied (\eqn{\Rightarrow} long).}
  \item{id_col, task_col, alt_col, choice_col}{Column names. For wide format,
    \code{task_col = NULL} synthesises 1..T per respondent;
    \code{choice_col} is the chosen alternative index (1..J).}
  \item{attribute_cols}{In long format, a character vector of generic
    attribute column names. In wide format, a named list mapping each
    attribute to a length-\eqn{J} character vector of column names (one
    per alternative). \code{NA} in a wide slot encodes a structural zero
    on that alternative.}
  \item{price_col}{In long format, the single price column name. In wide
    format, a length-\eqn{J} character vector of price column names.}
  \item{choice_format}{(long) \code{"indicator"} = 0/1 chosen-row flag;
    \code{"alt_index"} = chosen-alt integer repeated on every row of the
    (id, task) group.}
  \item{avail_col}{Optional availability column(s). Scalar in long, length-\eqn{J}
    in wide. Tasks where any alternative has \code{avail == 0} are dropped,
    with a drop-rate message.}
  \item{scalings}{Optional named list of per-attribute scalings; each
    column is divided by the corresponding scalar. Keys: in long, match
    the entries of \code{attribute_cols} / \code{price_col}; in wide,
    match \code{names(attribute_cols)} or the literal \code{"price"}.}
  \item{price_scaling}{(long, shortcut) Equivalent to
    \code{scalings = list(<price_col> = price_scaling)}.}
  \item{attributes, price, availability}{Deprecated aliases for
    \code{attribute_cols}, \code{price_col}, \code{avail_col} in wide
    format. Accepted silently for backward compatibility; please switch
    to the canonical names.}
  \item{verbose}{Print build progress, filter rates, and final dimensions.}
  \item{\ldots}{Passed to \code{klue_database_long} or
    \code{klue_database_wide}.}
}
\value{
  A data frame with columns \code{ID}, \code{TASK}, \code{x1_1..xN_J},
  \code{price_1..price_J}, \code{CHOICE}. Attributes set:
  \code{attr_labels} (display names for the betas table),
  \code{n_alternatives}, \code{n_generic}.
}
\examples{
\dontrun{
# Long format
db <- klue_database_long(
  data = "data.csv",
  id_col = "respondent_id", task_col = "task",
  alt_col = "alternative", choice_col = "chosen",
  attribute_cols = c("attr1", "attr2"),
  price_col = "price",
  scalings = list(price = 10)
)

# Wide format
db <- klue_database_wide(
  data = mydata,
  id_col = "ID", task_col = "task", choice_col = "choice",
  attribute_cols = list(
    time = c("time_a", "time_b", "time_c"),
    qual = c("qual_a", "qual_b", "qual_c")
  ),
  price_col = c("cost_a", "cost_b", "cost_c"),
  avail_col = c("av_a", "av_b", "av_c"),
  scalings = list(time = 60, price = 10)
)
}
}
\seealso{\code{\link{klue}}}
