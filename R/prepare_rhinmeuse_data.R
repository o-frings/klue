# ============================================================================ #
#   R/prepare_rhinmeuse_data.R                                                 #
#   Rebuilds the analysis file of the Rhin-Meuse water-quality DCE             #
#   (test_data_files/data.csv, loaded by R/empirical_application.R) from the   #
#   public deposit and the survey authors' choice-card design.                 #
#                                                                              #
#   Inputs (test_data_files/):                                                 #
#     MXDUB2/Non-use value of improving water quality and biodiversity.csv     #
#         public deposit doi:10.57745/MXDUB2 (Amiri, Abildtrup, Garcia and     #
#         Montagne-Huck 2024): 1,365 survey answers, one row per respondent    #
#     rhinmeuse_design.csv                                                     #
#         the choice-card design from the survey's authors: 3 versions x 8     #
#         cards x 3 alternatives (alt 1 = no project), with the eight policy   #
#         attributes, Waterbill and ASCsq                                      #
#     rhinmeuse_kept.csv                                                       #
#         the 930 respondents the authors kept after their Stata cleaning      #
#         (deposit files Stata_file_*.do), with the design version each saw    #
#                                                                              #
#   Output: test_data_files/data_rebuilt.csv, one row per respondent x card x  #
#   alternative (22,320 rows), with the columns the analysis reads.            #
#   With --check, the script compares the rebuilt file with data.csv on those  #
#   columns and stops on any difference.                                       #
#                                                                              #
#   Run from the project folder: Rscript R/prepare_rhinmeuse_data.R [--check]  #
# ============================================================================ #

dir_in   <- "test_data_files"
deposit  <- file.path(dir_in, "MXDUB2",
                      "Non-use value of improving water quality and biodiversity.csv")
design_f <- file.path(dir_in, "rhinmeuse_design.csv")
kept_f   <- file.path(dir_in, "rhinmeuse_kept.csv")
out_f    <- file.path(dir_in, "data_rebuilt.csv")

for (f in c(deposit, design_f, kept_f))
  if (!file.exists(f)) stop("missing input: ", f,
                            " (see README_REPRODUCE.md, section Data)")

dep    <- read.csv(deposit, check.names = FALSE, stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8")
design <- read.csv(design_f, stringsAsFactors = FALSE)
kept   <- read.csv(kept_f, stringsAsFactors = FALSE)

id_col    <- "id. ID de la réponse"   # survey response id
panel_col <- "ID."                     # hashed respondent id used by the analysis
answer_col <- function(cs) sprintf("Q%dbloc3[a]. Choice set %d", 15 + cs, cs)
answers   <- c("Absence de projet" = 1L, "Projet A" = 2L, "Projet B" = 3L)
attr_cols <- c("Ban_pesticides", "RenaturationYes", "Local_Negative_Effect_Ferti",
               "Mixedfertilizers", "Hedgesconservation", "Hedgesplantation",
               "Forest_Management_For_Water", "Forest_Managment_For_Biodiv",
               "Waterbill", "ASCsq")

stopifnot(all(c(id_col, panel_col, vapply(1:8, answer_col, "")) %in% names(dep)),
          all(c("version", "CS", "alt", attr_cols) %in% names(design)),
          nrow(design) == 3 * 8 * 3, anyDuplicated(kept$response_id) == 0)

dep <- dep[match(kept$response_id, dep[[id_col]]), ]
if (anyNA(dep[[id_col]])) stop("kept respondents missing from the deposit")

# chosen alternative per respondent and card
chosen <- sapply(1:8, function(cs) answers[dep[[answer_col(cs)]]])
if (anyNA(chosen)) stop("unrecognised or missing choice answers in kept respondents")

long <- expand.grid(alt = 1:3, CS = 1:8, r = seq_len(nrow(dep)))
long$idIDdelaréponse <- dep[[id_col]][long$r]
long$ID      <- dep[[panel_col]][long$r]
long$version <- kept$version[long$r]
long$choice  <- as.integer(chosen[cbind(long$r, long$CS)] == long$alt)
long <- merge(long, design, by = c("version", "CS", "alt"), sort = FALSE)
long <- long[order(long$r, long$CS, long$alt),
             c("ID", "idIDdelaréponse", "CS", "alt", "choice", attr_cols)]
stopifnot(nrow(long) == nrow(dep) * 8 * 3)

write.csv(long, out_f, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("wrote %s: %d respondents, %d rows\n", out_f, nrow(dep), nrow(long)))

if ("--check" %in% commandArgs(trailingOnly = TRUE)) {
  orig <- read.csv(file.path(dir_in, "data.csv"), stringsAsFactors = FALSE)
  key  <- c("ID", "CS", "alt")
  cmp_cols <- c("idIDdelaréponse", "choice", attr_cols)
  orig <- orig[do.call(order, orig[key]), c(key, cmp_cols)]
  new  <- long[do.call(order, long[key]), c(key, cmp_cols)]
  num  <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
  for (col in c("CS", "alt", cmp_cols)) { orig[[col]] <- num(orig[[col]]); new[[col]] <- num(new[[col]]) }
  rownames(orig) <- rownames(new) <- NULL
  same <- isTRUE(all.equal(orig, new, check.attributes = FALSE))
  if (!same) stop("rebuilt file differs from data.csv on the analysis columns")
  cat("check: identical to data.csv on", length(cmp_cols) + 3, "analysis columns\n")
}
