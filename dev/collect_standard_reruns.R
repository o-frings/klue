#!/usr/bin/env Rscript
# =============================================================================
# dev/collect_standard_reruns.R
#
# Purpose. Collect every manuscript number that the klue 0.10.0 standard-Apollo
# reruns (bash dev/rerun_standard_apollo.sh) replace, from the new result
# files, next to the value published in version12.tex, so the revision can be
# written from one file. Covers
#   1. tab:emp_summary and the appendix MMNL rows: BIC-best LCMNL per dataset,
#      the four MMNL benchmarks M1-M4 (BIC, LL, k, convergence, Apollo code,
#      nested_ok), the best MMNL with and without random constants, BIC
#      margins and the verdict under fixed constants only (M1, M2) and under
#      all four benchmarks;
#   2. tab:mmnl_blocked: winner counts (MNL / LCMNL / MMNL) and LCMNL
#      over-selection, K* = 1 and K* >= 2;
#   3. tab:stress_blocked, tab:stress (and fig:time): mean gap to the best LL
#      and the "stuck" count per rung x arm, mean seconds per arm;
#   4. tab:mmnl, tab:mmnl_corr: winner counts by DGP type;
#   5. output/h4_misspec_blocked.csv (an arm added after version12).
#
# Definitions, as in version12.tex and the drivers:
#   BIC        -2 LL + k log(N), N = respondents (klue's $BIC; Apollo's
#              observation-based BIC is $BIC_apollo and is not used).
#   converged  Apollo's successfulEstimation (klue_mmnl()$converged); an MMNL
#              that did not converge does not enter a comparison.
#   winner     lowest BIC (the simulation drivers store it as `winner` or
#              `bic_prefers`; a failed MMNL carries BIC = Inf there).
#   verdict    "discrete" if the BIC-best LCMNL has a lower BIC than every
#              converged MMNL benchmark in the set, else "continuous".
#   margin     best MMNL BIC - BIC-best LCMNL BIC (> 0 favours LCMNL).
#   gap, stuck LL - best LL found by any arm in the cell; stuck = gap < -0.5.
#
# Partial reruns. Nothing here fails on a missing or half-finished result.
# Missing fits, cells, arms and files are reported as "pending"; a fit, cell
# or step that has started is reported from the lane log; a verdict computed
# while a benchmark is still missing is "provisional" if it reads discrete
# (a further benchmark can only lower the best MMNL BIC, so a continuous
# reading is already final); a result file last written before the rerun was
# launched (2026-09-29 23:33, output/rerun_std_launch.log) is reported as the
# superseded original, not as a new number.
#
# Inputs (read-only; paths relative to the repository root):
#   version12.tex                                     published values
#   output/rerun_std_lane.log, the driver logs it names   progress
#   output/mmnl_bench_std/<dataset>_M{1..4}.rds       new MMNL benchmark fits
#   output/vittel_respec_lcmnl_ext.csv                Vittel LCMNL, identified coding
#   output/emp_cmax_extend_{Vittel,Electricity,Swissmetro}.csv,
#   output/mode_lcmnl_results.csv, output/swiss_lcmnl_results.csv   LCMNL ladders
#   output/k1_mmnl_blocked_full.csv, output/discrete_mmnl_blocked.csv,
#   output/h4_misspec_blocked.csv                     blocked H4 arms
#   output/stress_replicate_blocked.rds, output/stress_replicate.rds   ladders
#   output/mmnl_results.csv, output/mmnl_correlated_results.csv        tab:mmnl(_corr)
#   output/superseded_klue095/*                       the published originals
# Outputs:
#   <out>/standard_rerun_summary.csv   long format, one number per row:
#       section, table, row, column, metric, new, published, status,
#       new_source, published_source, note
#   <out>/standard_rerun_summary.md    the same as readable tables
#   <out> defaults to output/.
#
# Run from the repository root. Base R only: no estimation, no RNG, no
# package loading, a few MB of memory, so it is safe beside a running lane.
#   Rscript dev/collect_standard_reruns.R               # writes output/standard_rerun_summary.*
#   Rscript dev/collect_standard_reruns.R --out DIR     # writes DIR/standard_rerun_summary.*
#   KLUE_COLLECT_OUT=DIR Rscript dev/collect_standard_reruns.R
# Deterministic: the outputs depend only on the input files (no clock time is
# written; input modification times are listed).
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)

# ---- arguments and paths ------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
OUT_DIR <- Sys.getenv("KLUE_COLLECT_OUT", "output")
for (i in seq_along(args)) {
  if (identical(args[i], "--out")) {
    if (i == length(args)) stop("--out needs a directory")
    OUT_DIR <- args[i + 1L]
  } else if (startsWith(args[i], "--out=")) OUT_DIR <- sub("^--out=", "", args[i])
}
if (!file.exists("version12.tex") || !dir.exists("output"))
  stop("run from the repository root (version12.tex and output/ not found here)")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TEX    <- "version12.tex"
SUP    <- "output/superseded_klue095"
BENCH  <- "output/mmnl_bench_std"
LANE   <- "output/rerun_std_lane.log"   # the default single lane; with --two-lanes
                                        # the step states read 'not started' (see
                                        # output/rerun_std_lane{A,B}.log instead)
# Rerun launch: output/rerun_std_launch.log, first START in output/rerun_std_laneA.log.
CUTOFF <- as.POSIXct("2026-09-29 23:33:00", tz = "Europe/Paris")

# ---- small helpers ------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a
is_new  <- function(f) file.exists(f) && file.mtime(f) > CUTOFF
mtime   <- function(f) if (file.exists(f))
  format(file.mtime(f), "%Y-%m-%d %H:%M", tz = "Europe/Paris") else "missing"
file_state <- function(f) if (!file.exists(f)) "missing" else
  if (is_new(f)) "new" else "superseded original"

# one value -> string for the CSV (numbers to 4 decimals, trailing zeros cut)
vstr <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  x <- x[[1]]
  if (is.logical(x)) return(if (is.na(x)) NA_character_ else as.character(x))
  if (is.numeric(x)) {
    if (is.na(x)) return(NA_character_)
    if (!is.finite(x)) return(as.character(x))
    return(sub("\\.$", "", sub("0+$", "", sprintf("%.4f", x))))
  }
  as.character(x)
}
# numbers -> string for the markdown tables
fm <- function(x, d = 1) {
  if (is.null(x) || !length(x)) return("")
  vapply(as.list(x), function(v) {
    if (is.null(v) || !length(v) || is.na(v)) return("")
    if (!is.numeric(v)) return(as.character(v))
    if (!is.finite(v)) return(as.character(v))
    formatC(round(v, d) + 0, format = "f", digits = d, big.mark = ",")   # + 0 drops "-0.00"
  }, character(1))
}

ROWS <- list()
add <- function(section, table, row, column, metric, new = NA, published = NA,
                status = "", new_source = "", published_source = "", note = "") {
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    section = section, table = table, row = row, column = column, metric = metric,
    new = vstr(new), published = vstr(published), status = status,
    new_source = new_source, published_source = published_source, note = note)
}
MD <- character()
md <- function(...) MD <<- c(MD, ...)
md_table <- function(df) {
  df[] <- lapply(df, function(v) { v <- as.character(v); v[is.na(v)] <- ""; v })
  c(paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")), "")
}

# ---- version12.tex: table rows --------------------------------------------------
TEXL <- readLines(TEX, warn = FALSE)
# Rows of the table carrying \label{<label>}: list(line, cells) per row
# ending in "\\", from the label to the table's \bottomrule.
tex_rows <- function(label) {
  i <- grep(paste0("\\label{", label, "}"), TEXL, fixed = TRUE)[1]
  if (is.na(i)) return(list())
  j <- i + which(grepl("\\bottomrule", TEXL[(i + 1L):length(TEXL)], fixed = TRUE))[1]
  if (is.na(j)) return(list())
  out <- list()
  for (l in (i + 1L):(j - 1L)) {
    s <- TEXL[l]
    if (!grepl("&", s, fixed = TRUE) || !grepl("\\\\\\\\\\s*$", s)) next
    cells <- trimws(strsplit(sub("\\\\\\\\\\s*$", "", s), "&", fixed = TRUE)[[1]])
    out[[length(out) + 1L]] <- list(line = l, cells = cells)
  }
  out
}
# all numbers in a LaTeX cell ("$12{,}195$ ($C{=}7$)" -> 12195, 7)
tex_nums <- function(s) {
  s <- gsub("{,}", "", s, fixed = TRUE)
  as.numeric(regmatches(s, gregexpr("-?[0-9]+(\\.[0-9]+)?", s))[[1]])
}
tex_ref <- function(line) sprintf("%s:%d", TEX, line)
# first row of a table whose first cell matches `pattern`
tex_find <- function(rows, pattern) {
  for (r in rows) if (grepl(pattern, r$cells[1])) return(r)
  NULL
}

# ---- lane progress --------------------------------------------------------------
PLAN <- c("dev/rerun_mmnl_benchmark.R Vittel",
          "dev/rerun_mmnl_benchmark.R Swissmetro Mode SwissRoute Electricity",
          "dev/run_mmnl_second_start.R",
          "dev/rerun_mmnl_benchmark.R summary",
          "dev/run_k1_mmnl_blocked.R",
          "dev/run_discrete_mmnl_blocked.R",
          "dev/run_h4_misspec_blocked.R",
          "dev/stress_replicate_blocked.R",
          "dev/stress_replicate.R",
          "dev/run_mmnl_studies.R mmnl mmnl_corr")      # dev/rerun_standard_apollo.sh, laneA + laneB
lane <- data.frame(step = PLAN, state = "not started", started = "", ended = "",
                   exit = NA_integer_, log = "")
if (file.exists(LANE)) {
  open <- NA_integer_
  for (l in readLines(LANE, warn = FALSE)) {
    ts <- sub("^\\[([^]]+)\\].*$", "\\1", l)
    if (grepl("^\\[[^]]+\\] START ", l)) {
      cmd <- sub("^\\[[^]]+\\] START (.*) \\(log: (.*)\\)$", "\\1", l)
      lg  <- sub("^\\[[^]]+\\] START (.*) \\(log: (.*)\\)$", "\\2", l)
      open <- match(cmd, lane$step)
      if (!is.na(open)) {
        lane$state[open] <- "started, no EXIT yet"; lane$started[open] <- ts
        lane$log[open] <- lg; lane$ended[open] <- ""; lane$exit[open] <- NA
      }
    } else if (grepl("^\\[[^]]+\\] EXIT ", l) && !is.na(open)) {
      code <- suppressWarnings(as.integer(sub("^\\[[^]]+\\] EXIT (-?[0-9]+) .*$", "\\1", l)))
      lane$exit[open]  <- code; lane$ended[open] <- ts
      lane$state[open] <- if (isTRUE(code == 0L)) "done" else sprintf("exited %s", code)
      open <- NA_integer_
    }
  }
}
running_step <- lane$step[lane$state == "started, no EXIT yet"]
for (k in seq_len(nrow(lane)))
  add("lane", "rerun_standard_apollo.sh", lane$step[k], "", "state", lane$state[k],
      status = lane$state[k], new_source = LANE,
      note = paste0(if (nzchar(lane$started[k])) paste0("started ", lane$started[k]) else "",
                    if (nzchar(lane$ended[k])) paste0(", ended ", lane$ended[k]) else "",
                    if (nzchar(lane$log[k])) paste0("; log ", lane$log[k]) else ""))

# MMNL benchmark fit currently running: START without DONE in the driver log
# of a benchmark step that has started but not exited.
bench_running <- character()
for (k in which(lane$state == "started, no EXIT yet" &
                startsWith(lane$step, "dev/rerun_mmnl_benchmark.R"))) {
  if (!file.exists(lane$log[k])) next
  L  <- readLines(lane$log[k], warn = FALSE)
  st <- sub("^.*START (\\S+) (M[1-4]) .*$", "\\1_\\2",
            grep("\\] START \\S+ M[1-4] ", L, value = TRUE, perl = TRUE), perl = TRUE)
  dn <- sub("^.*DONE  (\\S+) (M[1-4]) .*$", "\\1_\\2",
            grep("\\] DONE  \\S+ M[1-4] ", L, value = TRUE, perl = TRUE), perl = TRUE)
  bench_running <- c(bench_running, setdiff(st, dn))
}

# =============================================================================
# 1. Empirical: tab:emp_summary (+ appendix MMNL rows)
# =============================================================================
DATASETS <- c("Vittel", "Mode", "SwissRoute", "Electricity", "Swissmetro")   # manuscript order
TEX_ROW  <- c(Vittel = "^Vittel", Mode = "^Apollo mode", SwissRoute = "^Apollo Swiss route",
              Electricity = "^Electricity", Swissmetro = "^Swissmetro")
TEX_APP  <- c(Vittel = "tab:emp_mmnl", Mode = "tab:mode_mmnl", SwissRoute = "tab:swiss_mmnl",
              Electricity = "tab:elec_mmnl", Swissmetro = "tab:sm_mmnl")
# LCMNL ladders. The LCMNL side is not Apollo-dependent and was not rerun,
# except Vittel, whose benchmark uses the identified coding (README_REPRODUCE.md,
# "Vittel, identified specification"); version12 used the rank-deficient one.
LC_NEW <- list(Vittel      = c("output/vittel_respec_lcmnl_ext.csv", "output/vittel_respec_lcmnl.csv"),
               Mode        = "output/mode_lcmnl_results.csv",
               SwissRoute  = "output/swiss_lcmnl_results.csv",
               Electricity = "output/emp_cmax_extend_Electricity.csv",
               Swissmetro  = "output/emp_cmax_extend_Swissmetro.csv")
LC_PUB <- list(Vittel = "output/emp_cmax_extend_Vittel.csv", Mode = LC_NEW$Mode,
               SwissRoute = LC_NEW$SwissRoute, Electricity = LC_NEW$Electricity,
               Swissmetro = LC_NEW$Swissmetro)
MODELS <- c(M1 = "independent, fixed constants", M2 = "correlated, fixed constants",
            M3 = "independent, random constants", M4 = "correlated, random constants")
PREREQ <- c(M2 = "M1", M4 = "M3")   # dev/rerun_mmnl_benchmark.R:96-100

lc_best <- function(files) {
  f <- files[file.exists(files)][1]
  if (is.na(f)) return(NULL)
  d <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  ok <- d$converged %in% TRUE & is.finite(d$BIC)
  if (!any(ok)) return(NULL)
  d  <- d[ok, ]; d <- d[order(d$BIC), ]
  b  <- d[1, ]
  list(file = f, C = b$C, BIC = b$BIC, LL = b$LL, k = b$k,
       singular = if ("singular" %in% names(b)) isTRUE(b$singular) else NA,
       Cmax = max(d$C), boundary = b$C == max(d$C),
       runner_C = if (nrow(d) > 1) d$C[2] else NA, runner_dBIC = if (nrow(d) > 1) d$BIC[2] - b$BIC else NA)
}
read_fit <- function(ds, m) {
  f <- file.path(BENCH, sprintf("%s_%s.rds", ds, m))
  if (!file.exists(f)) return(NULL)
  o <- tryCatch(readRDS(f), error = function(e) NULL)   # e.g. caught mid-write
  if (is.null(o)) return(list(file = f, unreadable = TRUE))
  st <- o$apollo_status %||% list()
  cv <- isTRUE(o$converged)   # a failed fit carries LL = -Inf, BIC = Inf, k = 0: report NA
  list(file = f, unreadable = FALSE, converged = cv, LL = if (cv) o$LL else NA,
       k = if (cv) o$k else NA, BIC = if (cv) o$BIC else NA, BIC_apollo = o$BIC_apollo %||% NA,
       nested_ok = if (is.null(o$nested_ok)) NA else o$nested_ok,
       code = st$code %||% NA, message = st$message %||% NA, nIter = st$nIter %||% NA,
       minutes = (o$secs %||% NA) / 60, vcov = !is.null(o$vcov) && !anyNA(o$vcov),
       N_implied = if (isTRUE(o$converged) && o$k > 0) exp((o$BIC + 2 * o$LL) / o$k) else NA,
       start_used = o$start_used %||% NA, reason = o$reason %||% NA,
       old = old_fit_check(ds, m, o))
}
# The superseded klue 0.9.5 fit of the same model (output/mmnl_bench/): with the
# same k, draws and N the two models are the same likelihood (0.9.5 used
# exp(sigma), which spans the same models), so a new LL more than 0.5 below the
# old one means the new single estimation stopped at a lower local maximum.
old_fit_check <- function(ds, m, o) {
  f <- file.path("output/mmnl_bench", sprintf("%s_%s.rds", ds, m))
  if (!file.exists(f) || !isTRUE(o$converged)) return(list(flag = NA, LL = NA))
  p <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(p) || !isTRUE(p$converged)) return(list(flag = NA, LL = NA))
  same <- identical(as.numeric(p$k), as.numeric(o$k)) &&
          isTRUE(all.equal(as.numeric(p$settings$n_draws), as.numeric(o$settings$n_draws)))
  list(flag = same && p$LL > o$LL + 0.5, LL = p$LL, same_spec = same)
}

emp_summary <- tex_rows("tab:emp_summary")
emp <- list()
for (ds in DATASETS) {
  # --- published (version12.tex) ---
  r   <- tex_find(emp_summary, TEX_ROW[[ds]])
  pub <- list(line = NA, lc_BIC = NA, lc_C = NA, ind = NA, corr = NA)
  if (!is.null(r) && length(r$cells) >= 5) {
    n3 <- tex_nums(r$cells[3])
    pub <- list(line = r$line, lc_BIC = n3[1], lc_C = n3[2],
                ind = tex_nums(r$cells[4])[1], corr = tex_nums(r$cells[5])[1])
  }
  app  <- tex_rows(TEX_APP[[ds]])
  pa   <- list(M1 = tex_find(app, "^MMNL \\(indep"), M2 = tex_find(app, "^MMNL \\(correlated"))
  pubm <- lapply(pa, function(x) if (is.null(x)) NULL else
    list(line = x$line, k = tex_nums(x$cells[2])[1], LL = tex_nums(x$cells[3])[1],
         BIC = tex_nums(x$cells[4])[1]))
  # --- new ---
  lc   <- lc_best(LC_NEW[[ds]])
  lcp  <- lc_best(LC_PUB[[ds]])
  fits <- lapply(names(MODELS), function(m) read_fit(ds, m)); names(fits) <- names(MODELS)
  mstat <- vapply(names(MODELS), function(m) {
    x <- fits[[m]]
    if (!is.null(x) && isTRUE(x$unreadable)) return("pending (file unreadable, being written?)")
    if (!is.null(x)) return(if (x$converged) "done" else "done, not converged")
    p <- PREREQ[m]
    if (!is.na(p) && !is.null(fits[[p]]) && !isTRUE(fits[[p]]$unreadable) && !fits[[p]]$converged)
      return(sprintf("not run (%s not converged)", p))
    if (paste0(ds, "_", m) %in% bench_running) return("running")
    "pending"
  }, character(1))
  n_done   <- sum(startsWith(mstat, "done"))
  n_expect <- sum(!startsWith(mstat, "not run"))
  ds_status <- if (n_done == 0) "pending" else if (n_done == n_expect) "done" else "partial"
  conv_bic <- vapply(names(MODELS), function(m) {
    x <- fits[[m]]; if (is.null(x) || isTRUE(x$unreadable) || !x$converged) NA_real_ else x$BIC
  }, numeric(1))
  waiting <- names(MODELS)[mstat %in% c("pending", "running") | startsWith(mstat, "pending")]

  judge <- function(set) {
    b <- conv_bic[set]; miss <- intersect(set, waiting)
    if (all(is.na(b)) || is.null(lc))
      return(list(best = NA, BIC = NA, margin = NA,
                  verdict = if (length(miss)) "pending" else "no converged MMNL"))
    best <- names(which.min(b)); margin <- min(b, na.rm = TRUE) - lc$BIC
    v <- if (lc$BIC < min(b, na.rm = TRUE)) "discrete" else "continuous"
    if (length(miss)) v <- if (v == "discrete")
      sprintf("discrete (provisional: %s pending)", paste(miss, collapse = ", ")) else
      sprintf("continuous (final although %s pending)", paste(miss, collapse = ", "))
    list(best = best, BIC = min(b, na.rm = TRUE), margin = margin, verdict = v)
  }
  jf <- judge(c("M1", "M2")); jr <- judge(c("M3", "M4")); ja <- judge(names(MODELS))
  pub_best <- suppressWarnings(min(c(pub$ind, pub$corr), na.rm = TRUE))
  pub_verdict <- if (is.finite(pub_best) && !is.na(pub$lc_BIC))
    (if (pub$lc_BIC < pub_best) "discrete" else "continuous") else NA
  emp[[ds]] <- list(pub = pub, pubm = pubm, lc = lc, lcp = lcp, fits = fits, mstat = mstat,
                    status = ds_status, jf = jf, jr = jr, ja = ja,
                    pub_best = pub_best, pub_verdict = pub_verdict)

  # --- CSV rows ---
  SEC <- "empirical"; TAB <- "tab:emp_summary"
  src_pub <- if (is.na(pub$line)) "version12.tex (row not found)" else tex_ref(pub$line)
  add(SEC, TAB, ds, "", "dataset_status", ds_status, status = ds_status, new_source = BENCH,
      note = paste(sprintf("%s %s", names(MODELS), mstat), collapse = "; "))
  if (!is.null(lc)) {
    add(SEC, TAB, ds, "LCMNL", "best_C", lc$C, pub$lc_C, "done", lc$file, src_pub,
        note = sprintf("C searched 1..%d%s; runner-up C=%s is %s BIC behind%s", lc$Cmax,
                       if (lc$boundary) ", minimum at the boundary" else "",
                       vstr(lc$runner_C), vstr(round(lc$runner_dBIC, 1)),
                       if (isTRUE(lc$singular)) "; best fit flagged singular" else ""))
    add(SEC, TAB, ds, "LCMNL", "best_BIC", lc$BIC, pub$lc_BIC, "done", lc$file, src_pub,
        note = if (!is.null(lcp) && lcp$file != lc$file)
          sprintf("published from %s (C=%d, BIC %.1f)", lcp$file, lcp$C, lcp$BIC) else "")
    add(SEC, TAB, ds, "LCMNL", "best_LL", lc$LL, if (!is.null(lcp)) lcp$LL else NA, "done", lc$file,
        if (!is.null(lcp)) lcp$file else "")
    add(SEC, TAB, ds, "LCMNL", "best_k", lc$k, if (!is.null(lcp)) lcp$k else NA, "done", lc$file,
        if (!is.null(lcp)) lcp$file else "")
  } else add(SEC, TAB, ds, "LCMNL", "best_BIC", NA, pub$lc_BIC, "pending", paste(LC_NEW[[ds]], collapse = " | "), src_pub)
  for (m in names(MODELS)) {
    x  <- fits[[m]]; pm <- pubm[[m]]
    pref <- if (!is.null(pm)) pm$BIC else if (m == "M1") pub$ind else if (m == "M2") pub$corr else NA
    psrc <- if (!is.null(pm)) tex_ref(pm$line) else if (!is.na(pref)) src_pub else "not in version12"
    nsrc <- file.path(BENCH, sprintf("%s_%s.rds", ds, m))
    add(SEC, TAB, ds, m, "status", mstat[[m]], status = mstat[[m]], new_source = nsrc, note = MODELS[[m]])
    ok <- !is.null(x) && !isTRUE(x$unreadable)
    add(SEC, TAB, ds, m, "BIC", if (ok) x$BIC else NA, pref, mstat[[m]], nsrc, psrc)
    add(SEC, TAB, ds, m, "LL", if (ok) x$LL else NA, if (!is.null(pm)) pm$LL else NA, mstat[[m]], nsrc, psrc)
    add(SEC, TAB, ds, m, "k", if (ok) x$k else NA, if (!is.null(pm)) pm$k else NA, mstat[[m]], nsrc, psrc)
    if (ok) {
      add(SEC, TAB, ds, m, "converged", x$converged, status = mstat[[m]], new_source = nsrc,
          note = "Apollo successfulEstimation")
      add(SEC, TAB, ds, m, "apollo_code", x$code, status = mstat[[m]], new_source = nsrc, note = x$message)
      add(SEC, TAB, ds, m, "nested_ok", x$nested_ok, status = mstat[[m]], new_source = nsrc,
          note = if (m %in% names(PREREQ)) sprintf("LL >= %s LL (same draws)", PREREQ[[m]]) else "not applicable")
      add(SEC, TAB, ds, m, "vcov_available", x$vcov, status = mstat[[m]], new_source = nsrc)
      add(SEC, TAB, ds, m, "minutes", x$minutes, status = mstat[[m]], new_source = nsrc)
      add(SEC, TAB, ds, m, "N_implied_by_BIC", x$N_implied, status = mstat[[m]], new_source = nsrc,
          note = "exp((BIC + 2 LL) / k): respondents if BIC is on respondents")
      if (isTRUE(x$old$same_spec))
        add(SEC, TAB, ds, m, "LL_superseded_0.9.5_fit", x$old$LL, status = mstat[[m]],
            new_source = file.path("output/mmnl_bench", sprintf("%s_%s.rds", ds, m)),
            note = if (isTRUE(x$old$flag))
              "WARNING: the 0.9.5 fit of the same model is > 0.5 LL higher; the new fit stopped at a lower local maximum"
            else "new fit within 0.5 LL of the 0.9.5 fit or better")
    }
  }
  add(SEC, TAB, ds, "MMNL", "correlation_gain_fixed_M2_minus_M1",
      conv_bic["M2"] - conv_bic["M1"], pub$corr - pub$ind,
      if (all(!is.na(conv_bic[c("M1", "M2")]))) "done" else "pending", BENCH, src_pub,
      "BIC; < 0 means the correlated model is better")
  add(SEC, TAB, ds, "MMNL", "correlation_gain_random_M4_minus_M3",
      conv_bic["M4"] - conv_bic["M3"], NA,
      if (all(!is.na(conv_bic[c("M3", "M4")]))) "done" else "pending", BENCH, "not in version12")
  for (j in list(list("fixed_constants_M1_M2", jf, pub_best, pub_verdict),
                 list("random_constants_M3_M4", jr, NA, NA),
                 list("all_four_M1_M4", ja, NA, NA))) {
    jj <- j[[2]]
    st <- if (startsWith(jj$verdict, "pending")) "pending" else
      if (grepl("pending", jj$verdict)) "partial" else "done"
    add(SEC, TAB, ds, j[[1]], "best_MMNL", jj$best, if (j[[1]] == "fixed_constants_M1_M2")
      (if (is.finite(pub_best)) (if (pub_best == pub$ind) "MMNL_ind" else "MMNL_corr") else NA) else NA,
      st, BENCH, if (j[[1]] == "fixed_constants_M1_M2") src_pub else "not in version12")
    add(SEC, TAB, ds, j[[1]], "best_MMNL_BIC", jj$BIC, j[[3]], st, BENCH,
        if (j[[1]] == "fixed_constants_M1_M2") src_pub else "not in version12")
    if (j[[1]] != "random_constants_M3_M4") {
      add(SEC, TAB, ds, j[[1]], "margin_bestMMNL_minus_LCMNL", jj$margin,
          if (is.finite(j[[3]] %||% NA)) j[[3]] - pub$lc_BIC else NA, st, BENCH,
          if (j[[1]] == "fixed_constants_M1_M2") src_pub else "not in version12",
          "> 0: LCMNL has the lower BIC")
      add(SEC, TAB, ds, j[[1]], "verdict", jj$verdict, j[[4]], st, BENCH,
          if (j[[1]] == "fixed_constants_M1_M2") src_pub else "not in version12",
          if (j[[1]] == "fixed_constants_M1_M2")
            "version12's benchmarks had fixed constants (MMNL_ind ~ M1, MMNL_corr ~ M2)" else "")
    }
  }
}

# =============================================================================
# 2. tab:mmnl_blocked (K* = 1 and K* >= 2 arms)
# =============================================================================
H4 <- list(
  K1   = list(file = "output/k1_mmnl_blocked_full.csv", key = c("sigma", "rep"), expected = 80L,
              label = "K* = 1 (continuous)", tex = "K\\^\\* = 1", driver = "dev/run_k1_mmnl_blocked.R"),
  Kge2 = list(file = "output/discrete_mmnl_blocked.csv", key = c("true_K", "kappa", "sigma", "rep"),
              expected = 36L, label = "K* >= 2 (discrete)", tex = "geq 2", driver = "dev/run_discrete_mmnl_blocked.R"))
WIN <- c("MNL", "LCMNL", "MMNL")
h4_read <- function(f) {
  if (!file.exists(f)) return(NULL)
  d <- tryCatch(read.csv(f), error = function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d$winner <- sub("\\..*$", "", d$winner)          # superseded files store "MMNL.model"
  d
}
h4_summ <- function(d) {
  if (is.null(d)) return(NULL)
  tb <- table(factor(d$winner, levels = WIN))
  list(n = nrow(d), MNL = tb[["MNL"]], LCMNL = tb[["LCMNL"]], MMNL = tb[["MMNL"]],
       over = sum(d$lc_bic < d$mnl_bic), mmnl_fail = sum(!is.finite(d$mmnl_bic)),
       lc_C = paste(sprintf("C=%s:%d", names(table(d$lc_C[d$lc_bic < d$mnl_bic])),
                            as.integer(table(d$lc_C[d$lc_bic < d$mnl_bic]))), collapse = " "))
}
mb_tex <- tex_rows("tab:mmnl_blocked")
# over-selection count published in the table notes (version12.tex:428)
over_line <- grep("over-selects to $C \\geq 2$ in $", TEXL, fixed = TRUE)[1]
over_pub  <- if (is.na(over_line)) NA else
  as.numeric(sub("^.*over-selects to \\$C \\\\geq 2\\$ in \\$([0-9]+)\\$ of the.*$", "\\1", TEXL[over_line]))
h4res <- list()
for (a in names(H4)) {
  h  <- H4[[a]]
  st_file <- file_state(h$file)
  dn <- if (st_file == "new") h4_read(h$file) else NULL
  dp <- h4_read(file.path(SUP, basename(h$file)))
  sn <- h4_summ(dn); sp <- h4_summ(dp)
  # published counts on the same cells as the new file (paired)
  spp <- NULL; changed <- NA
  if (!is.null(dn) && !is.null(dp)) {
    kn <- do.call(paste, dn[h$key]); kp <- do.call(paste, dp[h$key])
    dpp <- dp[kp %in% kn, ]; spp <- h4_summ(dpp)
    m <- match(kn, kp); changed <- sum(dn$winner != dp$winner[m], na.rm = TRUE)
  }
  r  <- tex_find(mb_tex, h$tex)
  pt <- if (is.null(r)) NULL else list(line = r$line, n = tail(tex_nums(r$cells[1]), 1),
                                       MNL = tex_nums(r$cells[2])[1], LCMNL = tex_nums(r$cells[3])[1],
                                       MMNL = tex_nums(r$cells[4])[1])
  n_new  <- if (is.null(sn)) 0L else sn$n
  status <- if (st_file == "missing") "pending (file not yet written)" else
    if (st_file != "new") "pending (superseded original still in place)" else
    if (n_new >= h$expected) "done" else sprintf("partial (%d/%d cells)", n_new, h$expected)
  if (any(lane$state[lane$step == h$driver] == "started, no EXIT yet")) status <- paste0(status, "; driver running")
  h4res[[a]] <- list(h = h, sn = sn, sp = sp, spp = spp, pt = pt, status = status,
                     changed = changed, st_file = st_file)
  SEC <- "simulation"; TAB <- "tab:mmnl_blocked"
  psrc <- if (is.null(pt)) "version12.tex (row not found)" else tex_ref(pt$line)
  add(SEC, TAB, h$label, "", "cells_done", n_new, if (is.null(pt)) NA else pt$n, status, h$file, psrc,
      sprintf("expected %d; file %s (mtime %s)", h$expected, st_file, mtime(h$file)))
  for (w in WIN)
    add(SEC, TAB, h$label, w, "winner_count", if (is.null(sn)) NA else sn[[w]],
        if (is.null(pt)) NA else pt[[w]], status, h$file, psrc,
        if (!is.null(spp)) sprintf("published on the same %d cells: %d (%s)", spp$n, spp[[w]],
                                   file.path(SUP, basename(h$file))) else "")
  if (a == "Kge2") {
    cc <- function(d) if (is.null(d)) NA else sum(d$lc_C == d$true_K)
    add(SEC, TAB, h$label, "LCMNL", "best_C_equals_true_K", cc(dn), cc(dp), status, h$file,
        file.path(SUP, basename(h$file)), "best LCMNL over C = 2..5")
  } else
  add(SEC, TAB, h$label, "LCMNL", "overselects_best_C2plus_below_MNL", if (is.null(sn)) NA else sn$over,
      if (a == "K1") (if (!is.na(over_pub)) over_pub else if (!is.null(sp)) sp$over else NA) else
        (if (!is.null(sp)) sp$over else NA), status, h$file,
      if (a == "K1" && !is.na(over_line)) tex_ref(over_line) else file.path(SUP, basename(h$file)),
      paste0(if (!is.null(sn)) paste0("new by C: ", sn$lc_C) else "",
             if (!is.null(spp)) sprintf("; published on the same cells: %d", spp$over) else ""))
  add(SEC, TAB, h$label, "MMNL", "mmnl_not_converged", if (is.null(sn)) NA else sn$mmnl_fail,
      if (!is.null(sp)) sp$mmnl_fail else NA, status, h$file, file.path(SUP, basename(h$file)),
      "BIC = Inf: MMNL excluded from the cell's comparison")
  add(SEC, TAB, h$label, "", "winner_changed_vs_published_same_cells", changed, NA, status, h$file,
      file.path(SUP, basename(h$file)))
  if (!is.null(sp) && !is.null(pt))
    add(SEC, TAB, h$label, "", "published_check_superseded_file_vs_tex",
        all(c(sp$n, sp$MNL, sp$LCMNL, sp$MMNL) == c(pt$n, pt$MNL, pt$LCMNL, pt$MMNL)),
        status = "check", new_source = file.path(SUP, basename(h$file)), published_source = psrc,
        note = sprintf("superseded file: n=%d MNL=%d LCMNL=%d MMNL=%d", sp$n, sp$MNL, sp$LCMNL, sp$MMNL))
}
# overall row
ov_new <- if (all(vapply(h4res, function(x) !is.null(x$sn), logical(1))))
  lapply(setNames(c("n", WIN), c("n", WIN)), function(w) sum(vapply(h4res, function(x) x$sn[[w]], numeric(1)))) else NULL
r <- tex_find(mb_tex, "^Overall")
ov_pub <- if (is.null(r)) NULL else list(line = r$line, n = tail(tex_nums(r$cells[1]), 1),
                                         MNL = tex_nums(r$cells[2])[1], LCMNL = tex_nums(r$cells[3])[1],
                                         MMNL = tex_nums(r$cells[4])[1])
ov_status <- if (is.null(ov_new)) "pending (an arm has no new file yet)" else
  if (all(vapply(h4res, function(x) x$status == "done", logical(1)))) "done" else "partial"
for (w in c("n", WIN))
  add("simulation", "tab:mmnl_blocked", "Overall", if (w == "n") "" else w,
      if (w == "n") "cells_done" else "winner_count",
      if (is.null(ov_new)) NA else ov_new[[w]], if (is.null(ov_pub)) NA else ov_pub[[w]], ov_status,
      paste(vapply(H4, `[[`, "", "file"), collapse = " + "),
      if (is.null(ov_pub)) "" else tex_ref(ov_pub$line))

# =============================================================================
# 3. Stress ladders: tab:stress_blocked, tab:stress (+ fig:time)
# =============================================================================
LADDERS <- list(
  blocked = list(table = "tab:stress_blocked", file = "output/stress_replicate_blocked.rds",
                 driver = "dev/stress_replicate_blocked.R", seeds = 10L,
                 rungs = c("moderate", "hard", "very_hard"),
                 arms = c("klue_ml", "klue_em", "klue_em_rp", "apollo_searchStart", "apollo_lcEM",
                          "apollo_clust", "gmnl_lc", "apollo_ss_published")),
  randomised = list(table = "tab:stress", file = "output/stress_replicate.rds",
                    driver = "dev/stress_replicate.R", seeds = 5L,
                    rungs = c("moderate", "hard", "very_hard", "hard_corr"),
                    arms = c("klue_ml", "klue_em", "apollo_searchStart", "apollo_lcEM",
                             "gmnl_lc", "apollo_ss_published")))
ARM_LABEL <- c(klue_ml = "Clustering+ML", klue_em = "Clustering+EM", klue_em_rp = "RP+EM",
               apollo_searchStart = "Apollo searchStart (defaults)", apollo_lcEM = "apollo_lcEM",
               apollo_clust = "Apollo from clustering starts", gmnl_lc = "gmnl",
               apollo_ss_published = "Apollo searchStart (published setting)")
TEX_ARM <- list(klue_ml = "Clustering\\$\\+\\$ML", klue_em = "Clustering\\$\\+\\$EM", klue_em_rp = "^RP",
                apollo_ss_published = "Apollo", gmnl_lc = "gmnl")
rung_of <- function(s) if (grepl("^Very hard", s)) "very_hard" else if (grepl("^Hard\\$\\+\\$corr", s))
  "hard_corr" else if (grepl("^Hard", s)) "hard" else if (grepl("^Moderate", s)) "moderate" else NA

# Name the arms as the drivers do after their migrate_cell(): a cell not yet
# topped up (no harness tag) holds the published non-default Apollo arm under
# "apollo_searchStart" (dev/stress_replicate_blocked.R:76-99, dev/stress_replicate.R:51-74).
migrate <- function(cell) {
  if (identical(cell$harness, "apollo_std")) return(cell)
  for (f in c("LL", "gap", "secs")) if (!is.null(cell[[f]]))
    names(cell[[f]])[names(cell[[f]]) == "apollo_searchStart"] <- "apollo_ss_published"
  cell
}
agg <- function(cells, rung, arm) {
  cs <- Filter(function(c) identical(c$rung, rung), cells)
  has <- Filter(function(c) arm %in% names(c$LL), cs)
  g <- vapply(has, function(c) {
    gg <- if (!is.null(c$gap) && arm %in% names(c$gap)) c$gap[[arm]] else c$LL[[arm]] - c$best
    as.numeric(gg) }, numeric(1))
  s <- vapply(has, function(c) if (arm %in% names(c$secs)) as.numeric(c$secs[[arm]]) else NA_real_, numeric(1))
  fin <- is.finite(g)
  list(n_cells = length(cs), n_arm = length(has), n_fin = sum(fin), n_err = sum(!fin),
       mean_gap = if (any(fin)) mean(g[fin]) else NA, stuck = sum(g[fin] < -0.5),
       secs = if (any(is.finite(s))) mean(s[is.finite(s)]) else NA)
}
stress <- list()
for (lad in names(LADDERS)) {
  L <- LADDERS[[lad]]
  cur  <- if (file.exists(L$file)) tryCatch(readRDS(L$file), error = function(e) NULL) else NULL
  pubc <- if (file.exists(file.path(SUP, basename(L$file))))
    lapply(readRDS(file.path(SUP, basename(L$file))), migrate) else NULL
  cur  <- if (is.null(cur)) list() else lapply(cur, migrate)
  newc <- Filter(function(c) identical(c$harness, "apollo_std"), cur)
  running <- any(lane$state[lane$step == L$driver] == "started, no EXIT yet")
  # published table rows from version12.tex
  tr <- tex_rows(L$table); hdr <- tr[[1]]$cells
  pub_tex <- list()
  for (r in tr[-1]) {
    rg <- rung_of(r$cells[1]); if (is.na(rg)) next
    for (a in names(TEX_ARM)) {
      col <- which(grepl(TEX_ARM[[a]], hdr))[1]
      if (is.na(col)) next
      v <- tex_nums(r$cells[col])
      pub_tex[[paste(rg, a)]] <- list(line = r$line, gap = v[1], stuck = v[2], n = v[3])
    }
  }
  res <- list(); check_bad <- character()
  for (rg in L$rungs) for (a in L$arms) {
    an <- agg(newc, rg, a)
    ap <- if (!is.null(pubc)) agg(pubc, rg, a) else NULL
    pt <- pub_tex[[paste(rg, a)]]
    status <- if (an$n_cells == 0) "pending" else
      if (an$n_arm == 0) "pending (arm not computed)" else
      if (an$n_cells < L$seeds) sprintf("partial (%d/%d seeds)", an$n_cells, L$seeds) else "done"
    if (running && status != "done") status <- paste0(status, "; driver running")
    res[[paste(rg, a)]] <- list(rung = rg, arm = a, new = an, pubfile = ap, pubtex = pt, status = status)
    if (!is.null(pt) && !is.null(ap) &&
        !(isTRUE(abs(round(ap$mean_gap, 2) - pt$gap) < 0.006) && ap$stuck == pt$stuck))
      check_bad <- c(check_bad, sprintf("%s/%s file %.2f (%d) vs tex %.2f (%d)", rg, a,
                                        ap$mean_gap, ap$stuck, pt$gap, pt$stuck))
    SEC <- "simulation"; TAB <- L$table
    nsrc <- L$file; psrc <- if (!is.null(pt)) tex_ref(pt$line) else
      if (!is.null(ap) && ap$n_arm > 0) file.path(SUP, basename(L$file)) else "not in version12"
    pub_gap <- if (!is.null(pt)) pt$gap else if (!is.null(ap) && ap$n_arm > 0) ap$mean_gap else NA
    pub_stk <- if (!is.null(pt)) pt$stuck else if (!is.null(ap) && ap$n_arm > 0) ap$stuck else NA
    add(SEC, TAB, rg, a, "mean_gap", an$mean_gap, pub_gap, status, nsrc, psrc,
        sprintf("%s; new over %d topped-up cells", ARM_LABEL[[a]], an$n_cells))
    add(SEC, TAB, rg, a, "stuck_gap_below_-0.5", an$stuck, pub_stk, status, nsrc, psrc,
        sprintf("of %d finite (%d errored)", an$n_fin, an$n_err))
    add(SEC, TAB, rg, a, "mean_seconds", an$secs, if (!is.null(ap) && ap$n_arm > 0) ap$secs else NA,
        status, nsrc, file.path(SUP, basename(L$file)),
        if (a == "apollo_ss_published") "published secs of this arm are stale (dev/stress_replicate_blocked.R:24-27)" else "")
  }
  # best LL moved by a new arm (shifts every existing gap)
  improved <- 0L; max_imp <- 0
  if (!is.null(pubc)) for (k in names(newc)) if (!is.null(pubc[[k]])) {
    d <- newc[[k]]$best - pubc[[k]]$best
    if (is.finite(d) && d > 1e-6) { improved <- improved + 1L; max_imp <- max(max_imp, d) }
  }
  add("simulation", L$table, "all rungs", "", "cells_topped_up", length(newc),
      length(L$rungs) * L$seeds, if (length(newc) == length(L$rungs) * L$seeds) "done" else
        if (length(newc)) "partial" else "pending", L$file, file.path(SUP, basename(L$file)),
      sprintf("cells with harness = apollo_std; file %s (mtime %s)", file_state(L$file), mtime(L$file)))
  add("simulation", L$table, "all rungs", "", "cells_where_new_arm_raised_best_LL", improved, NA,
      if (length(newc)) "done" else "pending", L$file, "", sprintf("largest rise %.3f LL", max_imp))
  add("simulation", L$table, "all rungs", "", "published_check_superseded_file_vs_tex",
      !length(check_bad), status = "check", new_source = file.path(SUP, basename(L$file)),
      published_source = paste0(TEX, " ", L$table),
      note = if (length(check_bad)) paste(check_bad, collapse = "; ") else "superseded file reproduces every published cell")
  stress[[lad]] <- list(L = L, res = res, n_new = length(newc), improved = improved,
                        max_imp = max_imp, check_bad = check_bad, running = running)
}

# =============================================================================
# 4. tab:mmnl and tab:mmnl_corr (randomised design)
# =============================================================================
MM <- list(
  mmnl = list(table = "tab:mmnl", file = "output/mmnl_results.csv", expected = c(K1 = 9L, Kge2 = 54L),
              cols = c("MNL", "LCMNL", "MMNL")),
  mmnl_corr = list(table = "tab:mmnl_corr", file = "output/mmnl_correlated_results.csv",
                   expected = c(K1 = 6L, Kge2 = 16L), cols = c("MNL", "LCMNL", "MMNL_indep", "MMNL_corr")))
mm_counts <- function(d, cols) {
  if (is.null(d)) return(NULL)
  w <- d$bic_prefers
  if ("lcmnl_C" %in% names(d) && "MMNL_corr" %in% cols)    # corr study: "LCMNL" includes C = 1
    w[w == "LCMNL" & d$lcmnl_C == 1] <- "MNL"
  out <- list()
  for (g in c("K1", "Kge2", "Overall")) {
    s <- if (g == "K1") d$true_K == 1 else if (g == "Kge2") d$true_K > 1 else rep(TRUE, nrow(d))
    tb <- table(factor(w[s], levels = cols))
    out[[g]] <- c(n = sum(s), setNames(as.integer(tb), cols))
  }
  out
}
mmres <- list()
for (st in names(MM)) {
  M <- MM[[st]]
  fs <- file_state(M$file)
  dn <- if (fs == "new") tryCatch(read.csv(M$file), error = function(e) NULL) else NULL
  dp <- if (file.exists(file.path(SUP, basename(M$file)))) read.csv(file.path(SUP, basename(M$file))) else
    if (fs == "superseded original") read.csv(M$file) else NULL
  cn <- mm_counts(dn, M$cols); cp <- mm_counts(dp, M$cols)
  tr <- tex_rows(M$table)
  pt <- list()
  for (g in c("K1", "Kge2", "Overall")) {
    r <- tex_find(tr, switch(g, K1 = "K\\^\\* = 1", Kge2 = "geq 2", Overall = "^Overall"))
    if (is.null(r)) next
    v <- vapply(r$cells[-1], function(s) tex_nums(s)[1], numeric(1))
    pt[[g]] <- list(line = r$line, n = tail(tex_nums(r$cells[1]), 1), v = v)
  }
  status <- if (fs == "missing") "pending (file missing)" else
    if (fs != "new") "pending (superseded original still in place)" else "done"
  if (any(lane$state[lane$step == "dev/run_mmnl_studies.R mmnl mmnl_corr"] == "started, no EXIT yet") &&
      status != "done") status <- paste0(status, "; driver running")
  mmres[[st]] <- list(M = M, cn = cn, cp = cp, pt = pt, status = status, fs = fs)
  SEC <- "simulation"; TAB <- M$table
  for (g in c("K1", "Kge2", "Overall")) {
    psrc <- if (!is.null(pt[[g]])) tex_ref(pt[[g]]$line) else ""
    exp_n <- if (g == "Overall") sum(M$expected) else M$expected[[g]]
    add(SEC, TAB, g, "", "n_conditions_kept", if (is.null(cn)) NA else cn[[g]][["n"]],
        if (is.null(pt[[g]])) NA else pt[[g]]$n, status, M$file, psrc,
        sprintf("of %d run; a condition is dropped when an MMNL fit does not converge; file %s (mtime %s)",
                exp_n, fs, mtime(M$file)))
    # published columns: tab:mmnl = MNL, LCMNL, MMNL; tab:mmnl_corr = LCMNL or MNL, MMNL indep, MMNL corr
    pubv <- if (is.null(pt[[g]])) NULL else pt[[g]]$v
    for (cc in M$cols) {
      pv <- if (is.null(pubv)) NA else if (st == "mmnl") pubv[match(cc, M$cols)] else
        switch(cc, MNL = NA, LCMNL = NA, MMNL_indep = pubv[2], MMNL_corr = pubv[3])
      add(SEC, TAB, g, cc, "winner_count", if (is.null(cn)) NA else cn[[g]][[cc]], pv, status, M$file, psrc,
          if (!is.null(cp)) sprintf("superseded file: %d", cp[[g]][[cc]]) else "")
    }
    if (st == "mmnl_corr")
      add(SEC, TAB, g, "LCMNL_or_MNL", "winner_count",
          if (is.null(cn)) NA else cn[[g]][["LCMNL"]] + cn[[g]][["MNL"]],
          if (is.null(pubv)) NA else pubv[1], status, M$file, psrc,
          "published first column; LCMNL with C = 1 is the MNL")
  }
  if (!is.null(cp) && length(pt))
    add(SEC, TAB, "all", "", "published_check_superseded_file_vs_tex",
        all(vapply(names(pt), function(g) {
          v <- if (st == "mmnl") cp[[g]][M$cols] else
            c(cp[[g]][["LCMNL"]] + cp[[g]][["MNL"]], cp[[g]][["MMNL_indep"]], cp[[g]][["MMNL_corr"]])
          isTRUE(all(v == pt[[g]]$v)) && cp[[g]][["n"]] == pt[[g]]$n }, logical(1))),
        status = "check", new_source = file.path(SUP, basename(M$file)), published_source = paste(TEX, M$table))
}

# =============================================================================
# 5. H4 misspecification arm (blocked; not in version12)
# =============================================================================
H4M_FILE <- "output/h4_misspec_blocked.csv"
h4m_state <- file_state(H4M_FILE)
h4m <- if (h4m_state == "new") tryCatch(read.csv(H4M_FILE), error = function(e) NULL) else NULL
h4m_cells <- data.frame(cell = c("skew", "skew", "overlap", "overlap"), sigma = c(0.25, 0.40, 0.30, 0.40),
                        correct = c("MNL or MMNL", "MNL or MMNL", "LCMNL", "LCMNL"))
h4m_running <- any(lane$state[lane$step == "dev/run_h4_misspec_blocked.R"] == "started, no EXIT yet")
h4m_tab <- h4m_cells
h4m_tab$n <- 0L; h4m_tab$MNL <- NA; h4m_tab$LCMNL <- NA; h4m_tab$MMNL <- NA
h4m_tab$over <- NA; h4m_tab$correct_n <- NA; h4m_tab$mmnl_fail <- NA
for (k in seq_len(nrow(h4m_cells))) {
  status <- "pending"
  if (!is.null(h4m)) {
    d <- h4m[h4m$cell == h4m_cells$cell[k] & abs(h4m$sigma - h4m_cells$sigma[k]) < 1e-9, ]
    d$winner <- sub("\\..*$", "", d$winner)
    tb <- table(factor(d$winner, levels = WIN))
    h4m_tab$n[k] <- nrow(d); h4m_tab$MNL[k] <- tb[["MNL"]]; h4m_tab$LCMNL[k] <- tb[["LCMNL"]]
    # over-selection: the BIC-best class count within the LCMNL family (C = 1
    # when the MNL beats every C >= 2) exceeds the true K* (1 for skew cells,
    # 2 for overlap cells); the driver's lcmnl_overselects column only says
    # whether some C >= 2 beats the MNL
    best_C <- ifelse(d$lc_bic < d$mnl_bic, d$lc_C, 1L)
    h4m_tab$MMNL[k] <- tb[["MMNL"]]; h4m_tab$over[k] <- sum(best_C > d$true_K)
    h4m_tab$mmnl_fail[k] <- sum(!is.finite(d$mmnl_bic))
    h4m_tab$correct_n[k] <- if (h4m_cells$cell[k] == "skew") tb[["MNL"]] + tb[["MMNL"]] else tb[["LCMNL"]]
    status <- if (nrow(d) >= 10) "done" else if (nrow(d)) sprintf("partial (%d/10)", nrow(d)) else "pending"
  } else if (h4m_state == "superseded original") status <- "pending (file predates the rerun)"
  if (h4m_running && status != "done") status <- paste0(status, "; driver running")
  h4m_tab$status[k] <- status
  lab <- sprintf("%s sigma=%.2f", h4m_cells$cell[k], h4m_cells$sigma[k])
  for (w in c("n", WIN))
    add("simulation", "h4_misspec_blocked", lab, if (w == "n") "" else w,
        if (w == "n") "replications_done" else "winner_count", h4m_tab[[w]][k], NA, status, H4M_FILE,
        "not in version12", if (w == "n") sprintf("of 10; correct reading: %s", h4m_cells$correct[k]) else "")
  add("simulation", "h4_misspec_blocked", lab, "", "correct_reading_count", h4m_tab$correct_n[k], NA,
      status, H4M_FILE, "not in version12", h4m_cells$correct[k])
  add("simulation", "h4_misspec_blocked", lab, "LCMNL", "overselects_best_C_above_true_K",
      h4m_tab$over[k], NA, status, H4M_FILE, "not in version12")
  add("simulation", "h4_misspec_blocked", lab, "MMNL", "mmnl_not_converged",
      h4m_tab$mmnl_fail[k], NA, status, H4M_FILE, "not in version12")
}

# =============================================================================
# Write CSV
# =============================================================================
csv <- do.call(rbind, ROWS)
csv_path <- file.path(OUT_DIR, "standard_rerun_summary.csv")
write.csv(csv, csv_path, row.names = FALSE, na = "")

# =============================================================================
# Write markdown
# =============================================================================
md("# Standard-Apollo rerun summary (klue 0.10.0)", "",
   "Written by `dev/collect_standard_reruns.R` from the files listed under Inputs; ",
   "each number sits next to its version12.tex value. Status words: done, partial, pending; ",
   "\"superseded original\" = a file last written before the rerun launch (2026-09-29 23:33).", "",
   "Definitions: BIC on respondents; an MMNL counts only if Apollo reports successfulEstimation; ",
   "winner = lowest BIC; margin = best MMNL BIC - BIC-best LCMNL BIC (> 0 favours LCMNL); ",
   "stuck = more than 0.5 LL below the best LL found in the cell.", "")

md("## Lane progress", "", paste0("From `", LANE, "` (step order as in `dev/rerun_standard_apollo.sh`)."), "")
md(md_table(data.frame(step = lane$step, state = lane$state, started = lane$started,
                       ended = lane$ended, log = lane$log)))

md("## 1. Empirical benchmark (tab:emp_summary, appendix MMNL tables)", "")
v1 <- do.call(rbind, lapply(DATASETS, function(ds) {
  e <- emp[[ds]]
  data.frame(dataset = ds, status = e$status,
             `LCMNL C` = if (is.null(e$lc)) "" else e$lc$C,
             `LCMNL BIC` = if (is.null(e$lc)) "" else fm(e$lc$BIC),
             `best fixed (M1/M2)` = if (is.na(e$jf$best)) "" else sprintf("%s %s", e$jf$best, fm(e$jf$BIC)),
             `best random (M3/M4)` = if (is.na(e$jr$best)) "" else sprintf("%s %s", e$jr$best, fm(e$jr$BIC)),
             `margin fixed` = fm(e$jf$margin), `margin all four` = fm(e$ja$margin),
             `verdict fixed` = e$jf$verdict, `verdict all four` = e$ja$verdict,
             check.names = FALSE)
}))
md("### New readings", "", md_table(v1))
v2 <- do.call(rbind, lapply(DATASETS, function(ds) {
  e <- emp[[ds]]
  data.frame(dataset = ds, `LCMNL BIC (C)` = sprintf("%s (C=%s)", fm(e$pub$lc_BIC, 0), e$pub$lc_C),
             `MMNL_ind` = fm(e$pub$ind, 0), `MMNL_corr` = fm(e$pub$corr, 0),
             margin = fm(e$pub_best - e$pub$lc_BIC, 0), verdict = e$pub_verdict,
             source = if (is.na(e$pub$line)) "" else tex_ref(e$pub$line), check.names = FALSE)
}))
md("### Published in version12 (tab:emp_summary)", "", md_table(v2),
   "version12's benchmarks had fixed constants: MMNL_ind corresponds to M1 and MMNL_corr to M2. ",
   "Published Vittel used the rank-deficient forest coding (k = 20 / 56); there the second, random forest ",
   "dummy equals the policy-alternative indicator minus the first, so it acted as a random ",
   "policy-versus-status-quo constant. The new Vittel numbers use the identified coding ",
   "(`output/vittel_respec_lcmnl_ext.csv`, M1 k = 18).", "")
v3 <- do.call(rbind, lapply(DATASETS, function(ds) {
  e <- emp[[ds]]
  do.call(rbind, lapply(names(MODELS), function(m) {
    x <- e$fits[[m]]; ok <- !is.null(x) && !isTRUE(x$unreadable)
    pm <- e$pubm[[m]]
    pb <- if (!is.null(pm)) pm$BIC else if (m == "M1") e$pub$ind else if (m == "M2") e$pub$corr else NA
    data.frame(dataset = ds, model = m, status = e$mstat[[m]],
               BIC = if (ok) fm(x$BIC) else "", LL = if (ok) fm(x$LL, 2) else "", k = if (ok) x$k else "",
               converged = if (ok) x$converged else "", `Apollo code` = if (ok) vstr(x$code) else "",
               nested_ok = if (ok) vstr(x$nested_ok) else "", minutes = if (ok) fm(x$minutes) else "",
               `pub BIC` = fm(pb, 1), `pub LL` = if (!is.null(pm)) fm(pm$LL) else "",
               `pub k` = if (!is.null(pm)) pm$k else "", check.names = FALSE)
  }))
}))
md("### MMNL fits (output/mmnl_bench_std/)", "",
   "M1 independent, fixed constants; M2 correlated, fixed constants; M3 independent, random constants; ",
   "M4 correlated, random constants. Apollo code 4 = relative function convergence. ",
   "Published BIC, LL, k from the appendix tables (one decimal); Mode and Swiss route M2 only from tab:emp_summary (rounded).", "", md_table(v3))
v4 <- do.call(rbind, lapply(DATASETS, function(ds) {
  e <- emp[[ds]]
  data.frame(dataset = ds, file = if (is.null(e$lc)) "missing" else e$lc$file,
             C = if (is.null(e$lc)) "" else e$lc$C, BIC = if (is.null(e$lc)) "" else fm(e$lc$BIC),
             `C range` = if (is.null(e$lc)) "" else sprintf("1..%d%s", e$lc$Cmax, if (e$lc$boundary) " (min at boundary)" else ""),
             `runner-up` = if (is.null(e$lc) || is.na(e$lc$runner_C)) "" else sprintf("C=%d, +%s", e$lc$runner_C, fm(e$lc$runner_dBIC)),
             singular = if (is.null(e$lc)) "" else vstr(e$lc$singular),
             `published file` = if (is.null(e$lcp)) "" else sprintf("%s (C=%d, %s)", e$lcp$file, e$lcp$C, fm(e$lcp$BIC)),
             check.names = FALSE)
}))
md("### BIC-best LCMNL sources", "", md_table(v4))

md("## 2. tab:mmnl_blocked", "")
v5 <- do.call(rbind, lapply(h4res, function(x) {
  data.frame(arm = x$h$label, status = x$status, file = x$h$file,
             n = if (is.null(x$sn)) "" else x$sn$n,
             `MNL / LCMNL / MMNL` = if (is.null(x$sn)) "" else sprintf("%d / %d / %d", x$sn$MNL, x$sn$LCMNL, x$sn$MMNL),
             `over-selects (K* = 1)` = if (is.null(x$sn) || x$h$expected != 80L) "" else x$sn$over,
             `MMNL failed` = if (is.null(x$sn)) "" else x$sn$mmnl_fail,
             `published same cells` = if (is.null(x$spp)) "" else
               sprintf("%d / %d / %d%s", x$spp$MNL, x$spp$LCMNL, x$spp$MMNL,
                       if (x$h$expected == 80L) sprintf(", over %d", x$spp$over) else ""),
             `winner changed` = vstr(x$changed),
             `published (tex)` = if (is.null(x$pt)) "" else sprintf("n=%d: %d / %d / %d", x$pt$n, x$pt$MNL, x$pt$LCMNL, x$pt$MMNL),
             check.names = FALSE)
}))
md(md_table(v5),
   sprintf("Published over-selection, K* = 1: %s of 80 (%s). Overall: new %s; published %s.",
           vstr(over_pub), if (is.na(over_line)) "not found" else tex_ref(over_line),
           if (is.null(ov_new)) "pending" else sprintf("n=%d: %d / %d / %d", ov_new$n, ov_new$MNL, ov_new$LCMNL, ov_new$MMNL),
           if (is.null(ov_pub)) "" else sprintf("n=%d: %d / %d / %d", ov_pub$n, ov_pub$MNL, ov_pub$LCMNL, ov_pub$MMNL)), "")

for (lad in names(stress)) {
  s <- stress[[lad]]; L <- s$L
  md(sprintf("## 3%s. %s (%s ladder, %d seeds per rung)", if (lad == "blocked") "a" else "b",
             L$table, lad, L$seeds), "",
     sprintf("`%s`: %d of %d cells topped up with the standard-Apollo arms%s; a new arm raised the best LL in %d of them (largest rise %.3f LL). Superseded file reproduces the published table: %s.",
             L$file, s$n_new, length(L$rungs) * L$seeds, if (s$running) " (driver running)" else "",
             s$improved, s$max_imp, if (length(s$check_bad)) paste("NO:", paste(s$check_bad, collapse = "; ")) else "yes"), "")
  cellstr <- function(a) if (is.null(a) || a$n_arm == 0) "pending" else
    sprintf("%s (%d/%d)%s", fm(a$mean_gap, 2), a$stuck, a$n_fin, if (a$n_err) sprintf(" [%d err]", a$n_err) else "")
  tabn <- do.call(rbind, lapply(L$rungs, function(rg) {
    v <- vapply(L$arms, function(a) cellstr(s$res[[paste(rg, a)]]$new), character(1))
    as.data.frame(as.list(c(rung = rg, v)), check.names = FALSE)
  }))
  names(tabn) <- c("rung", ARM_LABEL[L$arms])
  md("New: mean gap (stuck / seeds with a finite LL)", "", md_table(tabn))
  tabp <- do.call(rbind, lapply(L$rungs, function(rg) {
    v <- vapply(L$arms, function(a) {
      pt <- s$res[[paste(rg, a)]]$pubtex
      if (is.null(pt)) "" else sprintf("%s (%d/%d)", fm(pt$gap, 2), as.integer(pt$stuck), as.integer(pt$n))
    }, character(1))
    as.data.frame(as.list(c(rung = rg, v)), check.names = FALSE)
  }))
  names(tabp) <- c("rung", ARM_LABEL[L$arms])
  md("Published (version12.tex; the published Apollo column is apollo_ss_published)", "", md_table(tabp))
  tabs <- do.call(rbind, lapply(L$rungs, function(rg) {
    v <- vapply(L$arms, function(a) {
      n <- s$res[[paste(rg, a)]]$new; p <- s$res[[paste(rg, a)]]$pubfile
      sprintf("%s / %s", if (n$n_arm) fm(n$secs) else "pending", if (!is.null(p) && p$n_arm) fm(p$secs) else "-")
    }, character(1))
    as.data.frame(as.list(c(rung = rg, v)), check.names = FALSE)
  }))
  names(tabs) <- c("rung", ARM_LABEL[L$arms])
  md(paste0("Mean seconds per fit, new / published file", if (lad == "blocked")
    " (the published Apollo seconds in the blocked file are stale: dev/stress_replicate_blocked.R:24-27)" else ""), "",
     md_table(tabs))
}

md("## 4. tab:mmnl and tab:mmnl_corr (randomised design)", "")
for (st in names(mmres)) {
  x <- mmres[[st]]; M <- x$M
  md(sprintf("### %s: `%s` is the %s (mtime %s); status: %s", M$table, M$file, x$fs, mtime(M$file), x$status), "")
  tb <- do.call(rbind, lapply(c("K1", "Kge2", "Overall"), function(g) {
    newv <- if (is.null(x$cn)) "pending" else as.character(x$cn[[g]][c("n", M$cols)])
    pubv <- if (is.null(x$pt[[g]])) "" else paste(c(sprintf("n=%d", as.integer(x$pt[[g]]$n)), x$pt[[g]]$v), collapse = " / ")
    supv <- if (is.null(x$cp)) "" else paste(x$cp[[g]][c("n", M$cols)], collapse = " / ")
    data.frame(DGP = g, new = paste(newv, collapse = " / "), `superseded file` = supv,
               `published (tex)` = pubv, check.names = FALSE)
  }))
  md(sprintf("Columns: n / %s. Published columns: %s.", paste(M$cols, collapse = " / "),
             if (st == "mmnl") "n / MNL / LCMNL / MMNL" else "n / LCMNL or MNL / MMNL indep / MMNL corr"), "",
     md_table(tb))
}

md("## 5. H4 misspecification arm, blocked (output/h4_misspec_blocked.csv)", "",
   sprintf("File: %s (mtime %s). Not in version12.", h4m_state, mtime(H4M_FILE)), "")
md(md_table(data.frame(cell = h4m_tab$cell, sigma = h4m_tab$sigma, status = h4m_tab$status,
                       n = h4m_tab$n, MNL = vapply(h4m_tab$MNL, vstr, ""),
                       LCMNL = vapply(h4m_tab$LCMNL, vstr, ""), MMNL = vapply(h4m_tab$MMNL, vstr, ""),
                       `correct reading` = h4m_tab$correct, `n correct` = vapply(h4m_tab$correct_n, vstr, ""),
                       `over-selects` = vapply(h4m_tab$over, vstr, ""), check.names = FALSE)))

inputs <- c(TEX, LANE, file.path(BENCH, sprintf("%s_%s.rds", rep(DATASETS, each = 4), names(MODELS))),
            unique(unlist(LC_NEW)), LC_PUB$Vittel, vapply(H4, `[[`, "", "file"), H4M_FILE,
            vapply(LADDERS, `[[`, "", "file"), vapply(MM, `[[`, "", "file"),
            file.path(SUP, c("k1_mmnl_blocked_full.csv", "discrete_mmnl_blocked.csv",
                             "stress_replicate_blocked.rds", "stress_replicate.rds",
                             "mmnl_results.csv", "mmnl_correlated_results.csv")))
md("## Inputs", "", md_table(data.frame(file = inputs, mtime = vapply(inputs, mtime, ""),
                                         state = vapply(inputs, function(f)
                                           if (startsWith(f, SUP) || f == TEX) (if (file.exists(f)) "reference" else "missing")
                                           else if (f %in% c(unlist(LC_NEW), LC_PUB$Vittel) && !is_new(f))
                                             (if (file.exists(f)) "unchanged (LCMNL, not rerun)" else "missing")
                                           else file_state(f), ""))))
md_path <- file.path(OUT_DIR, "standard_rerun_summary.md")
writeLines(MD, md_path)
cat(sprintf("wrote %s (%d rows) and %s\n", csv_path, nrow(csv), md_path))
