# =============================================================================
# Gate 5: Extend the empirical class-count search until BIC turns (no fixed cap)
# for the three datasets whose BIC is monotone at the current C=6 boundary:
# Vittel, Electricity, Swissmetro.
#
# For each dataset, fit LCMNL via the clustering multistart for C = 1, 2, ...,
# stopping when BIC has risen for PATIENCE consecutive steps past its minimum
# (an interior optimum), OR estimation fails to converge (degeneracy), OR a
# safety cap is reached -- in which case BIC is "still falling at C = cap",
# reported rather than hidden. Resumable (writes after every fit and continues
# where an interrupted run stopped).
# Seeds: klue's clustering starts use seed 123. The set.seed(20260608 + C)
# before each fit has no effect on the fits: klue seeds its clustering itself
# and restores the caller's stream, and nothing else draws random numbers.
#
# HEAVY -- parked until a compute window. Swissmetro (N=1,004) and Vittel
# (N=930, 8 attributes) at high C are the long poles. The cap is a backstop,
# not a scientific limit; raise it if BIC has not turned. Run from the
# repository root, prefixing each run with
# OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 nice -n 15:
#   Run (all three, cap 12):       Rscript dev/run_emp_cmax_extend.R
#   Run (one dataset, custom cap):  Rscript dev/run_emp_cmax_extend.R Swissmetro 10
# Writes: output/emp_cmax_extend_<dataset>.csv  (one row per C; fresh filenames,
#         does NOT overwrite the existing *_lcmnl_results.csv).
# =============================================================================

SEED     <- 20260608L
PATIENCE <- 2L    # consecutive BIC rises past the minimum that confirm a turn
CAP_DEF  <- 12L   # safety backstop; user chose "no cap" -- raise if still falling
C_START  <- 1L
COOL     <- as.numeric(Sys.getenv("KLUE_SLEEP", "0"))  # cooling pause (s) between C fits

if (file.exists("klue/DESCRIPTION")) {
  suppressWarnings(suppressMessages(pkgload::load_all("klue", quiet = TRUE)))
} else {
  suppressWarnings(suppressMessages(library(klue)))
}
if (!dir.exists("output")) dir.create("output")

# Keep Apollo/MMNL parallelism modest (the multistart itself is the cost here).
options(klue.mmnl.n_cores = max(1L, as.integer(Sys.getenv("KLUE_CORES", "2"))))

SPECS <- list(
  Vittel      = list(script = "R/empirical_application.R", skip = "empirical.skip",
                     loader = "load_empirical_database",   dgp = "EMP_DGP"),
  Electricity = list(script = "R/empirical_electricity.R", skip = "elec.skip",
                     loader = "load_electricity_database", dgp = "ELEC_DGP"),
  Swissmetro  = list(script = "R/empirical_swissmetro.R",  skip = "sm.skip",
                     loader = "load_swissmetro_database",  dgp = "SM_DGP")
)

# Replay saved converged rows to restore (best_bic, best_c, rises) on resume.
replay_state <- function(df) {
  best_bic <- Inf; best_c <- NA_integer_; rises <- 0L
  for (i in order(df$C)) {
    if (!isTRUE(df$converged[i])) next
    if (df$BIC[i] < best_bic - 1e-6) {
      best_bic <- df$BIC[i]; best_c <- df$C[i]; rises <- 0L
    } else rises <- rises + 1L
  }
  list(best_bic = best_bic, best_c = best_c, rises = rises)
}

fit_path <- function(name, cap = CAP_DEF) {
  s <- SPECS[[name]]
  options(structure(list(TRUE), names = s$skip))   # disable that script's autorun
  source(s$script, local = FALSE)                  # defines loader + DGP globally
  db  <- get(s$loader)()
  dgp <- get(s$dgp)
  out <- file.path("output", sprintf("emp_cmax_extend_%s.csv", name))

  rows <- list(); start <- C_START
  best_bic <- Inf; best_c <- NA_integer_; rises <- 0L
  if (file.exists(out)) {                          # resume
    prev <- read.csv(out, stringsAsFactors = FALSE)
    if (any(!prev$converged)) {
      cat(sprintf("[%s] already terminated (non-convergence on record); skipping\n", name))
      return(invisible(prev))
    }
    rows <- lapply(seq_len(nrow(prev)), function(i) prev[i, , drop = FALSE])
    st <- replay_state(prev); best_bic <- st$best_bic
    best_c <- st$best_c; rises <- st$rises
    if (rises >= PATIENCE) {
      cat(sprintf("[%s] BIC already turned at C=%d; nothing to do\n", name, best_c))
      return(invisible(prev))
    }
    start <- max(prev$C) + 1L
    cat(sprintf("[%s] resuming from C=%d (best=%.1f@C%d)\n", name, start, best_bic, best_c))
  }

  stop_reason <- sprintf("hit safety cap C=%d (BIC still falling)", cap)
  for (C in start:cap) {
    set.seed(SEED + C)                             # no effect on the fits (see header)
    t0 <- Sys.time()
    m  <- tryCatch(klue_lcmnl(db, C, dgp = dgp),
                   error = function(e) NULL)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (is.null(m) || !isTRUE(m$converged)) {
      rows[[length(rows) + 1L]] <- data.frame(
        dataset = name, C = C, converged = FALSE, LL = NA_real_, k = NA_integer_,
        BIC = NA_real_, AIC = NA_real_, ICL = NA_real_, ICL_BIC = NA_real_,
        min_class_prob = NA_real_, best_method = NA_character_, secs = round(secs, 1))
      write.csv(do.call(rbind, rows), out, row.names = FALSE)
      stop_reason <- sprintf("non-convergence/degenerate at C=%d", C); break
    }

    mcp <- if (!is.null(m$class_probs)) min(m$class_probs) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = name, C = C, converged = TRUE, LL = m$LL, k = m$k,
      BIC = m$BIC, AIC = m$AIC, ICL = m$ICL,
      ICL_BIC = ifelse(is.null(m$ICL_BIC), NA_real_, m$ICL_BIC),
      min_class_prob = mcp,
      best_method = ifelse(is.null(m$best_method) || is.na(m$best_method),
                           "MNL", m$best_method),
      secs = round(secs, 1))
    write.csv(do.call(rbind, rows), out, row.names = FALSE)   # save after every fit
    if (COOL > 0) Sys.sleep(COOL)                              # cooling pause

    is_min <- m$BIC < best_bic - 1e-6
    cat(sprintf("[%s] C=%2d BIC=%.1f %s minP=%.3f %.0fs\n", name, C, m$BIC,
                if (is_min) "(new min)" else sprintf("(>min@C%d)", best_c),
                ifelse(is.na(mcp), -1, mcp), secs))
    if (is_min) { best_bic <- m$BIC; best_c <- C; rises <- 0L }
    else {
      rises <- rises + 1L
      if (rises >= PATIENCE) {
        stop_reason <- sprintf("BIC turned: interior minimum at C=%d", best_c); break
      }
    }
  }
  cat(sprintf(">> %s: %s ; best C=%s, BIC=%.1f  [%s]\n",
              name, stop_reason, as.character(best_c), best_bic, out))
  invisible(do.call(rbind, rows))
}

if (sys.nframe() == 0L) {
  args     <- commandArgs(trailingOnly = TRUE)
  cap      <- if (length(args) >= 2) as.integer(args[2]) else CAP_DEF
  datasets <- if (length(args) >= 1) args[1] else names(SPECS)
  for (nm in datasets) fit_path(nm, cap = cap)
}
