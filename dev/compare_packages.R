# dev/compare_packages.R
# ============================================================================
# Head-to-head comparison harness for klue's clustering-assisted LCMNL
# initialisation against the ACTUAL starting-value / estimation routines of the
# two dominant R choice packages: Apollo and gmnl.
#
# On a simulated LCMNL dataset (fixed C = true number of classes) we run, per
# method: final log-likelihood, gap to the best LL found, converged (T/F), and
# wall-clock seconds.
#
# METHODS:
#   1. klue clustering, direct-ML  : klue_lcmnl(db, C, dgp, estimator="ml")  [reference]
#   2. klue clustering, EM         : klue_lcmnl(db, C, dgp, estimator="em")
#   3. Apollo apollo_lcEM          : Apollo's latent-class EM, as in its EM example
#   4. Apollo searchStart + estimate: Apollo's scatter-search start finder at its
#                                     default settings + gradient ML (BGW)
# The Apollo methods follow Apollo's official example scripts; see the block
# "Apollo latent-class model code" below.
#   5. gmnl model="lc", Q=C        : gmnl with its DEFAULT starting values (gradient ML)
#
# All runs are SINGLE CORE and data is kept tiny. Do NOT raise nCores.
#
# Run from the repository root, one core:
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     nice -n 15 Rscript dev/compare_packages.R
# Run directly, it prints the validation gate (main() at the end of this file:
# klue_simulate seed 42, C = 3) and writes no result file. Sourced (the stress
# ladders use sys.source), it only defines the runners. Either way it stops
# unless apollo, idefix, mclust and cluster match dev/pinned_versions.R; the
# gmnl arm also checks gmnl and mlogit (KLUE_ALLOW_VERSION_DRIFT=1 turns the
# stop into a warning).
# Seeds: klue's clustering starts use seed 123; run_klue_em_rp draws partition
# p under set.seed(seed*100000 + 50000 + p); apollo_searchStart draws its
# candidates after set.seed(17) (Apollo's default seed 13, plus 4). klue's
# BFGS and EM, apollo_lcEM, apollo_estimate (BGW) and gmnl use no random
# numbers that affect the estimates.
# ============================================================================

# gmnl 1.1.3.2 needs the pre-dfidx mlogit.data structure (mlogit < 1.1); with
# mlogit >= 1.1 it dies in its own code ("argument is not a matrix") because
# mlogit.data now returns a dfidx object. If a compatible mlogit has been
# installed to an isolated library, put it first on the path so gmnl picks it
# up. Nothing else here uses mlogit (klue has its own engine; Apollo doesn't
# use mlogit), so shadowing only affects gmnl. Build it with:
#   remotes::install_version("mlogit","1.0-3.1", lib="~/R/oldmlogit_lib", dependencies=FALSE)
# KLUE_OLD_MLOGIT_LIB points elsewhere if needed. The gmnl arm stops when the
# side library is missing (run_gmnl), so a ladder cannot record a silent NA.
OLD_MLOGIT_LIB <- path.expand(Sys.getenv("KLUE_OLD_MLOGIT_LIB", "~/R/oldmlogit_lib"))
if (dir.exists(OLD_MLOGIT_LIB)) .libPaths(c(OLD_MLOGIT_LIB, .libPaths()))

suppressWarnings(suppressMessages({
  pkgload::load_all("klue", quiet = TRUE)   # NB: do NOT library(klue); installed copy is stale
  library(apollo)
}))
# Stop unless the library matches the versions behind the paper's tables.
source("dev/pinned_versions.R", local = TRUE)
check_pinned_versions(c("apollo", "idefix", "mclust", "cluster"))

options(mc.cores = 1L)        # thermally constrained laptop; a 2-core batch is running
HAVE_GMNL <- requireNamespace("gmnl",   quietly = TRUE)
HAVE_MLOGIT <- requireNamespace("mlogit", quietly = TRUE)

# ----------------------------------------------------------------------------
# Generic helpers
# ----------------------------------------------------------------------------

get_apollo_ll <- function(m) {
  if (is.null(m)) return(NA_real_)
  if (!is.null(m$maximum))  return(as.numeric(m$maximum))
  if (!is.null(m$LLout))    return(as.numeric(m$LLout[1]))
  am <- attr(m, "maximum")
  if (!is.null(am)) return(as.numeric(am))
  NA_real_
}

# Did an Apollo model converge? successfulEstimation flag if present, else
# the presence of a finite maximum.
apl_converged <- function(m) {
  if (is.null(m)) return(FALSE)
  if (!is.null(m$successfulEstimation)) return(isTRUE(m$successfulEstimation))
  isTRUE(is.finite(get_apollo_ll(m)))
}

# ----------------------------------------------------------------------------
# Apollo latent-class model code, written as Apollo's official examples write
# it (apollochoicemodelling.com, "6 Mixture models/LC_no_covariates.r" for
# apollo_estimate and apollo_searchStart, "8 Alternative estimation
# approaches/EM_LC_no_covariates.r" for apollo_lcEM), generated for any number
# of alternatives J, attributes and classes C:
#   - class-specific parameters asc_alt{j}_{c}, b_x{a}_{c}, b_price_{c}; the
#     reference alternative J keeps its constant in apollo_beta at 0 and fixed
#     (EM example); the last class is the allocation reference (delta_C = 0
#     fixed), as in both examples;
#   - allocation: apollo_classAlloc for apollo_estimate (LC example);
#     an explicit logit list for apollo_lcEM (EM example);
#   - the class loop for(s in 1:C) with the class count written as a literal
#     (LC example: 1:2), or for(s in 1:length(pi_values)) for apollo_lcEM,
#     which requires that form (EM example).
# The two forms are not interchangeable in Apollo 0.3.5: apollo_estimate with
# analytic gradients needs apollo_classAlloc, and the length(pi_values) loop
# combined with apollo_classAlloc gives zero analytic gradients for every
# class-specific parameter (Apollo cannot expand that loop), which validation
# reports as "parameters do not influence the log-likelihood".
# ----------------------------------------------------------------------------
.apollo_lc_check_dgp <- function(dgp) {
  if (!identical(as.integer(dgp$asc_map), c(seq_len(dgp$n_alternatives - 1L), 0L)))
    stop("the Apollo LC harness supports the default asc_map only", call. = FALSE)
}
.apollo_lc_pnames <- function(dgp)
  c(paste0("asc_alt", seq_len(dgp$n_alternatives)),
    paste0("b_x", seq_len(dgp$n_generic)), "b_price")
.apollo_lc_fixed <- function(C, dgp, ref = c("last", "first"))
  c(paste0("asc_alt", dgp$n_alternatives, "_", 1:C),
    paste0("delta_", if (match.arg(ref) == "last") C else 1L))

# Starting values for the Apollo LC arms.
#   "mnl_scaled" (default): Apollo's examples hand-set sign-informed values that
#     differ across classes (class b = half of class a in LC_no_covariates.r),
#     and Apollo's manual reads starts from a simpler model (apollo_readBeta).
#     Here class c starts at the pooled-MNL estimates divided by c (so class 2
#     is half of class 1, as in the example), constants included; deltas 0.
#     klue's MNL matches Apollo's to about 1e-5. Without a start search, a
#     gradient fit from this start can collapse a class (seen on test data:
#     'Singular convergence', a delta near -255); the runners report Apollo's
#     successfulEstimation, and every arm on this start searches or uses EM.
#   "generic": the start of the published Apollo arms (constants 0,
#     b_x = 0.1 * (c - 1), b_price = -0.5, deltas 0), kept to reproduce them.
# order: "parameter" lists apollo_beta parameter by parameter over the classes
#   (asc_alt1_1, asc_alt1_2, ..., delta_C), as the examples do; "class" lists
#   it class by class, the layout of the published arms. apollo_searchStart
#   assigns its candidate columns by position, so the order changes its
#   candidates (not the model).
# ref: the class whose delta is fixed at 0: "last" as in the examples, "first"
#   as in the published arms.
build_apollo_lc_beta <- function(C, dgp, start = c("mnl_scaled", "generic"),
                                 db = NULL, order = c("parameter", "class"),
                                 ref = c("last", "first")) {
  start <- match.arg(start); order <- match.arg(order); ref <- match.arg(ref)
  .apollo_lc_check_dgp(dgp)
  J <- dgp$n_alternatives; ng <- dgp$n_generic; nb <- dgp$n_beta
  mnl <- NULL
  if (start == "mnl_scaled") {
    if (is.null(db)) stop("start = \"mnl_scaled\" needs the database", call. = FALSE)
    fit <- estimate_lcmnl(db, C = 1, start_betas = matrix(0, 1, nb), dgp = dgp,
                          vcov = FALSE)
    if (is.null(fit$par) || !all(is.finite(fit$par)))
      stop("pooled MNL for the Apollo start failed", call. = FALSE)
    mnl <- fit$par
  }
  val <- function(pname, cc) {
    m <- 1 / cc
    if (grepl("^asc_alt", pname)) {
      j <- as.integer(sub("^asc_alt", "", pname))
      if (j == J || is.null(mnl)) 0 else m * mnl[[paste0("asc", j)]]
    } else if (pname == "b_price") {
      if (is.null(mnl)) -0.5 else m * mnl[[paste0("b", nb)]]
    } else {
      a <- as.integer(sub("^b_x", "", pname))
      if (is.null(mnl)) 0.1 * (cc - 1) else m * mnl[[paste0("b", a)]]
    }
  }
  pn <- .apollo_lc_pnames(dgp)
  grid <- if (order == "parameter") expand.grid(c = 1:C, p = pn, stringsAsFactors = FALSE)
          else expand.grid(p = pn, c = 1:C, stringsAsFactors = FALSE)
  beta <- setNames(mapply(val, grid$p, grid$c), paste0(grid$p, "_", grid$c))
  c(beta, setNames(rep(0, C), paste0("delta_", 1:C)))
}

# apollo_lcPars: one list per class-specific parameter plus pi_values, from
# apollo_classAlloc (LC example) or an explicit logit (EM example).
make_apollo_lcPars <- function(C, dgp = DGP_DEFAULT, alloc = c("classAlloc", "logit")) {
  alloc <- match.arg(alloc)
  lines <- c("function(apollo_beta, apollo_inputs) {", "  lcpars <- list()")
  for (pname in .apollo_lc_pnames(dgp))
    lines <- c(lines, sprintf('  lcpars[["%s"]] <- list(%s)', pname,
                              paste(paste0(pname, "_", 1:C), collapse = ", ")))
  if (alloc == "classAlloc") {
    lines <- c(lines,
      '  V <- list()',
      sprintf('  V[["class_%d"]] <- delta_%d', 1:C, 1:C),
      sprintf('  classAlloc_settings <- list(classes = c(%s), utilities = V)',
              paste(sprintf('class_%d = %d', 1:C, 1:C), collapse = ", ")),
      '  lcpars[["pi_values"]] <- apollo_classAlloc(classAlloc_settings)')
  } else {
    den <- paste(sprintf("exp(delta_%d)", 1:C), collapse = " + ")
    lines <- c(lines, sprintf('  lcpars[["pi_values"]] <- list(%s)',
                              paste(sprintf("class_%d = exp(delta_%d)/(%s)", 1:C, 1:C, den),
                                    collapse = ", ")))
  }
  lines <- c(lines, '  return(lcpars)', '}')
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- globalenv()
  fn
}

# apollo_probabilities for the LC model, with the examples' class loop.
make_apollo_probabilities_lc <- function(C, dgp, loop = c("literal", "length")) {
  loop <- match.arg(loop)
  J <- dgp$n_alternatives; ng <- dgp$n_generic
  util <- vapply(1:J, function(j)
    sprintf('    V[["alt%d"]] <- %s', j, paste(c(
      sprintf('asc_alt%d[[s]]', j),
      sprintf('b_x%d[[s]] * x%d_%d', 1:ng, 1:ng, j),
      sprintf('b_price[[s]] * price_%d', j)), collapse = " + ")), character(1))
  lines <- c(
    'function(apollo_beta, apollo_inputs, functionality = "estimate") {',
    '  apollo_attach(apollo_beta, apollo_inputs)',
    '  on.exit(apollo_detach(apollo_beta, apollo_inputs))',
    '  P <- list()',
    sprintf('  mnl_settings <- list(alternatives = c(%s), avail = list(%s), choiceVar = CHOICE)',
            paste(sprintf('alt%d = %d', 1:J, 1:J), collapse = ", "),
            paste(sprintf('alt%d = 1', 1:J), collapse = ", ")),
    sprintf('  for (s in %s) {', if (loop == "literal") paste0("1:", C) else "1:length(pi_values)"),
    '    V <- list()',
    util,
    '    mnl_settings$utilities <- V',
    '    P[[paste0("Class_", s)]] <- apollo_mnl(mnl_settings, functionality)',
    '    P[[paste0("Class_", s)]] <- apollo_panelProd(P[[paste0("Class_", s)]], apollo_inputs, functionality)',
    '  }',
    '  lc_settings <- list(inClassProb = P, classProb = pi_values)',
    '  P[["model"]] <- apollo_lc(lc_settings, apollo_inputs, functionality)',
    '  P <- apollo_prepareProb(P, apollo_inputs, functionality)',
    '  return(P)',
    '}'
  )
  fn <- eval(parse(text = paste(lines, collapse = "\n")))
  environment(fn) <- globalenv()
  fn
}

# ----------------------------------------------------------------------------
# METHOD 1 & 2: klue
# ----------------------------------------------------------------------------
run_klue <- function(db, C, dgp, estimator) {
  t0 <- Sys.time()
  res <- tryCatch(klue_lcmnl(db, C, dgp = dgp, estimator = estimator),
                  error = function(e) { message("klue ", estimator, " error: ",
                                                 conditionMessage(e)); NULL })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(res)) return(list(LL = NA_real_, converged = FALSE, seconds = secs,
                                note = "errored"))
  list(LL = as.numeric(res$LL), converged = isTRUE(res$converged),
       seconds = secs, note = "")
}

# ----------------------------------------------------------------------------
# METHOD 2b: EM from random-partition starts -- the combination Stata's lclogit
# and Latent GOLD ship by default (random assignment of respondents to classes,
# class-wise MNL fits as starts, EM). Best of `n_starts` partitions, mirroring
# the six-start clustering budget of run_klue; partitions are deterministic
# given `seed` (thread the dataset seed through for independent draws per cell).
# ----------------------------------------------------------------------------
run_klue_em_rp <- function(db, C, dgp, n_starts = 6L, seed = 1L) {
  t0 <- Sys.time()
  res <- tryCatch({
    starts <- get_random_partition_starts(db, C, n_starts = n_starts,
                                          seed = seed, dgp = dgp)
    fits <- lapply(starts, function(st)
      tryCatch(estimate_lcmnl_em(db, C, start_betas = st$betas,
                                 start_shares = st$shares, dgp = dgp),
               error = function(e) NULL))
    lls <- vapply(fits, function(f) if (is.null(f) || !is.finite(f$LL)) -Inf
                                    else f$LL, numeric(1))
    if (all(!is.finite(lls))) NULL else fits[[which.max(lls)]]
  }, error = function(e) { message("klue em_rp error: ", conditionMessage(e)); NULL })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(res)) return(list(LL = NA_real_, converged = FALSE, seconds = secs,
                                note = "errored"))
  list(LL = as.numeric(res$LL), converged = isTRUE(res$converged),
       seconds = secs, note = "")
}

# ----------------------------------------------------------------------------
# Apollo LC plumbing shared by the Apollo methods. Apollo requires the apollo_*
# objects in the global environment: .apollo_lc_setup assigns them (as the
# examples' top-level script does) and validates inputs with explicit arguments
# (officially supported); .apollo_result clears them and returns a runner's
# result list. form = "estimate" builds the LC-example model, form = "em" the
# EM-example model (explicit-logit allocation, length(pi_values) loop, a
# weights column, noValidation/noDiagnostics as in that example).
# estimate_settings: writeIter = FALSE (no iterations file) and silent = TRUE
# (log size) change no estimate; everything else is Apollo's default.
# ----------------------------------------------------------------------------
.apollo_lc_setup <- function(db, C, dgp, model_name, form = c("estimate", "em"),
                             start = "mnl_scaled", order = "parameter",
                             ref = "last") {
  form <- match.arg(form)
  cleanup_apollo()
  apollo_control <<- c(list(modelName = model_name, modelDescr = "LC",
                            indivID = "ID", nCores = 1L, outputDirectory = tempdir()),
                       if (form == "em") list(noValidation = TRUE, noDiagnostics = TRUE))
  if (form == "em") db$weights <- 1
  apollo_beta          <<- tryCatch(build_apollo_lc_beta(C, dgp, start, db = db,
                                                         order = order, ref = ref),
                                    error = function(e) { message(model_name, " start: ",
                                                                  conditionMessage(e)); NULL })
  if (is.null(apollo_beta)) return(FALSE)
  apollo_fixed         <<- .apollo_lc_fixed(C, dgp, ref)
  apollo_lcPars        <<- make_apollo_lcPars(C, dgp, if (form == "em") "logit" else "classAlloc")
  apollo_probabilities <<- make_apollo_probabilities_lc(C, dgp, if (form == "em") "length" else "literal")
  apollo_inputs        <<- tryCatch(
    apollo_validateInputs(apollo_beta = apollo_beta, apollo_fixed = apollo_fixed,
                          database = db, apollo_control = apollo_control),
    error = function(e) { message(model_name, " validateInputs error: ",
                                  conditionMessage(e)); NULL })
  !is.null(apollo_inputs)
}
.apollo_result <- function(m = NULL, t0 = NULL, note = "") {
  cleanup_apollo()
  list(LL = get_apollo_ll(m), converged = apl_converged(m),
       seconds = if (is.null(t0)) NA_real_ else as.numeric(difftime(Sys.time(), t0, units = "secs")),
       note = note)
}
# The start goes into the global apollo_beta first, as the examples do
# (apollo_beta = apollo_searchStart(...); apollo_estimate(apollo_beta, ...)):
# Apollo 0.3.5 checks each rewrite of the model code (loop expansion, scaling)
# by evaluating it at the GLOBAL apollo_beta, and when that differs from the
# start it keeps the unexpanded class loop, whose analytic gradients are zero
# ("parameters do not influence the log-likelihood").
.apollo_est <- function(beta, label = "apollo_estimate") {
  apollo_beta <<- beta
  tryCatch(apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                           estimate_settings = list(writeIter = FALSE, silent = TRUE)),
           error = function(e) { message(label, " error: ", conditionMessage(e)); NULL })
}

# ----------------------------------------------------------------------------
# METHOD 3: Apollo apollo_lcEM, as in EM_LC_no_covariates.r (EMmaxIterations =
# 100, the example's setting; Apollo then runs its post-EM ML step).
# ----------------------------------------------------------------------------
run_apollo_lcEM <- function(db, C, dgp, start = "mnl_scaled") {
  t0 <- Sys.time()
  # apollo_lcEM reads the data from a global `database`, which the example
  # script defines at top level; set it for this call and restore afterwards.
  had_db <- exists("database", envir = globalenv(), inherits = FALSE)
  if (had_db) old_db <- get("database", envir = globalenv())
  on.exit(if (had_db) assign("database", old_db, envir = globalenv())
          else if (exists("database", envir = globalenv(), inherits = FALSE))
            rm("database", envir = globalenv()), add = TRUE)
  db$weights <- 1
  assign("database", db, envir = globalenv())
  if (!.apollo_lc_setup(db, C, dgp, "cmp_lcEM", form = "em", start = start))
    return(.apollo_result(t0 = t0, note = "setup failed"))
  m <- tryCatch(
    apollo_lcEM(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                lcEM_settings     = list(EMmaxIterations = 100),
                estimate_settings = list(writeIter = FALSE, silent = TRUE)),
    error = function(e) { message("apollo_lcEM error: ", conditionMessage(e)); NULL })
  .apollo_result(m, t0, if (is.null(m)) "errored" else "")
}

# ----------------------------------------------------------------------------
# METHOD 4: Apollo apollo_searchStart + apollo_estimate, the examples' call
# pattern (apollo_beta <- apollo_searchStart(...); apollo_estimate(...)).
#   settings = NULL (default): Apollo 0.3.5's defaults, i.e. no
#     searchStart_settings passed: 100 candidates in apollo_beta +- 0.1,
#     smartStart = FALSE, 5 stages, dTest = 1, gTest = 1e-3, llTest = 3,
#     bfgsIter = 20.
#   SS_PUBLISHED_RANDOMISED / SS_PUBLISHED_BLOCKED: the non-default
#     configurations behind the published Apollo arms, run with the published
#     "generic" start, class-major parameter order and delta_1 fixed, which
#     reproduces them (moderate|seed2 of the blocked ladder: -2374.296231).
#     `window` sets apolloBetaMin/Max for every parameter.
# ----------------------------------------------------------------------------
SS_PUBLISHED_RANDOMISED <- list(nCandidates = 30, smartStart = TRUE, dTest = 1,
                                gTest = 10, maxStages = 3, bfgsIter = 20)
SS_PUBLISHED_BLOCKED    <- list(nCandidates = 30, smartStart = FALSE, dTest = 1,
                                gTest = 10, maxStages = 3, bfgsIter = 20,
                                window = c(-3, 3))
run_apollo_searchStart <- function(db, C, dgp, settings = NULL, start = "mnl_scaled",
                                   order = "parameter", ref = "last") {
  t0 <- Sys.time()
  if (!.apollo_lc_setup(db, C, dgp, "cmp_search", start = start, order = order, ref = ref))
    return(.apollo_result(t0 = t0, note = "setup failed"))
  if (!is.null(settings$window)) {
    nm <- names(apollo_beta)
    settings$apolloBetaMin <- setNames(rep(settings$window[1], length(nm)), nm)
    settings$apolloBetaMax <- setNames(rep(settings$window[2], length(nm)), nm)
    settings$window <- NULL
  }
  start_beta <- tryCatch(
    if (is.null(settings))
      apollo_searchStart(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)
    else
      apollo_searchStart(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                         searchStart_settings = settings),
    error = function(e) { message("apollo_searchStart error: ", conditionMessage(e)); NULL })
  m <- if (!is.null(start_beta)) .apollo_est(start_beta)
  .apollo_result(m, t0, if (is.null(start_beta)) "searchStart failed"
                        else if (is.null(m)) "estimate errored" else "")
}
run_apollo_ss_published_randomised <- function(db, C, dgp)
  run_apollo_searchStart(db, C, dgp, settings = SS_PUBLISHED_RANDOMISED, start = "generic",
                         order = "class", ref = "first")
run_apollo_ss_published_blocked <- function(db, C, dgp)
  run_apollo_searchStart(db, C, dgp, settings = SS_PUBLISHED_BLOCKED, start = "generic",
                         order = "class", ref = "first")

# ----------------------------------------------------------------------------
# METHOD 4c: Apollo's estimator from klue's six clustering starts, keeping the
# best. Same starts as METHOD 1, different optimiser (apollo_estimate, BGW):
# separates the contribution of the start from that of the estimator.
# Constants start at 0 (as in klue); shares map to delta_c = log(pi_c / pi_C)
# with delta_C fixed, the examples' normalisation.
# ----------------------------------------------------------------------------
run_apollo_from_clustering <- function(db, C, dgp) {
  t0 <- Sys.time()
  starts <- tryCatch(Filter(Negate(is.null), get_all_starts(db, C, dgp = dgp)),
                     error = function(e) { message("get_all_starts error: ",
                                                   conditionMessage(e)); list() })
  if (!length(starts) || !.apollo_lc_setup(db, C, dgp, "cmp_clust", start = "generic"))
    return(.apollo_result(t0 = t0, note = "no starts or validateInputs failed"))
  ng <- dgp$n_generic; nb <- dgp$n_beta; J <- dgp$n_alternatives
  fits <- lapply(seq_along(starts), function(i) {
    st <- starts[[i]]; b <- apollo_beta; sh <- pmax(st$shares, 0.01)
    for (cc in 1:C) {
      for (j in 1:J) b[paste0("asc_alt", j, "_", cc)] <- 0
      for (a in 1:ng) b[paste0("b_x", a, "_", cc)] <- st$betas[cc, a]
      b[paste0("b_price_", cc)] <- st$betas[cc, nb]
      b[paste0("delta_", cc)]   <- log(sh[cc] / sh[C])
    }
    .apollo_est(b, sprintf("apollo_estimate (start %d)", i))
  })
  lls <- vapply(fits, get_apollo_ll, numeric(1))
  if (all(!is.finite(lls))) return(.apollo_result(t0 = t0, note = "all starts errored"))
  .apollo_result(fits[[which.max(lls)]], t0)
}

# ----------------------------------------------------------------------------
# METHOD 5: gmnl model="lc"
#
# gmnl wants data in dfidx/mlogit long format: one row per (choice situation,
# alternative). We reshape the wide klue database into long form carrying, for
# each (ID, TASK, alt): the chosen indicator, the generic attributes
# (x1..x{ng}, price) for that alt, an alt id, and a choice-situation id (chid).
# The reference alternative carries no ASC, matching klue's spec, which gmnl
# handles automatically (it drops one ASC).
# ----------------------------------------------------------------------------
reshape_for_gmnl <- function(db, dgp) {
  J <- dgp$n_alternatives; ng <- dgp$n_generic
  n_sit <- nrow(db)                      # one row per (ID, TASK) = one choice situation
  chid  <- seq_len(n_sit)

  long <- data.frame(
    chid  = rep(chid, each = J),
    id    = rep(db$ID, each = J),
    alt   = rep(1:J, times = n_sit),
    choice = as.vector(t(outer(db$CHOICE, 1:J, FUN = `==`)))  # TRUE for chosen alt
  )
  # Attribute columns: gather alt-specific x{a}_{j} / price_{j} into generic cols.
  for (a in 1:ng) {
    mat <- as.matrix(db[, paste0("x", a, "_", 1:J), drop = FALSE])  # n_sit x J
    long[[paste0("x", a)]] <- as.vector(t(mat))
  }
  pmat <- as.matrix(db[, paste0("price_", 1:J), drop = FALSE])
  long$price <- as.vector(t(pmat))
  long
}

run_gmnl <- function(db, C, dgp) {
  if (!HAVE_GMNL || !HAVE_MLOGIT) {
    return(list(LL = NA_real_, converged = FALSE, seconds = NA_real_,
                note = paste0("not run: missing package(s) ",
                              paste(c(if (!HAVE_GMNL) "gmnl",
                                      if (!HAVE_MLOGIT) "mlogit"), collapse = ", "))))
  }
  # gmnl 1.1-3.2 errors on mlogit >= 1.1, which used to end as LL = NA with
  # note "errored"; stop instead, before any time is spent.
  check_pinned_versions(c("gmnl", "mlogit"))
  ng <- dgp$n_generic
  long <- reshape_for_gmnl(db, dgp)

  t0 <- Sys.time()
  m <- tryCatch({
    suppressWarnings(suppressMessages({
      mdat <- mlogit::mlogit.data(long, choice = "choice", shape = "long",
                                  alt.var = "alt", chid.var = "chid",
                                  id.var = "id")
      # gmnl LC 5-part Formula: part 1 = attributes (class-specific coefs in LC);
      # part 2 = "1" => class-specific ASCs (matches klue's per-class ASCs, ref
      # alt dropped automatically); parts 3-4 = 0; part 5 = "1" = class-membership
      # intercepts (constant class shares). Validated: this reaches klue's exact
      # global (LL identical, k=23) under mlogit 1.0-3.1.
      rhs_generic <- paste(c(paste0("x", 1:ng), "price"), collapse = " + ")
      fml <- stats::as.formula(paste0("choice ~ ", rhs_generic, " | 1 | 0 | 0 | 1"))
      gmnl::gmnl(fml, data = mdat, model = "lc", Q = C, panel = TRUE,
                 method = "bfgs", print.level = 0)
    }))
  }, error = function(e) { message("gmnl error: ", conditionMessage(e)); NULL })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(m)) return(list(LL = NA_real_, converged = FALSE, seconds = secs,
                              note = "errored"))
  ll   <- tryCatch(as.numeric(stats::logLik(m)), error = function(e) NA_real_)
  # gmnl 1.1-3.2 keeps maxLik's return code in m$logLik$code (0 = converged)
  conv <- tryCatch(isTRUE(m$logLik$code == 0) || isTRUE(m$logLik$convergence == 0) ||
                   isTRUE(m$convergence == 0), error = function(e) NA)
  if (is.na(conv)) conv <- is.finite(ll)
  list(LL = ll, converged = isTRUE(conv), seconds = secs, note = "")
}

# ----------------------------------------------------------------------------
# Driver: run all methods on one dataset, assemble the comparison table.
# ----------------------------------------------------------------------------
compare_on_dataset <- function(db, C, dgp, label = "") {
  cat(sprintf("\n=== Comparison on dataset: %s (C=%d, N=%d, rows=%d) ===\n",
              label, C, length(unique(db$ID)), nrow(db)))

  methods <- list(
    "klue_ml (ref)"        = function() run_klue(db, C, dgp, "ml"),
    "klue_em"              = function() run_klue(db, C, dgp, "em"),
    "apollo_lcEM"          = function() run_apollo_lcEM(db, C, dgp),
    "apollo_searchStart"   = function() run_apollo_searchStart(db, C, dgp),
    "gmnl_lc"              = function() run_gmnl(db, C, dgp)
  )

  rows <- list()
  for (nm in names(methods)) {
    cat(sprintf("  -> running %-20s ...\n", nm))
    r <- methods[[nm]]()
    rows[[nm]] <- data.frame(method = nm, LL = r$LL, converged = r$converged,
                             seconds = r$seconds, note = r$note,
                             stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  best_ll <- suppressWarnings(max(tab$LL[is.finite(tab$LL)], na.rm = TRUE))
  tab$gap_to_best <- ifelse(is.finite(tab$LL), tab$LL - best_ll, NA_real_)

  tab <- tab[, c("method", "LL", "gap_to_best", "converged", "seconds", "note")]
  tab
}

print_table <- function(tab) {
  out <- tab
  out$LL          <- sprintf("%.3f", out$LL)
  out$gap_to_best <- ifelse(is.na(tab$gap_to_best), "NA",
                            sprintf("%+.3f", tab$gap_to_best))
  out$seconds     <- ifelse(is.na(tab$seconds), "NA", sprintf("%.1f", tab$seconds))
  print(out, row.names = FALSE)
}

# ============================================================================
# HARD VALIDATION GATE: well-separated dataset -> global optimum unambiguous.
# Every correctly-implemented method must reach essentially the SAME LL.
# Apollo (both routines) and gmnl should land within ~2 LL of klue's best.
# ============================================================================
main <- function() {
  dgp <- klue_dgp(n_generic = 4, n_alternatives = 3)
  C   <- 3

  cat("\n############################################################\n")
  cat("# VALIDATION GATE: separation=1.5, true_K=3, N_per_class=100, T=12, seed=42\n")
  cat("############################################################\n")
  d_easy <- klue_simulate(N_per_class = 100, T_tasks = 12, true_K = C,
                          separation = 1.5, heterogeneity = 0.2, seed = 42, dgp = dgp)
  tab <- compare_on_dataset(d_easy$database, C, dgp, label = "EASY (separation=1.5)")

  cat("\n----- VALIDATION TABLE (easy dataset) -----\n")
  print_table(tab)

  # Gate verdict
  ref_ll  <- tab$LL[tab$method == "klue_ml (ref)"]
  best_ll <- suppressWarnings(max(tab$LL[is.finite(tab$LL)], na.rm = TRUE))
  tol <- 2.0
  check <- function(nm) {
    ll <- tab$LL[tab$method == nm]
    if (length(ll) == 0 || !is.finite(ll)) return(sprintf("%-20s : NOT RUN / NA", nm))
    gap <- ll - best_ll
    sprintf("%-20s : LL=%.3f  gap=%+.3f  %s", nm, ll, gap,
            if (abs(gap) <= tol) "PASS" else "FAIL (>2 LL below best)")
  }
  cat("\n----- GATE VERDICT (tolerance = 2 LL vs best) -----\n")
  cat(sprintf("best LL found = %.3f ; klue_ml ref = %.3f\n", best_ll, ref_ll))
  for (nm in c("klue_em", "apollo_lcEM", "apollo_searchStart", "gmnl_lc"))
    cat("  ", check(nm), "\n")

  invisible(tab)
}

# Autorun only when invoked directly (Rscript), not when sys.source'd for its
# functions (e.g. by dev/stress_replicate.R). sys.nframe() is 0 only at the
# Rscript top level.
# (The single-seed stress ladder that used to live here is superseded by
# dev/stress_replicate.R and dev/stress_replicate_blocked.R.)
if (sys.nframe() == 0L) main()
