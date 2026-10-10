# klue (development version)

- The citation of the methodology paper (`citation("klue")`, `CITATION.cff` and
  the message printed by `klue()`) uses its new title, *Classes or Continuum? A
  Pre-Registrable Specification Rule for Latent Class Logit Models*.

# klue 0.10.0 (2026-09-29)

`klue_mmnl()` now estimates the mixed logit the way Apollo's own example
scripts do (`MMNL_preference_space.r`, `MMNL_preference_space_correlated.r`),
and reproduces hand-written scripts in that form to machine precision (same
log-likelihood, estimates and robust standard errors to about 1e-10, same
convergence code and iteration count). MMNL results change; refit anything
estimated with 0.9.x.

- Random coefficients are `mu_p + sigma_p * draws_p` with `sigma_p`
  unconstrained, as in Apollo's examples; `sigma_p` is the standard deviation
  itself (its sign is arbitrary) and `$sigma` reports `|sigma_p|`. Up to 0.9.5
  the independent model used `mu_p + exp(sigma_p) * draws_p`. Apollo's manual
  advises against such transforms (they slow estimation and need not reach the
  unconstrained solution); here they could stop at a saddle point when a
  standard deviation is zero.
- One `apollo_estimate()` call at `n_draws` (default 3000 MLHS) with Apollo's
  default settings (BGW, `silent = FALSE`; `writeIter = FALSE` only skips the
  iterations file). The coarse first stage at `n_draws_stage1` draws and every
  automatic retry and fallback are gone; a failed estimation is reported, not
  retried. The draw count is passed as a double, as in the examples: an
  integer made Apollo 0.3.5's core-count hint overflow (`nrow * draws^2` in
  integer arithmetic), which is why 0.9.x ran Apollo with `silent = TRUE`.
- `converged` is Apollo's `successfulEstimation`, `LL` is Apollo's final
  log-likelihood, and `$apollo_status` records Apollo's `code`, `message` and
  `nIter`; `$BIC_apollo` is Apollo's own BIC (on observations; `$BIC` stays on
  respondents, as for the LCMNL). Up to 0.9.5 a fit counted as converged
  whenever its estimates were finite.
- The correlated model starts from the independent estimates (Apollo's
  `apollo_readBeta` step: same means, the independent standard deviations on
  the Cholesky diagonal). The off-diagonals start at 0 (the current example
  uses 0.01, an earlier revision 0), so that with the same draws the start
  reproduces the independent log-likelihood exactly. Without a `warm_start`,
  the informed start fits the independent model with the same draws first, so
  `nested_ok` is checked on this path too; `klue()` passes its independent fit
  under the informed start and skips the correlated fit when the independent
  one failed.
- The informed independent start takes every pooled-MNL estimate, constants
  included. A log-normal price starts at `mu = log(-b_MNL)`, or at Apollo's
  recommended `-3` when the pooled price coefficient is not negative; a pooled
  MNL that stopped at its iteration cap is still used. Standard deviations
  start at 0.01, as in the current examples (0.1 in 0.9.5).
- `warm_start` must come from klue >= 0.10.0 (its `sigma_*` were log standard
  deviations before) and match the call's `random` set, `price` distribution
  and parameters; `start` is ignored, with a warning, when a warm start is
  given. `nested_ok` also requires the same data (`N`, rows).
- Removed the price box (`mu_price_bounds`, `sigma_price_bounds`, options
  `klue.mmnl.*_bounds`): Apollo 0.3.5 has no `bounds` estimation setting, so
  the box was never applied. `n_draws_stage1`, `mu_price_bounds` and
  `sigma_price_bounds` are accepted with a deprecation warning and ignored;
  they moved to the end of the argument list, so a positional call that
  reached them binds differently.
- Apollo runs on 1 core unless `n_cores` or `options(klue.mmnl.n_cores)` says
  otherwise (it was the number of physical cores minus one). The generated
  model functions live in the global environment, as in the examples, and
  during estimation `apollo_beta`, `apollo_fixed` and `apollo_probabilities`
  are global as in the examples' script; `klue_mmnl()` restores any `apollo_*`
  objects the caller had and the caller's random-number stream (Apollo reseeds
  it for its draws). Model names are unique per call.
- Failure results keep Apollo's status, the estimates where Apollo returned
  them, the internal independent fit and its log; fits carry
  `sd_scale = "linear"`, `start_values` and `start_used`.
- New `memory_saver` argument (option `klue.mmnl.memory_saver`, environment
  variable `KLUE_MMNL_MEMORY_SAVER`): Apollo's `apollo_control$memorySaver`,
  which computes the analytic gradient in chunks of about two respondents.
  The estimates are bitwise identical; on the Vittel correlated MMNL (46
  parameters, 1000 draws) peak memory fell from 7.5 to 4.1 GB at the same
  speed. The default stays Apollo's (`FALSE`).
- Kept, because Apollo 0.3.5 needs it: on data with availability columns,
  utilities are measured from the chosen alternative's (see 0.9.5).
- `klue_design()` runs idefix's Modfed sequentially (`parallel = FALSE`)
  instead of on `detectCores() - 1` workers; the design is identical.
- `klue_demo(full = TRUE)` estimates its MMNL with 500 MLHS draws, the draw
  count of Apollo's MMNL example.

Additions for inference and re-estimation (5 October 2026, written for the
Markets vs Governments analysis). Default fits are unchanged: the reproduction
checks of 6 October (`dev/check_ladder_carryover.R`,
`dev/check_empirical_lcmnl_maxima.R`) ran on this tree and reproduced the
stored log-likelihoods.

- New `klue_refit()`: re-estimate a fitted MNL or latent class model on new
  data, starting from the fit's estimates matched by label (subsets, one arm
  of an experiment, added or dropped terms, a membership model).
- New `klue_predict()` (choice probabilities of every task row),
  `klue_draw_betas()` (respondent coefficients drawn from a fit, for
  `klue_simulate_choices()`) and `klue_derror()` (MNL D-error of a design,
  as Ngene reports it for one respondent).
- New `klue_nlcom()` (delta-method inference on any function of the
  parameters), with `klue_mean_fn()` (share-weighted mean of a coefficient)
  and `klue_wtp_fn()` (willingness to pay, per class or as the share-weighted
  mean); `klue_par()` gives the parameter name of a labelled coefficient.
- `klue_wald()` takes `fn` (a function of the parameters, tested with a
  generalised inverse) and `tol`; `klue_fieller()` and `klue_tost()` accept
  functions of the parameters as well as weight vectors.
- `klue_dgp()` takes `labels`, `asc_labels` and `common` (coefficients shared
  by all classes).
- `klue_lcmnl()` takes `clean`: fit every start at full precision with its
  covariance and return the best start that has not failed; also `maxit`,
  `reltol` and `cluster_cap`.
- Every fit reports `failed` (the optimiser did not converge, or its
  covariance is singular or has missing entries), `max_abs_par` and
  `entropy_norm` (1 - H / (N ln C)).
- `DESCRIPTION`: the Description field now says what the package does; the
  version history it carried is in this file (0.9.0–0.9.2).

# klue 0.9.5 (2026-09-26)

- `klue_mmnl()` can make any subset of parameters random, constants included.
  New `random` names them among `x1..xN`, `price` and `asc1..ascK`, where
  `ascK` indexes the groups of `dgp$asc_map` (constants now follow `asc_map`,
  so alternatives in one group share one constant); unnamed parameters are
  fixed. New `price = c("lognormal", "normal", "fixed")` sets the distribution
  of a random price. Parameter names: fixed `b_x{a}`, `b_price`, fixed
  constants `asc_alt{g}` under the default map (`asc{g}` otherwise); random
  `mu_{p}` with `sigma_{p}` (independent) or Cholesky rows `s_{p}_{l}`
  (correlated; `s_pr` for price).
- Fits return `par`, `vcov` and `robust_vcov` (Apollo's `varcov` and
  `robvarcov`), plus `random`, `price` and `correlation`, so `klue_lincom()`,
  `klue_wald()` and the other inference helpers run on an MMNL fit directly.
- New `warm_start`: restart from an earlier fit, skipping the stage-1 warm
  start. An independent fit given to a correlated call is mapped to the point
  where the two models coincide, so with the same draws the correlated fit
  cannot end below it; the result's `nested_ok` checks this.
- Fixed: the correlated model's warm start put the independent fit's
  log standard deviations on the Cholesky diagonal, which enters the utility
  linearly (an sd of 0.1 started at 2.3), and reset the constants to zero; the
  neutral start had the same diagonal error. It now starts at the standard
  deviations and carries every estimate over. Correlated fits can change; the
  default independent specification is bit-identical to 0.9.4.
- Fixed: on data with availability columns, MMNL likelihoods could turn NaN.
  Apollo 0.3.5's `apollo_mnl` sets unavailable utilities to 0 and then
  multiplies `exp(V_j - V_chosen)` by availability, which is `Inf * 0 = NaN`
  once a chosen utility falls below about -709 (extreme random-coefficient
  draws); Apollo then refuses to estimate. klue now measures utilities from the
  chosen alternative's, which leaves the probabilities unchanged and makes that
  term `exp(0) * 0 = 0`. Data without availability columns are unaffected.
- If estimation from a `warm_start` fails, `klue_mmnl()` falls back to the
  two-stage path from the same point (`warm_start_fallback = TRUE`).

# klue 0.9.4 (2026-09-21)

- Rank-deficient observed information is now detected and reported. Previously
  only an exactly computationally singular matrix was caught (`solve()`
  erroring); a design with an exact linear dependency leaves an eigenvalue at
  ~0, which `solve()` still inverts into huge negative variances and NaN
  standard errors with `singular = FALSE`. Fits now carry `aliased`, naming the
  parameters that are not separately identified, and set `vcov`/`robust_vcov`
  to `NA` with an explanatory message. Estimates, log-likelihood, class
  assignment and `k` are unchanged, and the covariance of a well-identified fit
  is bit-identical to 0.9.3.
- `klue()` documents its `asc_map` argument (`R CMD check` clean).

# klue 0.9.3 (2026-09-17)

API consolidation: one package, canonical names only.

- The 0.6.x alias exports are gone (`run_*`, `build_*`, `generate_*`,
  `estimate_*multistart*`, `make_dgp_config`); the `klue_*` names are the only
  API. 26 exports.
- The paper's Monte Carlo study drivers (`klue_study_*`) and the legacy
  fixed-design simulator moved out of the package into
  `studies/klue_studies.R` in the reproduction repository; `source()` it after
  loading the package.
- `klue_database_long()` orders respondents and tasks numerically, carries
  `covariate_cols`, and attaches respondent ids as `attr(db, "ids")`.
- New `klue_simulate_choices()`: simulate choices on an existing design.
- Fit objects carry a `singular` flag.
- `cleanup_apollo()` and the three start generators used by the initialisation
  ablations are exported.

# klue 0.9.2 (2026-07-09)

Hardening release. Principle applied throughout: no silent fallbacks — every
degraded path warns, messages, or stops.

- `apollo` moved from Depends to Suggests: only `klue_mmnl()` needs it (hard
  `stop()` with an install hint if absent); `library(klue)` no longer attaches
  apollo, so LCMNL-only pipelines load faster.
- Reproduction switch: `klue_lcmnl()`, `klue_starts()`, and the internal start
  generators gain `screen` (default `getOption("klue.screen", TRUE)`).
  `options(klue.screen = FALSE)` restores the pre-0.9.1 full-precision path
  end to end (tight cluster-MNL start fits, every start estimated at full
  precision, no polish) for exact reproduction of results produced with it.
  Verified against 0.9.0: equivalence harness 15/15 components
  (engine-equivalence mode exact; screened default within the LL/BIC
  contract).
- Deterministic seeding (clustering, subsampling, baseline start generators)
  saves and restores `.Random.seed`: klue no longer resets the caller's RNG
  stream.
- `C = 1` short-circuits to a single fit (the six clustering starts coincide
  there); `best_method` is reported as `"pooled"`.
- A dead `mclapply` worker now raises a warning naming the lost start instead
  of crashing the winner scan; `gc()` before forking shrinks the inherited
  copy-on-write heap.
- `compute_recovery()` refuses K > 8 loudly instead of attempting a K^K
  permutation enumeration.
- Loud fallbacks: cluster-start fallback coefficients, EM collapsed-class and
  failed M-step skips, and all MMNL warm-start fallbacks report what happened
  and what is used instead.
- `estimate_lcmnl_em()` returns `fit$par` like the ML path;
  `apollo_probabilities` cleanup is registered before the global assignment,
  so an apollo error cannot leak it into the global environment.

# klue 0.9.1 (2026-07-07)

Large-N fixes, motivated by an unbounded run at N = 6,000 (55+ minutes in the
O(N^2) clustering-starts phase).

- Clustering starts are computed on a fixed-seed respondent subsample capped
  at `cluster_cap = 2000` (pam, Mclust, hclust are O(N^2) or worse; starting
  values do not need every respondent). The LCMNL itself still fits the full
  sample.
- Two-stage multistart in `klue_lcmnl()`: starts are screened at
  `maxit = 200`, `reltol = 1e-6`, and only the best log-likelihood is
  polished at full precision, warm-started via the new
  `estimate_lcmnl(start_par =)`. Fits that hit the iteration cap keep their
  finite LL instead of being discarded. Tie-breaks between equal-LL starts
  can differ from 0.9.0; the polished optimum is the same.
- Fit objects return the full parameter vector as `fit$par` (betas, per-class
  ASCs, share deltas).
- Each start logs its LL and elapsed time via `message()`.

# klue 0.9.0 (2026-06-15)

Full rewrite of the 0.6.x engine: same exported API (all 0.6.x aliases kept)
at roughly 70% of the line count. One shared MNL likelihood/gradient kernel;
one simulator (`klue_simulate` with `covariates =` and `design =`, RNG order
preserved so seeded datasets are bit-identical to 0.6.x); clustering starts
behind a `klue_starts(method =, feature_type =)` registry; merged
`klue_mmnl(correlation =)`; opt-in parallel multistart (`n_cores`) and
explicit `n_cores` on all study drivers (output identical at any core
count). Verified equivalent to 0.6.3: 16/16 components (data generation
bit-exact; LCMNL LL/BIC to 1e-5; MMNL to 1e-2); testthat suite green.
