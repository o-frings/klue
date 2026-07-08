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
