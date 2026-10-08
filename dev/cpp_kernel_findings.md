# C++ (RcppArmadillo) kernel for klue — findings

Date: 2026-06-14. Built and verified the RcppArmadillo LCMNL kernel from the
rebuild plan (`klue/src/lcmnl.cpp`, wired behind `engine = c("cpp","r")`).

## Headline: C++ does NOT speed up this workload

| Workload | R engine | C++ engine | ratio |
|---|---|---|---|
| `estimate_lcmnl` C=3, T_total=9000, ×10 | 3.8 s | 5.0 s | **0.8× (slower)** |
| single C=3 fit, T_total=48000, ×5 | 12.5 s | 14.6 s | **0.85× (slower)** |
| `klue_lcmnl` C=3 multistart | 3.1 s | 3.6 s | **0.86× (slower)** |

The C++ is consistently ~15–30% **slower**. Reason: the R kernel
(`.mnl_eval`/`.mnl_score`) is already fully vectorised dense linear algebra —
the per-iteration cost is BLAS matrix products plus a vectorised `exp`, which
R dispatches to the same optimised BLAS the C++ uses. There is no R-interpreter
hot loop to remove, so RcppArmadillo only adds per-call marshalling overhead.
This is the expected outcome when porting already-vectorised R: Rcpp wins
against interpreted R *loops*, not against BLAS-bound array code.

Profiling `klue_lcmnl` C=3: 95% of time is the `.Call` to the kernel; of the
3.1 s total, clustering starts are 0.38 s and the six optimisations are 2.73 s.

## The kernel is numerically correct (verified, kept as opt-in)

- Fixed-parameter kernel agreement R vs C++: **≤ 1e-11** (nll, grad, weighted
  MNL, LCMNL C=2/3/4, panel log-lik).
- Exhaustive grid `dev/engine_grid.R`: **102/102 configs equivalent**
  across n_generic∈{2,4}, J∈{2,3,4}, true_K∈{2,3}, C∈1:4, {ml,em}, plus
  unbalanced multistart. Criterion: log-likelihood (the optimised objective)
  agrees in every config (worst dLL = 1.4e-5); class coefficients/posteriors
  agree to ~1e-8 wherever classes are identified (C ≤ true_K). On
  over-specified C the empty classes are non-identified — their coefficients
  differ between engines, but also between any two R runs from a 1e-9-nudged
  start (verified), so this is benign, not a kernel error.
- Adversarial static review of the C++ against the R reference: no math bug;
  one nit (the no-copy `panel_sum` view requires T_total divisible by N, which
  the invariant guarantees).
- testthat `test-engine.R` (52 assertions) locks this in; skips when the
  compiled code is absent.

The default engine is therefore **"r"**; "cpp" is an opt-in
(`options(klue.engine="cpp")`) independent cross-check, not a faster path.

## The genuine speed lever: parallel multistart

The six clustering starts inside `klue_lcmnl` run sequentially and are 88% of
the time. Parallelising them over forks:

| | sequential | parallel (6 cores, 12 available) |
|---|---|---|
| 6 starts, C=3 | 2.58 s | 0.69 s — **3.8×** |

Caveat: the study drivers already parallelise across conditions via
`mclapply`, so parallel-starts must be **opt-in / default-off** inside
`klue_lcmnl` (a `n_cores`/`parallel` argument) to avoid nested
oversubscription. It is the real win for interactive single-dataset use
(`klue()` / `klue_lcmnl()` on one dataset).

## Toolchain issue found on this machine (separate from klue)

The package's `src/` will not compile until the Command Line Tools are fixed.
Root cause: `/Library/Developer/CommandLineTools/usr/include/c++/v1/` is a
stale 11-file libc++ stub (Oct 2023) that sits first in clang's search path
and shadows the complete 193-file libc++ in the SDK, so `<cmath>` is not found.

- Workaround used for all builds here (local to the command, no global change):
  `CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1`
- Proper fix (recommended): reinstall CLT —
  `sudo rm -rf /Library/Developer/CommandLineTools && sudo xcode-select --install`

This breakage affects ALL R source-package installs on the machine, not just
klue. r-universe/CRAN build servers are unaffected.

## Decision taken (2026-06-14)

Parked the C++ and implemented parallel multistart:

- The C++ kernel and its tests moved to `dev/cpp_experiment/` (with a README on
  restoring it for a cross-check). The package is pure R again — no `src/`, no
  `Rcpp`/`RcppArmadillo` dependency, no compiler needed to install.
- `klue_lcmnl(..., n_cores = 1L)` and `klue(..., n_cores = 1L)` gained an
  opt-in parallel path: with `n_cores > 1` the six per-start fits run via
  `parallel::mclapply`. Default 1 (sequential) so the study drivers'
  across-condition `mclapply` is not nested/oversubscribed.
- Verified: parallel == sequential **exactly** (dLL = 0, dbetas = 0, identical
  best method) for C = 2,3,4; full C=1:5 sweep **19.6 s -> 7.1 s (2.8x)** on 6
  of 12 cores. Default path byte-identical to old klue (batch checks pass);
  57-test suite green.
