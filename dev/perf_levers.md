# klue performance levers — measured (2026-06-14)

Machine: 12 logical / 6 physical cores (Apple Silicon), R 4.3.1, reference BLAS.
Benchmarks are sequential (CPU timing can't be parallelised without corrupting
it). Numbers are for the standard DGP (N=300-450, T=20, n_beta=5, J=3).

## Profile: where time goes

- `klue_lcmnl` C=3 (one dataset): 3.1 s total = 0.38 s clustering starts +
  2.73 s for the six BFGS optimisations (88%).
- LCMNL study condition (BIC scan C=1:5): ~5 s.
- MMNL study condition (Apollo, 3000 draws): ~3 min — i.e. ~35x an LCMNL
  condition. **The MMNL study arms dominate the full study's wall-clock.**

## Lever scoreboard

| Lever | Speedup | Risk | Status |
|---|---|---|---|
| A. Raise study core cap 4 -> ~physical | ~1.45x (LCMNL studies) | none (output invariant to core count) | **DONE** |
| B. Parallel multistart `n_cores` | 2.8x (interactive C-sweep) | none (exact; default off) | **DONE** |
| C. Accelerate BLAS (vecLib) | LCMNL ~1x; Apollo MMNL TBD | reproducibility (multi-thread FP) + sudo | proposed, untested |
| D. Parallel MMNL across conditions | 0.9x (SLOWER) | n/a | rejected (memory-bound; see below) |
| E. Relax optim reltol 1e-10 -> 1e-8 | ~1.2-1.4x per optim | changes paper numbers slightly | not recommended for paper |
| F. RcppArmadillo kernel | 0.8x (SLOWER) | n/a | rejected (see cpp_kernel_findings.md) |

## Detail

### A — study core cap (DONE)
`.klue_cores()` capped at 4; raised default to `min(16, cores-2)`. Scaling on
24 LCMNL conditions: 2c 221 s, 4c 120 s, 8c 83 s, 10c 81 s. Sweet spot ~physical
cores; >physical gives diminishing returns. Output byte-identical at 4 vs 10
cores (per-condition seeds make mclapply order-independent).

Sequencing vs parallelisation is now an explicit per-call knob: every study
driver and the master `klue_study()` take `n_cores` (default `.klue_cores()`).
`n_cores = 1` runs fully sequentially (lowest memory, reproducible timing,
shared-machine-friendly); higher parallelises across conditions. Same output
either way. The global `options(klue.cores = n)` still sets the default. Choose
per run by time budget: e.g. `klue_study(n_cores = 1)` overnight on a laptop
you are using, `klue_study(n_cores = 10)` when the machine is free.

### B — parallel multistart (DONE)
`klue_lcmnl(..., n_cores=)` / `klue(..., n_cores=)` fork the six clustering
starts. Exact (dLL=0, dbetas=0, same method). Full C=1:5 sweep 19.6 s -> 7.1 s.
Default 1 so it never nests with the study drivers' across-condition mclapply.

### C — Accelerate BLAS (to test)
R currently uses the reference BLAS (2000x2000 matmul = 2.5 s; Accelerate
~0.1 s, ~25x). BUT the LCMNL hot path is skinny (T x 5) gemv, memory-bound, so
Accelerate barely helps it (and is why the C++ port didn't help either).
Apollo MMNL does larger draw x obs x param dense algebra and MAY benefit; this
is the cheapest thing to try for the MMNL-dominated study. Enable (reversible):
```sh
cd /Library/Frameworks/R.framework/Resources/lib
sudo ln -sf libRblas.vecLib.dylib libRblas.dylib   # revert: ln -sf libRblas.0.dylib libRblas.dylib
```
Caveat: Accelerate is multi-threaded -> (i) competes with mclapply in the study
drivers (set klue.cores lower if both on); (ii) floating-point results can
differ at ~1e-10 and are not bit-reproducible across runs, which matters for
the paper's exact-reproduction claim. Recommend: pin paper reproduction to
reference BLAS; offer Accelerate as an opt-in accelerator for interactive use.

### D — parallel MMNL across conditions (REJECTED, tested)
Tested: 4 MMNL conditions (1000 draws) sequential with Apollo nCores=4 vs
parallel via `mclapply(mc.cores=4)` with Apollo nCores=1. Result: parallel was
SLOWER (598 s -> 691 s, 0.9x) and bit-identical (max |dLL| = 0, all converged).
Forking is numerically safe (each child has its own globalenv, so the
"Apollo global state" serialization is not actually required), but each MMNL
fit allocates large draw x obs matrices; running several at once saturates
memory bandwidth and loses more than the concurrency gains. Apollo's own
internal multicore already parallelises each fit cache-efficiently within one
process. Conclusion: leave the study MMNL phase as-is (sequential conditions,
Apollo nCores internal). The only remaining MMNL lever is the BLAS (C).

Apollo internal nCores scaling, one MMNL fit (1000 draws, LL identical across
all): 1c 354 s, 2c 211 s (1.7x), 4c 140 s (2.5x), 6c 128 s (2.8x). Sub-linear,
flat past ~physical cores. The default `.klue_default_mmnl_cores()` =
`detectCores(logical=FALSE) - 1` (= 11 here) is already at/above diminishing
returns; trimming to ~physical-core count would save scheduler overhead at no
speed cost but is results-neutral (LL is nCores-invariant), so not changed.
