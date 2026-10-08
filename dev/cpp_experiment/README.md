# Parked: RcppArmadillo LCMNL kernel (negative result)

This folder holds the C++ kernel experiment. It is **not** part of the shipped
`klue` package: benchmarks showed it is ~15–30% *slower* than the pure-R
kernel, because that kernel is already BLAS-bound vectorised algebra with no
interpreter loop for C++ to remove (full write-up: `../cpp_kernel_findings.md`).
The genuine speed lever — opt-in parallel multistart (`n_cores`) — was adopted
in the package instead.

The kernel is kept because it is numerically exact (verified to 1e-11 at the
kernel level, 102/102 configs end-to-end) and is a useful independent
cross-check of the R estimator maths.

Files:
- `lcmnl.cpp` — the kernel (nll/grad/weighted-MNL/LCMNL/panel-loglik).
- `test-engine.R` — testthat case asserting cpp == r to 1e-8.
- `engine_grid.R` — exhaustive cpp-vs-r grid across DGP shapes / C / estimator.

## Restoring it for a cross-check run

The machine's Command Line Tools are broken (a stale libc++ stub shadows the
SDK's); compile with the SDK's libc++ on `CPATH`:

```sh
SDK_CXX=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/c++/v1
cp dev/cpp_experiment/lcmnl.cpp klue/src/lcmnl.cpp        # mkdir klue/src first
# re-add to klue/DESCRIPTION:  Imports: Rcpp ;  LinkingTo: Rcpp, RcppArmadillo
# re-add to klue/NAMESPACE:    useDynLib(klue, .registration = TRUE)
#                               importFrom(Rcpp, evalCpp)
# re-add the engine = c("cpp","r") dispatch to klue/R/03_estimate.R (see git
#   history of the commit that added it), then:
Rscript -e 'Rcpp::compileAttributes("klue")'
CPATH=$SDK_CXX Rscript -e 'devtools::load_all("klue"); devtools::test("klue")'
CPATH=$SDK_CXX Rscript dev/cpp_experiment/engine_grid.R
```

Proper fix for the toolchain (recommended regardless):
`sudo rm -rf /Library/Developer/CommandLineTools && sudo xcode-select --install`
