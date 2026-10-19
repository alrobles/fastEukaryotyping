# fastEukaryotyping

## fastEukaryotyping 0.1.0 ($today)

- Initial package skeleton.
- Pure-R reference implementation of `euk_heuristic_distance_r()`.
- C++17/RcppParallel implementation of `euk_heuristic_distance_cpp()`.
- Added `impute_data` flag to `euk_heuristic_distance()` and reference implementation.
- Two-pass missing-locus imputation ported from `Joel-Barratt/Eukaryotyping`.
- `euk_heuristic_distance()` matches the original Barratt heuristic output on the
  filtered 68-sample HDS fixture (`barratt_expected.rds`) to within 1e-10.
- Pure-R reference `euk_bayesian_distance_r()` porting Plucinski's
  `euk_bayesian_fulldataset_V2.r`, including epsilon fallbacks and the
  multiset/unique-overlap pair semantics of the original.
- C++17/RcppParallel implementation `.euk_bayesian_distance_cpp()`.
- `euk_bayesian_distance()` matches the original script output on the same
  68-sample fixture to within 1e-10 (R and C++ backends).
- `euk_ensemble_distance()`: rank-mapped mean of the heuristic and Bayesian
  matrices, per the original ensemble workflow.
- Added `inst/COPYRIGHTS` and `cran-comments.md`; `R CMD check --as-cran`
  reports 0 errors, 0 warnings and 2 expected notes (new submission,
  system compiler flag).
