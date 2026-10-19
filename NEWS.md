# fastEukaryotyping 0.1.0

- Initial package skeleton.
- Pure-R reference implementation of `euk_heuristic_distance_r()`.
- C++17/RcppParallel implementation of `euk_heuristic_distance_impl()`.
- Added `impute_data` flag to `euk_heuristic_distance()` and reference implementation.
- Two-pass missing-locus imputation ported from `Joel-Barratt/Eukaryotyping`.
- `euk_heuristic_distance()` matches the original Barratt heuristic output on the
  filtered 68-sample HDS fixture (`barratt_expected.rds`) to within 1e-10.
- Pure-R reference `euk_bayesian_distance_r()` porting Plucinski's
  `euk_bayesian_fulldataset_V2.r`, including epsilon fallbacks and the
  multiset/unique-overlap pair semantics of the original.
- C++17/RcppParallel implementation `euk_bayesian_distance_impl()`.
- `euk_bayesian_distance()` matches the original script output on the same
  68-sample fixture to within 1e-10 (R and C++ backends).
- `euk_ensemble_distance()`: rank-mapped mean of the heuristic and Bayesian
  matrices, per the original ensemble workflow.
- Added `inst/COPYRIGHTS` and `cran-comments.md`; `R CMD check --as-cran`
  reports 0 errors, 0 warnings and 2 expected notes (new submission,
  system compiler flag).
- Original authors moved out of `Authors@R` pending their consent;
  attribution retained in `inst/COPYRIGHTS` and the method citation.

### Migration: _cpp suffix removed, _r functions deprecated

- Internal C++ functions renamed from `.euk_*_distance_cpp` to
  `.euk_*_distance_impl` (the `_cpp` suffix is no longer part of the
  release API).
- `use_r = FALSE` is now the default; the C++ backend is the primary
  implementation path. `use_r = TRUE` emits a deprecation warning.
- Pure-R reference functions (`euk_heuristic_distance_r`,
  `euk_bayesian_distance_r`, `euk_ensemble_distance_r`) are marked
  deprecated and emit a deprecation warning when called. They remain
  available for cross-validation and development.

### CRAN and documentation

- Added `_pkgdown.yml` for pkgdown site generation.
- Added `LICENSE` (GPL >= 3).
- Added `inst/CITATION` for `citation("fasteukaryotyping")`.
- Added vignettes: `introduction.Rmd` (package overview and methods)
  and `mlst-typing-workflow.Rmd` (end-to-end analysis workflow).
- Updated man/ Rd files via roxygen2.
- R CMD check: 0 errors, 0 warnings, 1 note (Rd braces only).
- 16 tests pass, covering fixture equivalence, random-data parity
  between R and C++, symmetry, and range validation.