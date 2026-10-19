## R CMD check results

0 errors | 0 warnings | 1 note

## Notes explained

### Rd file formatting

- Benign "Lost braces" note in `dot-euk_heuristic_distance_impl.Rd` and
  `euk_bayesian_distance.Rd` (escaped backslashes inside Rd markup).
  Cosmetic only; does not affect documentation rendering.

## Attribution

This package is a reimplementation of algorithms from
Joel-Barratt/Eukaryotyping (GPL >= 3). Provenance is documented in
`inst/COPYRIGHTS` and the method is cited via doi:10.1017/S0031182019000581.
The original authors have been notified and will be added to `Authors@R`
only with their consent.

## Downstream dependencies

None (new submission).

## Test environments

- Local: R 4.5.2, Ubuntu 26.04, gcc 15.2.0

## History

- 0.0.0.9000: Initial development submission
- 0.1.0:    Release; C++ default, _r deprecated, vignettes, CRAN prep