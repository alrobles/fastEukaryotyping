## R CMD check results

0 errors | 0 warnings | 2 notes

## Notes explained

### CRAN incoming feasibility

- This is a **new submission** of a new package.
- `URL`/`BugReports` point to `github.com/alrobles/fastEukaryotyping`, the
  public release repository; it returns 404 until the repository is made
  public at release time.
- Version `0.0.0.9000` triggers the "large components" flag; it will be
  bumped to `0.1.0` for the release tarball.

### Compilation flag

`-mno-omit-leaf-frame-pointer` is injected by this system's R Makeconf
(R 4.5.2, Ubuntu 26.04), not by the package `Makevars`; the package only
requests `CXX_STD = CXX17` plus RcppParallel flags.

## Attribution

This package is a reimplementation of algorithms from
Joel-Barratt/Eukaryotyping (GPL >= 3). Original authors are credited as
contributors in DESCRIPTION and `inst/COPYRIGHTS`; the method is cited via
doi:10.1017/S0031182019000581.

## Downstream dependencies

None (new submission).

## Test environments

- Local: R 4.5.2, Ubuntu 26.04, gcc 15.2.0
