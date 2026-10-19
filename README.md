# fastEukaryotyping

Fast C++17 reimplementation of the **Eukaryotyping** ensemble for MLST haplotype data sheets.

This package migrates the original R scripts of Barratt, Plucinski and colleagues
(Barratt et al. 2019, *Parasitology*, [doi:10.1017/S0031182019000581](https://doi.org/10.1017/S0031182019000581))
into a CRAN-quality R package with Rcpp/RcppParallel backends.

## Features

- `euk_heuristic_distance()` — Barratt's unsupervised heuristic mixture-model distance.
- `euk_bayesian_distance()` — Plucinski's naive Bayes genetic distance.
- `euk_ensemble_distance()` — Rank-mapped mean of the two matrices.
- Handles missing data, heterozygosity, ploidy differences and mixed-genotype infections.
- Pure-R reference implementations (`*_r`) are kept for cross-language testing.

## Installation

```r
remotes::install_github("alrobles/fastEukaryotyping")
```

## Usage

```r
library(fastEukaryotyping)

# geno is a list of n x maxMOI integer matrices; NA for empty slots
# ploidy is an integer vector (1 = haploid, 2 = diploid) for each locus
d <- euk_heuristic_distance(geno, ploidy)
```

## Verification against the original Barratt script

The package ships a filtered HDS fixture taken from
`Joel-Barratt/Eukaryotyping` (68 samples, 14 loci) and the original
`Heuristic_pairwisedistancematrix` saved as
`tests/testthat/fixtures/barratt_expected.rds`.

```r
library(fastEukaryotyping)

geno   <- readRDS(system.file("testthat/fixtures", "barratt_geno.rds", package = "fastEukaryotyping"))
ploidy <- readRDS(system.file("testthat/fixtures", "barratt_ploidy.rds", package = "fastEukaryotyping"))
ids    <- readRDS(system.file("testthat/fixtures", "barratt_ids.rds", package = "fastEukaryotyping"))
expected <- readRDS(system.file("testthat/fixtures", "barratt_expected.rds", package = "fastEukaryotyping"))

D <- euk_heuristic_distance(geno, ploidy, ids = ids, impute_data = TRUE)
all.equal(D, expected, tolerance = 1e-10)  # TRUE
```

`R CMD check --as-cran` also runs `testthat::test()` which enforces this
identity for both the R reference and the C++ implementation.

## Development

The package follows the `_r` → `_cpp` migration pattern:

1. `R/*_r.R` — base-R reference implementations.
2. `src/*.cpp` — C++17/RcppParallel implementations.
3. `tests/testthat/test-*-cpp-vs-r.R` — numerical parity tests.

## License

GPL (>= 3). Original Eukaryotyping scripts are credited in `DESCRIPTION`.
