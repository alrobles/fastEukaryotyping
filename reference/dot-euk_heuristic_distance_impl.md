# Internal heuristic distance (C++ implementation)

RcppParallel reimplementation of \codeeuk_heuristic_distance_r(). Each
locus may have at most 64 distinct alleles.

## Usage

``` r
.euk_heuristic_distance_impl(geno, ploidy, m, impute_data)
```

## Arguments

- geno:

  List of integer matrices.

- ploidy:

  Integer vector of ploidy values.

- m:

  Integer vector of m values.

- impute_data:

  Integer (0/1). If 1, apply the two-pass Barratt missing-locus
  imputation.

## Value

A symmetric distance matrix.
