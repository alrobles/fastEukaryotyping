# Internal Bayesian distance (C++ implementation)

RcppParallel reimplementation of `euk_bayesian_distance_r()`
(Plucinski's unsupervised naive Bayes classifier).

## Usage

``` r
.euk_bayesian_distance_impl(geno, ploidy, epsilon)
```

## Arguments

- geno:

  List of integer matrices.

- ploidy:

  Integer vector of ploidy values.

- epsilon:

  Rate of missed alleles.

## Value

A symmetric distance matrix with values in \\\[0, 2\]\\.
