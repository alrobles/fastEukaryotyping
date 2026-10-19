# Bayesian Eukaryotyping distance

Computes Plucinski's unsupervised naive-Bayes distance matrix for MLST
haplotype data sheets.

## Usage

``` r
euk_bayesian_distance(
  geno,
  ploidy,
  ids = NULL,
  epsilon = 0.3072,
  use_r = FALSE
)
```

## Arguments

- geno:

  List of \`n x maxMOI\` integer matrices (one per locus). \`NA\` marks
  empty slots.

- ploidy:

  Integer vector (1 = haploid, 2 = diploid) for each locus.

- ids:

  Optional character vector of sample identifiers.

- epsilon:

  Rate of missed alleles (default \`0.3072\`).

- use_r:

  Logical. If \`TRUE\` use the pure-R reference implementation
  (deprecated); if \`FALSE\` (default) use the optimised C++
  implementation.

## Value

An \`n x n\` numeric distance matrix in \eqn\[0, 2\].

## Examples

``` r
set.seed(42)
n <- 10
geno <- list(
  matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
  matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
)
d <- euk_bayesian_distance(geno, ploidy = c(2L, 2L))
```
