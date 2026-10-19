# Heuristic Eukaryotyping distance

Computes Barratt's unsupervised heuristic distance matrix for MLST
haplotype data sheets.

## Usage

``` r
euk_heuristic_distance(
  geno,
  ploidy,
  ids = NULL,
  m = rep(1L, length(geno)),
  impute_data = TRUE,
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

- m:

  Integer vector of \`m\` values, one per locus. Default \`rep(1L,
  length(geno))\`.

- impute_data:

  Logical. If \`TRUE\` (default), missing loci are imputed with the
  two-pass Barratt mean procedure so that the returned matrix is fully
  populated. If \`FALSE\`, pairs sharing a missing locus return \`NA\`.

- use_r:

  Logical. If \`TRUE\` use the pure-R reference implementation
  (deprecated); if \`FALSE\` (default) use the optimised C++
  implementation.

## Value

An \`n x n\` numeric distance matrix.

## Examples

``` r
set.seed(42)
n <- 10
geno <- list(
  matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
  matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
)
d <- euk_heuristic_distance(geno, ploidy = c(2L, 2L))
```
