# Eukaryotyping ensemble distance

Rank-mapped mean of the heuristic and Bayesian distances, following
Barratt et al. (2019): each matrix is converted to fractional ranks over
the unordered pairs and the ensemble entry is the mean of the two ranks.

## Usage

``` r
euk_ensemble_distance(geno, ploidy, ids = NULL, epsilon = 0.3072, ...)
```

## Arguments

- geno:

  List of \`n x maxMOI\` integer matrices (one per locus).

- ploidy:

  Integer vector (1 = haploid, 2 = diploid) for each locus.

- ids:

  Optional character vector of sample identifiers.

- epsilon:

  Rate of missed alleles for the Bayesian component.

- ...:

  Passed to \[euk_heuristic_distance()\] (e.g. \`m\`, \`impute_data\`,
  \`use_r\`).

## Value

An \`n x n\` matrix of mean fractional ranks in \\\[0, 1\]\\.
