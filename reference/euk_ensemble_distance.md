# Eukaryotyping ensemble distance

Combines the heuristic and Bayesian distance matrices into a single
ensemble. Two rules are available: \`"original"\` (default) reproduces
the distribution-mapping procedure of the original \`run.r\` script
shipped with Barratt et al. (2019); \`"rank"\` is a package-specific
direction-consistent fractional-rank mean.

## Usage

``` r
euk_ensemble_distance(
  geno,
  ploidy,
  ids = NULL,
  method = c("original", "rank"),
  seed = NULL,
  epsilon = 0.3072,
  ...
)
```

## Arguments

- geno:

  List of \`n x maxMOI\` integer matrices (one per locus).

- ploidy:

  Integer vector (1 = haploid, 2 = diploid) for each locus.

- ids:

  Optional character vector of sample identifiers.

- method:

  Combination rule. \`"original"\` (default) reproduces the combining
  procedure of the original \`run.r\` script: the Bayesian matrix is
  reversed (\`2 - D_bayesian\`), both matrices are min-max normalised
  over all entries including the diagonal, the sorted Bayesian values
  are mapped onto the ranks of the heuristic matrix, the ensemble is the
  elementwise mean of the mapped and heuristic matrices, rescaled by its
  maximum. Random tie breaking makes this method non-deterministic
  without \`seed\` and can yield a non-symmetric matrix. \`"rank"\` is a
  package-specific fractional-rank mean over the strictly upper
  triangle, with the Bayesian matrix reversed so both components are
  divergence-oriented; its diagonal is \`NA\`.

- seed:

  Optional integer seed for the \`"original"\` method's random tie
  breaking. Ignored for \`method = "rank"\`.

- epsilon:

  Rate of missed alleles for the Bayesian component.

- ...:

  Passed to \[euk_heuristic_distance()\] (e.g. \`m\`, \`impute_data\`,
  \`use_r\`).

## Value

An \`n x n\` matrix of ensemble distances (in \\\[0, 1\]\\ for \`method
= "rank"\` off-diagonal entries). Larger values indicate pairs that both
components rank as more divergent.
