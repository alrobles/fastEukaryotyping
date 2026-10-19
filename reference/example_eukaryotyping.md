# Example 8-locus haplotype data for fasteukaryotyping

A small synthetic dataset used to test and demonstrate
[`euk_heuristic_distance`](https://alrobles.github.io/fastEukaryotyping/reference/euk_heuristic_distance.md)
and
[`euk_import_hds`](https://alrobles.github.io/fastEukaryotyping/reference/euk_import_hds.md).

## Usage

``` r
example_eukaryotyping
```

## Format

A named list with 5 elements:

- `geno`:

  A list of 3 integer matrices (one per locus) with 12 rows (samples)
  and 1-2 columns (ploidy). Missing values are coded as `NA`.

- `ploidy`:

  Integer vector of length 3 (c(2, 2, 1)).

- `ids`:

  Character vector of sample identifiers.

- `D`:

  The 12 x 12 heuristic distance matrix computed from
  `euk_heuristic_distance_r`.

- `hds`:

  A data frame with an X-marked Haplotype Data Sheet for
  `euk_import_hds`.

## Source

Synthetic data created for the fasteukaryotyping test suite.

## Details

The object was generated with `set.seed(42)` so that tests and examples
are deterministic. Use `data(example_eukaryotyping)` to load the list
into the working environment.
