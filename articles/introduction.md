# Introduction to fasteukaryotyping

## What is MLST and why distance matrices?

Multilocus Sequence Typing (MLST) characterises microbial isolates by
the alleles present at a set of housekeeping loci. Given a collection of
samples (isolates), each with potentially multiple alleles per locus
(due to heterozygosity, mixed infections, or polyploidy), the central
question is: *how genetically similar are any two samples?*

`fasteukaryotyping` implements the three-method ensemble described by
Barratt et al. (2019) (<doi:10.1017/S0031182019000581>) for computing
pairwise genetic distances from MLST haplotype data sheets (HDS).

``` r

library(fasteukaryotyping)
```

## The three distance methods

### Heuristic distance (`euk_heuristic_distance`)

Barratt’s unsupervised heuristic mixture-model distance computes, for
each locus, a delta value that depends on the number of shared alleles,
the allele-set sizes, a per-locus entropy term \\H\_{\nu}\\, and a
frequency-informed penalty \\P\\. The per-locus deltas are summed across
loci to produce the pairwise distance. Missing loci are imputed with a
two-pass mean procedure that matches the original Barratt script.

``` r

set.seed(42)
n <- 10
geno <- list(
  matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
  matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
)
d_h <- euk_heuristic_distance(geno, ploidy = c(2L, 2L))
dim(d_h)
```

    ## [1] 10 10

### Bayesian distance (`euk_bayesian_distance`)

Plucinski’s naive Bayes genetic distance evaluates, for each pair and
locus, three log-likelihoods:

- **ll0**: samples are unrelated – every allele drawn independently from
  the locus allele frequencies.
- **ll1**: samples share exactly one allele – the shared allele is the
  most likely one; remaining alleles drawn independently.
- **ll2** (diploid loci only): samples share both alleles – evaluated
  over all 2-subsets; missing allele draws get probability \\\epsilon
  \cdot \min(p_1, p_2)\\.

The per-locus values are summed, normalised to a posterior over the
three hypotheses, and the distance is the posterior expectation in
\\\[0, 2\]\\.

``` r

d_b <- euk_bayesian_distance(geno, ploidy = c(2L, 2L))
dim(d_b)
```

    ## [1] 10 10

### Ensemble distance (`euk_ensemble_distance`)

The ensemble converts each matrix to fractional ranks over the unordered
pairs and takes the mean of the two ranks. This copula-style average is
invariant to monotone transformations of either distance, producing a
consensus value in \\\[0, 1\]\\.

``` r

d_e <- euk_ensemble_distance(geno, ploidy = c(2L, 2L))
dim(d_e)
```

    ## [1] 10 10

## Input data format

The core input `geno` is a list of \\L\\ integer matrices (one per
locus), each with \\n\\ rows (samples) and up to `maxMOI` columns.
Missing alleles are coded as `NA` (or `0`). The ploidy vector specifies
whether each locus is haploid (`1`) or diploid (`2`).

MLST data is often stored as a *Haplotype Data Sheet* (HDS): a CSV with
a `Seq_ID` column and one column per haplotype per locus, where presence
is marked by `"X"`. Use
[`euk_import_hds()`](https://alrobles.github.io/fasteukaryotyping/reference/euk_import_hds.md)
to convert this format into the internal `geno` list:

``` r

hds <- read.csv(system.file("extdata", "example_hds.csv",
                             package = "fasteukaryotyping"),
                stringsAsFactors = FALSE)
out <- euk_import_hds(hds, locinames = "Hap", ploidy = 2L)
attr(out, "geno")
```

    ## [[1]]
    ##       1_Hap_Hap 2_Hap_Hap
    ##  [1,]         1         1
    ##  [2,]         2         2
    ##  [3,]         1         1
    ##  [4,]         2         2
    ##  [5,]         1         1
    ##  [6,]         1         1
    ##  [7,]         2         2
    ##  [8,]         2         2
    ##  [9,]         1         1
    ## [10,]         2         2
    ## [11,]         1         1
    ## [12,]         2         2

## Performance

The C++ implementation uses RcppParallel for multi-threaded pairwise
computation. For a dataset of \\n\\ samples and \\L\\ loci, the
heuristic distance processes \\\binom{n}{2}\\ pairs in parallel per
locus, while the Bayesian distance parallelises over the full \\n \times
n\\ matrix.

Compare the R reference and C++ implementations on the built-in example:

``` r

eg <- example_eukaryotyping
system.time({
  d_r <- euk_heuristic_distance(eg$geno, eg$ploidy, use_r = TRUE)
})
```

    ##    user  system elapsed 
    ##    0.02    0.00    0.02

``` r

system.time({
  d_cpp <- euk_heuristic_distance(eg$geno, eg$ploidy, use_r = FALSE)
})
```

    ##    user  system elapsed 
    ##       0       0       0
