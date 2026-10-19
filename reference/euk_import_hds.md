# Import a Haplotype Data Sheet (HDS) into the internal geno list

Import a Haplotype Data Sheet (HDS) into the internal geno list

## Usage

``` r
euk_import_hds(
  data,
  locinames,
  locinames_base = locinames,
  ploidy = rep(2L, length(locinames)),
  include_seq_id = TRUE
)
```

## Arguments

- data:

  A \`data.frame\` with a \`Seq_ID\` column and one column per haplotype
  per locus. Presence is marked by \`"X"\`.

- locinames:

  Character vector of marker names (as used in column names).

- locinames_base:

  Optional character vector of base marker names for markers split into
  sub-markers. Defaults to \`locinames\`.

- ploidy:

  Integer vector (1 = haploid, 2 = diploid) for each locus in
  \`locinames\`.

- include_seq_id:

  Logical. Should the \`Seq_ID\` column be kept in the output? Default
  \`TRUE\`.

## Value

A \`data.frame\` with \`Seq_ID\` and columns named
\`\<slot\>\_Hap\_\<lociname\>\`, plus a \`geno\` attribute containing
the list of per-locus integer matrices.

## Examples

``` r
# See inst/extdata/Example_haplotype_data_sheet.txt
```
