# How the Eukaryotyping distance ensemble works

This vignette explains the three distance methods that make up the
Eukaryotyping ensemble. The component methods (heuristic and Bayesian)
are faithful reimplementations of the original R scripts from Barratt,
Park, Nascimento et al. (2019;
[doi:10.1017/S0031182019000581](https://doi.org/10.1017/S0031182019000581)),
verified to within \\10^{-10}\\ on a 68-sample *C. cayetanensis* HDS
fixture. The ensemble `"original"` mode reproduces the published `run.r`
workflow verbatim; the `"rank"` mode is a package-specific variant.

## Notation

Let:

- \\n\\ be the number of isolates (samples)
- \\L\\ be the number of loci
- For locus \\j\\, \\\mathcal{A}\_j\\ is the set of observed alleles
  across all samples
- For isolate \\i\\ and locus \\j\\, \\\mathbf{g}\_{ij}\\ is the
  multiset of alleles carried by that isolate (size up to
  \\\text{MOI}\_j\\, the maximum multiplicity for that locus)
- \\p\_{ja}\\ is the empirical frequency of allele \\a\\ at locus \\j\\
- \\H\_{\nu j}\\ is the entropy of locus \\j\\

## 1. Heuristic distance (`euk_heuristic_distance`)

The heuristic distance is Barratt’s unsupervised mixture-model distance.
For each pair of isolates \\(i, k)\\ and each locus \\j\\, a per-locus
delta \\\delta\_{ik}^{(j)}\\ is computed; the overall distance is the
sum of per-locus deltas. A contribution is missing (no comparable locus)
when **either** isolate has no allele at that locus, not only when both
lack it.

### Per-locus delta

Let \\S\_{ij}\\ be the set of distinct alleles carried by isolate \\i\\
at locus \\j\\. Define:

- \\x = \|S\_{ij}\| + \|S\_{kj}\|\\ — total distinct allele count
- \\y = \|S\_{ij} \cap S\_{kj}\|\\ — number of shared distinct alleles
- \\n = \min(\|S\_{ij}\|, \|S\_{kj}\|)\\
- \\m_j\\ — a per-locus tuning parameter (default 1)

A frequency-informed penalty \\P\\ moderates the contribution of shared
alleles that are rare:

\\ P = \begin{cases} \Bigl(\frac{\\\\r : S\_{ij} \cap S\_{kj} \subseteq
S\_{rj}\\}{n\_{\text{nonmissing},j}}\Bigr)^2 & y \> 0 \\\[4pt\] 1 & y =
0 \end{cases} \\

where \\n\_{\text{nonmissing},j}\\ is the number of isolates with at
least one non-missing allele at locus \\j\\.

The per-locus delta for **diploid** loci (\\\text{ploidy} \> 1\\) is:

\\ \delta = H\_{\nu j} \cdot \Bigl( \delta\_{\text{raw}} \cdot
\mathbb{I}(\delta\_{\text{raw}} \> 0) + P \cdot
\mathbb{I}(\delta\_{\text{raw}} = 0) \Bigr) \cdot k \\

where \\k = 1\\ if \\y = 0\\ and \\k = P\\ otherwise, and
\\\delta\_{\text{raw}}\\ is a piecewise integer function of \\x\\,
\\y\\, \\n\\, and \\m_j\\ that implements the original Barratt lookup
table (see `src/eukaryotyping.cpp:sub_per_locus()` for the exact
mapping).

For **haploid** loci the expression retains the same structure:

\\ \delta = H\_{\nu j} \cdot \bigl( \delta\_{\text{ex,raw}} \cdot
\mathbb{I}(\delta\_{\text{ex,raw}} \> 0) + P \cdot
\mathbb{I}(\delta\_{\text{ex,raw}} = 0) \bigr) \cdot k \\

with \\\delta\_{\text{ex,raw}} = 2x\\ when \\y = 0\\ and
\\\delta\_{\text{ex,raw}} = 0\\ otherwise. When \\y \> 0\\ the \\P\\
replacement term contributes \\H\_{\nu j} \cdot P^2\\ — shared alleles
do **not** give a zero contribution. The entropy \\H\_{\nu j}\\ is the
Shannon entropy of the empirical allele distribution:

\\ H\_{\nu j} = -\sum\_{a \in \mathcal{A}\_j} p\_{ja} \log_2 p\_{ja} \\

### Preprocessing and diagonal

The R reference (`euk_heuristic_distance_r`) and the C++ implementation
preprocess each locus by converting `0` to `NA` and replicating the
single observed allele into both columns for diploid loci with MOI = 1.
Donor matching during imputation compares the **post-preprocessing**
empty-slot patterns, so two isolates that agree on their non-missing
alleles but differ in whether a singleton diploid was replicated are not
considered matches.

Self-pairs (the diagonal) are computed through the full per-locus delta
machinery and are **not** forced to zero — a heterozygous diploid
individual has a non-zero within-sample distance.

### Missing-locus imputation

Barratt’s two-pass mean procedure fills entries where a locus is missing
for an isolate:

1.  For each isolate with a missing locus, find all *complete* isolates
    whose post-preprocessing allele sets and fullness flags match at
    every non-missing locus. Impute the missing-locus delta as the mean
    of those matched isolates’ deltas.
2.  If no match exists, impute using the global mean of the locus.
3.  A second pass catches any rows that remain partially empty after the
    first pass, re-imputing against the first-pass results.

The final distance for pair \\(i, k)\\ is:

\\ D\_{ik}^{(H)} = \sum\_{j=1}^{L} \delta\_{ik}^{(j)} \\

## 2. Bayesian distance (`euk_bayesian_distance`)

Plucinski’s naive-Bayes distance evaluates, for each pair \\(i, k)\\,
three log-likelihoods per locus corresponding to three relatedness
hypotheses. If a pair has **no** comparable locus (each locus is
unobserved for at least one member) the method skips all loci (via
`colSums(..., na.rm = TRUE)`) and returns exactly \\1.0\\ — a uniform
posterior — rather than imputing.

The Bayesian functions do **not** preprocess by duplicating diploid
singleton alleles (unlike the heuristic path), so an isolate with a
single observed allele at a diploid locus contributes one allele instead
of two.

### Per-locus log-likelihoods

Let \\\mathbf{v}\_{ij}\\ be the allele values (not the set, but the
actual observed alleles, with multiplicities) carried by isolate \\i\\
at locus \\j\\, and let \\p\_{ja}\\ be the empirical frequency of allele
\\a\\ as above.

**Hypothesis 0 — unrelated:**

Both isolates draw all their alleles independently from the locus allele
frequencies:

\\ \ell_0 = \sum\_{a \in \mathbf{v}\_{ij}} \log p\_{ja} + \sum\_{a \in
\mathbf{v}\_{kj}} \log p\_{ja} \\

**Hypothesis 1 — share one allele:**

The pair shares exactly one allele; the shared allele is the most likely
one. For every pair of positions \\(u, v)\\ where \\v\_{iju} =
v\_{kjv}\\, compute:

\\ \ell\_{1,uv} = \ell_0 - \log p\_{j,v\_{iju}} \\

and take \\\ell_1 = \max\_{uv} \ell\_{1,uv}\\. If no shared allele
exists, \\\ell_1 = \log(\epsilon \cdot \min(p\_{\text{obs}}))\\, where
\\\epsilon\\ is the estimated allele-capture failure rate (default
0.3072, estimated by Joel Barratt in May 2019) and
\\\min(p\_{\text{obs}})\\ is the smallest frequency among the alleles
carried by the pair.

**Hypothesis 2 — share both alleles:**

Both isolates share a pair of alleles (multiset equality). For every
pair of 2-subsets \\(\mathbf{a}, \mathbf{b})\\ from \\\mathbf{v}\_{ij}\\
and \\\mathbf{v}\_{kj}\\ that are equal as multisets:

\\ \ell\_{2,\mathbf{a}\mathbf{b}} = \sum\_{a \in \mathbf{v}\_{ij}
\setminus \mathbf{a}} \log p\_{ja} + \sum\_{a \in \mathbf{v}\_{kj}} \log
p\_{ja} \\

which is equivalent to \\\ell_0\\ with the log-frequencies of the
**removed** pair elements subtracted (the original R code uses
[`match()`](https://rdrr.io/r/base/match.html) to locate the
first-position index of each pair value, then removes those indices via
negative indexing). Let

\\ n\_{\text{shared}} = \max\_{\mathbf{a},\mathbf{b}}
\bigl\|\operatorname{intersect}(\mathbf{a}, \mathbf{b})\bigr\| \\

be the maximum distinct-allele overlap over all 2-subset combinations.

If a multiset equality is found, \\\ell_2 = \max
\ell\_{2,\mathbf{a}\mathbf{b}}\\. Otherwise:

\\ \ell_2 = \log\bigl\[(\epsilon \cdot \min(p\_{\text{obs}}))^{2 -
n\_{\text{shared}}}\bigr\] \\

For ploidy 1 loci, and when **either** allele vector has length 1, the
package sets \\\ell_2 = \ell_1\\ (a package-level decision — the
original script crashes on these cases).

### Posterior and distance

The per-locus log-likelihoods are summed across loci:

\\ s_0 = \sum\_{j} \ell_0^{(j)}, \quad s_1 = \sum\_{j} \ell_1^{(j)},
\quad s_2 = \sum\_{j} \ell_2^{(j)} \\

These are normalised to a posterior over the three hypotheses:

\\ \Pr(h \mid \text{data}) = \frac{e^{s_h}}{e^{s_0} + e^{s_1} +
e^{s_2}}, \quad h \in \\0, 1, 2\\ \\

The distance is the posterior expected number of shared-allele
hypotheses:

\\ D\_{ik}^{(B)} = 0 \cdot \Pr(0) + 1 \cdot \Pr(1) + 2 \cdot \Pr(2) \in
\[0, 2\] \\

Values near 0 indicate the pair appears unrelated; values near 2 suggest
complete allelic sharing. Self-pair diagonals range approximately
\\1.46\\–\\2.0\\ (a homozygous diploid with common alleles scores near
2; a singleton diploid with a rare allele scores lower).

## 3. Ensemble distance (`euk_ensemble_distance`)

The heuristic and Bayesian distances operate on different scales and
have different sensitivities. The ensemble combines them via one of two
methods chosen with the `method` argument.

### `method = "original"` (default)

Reproduces the published `run.r` workflow verbatim:

1.  **Reverse** the Bayesian direction: \\D_B' = 2 - D_B\\ so that both
    matrices are divergence-oriented (larger = more dissimilar).
2.  **Min-max normalise** each matrix over **all** entries (diagonal,
    upper and lower triangles together) to \\\[0, 1\]\\.
3.  Let \\H\\ be the normalised heuristic and \\B\\ the reversed-
    normalised Bayesian matrix. Sort \\B\\’s values and map them onto
    the rank-order of \\H\\: \\ D\_{ik}^{(E)} =
    B\_{(\text{rank}(H\_{ik}))} \\ where ties in \\H\\ are broken
    **randomly** (`ties.method = "random"`), making the output
    **non-deterministic** unless a `seed` is provided.
4.  Take the elementwise mean of the distribution-mapped matrix and
    \\H\\, then rescale by the global maximum.

Because ties in the heuristic ranking are resolved at random, two calls
with the same data can produce different matrices. The output matrix can
also be **asymmetric** (the random tie-breaking and the full-matrix
(non-triangular) sorting do not guarantee symmetry).

### `method = "rank"`

A package-specific variant that avoids random tie-breaking and preserves
symmetry:

1.  **Reverse** the Bayesian direction: \\D_B' = 2 - D_B\\, so both
    matrices are divergence-oriented.
2.  For an \\n \times n\\ distance matrix \\D\\, let \\U\\ be the set of
    \\\binom{n}{2}\\ strictly upper-triangular entries. Assign
    fractional ranks: \\ r(D)\_{ik} = \frac{\text{rank of } D\_{ik}
    \text{ among } U} {\binom{n}{2} + 1} \\ with ties resolved by
    `ties.method = "average"`. The diagonal is set to `NA`.
3.  The ensemble entry is: \\ D\_{ik}^{(E)} = \frac{r(D_H)\_{ik} +
    r(D_B')\_{ik}}{2} \\

### Interpretation

After reversal both components are divergence-oriented:

- Values near **0** indicate both methods agree on high similarity (low
  divergence).
- Values near **1** indicate agreement on strong divergence.

For the `"rank"` method the output lies in \\(0, 1)\\; for `"original"`
the range depends on the global min-max scaling. The `"rank"` method is
a package variant and is **not** described in the original publication.

## References

Barratt JLN, Park S, Nascimento FS, et al. (2019) “Genotyping
genetically heterogeneous *Cyclospora cayetanensis* infections to
complement epidemiological case linkage”. *Parasitology*
**146**:1275–1283. <doi:10.1017/S0031182019000581>.

Component distances (heuristic and Bayesian) verified against the
original R scripts to within \\10^{-10}\\ on a 68-sample fixture. The
ensemble `"rank"` method is a package-level extension.
