# How the Eukaryotyping distance ensemble works

This vignette explains the three distance methods that make up the
Eukaryotyping ensemble. The methods were first published by Barratt et
al. (2019)
([doi:10.1017/S0031182019000581](https://doi.org/10.1017/S0031182019000581))
for genotyping *Cyclospora cayetanensis* from MLST data, and are
reimplemented here in C++17 with RcppParallel.

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
sum of per-locus deltas.

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

For **haploid** loci the expression simplifies to:

\\ \delta = H\_{\nu j} \cdot 2x \cdot \mathbb{I}(y = 0) \cdot k \\

with the same \\k\\ convention. The entropy \\H\_{\nu j}\\ is the
Shannon entropy of the empirical allele distribution:

\\ H\_{\nu j} = -\sum\_{a \in \mathcal{A}\_j} p\_{ja} \log_2 p\_{ja} \\

### Missing-locus imputation

A pair of isolates that share **no comparable locus** (both missing)
cannot contribute to the sum. Barratt’s two-pass mean procedure fills
these entries:

1.  For each isolate with a missing locus, find all *complete* isolates
    whose allele sets match at every non-missing locus. Impute the
    missing-locus delta as the mean of those matched isolates’ deltas.
2.  If no match exists, impute using the global mean of the locus.
3.  A second pass catches any rows that remain partially empty after the
    first pass, re-imputing against the first-pass results.

The final distance for pair \\(i, k)\\ is:

\\ D\_{ik}^{(H)} = \sum\_{j=1}^{L} \delta\_{ik}^{(j)} \\

## 2. Bayesian distance (`euk_bayesian_distance`)

Plucinski’s naive-Bayes distance evaluates, for each pair \\(i, k)\\,
three log-likelihoods per locus corresponding to three relatedness
hypotheses.

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

**Hypothesis 2 — share both alleles** (diploid loci only):

Both isolates share a pair of alleles (multiset equality). For every
pair of 2-subsets \\(\mathbf{a}, \mathbf{b})\\ of \\\mathbf{v}\_{ij}\\
and \\\mathbf{v}\_{kj}\\ that are equal as multisets:

\\ \ell\_{2,\mathbf{a}\mathbf{b}} = \bigl(\ell_0 - \log p\_{j,a_1} -
\log p\_{j,a_2}\bigr) + \sum\_{a \in \mathbf{v}\_{kj}} \log p\_{ja} \\

and \\\ell_2 = \max \ell\_{2,\mathbf{a}\mathbf{b}}\\. If no multiset
equality exists:

\\ \ell_2 = \log\bigl\[(\epsilon \cdot \min(p\_{\text{obs}}))^{2 -
n\_{\text{shared}}}\bigr\] \\

where \\n\_{\text{shared}}\\ is the number of distinct allele values
shared by the two subsets.

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
complete allelic sharing.

## 3. Ensemble distance (`euk_ensemble_distance`)

The heuristic and Bayesian distances operate on different scales and
have different sensitivities. The ensemble combines them via
**fractional rank averaging**, a copula-style summary that is invariant
to monotone transformations of either distance.

For an \\n \times n\\ distance matrix \\D\\, let \\U\\ be the set of
\\\binom{n}{2}\\ unordered upper-triangular entries. Assign fractional
ranks:

\\ r(D)\_{ik} = \frac{\text{rank of } D\_{ik} \text{ among } U} {\|U\| +
1} \\

so that \\r(D)\_{ik} \in (0, 1)\\. The ensemble entry is:

\\ D\_{ik}^{(E)} = \frac{r\bigl(D^{(H)}\bigr)\_{ik} +
r\bigl(D^{(B)}\bigr)\_{ik}}{2} \\

Values near 0 indicate both methods agree on high similarity; values
near 1 indicate strong divergence.

## References

Barratt JLN, Plucinski M, Nascimento F, et al. (2019) “A Genotyping
Approach to Determine the Geographic Origin and Relatedness of
Cyclospora cayetanensis” *Parasitology*,
<doi:10.1017/S0031182019000581>.
