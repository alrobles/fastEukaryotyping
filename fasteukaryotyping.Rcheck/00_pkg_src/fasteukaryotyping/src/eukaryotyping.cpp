// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <cstdint>

using namespace Rcpp;
using namespace RcppParallel;

// population count for uint64_t
static inline int popcount64(uint64_t x) {
#ifdef __GNUC__
  return __builtin_popcountll(x);
#else
  int c = 0;
  while (x) { c++; x &= x - 1; }
  return c;
#endif
}

static inline int single_bit_position(uint64_t x) {
  // x is assumed to have exactly one bit set.
#ifdef __GNUC__
  return __builtin_ctzll(x);
#else
  for (int i = 0; i < 64; ++i) {
    if ((x >> i) & 1ULL) return i;
  }
  return 0;
#endif
}

static inline double log2_d(double x) {
  return std::log(x) / std::log(2.0);
}

struct LocusData {
  int ploidy;
  int m;
  int max_moi;
  double H_nu;
  int nonmissing_total;
  std::vector<uint64_t> row_mask;       // distinct allele bitmasks per row
  std::vector<int> row_count;          // number of non-missing slots per row
  std::vector<uint8_t> is_full;        // effective count == max_moi (diploid MOI==1 replicated)
  std::vector<uint64_t> distinct_masks;
  std::vector<int> distinct_mask_counts;
};

static double compute_H_nu(const std::vector<int>& allele_counts, int total_slots) {
  if (total_slots == 0) return 0.0;
  double H = 0.0;
  for (size_t a = 0; a < allele_counts.size(); ++a) {
    if (allele_counts[a] > 0) {
      double p = static_cast<double>(allele_counts[a]) / total_slots;
      H -= p * log2_d(p);
    }
  }
  return H;
}

static int count_rows_with_mask(const std::vector<uint64_t>& distinct_masks,
                                const std::vector<int>& distinct_mask_counts,
                                uint64_t subset_mask) {
  int total = 0;
  for (size_t k = 0; k < distinct_masks.size(); ++k) {
    if ((distinct_masks[k] & subset_mask) == subset_mask) {
      total += distinct_mask_counts[k];
    }
  }
  return total;
}

static double sub_per_locus(const LocusData& loc, int i, int j) {
  int count_i = loc.row_count[i];
  int count_j = loc.row_count[j];
  if (count_i == 0 || count_j == 0) {
    return NA_REAL;
  }

  uint64_t v1_mask = loc.row_mask[i];
  uint64_t v2_mask = loc.row_mask[j];

  int pop1 = popcount64(v1_mask);
  int pop2 = popcount64(v2_mask);
  int x = pop1 + pop2;

  uint64_t shared_mask = v1_mask & v2_mask;
  int y = popcount64(shared_mask);

  double P = 1.0;
  if (y > 0) {
    int numerator = count_rows_with_mask(loc.distinct_masks, loc.distinct_mask_counts, shared_mask);
    int denom = loc.nonmissing_total;
    if (denom > 0) {
      double ratio = static_cast<double>(numerator) / denom;
      P = ratio * ratio;
    }
  }

  double delta;
  if (loc.ploidy > 1) {
    int n = std::min(pop1, pop2);
    int w = 0;
    if (n > 1) {
      w = x;
    } else if (n == 1 && x == 2) {
      w = 4;
    } else if (n == 1 && x > 2) {
      w = 1 + x;
    }

    int m = loc.m;
    int jj = (m == 1) ? 2 : m;

    int z = 0;
    int term1_cond = (n == 1 && y == 1 && x > 2) ? 1 : 0;
    z += 3 * term1_cond;

    int a = (n > 1 && jj >= y) ? y : 0;
    int b = (n > 1 && y > jj) ? jj : 0;
    int c = (n == 1 && x == 2 && y == 1) ? jj : 0;
    z += 2 * (a + b + c);

    int delta_nu_raw = 0;
    if (y == 0) {
      delta_nu_raw = w;
    } else {
      delta_nu_raw = 2 * jj;
      // subtract z if it lies in [jj, 2*jj]
      if (z >= jj && z <= 2 * jj) {
        delta_nu_raw -= z;
      }
    }

    double kk = (y == 0) ? 1.0 : P;
    double raw = (delta_nu_raw > 0) ? static_cast<double>(delta_nu_raw)
                                     : P;
    delta = loc.H_nu * raw * kk;
  } else {
    int delta_ex_raw = 2 * x * (y == 0 ? 1 : 0);
    double kk = (y == 0) ? 1.0 : P;
    double raw = (delta_ex_raw > 0) ? static_cast<double>(delta_ex_raw)
                                    : P;
    delta = loc.H_nu * raw * kk;
  }

  return delta;
}

struct HeuristicWorker : public Worker {
  const std::vector<LocusData>& loci;
  const int n;
  const int L;
  const size_t per_loc;
  std::vector<double>& deltas;

  HeuristicWorker(const std::vector<LocusData>& loci_, int n_, int L_, std::vector<double>& deltas_)
    : loci(loci_), n(n_), L(L_), per_loc(static_cast<size_t>(n_) * n_), deltas(deltas_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t pp = begin; pp < end; ++pp) {
      int64_t p = static_cast<int64_t>(pp);
      // invert p = i*(i+1)/2 + j, 0 <= j <= i < n
      int i = static_cast<int>((std::sqrt(1.0 + 8.0 * p) - 1.0) / 2.0);
      if (i >= n) i = n - 1;
      int64_t tri = static_cast<int64_t>(i) * (i + 1) / 2;
      while (tri > p) {
        --i;
        tri = static_cast<int64_t>(i) * (i + 1) / 2;
      }
      while (tri + i < p) {
        ++i;
        tri = static_cast<int64_t>(i) * (i + 1) / 2;
      }
      int j = static_cast<int>(p - tri);

      for (int l = 0; l < L; ++l) {
        double d = sub_per_locus(loci[l], i, j);
        size_t base = static_cast<size_t>(l) * per_loc;
        deltas[base + static_cast<size_t>(i) * n + j] = d;
        if (i != j) {
          deltas[base + static_cast<size_t>(j) * n + i] = d;
        }
      }
    }
  }
};

static void compute_matrix_from_deltas(const std::vector<double>& deltas,
                                       NumericMatrix& D,
                                       int n, int L) {
  size_t per_loc = static_cast<size_t>(n) * n;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double total = 0.0;
      bool has_na = false;
      size_t idx_base = static_cast<size_t>(i) * n + j;
      for (int l = 0; l < L; ++l) {
        double d = deltas[static_cast<size_t>(l) * per_loc + idx_base];
        if (std::isnan(d) || ISNA(d)) {
          has_na = true;
          break;
        }
        total += d;
      }
      D(i, j) = has_na ? NA_REAL : total;
    }
  }
}

static void impute_pass(const std::vector<LocusData>& loci,
                        const std::vector<double>& source,
                        std::vector<double>& target,
                        const std::vector<int>& whichna,
                        const std::vector<uint8_t>& is_whichna,
                        int n, int L) {
  size_t per_loc = static_cast<size_t>(n) * n;
  for (int isolate : whichna) {
    std::vector<int> missing_loci;
    std::vector<int> nonmissing_loci;
    missing_loci.reserve(L);
    nonmissing_loci.reserve(L);
    for (int l = 0; l < L; ++l) {
      if (loci[l].row_count[isolate] == 0) {
        missing_loci.push_back(l);
      } else {
        nonmissing_loci.push_back(l);
      }
    }

    std::vector<int> match;
    for (int s = 0; s < n; ++s) {
      if (is_whichna[s]) continue;
      bool ok = true;
      for (int l : nonmissing_loci) {
        if (loci[l].row_mask[isolate] != loci[l].row_mask[s] ||
            loci[l].is_full[isolate] != loci[l].is_full[s]) {
          ok = false;
          break;
        }
      }
      if (ok) match.push_back(s);
    }

    for (int l : missing_loci) {
      size_t base = static_cast<size_t>(l) * per_loc;
      if (!match.empty()) {
        for (int x = 0; x < n; ++x) {
          double sum = 0.0;
          int cnt = 0;
          size_t row_off = static_cast<size_t>(x) * n;
          for (int s : match) {
            double d = source[base + row_off + s];
            if (!std::isnan(d) && !ISNA(d)) {
              sum += d;
              ++cnt;
            }
          }
          double val = (cnt > 0) ? (sum / cnt) : NA_REAL;
          target[base + static_cast<size_t>(isolate) * n + x] = val;
          target[base + static_cast<size_t>(x) * n + isolate] = val;
        }
        double sumd = 0.0;
        int cntd = 0;
        for (int s : match) {
          double d = source[base + static_cast<size_t>(s) * n + s];
          if (!std::isnan(d) && !ISNA(d)) {
            sumd += d;
            ++cntd;
          }
        }
        target[base + static_cast<size_t>(isolate) * n + isolate] =
          (cntd > 0) ? (sumd / cntd) : NA_REAL;
      } else {
        double sum = 0.0;
        int cnt = 0;
        for (int a = 0; a < n; ++a) {
          for (int b = 0; b < n; ++b) {
            double d = source[base + static_cast<size_t>(a) * n + b];
            if (!std::isnan(d) && !ISNA(d)) {
              sum += d;
              ++cnt;
            }
          }
        }
        double val = (cnt > 0) ? (sum / cnt) : NA_REAL;
        for (int x = 0; x < n; ++x) {
          target[base + static_cast<size_t>(isolate) * n + x] = val;
          target[base + static_cast<size_t>(x) * n + isolate] = val;
        }
        double sumd = 0.0;
        int cntd = 0;
        for (int s = 0; s < n; ++s) {
          double d = target[base + static_cast<size_t>(s) * n + s];
          if (!std::isnan(d) && !ISNA(d)) {
            sumd += d;
            ++cntd;
          }
        }
        target[base + static_cast<size_t>(isolate) * n + isolate] =
          (cntd > 0) ? (sumd / cntd) : val;
      }
    }
  }
}

//' Internal heuristic distance (C++ implementation)
//'
//' RcppParallel reimplementation of \\code{euk_heuristic_distance_r()}.
//' Each locus may have at most 64 distinct alleles.
//'
//' @param geno List of integer matrices.
//' @param ploidy Integer vector of ploidy values.
//' @param m Integer vector of m values.
//' @param impute_data Integer (0/1). If 1, apply the two-pass Barratt
//'   missing-locus imputation.
//' @return A symmetric distance matrix.
// [[Rcpp::export(.euk_heuristic_distance_impl)]]
NumericMatrix euk_heuristic_distance_impl(List geno, IntegerVector ploidy, IntegerVector m, int impute_data) {
  if (geno.size() == 0) {
    stop("'geno' must contain at least one locus.");
  }
  if (geno.size() != ploidy.size()) {
    stop("'ploidy' must have length equal to 'geno'.");
  }
  if (geno.size() != m.size()) {
    stop("'m' must have length equal to 'geno'.");
  }

  IntegerMatrix first = geno[0];
  int n = first.nrow();
  int L = geno.size();
  size_t per_loc = static_cast<size_t>(n) * n;

  std::vector<LocusData> loci;
  loci.reserve(L);

  for (int l = 0; l < L; ++l) {
    IntegerMatrix mat = geno[l];
    if (mat.nrow() != n) {
      stop("All loci must have the same number of rows.");
    }

    LocusData loc;
    loc.ploidy = ploidy[l];
    loc.m = m[l];
    loc.max_moi = mat.ncol();
    loc.row_mask.assign(n, 0ULL);
    loc.row_count.assign(n, 0);
    loc.is_full.assign(n, 0);

    std::map<int, int> allele_to_bit;
    std::vector<int> allele_counts;
    int total_slots = 0;

    for (int r = 0; r < n; ++r) {
      uint64_t mask = 0ULL;
      int count = 0;
      for (int c = 0; c < loc.max_moi; ++c) {
        int val = mat(r, c);
        if (IntegerVector::is_na(val) || val == 0) continue;
        ++count;
        auto it = allele_to_bit.find(val);
        int bit;
        if (it == allele_to_bit.end()) {
          bit = allele_to_bit.size();
          if (bit >= 64) {
            stop("Locus %d has more than 64 distinct alleles; bitset fallback not yet implemented.", l + 1);
          }
          allele_to_bit[val] = bit;
          allele_counts.push_back(0);
        } else {
          bit = it->second;
        }
        mask |= (1ULL << bit);
        allele_counts[bit]++;
        total_slots++;
      }
      // Replicate diploid single-allele rows as import_data_V2.r does:
      // the single observed allele occupies two slots, so the frequency count
      // must be incremented once more.
      if (loc.ploidy > 1 && count == 1 && loc.max_moi >= 2) {
        int bit = single_bit_position(mask);
        allele_counts[bit]++;
        total_slots++;
        count++;
      }

      loc.row_mask[r] = mask;
      loc.row_count[r] = count;
      loc.is_full[r] = (count == loc.max_moi) ? 1 : 0;
    }

    loc.H_nu = compute_H_nu(allele_counts, total_slots);

    loc.nonmissing_total = 0;
    std::map<uint64_t, int> mask_counts;
    for (int r = 0; r < n; ++r) {
      if (loc.row_count[r] > 0) {
        ++loc.nonmissing_total;
        mask_counts[loc.row_mask[r]]++;
      }
    }
    for (const auto& kv : mask_counts) {
      loc.distinct_masks.push_back(kv.first);
      loc.distinct_mask_counts.push_back(kv.second);
    }

    loci.push_back(std::move(loc));
  }

  // Compute all per-locus pairwise delta matrices
  std::vector<double> deltas(static_cast<size_t>(L) * per_loc, NA_REAL);
  int64_t npairs = static_cast<int64_t>(n) * (n + 1) / 2;
  if (npairs > 0) {
    HeuristicWorker worker(loci, n, L, deltas);
    parallelFor(0, static_cast<size_t>(npairs), worker);
  }

  if (!impute_data) {
    NumericMatrix D(n, n);
    compute_matrix_from_deltas(deltas, D, n, L);
    return D;
  }

  // Samples with at least one missing locus (rows that are all NA before imputation)
  std::vector<int> whichna;
  std::vector<uint8_t> is_whichna(n, 0);
  for (int i = 0; i < n; ++i) {
    bool has_missing = false;
    for (int l = 0; l < L; ++l) {
      if (loci[l].row_count[i] == 0) {
        has_missing = true;
        break;
      }
    }
    if (has_missing) {
      whichna.push_back(i);
      is_whichna[i] = 1;
    }
  }

  // First imputation pass
  std::vector<double> deltas_imp = deltas;
  impute_pass(loci, deltas, deltas_imp, whichna, is_whichna, n, L);

  // Second pass: identify rows that still contain NA after first pass
  std::vector<int> whichna2;
  for (int i = 0; i < n; ++i) {
    bool any_na = false;
    for (int j = 0; j < n && !any_na; ++j) {
      for (int l = 0; l < L; ++l) {
        double d = deltas_imp[static_cast<size_t>(l) * per_loc + static_cast<size_t>(i) * n + j];
        if (std::isnan(d) || ISNA(d)) {
          any_na = true;
          break;
        }
      }
    }
    if (any_na) whichna2.push_back(i);
  }

  // Copy first-pass results into the original array and run second pass there
  std::copy(deltas_imp.begin(), deltas_imp.end(), deltas.begin());
  impute_pass(loci, deltas_imp, deltas, whichna2, is_whichna, n, L);

  NumericMatrix D(n, n);
  compute_matrix_from_deltas(deltas, D, n, L);
  return D;
}


// ---------------------------------------------------------------------------
// Bayesian (Plucinski unsupervised naive Bayes) distance
// ---------------------------------------------------------------------------

struct BayesLocus {
  int ploidy;
  int max_moi;
  std::vector<std::vector<int>> row_alleles;   // non-missing allele values per row
  std::vector<std::vector<int>> row_freqidx;   // index into allele log-freq per allele
  std::vector<double> log_freq;                // log allele frequencies
  double min_freq;
};

// Per-locus triple log-likelihood: port of calculate_loglikelihood2.
static inline void bayes_locus_ll(const BayesLocus& L,
                                  const std::vector<int>& v1,
                                  const std::vector<int>& f1,
                                  const std::vector<int>& v2,
                                  const std::vector<int>& f2,
                                  double epsilon, double out[3]) {
  int n1 = (int)v1.size(), n2 = (int)v2.size();
  double ll0 = 0.0;
  for (int i = 0; i < n1; ++i) ll0 += L.log_freq[f1[i]];
  for (int k = 0; k < n2; ++k) ll0 += L.log_freq[f2[k]];
  double sump2 = 0.0;
  for (int k = 0; k < n2; ++k) sump2 += L.log_freq[f2[k]];
  std::vector<double> lp1(n1);
  for (int i = 0; i < n1; ++i) lp1[i] = L.log_freq[f1[i]];
  // min(c(p1,p2)): minimum frequency over the alleles carried by the pair
  double minobs = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n1; ++i) minobs = std::min(minobs, std::exp(L.log_freq[f1[i]]));
  for (int k = 0; k < n2; ++k) minobs = std::min(minobs, std::exp(L.log_freq[f2[k]]));

  // ll1: best single shared allele
  double best1 = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n1; ++i) {
    for (int k = 0; k < n2; ++k) {
      if (v1[i] == v2[k]) {
        // exp(sum log p1[-i] + sum log p2) = exp(ll0 - lp1[i])
        double cand = std::exp(ll0 - lp1[i]);
        if (cand > best1) best1 = cand;
      }
    }
  }
  double ll1 = std::log(best1);
  if (!std::isfinite(ll1)) ll1 = std::log(epsilon * minobs);

  // ll2: best shared allele pair over 2-subsets (diploid only)
  double ll2 = ll1;
  if (L.ploidy > 1 && n1 >= 2 && n2 >= 2) {
    double best2 = -std::numeric_limits<double>::infinity();
    int nshared = 0;
    for (int a = 0; a < n1; ++a) {
      for (int b = a + 1; b < n1; ++b) {
        for (int c = 0; c < n2; ++c) {
          for (int d = c + 1; d < n2; ++d) {
            // nshared: |unique(pair1) intersect unique(pair2)| — R's
            // length(intersect(a, b)); duplicated values count once.
            int ov = (v1[a] == v2[c] || v1[a] == v2[d]) +
                     (v1[b] != v1[a] && (v1[b] == v2[c] || v1[b] == v2[d]));
            if (ov > nshared) nshared = ov;
            // cand test: multiset equality, R's
            // sum(sort(pair1) == sort(pair2)) == 2
            bool matched2 =
                (v1[a] == v2[c] && v1[b] == v2[d]) ||
                (v1[a] == v2[d] && v1[b] == v2[c]);
            if (matched2) {
              // original: p1[-match(pair, v1)] -- removes the FIRST-position
              // index of each pair VALUE in v1 (duplicates collapse)
              int ra = 0;
              while (v1[ra] != v1[a]) ++ra;
              int rb = 0;
              while (v1[rb] != v1[b]) ++rb;
              double s = 0.0;
              for (int i = 0; i < n1; ++i)
                if (i != ra && i != rb) s += lp1[i];
              double cand = std::exp(s + sump2);
              if (cand > best2) best2 = cand;
            }
          }
        }
      }
    }
    ll2 = std::log(best2);
    if (!std::isfinite(ll2)) {
      ll2 = std::log(std::pow(epsilon * minobs, 2.0 - nshared));
    }
  }
  out[0] = ll0; out[1] = ll1; out[2] = ll2;
}

struct BayesWorker : public Worker {
  const std::vector<BayesLocus>& loci;
  const double epsilon;
  const int n;
  RMatrix<double> out;

  BayesWorker(const std::vector<BayesLocus>& loci, double epsilon,
              int n, NumericMatrix out)
    : loci(loci), epsilon(epsilon), n(n), out(out) {}

  void operator()(std::size_t begin, std::size_t end) {
    int L = (int)loci.size();
    for (std::size_t idx = begin; idx < end; ++idx) {
      int i = (int)(idx / n);
      int k = (int)(idx % n);
      if (k <= i) continue;
      double s0 = 0.0, s1 = 0.0, s2 = 0.0;
      bool any_locus = false;
      for (int j = 0; j < L; ++j) {
        const std::vector<int>& v1 = loci[j].row_alleles[i];
        const std::vector<int>& v2 = loci[j].row_alleles[k];
        if (v1.empty() || v2.empty()) continue;
        double ll[3];
        bayes_locus_ll(loci[j], v1, loci[j].row_freqidx[i],
                       v2, loci[j].row_freqidx[k], epsilon, ll);
        s0 += ll[0]; s1 += ll[1]; s2 += ll[2];
        any_locus = true;
      }
      if (!any_locus) { out(i, k) = NA_REAL; out(k, i) = NA_REAL; continue; }
      double e0 = std::exp(s0), e1 = std::exp(s1), e2 = std::exp(s2);
      double den = e0 + e1 + e2;
      out(i, k) = out(k, i) = (e1 + 2.0 * e2) / den;
    }
  }
};

//' Internal Bayesian distance (C++ implementation)
//'
//' RcppParallel reimplementation of \code{euk_bayesian_distance_r()}
//' (Plucinski's unsupervised naive Bayes classifier).
//'
//' @param geno List of integer matrices.
//' @param ploidy Integer vector of ploidy values.
//' @param epsilon Rate of missed alleles.
//' @return A symmetric distance matrix with values in \eqn{[0, 2]}.
// [[Rcpp::export(.euk_bayesian_distance_impl)]]
NumericMatrix euk_bayesian_distance_impl(List geno, IntegerVector ploidy, double epsilon) {
  if (geno.size() == 0) {
    stop("'geno' must contain at least one locus.");
  }
  if (geno.size() != ploidy.size()) {
    stop("'ploidy' must have length equal to 'geno'.");
  }
  IntegerMatrix first = geno[0];
  int n = first.nrow();
  int L = geno.size();

  std::vector<BayesLocus> loci(L);
  for (int j = 0; j < L; ++j) {
    IntegerMatrix mat = geno[j];
    BayesLocus& B = loci[j];
    B.ploidy = ploidy[j];
    B.max_moi = mat.ncol();
    B.row_alleles.assign(n, std::vector<int>());
    B.row_freqidx.assign(n, std::vector<int>());
    // allele frequencies over all slots
    std::map<int, int> counts;
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < B.max_moi; ++c) {
        int a = mat(i, c);
        if (a != NA_INTEGER && a != 0) counts[a]++;
      }
    }
    if (counts.empty()) stop("Locus has no non-missing alleles.");
    std::vector<int> als;
    double tot = 0.0;
    for (auto& kv : counts) { als.push_back(kv.first); tot += kv.second; }
    B.min_freq = 1.0;
    for (size_t a = 0; a < als.size(); ++a) {
      double f = counts[als[a]] / tot;
      B.log_freq.push_back(std::log(f));
      if (f < B.min_freq) B.min_freq = f;
    }
    for (int i = 0; i < n; ++i) {
      for (int c = 0; c < B.max_moi; ++c) {
        int a = mat(i, c);
        if (a == NA_INTEGER || a == 0) continue;
        auto it = std::find(als.begin(), als.end(), a);
        B.row_alleles[i].push_back(a);
        B.row_freqidx[i].push_back((int)(it - als.begin()));
      }
    }
  }

  NumericMatrix D(n, n);
  for (int i = 0; i < n; ++i) {
    // diagonal through the same code path (i == k): compute like a pair
    double s0 = 0.0, s1 = 0.0, s2 = 0.0; bool any_locus = false;
    for (int j = 0; j < L; ++j) {
      const std::vector<int>& v1 = loci[j].row_alleles[i];
      if (v1.empty()) continue;
      double ll[3];
      bayes_locus_ll(loci[j], v1, loci[j].row_freqidx[i],
                     v1, loci[j].row_freqidx[i], epsilon, ll);
      s0 += ll[0]; s1 += ll[1]; s2 += ll[2];
      any_locus = true;
    }
    if (!any_locus) { D(i, i) = NA_REAL; }
    else {
      double e0 = std::exp(s0), e1 = std::exp(s1), e2 = std::exp(s2);
      D(i, i) = (e1 + 2.0 * e2) / (e0 + e1 + e2);
    }
  }
  BayesWorker w(loci, epsilon, n, D);
  parallelFor(0, static_cast<std::size_t>(n) * n, w);
  return D;
}
