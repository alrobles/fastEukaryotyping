#' Eukaryotyping Bayesian distance --- pure R reference
#'
#' \strong{Deprecated:} use \code{\link{euk_bayesian_distance}} with
#' \code{use_r = FALSE} (the default) for the optimised C++ implementation.
#' This reference implementation is retained for cross-validation and
#' development only.
#'
#' Faithful port of Plucinski's unsupervised naive Bayes distance
#' (`euk_bayesian_fulldataset_V2.r` in `Joel-Barratt/Eukaryotyping`).
#'
#' For each pair of isolates and locus, three log-likelihoods are computed:
#' \itemize{
#'   \item \code{ll0}: samples unrelated --- every allele drawn independently
#'     from the locus allele frequencies.
#'   \item \code{ll1}: samples share exactly one allele --- the shared allele
#'     is the most likely one; remaining alleles drawn independently.
#'   \item \code{ll2} (diploid loci only): samples share both alleles ---
#'     evaluated over all 2-subsets; missing allele draws get probability
#'     \code{epsilon * min(p1, p2)}.
#' }
#' The per-locus values are summed across loci, normalised to a posterior over
#' the three hypotheses, and the distance is the posterior expectation
#' \code{sum(lik * c(0, 1, 2))} in \eqn{[0, 2]}.
#'
#' @param geno List of `n x maxMOI` integer matrices (one per locus). `NA`
#'   marks empty slots; `0` is also treated as missing.
#' @param ploidy Integer vector of length `length(geno)`: 1 = haploid,
#'   2 = diploid.
#' @param epsilon Rate of missed alleles. Default `0.3072`, the value estimated
#'   by Joel Barratt (May 2019) in the original script.
#' @return An `n x n` distance matrix with values in \eqn{[0, 2]}; a pair
#'   sharing no comparable locus returns exactly 1 (uniform posterior), as in
#'   the original script's `colSums(na.rm = TRUE)` convention.
#' @importFrom utils combn
#' @keywords internal
#' @noRd
#' @deprecated Use \code{\link{euk_bayesian_distance}(..., use_r = FALSE)}
#'   instead.
euk_bayesian_distance_r <- function(geno, ploidy, epsilon = 0.3072) {
  message("DeprecationWarning: euk_bayesian_distance_r is deprecated. Use euk_bayesian_distance(..., use_r = FALSE) instead.")
  nids <- nrow(geno[[1L]])
  nloci <- length(geno)

  # --- preprocess: 0 -> NA, locus allele sets + frequencies -----------------
  obs <- vector("list", nloci)
  alleles <- vector("list", nloci)
  frequencies <- vector("list", nloci)
  for (j in seq_len(nloci)) {
    mat <- geno[[j]]
    mat[mat == 0L] <- NA_integer_
    obs[[j]] <- mat
    raw <- as.vector(mat)
    raw <- raw[!is.na(raw)]
    if (length(raw) == 0L) {
      stop(sprintf("Locus %d has no non-missing alleles.", j))
    }
    alleles[[j]] <- sort(unique(raw))
    fr <- vapply(alleles[[j]], function(a) sum(raw == a), numeric(1L))
    frequencies[[j]] <- fr / sum(fr)
  }

  # --- per-locus triple log-likelihood (calculate_loglikelihood2) -----------
  locus_ll <- function(v1, v2, p1, p2, ploid) {
    n1 <- length(v1); n2 <- length(v2)
    ll0 <- sum(log(p1)) + sum(log(p2))

    # ll1: best single shared allele (exact port: match removes FIRST
    # occurrence positions of the shared pair's values)
    best1 <- -Inf
    for (i in seq_len(n1)) for (k in seq_len(n2)) {
      if (v1[i] == v2[k]) {
        cand <- exp(sum(log(p1[-i])) + sum(log(p2)))
        if (cand > best1) best1 <- cand
      }
    }
    ll1 <- if (best1 > 0) log(best1) else log(epsilon * min(c(p1, p2)))

    # ll2: best pair of shared alleles over 2-subsets (diploid only)
    if (ploid <= 1L || n1 < 2L || n2 < 2L) {
      ll2 <- ll1
    } else {
      pairs1 <- combn(v1, 2L, simplify = FALSE)
      pairs2 <- combn(v2, 2L, simplify = FALSE)
      best2 <- -Inf; nshared <- 0L
      for (a in pairs1) for (b in pairs2) {
        # cand test: multiset equality  sum(sort(a) == sort(b)) == 2
        # nshared:   unique overlap      length(intersect(a, b))
        ov <- length(intersect(a, b))
        if (ov > nshared) nshared <- ov
        if (a[1L] == b[1L] && a[2L] == b[2L] ||
            a[1L] == b[2L] && a[2L] == b[1L]) {
          # original: p1[-match(pair, v1)] -- removes first-position matches
          rest1f <- p1[-match(a, v1)]
          cand <- exp(sum(log(rest1f)) + sum(log(p2)))
          if (cand > best2) best2 <- cand
        }
      }
      ll2 <- if (best2 > 0) log(best2) else {
        log((epsilon * min(c(p1, p2)))^(2L - nshared))
      }
    }
    c(ll0, ll1, ll2)
  }

  # --- pairwise distances ---------------------------------------------------
  D <- matrix(NA_real_, nids, nids)
  for (i in seq_len(nids)) {
    for (k in seq_len(nids)) {
      ll <- matrix(NA_real_, nloci, 3L)
      for (j in seq_len(nloci)) {
        v1 <- obs[[j]][i, ]; v1 <- v1[!is.na(v1)]
        v2 <- obs[[j]][k, ]; v2 <- v2[!is.na(v2)]
        if (length(v1) > 0L && length(v2) > 0L) {
          p1 <- frequencies[[j]][match(v1, alleles[[j]])]
          p2 <- frequencies[[j]][match(v2, alleles[[j]])]
          ll[j, ] <- locus_ll(v1, v2, p1, p2, ploidy[j])
        }
      }
      s <- colSums(ll, na.rm = TRUE)
      lik <- exp(s) / sum(exp(s))
      D[i, k] <- sum(lik * c(0, 1, 2))
    }
  }
  D
}


#' Eukaryotyping ensemble distance --- pure R engine
#'
#' Internal engine behind \code{\link{euk_ensemble_distance}}. Two
#' combination rules are provided.
#'
#' `method = "original"` is a verbatim port of the combining procedure in the
#' original `run.r` script: the Bayesian matrix is reversed (\eqn{2 - D_B}),
#' both matrices are min-max normalised over all entries (diagonal and both
#' triangles included), the sorted normalised Bayesian values are mapped onto
#' the ranks of the normalised heuristic matrix with random tie breaking, the
#' ensemble is the elementwise mean of that mapped matrix and the normalised
#' heuristic matrix, and the result is divided by its maximum. Random tie
#' breaking makes the original procedure non-deterministic and can produce a
#' non-symmetric matrix; pass `seed` for reproducibility.
#'
#' `method = "rank"` is a package-specific direction-consistent variant: each
#' matrix's strictly upper-triangular entries are converted to fractional
#' ranks (average ties, denominator `choose(n, 2) + 1`, missing pairs kept
#' `NA`), with the Bayesian matrix first reversed as `2 - D_bayesian` so that
#' both components are divergence-oriented. The ensemble entry is the mean of
#' the two fractional ranks; the diagonal is `NA`. Invariant to strictly
#' increasing transformations that preserve ties and missingness.
#'
#' @param D_heuristic `n x n` heuristic distance matrix.
#' @param D_bayesian `n x n` Bayesian distance matrix (same dimnames).
#' @param method `"original"` (default) or `"rank"`.
#' @param seed Optional integer; `set.seed(seed)` is called before the
#'   original method's random tie breaking. Ignored for `"rank"`.
#' @return `n x n` numeric matrix of ensemble distances.
#' @keywords internal
#' @noRd
euk_ensemble_distance_r <- function(D_heuristic, D_bayesian,
                                  method = c("original", "rank"),
                                  seed = NULL) {
  method <- match.arg(method)

  if (method == "original") {
    # verbatim port of Joel-Barratt/Eukaryotyping run.r lines 15-37
    B_temp <- 2 - D_bayesian
    B_temp <- B_temp - min(B_temp)
    B_norm <- B_temp / max(B_temp)
    H_temp <- D_heuristic - min(D_heuristic, na.rm = TRUE)
    H_norm <- H_temp / max(H_temp, na.rm = TRUE)
    if (!is.null(seed)) set.seed(seed)
    D <- sort(B_norm)[rank(H_norm, ties.method = "random")]
    arr <- array(NA, dim = c(nrow(D_heuristic), ncol(D_heuristic), 2))
    arr[, , 1] <- D
    arr[, , 2] <- H_norm
    out <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))
    out <- out / max(out)
    dimnames(out) <- dimnames(D_heuristic)
    return(out)
  }

  rank_upper <- function(D) {
    iu <- which(upper.tri(D), arr.ind = TRUE)
    r <- rank(D[iu], ties.method = "average", na.last = "keep")
    out <- matrix(NA_real_, nrow(D), ncol(D))
    out[iu] <- r
    out[cbind(iu[, 2L], iu[, 1L])] <- r
    out / (length(r) + 1)
  }
  # 2 - D_bayesian: reverse the similarity direction so both components
  # rank divergence in the same direction
  (rank_upper(D_heuristic) + rank_upper(2 - D_bayesian)) / 2
}
