#' Eukaryotyping heuristic distance --- pure R reference
#'
#' \strong{Deprecated:} use \code{\link{euk_heuristic_distance}} with
#' \code{use_r = FALSE} (the default) for the optimised C++ implementation.
#' This reference implementation is retained for cross-validation and
#' development only.
#'
#' Internal reference implementation that mirrors the original
#' `euk_heuristic_fulldataset.r` script line-by-line, including the
#' two-pass missing-locus imputation.
#'
#' @param geno List of `n x maxMOI` integer matrices (one per locus). `NA`
#'   marks empty slots; `0` is also treated as missing.
#' @param ploidy Integer vector of length `length(geno)`: 1 = haploid,
#'   2 = diploid.
#' @param m Integer vector of length `length(geno)`. Default `rep(1L, nloci)`,
#'   matching the original script.
#' @param impute_data Logical. If `TRUE`, missing loci are imputed with the
#'   Barratt two-step mean procedure.
#' @return An `n x n` distance matrix. Symmetric. The diagonal is the
#'   within-sample distance and is not forced to zero.
#' @importFrom utils combn
#' @keywords internal
#' @noRd
#' @deprecated Use \code{\link{euk_heuristic_distance}(..., use_r = FALSE)}
#'   instead.
euk_heuristic_distance_r <- function(geno, ploidy,
                                     m = rep(1L, length(geno)),
                                     impute_data = FALSE) {
  message("DeprecationWarning: euk_heuristic_distance_r is deprecated. Use euk_heuristic_distance(..., use_r = FALSE) instead.")
  nids <- nrow(geno[[1L]])
  nloci <- length(geno)

  # --- Pre-process geno: 0 -> NA and replicate diploid single-allele rows -----
  observeddatamatrix <- vector("list", nloci)
  for (j in seq_len(nloci)) {
    mat <- geno[[j]]
    if (is.null(mat) || nrow(mat) != nids) {
      stop("All elements of 'geno' must be matrices with the same number of rows.")
    }
    mat[mat == 0L] <- NA_integer_

    # Original `import_data_V2.r` duplicates the single observed allele into
    # the first two columns for diploid loci with MOI == 1.
    if (ploidy[j] > 1L && ncol(mat) >= 2L) {
      for (i in seq_len(nids)) {
        nonna <- which(!is.na(mat[i, ]))
        if (length(nonna) == 1L) {
          mat[i, 1L:2L] <- mat[i, nonna]
        }
      }
    }
    observeddatamatrix[[j]] <- mat
  }

  # --- Per-locus allele frequencies and entropy -----------------------------
  alleles <- vector("list", nloci)
  frequencies <- vector("list", nloci)
  H_nu <- numeric(nloci)
  nonmissing_total <- integer(nloci)
  row_sets <- vector("list", nloci)
  is_full <- vector("list", nloci)

  for (j in seq_len(nloci)) {
    mat <- observeddatamatrix[[j]]
    max_moi <- ncol(mat)
    raw_alleles <- as.vector(mat)
    raw_alleles <- raw_alleles[!is.na(raw_alleles)]
    if (length(raw_alleles) == 0L) {
      stop(sprintf("Locus %d has no non-missing alleles.", j))
    }
    alleles[[j]] <- sort(unique(raw_alleles))
    freq <- vapply(alleles[[j]], function(a) sum(raw_alleles == a), integer(1L))
    freq <- freq / sum(freq)
    frequencies[[j]] <- freq
    H_nu[j] <- -sum(freq * logb(freq, 2L))

    nonmissing_total[j] <- sum(apply(mat, 1L, function(r) any(!is.na(r))))

    # Store unique allele sets and "fullness" (no NA slot) for fast set-equal
    sets_j <- vector("list", nids)
    full_j <- logical(nids)
    for (i in seq_len(nids)) {
      vals <- mat[i, ]
      vals <- vals[!is.na(vals)]
      sets_j[[i]] <- unique(vals)
      # Effective count after diploid MOI == 1 replication
      eff_count <- length(vals) + (ploidy[j] > 1L && length(vals) == 1L)
      full_j[i] <- eff_count == max_moi
    }
    row_sets[[j]] <- sets_j
    is_full[[j]] <- full_j
  }

  # --- Per-locus delta function ----------------------------------------------
  sub_per_locus <- function(isolate1, isolate2, j) {
    s1 <- row_sets[[j]][[isolate1]]
    s2 <- row_sets[[j]][[isolate2]]
    if (length(s1) == 0L || length(s2) == 0L) {
      return(NA_real_)
    }

    x <- length(s1) + length(s2)
    y <- length(intersect(s1, s2))

    shared_alleles <- intersect(s1, s2)
    if (length(shared_alleles) > 0L) {
      numerator <- sum(vapply(seq_len(nids), function(r) {
        all(shared_alleles %in% row_sets[[j]][[r]])
      }, logical(1L)))
      P <- (numerator / nonmissing_total[j])^2
    } else {
      P <- 1.0
    }

    if (ploidy[j] > 1L) {
      n <- min(length(s1), length(s2))
      w <- x * (n > 1L) +
           4L * (n == 1L) * (x == 2L) +
           (1L + x) * (n == 1L) * (x > 2L)
      jj <- 2L * (m[j] == 1L) + m[j] * (m[j] > 1L)
      z <- 3L * (((2L * (n == 1L) + 1L * (y == 1L) * (x > 2L))) == 3L) +
           2L * (y * (n > 1L) * (jj >= y) +
                 jj * (n > 1L) * (y > jj) +
                 jj * (n == 1L) * (x == 2L) * (y == 1L))
      delta_nu_raw <- w * (y == 0L) +
                      2L * jj * (y > 0L) +
                      sum(sapply(jj:(2L * jj), function(ii) -ii * (z == ii)))
      k <- 1L * (y == 0L) + P * (y > 0L)
      delta <- H_nu[j] * (delta_nu_raw * (delta_nu_raw > 0L) +
                           P * (delta_nu_raw == 0L)) * k
    } else {
      delta_ex_raw <- 2L * x * (y == 0L)
      k <- 1L * (y == 0L) + P * (y > 0L)
      delta <- H_nu[j] * (delta_ex_raw * (delta_ex_raw > 0L) +
                           P * (delta_ex_raw == 0L)) * k
    }
    delta
  }

  # --- Compute per-locus distance matrices (all pairs, including diagonal) ---
  D_list <- lapply(seq_len(nloci), function(j) matrix(NA_real_, nids, nids))
  for (j in seq_len(nloci)) {
    for (i in seq_len(nids)) {
      for (k in seq_len(i)) {
        d <- sub_per_locus(i, k, j)
        D_list[[j]][i, k] <- d
        D_list[[j]][k, i] <- d
      }
    }
  }

  D_total <- Reduce("+", D_list)

  if (!impute_data) {
    return(D_total)
  }

  # --- Missing-locus imputation (two passes, matching original script) -------

  # Which samples have at least one completely missing locus?
  has_missing_locus <- vapply(seq_len(nids), function(i) {
    any(vapply(seq_len(nloci), function(j) length(row_sets[[j]][[i]]) == 0L, logical(1L)))
  }, logical(1L))
  whichna <- which(has_missing_locus)

  set_equal_rows <- function(isolate, j, candidates) {
    target_set <- row_sets[[j]][[isolate]]
    target_full <- is_full[[j]][[isolate]]
    vapply(candidates, function(s) {
      setequal(target_set, row_sets[[j]][[s]]) &&
        (target_full == is_full[[j]][[s]])
    }, logical(1L))
  }

  impute_pass <- function(source, target, isolate_indices) {
    if (length(isolate_indices) == 0L) {
      return(target)
    }
    for (i in isolate_indices) {
      missing_loci <- which(vapply(seq_len(nloci), function(j) {
        length(row_sets[[j]][[i]]) == 0L
      }, logical(1L)))
      nonmissing_loci <- setdiff(seq_len(nloci), missing_loci)

      if (length(nonmissing_loci) > 0L) {
        match_list <- lapply(nonmissing_loci, function(j) {
          set_equal_rows(i, j, seq_len(nids))
        })
        matches <- Reduce("&", match_list)
        matches[whichna] <- FALSE
        match_samples <- which(matches)
      } else {
        match_samples <- integer(0L)
      }

      for (j in missing_loci) {
        if (length(match_samples) > 0L) {
          for (x in seq_len(nids)) {
            val <- mean(source[[j]][x, match_samples], na.rm = TRUE)
            target[[j]][i, x] <- val
            target[[j]][x, i] <- val
          }
          target[[j]][i, i] <- mean(diag(source[[j]])[match_samples], na.rm = TRUE)
        } else {
          val <- mean(source[[j]], na.rm = TRUE)
          if (is.nan(val)) {
            val <- NA_real_
          }
          target[[j]][i, ] <- val
          target[[j]][, i] <- val
          target[[j]][i, i] <- mean(diag(target[[j]]), na.rm = TRUE)
        }
      }
    }
    target
  }

  # First pass: impute from the raw per-locus matrices
  D_imp <- lapply(D_list, function(x) x[,])
  D_imp <- impute_pass(D_list, D_imp, whichna)
  D_total1 <- Reduce("+", D_imp)
  whichna2 <- which(rowSums(is.na(D_total1)) > 0L)

  # Second pass: re-impute remaining rows using the first-pass matrices
  D_imp2 <- lapply(D_imp, function(x) x[,])
  D_imp2 <- impute_pass(D_imp, D_imp2, whichna2)
  D_total <- Reduce("+", D_imp2)

  D_total
}
