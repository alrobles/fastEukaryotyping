#' Import a Haplotype Data Sheet (HDS) into the internal geno list
#'
#' @param data A `data.frame` with a `Seq_ID` column and one column per
#'   haplotype per locus. Presence is marked by `"X"`.
#' @param locinames Character vector of marker names (as used in column names).
#' @param locinames_base Optional character vector of base marker names for
#'   markers split into sub-markers. Defaults to `locinames`.
#' @param ploidy Integer vector (1 = haploid, 2 = diploid) for each locus in
#'   `locinames`.
#' @param include_seq_id Logical. Should the `Seq_ID` column be kept in the
#'   output? Default `TRUE`.
#' @return A `data.frame` with `Seq_ID` and columns named
#'   `<slot>_Hap_<lociname>`, plus a `geno` attribute containing the list of
#'   per-locus integer matrices.
#' @export
#' @examples
#' # See inst/extdata/Example_haplotype_data_sheet.txt
euk_import_hds <- function(data, locinames, locinames_base = locinames,
                           ploidy = rep(2L, length(locinames)),
                           include_seq_id = TRUE) {
  checkmate::assert_data_frame(data)
  checkmate::assert_character(locinames, min.len = 1L, any.missing = FALSE)
  checkmate::assert_character(locinames_base, len = length(locinames))
  checkmate::assert_integerish(ploidy, len = length(locinames))
  checkmate::assert_flag(include_seq_id)

  if (!"Seq_ID" %in% names(data)) {
    stop("'data' must contain a 'Seq_ID' column.")
  }

  data <- data[!is.na(data$Seq_ID) & data$Seq_ID != "", ]
  ids <- data$Seq_ID
  nids <- length(ids)
  nloci <- length(locinames)

  newdata <- NULL
  datacompleteness <- NULL
  geno <- vector("list", nloci)

  for (j in seq_len(nloci)) {
    locicolumns <- grepl(paste0("^", locinames[j]), names(data)) |
                   grepl(locinames[j], names(data))
    if (!any(locicolumns)) {
      stop(sprintf("No columns matching locus '%s' found.", locinames[j]))
    }
    raw_alleles <- data[, locicolumns, drop = FALSE]
    moi <- rowSums(raw_alleles == "X", na.rm = TRUE)
    max_moi <- max(c(2L, moi))
    newdatatemp <- matrix(NA_integer_, nrow = nids, ncol = max_moi)

    for (x in seq_len(nids)) {
      if (moi[x] > 0L) {
        hap_idx <- which(raw_alleles[x, ] == "X")
        nh <- min(length(hap_idx), max_moi)
        newdatatemp[x, seq_len(nh)] <- hap_idx[seq_len(nh)]
      }
    }

    if (ploidy[j] > 1L) {
      for (x in seq_len(nids)) {
        if (moi[x] == 1L) {
          newdatatemp[x, 1L:2L] <- newdatatemp[x, 1L]
        }
      }
    }

    colnames(newdatatemp) <- paste0(seq_len(max_moi), "_Hap_", locinames[j])
    newdata <- cbind(newdata, newdatatemp)
    datacompleteness <- cbind(datacompleteness, moi)
    geno[[j]] <- newdatatemp
  }

  out <- data.frame(ids = ids, newdata, stringsAsFactors = FALSE)
  if (include_seq_id) {
    names(out)[1L] <- "Seq_ID"
  } else {
    out <- out[, -1L, drop = FALSE]
  }
  attr(out, "geno") <- geno
  out
}

#' Heuristic Eukaryotyping distance
#'
#' Computes Barratt's unsupervised heuristic distance matrix for MLST haplotype
#' data sheets.
#'
#' @param geno List of `n x maxMOI` integer matrices (one per locus). `NA` marks
#'   empty slots.
#' @param ploidy Integer vector (1 = haploid, 2 = diploid) for each locus.
#' @param ids Optional character vector of sample identifiers.
#' @param m Integer vector of `m` values, one per locus. Default `rep(1L, length(geno))`.
#' @param impute_data Logical. If `TRUE` (default), missing loci are imputed
#'   with the two-pass Barratt mean procedure so that the returned matrix is
#'   fully populated. If `FALSE`, pairs sharing a missing locus return `NA`.
#' @param use_r Logical. If `TRUE` use the pure-R reference implementation
#'   (deprecated); if `FALSE` (default) use the optimised C++ implementation.
#' @return An `n x n` numeric distance matrix.
#' @export
#' @examples
#' set.seed(42)
#' n <- 10
#' geno <- list(
#'   matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
#'   matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
#' )
#' d <- euk_heuristic_distance(geno, ploidy = c(2L, 2L))
#'
euk_heuristic_distance <- function(geno, ploidy, ids = NULL, m = rep(1L, length(geno)),
                                   impute_data = TRUE, use_r = FALSE) {
  checkmate::assert_list(geno, min.len = 1L)
  checkmate::assert_integerish(ploidy, len = length(geno))
  checkmate::assert_integerish(m, len = length(geno))
  checkmate::assert_character(ids, len = nrow(geno[[1L]]), null.ok = TRUE)
  checkmate::assert_flag(impute_data)
  checkmate::assert_flag(use_r)

  n <- nrow(geno[[1L]])
  for (j in seq_along(geno)) {
    checkmate::assert_matrix(geno[[j]], mode = "integerish", nrows = n)
  }

  if (use_r) {
    message("DeprecationWarning: use_r = TRUE is deprecated; the C++ backend is now the default.")
    D <- euk_heuristic_distance_r(geno, ploidy, m, impute_data)
  } else {
    D <- .euk_heuristic_distance_impl(geno, ploidy, as.integer(m), as.integer(impute_data))
  }

  if (!is.null(ids)) {
    rownames(D) <- ids
    colnames(D) <- ids
  }
  D
}


#' Bayesian Eukaryotyping distance
#'
#' Computes Plucinski's unsupervised naive-Bayes distance matrix for MLST
#' haplotype data sheets.
#'
#' @param geno List of `n x maxMOI` integer matrices (one per locus). `NA` marks
#'   empty slots.
#' @param ploidy Integer vector (1 = haploid, 2 = diploid) for each locus.
#' @param ids Optional character vector of sample identifiers.
#' @param epsilon Rate of missed alleles (default `0.3072`).
#' @param use_r Logical. If `TRUE` use the pure-R reference implementation
#'   (deprecated); if `FALSE` (default) use the optimised C++ implementation.
#' @return An `n x n` numeric distance matrix in \eqn{[0, 2]}.
#' @export
#' @examples
#' set.seed(42)
#' n <- 10
#' geno <- list(
#'   matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
#'   matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
#' )
#' d <- euk_bayesian_distance(geno, ploidy = c(2L, 2L))
#'
euk_bayesian_distance <- function(geno, ploidy, ids = NULL,
                                  epsilon = 0.3072, use_r = FALSE) {
  checkmate::assert_list(geno, min.len = 1L)
  checkmate::assert_integerish(ploidy, len = length(geno))
  checkmate::assert_character(ids, len = nrow(geno[[1L]]), null.ok = TRUE)
  checkmate::assert_number(epsilon, lower = 0, upper = 1)
  checkmate::assert_flag(use_r)

  n <- nrow(geno[[1L]])
  for (j in seq_along(geno)) {
    checkmate::assert_matrix(geno[[j]], mode = "integerish", nrows = n)
  }

  if (use_r) {
    message("DeprecationWarning: use_r = TRUE is deprecated; the C++ backend is now the default.")
    D <- euk_bayesian_distance_r(geno, ploidy, epsilon)
  } else {
    D <- .euk_bayesian_distance_impl(geno, ploidy, epsilon)
  }

  if (!is.null(ids)) {
    rownames(D) <- ids
    colnames(D) <- ids
  }
  D
}


#' Eukaryotyping ensemble distance
#'
#' Combines the heuristic and Bayesian distance matrices into a single
#' ensemble. Two rules are available: `"original"` (default) reproduces the
#' distribution-mapping procedure of the original `run.r` script shipped
#' with Barratt et al. (2019); `"rank"` is a package-specific
#' direction-consistent fractional-rank mean.
#'
#' @param geno List of `n x maxMOI` integer matrices (one per locus).
#' @param ploidy Integer vector (1 = haploid, 2 = diploid) for each locus.
#' @param ids Optional character vector of sample identifiers.
#' @param method Combination rule. `"original"` (default) reproduces the
#'   combining procedure of the original `run.r` script: the Bayesian matrix
#'   is reversed (`2 - D_bayesian`), both matrices are min-max normalised
#'   over all entries including the diagonal, the sorted Bayesian values are
#'   mapped onto the ranks of the heuristic matrix, the ensemble is the
#'   elementwise mean of the mapped and heuristic matrices, rescaled by its
#'   maximum. Random tie breaking makes this method non-deterministic
#'   without `seed` and can yield a non-symmetric matrix. `"rank"` is a
#'   package-specific fractional-rank mean over the strictly upper
#'   triangle, with the Bayesian matrix reversed so both components are
#'   divergence-oriented; its diagonal is `NA`.
#' @param seed Optional integer seed for the `"original"` method's random
#'   tie breaking. Ignored for `method = "rank"`.
#' @param epsilon Rate of missed alleles for the Bayesian component.
#' @param ... Passed to [euk_heuristic_distance()] (e.g. `m`, `impute_data`,
#'   `use_r`).
#' @return An `n x n` matrix of ensemble distances (in \eqn{[0, 1]} for
#'   `method = "rank"` off-diagonal entries). Larger values indicate pairs
#'   that both components rank as more divergent.
#' @export
euk_ensemble_distance <- function(geno, ploidy, ids = NULL,
                                  method = c("original", "rank"),
                                  seed = NULL, epsilon = 0.3072, ...) {
  method <- match.arg(method)
  Dh <- euk_heuristic_distance(geno, ploidy, ids = NULL, ...)
  Db <- euk_bayesian_distance(geno, ploidy, ids = NULL, epsilon = epsilon)
  D <- euk_ensemble_distance_r(Dh, Db, method = method, seed = seed)
  if (!is.null(ids)) {
    rownames(D) <- ids
    colnames(D) <- ids
  }
  D
}
