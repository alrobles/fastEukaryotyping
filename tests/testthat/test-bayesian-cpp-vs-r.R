TOL <- 1e-10

test_that("euk_bayesian_distance_cpp matches euk_bayesian_distance_r on fixture", {
  eg <- example_eukaryotyping
  D_r   <- euk_bayesian_distance(eg$geno, eg$ploidy, ids = eg$ids, use_r = TRUE)
  D_cpp <- euk_bayesian_distance(eg$geno, eg$ploidy, ids = eg$ids, use_r = FALSE)
  expect_equal(D_cpp, D_r, tolerance = TOL)
})

test_that("euk_bayesian_distance_cpp matches euk_bayesian_distance_r on random data", {
  set.seed(321L)
  for (trial in 1:5) {
    n      <- sample(8:25, 1)
    L      <- sample(2:4, 1)
    ploidy <- sample(c(1L, 2L), L, replace = TRUE)

    geno <- lapply(seq_len(L), function(j) {
      mois <- sample(1:(ploidy[j] + 1L), n, replace = TRUE)
      mat  <- matrix(0L, nrow = n, ncol = max(mois))
      for (i in seq_len(n)) {
        if (mois[i] > 0L) {
          mat[i, seq_len(mois[i])] <- sample(1:5, mois[i], replace = TRUE)
        }
      }
      mat
    })

    D_r   <- euk_bayesian_distance(geno, ploidy, use_r = TRUE)
    D_cpp <- euk_bayesian_distance(geno, ploidy, use_r = FALSE)
    expect_equal(D_cpp, D_r, tolerance = TOL)
  }
})

test_that("bayesian distances are symmetric and within [0, 2]", {
  eg <- example_eukaryotyping
  D <- euk_bayesian_distance(eg$geno, eg$ploidy, use_r = TRUE)
  expect_true(isSymmetric(unname(D), tol = TOL))
  expect_true(all(D[!is.na(D)] >= 0 & D[!is.na(D)] <= 2))
})

test_that("identical samples score higher similarity than divergent ones", {
  # distance is posterior E[#shared hypotheses] on {0,1,2}: larger = closer
  geno <- list(matrix(c(1L, 1L, 3L, 4L), nrow = 2, ncol = 2))
  D <- euk_bayesian_distance(geno, ploidy = 1L)
  expect_true(D[1, 1] >= D[1, 2])
})

test_that("euk_ensemble_distance rank method returns means in [0, 1]", {
  eg <- example_eukaryotyping
  D <- euk_ensemble_distance(eg$geno, eg$ploidy, ids = eg$ids,
                             method = "rank")
  expect_true(isSymmetric(unname(D), tol = TOL))
  off <- D[upper.tri(D)]
  expect_true(all(off >= 0 & off <= 1, na.rm = TRUE))
})

test_that("original method may be asymmetric (random ties, as in run.r)", {
  # documented quirk of the original procedure, not a bug: rank() with
  # ties.method = "random" gives H[i,j] and H[j,i] different ranks
  eg <- example_eukaryotyping
  D <- euk_ensemble_distance(eg$geno, eg$ploidy, method = "original",
                             seed = 42)
  expect_true(is.matrix(D))
  expect_true(all(D >= 0 & D <= 1))
})
