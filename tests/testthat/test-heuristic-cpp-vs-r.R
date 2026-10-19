TOL <- 1e-10

test_that("euk_heuristic_distance_cpp matches euk_heuristic_distance_r on fixture", {
  eg <- example_eukaryotyping
  D_r   <- euk_heuristic_distance(eg$geno, eg$ploidy, ids = eg$ids, use_r = TRUE)
  D_cpp <- euk_heuristic_distance(eg$geno, eg$ploidy, ids = eg$ids, use_r = FALSE)
  expect_equal(D_cpp, D_r, tolerance = TOL)
})

test_that("euk_heuristic_distance_cpp matches euk_heuristic_distance_r on random data", {
  set.seed(123L)
  for (trial in 1:5) {
    n      <- sample(8:30, 1)
    L      <- sample(2:5, 1)
    ploidy <- sample(c(1L, 2L), L, replace = TRUE)
    m      <- sample(1:2, L, replace = TRUE)

    geno <- lapply(seq_len(L), function(j) {
      mois <- sample(1:(ploidy[j] + 1L), n, replace = TRUE)
      mat  <- matrix(0L, nrow = n, ncol = max(mois))
      for (i in seq_len(n)) {
        if (mois[i] > 0L) {
          mat[i, seq_len(mois[i])] <- sample(1:6, mois[i], replace = TRUE)
        }
      }
      mat
    })

    D_r   <- euk_heuristic_distance(geno, ploidy, m = m, use_r = TRUE)
    D_cpp <- euk_heuristic_distance(geno, ploidy, m = m, use_r = FALSE)
    expect_equal(D_cpp, D_r, tolerance = TOL, info = paste("trial", trial))
  }
})
