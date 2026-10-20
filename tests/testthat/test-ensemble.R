test_that("original ensemble reproduces run.r on the Barratt fixture", {
  geno   <- readRDS(test_path("fixtures", "barratt_geno.rds"))
  ploidy <- readRDS(test_path("fixtures", "barratt_ploidy.rds"))
  expected <- readRDS(test_path("fixtures", "barratt_ensemble_expected.rds"))

  # the fixture was generated with the pure-R components (run.r verbatim);
  # feed the engine the same components -- rank mapping amplifies fp noise
  Dh <- euk_heuristic_distance(geno, ploidy, impute_data = TRUE, use_r = TRUE)
  Db <- euk_bayesian_distance(geno, ploidy, use_r = TRUE)
  D <- fasteukaryotyping:::euk_ensemble_distance_r(Dh, Db,
                                                 method = "original",
                                                 seed = 42)
  expect_equal(unname(D), unname(expected), tolerance = 1e-10)
})

test_that("original ensemble is deterministic under a fixed seed", {
  geno   <- readRDS(test_path("fixtures", "barratt_geno.rds"))
  ploidy <- readRDS(test_path("fixtures", "barratt_ploidy.rds"))
  D1 <- euk_ensemble_distance(geno, ploidy, method = "original", seed = 1)
  D2 <- euk_ensemble_distance(geno, ploidy, method = "original", seed = 1)
  expect_identical(D1, D2)
})

test_that("rank ensemble is divergence-oriented in both components", {
  # H and B perfectly opposed: pair12 least divergent by H but most similar
  # by B would collapse to 0.5 under the old direction bug; with the
  # corrected direction (ranking 2 - Db) ordering must follow divergence.
  Dh <- matrix(c(0, 1, 3, 1, 0, 2, 3, 2, 0), 3)
  Db <- matrix(c(2, 2.0, 0.1, 2.0, 2, 1.0, 0.1, 1.0, 2), 3)
  E <- fasteukaryotyping:::euk_ensemble_distance_r(Dh, Db, method = "rank")
  expect_true(E[1, 2] < E[2, 3])
  expect_true(E[2, 3] < E[1, 3])
  expect_true(is.na(E[1, 1]))
})

test_that("no comparable locus returns 1.0 in both backends", {
  geno <- list(matrix(c(NA_integer_, 1L, NA_integer_, 2L), nrow = 2))
  Dr   <- euk_bayesian_distance(geno, ploidy = 1L, use_r = TRUE)
  Dcpp <- euk_bayesian_distance(geno, ploidy = 1L, use_r = FALSE)
  expect_equal(Dr[1, 2], 1.0)
  expect_equal(Dcpp[1, 2], 1.0)
})

test_that("diploid singleton locus does not crash and matches backends", {
  # original script errors on if(NA); package extends with ll2 = ll1
  geno <- list(matrix(c(1L, 1L, 2L, NA), nrow = 2))
  Dr   <- euk_bayesian_distance(geno, ploidy = 2L, use_r = TRUE)
  Dcpp <- euk_bayesian_distance(geno, ploidy = 2L, use_r = FALSE)
  expect_equal(Dcpp, Dr, tolerance = 1e-10)
})
