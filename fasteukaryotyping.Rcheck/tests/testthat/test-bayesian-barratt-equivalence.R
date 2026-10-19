test_that("euk_bayesian_distance matches original Plucinski script output", {
  geno   <- readRDS(test_path("fixtures", "barratt_geno.rds"))
  ploidy <- readRDS(test_path("fixtures", "barratt_ploidy.rds"))
  ids    <- readRDS(test_path("fixtures", "barratt_ids.rds"))
  expected <- readRDS(test_path("fixtures", "barratt_bayesian_expected.rds"))

  D_r   <- euk_bayesian_distance(geno, ploidy, ids = ids, use_r = TRUE)
  D_cpp <- euk_bayesian_distance(geno, ploidy, ids = ids, use_r = FALSE)

  expect_equal(D_r,   expected, tolerance = 1e-10)
  expect_equal(D_cpp, expected, tolerance = 1e-10)
})
