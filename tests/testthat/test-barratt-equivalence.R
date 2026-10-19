test_that("euk_heuristic_distance matches original Barratt output", {
  geno   <- readRDS(test_path("fixtures", "barratt_geno.rds"))
  ploidy <- readRDS(test_path("fixtures", "barratt_ploidy.rds"))
  ids    <- readRDS(test_path("fixtures", "barratt_ids.rds"))
  expected <- readRDS(test_path("fixtures", "barratt_expected.rds"))

  D_r   <- euk_heuristic_distance(geno, ploidy, ids = ids, impute_data = TRUE, use_r = TRUE)
  D_cpp <- euk_heuristic_distance(geno, ploidy, ids = ids, impute_data = TRUE, use_r = FALSE)

  expect_equal(D_r,   expected, tolerance = 1e-10)
  expect_equal(D_cpp, expected, tolerance = 1e-10)
})
