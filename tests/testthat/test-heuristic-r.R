test_that("euk_heuristic_distance_r returns a symmetric non-negative matrix", {
  geno <- readRDS(test_path("fixtures", "geno_fixture.rds"))
  ploidy <- readRDS(test_path("fixtures", "ploidy_fixture.rds"))
  D <- fasteukaryotyping:::euk_heuristic_distance_r(geno, ploidy, impute_data = TRUE)

  n <- nrow(geno[[1L]])
  expect_equal(dim(D), c(n, n))
  expect_true(all(D >= 0, na.rm = TRUE))
  expect_equal(D, t(D))
})

test_that("euk_heuristic_distance wrapper uses the pure-R implementation", {
  eg <- example_eukaryotyping
  D <- euk_heuristic_distance(eg$geno, ploidy = eg$ploidy, ids = eg$ids)

  expect_equal(rownames(D), eg$ids)
  expect_equal(colnames(D), eg$ids)
  expect_true(all(D >= 0, na.rm = TRUE))
})

test_that("euk_heuristic_distance agrees with stored fixture", {
  eg <- example_eukaryotyping
  D <- euk_heuristic_distance(eg$geno, ploidy = eg$ploidy, ids = eg$ids, use_r = TRUE)
  expect_equal(D, eg$D, tolerance = 1e-14)
})

test_that("euk_import_hds converts an X-marked HDS correctly", {
  csv_path <- system.file("extdata", "example_hds.csv", package = "fasteukaryotyping")
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  # Only one locus called "Hap" for this smoke test
  out <- euk_import_hds(df, locinames = "Hap", ploidy = 2L)
  expect_true("Seq_ID" %in% names(out))
  expect_type(attr(out, "geno"), "list")
})
