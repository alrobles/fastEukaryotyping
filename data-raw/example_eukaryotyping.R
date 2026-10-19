## code to prepare `example_eukaryotyping` and test fixtures goes here
## Run from the package root.

devtools::load_all()

set.seed(42L)

n <- 12L
nloci <- 3L
ploidy <- c(2L, 2L, 1L)

# Geno list: 3 loci, diploid for the first two and haploid for the third.
geno <- list(
  matrix(c(
    1L, 1L, 2L, 2L, 3L, NA, 1L, 2L, 3L, 4L, 1L, 2L,
    1L, 2L, 2L, 1L, NA, 1L, 2L, 1L, NA, 2L, 2L, 1L
  ), nrow = n, ncol = 2L),
  matrix(c(
    1L, 2L, NA, 1L, 2L, 2L, 1L, 1L, NA, 2L, 1L, 2L,
    2L, 1L, 1L, NA, 2L, 1L, 2L, NA, 1L, 1L, 2L, 2L
  ), nrow = n, ncol = 2L),
  matrix(c(
    1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L
  ), nrow = n, ncol = 1L)
)

ids <- paste0("S", seq_len(n))

D <- euk_heuristic_distance_r(geno, ploidy, impute_data = TRUE)
rownames(D) <- ids
colnames(D) <- ids

# X-marked Haplotype Data Sheet example
hds <- data.frame(
  Seq_ID = ids,
  Hap_1 = c("X", "", "X", "", "X", "X", "", "", "X", "", "X", ""),
  Hap_2 = c("", "X", "", "X", "", "", "X", "X", "", "X", "", "X"),
  stringsAsFactors = FALSE
)

example_eukaryotyping <- list(
  geno    = geno,
  ploidy  = ploidy,
  ids     = ids,
  D       = D,
  hds     = hds
)

# Package example data (must be .rda for the data/ directory)
save(example_eukaryotyping, file = "data/example_eukaryotyping.rda")

# Test fixtures
saveRDS(geno,   "tests/testthat/fixtures/geno_fixture.rds")
saveRDS(ploidy, "tests/testthat/fixtures/ploidy_fixture.rds")
saveRDS(ids,    "tests/testthat/fixtures/ids_fixture.rds")
saveRDS(D,      "tests/testthat/fixtures/expected_fixture.rds")

# CSV example for euk_import_hds
write.csv(hds, "inst/extdata/example_hds.csv", row.names = FALSE)
