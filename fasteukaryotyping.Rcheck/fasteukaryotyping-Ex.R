pkgname <- "fasteukaryotyping"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fasteukaryotyping')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("euk_bayesian_distance")
### * euk_bayesian_distance

flush(stderr()); flush(stdout())

### Name: euk_bayesian_distance
### Title: Bayesian Eukaryotyping distance
### Aliases: euk_bayesian_distance

### ** Examples

set.seed(42)
n <- 10
geno <- list(
  matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
  matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
)
d <- euk_bayesian_distance(geno, ploidy = c(2L, 2L))




cleanEx()
nameEx("euk_heuristic_distance")
### * euk_heuristic_distance

flush(stderr()); flush(stdout())

### Name: euk_heuristic_distance
### Title: Heuristic Eukaryotyping distance
### Aliases: euk_heuristic_distance

### ** Examples

set.seed(42)
n <- 10
geno <- list(
  matrix(sample(c(1L, 2L, 3L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2),
  matrix(sample(c(1L, 2L, NA), n * 2, replace = TRUE), nrow = n, ncol = 2)
)
d <- euk_heuristic_distance(geno, ploidy = c(2L, 2L))




cleanEx()
nameEx("euk_import_hds")
### * euk_import_hds

flush(stderr()); flush(stdout())

### Name: euk_import_hds
### Title: Import a Haplotype Data Sheet (HDS) into the internal geno list
### Aliases: euk_import_hds

### ** Examples

# See inst/extdata/Example_haplotype_data_sheet.txt



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
