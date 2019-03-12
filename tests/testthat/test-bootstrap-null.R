## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-bootstrap-null.R")
context("bootstrap null")

test_that("competeModels works as expected", {
  boot_df <- bootstrapNull(simulated_cell_extract_df %>% 
                             filter(clustername == "tp1"),
                           B = 3)
  expect_identical(nrow(boot_df), 3)
  expect_gt(boot_df$F_statistic, 0)
})