## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-recompute-signal.R")

test_that("recomputeSignalFromRatios works as expected", {
  out <- recomputeSignalFromRatios(simulated_cell_extract_df)
  expect_equivalent(simulated_cell_extract_df, out)
})
