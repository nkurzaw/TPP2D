## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-recompute-signal.R")

test_that("recomputeSignalFromRatios works as expected", {
  out <- recomputeSignalFromRatios(jq1_thp1_lys_subset)
  expect_identical(jq1_thp1_lys_subset, out)
})