## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-recompute-signal.R")
context("plot functions")

test_that("plot2dTppProfile generates ggplot object", {
  gg_profile <- plot2dTppProfile(simulated_cell_extract_df, "protein1")
  expect_identical(class(gg_profile), c("gg", "ggplot"))
})

test_that("plot2dTppProfile generates ggplot object", {
  gg_rel_profile <- plot2dTppRelProfile(simulated_cell_extract_df, "protein1")
  expect_identical(class(gg_rel_profile), c("gg", "ggplot"))
})
