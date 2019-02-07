## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-compete-models.R")

test_that("competeModels works as expected", {
  cm_df <- competeModels(simulated_cell_extract_df %>% 
                           filter(clustername == "tp1"))
  expect_identical(round(cm_df$F_statistic, 2), 112.26)
})

test_that("competeModels throws error", {
  expect_error(competeModels(simulated_cell_extract_df %>% 
                               dplyr::select(-log2_value)))
})
