## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-fits.R")

test_that("fitH0Model works as expected", {
  temp_df <- simulated_cell_extract_df %>% 
    filter(clustername == "tp1") %>% 
    mutate(nObs = n()) 
  fit_df <- fitH0Model(temp_df)
  expect_identical(round(fit_df$rssH0, 2), 4.83)
})

test_that("fitH1Model works as expected", {
  temp_df <- simulated_cell_extract_df %>% 
    filter(clustername == "tp1") %>% 
    mutate(nObs = n()) 
  fit_df <- fitH1Model(temp_df)
  expect_identical(round(fit_df$rssH1, 2), 0.07)
})
