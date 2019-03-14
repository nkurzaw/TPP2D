## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-fits.R")
context("fit models")

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
  expect_identical(round(fit_df$rssH1, 1), 0.1)
})

test_that("computeFstat works as expected", {
  temp_df <- simulated_cell_extract_df %>% 
    filter(clustername == "tp1") %>% 
    mutate(nObs = n()) 
  h0_df <- fitH0Model(temp_df)
  h1_df <- fitH1Model(temp_df)
  f_df <- computeFstat(h0_df, h1_df)
  expect_identical(round(f_df$F_statistic, -2), 100)
})

test_that("fitAndEvalDataset works as expected", {
  temp_df <- simulated_cell_extract_df %>% 
    filter(clustername == "tp1") %>% 
    mutate(nObs = n()) 
  f_df <- fitAndEvalDataset(temp_df)
  expect_identical(round(f_df$F_statistic, -2), 100)
})
