## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-compute-F-stat-from-params.R")
context("compute F stat from params")

test_that("computeFStatFromParams works as expected", {
    params_df <- getModelParamsDf(simulated_cell_extract_df %>%
                                      filter(representative %in% 1:3))
    fstat_df <- computeFStatFromParams(params_df)
    expect_equal(round(fstat_df$F_statistic, 2), c(1.23, 0.36, 0.24))
})