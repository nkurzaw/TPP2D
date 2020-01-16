## Tests for getModelParamsDf()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-getModelParamsDf.R")
context("getModelParamsDf")

params_df <- getModelParamsDf(.minObsFilter(simulated_cell_extract_df) %>% 
                                  filter(representative %in% 1:3))

test_that("getModelParamsDf: rss behave as expected", {
    expect_true(all(params_df$rssH0 - params_df$rssH1 > 0))
})

test_that("getModelParamsDf: length of residuals corresponds to obs.", {
    expect_true(all(params_df$nObs == sapply(params_df$residualsH0, length)))
    expect_true(all(params_df$nObs == sapply(params_df$residualsH1, length)))
})

