## Tests for recomputeSignalFromRatios()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-recompute-signal.R")
context("plot functions")

test_that("plot2dTppProfile generates ggplot object", {
  gg_profile <- plot2dTppProfile(simulated_cell_extract_df, "protein1")
  expect_identical(class(gg_profile), c("ggplot2::ggplot", "ggplot", 
                                        "ggplot2::gg", "S7_object", "gg"))
})

test_that("plot2dTppProfile generates ggplot object", {
  gg_rel_profile <- plot2dTppRelProfile(simulated_cell_extract_df, "protein1")
  expect_identical(class(gg_rel_profile), c("ggplot2::ggplot", "ggplot", 
                                            "ggplot2::gg", "S7_object", "gg"))
})

test_that("plot2dTppFcHeatmap generates ggplot object", {
    gg_rel_heatmap <- plot2dTppFcHeatmap(
        simulated_cell_extract_df, "tp2", drug_name = "drug1")
    expect_identical(class(gg_rel_heatmap), c("ggplot2::ggplot", "ggplot", 
                                              "ggplot2::gg", "S7_object", "gg"))
})
