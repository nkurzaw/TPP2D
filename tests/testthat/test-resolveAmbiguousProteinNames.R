## Tests for resolveAmbiguousProteinNames()
## library(TPP2D); library(testthat); source("setup-2dtpp-dataset.R"); source("test-recompute-signal.R")
context("resolve ambiguous protein names")

tst_df <- bind_rows(
    simulated_cell_extract_df,
    filter(simulated_cell_extract_df, 
           clustername == "tp3",
           temperature == 62) %>% 
        mutate(representative = "99")
)

test_that("resolveAmbiguousProteinNames works as expected", {
    out_df <- resolveAmbiguousProteinNames(tst_df)
    expect_equal(nrow(filter(out_df, clustername == "tp3")), 50)
})

test_that(paste("resolveAmbiguousProteinNames works as expected", 
                "with includeIsoforms = FALSE"), {
    out_df <- resolveAmbiguousProteinNames(
        tst_df, includeIsoforms = TRUE)
    expect_equal(nrow(filter(out_df, clustername == "20_tp3")), 50)
})