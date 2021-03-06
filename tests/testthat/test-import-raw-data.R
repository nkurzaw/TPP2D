## Tests for import2dDataset()
## library(TPP2D); library(testthat); source("setup-raw-2dtpp-dataset.R"); source("test-import-raw_data.R")
context("import data")

test_that("import2dDataset works as expected", {
  import_df <- import2dDataset(configTable = config_tab, 
                               data = raw_dat_list,
                               idVar = "protein_id",
                               intensityStr = "signal_sum_",
                               fcStr = "rel_fc_",
                               nonZeroCols = "qusm",
                               geneNameVar = "gene_name",
                               addCol = NULL,
                               qualColName = "qupm",
                               naStrs = c("NA", "n/d", "NaN"),
                               concFactor = 1e6,
                               medianNormalizeFC = TRUE,
                               filterContaminants = TRUE)
  recomp_sig_df <- data.frame(recomputeSignalFromRatios(import_df))
  expect_equivalent(data.frame(simulated_cell_extract_df), recomp_sig_df)
})
