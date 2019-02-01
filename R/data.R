#' @title Example subset of a simulated 2D-TPP cell 
#' extract dataset 
#' @name simulated_cell_extract_df
#' @docType data
#' @description Simulated example dataset obtained by 2D-TPP 
#' experiments for analysis by the TPP2D-package. It contains 
#' a tidy data frame after import and recomputing of robust 
#' signal intensities with 200 simulated protein profiles 
#' (protein1-200) and 3 spiked-in true positives (TP1-3)
#' @format data frame with columns representative (protein id), 
#' clustername (gene name), temperature, log_conc, raw_value, 
#' rel_value, value and log2_value
NULL

#' @title Example raw data for a subset of a simulated 
#' 2D-TPP cell extract dataset 
#' @name raw_dat_list
#' @docType data
#' @description Simulated example dataset obtained by 2D-TPP 
#' experiments for analysis by the TPP2D-package. It contains 
#' a list of data frames resembling raw data files returned 
#' from a MS database search with 200 simulated protein profiles 
#' (protein1-200) and 3 spiked-in true positives (TP1-3).
#' @format list of data frames with columns representative 
#' (protein id), clustername (gene name), temperature, log_conc, 
#' raw_value, rel_value, value and log2_value
NULL

#' @title Example config table for a import of a simulated 
#' 2D-TPP cell extract dataset 
#' @name config_tab
#' @docType data
#' @description Config table fot import of simulated example 
#' dataset obtained by 2D-TPP experiments for analysis by 
#' the TPP2D-package. It's a data frame with the columns
#' "Compound" describing the compound used for the assay, 
#' "Experiment" listing MS experiment ids of the separate runs
#' (typically comprising two multiplexed adjacent temperature), 
#' "Temperature": the temperature used for a given sub-experimet, 
#' the respective TMT labels "126"-"131L", RefCol referring to
#' the label used as a reference label for computing relative
#' fold changes (usually the label used for the control treatment).
#' Please note that when the data is not supplied as a list of
#' already imported data frames the config table for the import 
#' function should be a path to an txt, csv or xlsx file containing
#' an additional column "Path" listing for each row the respective
#' path to a searched protein output file.
#' @format "Compound" describing the compound used for the assay, 
#' "Experiment" listing MS experiment ids of the separate runs
#' (typically comprising two multiplexed adjacent temperature), 
#' "Temperature": the temperature used for a given sub-experimet, 
#' the respective TMT labels "126"-"131L", RefCol referring to
#' the label used as a reference label for computing relative
#' fold changes (usually the label used for the control treatment).
NULL