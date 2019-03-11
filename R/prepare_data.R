#' Recompute robust signal intensities based on 
#' bootstrapped TMT channel ratios
#' 
#' @param df tidy data_frame retrieved after 
#' import of a 2D-TPP dataset
#' 
#' @return A data_frame with recomputed signal 
#' intensities (columname: value) and log2 
#' transformed signal intensities (columnanme: log2_value) 
#' that more reliably reflect relative ratios 
#' between the TMT channels
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' recomputeSignalFromRatios(simulated_cell_extract_df)
#' 
#' @import dplyr
recomputeSignalFromRatios <- function(df){
  clustername <- log2_value <- raw_value <- rel_value <- 
    representative <- sum_raw <- sum_rel <- temperature <- 
    value <- NULL
  
  df_recomp <- df %>%
    group_by(representative, clustername, temperature) %>%
    mutate(sum_raw = sum(raw_value, na.rm = TRUE),
           sum_rel = sum(rel_value, na.rm = TRUE)) %>%
    mutate(value = (rel_value/sum_rel)*sum_raw) %>%
    ungroup() %>%
    dplyr::select(-sum_raw, -sum_rel) %>%
    mutate(log2_value = log2(value)) %>%
    filter(!is.na(log2_value))
  
  return(df_recomp)
}