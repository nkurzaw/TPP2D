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


#' Resolve ambiguous protein names
#' 
#' @param df tidy data_frame retrieved after 
#' import of a 2D-TPP dataset
#' @param includeIsoforms logical indicating
#' whether protein isoform should be kept for
#' analysis
#' 
#' @return data frame with resolved protein name
#' ambiguity
#' 
#' @examples 
#' tst_df <- bind_rows(tibble(representative = rep(1:3, each = 3), 
#'                            clustername = rep(letters[1:3], each = 3)), 
#'                     tibble(representative = rep(c(4, 5), c(3, 2)), 
#'                            clustername = rep(c("a", "b"), c(3, 2))))
#'                            
#' resolveAmbigousProteinNames(tst_df)
#' @export
resolveAmbigousProteinNames <- function(df, includeIsoforms = FALSE){
    representative <- clustername <- nObs <- NULL
    lookUpDf <- df %>% 
        group_by(representative, clustername) %>% 
        summarize(nObs = n()) %>% 
        ungroup() 
    
    if(includeIsoforms){
        df_fil <- df %>% 
            rowwise() %>% 
            mutate(clustername = 
                       paste(representative,
                             clustername, sep = "_")) %>% 
            ungroup
    }else{
        lookUpDf <- lookUpDf %>% 
            filter(nObs == max(nObs)) %>% 
            filter(!duplicated(clustername))
        
        df_fil <- df %>% 
            filter(representative %in% 
                       lookUpDf$representative)
    }
    return(df_fil)
}