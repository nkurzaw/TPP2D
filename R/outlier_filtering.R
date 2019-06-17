#' @import dplyr
.removeAmbiguousMeasurements <- function(in_df, 
                                         qualColName = "qupm"){
    representative <- temperature <- conc <- value <- NULL
    
    out_df <- in_df %>% 
        group_by(representative, temperature, conc) %>% 
        filter_(paste(qualColName, " == max(", qualColName, ")",
                      collapse = "")) %>% 
        filter(value == max(value)) %>% 
        group_by(representative) %>% 
        filter_(paste("any(", qualColName, " > 1)", collapse = "")) %>% 
        arrange(temperature, conc) %>% 
        mutate(temp_id = dense_rank(temperature),  
               conc_id = dense_rank(conc)) %>% 
        ungroup
    
    return(out_df)
}

#' @import dplyr
#' @importFrom stats sd
.detectOutliers <- function(in_df){
    representative <- temp_id <- conc_id <- NULL
    
    out_df <- bind_rows(lapply(unique(in_df$representative), function(repr){
        temp_repr <- filter(in_df, representative == repr)
        outlier_score <- 
            lapply(seq_len(nrow(temp_repr)), function(i){
                temp2 <- filter(temp_repr, temp_id %in% 
                                    c(temp_repr$temp_id[i] - 1, 
                                      temp_repr$temp_id[i], 
                                      temp_repr$temp_id[i] + 1), 
                                conc_id %in% 
                                    c(temp_repr$conc_id[i] - 1, 
                                      temp_repr$conc_id[i], 
                                      temp_repr$conc_id[i] + 1), 
                                !(temp_id == temp_repr$temp_id[i] & 
                                      conc_id == temp_repr$conc_id[i]))
                factor_sd <- abs(log2(filter(temp_repr, temp_id == temp_repr$temp_id[i] & 
                                                 conc_id == temp_repr$conc_id[i])$rel_value) - 
                                     mean(log2(temp2$rel_value), na.rm = TRUE))/ 
                    sd(log2(temp2$rel_value), na.rm = TRUE)
                temp3 <- filter(temp_repr, temp_id %in% 
                                    c(temp_repr$temp_id[i] - 1, 
                                      temp_repr$temp_id[i], 
                                      temp_repr$temp_id[i] + 1), 
                                conc_id %in% 
                                    c(temp_repr$conc_id[i] - 1, 
                                      temp_repr$conc_id[i], 
                                      temp_repr$conc_id[i] + 1))
                return(list("factor" = factor_sd,
                            "shrinked_value" = mean(temp3$rel_value, na.rm = TRUE)))
            })
        if(!is.null(outlier_score)){
            temp_df <- temp_repr 
            temp_df$outlier_score = 
                sapply(outlier_score, function(x) x[["factor"]])
            temp_df$shrinked_value = 
                sapply(outlier_score, function(x) x[["shrinked_value"]])
        }else{
            temp_df <- temp_repr %>% 
                mutate(outlier_score = NA,
                       shrinked_value = NA)
        }
        
        return(temp_df)
    }))
    
    return(out_df)
}

#' Moderate outlier measurements
#' 
#' @param df tidy data frame retrieved after import of a 2D-TPP 
#' dataset
#' @param outlier_quantile quantile that should be used as a cutoff to 
#' define outliers, default: 0.98
#' @param qualColName column indicating how many unique peptides 
#' were used for protein quantification
#' 
#' @return data frame containing outlier moderated 2D-TPP data
#' 
#' @examples
#' data("simulated_cell_extract_df")
#' moderateOutliers(simulated_cell_extract_df) 
#' @export
#' @import dplyr
moderateOutliers <- function(df, 
                             outlier_quantile = 0.98,
                             qualColName = "qupm"){
    
    filtered_df <- .removeAmbiguousMeasurements(df, 
                                                qualColName = qualColName)
    
    out_detected_df <- .detectOutliers(in_df = filtered_df)
    
    moderated_df <- 
        mutate(out_detected_df, rel_value = 
                   case_when(outlier_score < 
                                 quantile(out_detected_df$outlier_score, 
                                          outlier_quantile, na.rm = TRUE) ~ rel_value,
                             TRUE ~ shrinked_value))
    out_df <- recomputeSignalFromRatios(moderated_df)
        
    return(out_df)
}