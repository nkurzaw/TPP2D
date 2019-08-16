#' @import dplyr
.removeAmbiguousMeasurements <- function(in_df, 
                                         qualColName = "qupm"){
    representative <- temperature <- conc <- value <- 
        raw_value <- NULL
    
    out_df <- in_df %>% 
        group_by(representative, temperature, conc) %>% 
        filter_(paste(qualColName, " == max(", qualColName, ")",
                      collapse = "")) %>% 
        filter(raw_value == max(raw_value)) %>% 
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
                factor_sd1 <- abs(log2(filter(temp_repr, temp_id == temp_repr$temp_id[i] & 
                                                  conc_id == temp_repr$conc_id[i])$rel_value) - 
                                      mean(log2(temp2$rel_value), na.rm = TRUE))/ 
                    sd(log2(temp2$rel_value), na.rm = TRUE)
                factor_sd2 <- abs(log2(filter(temp_repr, temp_id == temp_repr$temp_id[i] & 
                                                  conc_id == temp_repr$conc_id[i])$rel_value) - 
                                      mean(log2(temp2$rel_value), na.rm = TRUE))/ 
                    sd(log2(temp_repr$rel_value), na.rm = TRUE)
                temp3 <- filter(temp_repr, temp_id %in% 
                                    c(temp_repr$temp_id[i] - 1, 
                                      temp_repr$temp_id[i], 
                                      temp_repr$temp_id[i] + 1), 
                                conc_id %in% 
                                    c(temp_repr$conc_id[i] - 1, 
                                      temp_repr$conc_id[i], 
                                      temp_repr$conc_id[i] + 1))
                return(list("factor1" = factor_sd1,
                            "factor2" = factor_sd2,
                            "shrinked_value" = mean(temp3$rel_value, na.rm = TRUE),
                            "conc_edge" = (temp_repr$conc_id[i] == max(temp_repr$conc_id))))
            })
        if(!is.null(outlier_score)){
            temp_df <- temp_repr 
            temp_df$outlier_score_local = 
                sapply(outlier_score, function(x) x[["factor1"]])
            temp_df$outlier_score_global =
                sapply(outlier_score, function(x) x[["factor2"]])
            temp_df$shrinked_value = 
                sapply(outlier_score, function(x) x[["shrinked_value"]])
            temp_df$conc_edge = 
                sapply(outlier_score, function(x) x[["conc_edge"]])
        }else{
            temp_df <- temp_repr %>% 
                mutate(outlier_score_local = NA,
                       outlier_score_global = NA,
                       shrinked_value = NA,
                       conc_edge = NA)
        }
        
        return(temp_df)
    }))
    
    return(out_df)
}