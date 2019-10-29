#' Compute FDR for given F statistics based on true and
#' null dataset
#' 
#' @param df_out data frame containing results from analysis by
#' fitAndEvalDataset
#' @param df_null data frame containing results from analysis by
#' bootstrapNull
#' 
#' @return data frame annotating each protein with a FDR based on 
#' it's F statistic and number of observations
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:5)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' example_out <- fitAndEvalDataset(temp_df)
#' example_null <- bootstrapNull(temp_df, B = 1)
#' computeFdr(example_out, example_null)
#'  
#' @export
#'
#' @import dplyr
computeFdr <- function(df_out, df_null){

    dataset <- nObs <- nObsRound <- F_statistic <- 
    is_decoy <- max_rank <- true_cumsum <- 
    null_cumsum <- representative <- clustername <- 
    dataset <- fdr <- all_true <- all_null <- NULL
    
    B <- max(as.numeric(
    gsub("bootstrap_", "", unique(df_null$dataset))))
    
    out_df <- bind_rows(df_out %>% mutate(dataset = "true"),
                      df_null) %>%
    mutate(nObsRound = round(nObs/10)*10) %>%
    group_by(nObsRound) %>%
    arrange(desc(F_statistic)) %>%
    mutate(max_rank = n(),
           rank = dense_rank(desc(F_statistic)),
           is_decoy = ifelse(dataset != "true", 1, 0)) %>%
    mutate(all_null = sum(is_decoy),
           all_true = sum(!is_decoy),
           true_cumsum = cumsum(!is_decoy),
           null_cumsum = cumsum(is_decoy)) %>% 
    mutate(pi = (all_true-true_cumsum)/((all_null-null_cumsum)/B)) %>% 
    mutate(fdr = pi * (null_cumsum/B)/true_cumsum) %>% 
    ungroup()
    
    return(out_df)
}

#' Compute p-values for given F statistics based on true and
#' null dataset using a kernel density estimation of the 
#' bootstrapped null distribution
#' 
#' @param df_out data frame containing results from analysis by
#' fitAndEvalDataset
#' @param df_null data frame containing results from analysis by
#' bootstrapNull
#' 
#' @return data frame annotating each protein with a FDR based on 
#' it's F statistic and number of observations
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' example_out <- fitAndEvalDataset(temp_df)
#' example_null <- bootstrapNull(temp_df, B = 2)
#' computePvalFromKernelDensity(
#'   example_out, 
#'   example_null)
#'  
#' @export
#'
#' @import dplyr
#' @importFrom stats density
#' @importFrom stats p.adjust
#' @importFrom sfsmisc integrate.xy
computePvalFromKernelDensity <- function(df_out, df_null){
    
    dataset <- nObs <- nObsRound <- F_statistic <- 
        representative <- clustername <- dataset  <- NULL
    
    tmp_df <- bind_rows(df_out %>% mutate(dataset = "true"),
                        df_null) %>%
        mutate(nObsRound = round(nObs/10)*10) 
    
    out_df <- bind_rows(lapply(unique(tmp_df$nObsRound), function(nobs){
        nobs_df <- filter(tmp_df, nObsRound == nobs)
        density_est <- density(filter(nobs_df, dataset != "true")$F_statistic, 
                               n = 1e6)
        full_integral <- integrate.xy(density_est$x, density_est$y)
        nobs_df %>% filter(dataset == "true") %>% 
            group_by(representative) %>% 
            mutate(p_value = 1 - (
                integrate.xy(density_est$x[which(density_est$x < F_statistic)],
                             density_est$y[which(density_est$x < F_statistic)])/
                    full_integral)) %>% 
            ungroup
    })) %>% 
        ungroup %>% 
        within(p_value[p_value == 0] <- .Machine$double.eps) %>% 
        mutate(p_adj = p.adjust(p_value, method = "BH"))
    
    return(out_df)
}

#' Find hits according to FDR threshold
#' 
#' @param fdr_df data frame obtained from computeFdr
#' @param alpha significance threshold, default is set to 0.1
#' 
#' @return data frame of significant hits at FDR = alpha
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:5)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' example_out <- fitAndEvalDataset(temp_df)
#' example_null <- bootstrapNull(temp_df, B = 1)
#' fdr_df <- computeFdr(example_out, example_null)
#' findHits(fdr_df, 0.1)
#' 
#' @export
#' 
#' @import dplyr
findHits <- function(fdr_df, alpha){
  
  nObsRound <- fdr <- max_rank_fdr <- 
    dataset <- min_rank_true <- NULL
  
  hits_df <- fdr_df %>% 
    group_by(nObsRound) %>% 
    mutate(min_rank_true = min(rank[dataset == "true"])) %>% 
    filter(rank >= min_rank_true) %>% 
    mutate(max_rank_fdr = min(rank[fdr > alpha], na.rm = TRUE)) %>% 
    filter(rank < max_rank_fdr) %>% 
    filter(dataset == "true") %>% 
    ungroup()
  
  return(hits_df)
}