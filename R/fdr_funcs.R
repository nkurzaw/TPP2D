#' Get FDR for given F statistics based on true and
#' null dataset 
#' 
#' @param df_out data frame containing results from analysis by
#' fitAndEvalDataset
#' @param df_null data frame containing results from analysis by
#' bootstrapNull
#' @param squeezeDenominator logical indicating whether F statistic
#' denominator should be shrinked using limma::squeezeVar
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
#' getFDR(example_out, example_null)
#'  
#' @export
#'
#' @import dplyr
getFDR <- function(df_out, df_null, squeezeDenominator = TRUE){
    dataset <- nObs <- nObsRound <- F_statistic <- 
        is_decoy <- max_rank <- true_cumsum <- 
        null_cumsum <- representative <- clustername <- 
        dataset <- FDR <- all_true <- all_null <- NULL
    
    if(squeezeDenominator){
        df_out <- .shrinkFstat(df_out, trueOrNull = "true")
        df_null <- .shrinkFstat(df_null, trueOrNull = "null")
    }
    
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
        mutate(FDR = pi * (null_cumsum/B)/true_cumsum) %>% 
        ungroup() %>% 
        within(FDR[is.na(F_statistic)] <- NA)
    
    return(out_df)
}

#' @importFrom limma squeezeVar 
.shrinkFstat <- function(inDf, trueOrNull = "true"){
    rssH1 <- df2 <- rssH0 <- rssH1Squeezed <- df1 <- 
        nObs <- nObsRound <- dataset <- NULL
    if(trueOrNull == "true"){
        outDf <- inDf %>% 
            mutate(nObsRound = round(nObs/10)*10) %>% 
            group_by(nObsRound) %>% 
            mutate(rssH1Squeezed = limma::squeezeVar(
                rssH1, df = df2)$var.post) %>% 
            ungroup %>% 
            mutate(F_statistic = (rssH0 - rssH1)/
                       (rssH1Squeezed) * df2/df1)
    }else if(trueOrNull == "null"){
        outDf <- inDf %>% 
            mutate(nObsRound = round(nObs/10)*10) %>% 
            group_by(dataset, nObsRound) %>% 
            mutate(rssH1Squeezed = limma::squeezeVar(
                rssH1, df = df2)$var.post) %>% 
            ungroup %>% 
            mutate(F_statistic = (rssH0 - rssH1)/
                       (rssH1Squeezed) * df2/df1)
    }
    return(outDf)
}
#' Compute FDR for given F statistics based on true and
#' null dataset (old function)
#' 
#' @param df_out data frame containing results from analysis by
#' fitAndEvalDataset
#' @param df_null data frame containing results from analysis by
#' bootstrapNull
#' 
#' @return data frame annotating each protein with a FDR based on 
#' it's F statistic and number of observations
#' 
#' @name computeFdr-deprecated
#' @seealso \code{\link{TPP2D-deprecated}}
#' @keywords internal
NULL

#' @rdname TPP2D-deprecated
#' @section \code{computeFdr}:
#' For \code{computeFdr}, use \code{\link{getFDR}}.
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
    .Deprecated("getFDR")
    dataset <- nObs <- nObsRound <- F_statistic <- 
    is_decoy <- max_rank <- true_cumsum <- 
    null_cumsum <- representative <- clustername <- 
    dataset <- FDR <- all_true <- all_null <- NULL
    
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
    mutate(FDR = pi * (null_cumsum/B)/true_cumsum) %>% 
    ungroup()
    
    return(out_df)
}

#' Compute p-values for given F statistics based on true and
#' null dataset 
#' 
#' @param df_out data frame containing results from analysis by
#' fitAndEvalDataset
#' @param df_null data frame containing results from analysis by
#' bootstrapNull
#' @param squeezeDenominator logical indicating whether F statistic
#' denominator should be shrinked using limma::squeezeVar
#' 
#' @return data frame annotating each protein with a FDR based on 
#' it's F statistic and number of observations
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:3)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' example_out <- fitAndEvalDataset(temp_df)
#' example_null <- bootstrapNull(temp_df, B = 2)
#' getPvalues(
#'   example_out, 
#'   example_null)
#'  
#' @export
#'
#' @import dplyr
#' @importFrom stats density
#' @importFrom stats p.adjust
getPvalues <- function(df_out, df_null,
                       squeezeDenominator = FALSE){
    dataset <- nObs <- nObsRound <- F_statistic <- 
        is_decoy <- max_rank <- true_cumsum <- 
        null_cumsum <- representative <- clustername <- 
        dataset <- FDR <- all_true <- all_null <- NULL
    
    if(squeezeDenominator){
        df_out <- .shrinkFstat(df_out, trueOrNull = "true")
        df_null <- .shrinkFstat(df_null, trueOrNull = "null")
    }
    
    out_df <- bind_rows(df_out %>% mutate(dataset = "true"),
                        df_null) %>%
        mutate(nObsRound = round(nObs/10)*10) %>%
        group_by(nObsRound) %>%
        arrange(desc(F_statistic)) %>%
        mutate(max_rank = n(),
               rank = dense_rank(desc(F_statistic)),
               is_decoy = ifelse(dataset != "true", 1, 0)) %>%
        mutate(all_null = sum(is_decoy),
               null_cumsum = cumsum(is_decoy)) %>% 
        mutate(p_value = (null_cumsum)/all_null) %>% 
        ungroup() %>% 
        within(p_value[p_value == 0] <- .Machine$double.eps) %>% 
        filter(dataset == "true") %>% 
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
#' fdr_df <- getFDR(example_out, example_null)
#' findHits(fdr_df, 0.1)
#' 
#' @export
#' 
#' @import dplyr
findHits <- function(fdr_df, alpha){
  
  nObsRound <- FDR <- max_rank_fdr <- dataset <- min_rank_true <- 
      representative <- clustername <- nObs <- nCoeffsH0 <- 
      nCoeffsH1 <- rssH0 <- rssH1 <- df1 <- df2 <- F_statistic <- 
      pEC50H1 <- slopeH1 <- pEC50_slopeH1 <- detected_effectH1 <- 
      df_fil <- NULL
  
  hits_df <- fdr_df %>% 
    filter(!is.na(F_statistic)) %>% 
    group_by(nObsRound) %>% 
    mutate(min_rank_true = min(rank[dataset == "true"])) %>% 
    filter(rank >= min_rank_true) %>% 
    mutate(max_rank_fdr = min(rank[FDR > alpha], na.rm = TRUE)) %>% 
    filter(rank < max_rank_fdr) %>% 
    filter(dataset == "true") %>% 
    ungroup() %>% 
    dplyr::select(representative, clustername, nObs, matches("qupm"),
                  nCoeffsH0, nCoeffsH1, rssH0, rssH1, 
                  df1, df2, F_statistic, pEC50H1, slopeH1, pEC50_slopeH1, 
                  detected_effectH1, nObsRound, pi, FDR)
  
  return(hits_df)
}
