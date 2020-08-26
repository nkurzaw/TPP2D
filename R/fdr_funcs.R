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
        dataset <- FDR <- all_true <- all_null <- 
        remaining_null <- remaining_true <- NULL
    
    if(squeezeDenominator){
        df_out <- .shrinkFstat(df_out, trueOrNull = "true")
        df_null <- .shrinkFstat(df_null, trueOrNull = "null")
    }
    
    B <- max(as.numeric(
        gsub("bootstrap_", "", unique(df_null$dataset))))
    
    out_df <- bind_rows(df_out %>% mutate(dataset = "true"),
                        df_null) %>%
        mutate(nObsRound = round(nObs, digits = -1)) %>%
        group_by(nObsRound) %>%
        arrange(desc(F_statistic)) %>%
        mutate(max_rank = n(),
               rank = dense_rank(desc(F_statistic)),
               is_decoy = ifelse(dataset != "true", 1, 0)) %>%
        mutate(all_null = sum(is_decoy),
               all_true = sum(!is_decoy),
               true_cumsum = cumsum(!is_decoy),
               null_cumsum = cumsum(is_decoy)) %>% 
        mutate(remaining_null = all_null - null_cumsum,
               remaining_true = all_true - true_cumsum) %>% 
        mutate(pi = remaining_true/(remaining_null/B)) %>% 
        mutate(FDR = pi * (null_cumsum/B)/true_cumsum) %>% 
        ungroup() %>% 
        mutate(FDR = ifelse(is.na(F_statistic), NA, FDR))
    
    return(out_df)
}

#' @importFrom limma squeezeVar 
.shrinkFstat <- function(inDf, trueOrNull = "true"){
  . <- rssH1 <- df2 <- rssH0 <- rssH1Squeezed <- 
    df1 <- df0 <- nObs <- nObsRound <- dataset <- NULL
  
  outDf <- inDf %>% 
    mutate(nObsRound = round(nObs, digits = -1)) %>% 
    { 
      if(trueOrNull == "true")
        group_by(., nObsRound)    
      else 
        group_by(., dataset, nObsRound) 
    } %>%
    do({
      squeezeResult <- limma::squeezeVar(.$rssH1, df = .$df2)
      mutate(., rssH1Squeezed = squeezeResult$var.post, 
             df0 = squeezeResult$df.prior)
    }) %>%
    mutate(df0 = ifelse(is.finite(df0), df0, 0)) %>% 
    ungroup %>% 
    mutate(F_statistic = (rssH0 - rssH1)/
             (rssH1Squeezed) * (df0 + df2)/df1)
  return(outDf)
}
#' Compute FDR for given F statistics based on true and
#' null dataset (old function)
#' @name computeFdr-defunct
#' @seealso \code{\link{TPP2D-defunct}}
#' @keywords internal
NULL

#' @rdname TPP2D-defunct
#' @section \code{computeFdr}:
#' For \code{computeFdr}, use \code{\link{getFDR}}.
computeFdr <- function(df_out, df_null){
    .Defunct("getFDR")
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
#' @param pseudo_count numeric larger or equal to 0 added to both
#' counts of protein with an F-statistic higher than a threshold
#' theta of the true and bootstrapped datasets
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
getPvalues <- function(df_out, df_null, pseudo_count = 1,
                       squeezeDenominator = FALSE){
  dataset <- nObs <- nObsRound <- F_statistic <- 
    is_decoy <- max_rank <- true_cumsum <- p_value <- 
    null_cumsum <- representative <- clustername <- 
    dataset <- FDR <- all_true <- all_null <- NULL
  
  stopifnot(pseudo_count >= 0)
  
  if(squeezeDenominator){
    df_out <- .shrinkFstat(df_out, trueOrNull = "true")
    df_null <- .shrinkFstat(df_null, trueOrNull = "null")
  }
  
  out_df <- bind_rows(df_out %>% mutate(dataset = "true"),
                      df_null) %>%
    mutate(nObsRound = round(nObs, digits = -1)) %>%
    group_by(nObsRound) %>%
    arrange(desc(F_statistic)) %>%
    mutate(max_rank = n(),
           rank = dense_rank(dplyr::desc(F_statistic)),
           is_decoy = ifelse(dataset != "true", 1, 0)) %>%
    mutate(all_null = sum(is_decoy),
           null_cumsum = cumsum(is_decoy)) %>% 
    mutate(p_value = (null_cumsum + pseudo_count)/
             (all_null + pseudo_count)) %>% 
    ungroup() %>% 
    filter(dataset == "true") %>% 
    mutate(p_adj = p.adjust(p_value, method = "BH"))
  
  return(out_df)
}

#' Find hits according to FDR threshold
#' 
#' @param fdr_df data frame obtained from computeFdr
#' @param alpha significance threshold, default is set to 0.1
#' 
#' @return data frame of significant hits at FDR <= alpha
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
    mutate(max_rank_fdr = min(c(Inf, rank[FDR > alpha]), na.rm = TRUE)) %>% 
    filter(rank < max_rank_fdr) %>% 
    filter(dataset == "true") %>% 
    ungroup() %>% 
    dplyr::select(representative, clustername, nObs, matches("qupm"),
                  nCoeffsH0, nCoeffsH1, rssH0, rssH1, 
                  df1, df2, F_statistic, pEC50H1, slopeH1, pEC50_slopeH1, 
                  detected_effectH1, nObsRound, pi, FDR)
  
  return(hits_df)
}
