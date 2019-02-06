#' Fit H0 model and evaluate fit statistics
#' 
#' @param df tidy data_frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param maxit maximal number of iterations the optimization
#' should be given, default is set to 500
#' @param optim_fun optimization function that should be used
#' for fitting the H0 model
#' @param gr_fun optional gradient function for optim_fun,
#' default is NULL
#' 
#' @return data frame with H0 model characteristics for each
#' protein
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:10)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#'   
#' fitH0Model(temp_df)
#' 
#' @import dplyr
#' @importFrom stats optim
fitH0Model <- function(df, 
                       maxit = 500,
                       optim_fun = min_RSS_h0,
                       gr_fun = NULL){
  
  representative <- clustername <- nObs <- 
    temperature <- . <- NULL
  
  h0_df <- df %>%
    group_by(representative, clustername, nObs) %>%
    mutate(temp_i = dense_rank(temperature)) %>%
    do({
      unique_temp <- unique(.$temperature)
      len_temp <- length(unique_temp)
      start_par = vapply(unique_temp, function(x)
        mean(filter(., temperature == x)$log2_value), 1)
      h0_model = try(optim(par = start_par,
                       fn = optim_fun,
                       len_temp = len_temp,
                       data = .,
                       method = "L-BFGS-B",
                       gr = gr_fun,
                       control = list(maxit = maxit)))
      eval_optim_result(h0_model, hypothesis = "H0",
                        data = .)
    }) %>%
    group_by(representative, clustername) %>%
    ungroup()
  
  return(h0_df)
}

#' Fit H1 model and evaluate fit statistics
#' 
#' @param df tidy data_frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param maxit maximal number of iterations the optimization
#' should be given, default is set to 500
#' @param optim_fun optimization function that should be used
#' for fitting the H0 model
#' @param optim_fun_2 optional secound optimization function for
#' fitting the H1 model that should be used based on the fitted 
#' parameters of the optimizationfor based on optim_fun
#' @param gr_fun optional gradient function for optim_fun,
#' default is NULL
#' @param gr_fun_2 optional gradient function for optim_fun_2,
#' default is NULL
#' @param ec50_lower_limit lower limit of ec50 parameter
#' @param ec50_upper_limit lower limit of ec50 parameter
#' @param slopEC50 logical flag indicating whether the h1 model is
#' fitted with a linear model describing the shift od the pEC50 over 
#' temperatures
#' 
#' @return data frame with H1 model characteristics for each
#' protein
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup
#'   
#' fitH1Model(temp_df)
#' 
#' fitH1Model(temp_df,
#'            optim_fun_2 = 
#'              TPP2D:::min_RSS_h1_trim)
#' 
#' @import dplyr
#' @importFrom stats optim
#' @importFrom methods is
fitH1Model <- function(df, 
                       maxit = 500,
                       optim_fun = min_RSS_h1,
                       optim_fun_2 = NULL,
                       gr_fun = NULL,
                       gr_fun_2 = NULL,
                       ec50_lower_limit = NULL,
                       ec50_upper_limit = NULL,
                       slopEC50 = FALSE){
  
  representative <- clustername <- nObs <- 
    temperature <- . <- NULL
  
  if(is.null(ec50_lower_limit)){
    ec50_lower_limit <- min(unique(df$log_conc)[
      which(is.finite(unique(df$log_conc)))])
  }
  if(is.null(ec50_upper_limit)){
    ec50_upper_limit <- max(unique(df$log_conc)[
      which(is.finite(unique(df$log_conc)))])
  }
  
  h1_df <- df %>%
    group_by(representative, clustername, nObs) %>%
    mutate(temp_i = dense_rank(temperature)) %>%
    do({
      unique_temp <- unique(.$temperature)
      len_temp <- length(unique_temp)
      start_par = getStartParameters(
        df = ., unique_temp = unique_temp, 
        len_temp = len_temp, 
        slopEC50 = slopEC50)
      opt_limits <- getOptLimits(
        ec50Limits = c(ec50_lower_limit,
                       ec50_upper_limit),
        len_temp = len_temp, 
        slopEC50 = slopEC50)
      lower = opt_limits$lower 
      upper = opt_limits$upper
      h1_model = try(optim(par = start_par,
                       fn = optim_fun,
                       len_temp = len_temp,
                       data = .,
                       method = "L-BFGS-B",
                       upper = upper,
                       lower = lower,
                       gr = gr_fun,
                       control = list(maxit = maxit)))
      if(!is.null(optim_fun_2) & 
         !is(h1_model, "try-error")){
        h1_model = try(optim(par = h1_model$par,
                         fn = optim_fun_2,
                         len_temp = len_temp,
                         data = .,
                         method = "L-BFGS-B",
                         upper = upper,
                         lower = lower,
                         gr = gr_fun_2,
                         control = list(maxit = maxit)))
      }
        eval_optim_result(h1_model, hypothesis = "H1",
                          data = ., len_temp = len_temp)
    }) %>%
    group_by(representative, clustername) %>%
    ungroup()
  
  return(h1_df)
}

#' @importFrom methods is
eval_optim_result <- function(optim_result, hypothesis = "H1",
                              data, len_temp = NULL,
                              slopEC50 = FALSE){
  # evaluate optimization results for H0 or H1 models 
  
  if(!is(optim_result, "try-error")){
    if(hypothesis == "H1"){
      
      pEC50 = -optim_result$par[1]
      if(!slopEC50){
        slope = optim_result$par[2]
      }else{
        pEC50_slope = optim_result$par[2]
        slope = optim_result$par[3]
      }
      rss = optim_result$value
      nCoeffs = length(optim_result$par)
      fitStats <- data.frame(rss = rss,
                             nCoeffs = nCoeffs,
                             pEC50 = pEC50,
                             slope = slope)
      if(slopEC50){
        fitStats$pEC50_slope <- pEC50_slope
      }
      
      if(!is.null(len_temp)){
        if(!slopEC50){
          alpha <- optim_result$par[(4 + len_temp):(3 + len_temp*2)]
        }else{
          alpha <- optim_result$par[(5 + len_temp):(4 + len_temp*2)]
        }
        
        if(alpha[1] > (max(alpha[-1])/3)){
          fitStats$detected_effect <- 
            "expression/solubility"
        }else{
          fitStats$detected_effect <- 
            "stability"
        }
      }
      names(fitStats) <- paste0(names(fitStats), hypothesis)
      return(fitStats)
      
    }else{
      rss = optim_result$value
      nCoeffs = length(optim_result$par)
      fitStats <- data.frame(rss = rss,
                             nCoeffs = nCoeffs)
      names(fitStats) <- paste0(names(fitStats), hypothesis)
      return(fitStats)
    }
  }else{
    fitStats <- data.frame(rss = NA,
                           nCoeffs = NA,
                           pEC50 = NA,
                           slope = NA)
    names(fitStats) <- paste0(names(fitStats), hypothesis)
    return(fitStats)
  }
}

#' Compute F statistic from H1 and H0 model characteristics
#'
#' @param h0_df data frame with H0 model characteristics for each
#' protein
#' @param h1_df data frame with H1 model characteristics for each
#' protein
#' 
#' @return data frame with H0 and H1 model characteristics for each
#' protein and respectively computed F statistics
#'
#'
#' @examples
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:20)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#'   
#' h0_df <- fitH0Model(temp_df)
#' h1_df <- fitH1Model(temp_df)
#'   
#' computeFstat(h0_df, h1_df)
#' 
#' @export
#' 
#' @import dplyr
computeFstat <- function(h0_df, h1_df){
  
  nCoeffsH1 <- nCoeffsH0 <- nObs <- rssH0 <- 
    rssH1 <- df2 <- df1 <- NULL
  
  sum_df <- 
    left_join(h0_df, h1_df, by = c("representative",
                                   "clustername", "nObs")) %>%
    ungroup() %>%
    mutate(df1 = nCoeffsH1 - nCoeffsH0, df2 = nObs - nCoeffsH1) %>%
    mutate(F_statistic = ((rssH0 - rssH1) / rssH1) * (df2/df1))
  
  return(sum_df)
}

#' Fit H0 and H1 model to 2D thermal profiles of proteins
#' and compute F statistic
#'
#' @param df tidy data_frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param maxit maximal number of iterations the optimization
#' should be given, default is set to 500
#' @param optim_fun_h0 optimization function that should be used
#' for fitting the H0 model
#' @param optim_fun_h1 optimization function that should be used
#' for fitting the H1 model
#' @param optim_fun_h1_2 optional additional optimization function 
#' that will be run with paramters retrieved from optim_fun_h1 and 
#' should be used for fitting the H1 model with the trimmed sum
#' model, default is NULL
#' @param gr_fun_h0 optional gradient function for optim_fun_h0,
#' default is NULL
#' @param gr_fun_h1 optional gradient function for optim_fun_h1,
#' default is NULL
#' @param gr_fun_h1_2 optional gradient function for optim_fun_h1_2,
#' default is NULL
#' @param ec50_lower_limit lower limit of ec50 parameter
#' @param ec50_upper_limit lower limit of ec50 parameter
#' @param slopEC50 logical flag indicating whether the h1 model is
#' fitted with a linear model describing the shift od the pEC50 over 
#' temperatures
#' 
#' @return data frame with H0 and H1 model characteristics for each
#' protein and respectively computed F statistics
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' fitAndEvalDataset(temp_df)  
#' 
#' @export
fitAndEvalDataset <- function(df, maxit = 500,
                              optim_fun_h0 = min_RSS_h0,
                              optim_fun_h1 = min_RSS_h1,
                              optim_fun_h1_2 = NULL,
                              gr_fun_h0 = NULL,
                              gr_fun_h1 = NULL,
                              gr_fun_h1_2 = NULL,
                              ec50_lower_limit = NULL,
                              ec50_upper_limit = NULL,
                              slopEC50 = FALSE){
  
  h0_df <- fitH0Model(df = df,
                      maxit = maxit,
                      optim_fun = optim_fun_h0,
                      gr_fun = gr_fun_h0)
  
  h1_df <- fitH1Model(df = df,
                      maxit = maxit,
                      optim_fun = optim_fun_h1,
                      optim_fun_2 = optim_fun_h1_2,
                      gr_fun = gr_fun_h1,
                      gr_fun_2 = gr_fun_h1_2,
                      ec50_lower_limit = ec50_lower_limit,
                      ec50_upper_limit = ec50_upper_limit,
                      slopEC50 = slopEC50)
  
  sum_df <- computeFstat(h0_df, h1_df)
  
  return(sum_df)
}


minObsFilter <- function(df, minObs = 20){
  # Filter data frame for a minimal number of observations
  # per protein
  representative <- clustername <- rel_value <- 
    nObs <- NULL
  
  df_fil <- df %>%
    group_by(representative, clustername) %>%
    mutate(nObs = n()) %>%
    filter(nObs >= minObs) %>%
    ungroup()
  
  return(df_fil)
}

independentFilter <- function(df, fcThres = 1.5){
  # Filter data frame independently based on maximal 
  # fold change per protein
  representative <- clustername <- rel_value <- NULL
  
  df_fil <- df %>%
    group_by(representative, clustername) %>%
    filter(any(rel_value > fcThres) | any(rel_value < 1/fcThres)) %>%
    ungroup
  
  return(df_fil)
}


getEC50Limits <- function(df){
  log_conc <- NULL
  
  ec50_lower_limit <- min(unique(df$log_conc)[
    which(is.finite(unique(df$log_conc)))])
  ec50_upper_limit <- max(unique(df$log_conc)[
    which(is.finite(unique(df$log_conc)))])
  
  return(c(ec50_lower_limit, 
           ec50_upper_limit))
}

#' Compete H0 and H1 models per protein and obtain F statistic
#' 
#' @param df tidy data frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param fcThres numeric value of minimal fold change 
#' (or inverse fold change) a protein has to show to be kept 
#' upon independent filtering
#' @param independentFiltering boolean flag indicating whether
#' independent filtering should be performed based on minimal
#' fold changes per protein profile
#' @param minObs numeric value of minimal number of observations
#' that should be required per protein
#' @param maxit maximal number of iterations the optimization
#' should be given, default is set to 500
#' @param optim_fun_h0 optimization function that should be used
#' for fitting the H0 model
#' @param optim_fun_h1 optimization function that should be used
#' for fitting the H1 model
#' @param optim_fun_h1_2 optional additional optimization function 
#' that will be run with paramters retrieved from optim_fun_h1 and 
#' should be used for fitting the H1 model with the trimmed sum
#' model, default is NULL
#' @param gr_fun_h0 optional gradient function for optim_fun_h0,
#' default is NULL
#' @param gr_fun_h1 optional gradient function for optim_fun_h1,
#' default is NULL
#' @param gr_fun_h1_2 optional gradient function for optim_fun_h1_2,
#' default is NULL
#' @param slopEC50 logical flag indicating whether the h1 model is
#' fitted with a linear model describing the shift od the pEC50 over 
#' temperatures
#' 
#' @return data frame summarising the fit characteristics of H0 and
#' H1 models and therof resulting computed F statistics per protein
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:10)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' competeModels(temp_df)  
#' 
#' @export
competeModels <- function(df, fcThres = 1.5,
                          independentFiltering = FALSE,
                          minObs = 20,
                          optim_fun_h0 = min_RSS_h0_trim,
                          optim_fun_h1 = min_RSS_h1,
                          optim_fun_h1_2 = min_RSS_h1_trim,
                          gr_fun_h0 = NULL,
                          gr_fun_h1 = NULL,
                          gr_fun_h1_2 = NULL,
                          maxit = 750,
                          slopEC50 = FALSE){
  
  
  ec50_limits <- getEC50Limits(df)
  
  df_fil <- minObsFilter(df, minObs = minObs)
  
  if(independentFiltering){
    message(paste("Independent Filtering: removing proteins without", 
            "any values crossing the set threshold."))
    df_fil <- independentFilter(df_fil, fcThres = fcThres) 
  }

  sum_df <- fitAndEvalDataset(
    df_fil,
    maxit = maxit,
    optim_fun_h0 = optim_fun_h0,
    optim_fun_h1 = optim_fun_h1,
    optim_fun_h1_2 = optim_fun_h1_2,
    ec50_lower_limit = ec50_limits[1],
    ec50_upper_limit = ec50_limits[2],
    slopEC50 = slopEC50)
  
  return(sum_df)
}