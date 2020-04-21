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
#'   filter(clustername %in% paste0("protein", 1:5)) %>% 
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
                       optim_fun = .min_RSS_h0,
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
      .eval_optim_result(h0_model, hypothesis = "H0",
                        data = .)
    }) %>%
    group_by(representative, clustername) %>%
    ungroup()
  
  return(h0_df)
}

.fitEvalH1 <- function(df_fil, unique_temp, len_temp,
                      optim_fun = .min_RSS_h1_slopeEC50, 
                      optim_fun_2 = NULL,
                      gr_fun = NULL,
                      gr_fun_2 = NULL,
                      ec50_lower_limit = NULL,
                      ec50_upper_limit= NULL,
                      slopEC50 = TRUE, maxit = 500){
  .min_RSS_h1_slopeEC50 <- NULL
  if(is.null(ec50_lower_limit)){
    ec50_lower_limit <- min(unique(df_fil$log_conc)[
      which(is.finite(unique(df_fil$log_conc)))])
  }
  if(is.null(ec50_upper_limit)){
    ec50_upper_limit <- max(unique(df_fil$log_conc)[
      which(is.finite(unique(df_fil$log_conc)))])
  }
  start_par = .getStartParameters(
    df = df_fil, 
    unique_temp = unique_temp, 
    len_temp = len_temp, 
    slopEC50 = slopEC50)
  opt_limits <- .getOptLimits(
    ec50Limits = c(ec50_lower_limit, ec50_upper_limit),
    len_temp = len_temp, 
    slopEC50 = slopEC50)
  lower = opt_limits$lower 
  upper = opt_limits$upper
  h1_model = try(optim(par = start_par,
                       fn = optim_fun,
                       gr = gr_fun,
                       len_temp = len_temp,
                       data = df_fil,
                       method = "L-BFGS-B",
                       upper = upper,
                       lower = lower,
                       control = list(maxit = maxit)))
  if(!is.null(optim_fun_2) & 
     !is(h1_model, "try-error")){
    h1_model = try(optim(par = h1_model$par,
                         fn = optim_fun_2,
                         gr = gr_fun_2,
                         len_temp = len_temp,
                         data = df_fil,
                         method = "L-BFGS-B",
                         upper = upper,
                         lower = lower,
                         control = list(maxit = maxit)))
  }
  return(h1_model)
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
#'   filter(clustername %in% paste0("protein", 1:5)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup
#'   
#' fitH1Model(temp_df)
#' 
#' @import dplyr
#' @importFrom stats optim
#' @importFrom methods is
fitH1Model <- function(df, 
                       maxit = 500,
                       optim_fun = .min_RSS_h1_slope_pEC50,
                       optim_fun_2 = NULL,
                       gr_fun = NULL,
                       gr_fun_2 = NULL,
                       ec50_lower_limit = NULL,
                       ec50_upper_limit = NULL,
                       slopEC50 = TRUE){
  
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
      h1_model <- .fitEvalH1(df_fil = .,
                            unique_temp = unique_temp, 
                            len_temp = len_temp,
                            optim_fun = optim_fun, 
                            optim_fun_2 = optim_fun_2,
                            gr_fun = gr_fun,
                            gr_fun_2 = gr_fun_2,
                            slopEC50 = slopEC50, 
                            maxit = maxit)
      .eval_optim_result(h1_model, hypothesis = "H1",
                        data = ., len_temp = len_temp,
                        slopEC50 = slopEC50)
    }) %>%
    group_by(representative, clustername) %>%
    ungroup()
  
  return(h1_df)
}

#' @importFrom methods is
#' @importFrom stats fft
#' @importFrom MASS lda
.eval_optim_result <- function(optim_result, hypothesis = "H1",
                              data, len_temp = NULL,
                              slopEC50 = TRUE){
  # evaluate optimization results for H0 or H1 models 
  spec <- NULL
  if(!is(optim_result, "try-error")){
    if(hypothesis == "H1"){
      pEC50 <-  -optim_result$par[1]
      if(!slopEC50){
        slope <-  optim_result$par[2]
      }else{
        pEC50_slope <- optim_result$par[2]
        slope <- optim_result$par[3]
      }
      rss <- optim_result$value
      nCoeffs <- length(optim_result$par)
      fitStats <- 
        data.frame(rss = rss, nCoeffs = nCoeffs,
                   pEC50 = pEC50, slope = slope)
      if(slopEC50){
        fitStats$pEC50_slope <- pEC50_slope
      }
      if(!is.null(len_temp)){
        if(!slopEC50){
            alpha <- optim_result$par[(4 + len_temp):(3 + len_temp*2)]
        }else{
            alpha <- optim_result$par[(5 + len_temp):(4 + len_temp*2)]
        }
        # alpha_fft_df <- data.frame(freq = seq_len(length(alpha)),
        #                            strength = Mod(fft(alpha)))
        # alpha_fft_fit <- lm(strength ~ poly(freq, 2), 
        #                     data = alpha_fft_df[-1,])
        # freq <- Mod(fft(c(alpha, rep(0, 12-length(alpha)))))
        # alpha_fft_df <- tibble(freq) %>% 
        #     mutate(spec = paste0("T", 1:n())) %>% 
        #     arrange(spec) %>% 
        #     spread(spec, freq)
        # lda_pred <- predict(lda_model, alpha_fft_df)
        # if(lda_pred$class == "yes"){
        #     fitStats$detected_effect <- "carry-over"
        # }else 
        if(alpha[1] > (max(alpha[-1])/3)){
            fitStats$detected_effect <- "expression/solubility"
        }else{
            fitStats$detected_effect <- "stability"
        }
      }
      names(fitStats) <- paste0(names(fitStats), hypothesis)
      return(fitStats)
    }else{
      rss = optim_result$value
      nCoeffs = length(optim_result$par)
      fitStats <- data.frame(rss = rss,  nCoeffs = nCoeffs)
      names(fitStats) <- paste0(names(fitStats), hypothesis)
      return(fitStats)
    }
  }else{
    fitStats <- 
      data.frame(rss = NA,nCoeffs = NA,
                 pEC50 = NA, slope = NA)
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
                              optim_fun_h0 = .min_RSS_h0,
                              optim_fun_h1 = .min_RSS_h1_slope_pEC50,
                              optim_fun_h1_2 = NULL,
                              gr_fun_h0 = NULL,
                              gr_fun_h1 = NULL,
                              gr_fun_h1_2 = NULL,
                              ec50_lower_limit = NULL,
                              ec50_upper_limit = NULL,
                              slopEC50 = TRUE){
  
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


.minObsFilter <- function(df, minObs = 20){
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

.independentFilter <- function(df, fcThres = 1.5){
  # Filter data frame independently based on maximal 
  # fold change per protein
  representative <- clustername <- rel_value <- NULL
  
  df_fil <- df %>%
    group_by(representative, clustername) %>%
    filter(any(rel_value > fcThres) | any(rel_value < 1/fcThres)) %>%
    ungroup
  
  return(df_fil)
}


.getEC50Limits <- function(df){
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
                          optim_fun_h0 = .min_RSS_h0,
                          optim_fun_h1 = .min_RSS_h1_slope_pEC50,
                          optim_fun_h1_2 = NULL,
                          gr_fun_h0 = NULL,
                          gr_fun_h1 = NULL,
                          gr_fun_h1_2 = NULL,
                          maxit = 750){
  
  .checkDfColumns(df)
  
  if(identical(optim_fun_h1, .min_RSS_h1_slope_pEC50)){
    slopEC50 = TRUE
  }else{
    slopEC50 = FALSE
  }
  
  ec50_limits <- .getEC50Limits(df)
  
  df_fil <- .minObsFilter(df, minObs = minObs)
  
  if(independentFiltering){
    message(paste("Independent Filtering: removing proteins without", 
            "any values crossing the set threshold."))
    df_fil <- .independentFilter(df_fil, fcThres = fcThres) 
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

#' Compute F statistics from paramter data frame
#' 
#' @param params_df data frame listing all null and alternative
#' model parameters as obtained by 'getModelParamsDf'
#' @param in_df data frame of 2D-TPP profiles obtaine after 
#' data import
#' @param shrinkTo2ndHighestTemperatureFstat logical indicating whether
#' F statistic should be shrunked towards the 2nd highest per Temperature
#' F statistic
#' 
#' @return data frame of all proteins and computed F statistics
#' and parameters that were used for the computation
#' 
#' @examples
#' data("simulated_cell_extract_df")
#' params_df <- getModelParamsDf(simulated_cell_extract_df)
#' computeFStatFromParams(params_df)
#'  
#' @export
#' 
#' @import dplyr
computeFStatFromParams <- function(params_df, in_df = NULL, 
                                   shrinkTo2ndHighestTemperatureFstat = FALSE){
    nCoeffsH0 <- nCoeffsH1 <- nObs <- rssH0 <- rssH1 <- df1 <- 
        df2 <- representative <- clustername <- min_qupm <- 
        max_qupm <- pEC50H1 <- slopeH1 <- pEC50_slopeH1 <- 
        detected_effectH1 <- F_statistic <- NULL
    
    fstat_df <- params_df %>%
        mutate(df1 = nCoeffsH1 - nCoeffsH0, df2 = nObs - nCoeffsH1) %>%
        mutate(F_statistic = ((rssH0 - rssH1) / rssH1) * (df2/df1)) %>% 
        dplyr::select(representative, clustername, nObs, 
                      min_qupm, max_qupm,
                      nCoeffsH0, nCoeffsH1, rssH0, rssH1,
                      pEC50H1, slopeH1, pEC50_slopeH1, 
                      detected_effectH1,
                      df1, df2, F_statistic) 
    
    if(shrinkTo2ndHighestTemperatureFstat){
        if(is.null(in_df)){
            stop(paste("Please supply the data frame obtained",
                       "after data import as 'in_df' is using",
                       "the option 'shrinkFTo2ndTemperature = TRUE'!"))
        }
        fstat_2nd_temp_df <- bind_rows(
            lapply(params_df$clustername, 
                   .get2ndHighestTemperatureFstat,
                   #.getPerTemperatureFstats,
                   in_df = in_df, params_df = params_df)
        )
        
        fstat_df <- left_join(
            fstat_df %>% rename(F_statistic_orig = F_statistic), 
            fstat_2nd_temp_df,
            by = c("representative", "clustername")) %>% 
            #mutate(F_statistic = F_statistic_orig + F_statisticT)
            rowwise() %>% 
            mutate(F_statistic = mean(
                c(F_statisticT, F_statistic_orig))) %>% 
            ungroup
    }
    
    return(fstat_df)
}

#' Get pEC50 for a protein of interest at a specific temperatures
#' (optimally the melting point of the protein)
#' 
#' @param fstat_df data frame as obtained after calling 
#' \code{getModelParamsDf}, containing fitted null and 
#' alternative model parameters for each protein
#' @param protein character string referring to the protein 
#' of interest
#' @param temperaturePEC50 temperature (numeric) at which pEC50
#' should be inferred
#' 
#' @return numeric value specifying the pEC50 for the 
#' indicated protein and temperature
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' 
#' model_params_df <- getModelParamsDf(
#'    df = filter(simulated_cell_extract_df, 
#'            clustername == "tp1"))
#' 
#' getPEC504Temperature(
#'     fstat_df = model_params_df, 
#'     protein = "tp1", 
#'     temperaturePEC50 = 60)
#' @import dplyr
#' @export
getPEC504Temperature <- function(fstat_df, protein, 
                                 temperaturePEC50 = 60){
    clustername <- NULL
    fstat_fil <- filter(fstat_df, clustername == protein)
    pEC50 <- fstat_fil$pEC50H1 - 
        (temperaturePEC50 *fstat_fil$pEC50_slopeH1)
    return(pEC50)
}


.get2ndHighestTemperatureFstat <- function(in_df, params_df, gene_name){
    filtered_params_df <- filter(params_df, clustername == gene_name)
    temp_df <- filter(in_df, clustername == gene_name) %>% 
        arrange(temperature) %>% 
        mutate(residualsH0 = filtered_params_df$residualsH0[[1]], 
               residualsH1 = filtered_params_df$residualsH1[[1]],
               nObs = n(),
               min_qupm = min(qupm),
               max_qupm = max(qupm),
               nCoeffsH0 = length(filtered_params_df$parH0[[1]]),
               nCoeffsH1 = length(filtered_params_df$parH1[[1]]),
               pEC50H1 = filtered_params_df$pEC50H1[1],
               slopeH1 = filtered_params_df$slopeH1[1],
               pEC50_slopeH1 = filtered_params_df$pEC50_slopeH1[1]) %>% 
        group_by(temperature) %>% 
        mutate(rssH0T = sum(residualsH0^2), 
               rssH1T = sum(residualsH1^2),
               nObsTemperature = n()) %>% 
        mutate(df1T = (nObsTemperature/nObs) * nCoeffsH1 - 
                   (nObsTemperature/nObs) * nCoeffsH0,
               df2T = nObsTemperature - 
                   (nObsTemperature/nObs) * nCoeffsH1) %>% 
        ungroup %>% 
        mutate(F_statisticT = ((rssH0T - rssH1T)/rssH1T) * (df2T/df1T))
    fselect <- sort(unique(temp_df$F_statisticT), decreasing = TRUE)[2]
    out_df <- filter(temp_df, F_statisticT == fselect) %>% 
        filter(!duplicated(clustername)) %>% 
        dplyr::select(representative, clustername, #rssH0T, rssH1T, df1T, df2T, 
                      F_statisticT)
    return(out_df)
}

.getPerTemperatureFstats <- function(in_df, params_df, gene_name){
    filtered_params_df <- filter(params_df, clustername == gene_name)
    temp_df <- filter(in_df, clustername == gene_name) %>% 
        arrange(temperature) %>% 
        mutate(residualsH0 = filtered_params_df$residualsH0[[1]], 
               residualsH1 = filtered_params_df$residualsH1[[1]],
               nObs = n(),
               min_qupm = min(qupm),
               max_qupm = max(qupm),
               nCoeffsH0 = length(filtered_params_df$parH0[[1]]),
               nCoeffsH1 = length(filtered_params_df$parH1[[1]]),
               pEC50H1 = filtered_params_df$pEC50H1[1],
               slopeH1 = filtered_params_df$slopeH1[1],
               pEC50_slopeH1 = filtered_params_df$pEC50_slopeH1[1]) %>% 
        group_by(temperature) %>% 
        mutate(rssH0T = sum(residualsH0^2), 
               rssH1T = sum(residualsH1^2),
               nObsTemperature = n()) %>% 
        mutate(df1T = (nObsTemperature/nObs) * nCoeffsH1 - 
                   (nObsTemperature/nObs) * nCoeffsH0,
               df2T = nObsTemperature - 
                   (nObsTemperature/nObs) * nCoeffsH1) %>% 
        ungroup %>% 
        mutate(F_statisticT = ((rssH0T - rssH1T)/rssH1T) * (df2T/df1T)) %>% 
        filter(!duplicated(temperature))
    fselect <- sort(unique(temp_df$F_statisticT), decreasing = TRUE)[2]
    out_df <- filter(temp_df, !duplicated(clustername)) %>% 
        mutate(F_statisticT_list = 
                   list(rep(temp_df$F_statisticT, temp_df$nObsTemperature))) %>% 
        dplyr::select(representative, clustername, F_statisticT_list)
    return(out_df)
}
