.fitH0ParamDf <- function(df, 
                          optim_fun = .min_RSS_h0,
                          gr_fun = NULL,
                          slopEC50 = TRUE, 
                          maxit = 500){
    h0_param_df <- df %>% 
        group_by(representative, clustername, nObs) %>%
        mutate(temp_i = dense_rank(temperature)) %>%
        arrange(clustername, temp_i, log_conc) %>% 
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
                                 control = list(maxit = maxit)))
            fit_df <- .getFitDf(
                df_fil = ., 
                conc_vec = unique(.$log_conc),
                model_type = "H0",
                optim_model = h0_model,
                slopEC50 = slopEC50)
            
            fit_est_df <- left_join(
                ., fit_df,
                by = c("log_conc", "temp_i")) %>%
                mutate(residuals = log2_value - y_hat)
            
            tibble(min_qupm = min(df_fil$qupm),
                   max_qupm = max(df_fil$qupm),
                   nCoeffs = length(h0_model$par),
                   rss = h0_model$value,
                   par = list(h0_model$par),
                   estimate = list(fit_est_df$y_hat),
                   residuals = list(fit_est_df$residuals))

        }) %>% 
        ungroup()
    return(h0_param_df)
}

.fitH1ParamDf <- function(df, optim_fun = .min_RSS_h1_slope_pEC50, 
                         optim_fun_2 = NULL,
                         gr_fun = NULL,
                         gr_fun_2 = NULL,
                         slopEC50 = TRUE, 
                         maxit = 500){
    h1_param_df <- df %>% 
        group_by(representative, clustername, nObs) %>%
        mutate(temp_i = dense_rank(temperature)) %>%
        arrange(clustername, temp_i, log_conc) %>% 
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
            fit_df <- .getFitDf(
                df_fil = ., 
                conc_vec = unique(.$log_conc),
                model_type = "H1",
                optim_model = h1_model,
                slopEC50 = TRUE) 
            
            fit_est_df <- left_join(
                ., fit_df,
                by = c("log_conc", "temp_i")) %>%
                mutate(residuals = log2_value - y_hat)
            
            tibble(nCoeffs = length(h1_model$par),
                   rss = h1_model$value,
                   par = list(h1_model$par),
                   estimate = list(fit_est_df$y_hat),
                   residuals = list(fit_est_df$residuals))
            
        }) %>% 
        ungroup()
    return(h1_param_df)
}

#' Get H0 and H1 model parameters
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
#' @param slopEC50 logical flag indicating whether the h1 model is
#' fitted with a linear model describing the shift od the pEC50 over 
#' temperatures
#' 
#' @return a data.frame with fitted null and alternative model 
#' parameters
#' 
#' @export
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' getModelParamsDf(.minObsFilter(simulated_cell_extract_df))
getModelParamsDf <- function(df,
                             optim_fun_h0 = .min_RSS_h0,
                             optim_fun_h1 = .min_RSS_h1_slope_pEC50,
                             optim_fun_h1_2 = NULL,
                             gr_fun_h0 = NULL,
                             gr_fun_h1 = NULL,
                             gr_fun_h1_2 = NULL,
                             slopEC50 = TRUE,
                             maxit = 750){
    
    h0_param_df <- .fitH0ParamDf(
        df,
        optim_fun = optim_fun_h0,
        gr_fun = gr_fun_h0,
        slopEC50 = slopEC50, 
        maxit = maxit)
    
    h1_param_df <- .fitH1ParamDf(
        df,
        optim_fun = optim_fun_h1, 
        optim_fun_2 = optim_fun_h1_2,
        gr_fun = gr_fun_h1,
        gr_fun_2 = gr_fun_h1_2,
        slopEC50 = slopEC50, 
        maxit = maxit)
    
    combo_param_df <- left_join(
        h0_param_df, 
        h1_param_df, 
        by = c("representative", "clustername", "nObs"), 
        suffix = c("H0", "H1"))
    
    return(combo_param_df)
}