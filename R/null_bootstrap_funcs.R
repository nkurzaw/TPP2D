#' Bootstrap null distribution of F statistics for FDR estimation
#'
#' @param df tidy data_frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param fcThres numeric value of minimal fold change 
#' (or inverse fold change) a protein has to show to be kept 
#' upon independent filtering
#' @param minObs numeric value of minimal number of observations
#' that should be required per protein
#' @param independentFiltering boolean flag indicating whether
#' independent filtering should be performed based on minimal
#' fold changes per protein profile
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
#' @param ncores numeric value of numbers of cores that the function 
#' should use to parallelize
#' @param B numeric value of rounds of bootstrap, default: 3
#' 
#' @return data frame containing F statistics of proteins with 
#' permuted 2D thermal profiles that are informative on the Null
#' distribution of F statistics
#' 
#' @examples
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>% 
#'   filter(clustername %in% paste0("protein", 1:3)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' boot_df <- bootstrapNull(temp_df, B = 2)  
#' 
#' @export
#'
#' @importFrom stats lm
#' @importFrom stats residuals
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @import dplyr
bootstrapNull <- function(df, maxit = 500,
                          independentFiltering = FALSE,
                          fcThres = 1.5, minObs = 20,
                          optim_fun_h0 = min_RSS_h0,
                          optim_fun_h1 = min_RSS_h1_slope_pEC50,
                          optim_fun_h1_2 = NULL,
                          gr_fun_h0 = NULL,
                          gr_fun_h1 = NULL,
                          gr_fun_h1_2 = NULL,
                          ncores = 1,
                          B = 3){
  
  clustername <- prot <- log2_value <- NULL
  
  checkDfColumns(df)
  
  if(identical(optim_fun_h1, min_RSS_h1_slope_pEC50)){
    slopEC50 = TRUE
  }else{
    slopEC50 = FALSE
  }
  
  ec50_limits <- getEC50Limits(df)
  
  df_fil <- minObsFilter(df, minObs = minObs)
  
  if(independentFiltering){
    message("Independent Filtering: removing proteins without 
            any values crossing the threshold.")
    df_fil <- independentFilter(df_fil, fcThres = fcThres) 
  }
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  unique_names <- unique(df_fil$clustername)
  null_list <- foreach(prot = unique_names) %dopar% {
    df_prot <- filter(df_fil, clustername == prot)
    prot_h0 <- lm(log2_value ~ 1 + as.factor(temperature),
                  data = df_prot)
    len_res <- length(residuals(prot_h0))
    
    out_list <- lapply(seq_len(B), function(boot){
      df_resample_prot <- df_prot %>%
        mutate(log2_value = log2_value - residuals(prot_h0) +
                 sample(residuals(prot_h0), size = len_res, replace = TRUE))
      
      sum_df <- fitAndEvalDataset(df_resample_prot, 
                                  optim_fun_h0 = optim_fun_h0,
                                  optim_fun_h1 = optim_fun_h1,
                                  optim_fun_h1_2 = optim_fun_h1_2,
                                  gr_fun_h0 = gr_fun_h0,
                                  gr_fun_h1 = gr_fun_h1,
                                  gr_fun_h1_2 = gr_fun_h1_2,
                                  ec50_lower_limit = ec50_limits[1],
                                  ec50_upper_limit = ec50_limits[2],
                                  slopEC50 = slopEC50)
      
      return(sum_df)
    })
  }
  stopCluster(cl)
  null_df <- bind_rows(lapply(null_list, function(x){
    bind_rows(lapply(seq_len(length(x)), function(i){
      x[[i]] %>%
        mutate(dataset = paste("bootstrap", 
                               as.character(i), sep = "_"))
      }))
  }))
  
  return(null_df)
}