#' Bootstrap null distribution of F statistics for FDR estimation
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
#' @param seed seed to set, default NULL corresponds to no seed set
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
#'   filter(clustername %in% paste0("protein", 1:20)) %>% 
#'   group_by(representative) %>% 
#'   mutate(nObs = n()) %>% 
#'   ungroup 
#' bootstrapNull(temp_df, B = 1)  
#' 
#' @export
#'
#' @importFrom stats lm
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
bootstrapNull <- function(df, maxit = 500,
                          optim_fun_h0 = min_RSS_h0,
                          optim_fun_h1 = min_RSS_h1,
                          optim_fun_h1_2 = NULL,
                          gr_fun_h0 = NULL,
                          gr_fun_h1 = NULL,
                          gr_fun_h1_2 = NULL,
                          seed = NULL,
                          ncores = 1,
                          B = 3){
  ec50_lower_limit = 
    min(unique(df$log_conc)[
      which(is.finite(unique(df$log_conc)))])
  
  ec50_upper_limit = 
    max(unique(df$log_conc)[
      which(is.finite(unique(df$log_conc)))])
  
  registerDoParallel(cores = ncores)
  if(!is.null(seed)){
    set.seed(seed, kind = "L'Ecuyer-CMRG")
  }
  unique_names <- unique(df$clustername)
  null_list <- foreach(prot = unique_names) %dopar% {
    df_prot <- filter(df, clustername == prot)
    prot_h0 <- lm(log2_value ~ 1 + as.factor(temperature),
                  data = df_prot)
    len_res <- length(residuals(prot_h0))
    
    out_list <- lapply(seq_len(B), function(boot){
      df_prot <- df_prot %>%
        mutate(log2_value = log2_value - residuals(prot_h0) +
                 sample(residuals(prot_h0), size = len_res, replace = TRUE))
      
      sum_df <- fitAndEvalDataset(df_prot, 
                                  optim_fun_h0 = optim_fun_h0,
                                  optim_fun_h1 = optim_fun_h1,
                                  optim_fun_h1_2 = optim_fun_h1_2,
                                  gr_fun_h0 = gr_fun_h0,
                                  gr_fun_h1 = gr_fun_h1,
                                  gr_fun_h1_2 = gr_fun_h1_2)
      
      return(sum_df)
    })
  }
  null_df <- bind_rows(lapply(null_list, function(x){
    bind_rows(lapply(seq_len(length(x)), function(i){
      x[[i]] %>%
        mutate(dataset = paste("bootstrap", 
                               as.character(i), sep = "_"))
      }))
  }))
  
  return(null_df)
}