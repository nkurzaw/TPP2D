getStartParameters <- function(df, unique_temp, len_temp){
  # function to generate vector of start parameters
  # for h1 model optimization
  log_conc <- log2_value <- NULL
  
  start_par <- c(mean(unique(df$log_conc)[which(is.finite(unique(df$log_conc)))]),
    coef(lm(log2_value ~ log_conc, filter(df, is.finite(log_conc))))[[2]],
    max(unlist(lapply(unique_temp, function(temp){
      max(filter(df, temperature == temp)$log2_value) -
        min(filter(df, temperature == temp)$log2_value)
    }))),
    unlist(lapply(unique_temp, function(temp){
      mean(filter(df, temperature == temp)$log2_value)
    })),
    rep(0, len_temp))
    return(start_par)
}

paste_rmNA <- function(x, sep = "|"){
  # function that pastes non-dedundant and na-filtered
  # vector elements into a string
  return(paste(x[!is.na(x)], collapse = sep))
}