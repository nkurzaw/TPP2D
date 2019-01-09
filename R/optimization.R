## Internal functions for optimization of fits

trim_sum <- function(x) 
# Trimmed sum (sum up all but the highest values)
{
  sum(x[-which(x == max(x))])
}

#' @importFrom Rcpp cppFunction
compute_residuals_cpp <- Rcpp::cppFunction("
    NumericVector compute_residuals_cpp(
       NumericVector tempi,
       double zeta,
       double slope,
       double betamax,
       NumericVector betaz,
       NumericVector alpha,
       NumericVector logconc,
       NumericVector log2value) {
        int n = tempi.size();
        NumericVector out(n);
        for(int i = 0; i < n; ++i){
          int tempindex = tempi[i] - 1;
          double quot = alpha[tempindex] * betamax;
          double divi = 1 + exp( -slope * (logconc[i] - zeta));
          out[i] = pow((betaz[tempindex] + (quot / divi)) - log2value[i], 2.0);
        }
        return out;
      }
")


min_RSS_h0_trim <- function(data, par, unique_temperature, len_temp)
# Optimization function for fitting an intercept model to a protein's 2D 
# thermal profile by minimizing the trimmed sum of squared errors
{
  beta_0 <- par[1:len_temp]
  
  trim_sum(
    with(data, (beta_0[temp_i] - log2_value)^2)
  )
}

min_RSS_h1 <- function(data, par, unique_temperature, len_temp)
# Optimization function for fitting an dose-response model to a 
# protein's 2D thermal profile by minimizing the sum of squared errors
{
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp*2)]
  
  sum(
    with(data, (beta_0[temp_i] + (alpha[temp_i] * beta_max)/
                  (1 + exp(-slope * (log_conc - zeta))) -
                  log2_value)^2)
  )
}

min_RSS_h1_cpp <- function(data, par, unique_temperature, len_temp)
  # Optimization function for fitting an dose-response model to a 
  # protein's 2D thermal profile by minimizing the sum of squared errors
{
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp*2)]
  
  sum(
    with(data, compute_residuals_cpp(temp_i, zeta, slope, beta_max, 
                                     beta_0, alpha, log_conc, log2_value))
  )
}

min_RSS_h1_trim <- function(data, par, unique_temperature, len_temp)
# Optimization function for fitting an dose-response model to a 
# protein's 2D thermal profile by minimizing the trimmed sum of 
# squared errors
{
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp*2)]
  
  trim_sum(
    with(data, (beta_0[temp_i] + (alpha[temp_i] * beta_max)/
                  (1 + exp(-slope * (log_conc - zeta))) -
                  log2_value)^2)
  )
}