## Internal functions for optimization of fits

.trim_sum <- function(x, na.rm = TRUE)
  # Trimmed sum (sum up all but the highest values)
{
  sum(x[-which(x == max(x, na.rm = TRUE))],
      na.rm = TRUE)
}

.min_RSS_h0 <- function(data, par, len_temp) {
  # Optimization function for fitting an intercept model to a protein's 2D
  # thermal profile by minimizing the sum of squared errors
  temp_i <- log2_value <- NULL
  
  beta_0 <- par[seq_len(len_temp)]
  
  sum(with(data, (beta_0[temp_i] - log2_value) ^ 2))
}

.min_RSS_h0_trim <- function(data, par, len_temp) {
  # Optimization function for fitting an intercept model to a protein's 2D
  # thermal profile by minimizing the trimmed sum of squared errors
  temp_i <- log2_value <- NULL
  
  beta_0 <- par[seq_len(len_temp)]
  
  .trim_sum(with(data, (beta_0[temp_i] - log2_value) ^ 2))
}

.min_RSS_h1 <- function(data, par, len_temp) {
  # Optimization function for fitting an dose-response model to a
  # protein's 2D thermal profile by minimizing the sum of squared errors
  temp_i <- log2_value <- log_conc <- NULL
  
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp * 2)]
  
  sum(with(data, (
    beta_0[temp_i] + (alpha[temp_i] * beta_max) /
      (1 + exp(-slope * (log_conc - zeta))) -
      log2_value
  ) ^ 2))
}

.min_RSS_h1_slope_pEC50 <- function(data, par, len_temp) {
  # Optimization function for fitting an dose-response model to a
  # protein's 2D thermal profile by minimizing the sum of squared errors
  temp_i <- log2_value <- log_conc <- temperature <- NULL
  
  zeta <- par[1]
  zeta_slope <- par[2]
  slope <- par[3]
  beta_max <- par[4]
  beta_0 <- par[5:(len_temp + 4)]
  alpha <- par[(5 + len_temp):(4 + len_temp * 2)]
  
  sum(with(data,
           (
             beta_0[temp_i] + (alpha[temp_i] * beta_max) /
               (1 + exp(-slope * (
                 log_conc - (zeta + zeta_slope * temperature)
               ))) -
               log2_value
           )^ 2))
}


.min_RSS_h1_trim <- function(data, par, len_temp) {
  # Optimization function for fitting an dose-response model to a
  # protein's 2D thermal profile by minimizing the trimmed sum of
  # squared errors
  temp_i <- log2_value <- log_conc <- NULL
  
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp * 2)]
  
  .trim_sum(with(data, (
    beta_0[temp_i] + (alpha[temp_i] * beta_max) /
      (1 + exp(-slope * (log_conc - zeta))) -
      log2_value
  ) ^ 2))
}

.min_RSS_h0_gradient <- function(data, par, len_temp) {
  # Analytically solved gradient function for .min_RSS_h1
  temp_i <- log2_value <- NULL
  
  beta_0 <- par[seq_len(len_temp)]
  
  data <-
    full_join(data,
              expand.grid(
                log_conc = unique(data$log_conc),
                temp_i = seq_len(max(data$temp_i))
              ),
              by = c("log_conc", "temp_i"))
  
  outer_dev <-
    with(data, 2 * (beta_0[temp_i] - log2_value))
  
  d_beta_0 <- apply(matrix(outer_dev, nrow = 5, byrow = TRUE),
                    2, sum, na.rm = TRUE)
  
  return(d_beta_0)
  
}

.min_RSS_h0_gradient_trim <- function(data, par, len_temp) {
  # Analytically solved gradient function for .min_RSS_h1
  temp_i <- log2_value <- NULL
  
  beta_0 <- par[seq_len(len_temp)]
  
  data <-
    full_join(data,
              expand.grid(
                log_conc = unique(data$log_conc),
                temp_i = seq_len(max(data$temp_i))
              ),
              by = c("log_conc", "temp_i"))
  
  outer_dev <-
    with(data, 2 * (beta_0[temp_i] - log2_value))
  
  d_beta_0 <- apply(matrix(outer_dev, nrow = 5, byrow = TRUE),
                    2, .trim_sum, na.rm = TRUE)
  
  return(d_beta_0)
  
}

.min_RSS_h1_gradient <- function(data, par, len_temp) {
  # Analytically solved gradient function for .min_RSS_h1
  temp_i <- log2_value <- log_conc <- NULL
  
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp * 2)]
  
  data <-
    full_join(data,
              expand.grid(
                log_conc = unique(data$log_conc),
                temp_i = seq_len(max(data$temp_i))
              ),
              by = c("log_conc", "temp_i"))
  outer_dev <-
    with(data, 2 * (beta_0[temp_i] + (alpha[temp_i] * beta_max) /
                      (1 + exp(
                        -slope * (log_conc - zeta)
                      )) -
                      log2_value))
  d_zeta <- sum(outer_dev * with(data,-(alpha * beta_max) / (1 + exp(
    -slope * (log_conc - zeta)
  )) ^ 2 *
    exp(-slope * (log_conc - zeta)) * slope), na.rm = TRUE)
  d_slope <- sum(outer_dev * with(data,-(alpha * beta_max) / (1 + exp(
    -slope * (log_conc - zeta)
  )) ^ 2 *
    exp(-slope * (log_conc - zeta)) * (-(log_conc - zeta))), na.rm = TRUE)
  d_beta_max <- sum(outer_dev * with(data, alpha / (1 + exp(
    -slope * (log_conc - zeta)
  ))), na.rm = TRUE)
  d_beta_0 <-
    apply(matrix(outer_dev, nrow = 5, byrow = TRUE), 2, sum, na.rm = TRUE)
  d_alpha <- apply(matrix(
    outer_dev * with(data, beta_max / (1 + exp(
      -slope * (log_conc - zeta)
    ))),
    nrow = 5,
    byrow = TRUE
  ), 2, sum, na.rm = TRUE)
  return(c(d_zeta, d_slope, d_beta_max, d_beta_0, d_alpha))
  
}

.min_RSS_h1_gradient_trim <- function(data, par, len_temp) {
  # Analytically solved gradient function for .min_RSS_h1_trim
  temp_i <- log2_value <- log_conc <- NULL
  
  zeta <- par[1]
  slope <- par[2]
  beta_max <- par[3]
  beta_0 <- par[4:(len_temp + 3)]
  alpha <- par[(4 + len_temp):(3 + len_temp * 2)]
  
  data <-
    full_join(data,
              expand.grid(
                log_conc = unique(data$log_conc),
                temp_i = seq_len(max(data$temp_i))
              ),
              by = c("log_conc", "temp_i"))
  outer_dev <-
    with(data, 2 * (beta_0[temp_i] + (alpha[temp_i] * beta_max) /
                      (1 + exp(
                        -slope * (log_conc - zeta)
                      )) -
                      log2_value))
  d_zeta <- .trim_sum(outer_dev * with(data,-(alpha * beta_max) / (1 + exp(
    -slope * (log_conc - zeta)
  )) ^ 2 *
    exp(-slope * (log_conc - zeta)) * slope), na.rm = TRUE)
  d_slope <- .trim_sum(outer_dev * with(data,-(alpha * beta_max) / (1 + exp(
    -slope * (log_conc - zeta)
  )) ^ 2 *
    exp(-slope * (log_conc - zeta)) * (-(log_conc - zeta))), na.rm = TRUE)
  d_beta_max <- .trim_sum(outer_dev * with(data, alpha / (1 + exp(
    -slope * (log_conc - zeta)
  ))), na.rm = TRUE)
  d_beta_0 <-
    apply(matrix(outer_dev, nrow = 5, byrow = TRUE), 2, .trim_sum, na.rm = TRUE)
  d_alpha <- apply(matrix(
    outer_dev * with(data, beta_max / (1 + exp(
      -slope * (log_conc - zeta)
    ))),
    nrow = 5,
    byrow = TRUE
  ),
  2,
  .trim_sum,
  na.rm = TRUE)
  return(c(d_zeta, d_slope, d_beta_max, d_beta_0, d_alpha))
  
}

.min_RSS_h1_slope50_gradient <- function(data, par, len_temp) {
  # Analytically solved gradient function for .min_RSS_h1_slope_pEC50
  temp_i <- log2_value <- log_conc <- temperature <- NULL
  
  zeta <- par[1]
  zeta_slope <- par[2]
  slope <- par[3]
  beta_max <- par[4]
  beta_0 <- par[5:(len_temp + 4)]
  alpha <- par[(5 + len_temp):(4 + len_temp * 2)]
  
  data <-
    full_join(data,
              expand.grid(
                log_conc = unique(data$log_conc),
                temp_i = seq_len(max(data$temp_i))
              ),
              by = c("log_conc", "temp_i"))
  outer_dev <-
    with(data,
         2 * (beta_0[temp_i] + (alpha[temp_i] * beta_max) /
                (1 + exp(
                  -slope * (log_conc - (zeta + zeta_slope * temperature))
                )) -
                log2_value))
  
  d_zeta <- sum(outer_dev *
                  with(data,-(alpha * beta_max) /
                         (1 + exp(
                           -slope * (log_conc - (zeta + zeta_slope * temperature))
                         )) ^ 2 *
                         exp(-slope * (
                           log_conc - (zeta + zeta_slope * temperature)
                         )) * slope),
                na.rm = TRUE)
  
  d_zeta_slope <- sum(outer_dev * with(
    data,
    -(alpha * beta_max) / (1 + exp(-slope * (
      log_conc - (zeta + zeta_slope * temperature)
    ))) ^ 2 *
      exp(-slope * (
        log_conc - (zeta + zeta_slope * temperature)
      )) * (slope * temperature)
  ), na.rm = TRUE)
  
  d_slope <- sum(outer_dev * with(data,-(alpha * beta_max) / (1 + exp(
    -slope * (log_conc - (zeta + zeta_slope * temperature))
  )) ^ 2 *
    exp(-slope * (
      log_conc - (zeta + zeta_slope * temperature)
    )) * (-(
      log_conc - (zeta + zeta_slope * temperature)
    ))), na.rm = TRUE)
  d_beta_max <- sum(outer_dev * with(data, alpha / (1 + exp(
    -slope * (log_conc - (zeta + zeta_slope * temperature))
  ))), na.rm = TRUE)
  d_beta_0 <-
    apply(matrix(outer_dev, nrow = 5, byrow = TRUE), 2, sum, na.rm = TRUE)
  d_alpha <- apply(matrix(
    outer_dev * with(data, beta_max / (1 + exp(
      -slope * (log_conc - (zeta + zeta_slope * temperature))
    ))),
    nrow = 5,
    byrow = TRUE
  ), 2, sum, na.rm = TRUE)
  return(c(d_zeta, d_zeta_slope, d_slope, d_beta_max, d_beta_0, d_alpha))
  
}
