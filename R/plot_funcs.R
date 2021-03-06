#' Plot qq-plot of true data and bootstrapped null with ggplot
#' 
#' @param x vector containing values of values of first 
#' distribution to compare
#' @param y vector containing values of values of secound 
#' distribution to compare
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param alpha transparency paramenter between 0 and 1
#' @param gg_theme ggplot theme, default is theme_classic()
#' @param offset offset for x and y axis on top of maximal 
#' values
#' @param plot_diagonal logical parameter indicating whether
#' an identity line should be plotted
#' 
#' @return A ggplot displaying the qq-plot of a true and a
#' a bootstrapped null distribution
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' recomputeSignalFromRatios(simulated_cell_extract_df)
#'
#' @import ggplot2
#' @importFrom stats approx
gg_qq <- function(x, y, 
                  xlab = "F-statistics from sampled Null distr.",
                  ylab = "observed F-statistics", alpha = 0.25,
                  gg_theme = theme_classic(), offset = 1,
                  plot_diagonal = TRUE){
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  
  if (leny < lenx)
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)
    sy <- approx(1L:leny, sy, n = lenx)$y
  df <- data.frame(sx, sy)
  
  p <- ggplot(df, aes(sx, sy)) +
    geom_point(alpha = alpha) +
    coord_fixed(xlim = c(0, max(c(sx,sy)) + offset),
                ylim = c(0, max(c(sx,sy)) + offset),
                expand = FALSE) +
    xlab(xlab) +
    ylab(ylab) +
    gg_theme
  
  if(plot_diagonal){
      p <- p + 
          geom_line(aes(x, y),
                    linetype = "dashed",
                    color = "gray",
                    data = data.frame(
                        x = seq(0, max(c(sx,sy))),
                        y = seq(0, max(c(sx,sy))))) 
  }
  return(p)
}

#' Plot 2D thermal profile intensities of a protein 
#' of choice
#' 
#' @param df tidy data frame of a 2D-TPP dataset 
#' @param name gene name (clustername) of protein that 
#' should be visualized
#' 
#' @return A ggplot displaying the thermal profile of
#' a protein of choice in a datset of choice
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' plot2dTppProfile(simulated_cell_extract_df, "protein1")
#'
#' @import ggplot2
plot2dTppProfile <- function(df, name){
  clustername <- log_conc <- log2_value <- 
    temperature <- NULL
  ggplot(filter(df, clustername == name),
         aes(log_conc, log2_value)) +
    geom_point() +
    facet_wrap(~temperature)
}

#' Plot 2D thermal profile ratios of a protein of choice
#' 
#' @param df tidy data frame of a 2D-TPP dataset 
#' @param name gene name (clustername) of protein that 
#' should be visualized
#' 
#' @return A ggplot displaying the thermal profile ratios of
#' a protein of choice in a datset of choice
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' plot2dTppRelProfile(simulated_cell_extract_df, "protein1")
#'
#' @import ggplot2
plot2dTppRelProfile <- function(df, name){
  clustername <- log_conc <- rel_value <- 
    temperature <- NULL
  ggplot(filter(df, clustername == name),
         aes(log_conc, rel_value)) +
    geom_point() +
    facet_wrap(~temperature)
}

#' Plot H0 or H1 fit of 2D thermal profile intensities of 
#' a protein of choice
#' 
#' @param df tidy data frame of a 2D-TPP dataset 
#' @param name gene name (clustername) of protein that 
#' should be visualized
#' @param model_type character string indicating whether
#' the "H0" or the "H1" model should be fitted
#' @param optim_fun optimization function that should be used
#' for fitting either the H0 or H1 model
#' @param optim_fun_2 optional additional optimization function 
#' that will be run with paramters retrieved from optim_fun and 
#' should be used for fitting the H1 model with the trimmed sum
#' model, default is NULL
#' @param maxit maximal number of iterations the optimization
#' should be given, default is set to 500
#' @param xlab character string of x-axis label of plot
#' @param ylab character string of y-axis label of plot
#' @param dot_size numeric indicating the size of the data points
#' to plot
#' @param line_type character string defining the line type 
#' of the fitted curve, default "dashed"
#' @param fit_color character string defining the color
#' of the fitted curve
#' 
#' @return A ggplot displaying the thermal profile of
#' a protein of choice in a datset of choice
#' 
#' @export
#'
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' plot2dTppProfile(simulated_cell_extract_df, "protein1")
#'
#' @import ggplot2
plot2dTppFit <- function(df, name,
                         model_type = "H0",
                         optim_fun = .min_RSS_h0,
                         optim_fun_2 = NULL,
                         maxit = 500,
                         xlab = "-log10(conc.)",
                         ylab = "log2(summed intensities)",
                         dot_size = 1,
                         line_type = "solid",
                         fit_color = "gray30"){
  
  clustername <- temperature <- temp_i <- 
    log_conc <- y_hat <- log2_value <- NULL
  
  .checkDfColumns(df)
  
  if(model_type == "H1"){
    optim_fun <- .min_RSS_h1_slope_pEC50
    slopEC50 <-  TRUE
  }else{
    slopEC50 <-  FALSE
  }
  
  df_fil <- filter(df, clustername == name) %>% 
    mutate(temp_i = dense_rank(temperature))
  unique_temp <- unique(df_fil$temperature)
  len_temp <- length(unique_temp)
  
  if(model_type == "H0"){
      start_par = vapply(unique_temp, function(x)
        mean(filter(df_fil, temperature == x)$log2_value), 1)
      h0_model = try(optim(par = start_par,
                           fn = optim_fun,
                           len_temp = len_temp,
                           data = df_fil,
                           method = "L-BFGS-B",
                           control = list(maxit = maxit)))
        
      fit_df <-
        .getFitDf(df_fil, model_type = model_type,
                 optim_model = h0_model,
                 slopEC50 = slopEC50)
  }else if(model_type == "H1"){
    h1_model <- .fitEvalH1(df_fil = df_fil,
                          unique_temp = unique_temp, 
                          len_temp = len_temp,
                          optim_fun = optim_fun, 
                          optim_fun_2 = optim_fun_2,
                          slopEC50 = slopEC50, 
                          maxit = maxit)
    fit_df <-
      .getFitDf(df_fil, model_type = model_type,
               optim_model = h1_model,
               slopEC50 = slopEC50)
  }else{
    stop("Please specify a valid model_type! Either H0 or H1!")
  }
  ggplot(fit_df, aes(log_conc, y_hat)) +
    geom_point(aes(log_conc, log2_value), 
               size = dot_size,
               shape = 21,
               data = df_fil) +
    geom_line(color = fit_color,
              linetype = line_type) +
    facet_wrap(~temperature) +
    ggtitle(name) +
    labs(x = xlab, y = ylab)
}


.getFitDf <- function(df_fil, conc_vec = seq(-9, -3, by = 0.1),
                     model_type = "H0", optim_model, 
                     slopEC50 = TRUE){
  # internal function to retrieve data frame with data
  # and predicted model values
  temp_i <- len_temp <- log_conc <- temperature <- NULL
  unique_temp <- unique(df_fil$temperature)
  unique_temp_len <- length(unique_temp)
  conc_len = length(conc_vec)
  
  if(model_type == "H0"){
    fit_df <-
      tibble(
        log_conc = rep(conc_vec, unique_temp_len),
        temperature = rep(unique_temp,
                          each = conc_len),
        temp_i = rep(seq_len(unique_temp_len), 
                     each = conc_len),
        len_temp = unique_temp_len) %>%
      mutate(y_hat = optim_model$par[temp_i])
    
  }else if(model_type == "H1"){
    if(slopEC50){
      fit_df <-
        tibble(
          log_conc = rep(conc_vec, unique_temp_len),
          temperature = rep(unique_temp, each = conc_len),
          temp_i = rep(seq_len(unique_temp_len), 
                       each = conc_len),
          len_temp = unique_temp_len) %>%
        mutate(y_hat = .evalH1SlopeEC50Model(
          optim_model = optim_model,
          temp_i = temp_i, len_temp = len_temp,
          log_conc = log_conc, temperature = temperature))
    }else{
      fit_df <-
        tibble(
          log_conc = rep(conc_vec, unique_temp_len),
          temperature = rep(unique_temp, each = conc_len),
          temp_i = rep(seq_len(unique_temp_len), 
                     each = conc_len),
          len_temp = unique_temp_len) %>%
        mutate(y_hat = .evalH1Model(
          optim_model = optim_model,
          temp_i = temp_i, len_temp = len_temp,
          log_conc = log_conc))
    }
  }else{
    stop("'model type' has to be either 'H0' or 'H1'!")
  }
  return(fit_df)
}

.evalH1Model <- function(optim_model, temp_i, log_conc, len_temp){
  optim_model$par[3 + temp_i] + 
    (optim_model$par[3 + len_temp + temp_i] * optim_model$par[3])/
    (1 + exp(-optim_model$par[2] * (log_conc - optim_model$par[1])))
}

.evalH1SlopeEC50Model <- function(optim_model, temp_i, len_temp,
                                 log_conc, temperature){
  optim_model$par[4 + temp_i] + 
    (optim_model$par[4 + len_temp + temp_i] * optim_model$par[4])/
    (1 + exp(-optim_model$par[3] * 
               (log_conc - (optim_model$par[1] + 
                  optim_model$par[2] * temperature))))
}

#' Plot heatmap of 2D thermal profile fold changes of 
#' a protein of choice
#' 
#' @param df tidy data frame of a 2D-TPP dataset 
#' @param name gene name (clustername) of protein that 
#' should be visualized
#' @param drug_name character string of profiled drug name
#' @param midpoint midpoint of fold changes for color
#' scaling, default: 1
#' 
#' @return A ggplot displaying the thermal profile 
#' as a heatmap of fold changes of
#' a protein of choice in a dataset of choice
#' 
#' @export
#' @examples
#' 
#' data("simulated_cell_extract_df")
#' plot2dTppFcHeatmap(simulated_cell_extract_df, 
#'  "tp2", drug_name = "drug1")
#'
#' @import ggplot2
#' @import dplyr
plot2dTppFcHeatmap <- function(df, name, 
                               drug_name = "",
                               midpoint = 1){
  clustername <- temperature <- conc <- 
    rel_value <- fc <- NULL
  
  heat_df <- df %>% 
    filter(clustername == name) %>% 
    dplyr::select(clustername, temperature, 
                  conc, rel_value) %>% 
    mutate(temperature = factor(
      temperature, levels = rev(sort(unique(temperature)))),
      conc = as.factor(conc)) %>% 
    group_by(clustername, temperature, conc) %>% 
    summarise(fc = mean(rel_value, na.rm = TRUE)) %>% 
    ungroup
  
  ggplot(heat_df, aes(conc, temperature)) +
    geom_tile(aes(fill = fc)) +
    scale_fill_gradient2(
      "Fold change", 
      low = "firebrick", mid = "khaki", 
      high = "darkgreen", midpoint = midpoint) +
    labs(x = paste(drug_name, "conc."),
         y = expression("Temperature ("*~degree*C*")")) +
    ggtitle(name) +
    theme_minimal()
}

#' Plot Volcano plot of TPP2D results
#' @param fdr_df data frame obtained from `getFDR`
#' @param hits_df hits_df data frame obtained from `findHits`
#' @param alpha transparency level of plotted points
#' @param title_string character argument handed over to ggtitle
#' @param x_lim vector with two numerics indicating the x axis limits
#' @param y_lim vector with two numerics indicating the y axis limits
#' @param facet_by_obs logical indicating whether plot should be facetted
#' by number of observations, default: FALSE
#'
#' @return a ggplot displaying a volcano plot of the results
#' obtained after a TPP2D analysis
#'
#' @import ggplot2
#' @import dplyr
#' 
#' @export
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' temp_df <- simulated_cell_extract_df %>%
#'   filter(clustername %in% paste0("protein", 1:5)) %>%
#'   group_by(representative) %>%
#'   mutate(nObs = n()) %>%
#'   ungroup
#' example_params <- getModelParamsDf(temp_df)
#' example_fstat <- computeFStatFromParams(example_params)
#' example_null <- bootstrapNullAlternativeModel(
#'    df = temp_df, params_df = example_params,
#'    B = 2)
#' fdr_df <- getFDR(example_fstat, example_null)
#' hits_df <- findHits(fdr_df, 0.1)
#' plot2dTppVolcano(fdr_df = fdr_df, hits_df = hits_df)
plot2dTppVolcano <- function(fdr_df, hits_df, 
                             alpha = 0.5,
                             title_string = "",
                             x_lim = NULL,
                             y_lim = NULL,
                             facet_by_obs = FALSE){
  dataset <- rssH0 <- rssH1 <- F_statistic <- group <- 
    clustername <- slopeH1 <- NULL
  
  stab_colors <- c("steelblue", "orange")
  names(stab_colors) <- c(
    "stabilized",
    "destabilized")
  
  fdr_true_df <- fdr_df %>% 
    filter(dataset == "true")
  
  if(is.null(x_lim)){
    x_lim <- c(
      floor(min(sign(fdr_true_df$slopeH1)*
            sqrt(fdr_true_df$rssH0 - fdr_true_df$rssH1)) - 0.5),
      round(max(sign(fdr_true_df$slopeH1)*
            sqrt(fdr_true_df$rssH0 - fdr_true_df$rssH1)) + 0.5)
    )
  }
  if(is.null(y_lim)){
    y_lim <- c(
      0,
      round(max(log2(fdr_true_df$F_statistic + 1)) + 0.5) 
    )
  }
  
  p <- ggplot(fdr_true_df %>% 
           mutate(group = case_when(slopeH1 > 0 ~ "stabilized",
                                    slopeH1 < 0 ~ "destabilized")), 
         aes(sign(slopeH1)*sqrt(rssH0 - rssH1), log2(F_statistic + 1))) +
    geom_point(color = "gray", alpha = alpha) + 
    geom_point(aes(color = group), alpha = alpha, 
               data = hits_df %>% 
                 mutate(group = case_when(
                   slopeH1 > 0 ~ "stabilized",
                   slopeH1 < 0 ~ "destabilized"))) + 
    geom_text(
      aes(label = clustername),
      data = hits_df, nudge_x = 0.5, nudge_y = 0.5) +
    scale_color_manual("", values = stab_colors) +
    labs(x = bquote(sign(kappa) %.% sqrt(~'RSS'^0~' - RSS'^1~'')),
         y = expression('log'[2]~'('*italic(F)*'-statistic + 1)')) +
    ggtitle(title_string) +
    coord_cartesian(xlim = x_lim,
                    ylim = y_lim) +
    theme(legend.position = "bottom")
  if(facet_by_obs){
    p <- p +
      facet_wrap(~nObsRound)
  }
  return(p)
}
