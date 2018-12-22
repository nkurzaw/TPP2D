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
                  gg_theme = theme_classic(), offset = 1){
  sx <- sort(x)
  sy <- sort(y)
  lenx <- length(sx)
  leny <- length(sy)
  
  if (leny < lenx)
    sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)
    sy <- approx(1L:leny, sy, n = lenx)$y
  df <- data.frame(sx, sy)
  
  ggplot(df, aes(sx, sy)) +
    geom_point(alpha = alpha) +
    geom_line(aes(x, y),
              linetype = "dashed",
              color = "gray",
              data = data.frame(x = seq(0, max(c(sx,sy))),
                                y = seq(0, max(c(sx,sy))))) +
    coord_fixed(xlim = c(0, max(c(sx,sy)) + offset),
                ylim = c(0, max(c(sx,sy)) + offset),
                expand = FALSE) +
    xlab(xlab) +
    ylab(ylab) +
    gg_theme
}

#' Plot 2D thermal profile intensities of a protein
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

#' Plot 2D thermal profile ratios of a protein
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