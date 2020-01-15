#' S4 TPP2D Experiment Class
#'
#' @slot configTable data.frame.
#' @slot tidyDataTable data.frame.
#' @slot modelParamsDf data.frame
#' @slot resultTable data.frame
#' @slot bootstrapNullDf data.frame
#' @slot hitTable data.frame
#'
#' @return an object of class tpp2dExperiment
#' with the following slots:
#' 1) configTable: 
#'
#' @export
#'
#' @examples
#' tpp2dObj <- new("tpp2dExperiment", 
#'    configFilePath = system.file(package = "TPP2D"))
tpcaResult <- setClass("tpp2dExperiment",
                       slots = list(
                           configTable = "data.frame",
                           tidyDataTable = "data.frame",
                           modelParamsDf = "data.frame",
                           resultTable = "data.frame",
                           bootstrapNullDf = "data.frame",
                           hitTable = "data.frame"
                           ))
