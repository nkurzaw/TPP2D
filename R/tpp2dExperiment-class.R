#' S4 TPP2D Experiment Class
#'
#' @slot configFilePath character.
#' @slot configTable data.frame.
#' @slot tidyDataTable data.frame.
#' @slot modelParamsDf data.frame
#' @slot resultTable
#'
#' @return an object of class tpp2dExperiment
#' with the following slots:
#' 1) configFilePath: 
#'
#' @export
#'
#' @examples
#' tpp2dObj <- new("tpp2dExperiment", 
#'    configFilePath = system.file(package = "TPP2D"))
tpcaResult <- setClass("tpp2dExperiment",
                       slots = list(
                           configFilePath = "character",
                           configTable = "data.frame",
                           tidyDataTable = "data.frame",
                           modelParamsDf = "data.frame",
                           resultTable = "data.frame"
                           ))
