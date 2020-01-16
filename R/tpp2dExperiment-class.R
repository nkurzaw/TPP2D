#' S4 TPP2D Experiment Class
#'
#' @slot configTable data.frame.
#' @slot idVar character.
#' @slot intensityStr character.
#' @slot fcStr character.
#' @slot nonZeroCols character.
#' @slot geneNameVar character.
#' @slot qualColName character.
#' @slot naStrs character.
#' @slot concFactor numeric.
#' @slot medianNormalizeFC logical.
#' @slot filterContaminants logical.
#' @slot minObs numeric.
#' @slot independentFiltering logical.
#' @slot fcThres numeric.
#' @slot optim_fun_h0 function.
#' @slot optim_fun_h1 function.
#' @slot slopEC50 logical.
#' @slot maxit numeric.
#' @slot BPPARAM character.
#' @slot B numeric
#' @slot byMsExp logical.
#' @slot alpha numeric.
#' @slot tidyDataTable data.frame.
#' @slot modelParamsDf data.frame
#' @slot resultTable data.frame
#' @slot bootstrapNullDf data.frame
#' @slot hitTable data.frame
#'
#' @return an object of class tpp2dExperiment
#'
#' @export
#'
#' @examples
#' tpp2dObj <- new("tpp2dExperiment")
tpcaResult <- setClass("tpp2dExperiment",
                       slots = list(
                           configTable = "data.frame",
                           idVar = "character",
                           intensityStr = "character",
                           fcStr = "character",
                           nonZeroCols = "character",
                           geneNameVar = "character",
                           qualColName = "character",
                           naStrs = "character",
                           concFactor = "numeric",
                           medianNormalizeFC = "logical",
                           filterContaminants = "logical",
                           minObs = "numeric",
                           independentFiltering = "logical",
                           fcThres = "numeric",
                           optim_fun_h0 = "function",
                           optim_fun_h1 = "function",
                           slopEC50 = "logical",
                           maxit = "numeric",
                           BPPARAM = "character",
                           B = "numeric",
                           byMsExp = "logical",
                           alpha = "numeric",
                           tidyDataTable = "data.frame",
                           modelParamsDf = "data.frame",
                           resultTable = "data.frame",
                           bootstrapNullDf = "data.frame",
                           hitTable = "data.frame"
                           ))
