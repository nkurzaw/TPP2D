#' Run complete TPP2D analysis
#'
#' @param df tidy data_frame retrieved after import of a 2D-TPP 
#' dataset, potential filtering and addition of a column "nObs"
#' containing the number of observations per protein
#' @param configTable character string of a file path to a config table
#' @param data possible list of datasets from different MS runs 
#' corresponding to a 2D-TPP dataset, circumvents loading datasets 
#' referencend in config table, default is NULL
#' @param idVar character string indicating which data column provides the 
#'   unique identifiers for each protein.
#' @param intensityStr character string indicating which columns contain 
#'   raw intensities measurements
#' @param fcStr character string indicating which columns contain the actual 
#'   fold change values. Those column names containing the suffix \code{fcStr} 
#'   will be regarded as containing fold change values.
#' @param naStrs character vector indicating missing values in the data table. 
#'   When reading data from file, this value will be passed on to the argument 
#'   \code{na.strings} in function \code{read.delim}.
#' @param qualColName character string indicating which column can be used for 
#'   additional quality criteria when deciding between different non-unique 
#'   protein identifiers.
#' @param medianNormalizeFC perform median normalization (default: TRUE).
#' @param addCol character string indicating additional column to import
#' @param filterContaminants logical variable indicating whether data 
#' should be filtered to exclude contaminants (default: TRUE).
#' @param recomputeSignalRatios logical variable indicaiting whether 
#' signals should be recomputed from relative fold changes, recommended
#' if Isobarquant was used for protein quantification
#' @param minObs number of minimal observations per protein to include it
#' in the analysis
#' @param independentFiltering logical variable indicating whether
#' independent filtering should be performed based on minimal
#' fold changes per protein profile
#' @param fcThres numeric value of minimal fold change 
#' (or inverse fold change) a protein has to show to be kept 
#' upon independent filtering
#' @param nonZeroCols column like default qssm that should be imported and
#' requested to be non-zero in analyzed data
#' @param geneNameVar character string of the column name that describes
#' the gene name of a given protein in the raw data files
#' @param concFactor numeric value that indicates how concentrations need to 
#' be adjusted to yield total unit e.g. default mmol - 1e6
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
#' @param slopEC50 logical flag indicating whether the h1 model is
#' fitted with a linear model describing the shift od the pEC50 over 
#' temperatures
#' @param BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
#' @param B numeric value indicating number of rounds of bootstraps
#' that should be performed to estimate the null distribution
#' @param byMsExp logical indicating whether bootstrapping should be
#' performed within MS experiments
#' @param alpha FDR level that should be controlled
#' 
#' @return a tpp2dExperiment object
#' 
#' @examples 
#' data("simulated_cell_extract_df")
#' runTPP2D(df = simulated_cell_extract_df %>% 
#'    filter(representative %in% 1:3),
#'    B = 1)
#' 
#' @export
#' 
#' @importFrom methods new
#' 
runTPP2D <- function(df = NULL,
                     configTable = NULL, 
                     data = NULL,
                     idVar = "protein_id",
                     intensityStr = "signal_sum_",
                     fcStr = "rel_fc_",
                     nonZeroCols = "qusm",
                     geneNameVar = "gene_name",
                     addCol = NULL,
                     qualColName = "qupm",
                     naStrs = c("NA", "n/d", "NaN"),
                     concFactor = 1e6,
                     medianNormalizeFC = TRUE,
                     filterContaminants = TRUE,
                     recomputeSignalRatios = FALSE,
                     minObs = 20,
                     independentFiltering = FALSE,
                     fcThres = 1.5,
                     optim_fun_h0 = .min_RSS_h0,
                     optim_fun_h1 = .min_RSS_h1_slope_pEC50,
                     optim_fun_h1_2 = NULL,
                     gr_fun_h0 = NULL,
                     gr_fun_h1 = NULL,
                     gr_fun_h1_2 = NULL,
                     slopEC50 = TRUE,
                     maxit = 750,
                     BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
                     B = 20,
                     byMsExp = TRUE,
                     alpha = 0.1){
    
    if(is.null(df) & !is.null(configTable)){
        message("Importing data")
        import_df <- import2dDataset(configTable = configTable, 
                                     data = data,
                                     idVar = idVar,
                                     intensityStr = intensityStr,
                                     fcStr = fcStr,
                                     nonZeroCols = nonZeroCols,
                                     geneNameVar = geneNameVar,
                                     addCol = addCol,
                                     qualColName = qualColName,
                                     naStrs = naStrs,
                                     concFactor = concFactor,
                                     medianNormalizeFC = medianNormalizeFC,
                                     filterContaminants = filterContaminants)
        
        cfgTab <- TPP_importCheckConfigTable(
            infoTable = configTable, type = "2D")
        .checkDfColumns(import_df)
        df <- import_df
    }else if(!is.null(df)){
        cfgTab <- data.frame()
    }else if(is.null(configTable)){
        stop(paste("Either a data.frame (df) or a configTable",
                   "need to be supplied in order to perform",
                   "the analysis!\n"))
    }
    
    if(identical(optim_fun_h1, .min_RSS_h1_slope_pEC50)){
        slopEC50 = TRUE
    }else{
        slopEC50 = FALSE
    }
    
    if(recomputeSignalRatios){
        df <- recomputeSignalFromRatios(df)
    }
    
    if(independentFiltering){
        message(paste("Independent Filtering: removing proteins without", 
                      "any values crossing the set threshold.\n"))
        df <- .independentFilter(df, fcThres = fcThres) 
    }
    
    message("Computing null and alternative model parameters\n")
    
    params_df <- getModelParamsDf(df,
                                  optim_fun_h0 = optim_fun_h0,
                                  optim_fun_h1 = optim_fun_h1,
                                  optim_fun_h1_2 = optim_fun_h1_2,
                                  gr_fun_h0 = gr_fun_h0,
                                  gr_fun_h1 = gr_fun_h1,
                                  gr_fun_h1_2 = gr_fun_h1_2,
                                  slopEC50 = slopEC50,
                                  maxit = maxit,
                                  qualColName = qualColName)
    
    message("Computing F statistics\n")
    
    fstat_df <- computeFStatFromParams(params_df)
    
    message("Bootstrapping null distribution\n")
    
    null_df <- bootstrapNullAlternativeModel(df, 
                 params_df,
                 maxit = maxit,
                 minObs = minObs,
                 optim_fun_h0 = optim_fun_h0,
                 optim_fun_h1 = optim_fun_h1,
                 optim_fun_h1_2 = optim_fun_h1_2,
                 gr_fun_h0 = gr_fun_h0,
                 gr_fun_h1 = gr_fun_h1,
                 gr_fun_h1_2 = gr_fun_h1_2,
                 BPPARAM = BPPARAM,
                 B = B,
                 byMsExp = byMsExp)
    
    message("Computing FDR\n")
    
    fdr_df <- computeFdr(df_out = fstat_df,
                         df_null = null_df)
    hits_df <- findHits(fdr_df = fdr_df,
                        alpha = alpha)
    
    tpp2dObj <- new("tpp2dExperiment",
                    configTable = cfgTab,
                    idVar = idVar,
                    intensityStr = intensityStr,
                    fcStr = fcStr,
                    nonZeroCols = nonZeroCols,
                    geneNameVar = geneNameVar,
                    qualColName = qualColName,
                    naStrs = naStrs,
                    concFactor = concFactor,
                    medianNormalizeFC = medianNormalizeFC,
                    filterContaminants = filterContaminants,
                    minObs = minObs,
                    independentFiltering = independentFiltering,
                    fcThres = fcThres,
                    optim_fun_h0 = optim_fun_h0,
                    optim_fun_h1 = optim_fun_h1,
                    slopEC50 = slopEC50,
                    maxit = maxit,
                    BPPARAM = class(BPPARAM)[1],
                    B = B,
                    byMsExp = byMsExp,
                    alpha = alpha,
                    tidyDataTable = df,
                    modelParamsDf = params_df,
                    resultTable = fstat_df,
                    bootstrapNullDf = null_df,
                    hitTable = hits_df)
    return(tpp2dObj)
}