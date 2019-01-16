#' @import TPP
import2dMain <- function(configTable, data, idVar, fcStr,
                         addCol, naStrs, intensityStr,
                         qualColName, nonZeroCols){
  
  # internal import main function, adapted from TPP package
  files <- configTable$Path
  if (!is.null(files)) {
    if (any(files == "")) {
      files <- NULL
    }
  }
  Experiment = Compound = Temperature = RefCol <- NULL
  expNames <- configTable$Experiment
  argList <- TPP:::importFct_CheckDataFormat(dataframes = data, files = files,
                                             expNames = expNames)
  data <- argList[["dataframes"]]
  files <- argList[["files"]]
  if (!is.null(files)) {
    files2 <- files[!duplicated(names(files))]
    data <- TPP:::importFct_readFiles(files = files2, naStrs = naStrs)
  }
  configTable %>% group_by(Experiment, Compound,
                           Temperature, RefCol)
  iVec <- 1:nrow(configTable)
  dataList <- lapply(iVec, function(iTmp) {
    rowTmp <- configTable[iTmp, ]
    expTmp <- rowTmp$Experiment
    message("Importing 2D-TPP dataset: ", expTmp)
    tTmp <- rowTmp$Temperature
    dataTmp <- data[[expTmp]]
    noFCCols <- c("Compound", "Experiment", "Temperature",
                  "RefCol", "Path", "Condition")
    allCols <- colnames(rowTmp)
    labelCols <- setdiff(allCols, noFCCols)
    labelValues <- suppressMessages(rowTmp[, labelCols] %>%
                                      as.numeric)
    labelColsNum <- labelCols[!is.na(labelValues)]
    signalCols <- paste(intensityStr, labelColsNum, sep = "")
    relevant.cols <- c(idVar, qualColName, nonZeroCols, addCol,
                       signalCols) %>% unique
    if (!is.null(fcStr)) {
      fcCols <- paste(fcStr, labelColsNum, sep = "")
      relevant.cols <- c(relevant.cols, fcCols)
      dataCols <- fcCols
    }
    else {
      dataCols <- signalCols
    }
    if (!all(relevant.cols %in% colnames(dataTmp))) {
      notFound <- paste(setdiff(relevant.cols, colnames(dataTmp)),
                        collapse = "', '")
      stop("The following columns could not be found: '",
           notFound, "'. Please check the suffices and the additional column names you have specified.")
    }
    dataFiltered <- TPP:::importFct_removeDuplicates(
      inDF = dataTmp,refColName = idVar, nonNAColNames = dataCols, qualColName = qualColName[1])
    idsTmp <- as.character(dataFiltered[, idVar])
    idsAnnotated <- paste(expTmp, tTmp, idsTmp, sep = "_")
    dataFinal <- dataFiltered %>% subset(select = relevant.cols) %>%
      mutate(temperature = tTmp, experiment = expTmp, unique_ID = idsAnnotated)
    return(dataFinal)
  })
  newNames <- sapply(seq(nrow(configTable)), function(iTmp) {
    rowTmp <- configTable[iTmp, ]
    tTmp <- rowTmp$Temperature
    expTmp <- rowTmp$Experiment
    newName <- paste(expTmp, tTmp, sep = "_")
    return(newName)
  })
  names(dataList) <- newNames
  out <- TPP:::importFct_rmZeroSias(configTable = configTable, data.list = dataList,
                                    intensityStr = intensityStr)
  return(out)
}

configWide2Long <- function(configWide){
  # internal function to tranform config table into long format
  if(any(grepl("Path", colnames(configWide)))){
    configLong <- configWide %>%
      dplyr::select(-Path) %>%
      gather(label, conc, -Compound, -Experiment, -Temperature, -RefCol) %>%
      filter(conc != "-")
  }else{
    configLong <- configWide %>%
      gather(label, conc, -Compound, -Experiment, -Temperature, -RefCol) %>%
      filter(conc != "-")
  }
}

annotateDataList <- function(dataList, geneNameVar, configLong,
                             intensityStr, fcStr){
  # internal function to annotate list of 2D-TPP data subtables with
  # information from config table
  combinedTab <- bind_rows(lapply(dataList, function(dat){
    datLong <- dat %>% tbl_df() %>%
      gather(channel, signal, matches(intensityStr), matches(fcStr)) %>%
      mutate(label = gsub(fcStr, "", gsub(intensityStr, "", channel))) %>%
      left_join(configLong %>% dplyr::select(Temperature, RefCol, label, conc),
                by = c("temperature" = "Temperature", "label")) %>%
      mutate(var = ifelse(grepl(fcStr, channel), "rel_value", "raw_value")) %>%
      dplyr::select(-channel, -unique_ID) %>%
      spread(var, signal)
  }))
  return(combinedTab)
}

filterOutContaminants <- function(dataLong){
  # internal function to filter out contaminants
  filter(dataLong, !grepl("##", representative))
}

checkRatioRef <- function(dataLong, idVar, concFactor = 1e6){
  # internal function to check that protein fold changes are computed
  # relative to the correct TMT channel
  if(!all(filter(dataLong, label == RefCol)$rel_value == 1, na.rm = TRUE)){
    message("Recomputing ratios!")
    dataOut <- dataLong %>%
      dplyr::group_by_(idVar, "temperature") %>%
      mutate(rel_value = rel_value/rel_value[label == RefCol]) %>%
      ungroup %>%
      filter(!is.na(raw_value)) %>%
      mutate(conc = as.numeric(conc)) %>%
      mutate(log_conc = log10(conc/concFactor))
    
    return(dataOut)
    
  }else{
    message("Ratios were correctly computed!")
    return(dataLong %>%
             filter(!is.na(raw_value)) %>%
             mutate(conc = as.numeric(conc)) %>%
             mutate(log_conc = log10(conc/concFactor)))
  }
}

medianNormalizeRatios <- function(dataLong){
  # internal function to perform median normalization of ratios
  dataOut <- dataLong %>%
    rename(raw_rel_value = rel_value) %>%
    group_by(temperature, conc) %>%
    mutate(rel_value = raw_rel_value / median(raw_rel_value, na.rm = TRUE)) %>%
    ungroup()
  
  return(dataOut)
}

renameColumns <- function(dataLong, idVar, geneNameVar){
  # internal function to rename column names to match lazyeval variable
  # names of main function
  dplyr::rename_(dataLong, "representative" = idVar, 
                 "clustername" = geneNameVar) %>%
    group_by(clustername) %>%
    mutate(representative =
             paste_rmNA(unique(unlist(strsplit(representative, split = "\\|"))), 
                        sep = "|")) %>%
    ungroup()
}

#' Import 2D-TPP dataset using a config table
#' 
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
#' @param filterContaminants boolean variable indicating whether data 
#' should be filtered to exclude contaminants (default: TRUE).
#' @param nonZeroCols column like default qssm that should be imported and
#' requested to be non-zero in analyzed data
#' @param geneNameVar character string of the column name that describes
#' the gene name of a given protein in the raw data files
#' @param concFactor numeric value that indicates how concentrations need to 
#' be adjusted to yield total unit e.g. default mmol - 1e6
#' 
#' @export
#'
#' @import TPP
import2dDataset <- function(configTable, data,
                            idVar = "representative",
                            intensityStr = "sumionarea_protein_",
                            fcStr = "rel_fc_protein_",
                            nonZeroCols = "qssm",
                            geneNameVar = "clustername",
                            addCol = "",
                            qualColName = "qupm",
                            naStrs = c("NA", "n/d", "NaN"),
                            concFactor = 1e6,
                            medianNormalizeFC = TRUE,
                            filterContaminants = TRUE){
  
  configWide <- TPP:::importCheckConfigTable(infoTable = configTable, type = "2D")
  configLong <- configWide2Long(configWide = configWide)
  
  dataList <- import2dMain(configTable = configWide,
                           data = data,
                           idVar = idVar,
                           fcStr = fcStr,
                           addCol = c(geneNameVar, addCol),
                           naStrs = naStrs,
                           intensityStr = intensityStr,
                           nonZeroCols = nonZeroCols,
                           qualColName = qualColName)
  
  dataLong <- annotateDataList(dataList = dataList,
                               geneNameVar = geneNameVar,
                               configLong = configLong,
                               intensityStr = intensityStr,
                               fcStr = fcStr)
  
  dataRatioChecked <- checkRatioRef(dataLong, idVar = idVar,
                                    concFactor = concFactor)
  
  if(medianNormalizeFC){
    message("Median normalizing fold changes...")
    dataNorm <- medianNormalizeRatios(dataRatioChecked)
  }else{
    dataNorm <- dataRatioChecked
  }
  
  dataOut <- renameColumns(dataNorm,
                           idVar = idVar,
                           geneNameVar = geneNameVar)
  
  if(filterContaminants){
    dataOut <- filterOutContaminants(dataOut)
  }
  
  return(dataOut)
}