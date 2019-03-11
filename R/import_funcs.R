TPP_importFct_CheckDataFormat <- function (files, dataframes, expNames){
    # internal function copied from TPP package to avoid 
    # import of non-exported package functions
    . <- NULL
    isDF <- !is.null(dataframes)
    isF <- !is.null(files)
    isBoth <- isDF & isF
    isNone <- !(isDF | isF)
    if (isBoth) {
        stop("Data import function received a",
             " filename AND a dataframe object. \n",
             "Please specify only one.")
    }
    else if (isNone) {
        stop("Data import function requires a", 
             " filename or a dataframe object. \n",
             "Please specify one.")
    }
    if (isDF) {
        isClassList <- is.list(dataframes) && !is.data.frame(dataframes)
        isClassDF <- is.data.frame(dataframes)
        if (isClassList) {
            classesInList <- dataframes %>% 
            vapply(. %>% inherits(., "data.frame"), TRUE)
            if (!all(classesInList)) {
                stop(paste("Argument 'dataframes' contains", 
                           "elements that are not of type", 
                           "'data.frame' at the following positions: "), 
                     which(!classesInList) %>% paste(collapse = ", "), ".")
            }
          }
          else if (isClassDF) {
              dataframes <- list(dataframes)
              names(dataframes) <- expNames
          }
          else {
              stop("Argument 'dataframes' must be either an object of class \n
                   'data.frame', or a list of such objects.")
          }
    }
    if (isF) {
        files <- as.character(files)
        names(files) <- expNames
    }
    return(list(files = files, dataframes = dataframes))
}

#' @importFrom utils read.delim
#' @importFrom RCurl url.exists
TPP_importFct_readFiles <- function (files, naStrs){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  expNames <- names(files)
  data <- vector("list", length(files))
  names(data) <- expNames
  for (expName in expNames) {
    fTmp <- files[[expName]]
    if (file.exists(fTmp) || url.exists(fTmp)) {
      data[[expName]] <- read.delim(fTmp, as.is = TRUE, 
                                    na.strings = naStrs, quote = "")
    }
    else {
      stop("File ", fTmp, " could not be found.")
    }
  }
  return(data)
}

TPP_importFct_removeDuplicates <- function(inDF, refColName, 
                                           nonNAColNames, qualColName){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  message("Removing duplicate identifiers using quality column '", 
          qualColName, "'...")
  nonUniques = unique(inDF[duplicated(inDF[[refColName]]), 
                           refColName])
  retDF = subset(inDF, !(get(refColName) %in% nonUniques))
  for (nU in nonUniques) {
    tmpDF = subset(inDF, get(refColName) == nU)
    nonNArows = NULL
    for (r in seq_len(nrow(tmpDF))) {
      if (any(!is.na(tmpDF[r, nonNAColNames]))) {
        nonNArows = c(nonNArows, r)
      }
    }
    if (length(nonNArows) > 1) {
      if (is.null(qualColName)) {
        useRow = 1
      }
      else {
        qualVals = tmpDF[nonNArows, qualColName]
        useRow = match(max(qualVals), qualVals)
      }
    }
    else {
      useRow = nonNArows[1]
    }
    retDF = rbind(retDF, tmpDF[useRow, ])
  }
  message(nrow(retDF), " out of ", nrow(inDF), 
          " rows kept for further analysis.")
  return(retDF)
}

TPP_replaceZeros <- function(x){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  x[which(x == 0)] <- NA
  return(x)
}

TPP_importFct_rmZeroSias <- function(configTable, data.list, 
                                     intensityStr){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  out <- lapply(names(data.list), function(l.name) {
    datTmp <- data.list[[l.name]]
    colsTmp <- colnames(datTmp)
    intensity.cols <- grep(intensityStr, colsTmp, value = TRUE)
    intensity.df <- subset(datTmp, select = intensity.cols) %>% 
      mutate_all(as.character) %>% mutate_all(as.numeric)
    new.intensity.df <- intensity.df %>% mutate_all(TPP_replaceZeros)
    datTmp[, intensity.cols] <- new.intensity.df
    return(datTmp)
  })
  names(out) <- names(data.list)
  return(out)
}

TPP_importFct_checkExperimentCol <- function(expCol){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  if (is.null(expCol)) {
    m <- paste("Config table needs an 'Experiment'", 
               "column with unique experiment IDs.")
    stop(m, "\n")
  }
  oldExpNames <- expCol
  newExpNames <- gsub("([^[:alnum:]])", "_", expCol)
  iChanged <- oldExpNames != newExpNames
  if (any(iChanged)) {
    m1 <- paste("Replaced non-alphanumeric characters", 
                "in the 'Experiment' column entries:")
    m2 <- paste("'", paste(oldExpNames[iChanged], collapse = "', '"), 
                "'\nby\n'", paste(newExpNames[iChanged], collapse = "', '"), 
                sep = "")
    message(m1, "\n", m2, "\n")
  }
  return(newExpNames)
}

TPP_importFct_checkComparisons <- function(confgTable){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  expConds <- confgTable$Condition
  expNames <- confgTable$Experiment
  compCols <- grep("Comparison", colnames(confgTable), ignore.case = TRUE, 
                   value = TRUE)
  compChars <- apply(confgTable[compCols], 2, function(x) {
    length(grep("[[:alnum:]]", x, value = TRUE))
  })
  comp_unequal_two <- compChars != 2
  if (any(comp_unequal_two)) {
    warning(paste("\nThe following comparison columns could not be evaluated", 
                  "because they did not contain exactly two entries:\n\t\t"), 
            paste(compCols[comp_unequal_two], collapse = ",\n\t\t"))
  }
  validCompCols <- compCols[!comp_unequal_two]
  allCompStrs <- c()
  if (length(validCompCols) > 0) {
    message("Comparisons will be performed between the following experiments:")
    for (colName in validCompCols) {
      current_compEntries <- confgTable[[colName]]
      current_compRows <- grep("[[:alnum:]]", current_compEntries)
      current_compExps <- expNames[current_compRows]
      compRef <- current_compExps[1]
      compTreatm <- current_compExps[2]
      if ("Condition" %in% names(confgTable)) {
        current_compConds <- expConds[current_compRows]
        if ("Vehicle" %in% current_compConds && "Treatment" %in% 
            current_compConds) {
          compRef <- current_compExps[current_compConds == 
                                        "Vehicle"]
          compTreatm <- current_compExps[current_compConds == 
                                           "Treatment"]
        }
      }
      compStr <- paste(compTreatm, "_vs_", compRef, sep = "")
      names(compStr) <- colName
      message(compStr)
      allCompStrs <- c(allCompStrs, compStr)
    }
    message("\n")
  }
  return(allCompStrs)
}

#' @importFrom stringr str_to_title
TPP_importFct_checkConditions <- function(condInfo, 
                                          expectedLength){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  flagGenerateConds <- FALSE
  if (is.null(condInfo)) {
    message("No information about experimental conditions given.", 
            "Assigning NA instead.\n",
            "Reminder: recognition of Vehicle and Treatment groups", 
            "during pairwise \n",
            "comparisons is only possible when they are specified ",
            "in the config table.\n")
    condInfo <- rep(NA_character_, expectedLength)
  }
  else {
    condInfo <- as.character(condInfo) %>% 
      stringr::str_to_title()
    condLevels <- unique(condInfo)
    invalidLevels = 
      setdiff(condLevels, c("Treatment", "Vehicle"))
    if (length(invalidLevels) > 0) {
      stop("The entry '", invalidLevels, 
           paste("' in the condition column is invalid.", 
                 "Only the values 'Treatment' and", 
                 "'Vehicle' are allowed. Please correct", 
                 "this and start again."))
    }
  }
  return(condInfo)
}

TPP_checkFunctionArgs <- 
  function(functionCall, expectedArguments){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  myArgs <- names(functionCall)
  lapply(expectedArguments, function(arg) {
    if (!arg %in% myArgs) {
      stop("Error in ", paste(functionCall)[1], 
           ": argument '", 
           arg, "' is missing, with no default", 
           call. = FALSE)
    }
  })
}

TPP_nonLabelColumns <- function(){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  out <- data.frame(
    column = c("Experiment", "Experiment",
               "Experiment", "Path", "Path", 
               "Path", "Condition", "Replicate", 
               "Compound", "Temperature", "RefCol"), 
    type = c("TR", "CCR", "2D", "TR", "CCR", "2D", 
             "TR", "TR", "2D", "2D", "2D"), 
    obligatory = c(TRUE, TRUE, TRUE, FALSE, FALSE, 
                   FALSE, TRUE, FALSE, TRUE, TRUE, TRUE), 
    exclusive = c(FALSE, FALSE, FALSE, FALSE, FALSE, 
                  FALSE, TRUE, TRUE, TRUE, TRUE, TRUE))
  return(out)
}

TPP_detectLabelColumnsInConfigTable <- 
  function(allColumns){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  TPP_checkFunctionArgs(match.call(), c("allColumns"))
  noLabelCols <- TPP_nonLabelColumns()$column %>% 
    as.character %>% 
    unique
  compCols <- grep("comparison", allColumns, value = TRUE, 
                   ignore.case = TRUE)
  noLabelCols <- c(noLabelCols, compCols)
  labelCols <- setdiff(allColumns, noLabelCols)
  return(labelCols)
}

TPP_importCheckTemperatures <- function(temp){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  tempMatrix <- as.matrix(temp)
  rownames(tempMatrix) <- NULL
  naRows <- apply(is.na(tempMatrix), 1, all)
  if (any(naRows)) {
    stop("Row(s) ", paste(which(naRows), collapse = ", "), 
         " in the configuration table contain", 
         " only missing temperature values.")
  }
  return(tempMatrix)
}

#' @importFrom openxlsx read.xlsx
#' @importFrom utils read.table
TPP_importFct_readConfigTable <- function(cfg){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  if (is.character(cfg)) {
    if (file.exists(cfg)) {
      strChunks <- strsplit(cfg, "\\.")[[1]]
      fileExtension <- strChunks[length(strChunks)]
      if (fileExtension == "txt") {
        tab <- read.table(
          file = cfg, header = TRUE, 
          check.names = FALSE, stringsAsFactors = FALSE, 
          sep = "\t")
      }
      else if (fileExtension == "csv") {
        tab <- read.table(
          file = cfg, header = TRUE, 
          check.names = FALSE, stringsAsFactors = FALSE, 
          sep = ",")
      }
      else if (fileExtension == "xlsx") {
        tab <- openxlsx::read.xlsx(cfg)
      }
      else {
        stop("Error during data import: ", cfg, 
             " does not belong to a valid configuration file.")
      }
    }
    else {
      stop("Error during data import: ", cfg, 
           " does not belong to a valid configuration file.")
    }
    cfg <- tab
  }
  return(cfg)
}

TPP_importCheckConfigTable <- function (infoTable, type = "2D"){
  # internal function copied from TPP package to avoid 
  # import of non-exported package functions
  TPP_checkFunctionArgs(match.call(), c("infoTable", "type"))
  Experiment = Path = Compound <- NULL
  isValidDF <- FALSE
  if (is.data.frame(infoTable)) {
    if ((ncol(infoTable) > 1) & 
        ("Experiment" %in% colnames(infoTable))) {
      isValidDF <- TRUE
    }
  }
  if (!is.character(infoTable) & !isValidDF) {
    stop("'infoTable' must either be a data frame", 
         " with an 'Experiment' column \n",
         "and at least one isobaric label column,", 
         "or a filename pointing at a \n",
         "table that fulfills the same criteria")
  }
  isValidType <- type %in% c("2D")
  if (!isValidType) {
    stop("'type' must have this value: '2D'")
  }
  infoTable <- TPP_importFct_readConfigTable(cfg = infoTable)
  infoTable$Experiment <- 
    TPP_importFct_checkExperimentCol(infoTable$Experiment)
  infoTable <- subset(infoTable, Experiment != "")
  givenPaths <- NULL
  if (any("Path" %in% colnames(infoTable))) {
    if (all(infoTable$Path == "") || all(is.na(infoTable$Path))) {
      message("Removing empty 'Path' column from config table")
      infoTable <- infoTable %>% select(-Path)
    }
    else {
      givenPaths <- infoTable$Path
    }
  }
  compStrs <- NA
  infoTable$Condition <- NULL
  allCols <- colnames(infoTable)
  labelCols <- TPP_detectLabelColumnsInConfigTable(allColumns = allCols)
  labelValues <- infoTable[, labelCols]
  labelValuesNum <- suppressWarnings(labelValues %>% apply(2, 
                                                           as.numeric))
  if (is.matrix(labelValuesNum)) {
    isInvalid <- labelValuesNum %>% apply(2, is.na) %>% apply(2, 
                                                              all)
  }
  else if (is.vector(labelValuesNum)) {
    isInvalid <- is.na(labelValuesNum)
  }
  invalidLabels <- labelCols[isInvalid]
  infoTable[, invalidLabels] <- NULL
  labelColsNew <- labelCols[!isInvalid]
  labelStr <- paste(labelColsNew, collapse = ", ")
  message("The following valid label columns were detected:\n", 
          labelStr, ".")
  if (type == "2D") {
    temperatures <- infoTable$Temperature
    if (is.null(temperatures) | length(temperatures) < 2) {
      m1 <- paste("Insufficient temperatures (<2)", 
                  "specified in config file.")
      m2 <- paste("Does your configuration table", 
                  "have the correct column names?")
      stop(m1, "\n", m2)
    }
    else if (length(which(!infoTable$RefCol %in% labelColsNew)) != 
             0) {
      stop("Labels in reference column not found", 
           "in any of teh label columns.")
    }
    hasCompoundCol <- any(allCols == "Compound")
    if (!hasCompoundCol) {
      m <- paste("Config table of a 2D-TPP experiment", 
                 "needs a 'Compound' column.")
      stop(m, "\n")
    }
    else {
      infoTable <- infoTable %>% 
        mutate(Compound = 
                 gsub("([^[:alnum:]])", "_", Compound))
    }
    out <- infoTable
  }
  else {
    temperatures <- subset(infoTable, select = labelColsNew)
    tempMatrix <- TPP_importCheckTemperatures(temp = temperatures)
    infoList <- list(
      expNames = as.character(infoTable$Experiment), 
      expCond = infoTable$Condition, files = givenPaths, 
      compStrs = compStrs, labels = labelColsNew, 
      tempMatrix = tempMatrix)
    out <- infoList
  }
  return(out)
}

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
  argList <- TPP_importFct_CheckDataFormat(dataframes = data, 
                                           files = files,
                                           expNames = expNames)
  data <- argList[["dataframes"]]
  files <- argList[["files"]]
  if (!is.null(files)) {
    files2 <- files[!duplicated(names(files))]
    data <- TPP_importFct_readFiles(files = files2, 
                                    naStrs = naStrs)
  }
  configTable %>% group_by(Experiment, Compound,
                           Temperature, RefCol)
  iVec <- seq_len(nrow(configTable))
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
           notFound, paste("'. Please check the suffices and the", 
                           "additional column names you have specified."))
    }
    dataFiltered <- TPP_importFct_removeDuplicates(
      inDF = dataTmp,refColName = idVar, 
      nonNAColNames = dataCols, 
      qualColName = qualColName[1])
    idsTmp <- as.character(dataFiltered[, idVar])
    idsAnnotated <- paste(expTmp, tTmp, idsTmp, sep = "_")
    dataFinal <- dataFiltered %>% subset(select = relevant.cols) %>%
      mutate(temperature = tTmp, experiment = expTmp, unique_ID = idsAnnotated)
    return(dataFinal)
  })
  newNames <- vapply(seq(nrow(configTable)), function(iTmp) {
    rowTmp <- configTable[iTmp, ]
    tTmp <- rowTmp$Temperature
    expTmp <- rowTmp$Experiment
    newName <- paste(expTmp, tTmp, sep = "_")
    return(newName)
  }, "")
  names(dataList) <- newNames
  out <- TPP_importFct_rmZeroSias(configTable = configTable, 
                                  data.list = dataList,
                                  intensityStr = intensityStr)
  return(out)
}

#' @importFrom tidyr gather
configWide2Long <- function(configWide){
  # internal function to tranform config table into long format
  
  Path <- label <- conc <- Compound <- Experiment <- 
    Temperature <- RefCol <- NULL
  
  if(any(grepl("Path", colnames(configWide)))){
    configLong <- configWide %>%
      dplyr::select(-Path) %>%
      gather(label, conc, -Compound, 
             -Experiment, -Temperature, -RefCol) %>%
      filter(conc != "-")
  }else{
    configLong <- configWide %>%
      gather(label, conc, -Compound, 
             -Experiment, -Temperature, -RefCol) %>%
      filter(conc != "-")
  }
}

#' @importFrom tidyr spread
annotateDataList <- function(dataList, geneNameVar, configLong,
                             intensityStr, fcStr){
  # internal function to annotate list of 2D-TPP data subtables with
  # information from config table
  channel <- signal <- Temperature <- RefCol <- label <- 
    conc <- unique_ID <- spread_var <- NULL
  
  combinedTab <- bind_rows(lapply(dataList, function(dat){
    datLong <- dat %>% tbl_df() %>%
      gather(channel, signal, matches(intensityStr), matches(fcStr)) %>%
      mutate(label = gsub(fcStr, "", gsub(intensityStr, "", channel))) %>%
      left_join(configLong %>% 
                  dplyr::select(Temperature, RefCol, label, conc),
                by = c("temperature" = "Temperature", "label")) %>%
      mutate(spread_var = 
               ifelse(grepl(fcStr, channel), "rel_value", "raw_value")) %>%
      dplyr::select(-channel, -unique_ID) %>%
      spread(spread_var, signal)
  }))
  return(combinedTab)
}

filterOutContaminants <- function(dataLong){
  # internal function to filter out contaminants
  representative <- NULL
  filter(dataLong, !grepl("##", representative))
}

checkRatioRef <- function(dataLong, idVar, concFactor = 1e6){
  # internal function to check that protein 
  # fold changes are computed
  # relative to the correct TMT channel
  label <- RefCol <- rel_value <- raw_value <- conc <- NULL
  
  if(!all(filter(dataLong, label == RefCol)$rel_value == 1, 
          na.rm = TRUE)){
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

#' @importFrom stats median
medianNormalizeRatios <- function(dataLong){
  # internal function to perform median normalization 
  # of ratios
  rel_value <- temperature <- conc <- 
    raw_rel_value <- NULL
  
  dataOut <- dataLong %>%
    rename(raw_rel_value = rel_value) %>%
    group_by(temperature, conc) %>%
    mutate(rel_value = raw_rel_value / 
             median(raw_rel_value, na.rm = TRUE)) %>%
    ungroup()
  
  return(dataOut)
}

renameColumns <- function(dataLong, idVar, geneNameVar){
  # internal function to rename column names to 
  # match lazyeval variable
  # names of main function
  clustername <- representative <- NULL
  
  dplyr::rename_(dataLong, "representative" = idVar, 
                 "clustername" = geneNameVar) %>%
    group_by(clustername) %>%
    mutate(representative =
             paste_rmNA(unique(unlist(strsplit(representative, 
                                               split = "\\|"))), 
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
#' @return tidy data frame representing a 2D-TPP dataset
#' 
#' @examples 
#' data("config_tab")
#' data("raw_dat_list")
#' import_df <- import2dDataset(configTable = config_tab, 
#'                              data = raw_dat_list,
#'                              idVar = "protein_id",
#'                              intensityStr = "signal_sum_",
#'                              fcStr = "rel_fc_",
#'                              nonZeroCols = "qusm",
#'                              geneNameVar = "gene_name",
#'                              addCol = NULL,
#'                              qualColName = "qupm",
#'                              naStrs = c("NA", "n/d", "NaN"),
#'                              concFactor = 1e6,
#'                              medianNormalizeFC = TRUE,
#'                              filterContaminants = TRUE)
#' 
#' @export
import2dDataset <- function(configTable, data,
                            idVar = "representative",
                            intensityStr = "sumionarea_protein_",
                            fcStr = "rel_fc_protein_",
                            nonZeroCols = "qssm",
                            geneNameVar = "clustername",
                            addCol = NULL,
                            qualColName = "qupm",
                            naStrs = c("NA", "n/d", "NaN"),
                            concFactor = 1e6,
                            medianNormalizeFC = TRUE,
                            filterContaminants = TRUE){
  
  configWide <- TPP_importCheckConfigTable(
    infoTable = configTable, type = "2D")
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