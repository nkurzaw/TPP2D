% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_funcs.R
\name{renameColumns}
\alias{renameColumns}
\title{Rename columns of imported data frame}
\usage{
renameColumns(dataLong, idVar, geneNameVar)
}
\arguments{
\item{dataLong}{long format data frame of imported dataset}

\item{idVar}{character string indicating which data column provides the 
unique identifiers for each protein.}

\item{geneNameVar}{character string of the column name that describes
the gene name of a given protein in the raw data files}
}
\value{
data frame containing imported data with renamed
columns
}
\description{
Rename columns of imported data frame
}
\examples{
data("config_tab")
data("raw_dat_list")

dataList <- import2dMain(configTable = config_tab,
                         data = raw_dat_list,
                         idVar = "protein_id",
                         fcStr = "rel_fc_",
                         addCol = "gene_name",
                         naStrs = NA,
                         intensityStr = "signal_sum_",
                         nonZeroCols = "qusm",
                         qualColName = "qupm")
configLong <- configWide2Long(configWide = config_tab)
annoDat <- annotateDataList(dataList = dataList,
                            geneNameVar = "gene_name",
                            configLong = configLong,
                            intensityStr = "signal_sum_",
                            fcStr = "rel_fc_")
renameColumns(annoDat, 
              idVar = "protein_id", 
              geneNameVar = "gene_name")
}
