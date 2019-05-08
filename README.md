# TPP2D: FDR-controlled analysis of 2D-TPP datasets with R

This package contains functions to analyze 2D-thermal proteome profiles using a functional analysis approach. This is done by fitting two competing models (H0 and H1) to the thermal profile of each protein and asking whether the H1 model explains the variance in the data significantly better than H0.

## Installation

### Installation from Bioconductor (recommended)

The `TPP2D` package has been accepted into Bioconductor! Therefore, the most reliable way to install it is via BiocManager:

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("TPP2D")
```

### Installation of the development version

```{R}
BiocManager::install(“nkurzaw/TPP2D”)
```

## Getting started

The easiest way of learning how to use the `TPP2D` package is to browse it's vignette:
```{R}
library(TPP2D)
browseVignettes("TPP2D")
```