% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_analysis.R
\name{fitAndEvalDataset}
\alias{fitAndEvalDataset}
\title{Fit H0 and H1 model to 2D thermal profiles of proteins
and compute F statistic}
\usage{
fitAndEvalDataset(
  df,
  maxit = 500,
  optim_fun_h0 = .min_RSS_h0,
  optim_fun_h1 = .min_RSS_h1_slope_pEC50,
  optim_fun_h1_2 = NULL,
  gr_fun_h0 = NULL,
  gr_fun_h1 = NULL,
  gr_fun_h1_2 = NULL,
  ec50_lower_limit = NULL,
  ec50_upper_limit = NULL,
  slopEC50 = TRUE
)
}
\arguments{
\item{df}{tidy data_frame retrieved after import of a 2D-TPP 
dataset, potential filtering and addition of a column "nObs"
containing the number of observations per protein}

\item{maxit}{maximal number of iterations the optimization
should be given, default is set to 500}

\item{optim_fun_h0}{optimization function that should be used
for fitting the H0 model}

\item{optim_fun_h1}{optimization function that should be used
for fitting the H1 model}

\item{optim_fun_h1_2}{optional additional optimization function 
that will be run with paramters retrieved from optim_fun_h1 and 
should be used for fitting the H1 model with the trimmed sum
model, default is NULL}

\item{gr_fun_h0}{optional gradient function for optim_fun_h0,
default is NULL}

\item{gr_fun_h1}{optional gradient function for optim_fun_h1,
default is NULL}

\item{gr_fun_h1_2}{optional gradient function for optim_fun_h1_2,
default is NULL}

\item{ec50_lower_limit}{lower limit of ec50 parameter}

\item{ec50_upper_limit}{lower limit of ec50 parameter}

\item{slopEC50}{logical flag indicating whether the h1 model is
fitted with a linear model describing the shift od the pEC50 over 
temperatures}
}
\value{
data frame with H0 and H1 model characteristics for each
protein and respectively computed F statistics
}
\description{
Fit H0 and H1 model to 2D thermal profiles of proteins
and compute F statistic
}
\examples{
data("simulated_cell_extract_df")
temp_df <- simulated_cell_extract_df \%>\% 
  group_by(representative) \%>\% 
  mutate(nObs = n()) \%>\% 
  ungroup 
fitAndEvalDataset(temp_df)  

}
