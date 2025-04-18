% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdSuSiE.R
\name{sdSuSiE_Plot}
\alias{sdSuSiE_Plot}
\title{The function creates GWAS plots from sex-stratified GWASs and univariate sex-dimorphic analysis, and PIP plots from sdSuSiE.}
\usage{
sdSuSiE_Plot(res, R_mat, summary_data)
}
\arguments{
\item{res}{The result object from sdSuSiE.}

\item{R_mat}{A list of length 2 containing LD matrices for each sex.}

\item{summary_data}{A list of data frames for each sex containing at least SNP, Z, and possibly POS (position) columns.}
}
\value{
A combined ggplot object visualizing the results of sex-stratified GWASs, univariate sex-dimorphic analysis and sdSuSiE.
}
\description{
The function creates GWAS plots from sex-stratified GWASs and univariate sex-dimorphic analysis, and PIP plots from sdSuSiE.
}
