% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdSuSiE.R
\name{compute_softmax}
\alias{compute_softmax}
\title{The function computes the posterior inclusion probability for the l'th effect using log Bayes Factors.}
\usage{
compute_softmax(lbf, prior_weights)
}
\arguments{
\item{lbf}{Log Bayes Factors.}

\item{prior_weights}{A vector of prior weights.}
}
\value{
A list containing the posterior inclusion probability (alpha_wmulti) and the logarithm of the likelihood (loglik).
}
\description{
The function computes the posterior inclusion probability for the l'th effect using log Bayes Factors.
}
