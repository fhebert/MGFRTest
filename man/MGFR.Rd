\name{MGFR}
\alias{MGFR}
\title{Moment-Generating Function Ratio Test}
\usage{
MGFR(Z, Z0, Sigma = NULL, eigSigma = NULL, vtt = NULL)
}
\arguments{
\item{Z}{the vector of score test statistics between the gene and the phenotype}

\item{Z0}{a matrix of score test statistics between the gene and permuted versions of the phenotype}

\item{Sigma}{the correlation matrix corresponding to \code{Z} (optional)}

\item{eigSigma}{the eigendecomposition of \code{Sigma} (optional). One of \code{Sigma} or \code{eigSigma} must be given}

\item{vtt}{an optional sequence of values for the parameter \code{t}}
}
\value{
The p-value of the MGF-R test.
}
\description{
\code{MGFRTest} is used to compute the p-value of the MGF-R test to test the presence of 
an association between a gene and a binary phenotype.
}
