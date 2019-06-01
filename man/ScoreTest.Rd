\name{ScoreTest}
\alias{ScoreTest}
\title{Vector of score test statistics between a gene and a binary phenotype}
\usage{
ScoreTest(X, Y, U = NULL, Y0 = NULL, N = 1000)
}
\arguments{
\item{X}{a n x p matrix of genotypes (coded as 0/1/2)}

\item{Y}{a vector of length n containing the phenotype values (0/1)}

\item{U}{an optional n x q matrix of covariates}

\item{Y0}{an optional matrix of permuted phenotypes}

\item{N}{the desired number of permutations of the phenotype.}
}
\value{
A list with 3 elements:
\item{Z}{the vector of test statistics}
\item{Z0}{the N x p matrix of vector computed on permuted phenotypes}
\item{Sigma}{the correlation matrix corresponding to Z}
}
\description{
\code{ScoreTest} is used to compute the vector of score test statistics between a gene and a binary 
phenotype.
}
