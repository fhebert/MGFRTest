\name{MGFRTest-package}
\alias{MGFRTest-package}
\docType{package}
\encoding{UTF-8}
\title{
  Moment Generating Function Ratio Test
}
\description{
The MGFRTest package provides functions to compute the global p-value of the MGF-R Test to test the presence
of an association between a gene and a binary phenotype. More generally, it can be used as a signal detection method for testing the global nullity of the mean vector of a multivariate normal distribution.
}
\author{
Florian Hébert, Mathieu Emily, David Causeur

Maintainer: Florian Hébert <florian.hebert@agrocampus-ouest.fr>
}
\references{Hébert, F., Emily, M., and Causeur, D. (2019). An Adaptive Decorrelation Procedure for Signal Detection. Submitted.}
\keyword{Genome-wide association studies, MGF-R Test}
\seealso{
  \code{\link{ScoreTest}}
  \code{\link{MGFR}}
}
\examples{
  #Matrix of 10 independent SNPs and 2000 individuals
  X = matrix(rbinom(10*2000,2,0.4),ncol=10)  
  
  #Phenotype: 1000 cases, 1000 controls
  Y = sample(c(rep(0,1000),rep(1,1000)))    
  
  #Compute the vector of test statistics, the permuted test statistics and the correlation matrix
  Z = ScoreTest(X,Y)
  
  #Compute the p-value of the MGF-R test
  MGFR(Z$Z,Z$Z0,Z$Sigma)
}