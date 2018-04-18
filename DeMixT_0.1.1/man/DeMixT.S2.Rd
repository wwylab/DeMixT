\name{DeMixT.S2}
\alias{DeMixT.S2}

\title{
Deconvolves expressions of each individual sample for unknown component
}
\description{
This function is designed to estimate the deconvolved expressions of individual mixed tumor samples for unknown component for each gene.
}
\usage{
DeMixT.S2(data.Y, data.comp1, data.comp2 = NULL, givenpi, nbin = 50, nthread = detectCores() - 1)
}

\arguments{
  \item{data.Y}{A matrix of expression data from mixed tumor samples. It is a \eqn{G} by \eqn{Sy} matrix where \eqn{G} is the number of genes and \eqn{Sy} is the number of mixed samples. Samples with the same tissue type should be placed together in columns.}

  \item{data.comp1}{A matrix of expression data from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{S1} matrix where \eqn{G} is the number of genes and \eqn{S1} is the number of samples for component 1. }

  \item{data.comp2}{A matrix of expression data from additional reference samples. It is a \eqn{G} by \eqn{S2} matrix where \eqn{G} is the number of genes and \eqn{S2} is the number of samples for component 2. Component 2 is needed only for running a three-component model.}
  
 \item{givenpi}{A vector of proportions for all mixed tumor samples.
 In two-component analysis, it gives the proportions of the known reference component, and in three-component analysis, it gives the proportions for the two known components.}
 
 \item{nbin}{The number of bins used in numerical integration for computing complete likelihood. A larger value increases accuracy in estimation but increases the running time, especially in a three-component deconvolution problem. The default is 50.}
  
  \item{nthread}{The number of threads used for deconvolution when OpenMP is availble in the system. The default is the number of whole threads minus one. In our no-OpenMP version, it is set to 1.}
}

\value{
  \item{decovExprT}{Matrix of deconvolved expression profiles corresponding to T-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.}  
  
  \item{decovExprN1}{Matrix of deconvolved expression profiles corresponding to N1-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.} 
  
  \item{decovExprN2}{Matrix of deconvolved expression profiles corresponding to N2-component in mixed samples for a given subset of genes in a three-component setting. Each row corresponds to one gene and each column corresponds to one sample.} 
 
  \item{decovMu}{Estimated \eqn{\mu} of log2-normal distribution for both known (\eqn{MuN1, MuN2}) and unknown component (\eqn{MuT}).}
  
 \item{decovSigma}{Estimated \eqn{\sigma} of log2-normal distribution for both known (\eqn{SigmaN1, SigmaN2}) and unknown component (\eqn{SigmaT}).}
 }

\author{Zeya Wang, Wenyi Wang}

\seealso{http://bioinformatics.mdanderson.org/main/DeMix:Overview}
 
\examples{
# Example 1: two-component deconvolution given proportions 
  data(test.data1)
  givenpi <- c(t(as.matrix(test.data1$truth[-2,])))
  res <- DeMixT.S2(data.Y = test.data1$y, data.comp1 = test.data1$comp1, givenpi = givenpi)
 
# Example 2: three-component deconvolution given proportions 
# This example takes 10 minutes to finish running
#  data(test.data2)
#  givenpi <- c(t(test.data2$truth[-3,]))
#  res <- DeMixT.S2(data.Y = test.data2$y, data.comp1 = test.data2$comp1, data.comp2 = test.data2$comp2, givenpi = givenpi)
}

\references{J. Besag. “On the statistical analysis of dirty pictures”. In: Journal of the Royal Statistical Society. Series B (Methodological) (1986), pp. 259–302.
}

\keyword{DeMixT.S2}


