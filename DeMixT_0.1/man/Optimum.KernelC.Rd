\name{Optimum.KernelC}
\alias{Optimum.KernelC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel function to call the C code used for parameter estimation in DeMixT
}
\description{
This function is used by DeMixT.S1 and DeMixT.S2 to perform parameter estimation and expression deconvolution and invokes a set of C codes to perform the majority of the computation.
}
\usage{
Optimum.KernelC(inputdata, groupid, nhavepi, givenpi, givenpiT, niter, ninteg, tol, 
sg0 = 0.5^2, mu0= 0.0, nthread = detectCores() - 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputdata}{
A matrix of expression data (e.g gene expressions) from reference (e.g. normal) and mixed samples (e.g. mixed tumor samples). It is a \eqn{G \times S} matrix where \eqn{G} is the number of genes and \eqn{S} is the number of samples including reference and mixed samples. Samples with the same tissue type should be placed together in columns (e.g. cbind(normal samples, mixed tumor samples).}
  \item{groupid}{
A vector of indicators to denote if the corresponding samples are reference samples or mixed tumor samples. DeMixT is able to deconvolve mixed tumor samples with at most three components. We use 1 and 2 to denote the samples referencing the first and the second known component in mixed tumor samples. We use 3 to indicate mixed tumor samples prepared to be deconvolved. For example, in two-component deconvolution, we have c(1,1,...,3,3) and in three-component deconvolution, we have c(1,1,...,2,2,...,3,3).}
  \item{nhavepi}{If it is set to 0, then deconvolution is performed without any given proportions; if set to 1, deconvolution with given proportions for the first and the second known component is run; if set to 2, deconvolution is run with given proportions of T-component for component merging. This option helps to perform deconvolution in different settings.}
 \item{givenpi}{A vector of proportions for all admixed samples, with size of the number of admixed samples for two component and two times the number of admixed samples for three component. It is fixed with given proportions for the first and the second known component of mixed tumor samples, or just for one known component when there is just one type of reference tissues. It has the form of Vector\eqn{(\pi^1_{N_1},\pi^2_{N_1},...,\pi^{S_T}_{N_1},\pi^1_{N_2},\pi^2_{N_2},...,\pi^{S_T}_{N_2})}. \eqn{S_T(S_T<S)} is the number of mixed tumor samples.}
\item{givenpiT}{A vector of proportions for all admixed samples used for component merging, with size of the number of admixed samples. When nhavepi is set to 2, givenpiT is fixed with given proportions for unknown component of mixed tumor samples. This option is used when we adopt a two-step estimation strategy in deconvolution. It has the form of Vector\eqn{(\pi^1_T, \pi^2_T, \cdots, \pi^{S_T}_T)}. If option is not 2, this vector can be given with any element.}
\item{niter}{The maximum number of iterations used in the algorithm of iterated conditional modes. A larger value can better guarantee the convergence in estimation.}
\item{ninteg}{The number of bins used in numerical integration for computing complete likelihood. A larger value can increase accuracy in estimation but also increase the running time. Especially in three-component deconvolution, the increase of number of bins can lengthen the running time.}
\item{tol}{The convergence criterion.}
\item{sg0}{A single initial value of \eqn{\sigma_T^2} in log2-normal distribution for each gene. The default value is 0.5^2.}
\item{mu0}{A single initial value of \eqn{\mu_T} in log2-normal distribution for each gene. It should be a postive value. If it is assigned with any non-positive value, the initial value will use an estimator of method of moments by assuming \eqn{\pi} satisfying a Dirichlet distribution. The default value is 0.0.}
\item{nthread}{The number of threads used for deconvolution. The default is the number of whole CPUs minus one.}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{obj_val}{Final negative log-likelihood value at convergence.}

 \item{pi}{Matrix of estimated proportions for each known component. The first row corresponds to the proportion estimate of each sample for the first known component (groupid = 1) and the second row corresponds to that for the second known component (groupid = 2).} 
 \item{decovExprT}{A matrix of deconvolved expression profiles corresponding to unknown T-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.}  
 \item{decovExprN1}{A matrix of deconvolved expression profiles corresponding to known N1-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.} 
 \item{decovExprN2}{A matrix of deconvolved expression profiles corresponding to known N2-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.} 
 \item{decovMu}{Estimated \eqn{\mu} of log2-normal distribution for tumor component.}
 \item{decovSigma}{Estimated \eqn{\sigma} of log2-normal distribution for tumor component}
 \item{pi1}{An \eqn{S_T \times I} Matrix of estimated proportions for each iteration \eqn{i \in \{1, \cdots, I\}} for the first known component. \eqn{I} is the number of iterations. \eqn{S_T} is the number of admixed samples.}
 \item{pi2}{An \eqn{S_T \times I} Matrix of estimated proportions for each iteration \eqn{i \in \{1, \cdots, I\}} for the second known component. \eqn{I} is the number of iterations. \eqn{S_T} is the number of admixed samples.}
   }
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\author{
Zeya Wang
}
\references{
}
\keyword{DeMixT.Kernel}
\seealso{
		http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DeMixT}% use one of  RShowDoc("KEYWORDS")
\keyword{DeMixT}% __ONLY ONE__ keyword per line

