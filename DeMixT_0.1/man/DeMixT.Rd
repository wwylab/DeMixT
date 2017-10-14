\name{DeMixT}
\alias{DeMixT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deconvolution of heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms
}
\description{
This is the main function to run the deconvolution analysis in a setting of two or three components.
}
\usage{
DeMixT(inputdata, groupid, niter = 10, ninteg1 = 50,ninteg2 = 50, filter.out = TRUE, 
filter.option = 1, filter.criteria1 = c(0.5,0.5), filter.criteria2 = c(250,250), filter.criteria3 = 0.25, 
if.filter = FALSE, tol = 10^(-5), sg0 = 0.5^2, mu0= 0.0, nthread = detectCore() - 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputdata}{
A matrix of expression data (e.g gene expressions) from reference (e.g. normal) and mixed samples (e.g. mixed tumor samples). It is a \eqn{G \times S} matrix where \eqn{G} is the number of genes and \eqn{S} is the number of samples including reference and mixed samples. Samples with the same tissue type should be placed together in columns (e.g. cbind(normal samples, mixed tumor samples).}
  \item{groupid}{
A vector of indicators to denote if the corresponding samples are reference samples or mixed tumor samples. DeMixT is able to deconvolve mixed tumor samples with at most three components. We use 1 and 2 to denote the samples referencing the first and the second known component in mixed tumor samples. We use 3 to indicate mixed tumor samples prepared to be deconvolved. For example, in two-component deconvolution, we have c(1,1,...,3,3) and in three-component deconvolution, we have c(1,1,...,2,2,...,3,3).}
 \item{niter}{The maximum number of iterations used in the algorithm of iterated conditional modes. A larger value can better guarantee the convergence in estimation. The default number is 10.}
\item{ninteg1}{The number of bins used in numerical integration for computing complete likelihood when proportions are estimated. A larger value can increase accuracy in estimation but also increase the running time. Especially in three-component deconvolution, the increase of number of bins can lengthen the running time. The default value is 50.}
\item{ninteg2}{The number of bins used in numerical integration for computing complete likelihood when individual expressions are deconvolved. The default value is 50.}
\item{filter.out}{Option to control if only genes with low-biased estimates are output for two-component deconvolution.}
\item{filter.option}{The option is used to process or filter zero count from input expression matrix. If it is set to 1, all those genes containing zero count in sample will be removed. If it is set to 2, input data matrix will be added 1 to vanish zero count. This option is used for RNA-seq data. Default is 1.}
\item{filter.criteria1}{Filtering threshold used to select genes with small sample standard deviations in reference samples. This threshold is set on the level of log2-transformed data. A vector of length 1 is given for two-component deconvolution and length 2 for three-component deconvolution. A larger value will relax the gene filtering for proportion estimation but might introduce bias. A smaller value will take fewer genes for proportion estimation. User should take care of the choice of this value because too stringent criteria can cause no gene to be input for estimation. It is only enabled when if.filter = TRUE. Default is 0.5 for each component.}
\item{filter.criteria2}{The percentage or the number of genes that will be used for proportion estimation after gene selection. It is a vector of length 2 for three-component setting and length 1 for two-component setting. In the two-component setting, a rank of fold change of mixed to reference is used to select most differentially expressed genes at the average level. In the three-component setting, this is used within our component merging strategy. The first element is set to give the number of genes selected in the first step when we merge three-component into two-component, and the second element of this vector is set to give the number of genes selected in the second step, where all three components are being estimated. It is enabled when if.filter = TRUE. Default value is 250.}
\item{filter.criteria3}{Threshold of mean expression difference given to select genes with close expressions between two known components. Genes with mean expression difference between two known components below this value will be included for deconvolution in the first step. This option is set in our first step of component merging strategy when we merge three-component to two-component by selecting genes with very close expressions for the two known components. More stringent value can merge from three-component to two-component more strictly but also cause fewer genes or no gene left for estimation. It is only enabled for the case of three-component. Default is 0.25.}
\item{if.filter}{Option to control if a pre-determined filtering rule is used to select genes for proportion estimation.}
\item{tol}{The convergence criterion. The default is 10^(-5).}
\item{sg0}{A single initial value of \eqn{\sigma_T^2} in log2-normal distribution for each gene. The default value is 0.5^2.}
\item{mu0}{A single initial value of \eqn{\mu_T} in log2-normal distribution for each gene. It should be a postive value. If it is assigned with any non-positive value, the initial value will use an estimator of method of moments by assuming \eqn{\pi} satisfying a Dirichlet distribution. The default value is 0.0.}
\item{nthread}{The number of threads used for deconvolution. The default is the number of whole CPUs minus one.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{pi}{Matrix of estimated proportions for each known component. The first row corresponds to the proportion estimate of each sample for the first known component (groupid = 1) and the second row corresponds to that for the second known component (groupid = 2), and the third component would be \eqn{1 - \pi_1 - \pi_2}.} 
 \item{pi_iteration}{A list of estimated proportions for the first and second known component in each iteration. The first element in the list provides estimates in each iteration for the first component and the second element provides values for the second component. In each of these two elements, each row corresponds to each gene and each column corresponds to each iteratoin, and the third component would be \eqn{1 - \pi_1 - \pi_2}.}
  \item{decovExprT}{A matrix of deconvolved expression profiles corresponding to unknown T-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample. Gene names and sample names are attached in the output.}  
 \item{decovExprN1}{A matrix of deconvolved expression profiles corresponding to known N1-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample. Gene names and sample names are attached in the ouptut.} 
 \item{decovExprN2}{A matrix of deconvolved expression profiles corresponding to known N2-component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample. Gene names and sample names are attached in the output.}   
  \item{decovMu}{Estimated \eqn{\mu} of log2-normal distribution for both known (\eqn{MuN1, MuN2}) and unknown component (\eqn{MuT}).}
 \item{decovSigma}{Estimated \eqn{\sigma} of log2-normal distribution for both known (\eqn{SigmaN1, SigmaN2}) and unknown component (\eqn{SigmaT}).}
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\author{
Zeya Wang, Wenyi Wang
}
\references{
}
\keyword{DeMixT}
\seealso{
	http://bioinformatics.mdanderson.org/main/DeMixT:Overview
 }
\examples{
	
#first example for simulated two-component data 
 data(simul2)
 inputdata <- as.matrix(simul2)
 groupid <- c(rep(1, 5), rep(3, 10));	
 res <- DeMixT(inputdata, groupid)
	
## Those two examples below will take around hours for finishing the running. 
## Be cautious if you want to test those examples.		
#second example for simulated three-component data 
 data(simul3)
 inputdata <- as.matrix(simul3)
 groupid <- c(rep(1, 5), rep(2,5), rep(3, 5));
 res <- DeMixT(inputdata, groupid)
	
#third example of applying component merging strategy for three-component 
#in proportion estimation to a mixed cell line data example (this example 
# takes a longer time to be finished.)
 data(cell_line_mix)
 inputdata <- as.matrix(cell_line_mix)
 groupid <- c(rep(1, 6), rep(2,6), rep(3, 8)); # input as 2-component
 res <- DeMixT(inputdata, groupid, if.filter = TRUE)
	

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DeMixT}% use one of  RShowDoc("KEYWORDS")
\keyword{DeMixT}% __ONLY ONE__ keyword per line

