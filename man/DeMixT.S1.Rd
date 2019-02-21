\name{DeMixT.S1}
\alias{DeMixT.S1}

\title{
Estimates the proportions of mixed samples for each mixing component 
}
\description{
This function is designed to estimate the proportions of all mixed samples for each mixing component with or without component merging.
}
\usage{
DeMixT(data.Y, data.comp1, data.comp2 = NULL, niter = 10, nbin = 50, 
       if.filter = TRUE, ngene.selected.for.pi = 250, mean.diff.in.CM = 0.25, 
       tol = 10^(-5), nthread = detectCores() - 1) 
}

\arguments{
   \item{data.Y}{A SummarizedExperiment object of expression data from mixed tumor samples. It is a \eqn{G} by \eqn{Sy} matrix where \eqn{G} is the number of genes and \eqn{Sy} is the number of mixed samples. Samples with the same tissue type should be placed together in columns.}

  \item{data.comp1}{A SummarizedExperiment object of expression data from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{S1} matrix where \eqn{G} is the number of genes and \eqn{S1} is the number of samples for component 1. }

  \item{data.comp2}{A SummarizedExperiment object of expression data from additional reference samples. It is a \eqn{G} by \eqn{S2} matrix where \eqn{G} is the number of genes and \eqn{S2} is the number of samples for component 2. Component 2 is needed only for running a three-component model.}
  
  \item{niter}{The maximum number of iterations used in the algorithm of iterated conditional modes (ICM, Ref[1]). A larger value better guarantees the convergence in estimation but increases the running time. The default is 10.}
 
  \item{nbin}{The number of bins used in numerical integration for computing complete likelihood. A larger value increases accuracy in estimation but increases the running time, especially in a three-component deconvolution problem. The default is 50.}
  
  \item{if.filter}{The logical flag indicating whether a predetermined filter rule is used to select genes for proportion estimation. The default is TRUE.}
  
  \item{ngene.selected.for.pi}{The percentage or the number of genes used for proportion estimation. The difference between the expression levels from mixted tumor samples and the known component(s) are evaluated, and the most differentially expressed genes are selected. It is enabled when if.filter = TRUE. The default is 250.}
  
  \item{mean.diff.in.CM}{Threshold of expression difference in selecting genes in the component merging strategy. We merge three-component to two-component by selecting genes with similar expressions for the two known components. Genes with the mean differences less than the threshold will be selected for component merging. It is used in the three-component setting, and is enabled when if.filter = TRUE. The default is 0.25.}
  
  \item{tol}{The convergence criterion. The default is 10^(-5).}

  \item{nthread}{The number of threads used for deconvolution when OpenMP is availble in the system. The default is the number of whole threads minus one. In our no-OpenMP version, it is set to 1.}
}

\value{
  \item{pi}{Matrix of estimated proportions for each known component. \eqn{\pi1} corresponds to the proportion estimate for the first known component. \eqn{\pi2} corresponds to the second known component.} 
  
 \item{pi.iter}{Estimated proportions in each iteration. It is a \eqn{number of iteration X Sy X 1} array in two-component setting, and a \eqn{number_of_iteration X Sy X 2} array in three-component setting. This is enabled only when output.iter = TRUE.}
 
 \item{gene.name}{The names of genes used in estimating the proportions. If no gene names are rpovided in the original data set, the genes will be automatically indexed.}
}

\author{Zeya Wang, Wenyi Wang}

\seealso{http://bioinformatics.mdanderson.org/main/DeMixT}
 
\examples{
# Example 1: estimate proportions for simulated two-component data 
 data(test.data1.y)
 data(test.data1.comp1)
 res <- DeMixT.S1(data.Y = test.data1.y, data.comp1 = test.data1.comp1, if.filter = FALSE)

# Example 2: estimate proportions for simulated three-component data
# This example takes 10 minutes to finish running
# data(test.data2.y)
# data(test.data2.comp1)
# data(test.data2.comp2)
# res <- DeMixT.S1(data.Y = test.data2.y, data.comp1 = test.data2.comp1, data.comp2 = test.data2.comp2, if.filter = FALSE)
	
# Example 3: estimate proportions for simulated three-component  mixed cell line data 
# This example takes 1 hour to finish running
# data(test.data2.y)
# data(test.data3.comp1)
# data(test.data3.comp2)
# res <- DeMixT.S1(data.Y = test.data3.y, data.comp1 = test.data3.comp1, data.comp2 = test.data3.comp2)
}

\references{J. Besag. “On the statistical analysis of dirty pictures”. In: Journal of the Royal Statistical Society. Series B (Methodological) (1986), pp. 259–302.
}

\keyword{DeMixT.S1}

