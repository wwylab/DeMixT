#' @title Simulated three-component mixed cell line test data
#' 
#' @description A list of simulated three-component mixed cell line
#' test data used in DeMixT function. Expression data with 700 genes 
#' and 20 samples are simulated, where 675 genes' \eqn{MuN1} is close 
#' to \eqn{MuN2}.
#' 
#' @return A list with 6 elements (3 more elements when output.more.info = 
#' TRUE), which are
#' \item{pi}{A matrix of estimated proportion. First row and second row 
#' corresponds to the proportion estimate for the known components and unkown 
#' component respectively for two or three component settings. Each column 
#' corresponds to one sample.}
#' \item{Mu}{Simulated \eqn{Mu} of log2-normal distribution for both known
#' (\eqn{MuN1, MuN2}) and unknown component (\eqn{MuT}).}
#' \item{Sigma}{Simulated \eqn{Sigma} of log2-normal distribution for both 
#' known (\eqn{SigmaN1, SigmaN2}) and unknown component
#'  (\eqn{SigmaT}).}
#' \item{data.Y}{A SummarizedExperiment object of simulated expression data 
#' from mixed tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G}
#' is the number of genes and \eqn{My} is the number of mixed samples. 
#' Samples with the same tissue type should be placed together in columns.}
#' \item{data.N1}{A SummarizedExperiment object of simulated expression data 
#' from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix 
#' where \eqn{G} is the number of genes and \eqn{M1} is the number of samples 
#' for component 1.}
#' \item{data.N2}{A SummarizedExperiment object of expression data from
#' additional reference samples. It is a \eqn{G} by \eqn{M2} matrix where 
#' \eqn{G} is the number of genes and \eqn{M2} is the number of samples for
#' component 2.}
#' \item{True.data.T}{A SummarizedExperiment object of simulated tumor expression 
#' data. It is a \eqn{G} by \eqn{My} matrix, where \eqn{G} is the number of 
#' genes and \eqn{My} is the number of mixed samples.This is shown only when 
#' output.more.info = TRUE.}
#' \item{True.data.N1}{A SummarizedExperiment object of simulated true 
#' expression data for reference component 1 (e.g., stroma). It is a \eqn{G} 
#' by \eqn{M1} matrix where \eqn{G} is the number of genes and \eqn{M1} is the 
#' number of samples for component 1. This is shown only when 
#' output.more.info = TRUE.}
#' \item{True.data.N2}{A SummarizedExperiment object of simulated true 
#' expression data for reference component 2 (e.g., immue). It is a \eqn{G} 
#' by \eqn{M2} matrix where \eqn{G} is the number of genes and \eqn{M2} is the 
#' number of samples for component 2. This is shown only when 
#' output.more.info = TRUE.}
'test.data.3comp'