#' @title Deconvolves expressions of each individual sample for unknown component
#'
#' @description This function is designed to estimate the deconvolved expressions 
#' of individual mixed tumor samples for unknown component for each gene.
#'
#' @param data.Y A SummarizedExperiment object of expression data from mixed 
#' tumor samples. It is a \eqn{G} by \eqn{My} matrix where \eqn{G} is the number
#' of genes and \eqn{My} is the number of mixed samples. Samples with the same
#' tissue type should be placed together in columns.
#' @param data.N1 A SummarizedExperiment object of expression data 
#' from reference component 1 (e.g., normal). It is a \eqn{G} by \eqn{M1} matrix 
#' where \eqn{G} is the number of genes and \eqn{M1} is the number of samples 
#' for component 1. 
#' @param data.N2 A SummarizedExperiment object of expression data from
#' additional reference samples. It is a \eqn{G} by \eqn{M2} matrix where 
#' \eqn{G} is the number of genes and \eqn{M2} is the number of samples for
#' component 2. Component 2 is needed only for running a three-component model.
#' @param givenpi A vector of proportions for all mixed tumor samples.
#' In two-component analysis, it gives the proportions of the unknown reference
#' component, and in three-component analysis, it gives the proportions for the 
#' two known components.
#' @param nbin Number of bins used in numerical integration for computing 
#' complete likelihood. A larger value increases accuracy in estimation but
#' increases the running time, especially in a three-component deconvolution 
#' problem. The default is 50.
#' @param nthread The number of threads used for deconvolution when OpenMP is
#' available in the system. The default is the number of whole threads minus one.
#' In our no-OpenMP version, it is set to 1.
#'
#' @return 
#' \item{decovExprT}{A matrix of deconvolved expression profiles corresponding to 
#' T-component in mixed samples for a given subset of genes. Each row 
#' corresponds to one gene and each column corresponds to one sample.}  
#' \item{decovExprN1}{A matrix of deconvolved expression profiles corresponding to 
#' N1-component in mixed samples for a given subset of genes. Each row 
#' corresponds to one gene and each column corresponds to one sample.} 
#' \item{decovExprN2}{A matrix of deconvolved expression profiles corresponding to 
#' N2-component in mixed samples for a given subset of genes in a 
#' three-component setting. Each row corresponds to one gene and each 
#' column corresponds to one sample.} 
#' \item{decovMu}{A matrix of estimated \eqn{Mu} of log2-normal distribution for 
#' both known (\eqn{MuN1, MuN2}) and unknown component (\eqn{MuT}). Each row 
#' corresponds to one gene.} 
#' \item{decovSigma}{Estimated \eqn{Sigma} of log2-normal distribution for both 
#' known (\eqn{SigmaN1, SigmaN2}) and unknown component (\eqn{SigmaT}). Each 
#' row corresponds to one gene.}
#' 
#' @author Zeya Wang, Wenyi Wang
#' 
#' @seealso http://bioinformatics.mdanderson.org/main/DeMixT
#' 
#' @examples
#' # Example 1: two-component deconvolution given proportions 
#'   data(test.data.2comp)
#'   givenpi <- c(t(as.matrix(test.data.2comp$pi[-2,])))
#'   res.S2 <- DeMixT_S2(data.Y = test.data.2comp$data.Y, 
#'                       data.N1 = test.data.2comp$data.N1,
#'                       data.N2 = NULL, 
#'                       givenpi = givenpi, 
#'                       nbin = 50)
#' #                  
#' # Example 2: three-component deconvolution given proportions 
#' # data(test.data.3comp)
#' # givenpi = c(t(test.data.3comp$pi[-3,])) 
#' # res <- DeMixT_S2(data.Y = test.data.3comp$data.Y, 
#' #                  data.N1 = test.data.3comp$data.N1,
#' #                  data.N2 = test.data.3comp$data.N2, 
#' #                  givenpi = givenpi, 
#' #                  nbin = 50)
#' 
#' @references Wang Z, Cao S, Morris J S, et al. Transcriptome Deconvolution of 
#' Heterogeneous Tumor Samples with Immune Infiltration. iScience, 2018, 9: 451-460.
#' 
#' @keywords DeMixT_S2
#' 
#' @export 
DeMixT_S2 <- function(data.Y, data.N1, data.N2 = NULL, 
                      givenpi, nbin = 50, 
                      nthread = parallel::detectCores() - 1) 
{
  filter.out = TRUE
  filter.option = 1
  data.Y <- SummarizedExperiment::assays(data.Y)[[1]]
  data.N1 <- SummarizedExperiment::assays(data.N1)[[1]]
  if (!is.null(data.N2)){
    data.N2 <- SummarizedExperiment::assays(data.N2)[[1]]
  }
  if (is.null(rownames(data.Y))) {
    rownames(data.Y) <- paste('Gene', seq = seq(1, nrow(data.Y)))

    rownames(data.N1) <- paste('Gene', seq = seq(1, nrow(data.N1)))
    if(!is.null(data.N2)){
      rownames(data.N2) <- paste('Gene', seq = seq(1, nrow(data.N2)))
    }
  }
  if (is.null(colnames(data.Y))) {
    colnames(data.Y) <- paste('Sample', seq = seq(1, ncol(data.Y)))
    colnames(data.N1) <- paste('Sample', seq = seq(1, ncol(data.N1)))
    if(!is.null(data.N2)){
      colnames(data.N2) <- paste('Sample', seq = seq(1, ncol(data.N2)))
    }
  }
  sample.id <- colnames(data.Y)
  if (is.null(data.N2)) {
    inputdata <- cbind(data.N1, data.Y)
    groupid <- c(rep(1, ncol(data.N1)), rep(3, ncol(data.Y)))
  }else {
    inputdata <- cbind(data.N1, data.N2, data.Y)
    groupid <- c(rep(1, ncol(data.N1)), rep(2, ncol(data.N2)), 
                 rep(3, ncol(data.Y)))
  }
  if (filter.option == 1) {
    index <- apply(inputdata, 1, function(x) sum(x <= 0) == 
                     0)
    inputdata <- inputdata[index, ]
    data.N1 <- data.N1[index, ]
    if (!is.null(data.N2)) 
      data.N2 <- data.N2[index, ]
    data.Y <- data.Y[index, ]
  }else if (filter.option == 2) {
    inputdata <- inputdata + 1
    data.N1 <- data.N1 + 1
    if (!is.null(data.N2)) 
      data.N2 <- data.N2 + 1
    data.Y <- data.Y + 1
  }else {
    stop("The argument filter.option can only be 1 or 2")
  }
  
  if (is.null(data.N2)) {
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.N1, 
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.N1, 1, function(x) 
        length(unique(x)) > 1), ]
    }
  }else {
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.N1, 
        1, function(x) length(unique(x)) > 1) & apply(data.N2, 
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.N1, 1, function(x) 
        length(unique(x)) > 1) & apply(data.N2, 
        1, function(x) length(unique(x)) > 1), ]
    }
  }

  gene.id <- rownames(inputdata)
  inputdata <- as.matrix(inputdata)
  res <- Optimum_KernelC(inputdata, groupid, setting.pi = 1,  nspikein = 0, 
    givenpi = givenpi, givenpiT = rep(0, ncol(data.Y)), niter = 1, ninteg = nbin, 
    tol = 1e-05, nthread = nthread)

  
  decovExprT <- res$decovExprT
  decovExprT[which(decovExprT < 0,arr.ind = TRUE)] <- 0
  colnames(decovExprT) <- sample.id
  row.names(decovExprT) <- gene.id
  decovExprN1 <- res$decovExprN1
  decovExprN1[which(decovExprN1 < 0,arr.ind = TRUE)] <- 0
  colnames(decovExprN1) <- sample.id
  row.names(decovExprN1) <- gene.id
  decovExprN2 <- res$decovExprN2
  decovExprN2[which(decovExprN2 < 0,arr.ind = TRUE)] <- 0
  colnames(decovExprN2) <- sample.id
  row.names(decovExprN2) <- gene.id
  
  
  if (is.null(data.N2)) {
    
    decovMuN <- rowMeans(log2(data.N1[gene.id, ]))
    decovSigmaN <- apply(log2(data.N1[gene.id, ]), 1 , sd)
    
    decovMu <- cbind(decovMuN, res$decovMu)
    colnames(decovMu) <- c("MuN1", "MuT")
    row.names(decovMu) <- gene.id
    decovSigma <- cbind(decovSigmaN, res$decovSigma)
    colnames(decovSigma) <- c("SigmaN1", "SigmaT")
    row.names(decovSigma) <- gene.id
    filter.in <- rep(TRUE, nrow(inputdata))
    if (filter.out == TRUE) filter.in <- (decovMuN - res$decovMu <= 4)
    return(list(decovExprT = decovExprT[filter.in, ], decovExprN1 = 
      decovExprN1[filter.in, ], decovExprN2 = decovExprN2[filter.in, ], 
      decovMu = decovMu[filter.in, ], decovSigma = decovSigma[filter.in, ]))
  }else {
    
    decovMuN1 <- rowMeans(log2(data.N1[gene.id,]))
    
    decovSigmaN1 <- apply(log2(data.N1[gene.id,]), 1, sd)
    
    decovMuN2 <- rowMeans(log2(data.N2[gene.id,]))
    
    decovSigmaN2 <- apply(log2(data.N2[gene.id,]), 1, sd)
    decovMu <- cbind(decovMuN1, decovMuN2, res$decovMu)
    colnames(decovMu) <- c("MuN1", "MuN2", "MuT")
    row.names(decovMu) <- gene.id
    decovSigma <- cbind(decovSigmaN1, decovSigmaN2, res$decovSigma)
    colnames(decovSigma) <- c("SigmaN1", "SigmaN2", "SigmaT")
    row.names(decovSigma) <- gene.id
    return(list(decovExprT = decovExprT, decovExprN1 = decovExprN1, 
                decovExprN2 = decovExprN2, decovMu = decovMu, 
                decovSigma = decovSigma))
  }
}