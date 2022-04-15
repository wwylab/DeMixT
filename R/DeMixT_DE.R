#' @title Estimates the proportions of mixed samples for each mixing component 
#'
#' @description This function is designed to estimate the deconvolved 
#' expressions of individual mixed tumor samples for unknown component 
#' for each gene.
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
#' @param niter The maximum number of iterations used in the algorithm of 
#' iterated conditional modes. A larger value better guarantees 
#' the convergence in estimation but increases the running time. The default is 
#' 10.
#' @param nbin The number of bins used in numerical integration for computing
#' complete likelihood. A larger value increases accuracy in estimation but
#' increases the running time, especially in a three-component deconvolution
#' problem. The default is 50.
#' @param if.filter The logical flag indicating whether a predetermined filter
#' rule is used to select genes for proportion estimation. The default is TRUE.
#' @param filter.sd The cut-off for the standard deviation of lognormal 
#' distribution. Genes whose log transferred standard deviation smaller than
#' the cut-off will be selected into the model. The default is 0.5.
#' @param ngene.selected.for.pi The percentage or the number of genes used for
#' proportion estimation. The difference between the expression levels from
#' mixed tumor samples and the known component(s) are evaluated, and the most
#' differential expressed genes are selected, which is called DE. It is enabled
#' when if.filter = TRUE. The default is \eqn{min(1500, 0.3*G)}, where
#' \eqn{G} is the number of genes. Users can also try using more genes,
#' ranging from \eqn{0.3*G} to \eqn{0.5*G}, and evaluate the outcome.
#' @param nspikein The number of spikes in normal reference used for proportion
#' estimation. The default value is \eqn{ min(200, 0.3*My)}, where 
#' \eqn{My} the number of mixed samples. If it is set to 0, proportion 
#' estimation is performed without any spike in normal reference.
#' @param mean.diff.in.CM Threshold of expression difference for selecting genes
#' in the component merging strategy. We merge three-component to two-component
#' by selecting genes with similar expressions for the two known components.
#' Genes with the mean differences less than the threshold will be selected for
#' component merging. It is used in the three-component setting, and is enabled
#' when if.filter = TRUE. The default is 0.25.
#' @param tol The convergence criterion. The default is 10^(-5).
#' @param pi01 Initialized proportion for first kown component. The default is 
#' \eqn{Null} and pi01 will be generated randomly from uniform distribution.
#' @param pi02 Initialized proportion for second kown component. pi02 is needed 
#' only for running a three-component model. The default is \eqn{Null} and pi02 
#' will be generated randomly from uniform distribution.
#' @param nthread The number of threads used for deconvolution when OpenMP is
#' available in the system. The default is the number of whole threads minus one.
#' In our no-OpenMP version, it is set to 1.
#'
#' @return 
#' \item{pi}{A matrix of estimated proportion. First row and second row 
#' corresponds to the proportion estimate for the known components and unkown 
#' component respectively for two or three component settings, and each column 
#' corresponds to one sample.}
#' \item{pi.iter}{Estimated proportions in each iteration. It is a \eqn{niter 
#' *Ny*p} array, where \eqn{p} is the number of components. This is 
#' enabled only when output.more.info = TRUE.}
#' \item{gene.name}{The names of genes used in estimating the proportions. 
#' If no gene names are provided in the original data set, the genes will
#' be automatically indexed.}
#' 
#' @author Zeya Wang, Wenyi Wang
#' 
#' @seealso http://bioinformatics.mdanderson.org/main/DeMixT
#'
#' @examples
#' # Example 1: estimate proportions for simulated two-component data 
#' # with spike-in normal reference
#'   data(test.data.2comp)
#' # res.DE = DeMixT_DE(data.Y = test.data.2comp$data.Y, 
#' #                    data.N1 = test.data.2comp$data.N1,
#' #                    niter = 10, nbin = 50, nspikein = 50,
#' #                    if.filter = TRUE, 
#' #                    mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#' #                    tol = 10^(-5))
#' #
#' # Example 2: estimate proportions for simulated two-component data 
#' # without spike-in normal reference
#' # data(test.data.2comp)
#' # res.DE = DeMixT_DE(data.Y = test.data.2comp$data.Y, 
#' #                    data.N1 = test.data.2comp$data.N1,
#' #                    niter = 10, nbin = 50, nspikein = 0,
#' #                    if.filter = TRUE, 
#' #                    mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#' #                    tol = 10^(-5))
#' #
#' # Example 3: estimate proportions for simulated three-component 
#' # mixed cell line data 
#' # data(test.data.3comp)
#' # res.DE <- DeMixT_DE(data.Y = test.data.3comp$data.Y,
#' #                     data.N1 = test.data.3comp$data.N1,
#' #                     data.N2 = test.data.3comp$data.N2, 
#' #                     if.filter = TRUE)
#' 
#' @references Wang Z, Cao S, Morris J S, et al. Transcriptome Deconvolution of 
#' Heterogeneous Tumor Samples with Immune Infiltration. iScience, 2018, 9: 451-460.
#' 
#' @keywords DeMixT_DE
#' 
#' @export 
DeMixT_DE <- function (data.Y, data.N1, data.N2 = NULL, niter = 10, 
                       nbin = 50, if.filter = TRUE, filter.sd = 0.5,
                       ngene.selected.for.pi = NA, nspikein = NULL,
                       mean.diff.in.CM = 0.25, tol = 10^(-5),
                       pi01 = NULL, pi02 = NULL,
                       nthread = parallel::detectCores() - 1) 
{
  filter.option = 1
  ## Transfering data format
  data.Y <- SummarizedExperiment::assays(data.Y)[[1]]
  data.N1 <- SummarizedExperiment::assays(data.N1)[[1]]
  if (is.na(ngene.selected.for.pi)){
    ngene.selected.for.pi = min(1500, round(0.3*nrow(data.Y)))
  } 
  nS = ncol(data.Y)
  ## Generate Spike-in normal reference 
  if (!is.null(data.N2)) nspikein = 0
  if (is.null(nspikein)) nspikein = min(200, ceiling(ncol(data.Y)*0.3))
  if (nspikein > 0){
    MuN = apply(log2(data.N1), 1, mean)
    MuN[which(!is.finite(MuN))] = 
      apply(log2(data.N1[which(!is.finite(MuN)), ] + 0.1), 1, mean)
    SigmaN = apply(log2(data.N1), 1, sd)
    SigmaN[which(!is.finite(SigmaN))] = 
      apply(log2(data.N1[which(!is.finite(SigmaN)), ] + 0.1), 1, sd)
    
    Spikein.normal = array(0, c(nrow(data.Y), nspikein))
    for(k in 1:nrow(data.Y)) {
      Spikein.normal[k, ] = 2^rnorm(n = nspikein, mean = MuN[k], sd = SigmaN[k])
      }
    data.Y = cbind(data.Y, Spikein.normal)
  }
  if (!is.null(data.N2)) 
    data.N2 <- SummarizedExperiment::assays(data.N2)[[1]]
  ## Creat row and column names for input data
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
  ## Merge the inputdata and creat group id
  if (is.null(data.N2)) {
    inputdata <- cbind(data.N1, data.Y)
    groupid <- c(rep(1, ncol(data.N1)), rep(3, ncol(data.Y)))
  }else {
    inputdata <- cbind(data.N1, data.N2, data.Y)
    groupid <- c(rep(1, ncol(data.N1)), rep(2, ncol(data.N2)), 
                 rep(3, ncol(data.Y)))
  }
  
  ## filter.option: 1 - remove genes containing zero; 2 - add 1 to to kill zeros
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
  
  ## filter out genes with constant value across all samples
  if (is.null(data.N2)) {
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.N1, 
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.N1, 1, function(x) 
        length(unique(x)) > 1), ]
    }
  }else{
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.N1, 
        1, function(x) length(unique(x)) > 1) & apply(data.N2,
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.N1, 1, function(x) 
        length(unique(x)) > 1) & apply(data.N2, 1, function(x) 
          length(unique(x)) > 1), ]
    }
  }
  
  
  filter2 <- function(inputdata1r, ngene.selected.for.pi, n = 1) {
    if ((ngene.selected.for.pi > 1) & (ngene.selected.for.pi%%1 == 
                                       0)) {
      id2 <- order(inputdata1r, decreasing = TRUE)
      id2 <- id2[seq(1, min(n * ngene.selected.for.pi, 
                            length(inputdata1r)))]
    }
    else if ((ngene.selected.for.pi < 1) & (ngene.selected.for.pi > 
                                            0)) {
      id2 <- (inputdata1r > quantile(inputdata1r, probs = 1 - 
                                       n * ngene.selected.for.pi))
    }
    else {
      stop("The argument ngene.selected.for.pi can only be 
           an integer or a percentage between 0 and 1")
    }
    if (sum(id2) < 20) 
      stop("there are too few genes for filtering stage 1.\n
           Please relax threshold for filtering ")
    return(inputdatamat1[id2, ])
  }
  
  
  inputdata = as.matrix(inputdata)
 
  
  if (if.filter == FALSE) {
    gene.name <- rownames(inputdata)                       
    res <- Optimum_KernelC(inputdata, groupid, 
                           setting.pi = 0, nspikein = nspikein,
      givenpi = rep(0, 2 * ncol(data.Y)), givenpiT = rep(0, ncol(data.Y)),
      niter = niter, ninteg = nbin, tol = tol, pi01 = pi01, pi02 = pi02,
      nthread = nthread)
  }
  if (if.filter == TRUE & is.null(data.N2)) {
    #inputdatans <- rowSds(log2(data.N1))
    inputdatans <- apply(log2(data.N1), 1, sd)
    id1 <- (inputdatans < filter.sd)
    if (sum(id1) < 20) 
      stop("The threshold of standard variation is too stringent. \n
           Please provide a larger threshold. ")
    inputdatamat1 <- inputdata[id1, ]
    inputdatamat1nm <- rowMeans(inputdatamat1[, groupid == 
                                                1])
    inputdatamat1ym <- rowMeans(inputdatamat1[, groupid == 
                                                3])
    inputdata1r <- inputdatamat1ym/inputdatamat1nm
    inputdata2 <- filter2(inputdata1r, ngene.selected.for.pi)
    gene.name <- rownames(inputdata2)
    res <- Optimum_KernelC(inputdata2, groupid, setting.pi = 0, nspikein = nspikein,
      givenpi = rep(0, 2 * ncol(data.Y)), givenpiT = rep(0, ncol(data.Y)),
      niter = niter, ninteg = nbin, tol = tol, pi01 = pi01, pi02 = pi02,
      nthread = nthread)
  }
  if (if.filter == TRUE & !is.null(data.N2)) {
    message("Fitering stage 1 starts\n")
    inputdatan1m <- rowMeans(log2(inputdata[, groupid == 
                                              1]))
    inputdatan2m <- rowMeans(log2(inputdata[, groupid == 
                                              2]))
    #inputdatan1s <- rowSds(log2(inputdata[, groupid == 1]))
    inputdatan1s <- apply(log2(inputdata[, groupid == 1]), 1, sd)
    #inputdatan2s <- rowSds(log2(inputdata[, groupid == 2]))
    inputdatan2s <- apply(log2(inputdata[, groupid == 2]), 1, sd)
    id1 <- ((abs(inputdatan1m - inputdatan2m) < mean.diff.in.CM) & 
              (inputdatan1s < filter.sd) & (inputdatan2s < filter.sd))
    if (sum(id1) < 20) 
      stop("The thresholds of standard variation and \n            
           mean difference are too stringent. \n            
           Please provide larger thresholds")
    inputdatamat1 <- inputdata[id1, ]
    inputdatamat1nm <- rowMeans(inputdatamat1[, (groupid == 
                                                   1) | (groupid == 2)])
    inputdatamat1ym <- rowMeans(inputdatamat1[, (groupid == 
                                                   3)])
    #inputdata1r <- rowSds(inputdatamat1[, groupid == 3])
    inputdata1r <- apply(inputdatamat1[, groupid == 3], 1, sd)
    inputdatamat2 <- filter2(inputdata1r, ngene.selected.for.pi)
    cnvgroup <- groupid
    cnvgroup[groupid == 2] <- 1
    res1 <- Optimum_KernelC(inputdatamat2, cnvgroup, setting.pi = 0, nspikein = nspikein,
      givenpi = rep(0, ncol(data.Y)), givenpiT = rep(0, ncol(data.Y)),
      niter = niter, ninteg = nbin, tol = tol, pi01 = pi01, pi02 = pi02,
      nthread = nthread)
    fixed.piT <- 1 - as.numeric(res1$pi[1, ])
    message("Filtering stage 1 is finished\n")
    message("Filtering stage 2 starts\n")
    id3 <- ((inputdatan1s < filter.sd) & (inputdatan2s < 
                                            filter.sd))
    if (sum(id3) < 20) 
      stop("The thresholds of standard variation and \n 
           mean difference are too stringent. \n
           Please provide larger thresholds")
    inputdatamat1 <- inputdata[id3, ]
    inputdatan1m <- rowMeans(log2(inputdatamat1[, groupid == 
                                                  1]))
    inputdatan2m <- rowMeans(log2(inputdatamat1[, groupid == 
                                                  2]))
    inputdata1d <- abs(inputdatan1m - inputdatan2m)
    inputdatamat2 <- filter2(inputdata1d, ngene.selected.for.pi, 
                             2)
    #inputdata1s <- rowSds(inputdatamat2[, groupid == 3])
    inputdata1s <- apply(inputdatamat2[, groupid == 3], 1, sd)
    id5 <- (inputdata1s > quantile(inputdata1s, probs = 0.5))
    inputdatamat3 <- inputdatamat2[id5, ]
    gene.name <- rownames(inputdatamat3)
    res <- Optimum_KernelC(inputdatamat3, groupid, 
        nspikein = nspikein, setting.pi = 2,  
      givenpi = rep(0, 2 * ncol(data.Y)), givenpiT = fixed.piT, 
      niter = niter, ninteg = nbin, tol = tol, pi01 = pi01, pi02 = pi02,
      nthread = nthread)
    message("Filtering stage 2 is finished")
  }
  
  
  if(nspikein == 0){
    pi <- rbind(t(as.matrix(res$pi[1, ])), 
                1 - t(as.matrix(res$pi[1, ])))
    row.names(pi) <- c("PiN1", "PiT")
    colnames(pi) <- colnames(data.Y)
    pi1 <- as.matrix(res$pi1)
    colnames(pi1) <- as.character(seq(1, ncol(pi1)))
    row.names(pi1) = colnames(data.Y)
    ## Save pi.iter
    pi.iter <- array(0, c(ncol(res$pi1), nrow(res$pi1), 2))
    pi.iter[ , , 1] <- t(res$pi1)
    pi.iter[ , , 2] <- 1 - pi.iter[ , , 1]
    colnames(pi.iter) <- colnames(data.Y)
    row.names(pi.iter) <- paste('iter ', seq = c(1:nrow(pi.iter)))
  }else{
    pi <- rbind(t(as.matrix(res$pi[1, c(1:nS)])),
                1 - t(as.matrix(res$pi[1, c(1:nS)])))
    row.names(pi) <- c("PiN1", "PiT")
    colnames(pi) <- colnames(data.Y)[c(1:nS)]
    pi1 <- as.matrix(res$pi1[c(1:nS),])
    colnames(pi1) <- as.character(seq(1, ncol(pi1)))
    row.names(pi1) = colnames(data.Y[,c(1:nS)])
    ## Save pi.iter
    pi.iter <- array(0, c(ncol(res$pi1), nrow(res$pi1)-nspikein, 2))
    pi.iter[ , , 1] <- t(res$pi1)[,c(1:nS)]
    pi.iter[ , , 2] <- 1 - pi.iter[ , , 1]
    colnames(pi.iter) <- colnames(data.Y)[1:nS]
    row.names(pi.iter) <- paste('iter ', seq = c(1:nrow(pi.iter)))
  }
  
  if(IQR(pi[2,]) < 0.01){
    stop('DeMixT Estimation Failed, try different initial values for pi')
  }
  
  if (!is.null(data.N2)) {
    pi <- rbind(res$pi, 
                t(as.matrix(apply(res$pi, 2, function(x) 1-sum(x)))))
    row.names(pi) <- c("PiN1" ,"PiN2", "PiT")
    colnames(pi) <- colnames(data.Y)
    pi2 <- as.matrix(res$pi2)
    colnames(pi2) <- as.character(seq(1, ncol(pi2)))
    row.names(pi2) = colnames(data.Y)
    ## Save pi.iter
    pi.iter <- array(0, c(ncol(res$pi1), nrow(res$pi1), 3))
    pi.iter[ , , 1] <- t(res$pi1)
    pi.iter[ , , 2] <- t(res$pi2)
    pi.iter[ , , 3] <- 1 - (pi.iter[ , , 1] + pi.iter[ , , 2])
    # pi.iter <- array(t(rbind(res$pi1, res$pi2)), 
    #                  dim <- c(ncol(res$pi1), nrow(res$pi1), 2))
    colnames(pi.iter) <- colnames(data.Y)
    row.names(pi.iter) <- paste('iter ', seq = c(1:nrow(pi.iter)))
  }

  cat("Estimation of Proportions finished:\n")
  print(round(t(pi),4))
  return(list(pi = pi, pi.iter = pi.iter, gene.name = gene.name))
}
