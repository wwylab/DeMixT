#' @title Estimates the proportions of mixed samples for each mixing component 
#' using profile likelihood gene selection
#'
#' @description This function is designed to estimate the proportions of all 
#' mixed samples for each mixing component with a new proposed profile likelihood 
#' based gene selection, which can select most identifiable genes as reference 
#' gene sets to achieve better model fitting quality. We first calculated the
#' Hessian matrix of the parameter spaces and then derive the confidence interval 
#' of the profile likelihood of each gene. We then utilized the length of 
#' confidence interval as a metric to rank the identifiability of genes. As 
#' a result, the proposed gene selection approach can improve the tumor-specific 
#' transcripts proportion estimation.
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
#' the cut-off will be selected into the model. The default is TRUE.
#' @param ngene.Profile.selected The number of genes used for proportion
#' estimation ranked by profile likelihood. The default is 
#' \eqn{min(1500,0.1*G)}, where \eqn{G} is the number of genes. 
#' @param ngene.selected.for.pi The percentage or the number of genes used for
#' proportion estimation. The difference between the expression levels from
#' mixed tumor samples and the known component(s) are evaluated, and the most
#' differential expressed genes are selected, which is called DE. It is enabled
#' when if.filter = TRUE. The default is \eqn{min(1500, 0.3*G)}, where
#' \eqn{G} is the number of genes. Users can also try using more genes,
#' ranging from \eqn{0.3*G} to \eqn{0.5*G}, and evaluate the outcome.
#' @param mean.diff.in.CM Threshold of expression difference for selecting
#' genes in the component merging strategy. We merge three-component to 
#' two-component by selecting genes with similar expressions for the two known
#' components. Genes with the mean differences less than the threshold will 
#' be selected for component merging. It is used in the three-component 
#' setting, and is enabled when if.filter = TRUE. The default is 0.25.
#' @param nspikein The number of spikes in normal reference used for proportion
#' estimation. The default value is \eqn{ min(200, 0.3*My)}, where 
#' \eqn{My} the number of mixed samples. If it is set to 0, proportion 
#' estimation is performed without any spike in normal reference.
#' @param tol The convergence criterion. The default is 10^(-5).
#' @param pi01 Initialized proportion for first kown component. The default is 
#' \eqn{Null} and pi01 will be generated randomly from uniform distribution.
#' @param pi02 Initialized proportion for second kown component. pi02 is needed 
#' only for running a three-component model. The default is \eqn{Null} and pi02 
#' will be generated randomly from uniform distribution.
#' @param nthread The number of threads used for deconvolution when OpenMP 
#' is available in the system. The default is the number of whole threads
#' minus one. In our no-OpenMP version, it is set to 1.
#'
#' @return
#' \item{pi}{A matrix of estimated proportion. First row and second row 
#' corresponds to the proportion estimate for the known components and unkown 
#' component respectively for two or three component settings, and each column 
#' corresponds to one sample.}
#' \item{pi.iter}{Estimated proportions in each iteration. It is a \eqn{niter 
#' *My*p} array, where \eqn{p} is the number of components. This is 
#' enabled only when output.more.info = TRUE.}
#' \item{gene.name}{The names of genes used in estimating the proportions. 
#' If no gene names are provided in the original data set, the genes will
#' be automatically indexed.}
#'
#'@note A Hessian matrix file will be created in the working directory and the
#' corresponding Hessian matrix with an encoded name from the mixed tumor
#' sample data will be saved under this file. If a user reruns this function 
#' with the same dataset, this Hessian matrix will be loaded to in place of 
#' running the profile likelihood method and reduce running time.
#'
#'@author Shaolong Cao, Zeya Wang, Wenyi Wang
#' 
#' @seealso http://bioinformatics.mdanderson.org/main/DeMixT
#'
#' @examples
#' 
#' # Example 1: estimate proportions for simulated two-component data 
#' # with spike-in normal reference
#'   data(test.data.2comp)
#' # res.GS = DeMixT_GS(data.Y = test.data.2comp$data.Y, 
#' #                    data.N1 = test.data.2comp$data.N1,
#' #                    niter = 10, nbin = 50, nspikein = 50,
#' #                    if.filter = TRUE, ngene.Profile.selected = 150,
#' #                    mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#' #                    tol = 10^(-5))
#' #
#' # Example 2: estimate proportions for simulated two-component data 
#' # without spike-in normal reference
#' # data(test.dtat.2comp)
#' # res.GS = DeMixT_GS(data.Y = test.data.2comp$data.Y, 
#' #                    data.N1 = test.data.2comp$data.N1,
#' #                    niter = 10, nbin = 50, nspikein = 0,
#' #                    if.filter = TRUE, ngene.Profile.selected = 150,
#' #                    mean.diff.in.CM = 0.25, ngene.selected.for.pi = 150,
#' #                    tol = 10^(-5))
#' 
#' 
#' @references Gene Selection and Identifiability Analysis of RNA 
#' Deconvolution Models using Profile Likelihood. Manuscript in 
#' preparation.
#' 
#' @keywords DeMixT_GS
#' 
#' @export
DeMixT_GS <- function(data.Y, data.N1, data.N2 = NULL, 
                      niter = 10, nbin = 50, if.filter = TRUE,
                      filter.sd = 0.5, ngene.Profile.selected = NA, 
                      ngene.selected.for.pi = NA, 
                      mean.diff.in.CM = 0.25, nspikein = NULL,
                      tol = 10^(-5), pi01 = NULL, pi02 = NULL,
                      nthread = parallel::detectCores() - 1) {
  
  IF_inverse <- function(m){
    return(class(try(solve(m),silent=T))=="matrix")
  } 
  ## Creat a folder for saving hessian matrix
  path = getwd()
  if (!file.exists(paste0(path,'/Hessian_Matrix'))){
    dir.create(path = paste0(path, '/Hessian_Matrix'))
  } 
  ## Transfering data format
  data.Y <- SummarizedExperiment::assays(data.Y)[[1]]
  data.N1 <- SummarizedExperiment::assays(data.N1)[[1]]
  nS = ncol(data.Y)
  ## Encode data.Y for hessian matrix
  if(nS > 10){
    data.Y_encode = base64encode(log2(as.matrix(data.Y)[1, 1:10] + 1))
  }else{
    data.Y_encode = base64encode(log2(as.matrix(data.Y)[1, 1:nS] + 1))
  }
  if (!is.null(data.N2)){
    message('File for Hessian matrix has existed')
    data.N2 <- SummarizedExperiment::assays(data.N2)[[1]]
    nspikein = 0
  }
  if (is.na(ngene.selected.for.pi)){
    ngene.selected.for.pi = min(1500, round(0.3*nrow(data.Y)))
  } 
  ## Generate Spike-in normal reference 
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

  ## Create the default value for ngene.Profile.selected
  filter.option = 1
  if(is.na(ngene.Profile.selected)){
    ngene.Profile.selected<-min(1500, round(0.3*nrow(data.Y)))
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
  
  
  ## case 1
  if (if.filter == FALSE) {
    gene.name <- rownames(inputdata)
    res <- Optimum_KernelC(inputdata, groupid, setting.pi = 0, 
                           nspikein = nspikein,
      givenpi = rep(0, 2 * ncol(data.Y)), 
      givenpiT = rep(0, ncol(data.Y)),
      niter = niter, ninteg = nbin, tol = tol, 
      pi01 = pi01, pi02 = pi02,
      nthread = nthread)
  }
  
  ## Gene selection based on profile likelihood method
  if (if.filter == TRUE & is.null(data.N2)) {
    # step 2 identifiability selection
    cat("Gene selection starts\n")
    inputdatans <- apply(log2(data.N1), 1, sd)
    id1 <- (inputdatans < filter.sd)
    if (sum(id1) < 20) 
      stop("The threshold of standard variation is too stringent. \n
           Please provide a larger threshold. ")
    inputdatamat1 <- inputdata[id1, ]

    inputdatamat1nm <- rowMeans(inputdatamat1[, groupid == 1])
    inputdatamat1ym <- rowMeans(inputdatamat1[, groupid == 3])
    inputdata1r <- inputdatamat1ym/inputdatamat1nm
    inputdata2 <- filter2(inputdata1r, ngene.selected.for.pi)
    
    gene.name.S1 <- rownames(inputdata2)
    res.S1 <- Optimum_KernelC(inputdata2, groupid, setting.pi = 0, 
                              nspikein = nspikein,
          givenpi = rep(0, 2 * ncol(data.Y)), 
          givenpiT = rep(0, ncol(data.Y)),
          niter = niter, ninteg = nbin, tol = tol, 
          pi01 = pi01, pi02 = pi02,
          nthread = nthread)
    
    Pi.all<-colSums(res.S1$pi)
    PiT.all <- 1-colSums(res.S1$pi)
    MuT.S1 <- as.numeric(res.S1$decovMu[,dim(res.S1$decovMu)[2]])
    SigmaT.S1 <- as.numeric(res.S1$decovSigma[,dim(res.S1$decovMu)[2]])
    sdn.obs<-apply(log2(inputdata2[, groupid == 1]+0.001),1,sd)
    mun.obs<-rowMeans(log2(inputdata2[, groupid == 1]+0.001))
    
    n.normal<-sum(groupid == 1)
    n.mix<-sum(groupid == 3)
    
    res.all <- Optimum_KernelC(inputdatamat1, groupid, setting.pi = 1, 
                               nspikein = nspikein, givenpi = Pi.all, 
                               givenpiT = PiT.all, niter = 1, ninteg = nbin, 
                               pi01 = pi01, pi02 = pi02,
                               tol = tol, nthread = nthread)
    
    MuT.all <- as.numeric(res.all$decovMu[,dim(res.all$decovMu)[2]])
    SigmaT.all <- as.numeric(res.all$decovSigma[,dim(res.all$decovMu)[2]])
    sdn.obs<-apply(log2(inputdatamat1[, groupid == 1]+0.001),1,sd)
    mun.obs<-rowMeans(log2(inputdatamat1[, groupid == 1]+0.001))
    
    n.normal<-sum(groupid == 1)
    n.mix<-sum(groupid == 3)
    
    if(file.exists(paste0(path, '/Hessian_Matrix/', 
                          data.Y_encode,'_Hessian.RData'))){
      load(paste0(path, '/Hessian_Matrix/', 
                  data.Y_encode,'_Hessian.RData'))
      message(paste0('Loading Hessian matrix from ',
                     paste0(path, '/Hessian_Matrix/', 
                            data.Y_encode,'_Hessian.RData \n')))
    }else{
      if(nspikein > 0){
        Hessian <- D2Loglikelihood_2D(y=inputdatamat1[,(n.normal+1):
                                                        (n.normal+n.mix-nspikein)], 
                                      Pi=Pi.all, MuN=mun.obs, MuT=MuT.all, 
                                      SigmaN=sdn.obs, SigmaT=SigmaT.all)
      }else{
        Hessian <- D2Loglikelihood_2D(y=inputdatamat1[,(n.normal+1):
                                                        (n.normal+n.mix)], 
                                      Pi=Pi.all, MuN=mun.obs, MuT=MuT.all, 
                                      SigmaN=sdn.obs, SigmaT=SigmaT.all)
      }
      save(Hessian, file = paste0(path, '/Hessian_Matrix/',
                                  data.Y_encode,'_Hessian.RData'))
      message(paste0('Hessian matrix has been saved in ',
                     paste0(path, '/Hessian_Matrix/', 
                            data.Y_encode,'_Hessian.RData \n')))
    }
    
    ## Check if the Hessian matrix contains infinity
    if(sum(!is.finite(Hessian)) == 0){
      
      if(nspikein > 0 ){
        S <- n.mix - nspikein
      }else{
        S <- n.mix
      }
      G <- nrow(inputdatamat1)
      if(IF_inverse(Hessian)){
        C.inverse <- 2*diag(solve(Hessian))
      }else{
        C.inverse <- rep(NA, S)
      }
      C.marginal <- abs(2/diag(Hessian))
      if(sum(is.na(C.inverse))>1){
        C <- C.marginal
        Non.inverse<-TRUE
      }else{
        C <- C.inverse
        Non.inverse<-FALSE
      }
      
      if(nspikein > 0){
        theta.valid <- x_update_2D(Pi.all[c(1:(length(Pi.all) - nspikein))], 
                                   MuT.all, SigmaT.all, S, G)
      }else{
        theta.valid <- x_update_2D(Pi.all, MuT.all, SigmaT.all, S, G)
      }
      CI.upper <- theta.valid + sqrt(abs(qchisq(0.95,1)*C))
      CI.upper.MuT <- x_update_inv_2D(CI.upper,S,G)$MuT
      CI.upper.SigmaT <- x_update_inv_2D(CI.upper,S,G)$SigmaT
      CI.lower <- theta.valid - sqrt(abs(qchisq(0.95,1)*C))
      CI.lower.MuT <- x_update_inv_2D(CI.lower,S,G)$MuT
      CI.lower.SigmaT <- x_update_inv_2D(CI.lower,S,G)$SigmaT
      
      id1.ID <- order(CI.upper.MuT-CI.lower.MuT)
      
      id2 <- id1.ID[1:min(ngene.Profile.selected, length(id1.ID))]
      inputmat<-inputdatamat1[id2,]
      
      pi01 = Pi.all; pi02 = PiT.all;
      res <- Optimum_KernelC(inputmat, groupid, setting.pi = 0, 
                             nspikein = nspikein, 
                             givenpi = rep(0, 2 * n.mix), 
                             givenpiT = rep(0, n.mix), 
                             niter = niter, ninteg = nbin, 
                             tol = tol, pi01 = pi01, pi02 = pi02,
                             nthread=nthread)
      gene.name <- rownames(inputmat)
      
    }else{
        message('Hessian matrix has infinity value, switch to S1 \n') 
        res = res.S1
        gene.name <- rownames(inputdata2)
    }
  }
  
  
  ## case 3: two-stage filtering, three components
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
    res1 <- Optimum_KernelC(inputdatamat2, cnvgroup, setting.pi = 0, 
                            nspikein = nspikein,
                            givenpi = rep(0, ncol(data.Y)), 
                            givenpiT = rep(0, ncol(data.Y)),
                            niter = niter, ninteg = nbin, tol = tol,
                            pi01 = pi01, pi02 = pi02,
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
    res <- Optimum_KernelC(inputdatamat3, groupid, setting.pi = 2, 
                           nspikein = nspikein,
                           givenpi = rep(0, 2 * ncol(data.Y)), 
                           givenpiT = fixed.piT, 
                           niter = niter, ninteg = nbin, tol = tol, 
                           pi01 = pi01, pi02 = pi02,
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
  return(list(pi = pi,  
  pi.iter = pi.iter, gene.name = gene.name))
}
