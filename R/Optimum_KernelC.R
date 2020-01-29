#' @title Kernel function for optimizing parameters and hidden variables in DeMixT
#' 
#' @description This function is invoked by DeMixT_GS or DeMixT_S1 and DeMixT_S2 to 
#' finish parameter estimation by iterated conditional mode algorithm and reconstitute
#' gene expression profile of all components.
#'
#' @param inputdata A matrix of expression data (e.g gene expressions) from
#' reference (e.g. normal) and mixed samples (e.g. mixed tumor samples). It is a
#' \eqn{G*M} matrix where \eqn{G} is the number of genes and \eqn{M} is the 
#' number of samples including reference and mixed samples. Samples with the 
#' same tissue type should be placed together in columns (e.g. cbind(normal
#' amples, mixed tumor samples).
#' @param groupid A vector of indicators to denote if the corresponding samples
#' are reference samples or mixed tumor samples. DeMixT is able to deconvolve
#' mixed tumor samples with at most three components. We use 1 and 2 to denote
#' the samples referencing the first and the second known component in mixed 
#' tumor samples. We use 3 to indicate mixed tumor samples prepared to be
#' deconvolved. For example, in two-component deconvolution, we have 
#' c(1,1,...,3,3) and in three-component deconvolution, we have 
#' c(1,1,...,2,2,...,3,3).
#' @param nspikein The number of spikes in normal reference used for proportion
#' estimation. The default value is \eqn{ min(200, 0.3*My)}, where 
#' \eqn{My} the number of mixed tumor samples. If it is set to 0, proportion 
#' estimation is performed without any spike in normal reference.
#' @param setting.pi If it is set to 0, then deconvolution is performed without 
#' any given proportions; if set to 1, deconvolution with given proportions 
#' for the first and the second known component is run; if set to 2, 
#' deconvolution is run with given tumor proportions. This option helps to 
#' perform deconvolution in different settings. In estimation of 
#' component-specific proportions, we use a subset of genes ; so when 
#' it is required to deconvolve another subset of genes, we just easily plug 
#' back our estimated proportions by setting this option to 1. In our two-step
#' estimation strategy in a three-component setting, this option is set to 2 to
#' implement the second step.
#' @param givenpi \eqn{ST}-Vector of proportions. Given the number of mixed
#' tumor samples is \eqn{My(My<M)}, \eqn{ST} is set to \eqn{2*My} in a
#' three-component setting and \eqn{My} in a two-component setting. When 
#' setting.pi is 1, it is fixed with the given proportions for the first and the
#' second known component of mixed tumor samples, or for one unknown 
#' component when there is just one type of reference tissues. It has the form 
#' of Vector \eqn{PiN1-1}, \eqn{PiN1-2}, ..., \eqn{PiN1-My}, \eqn{PiN2-1}, 
#' \eqn{PiN2-2}, ..., \eqn{PiN2-My}.
#' @param givenpiT \eqn{ST}-Vector of proportions. When setting.pi is set to 2,
#' givenpiT is fixed with given proportions for unknown component of mixed 
#' tumor samples. This option is used when we adopt a two-step estimation
#' strategy in deconvolution. It has the form of Vector  \eqn{PiT-1}, 
#' \eqn{PiT-2}, ..., \eqn{PiT-My}. If option is not 2, 
#' this vector can be given with any element.
#' @param niter The number of iterations used in the algorithm of iterated
#' conditional modes. A larger value can better guarantee the convergence 
#' in estimation but increase the computation time.
#' @param ninteg The number of bins used in numerical integration for computing
#' complete likelihood. A larger value can increase accuracy in estimation but
#' also increase the running time. Especially in three-component deconvolution,
#' the increase of number of bins can greatly lengthen the running time.
#' @param tol The convergence criterion. The default is 10^(-5).
#' @param sg0 Initial value for \eqn{\sigma^2}. The default is 0.5^2.
#' @param mu0 Initial value for \eqn{\mu}. The default is 0.
#' @param pi01 Initialized proportion for first kown component. The default is 
#' \eqn{Null} and pi01 will be generated randomly from uniform distribution.
#' @param pi02 Initialized proportion for second kown component. pi02 is needed 
#' only for running a three-component model. The default is \eqn{Null} and pi02 
#' will be generated randomly from uniform distribution.
#' @param nthread The number of threads used for deconvolution when OpenMP 
#' is available in the system. The default is the number of whole threads minus
#' one. In our no-OpenMP version, it is set to 1.
#'
#' @return
#' \item{pi}{Matrix of estimated proportions for each known component. The first
#' row corresponds to the proportion estimate of each sample for the first known
#' component (groupid = 1) and the second row corresponds to that for the second
#' known component (groupid = 2).}
#' \item{decovExpr}{A matrix of deconvolved expression profiles corresponding to
#' unknown (e.g tumor) component in mixed samples for a given subset of genes.
#' Each row corresponds to one gene and each column corresponds to one sample.} 
#' \item{decovMu}{Estimated \eqn{Mu} of log2-normal distribution for tumor
#' component.}
#' \item{decovSigma}{Estimated \eqn{Sigma} of log2-normal distribution for 
#' tumor component.}
#' \item{pi1}{An \eqn{My*I} matrix of estimated proportions for each 
#' iteration, where \eqn{I} is the number of iteration, for the first 
#' known component.}
#' \item{pi2}{An \eqn{My*I} matrix of estimated proportions for each 
#' iteration, where \eqn{I} is the number of iteration, for the second 
#' known component.}
#'
#' @author Zeya Wang, Wenyi Wang
#' 
#' @seealso http://bioinformatics.mdanderson.org/main/DeMixT
#'
#' @examples
#' # Example 1: simulated two-component data
#'   data(test.data.2comp)
#' # data.N1 <- SummarizedExperiment::assays(test.data.2comp$data.N1)[[1]]
#' # data.Y <- SummarizedExperiment::assays(test.data.2comp$data.Y)[[1]]
#' # inputdata <- cbind(data.N1, data.Y)
#' # groupid <- c(rep(1, ncol(data.N1)), rep(3, ncol(data.Y)))
#' # nspikein <- 0
#' # Optimum_KernelC(inputdata, groupid, 
#' #                 nspikein = nspikein, setting.pi = 0, 
#' #                 givenpi = rep(0, 2 * ncol(data.y)), 
#' #                 niter = 10, ninteg = 30, tol = 10^(-4))
#'                 
#' @references Wang Z, Cao S, Morris J S, et al. Transcriptome Deconvolution of 
#' Heterogeneous Tumor Samples with Immune Infiltration. iScience, 2018, 9: 451-460.
#' 
#' @keywords Optimum_KernelC
#' 
#' @export              
Optimum_KernelC <- function(
    inputdata, groupid, nspikein, setting.pi, givenpi, givenpiT, 
    niter, ninteg, tol, sg0 = 0.5^2, mu0 = 0.0, pi01 = NULL, pi02 = NULL,
    nthread = 1){

    #pi01 = NULL;pi02 = NULL
    if(!is.matrix(inputdata)) stop("argument inputdata must be a matrix");
    if(!is.vector(groupid)) stop("argument groupid must be a vector");
    if((FALSE %in% (unique(groupid) %in% c(1,2,3))) == TRUE) 
        stop("argument groupid must be a vector of 1, 2, 3")
    if((ninteg %% 1 != 0)||(ninteg <= 0)) 
        stop(print("ninteg can only be an positive integer"))
    
    if((niter %% 1 != 0)||(niter <= 0)) 
        stop(print("niter can only be an positive integer"))
    nsub <- as.integer(dim(inputdata)[2])   
    ## number of all the samples(normal + mixed)
    wgenes <- as.integer(dim(inputdata)[1]) 
    ## numebr of genes
    intx=sum(groupid == 3);                
    ## number of mixed samples
    intn=nsub-intx;                        
    ## number of normal samples
    nCid = length(unique(groupid))-1       
    ## number of known components
    groupid<-(as.array(groupid))
    
    ######################process for reorganizing groupid###########
    
    groupid1 = which(groupid == 1)
    if(nCid > 1) groupid2 = which(groupid == 2)
    groupidT = which(groupid == 3)
    
    sample.id = colnames(inputdata[,groupidT])
    gene.id = row.names(inputdata)
    
    if(nCid == 1){
        inputdata = cbind(inputdata[, groupid1], inputdata[, groupidT])
        groupid = c(rep(1, length(groupid1)), rep(3, length(groupidT)))
        }else if(nCid == 2){
        inputdata = cbind(inputdata[, groupid1], 
                            inputdata[, groupid2], inputdata[, groupidT])
        groupid = c(rep(1, length(groupid1)), 
                    rep(2, length(groupid2)), rep(3, length(groupidT)))
        }
    
    dataarray1 <- (as.array(matrix(inputdata, nrow = 1, byrow = FALSE))) 
    ## expression profile of the normal sample as an array
    ## Set given proportions
    if(setting.pi == 1){
        if(!is.vector(givenpi)) 
        stop("argument option must be a vector if pi is known");
        
        givenpi=(as.array(givenpi))
        givenpi3=(as.array(rep(0, intx)))
        }else if(setting.pi == 2){
        givenpi=(as.array(rep(0, nCid*intx)))
        givenpi3 = as.array(givenpiT)    
        }else if(setting.pi==0){
            givenpi=(as.array(rep(0, nCid*intx)))
            givenpi3=(as.array(rep(0, intx)))
            }else{
            stop("setting.pi argument must be set 0, 1 or 2")
            }
    
    if(nCid == 1){
        givenpi1 = givenpi[seq(1,intx)]
        givenpi2 = (as.array(rep(0, intx)))
        }else{
        givenpi1 = givenpi[seq(1,intx)]
        givenpi2 = givenpi[(intx+1):(2*intx)]
        }
    
    s0 = sg0 + 0.0
    m0 = mu0 + 0.0
    
    
    ## Initilize pi
    if(nCid == 1){
      if(is.null(pi01)){
        pi01 = pi02 = array(0 , c(1, intx))
        for(s in 1:intx){
          pi01[s] = runif(1, min = 0.01, max = 0.98)
        }
        colnames(pi01) <- paste('Sample', seq = seq(1:intx))
        rownames(pi01) <- paste('PiN1')
        cat("Initial of Proportions:\n")
        print(round(t(pi01),4))
      }
    }else if(nCid == 2){
      if(is.null(pi01) | is.null(pi02)){
        pi01 = pi02 =  array(0 , c(1, intx))
        for(s in 1:intx){
          pi01[s] = runif(1, min = 0.01, max = 0.97)
          pi02[s] = runif(1, min = 0.01, max = 0.98 - pi01[s])
        }
        # cat('Initial pi1 is \n', pi01, '\n')
        # cat('Initial pi2 is \n', pi02, '\n')
      }
    }
    
    rres <- .C("Tdemix", dataarray1, as.integer(groupid), as.integer(nsub), 
                as.integer(wgenes), as.integer(nspikein), as.integer(setting.pi), 
                pi01, pi02, givenpi1, givenpi2, 
                givenpi3, as.integer(nCid), as.integer(niter), 
                as.integer(ninteg), tol, as.integer(nthread), 
                s0, m0, rep(0, 2 * intx), 
                rep(0, intx * wgenes), rep(0, niter * wgenes), 
                rep(0, niter * wgenes), rep(0, niter * intx), 
                rep(0, niter * intx), rep(0,niter), 
                rep(0, intx * wgenes), rep(0, intx * wgenes))
    
    obj<-rres[[25]]
    
    if(sum(obj == 0)>1){
        niter1 <- which(obj==0)[1]-1
        }else{
        niter1 <- length(obj)
        }

    if (setting.pi != 1) {
        message('Objective function in each step: ')
        message(paste(obj[seq(1,niter1)], " "))
        message(' \n')
        }
    
    outcome1<-matrix(rres[[19]], ncol=intx, nrow=2, byrow=TRUE)
    outcome2<-matrix(rres[[20]], ncol=(intx), nrow=wgenes, byrow = TRUE)
    outcome3<-matrix(rres[[21]], ncol=niter, nrow=wgenes,byrow= TRUE)
    outcome4<-matrix(rres[[22]], ncol=niter, nrow=wgenes,byrow=TRUE)
    outcome5<-matrix(rres[[23]], ncol=niter, nrow=intx,byrow=TRUE)
    outcome6<-matrix(rres[[24]], ncol=niter, nrow=intx,byrow=TRUE)
    outcome21<-matrix(rres[[26]], ncol=(intx), nrow=wgenes, byrow = TRUE)
    outcome22<-matrix(rres[[27]], ncol=(intx), nrow=wgenes, byrow = TRUE)

#                 if(setting.pi == 1){
# 	            print(paste0('Totally ', round(100*sum(outcome4[,niter]>99.99)/wgenes,3),'% genes estimated touch the optimization bound'))
#                 }else{
#                  print(paste0('Totally ', round(100*sum(outcome4[,niter]>24.99)/wgenes,3),'% genes estimated touch the optimization bound'))
#                 }
                    
    outcome4 <- sqrt(outcome4)
    return(list(obj_val = obj[niter1], obj_iter=obj, pi = outcome1, decovExprT = outcome2, 
                decovExprN1 = outcome21, decovExprN2 = outcome22, 
                decovMu = as.matrix(outcome3[,seq(1,niter1)]), 
                decovSigma = as.matrix(outcome4[,seq(1,niter1)]), 
                pi1 = as.matrix(outcome5[,seq(1,niter1)]), 
                pi2 = as.matrix(outcome6[,seq(1,niter1)])))
}
