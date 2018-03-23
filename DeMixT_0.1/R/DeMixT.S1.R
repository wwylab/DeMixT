require(parallel)

DeMixT.S1 <- function(data.Y, data.comp1, data.comp2 = NULL, 
                      niter = 10, nbin = 50, 
                      if.filter = FALSE, ngene.selected.for.pi = 250, mean.diff.in.CM = 0.25,
                      tol = 10^(-5), nthread = detectCores() - 1) {
  
  filter.option = 1
  filter.sd = 0.5
  
  ## index gene and sample names
  if (is.null(rownames(data.Y))) {
    rownames(data.Y) <- as.character(1:nrow(data.Y))
    
  }
  if (is.null(colnames(data.Y))) {
    colnames(data.Y) <- as.character(1:ncol(data.Y))
  }
  
  ## combine datasets
  if (is.null(data.comp2)) { # two-component
    inputdata <- cbind(data.comp1, data.Y)
    groupid <- c(rep(1, ncol(data.comp1)), rep(3, ncol(data.Y)))
  } else { # three-component
    inputdata<- cbind(data.comp1, data.comp2, data.Y)
    groupid <- c(rep(1, ncol(data.comp1)), rep(2, ncol(data.comp2)), rep(3, ncol(data.Y)))
  }
  
  ## filter.option: 1 - remove genes containing zero; 2 - add 1 to to kill zeros
  if (filter.option == 1) {
    inputdata <- inputdata[apply(inputdata, 1, function(x) sum(x <= 0) == 0), ]
  } else if (filter.option == 2) {
    inputdata <- inputdata + 1
  } else {
    stop("The argument filter.option can only be 1 or 2")
  }
  
  ## filter out genes with constant value across all samples
  inputdata < ifelse(is.null(data.comp2), inputdata[apply(data.comp1, 1, function(x) length(unique(x)) > 1), ], 
                     inputdata[apply(data.comp1, 1, function(x) length(unique(x)) > 1) & apply(data.comp2, 1, function(x) length(unique(x)) > 1), ])
  
  #
  filter2 <- function(inputdata1r, ngene.selected.for.pi, n = 1) {
    if ((ngene.selected.for.pi > 1) & (ngene.selected.for.pi %% 1 == 0)) { # ngene.selected.for.pi can be integer or quantile
      id2 <- order(inputdata1r, decreasing = TRUE)
      id2 <- id2[1:min(n * ngene.selected.for.pi, length(inputdata1r))]
    } else if ((ngene.selected.for.pi < 1) & (ngene.selected.for.pi > 0)) {
      id2 <- (inputdata1r > quantile(inputdata1r, probs = 1 - n * ngene.selected.for.pi))
    } else {
      stop("The argument ngene.selected.for.pi can only be an integer or a percentage between 0 and 1")
    }
    if (sum(id2) < 20) stop("there are too few genes for filtering stage 1 please relax threshold for filtering ")		
    return(inputdatamat1[id2,])
  }
  
  ## case 1
  if(if.filter == FALSE){
    gene.name <- rownames(inputdata)
    res <- Optimum.KernelC(inputdata, groupid, nhavepi = 0, givenpi = rep(0, 2 * ncol(data.y)), givenpiT = rep(0, ncol(data.y)), 
                           niter = niter, ninteg = nbin, tol = tol)
  } else {
    
  }
    
  ## case 2
  if (if.filter == TRUE & is.null(data.comp2)) {
    # step 1
    inputdatans <- apply(log2(data.comp1), 1, sd)
    id1 <- (inputdatans < filter.sd)
    if (sum(id1) < 20) stop("The threshold of standard variation is too stringent. Please provide a larger threshold. ")
    inputdatamat1 <- inputdata[id1, ]
    # step 2
    inputdatamat1nm <- apply(inputdatamat1[, groupid == 1], 1, mean)
    inputdatamat1ym <- apply(inputdatamat1[, groupid == 3], 1, mean)	 
    inputdata1r <- inputdatamat1ym / inputdatamat1nm
    inputdata2 <- filter2(inputdata1r, ngene.selected.for.pi)
    # 
    gene.name <- rownames(inputdata2)
    res <- Optimum.KernelC(inputdata2, groupid, nhavepi = 0, givenpi = rep(0, 2 * ncol(data.y)), givenpiT = rep(0, ncol(data.y)), 
                           niter = niter, ninteg = nbin, tol = tol)
  }
  
  ## case 3: two-stage filtering
  if(if.filter == TRUE & !is.null(data.comp2)) {

    cat("Fitering stage 1 starts\n")
    # step 1
    inputdatan1m <- apply(log2(inputdata[,groupid == 1]), 1, mean)
    inputdatan2m <- apply(log2(inputdata[,groupid == 2]), 1, mean)
    inputdatan1s <- apply(log2(inputdata[,groupid == 1]), 1, sd)
    inputdatan2s <- apply(log2(inputdata[,groupid == 2]), 1, sd)
    id1 <- ((abs(inputdatan1m - inputdatan2m) < mean.diff.in.CM) & (inputdatan1s < filter.sd) & (inputdatan2s < filter.sd))
    if(sum(id1) < 20) stop("The thresholds of standard variation and mean difference are too stringent. Please provide larger thresholds")
    inputdatamat1 <- inputdata[id1, ]
    
    # step 2
    inputdatamat1nm <- apply(inputdatamat1[ ,(groupid == 1) | (groupid == 2)], 1, mean)
    inputdatamat1ym <- apply(inputdatamat1[ ,(groupid == 3)], 1, mean) # use "inputdata1r <- inputdatamat1ym/inputdatamat1nm"
    inputdata1r <- apply(inputdatamat1[, groupid == 3], 1, sd)
    inputdatamat2 <- filter2(inputdata1r, ngene.selected.for.pi)
    # step 3
    cnvgroup <- groupid; cnvgroup[groupid == 2] <- 1 # combine 3-component to 2-component
    #
    res1 <- Optimum.KernelC(inputdatamat2, cnvgroup, nhavepi = 0, givenpi = rep(0, ncol(data.y)), givenpiT = rep(0,ncol(data.y)), 
                            niter = niter, ninteg = nbin, tol = tol)
    fixed.piT <- 1 - as.numeric(res1$pi[1, ])
    cat("Filtering stage 1 is finished\n")

    cat("Filtering stage 2 starts\n")
    # step 1
    id3 <- ((inputdatan1s < filter.sd) & (inputdatan2s < filter.sd))
    if(sum(id3) < 20) stop("The thresholds of standard variation and mean difference are too stringent. Please provide larger thresholds")
    # step 2
    inputdatamat1 <- inputdata[id3, ]
    inputdatan1m <- apply(log2(inputdatamat1[, groupid == 1]), 1, mean)
    inputdatan2m <- apply(log2(inputdatamat1[, groupid == 2]), 1, mean)
    inputdata1d <- abs(inputdatan1m - inputdatan2m)
    inputdatamat2 <- filter2(inputdata1d, ngene.selected.for.pi, 2)
    # step 3
    inputdata1s <- apply(inputdatamat2[, groupid == 3], 1, sd)
    id5 <- (inputdata1s > quantile(inputdata1s, probs = 0.5))
    inputdatamat3 <- inputdatamat2[id5, ]
    #
    gene.name <- rownames(inputdatamat3)
    res <- Optimum.KernelC(inputdatamat3, groupid, nhavepi = 2, givenpi = rep(0, 2 * ncol(data.Y)), givenpiT = fixed.piT, 
                           niter = niter, ninteg = nbin, tol = tol)
    cat("Filtering stage 2 is finished")
  }
  
  pi <- t(as.matrix(res$pi[1, ])); row.names(pi)[1] <-  "pi1"
  pi1 <- as.matrix(res$pi1); colnames(pi1) <- as.character(1:ncol(pi1)); row.names(pi1) = colnames(data.Y)
  pi.iter <- array(t(res$pi1), dim <- c(ncol(res$pi1), nrow(res$pi1), 1))
  
  if(!is.null(data.comp2)) {
    pi <- rbind(pi, res$pi[2,]); row.names(pi)[2] <-  "pi2"
    pi2 <- as.matrix(res$pi2); colnames(pi2) <- as.character(1:ncol(pi2)); row.names(pi2) = colnames(data.Y)
    pi.iter <- array(t(rbind(res$pi1, res$pi2)), dim <- c(ncol(res$pi1), nrow(res$pi1), 2))
  }
  
  colnames(pi) <- colnames(data.Y)
  
  return(list(pi = pi, pi.iter = pi.iter, gene.name = gene.name))
}
