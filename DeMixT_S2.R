DeMixT_S2 <- function(data.Y, data.comp1, data.comp2 = NULL, givenpi, nbin = 50, 
          nthread = parallel::detectCores() - 1) 
{
  filter.out = TRUE
  filter.option = 1
  data.Y <- SummarizedExperiment::assays(data.Y)[[1]]
  data.comp1 <- SummarizedExperiment::assays(data.comp1)[[1]]
  if (!is.null(data.comp2)){
    data.comp2 <- SummarizedExperiment::assays(data.comp2)[[1]]
  }  
    
  if (is.null(rownames(data.Y))) {
    rownames(data.Y) <- as.character(seq(1, nrow(data.Y)))

    rownames(data.comp1) <- as.character(seq(1, nrow(data.comp1)))
    if(!is.null(data.comp2)){
      rownames(data.comp2) <- as.character(seq(1, nrow(data.comp2)))
    }
  }
  if (is.null(colnames(data.Y))) {
    colnames(data.Y) <- as.character(seq(1, ncol(data.Y)))
    colnames(data.comp1) <- as.character(seq(1, ncol(data.comp1)))
    if(!is.null(data.comp2)){
      colnames(data.comp2) <- as.character(seq(1, ncol(data.comp2)))
    }
  }
  sample.id <- colnames(data.Y)
  if (is.null(data.comp2)) {
    inputdata <- cbind(data.comp1, data.Y)
    groupid <- c(rep(1, ncol(data.comp1)), rep(3, ncol(data.Y)))
  }else {
    inputdata <- cbind(data.comp1, data.comp2, data.Y)
    groupid <- c(rep(1, ncol(data.comp1)), rep(2, ncol(data.comp2)), 
                 rep(3, ncol(data.Y)))
  }
  if (filter.option == 1) {
    index <- apply(inputdata, 1, function(x) sum(x <= 0) == 
                     0)
    inputdata <- inputdata[index, ]
    data.comp1 <- data.comp1[index, ]
    if (!is.null(data.comp2)) 
      data.comp2 <- data.comp2[index, ]
    data.Y <- data.Y[index, ]
  }else if (filter.option == 2) {
    inputdata <- inputdata + 1
    data.comp1 <- data.comp1 + 1
    if (!is.null(data.comp2)) 
      data.comp2 <- data.comp2 + 1
    data.Y <- data.Y + 1
  }else {
    stop("The argument filter.option can only be 1 or 2")
  }
  
  if (is.null(data.comp2)) {
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.comp1, 
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.comp1, 1, function(x) 
        length(unique(x)) > 1), ]
    }
  }else {
    if (dim(inputdata)[1] == 1) {
      inputdata <- t(as.matrix(inputdata[apply(data.comp1, 
        1, function(x) length(unique(x)) > 1) & apply(data.comp2, 
        1, function(x) length(unique(x)) > 1), ]))
    }
    else {
      inputdata <- inputdata[apply(data.comp1, 1, function(x) 
        length(unique(x)) > 1) & apply(data.comp2, 
        1, function(x) length(unique(x)) > 1), ]
    }
  }

  gene.id <- rownames(inputdata)
  inputdata <- as.matrix(inputdata)
  res <- Optimum_KernelC(inputdata, groupid, nhavepi = 1,  nspikein = 0, 
    givenpi = givenpi, givenpiT = rep(0, ncol(data.Y)), niter = 1, ninteg = nbin, 
    tol = 1e-05, nthread = nthread)

  
  decovExprT <- res$decovExprT
  decovExprT[which(decovExprT < 0,arr.ind = T)] <- 0
  colnames(decovExprT) <- sample.id
  row.names(decovExprT) <- gene.id
  decovExprN1 <- res$decovExprN1
  decovExprN1[which(decovExprN1 < 0,arr.ind = T)] <- 0
  colnames(decovExprN1) <- sample.id
  row.names(decovExprN1) <- gene.id
  decovExprN2 <- res$decovExprN2
  decovExprN2[which(decovExprN2 < 0,arr.ind = T)] <- 0
  colnames(decovExprN2) <- sample.id
  row.names(decovExprN2) <- gene.id
  
  
  if (is.null(data.comp2)) {
    
    decovMuN <- rowMeans(log2(data.comp1[gene.id, ]))
    decovSigmaN <- apply(log2(data.comp1[gene.id, ]), 1 , sd)
    
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
    
    decovMuN1 <- rowMeans(log2(data.comp1[gene.id,]))
    
    decovSigmaN1 <- apply(log2(data.comp1[gene.id,]), 1, sd)
    
    decovMuN2 <- rowMeans(log2(data.comp2[gene.id,]))
    
    decovSigmaN2 <- apply(log2(data.comp2[gene.id,]), 1, sd)
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