###DeMixT step II function
DeMixT.S2 <- function(inputdata, groupid, givenpi,ninteg = 50, filter.out = TRUE, filter.option = 1, nthread=-1){
    core.num <- round(detectCores())-1
    if(nthread == -1) nthread = core.num
	sample.id <- colnames(inputdata[, groupid == 3])
    if(is.null(row.names(inputdata))) row.names(inputdata) <- 1:nrow(inputdata)
	if(filter.option == 1){
		inputdata <- inputdata[apply(inputdata, 1, FUN = function(x) sum(x <= 0) == 0), ]
        gene.id <- row.names(inputdata)

	}else if(filter.option == 2){
	    inputdata <- inputdata + 1
		gene.id <- row.names(inputdata)
	}else{
		stop("the argument filter.option can only be 1 or 2")
	}
	
    res <- Optimum.KernelC(inputdata, groupid, nhavepi = 1, givenpi=givenpi, givenpiT = rep(0,sum(groupid==3)), niter=1, ninteg=ninteg, tol=1e-5, nthread=nthread)
    
    decovExprT <- res$decovExprT; colnames(decovExprT) = sample.id; row.names(decovExprT) = gene.id
    decovExprN1 <- res$decovExprN1; colnames(decovExprN1) = sample.id; row.names(decovExprN1) = gene.id
    decovExprN2 <- res$decovExprN2; colnames(decovExprN2) = sample.id; row.names(decovExprN2) = gene.id
    
    if(length(unique(groupid)) == 2){
    decovMuN <- apply(log2(inputdata[,groupid == 1]), 1, mean)
    decovSigmaN <- apply(log2(inputdata[,groupid == 1]), 1, sd)
    decovMu <- cbind(decovMuN, res$decovMu); colnames(decovMu) = c('MuN1', 'MuT'); row.names(decovMu) = gene.id
    decovSigma <- cbind(decovSigmaN, res$decovSigma); colnames(decovSigma) = c('SigmaN1', 'SigmaT'); row.names(decovSigma) = gene.id
    ## filter out genes with muN - muT > 4
    filter.in <- (decovMu[,1] - decovMu[,2] <= 4.0)
    if(filter.out = FALSE)  filter.in <- rep(TRUE, dim(inputdata)[1])
    }else{
    decovMuN1 <- apply(log2(inputdata[,groupid == 1]), 1, mean)
    decovSigmaN1 <- apply(log2(inputdata[,groupid == 1]), 1, sd)
    decovMuN2 <- apply(log2(inputdata[,groupid == 2]), 1, mean)
    decovSigmaN2 <- apply(log2(inputdata[,groupid == 2]), 1, sd)
    decovMu <- cbind(decovMuN1, decovMuN2, res$decovMu); colnames(decovMu) = c('MuN1', 'MuN2','MuT'); row.names(decovMu) = gene.id
    decovSigma <- cbind(decovSigmaN1, decovSigmaN2, res$decovSigma); colnames(decovSigma) = c('SigmaN1','SigmaN2','SigmaT'); row.names(decovSigma) = gene.id
    filter.in <- rep(TRUE, dim(inputdata)[1])
    }
    
    
    return(list(decovExprT = decovExprT[filter.in, ], decovExprN1 = decovExprN1[filter.in, ], decovExprN2 = decovExprN2[filter.in, ], decovMu = decovMu[filter.in, ], decovSigma = decovSigma[filter.in, ]))
    
}



