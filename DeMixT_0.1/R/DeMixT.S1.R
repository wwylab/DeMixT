###DeMix step I function###
DeMixT.S1 <- function(inputdata, groupid, niter = 10, ninteg = 50, filter.option = 1, filter.criteria1 = c(0.5,0.5), filter.criteria2 = c(250,250), filter.criteria3 = 0.25, if.filter = FALSE,tol=10^(-5), sg0 = 0.5^2, mu0 = 0.0, nthread=-1){
##filter.option:1, we remove all the genes containing zero count in any sample; 2, we add 1 to the original inputdatamat to kill all the zero
    core.num <- round(detectCores())-1
    if(nthread == -1) nthread = core.num
	if(filter.option == 1){
		inputdata <- inputdata[apply(inputdata, 1, FUN = function(x) sum(x <= 0) == 0), ]
	}else if(filter.option == 2){
	    inputdata <- inputdata + 1
	}else{
	stop("the argument filter.option can only be 1 or 2")
	}
	nCid = length(unique(groupid))-1       ## number of known components
	sample.id = colnames(inputdata[,groupid == 3])
	
	if (nCid == 1){
	  if.stage = FALSE
	}else{
		if(if.filter == FALSE){
			if.stage = FALSE
		}else{
			if.stage = TRUE
		}
	}
						  
	if(if.stage == FALSE){
	if(if.filter == TRUE){
	if(nCid == 1){
	   inputdatans = apply(log2(inputdata[, groupid == 1]), 1, sd)
	   id1 <- (inputdatans < filter.criteria1[1])
		if(sum(id1)<20) stop("the threshold of standard variation for filtering gene is too stringent please provide a larger threshold ")
	   inputdatamat1 = inputdata[id1,]
	   inputdatamat1nm <- apply(inputdatamat1[,groupid == 1], 1, mean)	 
	   inputdatamat1ym <- apply(inputdatamat1[,groupid == 3], 1, mean)	 
	   inputdata1r <- inputdatamat1ym/inputdatamat1nm
		
#filter.criteria2 can be an integer or quantile
		if((filter.criteria2[1] > 1)&(filter.criteria2[1]%%1 == 0)){
			id2 <- order(inputdata1r, decreasing = T)
			id2 <- id2[1:filter.criteria2[1]]
		}else if((filter.criteria2[1] < 1)&(filter.criteria2[1] > 0)){
			id2 <- (inputdata1r > quantile(inputdata1r, probs = 1 - filter.criteria2[1])) ##subset the final genes
		}else{
			stop("filter.criteria has too be an integer or a percentage between 0 and 1")
		}
	   
	   if(sum(id2)<20) stop("there are too few genes please relax threshold for filtering ")
	   inputdatamat2 <- inputdatamat1[id2,] 
	}else if(nCid == 2){
        if(length(filter.criteria1) < 2) stop("Filter criteria 1 should be a vector of length 2 for three component deconvovlution")
		inputdatans1 = apply(log2(inputdata[, groupid == 1]), 1, sd)
		inputdatans2 = apply(log2(inputdata[, groupid == 2]), 1, sd)
		id1 <- ((inputdatans1 < filter.criteria1[1])&(inputdatans2 < filter.criteria1[2]))
		if(sum(id1)<20) stop("the threshold of standard variation for filtering gene is too stringent please provide a larger threshold ")
		inputdatamat1 = inputdata[id1,]
		inputdata1s <- apply(inputdatamat1[,groupid == 3], 1, sd)	 
		id2 <- (inputdata1s > quantile(inputdata1s, probs = 1 - filter.criteria2[1])) ##subset the final genes
		if(sum(id2)<20) stop("there are too few genes please relax threshold for filtering ")
		inputdatamat2 <- inputdatamat1[id2,] 
	}else{
	 stop("The expected number of components should be 2 or 3")
	}
	}else{
	inputdatamat2 = inputdata
	}
	   res <- Optimum.KernelC(inputdatamat2, groupid, nhavepi=0, givenpi=rep(0, 2*sum(groupid==3)), givenpiT=rep(0,sum(groupid==3)), niter=niter, ninteg=ninteg, tol=tol, sg0=sg0,mu0=mu0, nthread=nthread)
	}else{
		if(nCid==1) stop("if.stage can only be applied to a three component deconvolution")
        ############two-step filtering for three component#########
		print("Fitering stage 1 starts")
        # step 1
		inputdatan1m <- apply(log2(inputdata[,groupid == 1]), 1, mean)
		inputdatan2m <- apply(log2(inputdata[,groupid == 2]), 1, mean)
		inputdatan1s <- apply(log2(inputdata[,groupid == 1]), 1, sd)
		inputdatan2s <- apply(log2(inputdata[,groupid == 2]), 1, sd)
		id1 <- ((abs(inputdatan1m - inputdatan2m) < filter.criteria3)
			   &(inputdatan1s < filter.criteria1[1])&(inputdatan2s < filter.criteria1[2]))
		if(sum(id1)<20) stop("the threshold of filer.criteria1 and filter.criteria3 for filtering gene is too stringent please provide a larger threshold ")
        # step 2
		inputdatamat1 <- inputdata[id1,]
		inputdatamat1nm <- apply(inputdatamat1[,(groupid == 1)|(groupid == 2)], 1, mean)
		inputdatamat1ym <- apply(inputdatamat1[,(groupid == 3)], 1, mean)
#inputdata1r <- inputdatamat1ym/inputdatamat1nm
		inputdata1r <- apply(inputdatamat1[,groupid == 3], 1, sd)
		
#filter.criteria2 can be an integer or quantile
		if((filter.criteria2[1] > 1)&(filter.criteria2[1]%%1 == 0)){
			id2 <- order(inputdata1r, decreasing = T)
			id2 <- id2[1:filter.criteria2[1]]
		}else if((filter.criteria2[1] < 1)&(filter.criteria2[1] > 0)){
			id2 <- (inputdata1r > quantile(inputdata1r, probs = 1 - filter.criteria2[1])) ##subset the final genes
		}else{
		  stop("filter.criteria has too be an integer or a percentage between 0 and 1")
		}
		
		if(sum(id2)<20) stop("there are too few genes for filtering stage 1 please relax threshold for filtering ")		
		inputdatamat2 <- inputdatamat1[id2,]
        # step 3
		cnvgroup <- groupid; cnvgroup[groupid == 2] = 1; # inputdata as 2-component
        # Run this function
		res1 <- Optimum.KernelC(inputdatamat2, cnvgroup, nhavepi=0, givenpi=rep(0, sum(groupid==3)), givenpiT=rep(0,sum(groupid==3)), niter=niter, ninteg=ninteg, tol=tol, mu0=mu0, sg0=sg0, nthread=nthread)
		print("Filtering stage 1 over")
		fixed.piT <- 1 - as.numeric(res1$pi[1,])
        #################################second stage##########################################
#################################second stage##########################################
 		print("Filtering stage 2 starts")
        #step 1
		id3 <- ((inputdatan1s < filter.criteria1[1])&(inputdatan2s < filter.criteria1[2]))
        if(sum(id3)<20) stop("there are too few genes for filtering stage 2 please relax threshold for filtering ")
        inputdatamat1 <- inputdata[id3,]
		inputdatan1m <- apply(log2(inputdatamat1[,groupid == 1]), 1, mean)
		inputdatan2m <- apply(log2(inputdatamat1[,groupid == 2]), 1, mean)
		inputdata1d <- abs(inputdatan1m - inputdatan2m)
#filter.criteria2 can be an integer or quantile
		if((filter.criteria2[2] > 1)&(filter.criteria2[2]%%1 == 0)){
			id4 <- order(inputdata1d, decreasing = T)
			id4 <- id4[1:(2*filter.criteria2[2])]
		}else if((filter.criteria2[2] < 1)&(filter.criteria2[2] > 0)){
			id4 <- (inputdata1d > quantile(inputdata1d, probs = 1 - 2*filter.criteria2[2]))
		}else{
			stop("filter.criteria has too be an integer or a percentage between 0 and 1")
			
		}
	
		inputdatamat2 <- inputdatamat1[id4,]
		if(sum(id4)<40) stop("there are too few genes please relax threshold for filtering ")
		inputdata1s <- apply(inputdatamat2[,groupid == 3], 1, sd)
		id5 <- (inputdata1s > quantile(inputdata1s, probs = 0.5))
		inputdatamat3 <- inputdatamat2[id5,]
        #step 2
		res <- Optimum.KernelC(inputdatamat3, groupid, nhavepi = 2, givenpi=rep(0, 2*sum(groupid==3)), givenpiT = fixed.piT, niter=niter, ninteg=ninteg, tol=tol, mu0=mu0, sg0=sg0, nthread=nthread)
		print("Filtering stage 2 over")
	}

	  pi <- t(as.matrix(res$pi[1,])); row.names(pi)[1] = "pi1";
      
	  pi1 <- as.matrix(res$pi1); colnames(pi1) = as.character(1:ncol(pi1)); row.names(pi1) = sample.id
	  pi_iteration <- list(); pi_iteration[[1]] <- pi1
	  if(nCid == 2) {
          pi <- rbind(pi, res$pi[2,])
          row.names(pi)[2] = "pi2"
		  pi2 <- as.matrix(res$pi2);colnames(pi2) = as.character(1:ncol(pi2)); row.names(pi2) = sample.id
          pi_iteration[[2]] <- pi2
				}
     colnames(pi) = sample.id
	
	 return(list(pi = pi, pi_iteration = pi_iteration))
      						  
						  
}
