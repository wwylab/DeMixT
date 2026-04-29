

#' @description Detect suspicious samples by hierarchical clustering
#' @param count.matrix raw expression count matrix: gene x sample, rownames are genes, columnnames are sample ids
#' @param normal.id normal sample ids
#' @param tumor.id tumor sample ids
#' @return list object
#' 

detect_suspicious_sample_by_hierarchical_clustering <- function(count.matrix, normal.id, tumor.id){
    normal.exp <- count.matrix[, match(normal.id, colnames(count.matrix))]
    tumor.exp <- count.matrix[, match(tumor.id, colnames(count.matrix))]
    gene.normal.mean <- apply(normal.exp, 1, function(x) mean(x))
    gene.tumor.mean <- apply(tumor.exp, 1, function(x) mean(x))
    
    gene.normal.index <- gene.normal.mean > quantile(gene.normal.mean, 0.25)
    gene.tumor.index <- gene.tumor.mean > quantile(gene.tumor.mean, 0.25)
    sTable = count.matrix[gene.normal.index | gene.tumor.index, ]
    rank.test.pvalue = NULL
    
    labels = c(rep("o", length(normal.id)), rep(".", length(tumor.id)))
    pch = c(rep(1, length(normal.id)), rep(20, length(tumor.id)))
    
    y <- rep("Normal", dim(count.matrix)[2])
    y[colnames(count.matrix) %in% tumor.id] <- "Tumor"
    names(y) <- colnames(count.matrix)
    
    test.pvalue <- apply(sTable, 1, function(x) wilcox.test(x ~ y)$p.value)
    sorted.pvalue = sort.int(test.pvalue, index.return = TRUE)
    #select the top 1000 genes with smallest p-value
    top.gene.index = sorted.pvalue$ix[1:1000]
    
    #pca analysis
    principal.res = prcomp(t(sTable[top.gene.index, ]),
                           retx = T,
                           center = T,
                           scale = T)
    
    top2.pcs <- principal.res$x[, 1:2]
    top2.pcs.dis <-  dist(top2.pcs, method = "euclidean")
    top2.pcs.hclust <- hclust(top2.pcs.dis, method = "ward.D2")
    top2.pcs.two.clusters <- cutree(top2.pcs.hclust, k = 2)
    
    labels <- y[top2.pcs.hclust$order]
    hc_labels <- list()
    hc_labels[["label"]] <- labels
    hc_labels[["cluster"]] <- top2.pcs.two.clusters
    hc_labels[["hc"]] <- top2.pcs.hclust
    
    return(hc_labels)
}

plot_sd  <- function(inputdata, label){
    n.normal <- table(label)[1]
    sdn.obs <- apply(log2(inputdata[,1:n.normal]), 1, sd)
    sdm.obs <- apply(log2(inputdata[,-c(1:n.normal)]), 1, sd)
    par(mfrow = c(1,2))
    plot(sdn.obs, 
         xlab = 'Genes', ylab = 'Standard Deviation', main = 'Normal Reference')
    plot(sdm.obs, 
         xlab = 'Genes', ylab = 'Standard Deviation', main = 'Mixed Tumor')
    par(mfrow = c(1,1))
}

subset_sd <- function(inputdata, label, 
                      cutoff_normal = c(0.1, 0.6), 
                      cutoff_tumor = c(0.2, 0.8)){
    n.normal <- table(label)[1]
    sdn.obs <- apply(log2(inputdata[,1:n.normal]), 1, sd)
    sdm.obs <- apply(log2(inputdata[,-c(1:n.normal)]), 1, sd)
    indx <- which(sdn.obs > cutoff_normal[1] & sdn.obs < cutoff_normal[2] &
                      sdm.obs > cutoff_tumor[1] & sdm.obs < cutoff_tumor[2])
    inputdata <- inputdata[indx, ]
}


# PRAD.filter = subset_sd(PRAD.filter, label = label, 
#                         cutoff_normal = c(0.1, 0.7),
#                         cutoff_tumor = c(0, 0.8))
# cat('Number of genes: ', dim(PRAD.filter)[1], '\n')

#' @description Make a PCA plot for the expression data and check if normal samples and mixed tumor samples are well separated based on PCA


plot_dim <- function(inputdata, label, reduction = 'pca',
                     legend.position = 'bottomleft',
                     legend.cex = 1.2){
    N.label = length(unique(label))
    inputdata[which(inputdata == 0, arr.ind = T)] = 1
    ## PCA dimension reduction
    res = 0
    if(reduction == 'pca'){
        res <- prcomp(t(log2(inputdata)), center = T, scale. = T)
        ## PC variance
        pct1 = round(res$sdev[1]/sum(res$sdev), 3)*100
        pct2 = round(res$sdev[2]/sum(res$sdev), 3)*100
        ## Dimension
        dim1 = res$x[,1]; dim2 = res$x[,2]
        ## x/ylab text
        xlab.text = paste("First Comp: ", as.character(pct1), "% variance", sep="")
        ylab.text = paste("Second Comp: ", as.character(pct2), "% variance", sep="")   
        main.text = 'PCA'
    }
    ## Plot
    plot(dim1, dim2,  cex = 1, 
         xlab=xlab.text, ylab=ylab.text, main = main.text,
         col=rainbow(N.label)[as.numeric(factor(label))],
         pch = as.numeric(factor((label))), lwd=1.5)
    abline(h=0, v=0, col="brown", lty=2)
    abline(h=0, v=0, col="brown", lty=2)
    center1<-tapply(dim1, factor(label), mean)
    center2<-tapply(dim2, factor(label), mean)
    ## add lines and centers
    for (ii in 1:length(center1)) {
        groupi<-cbind(dim1,dim2)[as.numeric(factor(label))==ii, 1:2]
        if (class(groupi)[1] =="matrix") {
            for (j in (1:nrow(groupi))) {
                segments(groupi[j,1], groupi[j,2], center1[ii], center2[ii], col = rainbow(N.label)[ii], lwd = 0.5)
            }
        }else {
            segments(groupi[1], groupi[2], center1[ii], center2[ii], col = rainbow(N.label)[ii], lwd = 0.5)
        }
    }
    points(center1, center2, pch = 7, lwd = 1.5, col = rainbow(N.label))
    legend(legend.position,legend=names(table(factor(label))), 
           text.col=rainbow(N.label), pch=1:N.label, 
           col=rainbow(N.label), lty=1, cex = legend.cex)
}

#' @title Quantile_Normalization_Scale
#' @description perform 75th quantile scale normalization on the input data

Quantile_Normalization_Scale <- function(Count.matrix){
    newt <- Count.matrix
    colnames(newt)=NULL
    rownames(newt)=NULL
    
    designs=c(rep("0", dim(Count.matrix)[2]))
    seqData=newSeqCountSet(as.matrix(newt), designs)
    
    # Quantile normaliszation/total/median   ###try different normalization method###
    seqData=estNormFactors(seqData, "quantile")
    k3=seqData@normalizationFactor
    mk3=median(k3)
    k3=k3/mk3
    
    temp<-newt
    
    for(i in 1:ncol(newt)){
        temp[,i] = temp[,i]/k3[i]
    }
    Count.matrix.normalized<-temp
    colnames(Count.matrix.normalized)<-colnames(Count.matrix)
    rownames(Count.matrix.normalized)<-rownames(Count.matrix)
    
    return(Count.matrix.normalized)
}

#' @title subset_sd_gene_remaining
#' @description find the cutoffs to filter out genes with large standard deviations of expressions in both normal and tumor samples
subset_sd_gene_remaining <- function(inputdata, label, 
                                     cutoff_normal_range = c(0.2, 0.6), 
                                     cutoff_tumor_range = c(0.2, 0.8),
                                     cutoff_step = 0.2){
    
    n.normal <- table(label)[1]
    sdn.obs <- apply(log2(inputdata[,1:n.normal]), 1, sd)
    sdm.obs <- apply(log2(inputdata[,-c(1:n.normal)]), 1, sd)
    
    cutoff_normal_range <- seq(cutoff_normal_range[1], cutoff_normal_range[2], by = cutoff_step)
    cutoff_tumor_range <- seq(cutoff_tumor_range[1], cutoff_tumor_range[2], by = cutoff_step)
    
    num_gene_remaining_different_cutoffs <- NULL
    
    for(cutoff_normal in cutoff_normal_range){
        if(cutoff_normal == cutoff_normal_range[1]) next
        
        for(cutoff_tumor in cutoff_tumor_range){
            if(cutoff_tumor == cutoff_tumor_range[1]) next
            
            indx <- which(sdn.obs > cutoff_normal_range[1] & sdn.obs < cutoff_normal &
                              sdm.obs > cutoff_tumor_range[1] & sdm.obs < cutoff_tumor)
            
            num_gene_remaining_different_cutoffs <- rbind(num_gene_remaining_different_cutoffs, 
                                                          data.frame(normal.cutoff.low=cutoff_normal_range[1],
                                                                     normal.cutoff.high=cutoff_normal,
                                                                     tumor.cutoff.low=cutoff_tumor_range[1],
                                                                     tumor.cutoff.high=cutoff_tumor,
                                                                     num.gene.remaining=length(indx)))
            
        }
    }
    return(num_gene_remaining_different_cutoffs)
}


DeMixT_preprocessing <- function(count.matrix, label, 
                                 cutoff_normal_range=c(0.1, 1.0), 
                                 cutoff_tumor_range=c(0, 2.5), 
                                 cutoff_step=0.2){
    
    num_gene_remaining_different_cutoffs <- subset_sd_gene_remaining(count.matrix, label, 
                                                                     cutoff_normal_range, 
                                                                     cutoff_tumor_range,
                                                                     cutoff_step)
    
    
    num_gene_remaining_different_cutoffs_filter <- num_gene_remaining_different_cutoffs[which(num_gene_remaining_different_cutoffs$num.gene.remaining - 9000 > 0), ]
    num_gene_remaining_different_cutoffs_filter <- num_gene_remaining_different_cutoffs_filter[order(num_gene_remaining_different_cutoffs_filter$num.gene.remaining), ]
    sd_cutoff_normal <- c(num_gene_remaining_different_cutoffs_filter$normal.cutoff.low[1], num_gene_remaining_different_cutoffs_filter$normal.cutoff.high[1])
    sd_cutoff_tumor <- c(num_gene_remaining_different_cutoffs_filter$tumor.cutoff.low[1], num_gene_remaining_different_cutoffs_filter$tumor.cutoff.high[1])
    
    count.matrix <- subset_sd(inputdata = count.matrix, label = label, cutoff_normal = sd_cutoff_normal, cutoff_tumor = sd_cutoff_tumor)
    
    count.matrix = Quantile_Normalization_Scale(count.matrix)
    
    preprocessing_output <- list()
    
    preprocessing_output[['count.matrix']] <- count.matrix
    
    preprocessing_output[["sd_cutoff_normal"]] <- sd_cutoff_normal
    preprocessing_output[["sd_cutoff_tumor"]] <- sd_cutoff_tumor
    
    return(preprocessing_output)
    
}

