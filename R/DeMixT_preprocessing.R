#' @title detect_suspicious_sample_by_hierarchical_clustering_2comp
#' @description Detect suspicious samples by a hierarchical clustering
#' 
#' This function is designed to evaluate the separation of tumor samples and normal samples in a PCA space. 
#' If some normal samples appear in the tumor-sample dominated cluster, these normal samples are likely to 
#' be tumor samples and they are supposed to be filtered out before downstream analysis. But for those tumor 
#' samples appearing in the normal-sample dominated cluster, we do not remove them since they may be the ones 
#' with low tumor purity.
#' @name detect_suspicious_sample_by_hierarchical_clustering_2comp
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @return list object
#' 
#' @export detect_suspicious_sample_by_hierarchical_clustering_2comp
detect_suspicious_sample_by_hierarchical_clustering_2comp <- function(count.matrix, normal.id, tumor.id){
  if(length(normal.id) + length(tumor.id) != ncol(count.matrix)){
    stop("Total number of normal and tumor samples in normal.id and tumor.id must be the same with the numbe of columns in count.matrix.")
  }
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
  
  labels <- y[top2.pcs.hclust$order]
  labels[labels == "Tumor"] <- "+"
  labels[labels == "Normal"] <- "."
  labels_colors <- labels
  labels_colors[labels_colors == '+'] <- "blue"
  labels_colors[labels_colors == '.'] <- "red"
  
  
  par(mar = c(1, 2.5, 2, 0))
  as.dendrogram(top2.pcs.hclust) %>% set("labels", labels) %>% set("labels_col", labels_colors) %>% plot(main = "Hierarchical clustering of tumor and normal samples")
  as.dendrogram(top2.pcs.hclust) %>% rect.dendrogram(k=2, border = 8, lty = 2, lwd = 1)
  
  legend(
    "topright",
    legend = c("Normal", "Tumor"),
    pch = c(20, 3),
    col = c("red", "blue"),
    bty = "n"
  )
  
  return(hc_labels)
}

#' @title plot_sd
#' @description Plot the standard deviation of log2 raw expression
#' @name plot_sd
#' @rdname detect_suspicious_sample_by_hierarchical_clustering_2comp
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @return 
#' 
#' @export plot_sd
plot_sd  <- function(count.matrix, normal.id, tumor.id){
  if(length(normal.id) + length(tumor.id) != ncol(count.matrix)){
    stop("Total number of normal and tumor samples in normal.id and tumor.id must be the same with the numbe of columns in count.matrix.")
  }
  count.matrix[which(count.matrix == 0, arr.ind = T)] = 1
  sdn.obs <- apply(log2(count.matrix[, match(normal.id, colnames(count.matrix))]), 1, sd)
  sdm.obs <- apply(log2(count.matrix[, match(tumor.id, colnames(count.matrix))]), 1, sd)
  par(mfrow = c(1, 2))
  plot(sdn.obs, 
       xlab = 'Genes', ylab = 'Standard Deviation', main = 'Normal Reference')
  plot(sdm.obs, 
       xlab = 'Genes', ylab = 'Standard Deviation', main = 'Mixed Tumor')
  par(mfrow = c(1,1))
}

#' @title subset_sd
#' @description Subset a count matrix given the the ranges of the standard deviations of the 
#' log2 expressions from the tumor and normal samples 
#' @name subset_sd
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param cutoff_normal A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.1, 0.6)
#' @param cutoff_tumor A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the tumor samples to subset. Default is c(0.2, 0.8)
#' @return A subset of the count matrix
#' 
#' @export subset_sd
subset_sd <- function(count.matrix, normal.id, tumor.id, 
                      cutoff_normal = c(0.1, 0.6), 
                      cutoff_tumor = c(0.2, 0.8)){
  
  count.matrix[which(count.matrix == 0, arr.ind = T)] = 1
  sdn.obs <- apply(log2(count.matrix[, match(normal.id, colnames(count.matrix))]), 1, sd)
  sdm.obs <- apply(log2(count.matrix[, match(tumor.id, colnames(count.matrix))]), 1, sd)
  indx <- which(sdn.obs > cutoff_normal[1] & sdn.obs < cutoff_normal[2] & sdm.obs > cutoff_tumor[1] & sdm.obs < cutoff_tumor[2])
  count.matrix.subset <- count.matrix[indx, ]
  return(count.matrix.subset)
}

#' @title plot_dim
#' @description Plot the distribution of tumor and normal samples in a 2D PCA space based on their expressions
#' @name plot_dim
#' @rdname detect_suspicious_sample_by_hierarchical_clustering_2comp
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param legend.position Position of legend in the plot. Default is bottomleft.
#' @param legend.cex Character expansion factor relative to current par("cex"). Default = 1.2
#' @return
#' 
#' @export plot_dim
plot_dim <- function(count.matrix, labels, 
                     legend.position = 'bottomleft',
                     legend.cex = 1.2){
  
  if(length(labels) != ncol(count.matrix)){
    stop("The length of labels must be the same with the numbe of columns in count.matrix.")
  }
  
  N.label <- length(unique(labels))
  
  count.matrix[which(count.matrix == 0, arr.ind = T)] = 1
  ## PCA dimension reduction
  res = 0
  res <- prcomp(t(log2(count.matrix)), center = T, scale. = T)
  ## PC variance
  pct1 = round(res$sdev[1]/sum(res$sdev), 3)*100
  pct2 = round(res$sdev[2]/sum(res$sdev), 3)*100
  ## Dimension
  dim1 = res$x[,1]; dim2 = res$x[,2]
  ## x/ylab text
  xlab.text = paste("First Comp: ", as.character(pct1), "% variance", sep="")
  ylab.text = paste("Second Comp: ", as.character(pct2), "% variance", sep="")   
  main.text = 'PCA'
  
  ## Plot
  plot(dim1, dim2,  cex = 1, 
       xlab=xlab.text, ylab=ylab.text, main = main.text,
       col=rainbow(N.label)[as.numeric(factor(labels))],
       pch = as.numeric(factor((labels))), lwd=1.5)
  abline(h=0, v=0, col="brown", lty=2)
  abline(h=0, v=0, col="brown", lty=2)
  center1<-tapply(dim1, factor(labels), mean)
  center2<-tapply(dim2, factor(labels), mean)
  ## add lines and centers
  for (ii in 1:length(center1)) {
    groupi<-cbind(dim1,dim2)[as.numeric(factor(labels))==ii, 1:2]
    if (class(groupi)[1] =="matrix") {
      for (j in (1:nrow(groupi))) {
        segments(groupi[j,1], groupi[j,2], center1[ii], center2[ii], col = rainbow(N.label)[ii], lwd = 0.5)
      }
    }else {
      segments(groupi[1], groupi[2], center1[ii], center2[ii], col = rainbow(N.label)[ii], lwd = 0.5)
    }
  }
  points(center1, center2, pch = 7, lwd = 1.5, col = rainbow(N.label))
  legend(legend.position, legend=names(table(factor(labels))), 
         text.col=rainbow(N.label), pch=1:N.label, 
         col=rainbow(N.label), lty=1, cex = legend.cex)
}


#' @title scale_normalization_75th_percentile
#' @description Quantile normalization for the raw count matrix of tumor and normal reference using the 0.75 quantile scale normalization
#' @name scale_normalization_75th_percentile
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @return the scale normalized count matrix
#' 
#' @export scale_normalization_75th_percentile
scale_normalization_75th_percentile <- function(count.matrix){
  newt <- count.matrix
  colnames(newt) = NULL
  rownames(newt) = NULL
  
  designs=c(rep("0", dim(count.matrix)[2]))
  seqData=newSeqCountSet(as.matrix(newt), designs)
  
  # quantile normalization/total/median   ###try different normalization method###
  seqData=estNormFactors(seqData, "quantile")
  k3=seqData@normalizationFactor
  mk3=median(k3)
  k3=k3/mk3
  
  temp<-newt
  
  for(i in 1:ncol(newt)){
    temp[,i] = temp[,i]/k3[i]
  }
  count.matrix.normalized<-temp
  colnames(count.matrix.normalized)<-colnames(count.matrix)
  rownames(count.matrix.normalized)<-rownames(count.matrix)
  
  return(count.matrix.normalized)
}


#' @title subset_sd_gene_remaining
#' @description Find the cutoffs to filter out genes with large standard deviations of log2 expressions in both normal and tumor samples
#' @name subset_sd_gene_remaining
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param cutoff_normal_range A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_tumor_range A vector of two numeric values, indicating the lower and upper bounds to search standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_step A scatter value indicating the step size of changing cutoff_normal_range and cutoff_tumor_range to find a 
#' suitable subset of count matrix for downstream analysis
#' 
#' @export subset_sd_gene_remaining
subset_sd_gene_remaining <- function(count.matrix, normal.id, tumor.id, 
                                     cutoff_normal_range = c(0.2, 0.6), 
                                     cutoff_tumor_range = c(0.2, 0.8),
                                     cutoff_step = 0.2){
  
  if(length(normal.id) + length(tumor.id) != ncol(count.matrix)){
    stop("Total number of normal and tumor samples in normal.id and tumor.id must be the same with the numbe of columns in count.matrix.")
  }
  count.matrix[which(count.matrix == 0, arr.ind = T)] = 1
  sdn.obs <- apply(log2(count.matrix[, match(normal.id, colnames(count.matrix))]), 1, sd)
  sdm.obs <- apply(log2(count.matrix[, match(tumor.id, colnames(count.matrix))]), 1, sd)
  
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

#' @title DeMixT_preprocessing
#' @description DeMixT preprocessing in one go
#' @name DeMixT_preprocessing
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param selected.genes A integer number indicating the number of genes selected before running DeMixT with the GS (Gene Selection) method
#' @param cutoff_normal_range A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_tumor_range A vector of two numeric values, indicating the lower and upper bounds to search standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_step A scatter value indicating the step size of changing cutoff_normal_range and cutoff_tumor_range to find a 
#' suitable subset of count matrix for downstream analysis
#' 
#' @return processed count matrix 
#' 
#' @export DeMixT_preprocessing
DeMixT_preprocessing <- function(count.matrix, normal.id, tumor.id, 
                                 selected.genes = 9000,
                                 cutoff_normal_range=c(0.1, 1.0), 
                                 cutoff_tumor_range=c(0, 2.5), 
                                 cutoff_step=0.2){
  

  stopifnot(cutoff_normal_range[1] >= 0)
  stopifnot(cutoff_tumor_range[1] >= 0)
  stopifnot(cutoff_normal_range[2] >= cutoff_normal_range[0])
  stopifnot(cutoff_tumor_range[2] >= cutoff_tumor_range[0])

  stopifnot(cutoff_step > 0)
  stopifnot(selected.genes > 0)
  ##stopifnot(is.integer(selected.genes))

  num_gene_remaining_different_cutoffs <- subset_sd_gene_remaining(count.matrix, normal.id, tumor.id,
                                                                   cutoff_normal_range, 
                                                                   cutoff_tumor_range,
                                                                   cutoff_step)
  
  
  num_gene_remaining_different_cutoffs_filter <- num_gene_remaining_different_cutoffs[which(num_gene_remaining_different_cutoffs$num.gene.remaining - selected.genes > 0), ]
  num_gene_remaining_different_cutoffs_filter <- num_gene_remaining_different_cutoffs_filter[order(num_gene_remaining_different_cutoffs_filter$num.gene.remaining), ]
  sd_cutoff_normal <- c(num_gene_remaining_different_cutoffs_filter$normal.cutoff.low[1], num_gene_remaining_different_cutoffs_filter$normal.cutoff.high[1])
  sd_cutoff_tumor <- c(num_gene_remaining_different_cutoffs_filter$tumor.cutoff.low[1], num_gene_remaining_different_cutoffs_filter$tumor.cutoff.high[1])
  
  count.matrix <- subset_sd(count.matrix, normal.id, tumor.id, cutoff_normal = sd_cutoff_normal, cutoff_tumor = sd_cutoff_tumor)
  
  count.matrix = scale_normalization_75th_percentile(count.matrix)
  
  preprocessing_output <- list()
  
  preprocessing_output[['count.matrix']] <- count.matrix
  
  preprocessing_output[["sd_cutoff_normal"]] <- sd_cutoff_normal
  preprocessing_output[["sd_cutoff_tumor"]] <- sd_cutoff_tumor
  
  return(preprocessing_output)
  
}

#' @title batch_correction
#' @description Batch effect correction for multiple batches of tumor samples using ComBat
#' @name batch_correction
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed tumor samples. Row names are genes
#'  column names are tumor sample ids. 
#' @param batch_labels Factor of tumor samples from different batches
#' @return Batch effect corrected count matrix for tumor samples
#' 
#' @export batch_correction
batch_correction <- function(count.matrix, batch_labels){
  
  if(length(batch_labels) != ncol(count.matrix)){
    stop("Total number of normal and tumor samples in normal.id and tumor.id must be the same with the numbe of columns in count.matrix.")
  }
  count.matrix[which(count.matrix == 0, arr.ind = T)] = 1
  count.matrix.combat.log2 = ComBat(dat = log2(count.matrix),
                                             batch = batch_labels, mod = NULL,
                                             par.prior=TRUE, mean.only = F)
  
  count.matrix.combat.combat <- 2^count.matrix.combat.log2
  return(count.matrix.combat.combat)
}


