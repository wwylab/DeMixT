#' @title subset_sd
#' @description Subset a count matrix given the the ranges of the standard deviations of the 
#' log2 expressions from the tumor and normal samples. 
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
  indx <- which(sdn.obs > cutoff_normal[1] & sdn.obs < cutoff_normal[2] & 
                  sdm.obs > cutoff_tumor[1] & sdm.obs < cutoff_tumor[2])
  count.matrix.subset <- count.matrix[indx, ]
  return(count.matrix.subset)
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
  seqData=DSS::newSeqCountSet(as.matrix(newt), designs)
  
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

#' @title DeMixNB_preprocessing
#' @description DeMixNB preprocessing in one go
#' @name DeMixNB_preprocessing
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#'  column names are sample ids. 
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param selected.genes A integer number indicating the number of genes selected before running DeMixNB
#' @param cutoff_normal_range A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_tumor_range A vector of two numeric values, indicating the lower and upper bounds to search standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.2, 0.6)
#' @param cutoff_step A scatter value indicating the step size of changing cutoff_normal_range and cutoff_tumor_range to find a 
#' suitable subset of count matrix for downstream analysis
#' 
#' @return processed count matrix 
#' 
#' @export DeMixNB_preprocessing
DeMixNB_preprocessing <- function(count.matrix, normal.id, tumor.id, 
                                 selected.genes = 9000,
                                 cutoff_normal_range=c(0.1, 1.0), 
                                 cutoff_tumor_range=c(0, 2.5), 
                                 cutoff_step=0.2){
  

  stopifnot(cutoff_normal_range[1] >= 0)
  stopifnot(cutoff_tumor_range[1] >= 0)
  stopifnot(cutoff_normal_range[2] >= cutoff_normal_range[1])
  stopifnot(cutoff_tumor_range[2] >= cutoff_tumor_range[1])

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
