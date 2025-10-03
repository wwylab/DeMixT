#' @title gene_selection_DE
#' @description Perform gene selection by first normalizing the log-transformed data and 
#' doing the differential expressions analysis.
#' @name gene_selection_DE
#' @param count.matrix A matrix of raw expression count with \eqn{G} by \eqn{(My + M1)}, where \eqn{G} is the number
#' of genes, \eqn{My} is the number of mixed samples and \eqn{M1} is the number of normal samples. Row names are genes
#' column names are sample ids.
#' @param normal.id A vector of normal sample ids
#' @param tumor.id A vector of tumor sample ids
#' @param ngene.selected.for.normalization A integer number indicating the number of genes selected before running DeMixNB
#' @param cutoff_normal_range A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the normal samples to subset. Default is c(0.1, 0.6)
#' @param cutoff_tumor_range A vector of two numeric values, indicating the lower and upper bounds of standard deviation of 
#' log2 count matrix from the tumor samples to subset. Default is c(0.2, 0.8)
#' @param cutoff_step A scatter value indicating the step size of changing cutoff_normal_range and cutoff_tumor_range to find a 
#' suitable subset of count matrix for downstream analysis
#' @param filter.sd The cut-off for the standard deviation of lognormal 
#' distribution. Genes whose log transferred standard deviation smaller than
#' the cut-off will be selected into the model. The default is TRUE.
#' @param ngene.selected.for.pi The percentage or the number of genes used for
#' proportion estimation. The difference between the expression levels from
#' mixed tumor samples and the known component(s) are evaluated, and the most
#' differential expressed genes are selected, which is called DE. It is enabled
#' when if.filter = TRUE. The default is \eqn{min(1500, 0.3*G)}, where
#' \eqn{G} is the number of genes. Users can also try using more genes,
#' ranging from \eqn{0.3*G} to \eqn{0.5*G}, and evaluate the outcome.
#' @return The selected genes
#' @export gene_selection_DE
gene_selection_DE <- function(count.matrix, normal.id, tumor.id, 
                                  ngene.selected.for.normalization = 9000,
                                  cutoff_normal_range=c(0.1, 1.0), 
                                  cutoff_tumor_range=c(0, 2.5), 
                                  cutoff_step=0.2,
                                  filter.sd=NA,
                                  ngene.selected.for.pi=NA) {
  preprocessed_data = DeMixNB_preprocessing(count.matrix, 
                                            normal.id, 
                                            tumor.id, 
                                            ngene.selected.for.normalization,
                                            cutoff_normal_range, 
                                            cutoff_tumor_range, 
                                            cutoff_step)
  count_filter = preprocessed_data$count.matrix
  sd_cutoff_normal = preprocessed_data$sd_cutoff_normal
  sd_cutoff_tumor = preprocessed_data$sd_cutoff_tumor
  data.N1_filter = count_filter[,normal.id]
  if (is.na(filter.sd))
  {
    filter.sd = sd_cutoff_normal[2]
  }
  
  ## filter out genes with constant value across all samples
  if (dim(count_filter)[1] == 1) {
    count_filter <- t(as.matrix(count_filter[apply(data.N1_filter,
                                             1, function(x) length(unique(x)) > 1), ]))
    data.N1_filter <- t(as.matrix(data.N1_filter[apply(data.N1_filter,
                                         1, function(x) length(unique(x)) > 1), ]))
  }
  else {
    count_filter <- count_filter[apply(data.N1_filter, 1, function(x) 
      length(unique(x)) > 1), ]
    data.N1_filter <- data.N1_filter[apply(data.N1_filter, 1, function(x) 
      length(unique(x)) > 1), ]
  }
  
  # Helper function for filtering
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
  
  inputdatans <- apply(log2(data.N1_filter), 1, sd)
  id1 <- (inputdatans < filter.sd)
  if (sum(id1) < 20) 
    stop("The threshold of standard variation is too stringent. \n
           Please provide a larger threshold. ")
  
  inputdatamat1 <- count_filter[id1, ]
  
  inputdatamat1nm <- rowMeans(inputdatamat1[, normal.id])
  inputdatamat1ym <- rowMeans(inputdatamat1[, tumor.id])
  inputdata1r <- inputdatamat1ym/inputdatamat1nm
  inputdata2 <- filter2(inputdata1r, ngene.selected.for.pi)
  gene.name <- rownames(inputdata2)
  
  return(gene.name)
}