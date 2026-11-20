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
  if (is.na(ngene.selected.for.pi)) {
    ngene.selected.for.pi <- min(1500, round(0.3*nrow(count.matrix)))
  }
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

#' @title gene_selection_MoransI
#' @description Select spatially variable genes for DeMixNB from spatial
#' transcriptomics data. The counts of the reference and mixed spots are
#' normalized with SCTransform, Moran's I is computed for the top variable
#' genes, and the genes with the largest observed Moran's I are returned.
#' This follows the gene selection in the DeMixNB tutorial.
#' @name gene_selection_MoransI
#' @param seurat.obj A Seurat object with a spatial assay of raw counts and a
#' tissue image holding the spot coordinates, e.g. from
#' \code{Seurat::Load10X_Spatial}. Quality control of genes and spots should be
#' done beforehand.
#' @param normal.id A vector of reference (non-tumor) spot ids
#' @param tumor.id A vector of mixed (tumor) spot ids
#' @param ngene.selected An integer giving the number of genes to return,
#' ranked by observed Moran's I. Default is 1000.
#' @param n.features.tested An integer giving the number of top SCTransform
#' variable genes for which Moran's I is computed. Moran's I only re-ranks
#' these genes, so with the default n.features.tested = ngene.selected = 1000,
#' as in the tutorial, the returned genes are the top 1000 variable genes
#' ordered by Moran's I. Set n.features.tested larger than ngene.selected to
#' select among more genes. Default is 1000.
#' @param assay The name of the assay with raw counts. Default is "Spatial".
#' @param image The name of the tissue image used for the spot coordinates.
#' Default is NULL, i.e. the first image of the object.
#' @param sct.args A list of further arguments passed to
#' \code{Seurat::SCTransform}. Default is an empty list.
#' @param return.stats A boolean flag. If TRUE, also return the Moran's I
#' statistics of all tested genes. Default is FALSE.
#' @param verbose A boolean flag controlling the messages from Seurat.
#' Default is FALSE.
#' @return The selected gene names, ordered by decreasing observed Moran's I.
#' If return.stats = TRUE, a list with gene.name and moransi, a data frame of
#' the observed Moran's I and its p-value for all tested genes.
#' @details Requires the Seurat package. Installing Rfast2 makes the Moran's I
#' computation much faster, and glmGamPoi speeds up SCTransform.
#' SCTransform is run on the reference and mixed spots together. The selected
#' genes are meant to subset the raw count matrices passed to DeMixNB.
#' @export gene_selection_MoransI
gene_selection_MoransI <- function(seurat.obj, normal.id, tumor.id,
                                   ngene.selected = 1000,
                                   n.features.tested = 1000,
                                   assay = "Spatial",
                                   image = NULL,
                                   sct.args = list(),
                                   return.stats = FALSE,
                                   verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("gene_selection_MoransI requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat"))
    stop("seurat.obj must be a Seurat object.")
  if (!assay %in% names(seurat.obj@assays))
    stop("The assay '", assay, "' is not in seurat.obj.")
  if (length(intersect(normal.id, tumor.id)) > 0)
    stop("normal.id and tumor.id must not overlap.")
  spots <- c(normal.id, tumor.id)
  missing.spots <- setdiff(spots, colnames(seurat.obj))
  if (length(missing.spots) > 0)
    stop(length(missing.spots), " spot ids in normal.id or tumor.id are not in seurat.obj.")
  for (arg in c("ngene.selected", "n.features.tested")) {
    v <- get(arg)
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < 1)
      stop(arg, " must be a single number of at least 1.")
  }
  if (!is.list(sct.args)) stop("sct.args must be a list.")

  ## Step 1: SCTransform on the reference and mixed spots
  obj <- seurat.obj[, colnames(seurat.obj) %in% spots]
  obj <- do.call(Seurat::SCTransform,
                 c(list(object = obj, assay = assay, verbose = verbose), sct.args))

  ## Step 2: Moran's I for the top variable genes
  var.features <- Seurat::VariableFeatures(obj)
  features <- var.features[seq_len(min(n.features.tested, length(var.features)))]
  obj <- Seurat::FindSpatiallyVariableFeatures(obj, assay = "SCT",
                                               features = features,
                                               image = image,
                                               selection.method = "moransi",
                                               verbose = verbose)

  ## Step 3: Rank the tested genes by observed Moran's I
  meta <- obj[["SCT"]][[]]
  if (!"MoransI_observed" %in% colnames(meta))
    stop("Moran's I results were not found in the SCT assay.")
  meta <- meta[!is.na(meta$MoransI_observed), , drop = FALSE]
  meta <- meta[order(meta$MoransI_observed, decreasing = TRUE), , drop = FALSE]
  if (nrow(meta) < ngene.selected)
    warning("Only ", nrow(meta), " genes have Moran's I, fewer than ngene.selected = ",
            ngene.selected, ".", call. = FALSE)
  gene.name <- rownames(meta)[seq_len(min(ngene.selected, nrow(meta)))]

  if (return.stats) {
    stats.cols <- intersect(c("MoransI_observed", "MoransI_p.value"), colnames(meta))
    return(list(gene.name = gene.name, moransi = meta[, stats.cols, drop = FALSE]))
  }
  return(gene.name)
}
