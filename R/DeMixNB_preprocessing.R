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


##----------------------------------------
## Preprocessing of spatial transcriptomics data (DeMixNB tutorial)
##----------------------------------------

#' @title filter_genes_ST
#' @description Remove unwanted genes (e.g., mitochondrial genes and pseudogenes) 
#' by name, and genes detected in fewer than min.spots spots.
#' @name filter_genes_ST
#' @param seurat.obj A Seurat object with spatial counts.
#' @param assay The name of the assay with raw counts. Default is "Spatial".
#' @param unwanted.pattern A regular expression for the names of genes to
#' remove. Set to NULL to skip the name filter. Default is the pattern used in
#' the DeMixNB tutorial.
#' @param min.spots The minimum number of spots in which a gene must be
#' detected (count > 0). Default is 10.
#' @return The Seurat object with the remaining genes.
#' @export filter_genes_ST
filter_genes_ST <- function(seurat.obj, assay = "Spatial",
                            unwanted.pattern = "^(MT-|MIR-|LINC|A[A-Z]|B[A-Z]|\n                        C[A-Z]|Rik|Gm|Vmn|Olfr|Trav)",
                            min.spots = 10) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("This function requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat")) stop("seurat.obj must be a Seurat object.")
  if (!is.numeric(min.spots) || length(min.spots) != 1 || min.spots < 0)
    stop("min.spots must be a single non-negative number.")
  counts <- Seurat::GetAssayData(seurat.obj, assay = assay, layer = "counts")
  genes <- rownames(counts)
  unwanted <- if (is.null(unwanted.pattern)) rep(FALSE, length(genes)) else grepl(unwanted.pattern, genes)
  n.detected <- Matrix::rowSums(counts > 0)
  keep <- genes[!unwanted & n.detected >= min.spots]
  message(sum(unwanted), " genes removed by name, ",
          sum(!unwanted & n.detected < min.spots), " genes detected in fewer than ",
          min.spots, " spots removed, ", length(keep), " genes kept.")
  return(seurat.obj[keep, ])
}

#' @title summarize_spot_annotation
#' @description Summarize a cell-level annotation table (one row per cell, 
#' with the spot barcode and the cell class) into spot-level cell counts and 
#' proportions of tumor, stromal and immune cells.
#' @name summarize_spot_annotation
#' @param cell.annotation A data frame with one row per cell.
#' @param barcodes A vector of spot barcodes to keep, e.g.
#' colnames(seurat.obj). Default is NULL, i.e. all barcodes.
#' @param barcode.col The column with the spot barcode. Default is "Barcode".
#' @param class.col The column with the cell class. Default is "class2".
#' @param tumor.label The class label of tumor cells. Default is "t".
#' @param stromal.label The class label of stromal cells. Default is "f".
#' @param immune.label The class label of immune cells. Default is "l".
#' @param x.col The column with the x coordinate of the cells, or NULL.
#' Default is "x_fullres".
#' @param y.col The column with the y coordinate of the cells, or NULL.
#' Default is "y_fullres".
#' @return A data frame with one row per spot, sorted by barcode, with columns
#' Barcode, x_coor and y_coor (mean cell coordinates, if x.col and y.col are
#' given), tumor_proportion, stromal_proportion, immune_proportion,
#' non_tumor_proportion, tumor_cell and number_of_cells.
#' @details Cells with a missing or empty barcode are removed. Cells with a
#' missing class count as non-tumor cells.
#' @export summarize_spot_annotation
summarize_spot_annotation <- function(cell.annotation, barcodes = NULL,
                                      barcode.col = "Barcode", class.col = "class2",
                                      tumor.label = "t", stromal.label = "f",
                                      immune.label = "l",
                                      x.col = "x_fullres", y.col = "y_fullres") {
  cell.annotation <- as.data.frame(cell.annotation)
  needed <- c(barcode.col, class.col, x.col, y.col)
  missing.cols <- setdiff(needed, colnames(cell.annotation))
  if (length(missing.cols) > 0)
    stop("Columns not found in cell.annotation: ", paste(missing.cols, collapse = ", "))
  bc <- as.character(cell.annotation[[barcode.col]])
  keep <- !is.na(bc) & bc != ""
  if (!is.null(barcodes)) keep <- keep & bc %in% barcodes
  if (!any(keep)) stop("No cells left after matching the barcodes.")
  cell.annotation <- cell.annotation[keep, , drop = FALSE]
  bc <- bc[keep]
  cls <- as.character(cell.annotation[[class.col]])
  
  n.cells <- tapply(rep(1, length(bc)), bc, sum)
  spots <- names(n.cells)
  count_class <- function(label) as.numeric(tapply(cls %in% label, bc, sum)[spots])
  tumor.cells <- count_class(tumor.label)
  spot.summary <- data.frame(Barcode = spots, stringsAsFactors = FALSE)
  if (!is.null(x.col) && !is.null(y.col)) {
    spot.summary$x_coor <- as.numeric(tapply(cell.annotation[[x.col]], bc, mean)[spots])
    spot.summary$y_coor <- as.numeric(tapply(cell.annotation[[y.col]], bc, mean)[spots])
  }
  n.cells <- as.numeric(n.cells)
  spot.summary$tumor_proportion <- tumor.cells / n.cells
  spot.summary$stromal_proportion <- count_class(stromal.label) / n.cells
  spot.summary$immune_proportion <- count_class(immune.label) / n.cells
  spot.summary$non_tumor_proportion <- 1 - spot.summary$tumor_proportion
  spot.summary$tumor_cell <- tumor.cells
  spot.summary$number_of_cells <- n.cells
  return(spot.summary)
}

#' @title qc_spots_ST
#' @description Keep the spots with enough total counts and detected genes.
#' @name qc_spots_ST
#' @param seurat.obj A Seurat object with spatial counts.
#' @param assay The name of the assay with raw counts. Default is "Spatial".
#' @param min.counts The minimum total count of a spot. Default is 600.
#' @param min.features The minimum number of detected genes of a spot. Default
#' is 300.
#' @return The Seurat object with the spots that pass both thresholds.
#' @details The thresholds are applied to the nCount and nFeature columns of
#' the spot metadata, e.g. nCount_Spatial and nFeature_Spatial.
#' @export qc_spots_ST
qc_spots_ST <- function(seurat.obj, assay = "Spatial",
                        min.counts = 600, min.features = 300) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("This function requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat")) stop("seurat.obj must be a Seurat object.")
  meta <- seurat.obj@meta.data
  count.col <- paste0("nCount_", assay)
  feature.col <- paste0("nFeature_", assay)
  if (!all(c(count.col, feature.col) %in% colnames(meta)))
    stop("The spot metadata must contain ", count.col, " and ", feature.col, ".")
  keep <- rownames(meta)[meta[[count.col]] >= min.counts & meta[[feature.col]] >= min.features]
  keep <- keep[!is.na(keep)]
  message(length(keep), " of ", nrow(meta), " spots pass the QC thresholds.")
  if (length(keep) == 0) stop("No spots pass the QC thresholds.")
  return(seurat.obj[, keep])
}

#' @title filter_spots_by_annotation
#' @description Keep the spots with enough annotated cells whose tumor 
#' proportion is either 0 (reference spots) or within tumor.prop.range (mixed spots).
#' @name filter_spots_by_annotation
#' @param spot.summary A data frame from \code{summarize_spot_annotation}.
#' @param min.cells The minimum number of annotated cells in a spot. Default
#' is 5.
#' @param tumor.prop.range The lower and upper bounds of the tumor proportion
#' of mixed spots. Default is c(0.05, 0.95).
#' @return The rows of spot.summary that pass the filter.
#' @export filter_spots_by_annotation
filter_spots_by_annotation <- function(spot.summary, min.cells = 5,
                                       tumor.prop.range = c(0.05, 0.95)) {
  if (!all(c("tumor_proportion", "number_of_cells") %in% colnames(spot.summary)))
    stop("spot.summary must contain tumor_proportion and number_of_cells.")
  if (!is.numeric(tumor.prop.range) || length(tumor.prop.range) != 2 ||
      tumor.prop.range[1] > tumor.prop.range[2])
    stop("tumor.prop.range must be c(lower, upper) with lower <= upper.")
  tp <- spot.summary$tumor_proportion
  keep <- spot.summary$number_of_cells >= min.cells &
    (tp == 0 | (tp >= tumor.prop.range[1] & tp <= tumor.prop.range[2]))
  keep[is.na(keep)] <- FALSE
  return(spot.summary[keep, , drop = FALSE])
}

#' @title classify_spots
#' @description Split the spots into reference spots (tumor proportion 0) and 
#' mixed spots (tumor proportion above 0).
#' @name classify_spots
#' @param spot.summary A data frame from \code{filter_spots_by_annotation} or
#' \code{summarize_spot_annotation}.
#' @return A list with normal.id, the barcodes of the reference spots, and
#' tumor.id, the barcodes of the mixed spots.
#' @export classify_spots
classify_spots <- function(spot.summary) {
  if (!all(c("Barcode", "tumor_proportion") %in% colnames(spot.summary)))
    stop("spot.summary must contain Barcode and tumor_proportion.")
  tp <- spot.summary$tumor_proportion
  normal.id <- spot.summary$Barcode[!is.na(tp) & tp == 0]
  tumor.id <- spot.summary$Barcode[!is.na(tp) & tp > 0]
  message(length(normal.id), " reference spots and ", length(tumor.id), " mixed spots.")
  return(list(normal.id = normal.id, tumor.id = tumor.id))
}

#' @title subset_spots_ST
#' @description Subset the Seurat object to the reference and mixed spots and 
#' label them in the spot metadata.
#' @name subset_spots_ST
#' @param seurat.obj A Seurat object with spatial counts.
#' @param normal.id A vector of reference spot barcodes.
#' @param tumor.id A vector of mixed spot barcodes.
#' @return The Seurat object with the reference and mixed spots, with a
#' sample_type column in the spot metadata ("Normal" or "Tumor Mixed").
#' @export subset_spots_ST
subset_spots_ST <- function(seurat.obj, normal.id, tumor.id) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("This function requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat")) stop("seurat.obj must be a Seurat object.")
  if (length(intersect(normal.id, tumor.id)) > 0)
    stop("normal.id and tumor.id must not overlap.")
  missing.spots <- setdiff(c(normal.id, tumor.id), colnames(seurat.obj))
  if (length(missing.spots) > 0)
    stop(length(missing.spots), " spot ids in normal.id or tumor.id are not in seurat.obj.")
  obj <- seurat.obj[, colnames(seurat.obj) %in% c(normal.id, tumor.id)]
  sample.type <- ifelse(colnames(obj) %in% normal.id, "Normal", "Tumor Mixed")
  names(sample.type) <- colnames(obj)
  obj <- Seurat::AddMetaData(obj, metadata = sample.type, col.name = "sample_type")
  return(obj)
}

#' @title plot_spot_classification
#' @description Plot the reference and mixed spots on the tissue image.
#' @name plot_spot_classification
#' @param seurat.obj A Seurat object from \code{subset_spots_ST}.
#' @param group.by The metadata column to color the spots by. Default is
#' "sample_type".
#' @param cols The colors of the groups. Default is c(Normal = "blue", "Tumor
#' Mixed" = "red").
#' @param pt.size.factor The scale factor of the spot size. Default is 1.5.
#' @param ... Further arguments passed to \code{Seurat::SpatialDimPlot}.
#' @return A ggplot object.
#' @export plot_spot_classification
plot_spot_classification <- function(seurat.obj, group.by = "sample_type",
                                     cols = c(Normal = "blue", "Tumor Mixed" = "red"),
                                     pt.size.factor = 1.5, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("This function requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat")) stop("seurat.obj must be a Seurat object.")
  Seurat::SpatialDimPlot(seurat.obj, group.by = group.by, cols = cols,
                         pt.size.factor = pt.size.factor, interactive = FALSE,
                         crop = FALSE, ...)
}

#' @title prepare_DeMixNB_input
#' @description Extract the raw counts of the mixed and reference spots as the 
#' inputs of \code{DeMixNB}.
#' @name prepare_DeMixNB_input
#' @param seurat.obj A Seurat object with spatial counts.
#' @param normal.id A vector of reference spot barcodes.
#' @param tumor.id A vector of mixed spot barcodes.
#' @param genes A vector of genes to keep, e.g. from
#' \code{gene_selection_MoransI}. Default is NULL, i.e. all genes.
#' @param assay The name of the assay with raw counts. Default is "Spatial".
#' @return A list with data.Y, the counts of the mixed spots, and data.N1, the
#' counts of the reference spots, each a SummarizedExperiment object with
#' genes in rows.
#' @examples
#' \dontrun{
#' seurat.obj <- filter_genes_ST(seurat.obj)
#' spot.summary <- summarize_spot_annotation(cell.annotation,
#'                                           barcodes = colnames(seurat.obj))
#' seurat.obj <- qc_spots_ST(seurat.obj)
#' spot.summary <- summarize_spot_annotation(cell.annotation,
#'                                           barcodes = colnames(seurat.obj))
#' spot.summary <- filter_spots_by_annotation(spot.summary)
#' ids <- classify_spots(spot.summary)
#' combined <- subset_spots_ST(seurat.obj, ids$normal.id, ids$tumor.id)
#' plot_spot_classification(combined)
#' input <- prepare_DeMixNB_input(combined, ids$normal.id, ids$tumor.id)
#' res <- DeMixNB(input$data.Y, input$data.N1, gene.selection.method = "ST",
#'                seurat.obj = combined)
#' }
#' @export prepare_DeMixNB_input
prepare_DeMixNB_input <- function(seurat.obj, normal.id, tumor.id,
                                  genes = NULL, assay = "Spatial") {
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("This function requires the Seurat package. Please install it first.")
  if (!inherits(seurat.obj, "Seurat")) stop("seurat.obj must be a Seurat object.")
  counts <- Seurat::GetAssayData(seurat.obj, assay = assay, layer = "counts")
  if (is.null(genes)) genes <- rownames(counts)
  missing.genes <- setdiff(genes, rownames(counts))
  if (length(missing.genes) > 0)
    stop(length(missing.genes), " genes are not in the counts of seurat.obj.")
  missing.spots <- setdiff(c(normal.id, tumor.id), colnames(counts))
  if (length(missing.spots) > 0)
    stop(length(missing.spots), " spot ids in normal.id or tumor.id are not in seurat.obj.")
  data.Y <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(counts[genes, tumor.id, drop = FALSE])))
  data.N1 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(counts[genes, normal.id, drop = FALSE])))
  return(list(data.Y = data.Y, data.N1 = data.N1))
}
