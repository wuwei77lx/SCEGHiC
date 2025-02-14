#' Create aggregated data for SCEG-HiC
#'
#' Function to generate aggregated inputs for SCEG-HiC.
#' \code{Aggregate_datasets}
#'
#' @param object Seurat object.
#' @param k_neigh Number of cells to be aggregated per group.
#' @param rna_assay Name of assay containing gene expression information.
#' @param atac_assay Name of assay containing peak information.
#' @param cellnames Name of one or more metadata columns to group (color) the
#' cells by Default is the current cell identities.
#' @param atacbinary Logical,whether the aggregated scATAC-seq data need binary?
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param reduction.name The reduction name of extracting the cell coordinates
#' used for aggregating.
#' @param size_factor_normalize Logical, should accessibility values be
#' normalized by size factor?
#' @param seed Random seed.
#' @param verbose Logical, should warning and info messages be printed?
#'
#' @import Seurat
#' @return Aggregated Seurat object.
#' @export
#'
Aggregate_datasets <- function(object, k_neigh = 50, rna_assay = "RNA",
                               atac_assay = "peaks", cellnames = NULL,
                               atacbinary = TRUE, max_overlap = 0.8,
                               reduction.name = NULL,
                               size_factor_normalize = TRUE,
                               seed = 123, verbose = TRUE) {
  if (!is.null(reduction.name)) {
    cell_coord <- object@reductions[[reduction.name]]
  } else {
    cell_coord <- object@reductions$umap@cell.embeddings
  }

  if (!(is.null(cellnames))) {
    group <- as.character(object@meta.data[[cellnames]])
  } else {
    group <- as.character(Idents(object))
  }
  uniqgroup <- unique(group)
  if (rna_assay %in% names(object@assays)) {
    if (inherits(object@assays[[rna_assay]], "Assay5")) {
      rna_new <- matrix(0, nrow = nrow(object@assays[[rna_assay]]@layers[["counts"]]), ncol = 1)
    } else {
      rna_new <- matrix(0, nrow = nrow(object@assays[[rna_assay]]@counts), ncol = 1)
    }
  }
  atac_new <- matrix(0, nrow = nrow(object@assays[[atac_assay]]@counts), ncol = 1)
  cell_sample <- matrix(0, nrow = 1, ncol = k_neigh)
  for (i in 1:length(uniqgroup)) {
    if (verbose) {
      message(paste0("Aggregating cluster ", uniqgroup[i]))
    }
    if (!(is.null(cellnames))) {
      keep <- which(object@meta.data[[cellnames]] %in% uniqgroup[i])
      subobject <- subset(object, cells = keep)
    } else {
      subobject <- subset(object, idents = uniqgroup[i]) # Seurat extracts the required data
    }
    sub_index <- which(group %in% uniqgroup[i])
    cell_coord_i <- cell_coord[sub_index, ]
    sub_aggregated_datasets <- generate_aggregated_datastes(subobject, cell_coord_i, rna_assay, atac_assay, k_neigh, atacbinary, max_overlap, seed, verbose)

    sub_cell_sample <- sub_aggregated_datasets$cell_sample
    message(paste0("have ", dim(sub_cell_sample)[1], " samples"))
    if (rna_assay %in% names(object@assays)) {
      rna_new <- cbind(rna_new, sub_aggregated_datasets$rna)
    }
    atac_new <- cbind(atac_new, sub_aggregated_datasets$atac)
    if (ncol(sub_cell_sample) < k_neigh) {
      sub_cell_sample_new <- as.matrix(sub_cell_sample)
      sub_cell_sample_new <- cbind(sub_cell_sample_new, matrix(0, nrow = 1, ncol = k_neigh - ncol(sub_cell_sample_new)))
    } else {
      sub_cell_sample_new <- apply(sub_cell_sample, 2, function(x) {
        sub_index[x] # for each column return original index
      })
      sub_cell_sample_new <- as.data.frame(sub_cell_sample_new)
      sub_cell_sample_new <- as.matrix(sub_cell_sample_new)
    }
    cell_sample <- rbind(cell_sample, sub_cell_sample_new)
  }
  if (rna_assay %in% names(object@assays)) {
    rna_new <- rna_new[, -1]
  }
  atac_new <- atac_new[, -1]
  cell_sample <- cell_sample[-1, ]
  # normalization
  if (size_factor_normalize) {
    if (rna_assay %in% names(object@assays)) {
      rna_new <- t(t(log(rna_new + 1)) / estimateSizeFactorsForMatrix(rna_new))
    }
    atac_new <- t(t(log(atac_new + 1)) / estimateSizeFactorsForMatrix(atac_new))
  }
  new_data <- list()
  if (rna_assay %in% names(object@assays)) {
    new_data$rna <- rna_new
  }
  new_data$atac <- atac_new
  new_data$cell_sample <- cell_sample
  return(new_data)
}

#' Create aggregated data for a certain cluster
#'
#' Function to generate aggregated inputs of a cetrain cluster. \code{generate_aggregated_datastes}
#' takes as input sparse data. This function will aggregate binary accessibility scores (or gene expression)
#' per cell cluster, if they do not overlap any existing group with more than 50% cells.
#'
#' @param object object Seurat object.
#' @param cell_coord cell_coord similarity matrix or dimiension reductions.
#' @param rna_assay Name of assay containing gene expression information.
#' @param atac_assay Name of assay containing peak information.
#' @param k_neigh Number of cells to aggregate per group.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary.
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param seed Random seed.
#' @param verbose verbose Logical, should warning and info messages be printed?
#'
#' @importFrom FNN knn.index
#' @import Matrix
#'
#' @return Aggregated data.
#' @export
#'
generate_aggregated_datastes <- function(object, cell_coord, rna_assay = "RNA", atac_assay = "peaks", k_neigh = 50, atacbinary = TRUE, max_overlap = 0.8, seed = 123, verbose = TRUE) {
  if (nrow(cell_coord) > k_neigh) {
    # Create a k-nearest neighbors map
    nn_map <- as.data.frame(FNN::knn.index(cell_coord,
      k = (k_neigh - 1)
    ))
    row.names(nn_map) <- row.names(cell_coord)
    nn_map$agg_cell <- 1:nrow(nn_map)
    good_choices <- 1:nrow(nn_map)

    if (verbose) {
      message("Sample cells randomly.")
    }

    # Sample cells randomly
    set.seed(seed)
    choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
    chosen <- good_choices[choice]
    good_choices <- good_choices[good_choices != good_choices[choice]]

    it <- 0
    ## Slow (contain calculating of overlapping between cell groups)
    while (length(good_choices) > 0 & it < nrow(cell_coord) / ((1 - max_overlap) * k_neigh)) {
      it <- it + 1
      choice <- sample(1:length(good_choices), size = 1, replace = FALSE)
      new_chosen <- c(chosen, good_choices[choice])
      good_choices <- good_choices[good_choices != good_choices[choice]]
      cell_sample <- nn_map[new_chosen, ]

      # calculate overlapping between cell groups
      combs <- data.frame(1:(nrow(cell_sample) - 1), nrow(cell_sample))
      shared <- apply(combs, 1, function(x) { # Slow
        (k_neigh * 2) - length(unique(as.vector(as.matrix(cell_sample[x, ]))))
      })
      if (max(shared) < max_overlap * k_neigh) {
        chosen <- new_chosen
      }
    }

    # aggregating both scRNA-seq and scATAC-seq counts of cells within one group
    if (rna_assay %in% names(object@assays)) {
      if (inherits(object@assays[[rna_assay]], "Assay5")) {
        rna_old <- as.matrix(object@assays[[rna_assay]]@layers[["counts"]])
        rownames(rna_old) <- rownames(object@assays[[rna_assay]])
      } else {
        rna_old <- as.matrix(object@assays[[rna_assay]]@counts)
      }
      rna_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(rna_old)) %in% cell_sample[x, , drop = FALSE])
      rna_mask <- Matrix::Matrix(rna_mask)
      rna_new <- rna_old %*% rna_mask
      rna_new <- as.matrix(rna_new)
    }

    atac_old <- object@assays[[atac_assay]]@counts
    # binarize
    if (atacbinary) {
      atac_old <- atac_old > 0
    }

    atac_mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(atac_old)) %in% cell_sample[x, , drop = FALSE])
    atac_mask <- Matrix::Matrix(atac_mask)
    atac_new <- atac_old %*% atac_mask
    atac_new <- as.matrix(atac_new)
  } else {
    if (rna_assay %in% names(object@assays)) {
      if (inherits(object@assays[[rna_assay]], "Assay5")) {
        rna_old <- as.matrix(object@assays[[rna_assay]]@layers[["counts"]])
        rownames(rna_old) <- rownames(object@assays[[rna_assay]])
      } else {
        rna_old <- as.matrix(object@assays[[rna_assay]]@counts)
      }
      rna_new <- rowSums(rna_old)
      rna_new <- as.matrix(rna_new)
    }

    atac_old <- object@assays[[atac_assay]]@counts
    # binarize
    if (atacbinary) {
      atac_old <- atac_old > 0
    }
    atac_new <- rowSums(atac_old)
    atac_new <- as.matrix(atac_new)
    cell_sample <- as.data.frame(t(matrix(seq(from = 1, to = nrow(cell_coord)))))
  }
  new_data <- list()
  if (rna_assay %in% names(object@assays)) {
    new_data$rna <- rna_new
  }
  new_data$atac <- atac_new
  new_data$cell_sample <- cell_sample
  return(new_data)
}

#' Function to calculate the size factor for the single-cell data
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts.
#' @param locfunc The location function used to find the representive value.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded.
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
#' @importFrom stats median
#' @import slam
#' @import Matrix
#' @export
#'
estimateSizeFactorsForMatrix <- function(counts, locfunc = median, round_exprs = TRUE, method = "mean-geometric-mean-total") {
  # library("slam")
  if (isSparseMatrix(counts)) {
    estimateSizeFactorsForSparseMatrix(counts, locfunc = locfunc, round_exprs = round_exprs, method = method)
  } else {
    estimateSizeFactorsForDenseMatrix(counts, locfunc = locfunc, round_exprs = round_exprs, method = method)
  }
}

#' Convert a sparseMatrix from Matrix package to a slam matrix
#'
#' @param sp_mat The matrix for the aggregated single cell data.
#'
#' @import slam
#' @export
#'
asSlamMatrix <- function(sp_mat) {
  sp <- Matrix::summary(sp_mat)
  simple_triplet_matrix(sp[, "i"], sp[, "j"], sp[, "x"], ncol = ncol(sp_mat), nrow = nrow(sp_mat), dimnames = dimnames(sp_mat))
}

#' Convert a sparseMatrix from Matrix package to a slam matrix
#'
#' @param x The sparseMatrix data.
#'
#' @import Matrix
#' @export
#'
isSparseMatrix <- function(x) {
  sum(class(x) %in% c("dgCMatrix", "dgTMatrix"))
}

#' Estimate size factors for each column, given a sparseMatrix from the Matrix package
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts.
#' @param locfunc The location function used to find the representive value.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded.
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
#' @import slam
#' @importFrom stats median
#' @export
#'
estimateSizeFactorsForSparseMatrix <- function(counts,
                                               locfunc = median,
                                               round_exprs = TRUE,
                                               method = "mean-geometric-mean-total") {
  CM <- counts
  if (round_exprs) {
    CM <- round(CM)
  }
  CM <- asSlamMatrix(CM)

  if (method == "weighted-median") {
    log_medians <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- rowapply_simple_triplet_matrix(CM, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <- weights * (log(cnts) - log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      # print (head(norm_cnts))
      exp(mean(norm_cnts))
    })
  } else if (method == "median-geometric-mean") {
    log_geo_means <- rowapply_simple_triplet_matrix(CM, function(x) {
      mean(log(CM))
    })

    sfs <- colapply_simple_triplet_matrix(CM, function(cnts) {
      norm_cnts <- log(cnts) - log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      # print (head(norm_cnts))
      exp(locfunc(norm_cnts))
    })
  } else if (method == "median") {
    stop("Error: method 'median' not yet supported for sparse matrices")
  } else if (method == "mode") {
    stop("Error: method 'mode' not yet supported for sparse matrices")
  } else if (method == "geometric-mean-total") {
    cell_total <- col_sums(CM)
    sfs <- log(cell_total) / mean(log(cell_total))
  } else if (method == "mean-geometric-mean-total") {
    cell_total <- col_sums(CM)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

#' Estimate size factors dense matrix
#'
#' @param counts The matrix for the gene expression data, either read counts or FPKM values or transcript counts.
#' @param locfunc The location function used to find the representive value.
#' @param round_exprs A logic flag to determine whether or not the expression value should be rounded.
#' @param method A character to specify the size factor calculation appraoches. It can be either "mean-geometric-mean-total" (default),
#' "weighted-median", "median-geometric-mean", "median", "mode", "geometric-mean-total".
#'
#' @importFrom stats median
#' @export
#'
estimateSizeFactorsForDenseMatrix <- function(counts, locfunc = median, round_exprs = TRUE, method = "mean-geometric-mean-total") {
  CM <- counts
  if (round_exprs) {
    CM <- round(CM)
  }
  if (method == "weighted-median") {
    log_medians <- apply(CM, 1, function(cell_expr) {
      log(locfunc(cell_expr))
    })

    weights <- apply(CM, 1, function(cell_expr) {
      num_pos <- sum(cell_expr > 0)
      num_pos / length(cell_expr)
    })

    sfs <- apply(CM, 2, function(cnts) {
      norm_cnts <- weights * (log(cnts) - log_medians)
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      # print (head(norm_cnts))
      exp(mean(norm_cnts))
    })
  } else if (method == "median-geometric-mean") {
    log_geo_means <- rowMeans(log(CM))

    sfs <- apply(CM, 2, function(cnts) {
      norm_cnts <- log(cnts) - log_geo_means
      norm_cnts <- norm_cnts[is.nan(norm_cnts) == FALSE]
      norm_cnts <- norm_cnts[is.finite(norm_cnts)]
      # print (head(norm_cnts))
      exp(locfunc(norm_cnts))
    })
  } else if (method == "median") {
    row_median <- apply(CM, 1, median)
    sfs <- apply(Matrix::t(Matrix::t(CM) - row_median), 2, median)
  } else if (method == "mode") {
    sfs <- estimate_t(CM)
  } else if (method == "geometric-mean-total") {
    cell_total <- apply(CM, 2, sum)
    sfs <- log(cell_total) / mean(log(cell_total))
  } else if (method == "mean-geometric-mean-total") {
    cell_total <- apply(CM, 2, sum)
    sfs <- cell_total / exp(mean(log(cell_total)))
  }

  sfs[is.na(sfs)] <- 1
  sfs
}

#' Find the most commonly occuring relative expression value in each cell
#'
#' Converting relative expression values to mRNA copies per cell requires
#' knowing the most commonly occuring relative expression value in each cell.
#' This value typically corresponds to an RPC value of 1. This function
#' finds the most commonly occuring (log-transformed) relative expression value
#' for each column in the provided expression matrix.
#'
#' @param relative_expr_matrix a matrix of relative expression values for
#' values with each row and column representing genes/isoforms and cells,
#' respectively. Row and column names should be included.
#' Expression values should not be log-transformed.
#' @param relative_expr_thresh Relative expression values below this threshold
#' are considered zero.
#'
#' @return an vector of most abundant relative_expr value corresponding to the RPC 1.
#' @details This function estimates the most abundant relative expression value
#' (t^*) using a gaussian kernel density function. It can also optionally
#' output the t^* based on a two gaussian mixture model based on the smsn.mixture from mixsmsn package.
#' @export
#'
#' @examples
#' \dontrun{
#' HSMM_fpkm_matrix <- exprs(HSMM)
#' t_estimate <- estimate_t(HSMM_fpkm_matrix)
#' }
estimate_t <- function(relative_expr_matrix, relative_expr_thresh = 0.1) {
  # apply each column
  unlist(apply(relative_expr_matrix, 2, function(relative_expr) 10^mean(dmode(log10(relative_expr[relative_expr > relative_expr_thresh]))))) # avoid multiple output
}

#' use gaussian kernel to calculate the mode of transcript counts
#'
#' @param x log tranformed relative expression
#' @param breaks breaks control parameter
#'
#' @importFrom stats density
dmode <- function(x, breaks = "Sturges") {
  if (length(x) < 2) {
    return(0)
  }
  den <- stats::density(x, kernel = c("gaussian"))
  (den$x[den$y == max(den$y)])
}

#' Create SCEG-HiC preprocessed input data
#'
#' Function to create preprocessed input data for SCEG-HIC.\code{calculateHiCWeights} and
#' \code{Run_SCEG_HiC}takes as input the generated data.This function can generate clustered data
#'  or non clustered data for a specific cell type.
#'
#' @param object Seurat object.
#' @param aggregate Logically,  should the data be clustered ?If FALSE, \code{celltype} must be specified as one cell type.
#' @param celltype Which cells to generate data
#' @param rna_assay Name of assay containing gene expression information.If \code{aggregate} is TRUE, Name of assay containing raw gene expression
#' information.If \code{aggregate} is FALSE, Name of assay containing normalize gene expression information.
#' @param atac_assay Name of assay containing peak information.If \code{aggregate} is TRUE, Name of assay containing raw peak information.
#'  If \code{aggregate} is FALSE, Name of assay containing normalize peak information.
#' @param cellnames Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities.
#' @param atacbinary Logical, whether the aggregated scATAC-seq data need binary?
#' @param k_neigh Number of cells to be aggregated per group.
#' @param max_overlap The maximum overlapping ratio of two groups.
#' @param size_factor_normalize Logical, should accessibility values be normalized by size factor?
#' @param reduction.name The reduction name of extracting the cell coordinates used for aggregating.
#' @param seed Random seed.
#' @param verbose Logical, should warning and info messages be printed?
#'
#' @import Signac
#' @import Seurat
#' @import Matrix
#' @return Preprocessed data for SCEG-HiC
#' @export
#'
#' @examples
#' data(multiomic_small)
#' # aggregated with paired scRNA-seq and scATAC-seq data
#' SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#'
#' # select "1" cells with paired scRNA-seq and scATAC-seq data
#' SCEGdata <- process_data(multiomic_small, aggregate = FALSE, celltype = "1")
#'
process_data <- function(object, aggregate = TRUE, celltype = NULL, rna_assay = "RNA", atac_assay = "peaks", cellnames = NULL, atacbinary = TRUE, k_neigh = 50, max_overlap = 0.8, size_factor_normalize = TRUE, reduction.name = NULL, seed = 123, verbose = TRUE) {
  if (aggregate) {
    message("Generating aggregated data")
    if ("aggregated_data" %in% names(Misc(object))) {
      agg.data <- Misc(object, slot = "aggregated_data")
    } else {
      agg.data <- Aggregate_datasets(object, k_neigh = k_neigh, rna_assay = rna_assay, atac_assay = atac_assay, atacbinary = atacbinary, cellnames = cellnames, max_overlap = max_overlap, reduction.name = reduction.name, size_factor_normalize = size_factor_normalize, seed = seed, verbose = verbose)
      Misc(object, slot = "aggregated_data") <- agg.data
      return(agg.data)
    }
  } else {
    if (is.null(celltype)) {
      stop("Please provide a specific cell type")
    }
    if (length(celltype) > 1) {
      stop("Please provide only one cell type")
    }
    if (!(is.null(cellnames))) {
      keep <- which(object@meta.data[[cellnames]] %in% celltype)
      dataset <- subset(object, cells = keep)
    } else {
      dataset <- subset(object, idents = celltype)
    }

    if (rna_assay %in% names(object@assays)) {
      rna <- dataset@assays[[rna_assay]]@data
      rna <- rna[rowSums(rna) != 0, ]
    }
    atac <- dataset@assays[[atac_assay]]@data
    atac <- atac[rowSums(atac) != 0, ]

    # binarize
    if (atacbinary) {
      atac@x[atac@x > 0] <- 1
    }

    singledata <- list()
    if (rna_assay %in% names(object@assays)) {
      singledata$rna <- rna
    }
    singledata$atac <- atac
    return(singledata)
  }
}
