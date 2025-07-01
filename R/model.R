#' Normalize the bulk average Hi-C matrix
#'
#' This function normalizes the raw bulk average Hi-C contact scores using different methods.
#' @param hicweight A data.frame containing raw bulk average Hi-C contacts with columns: element1, element2, and score.
#' @param normalizeMethod Method used to normalize the average Hi-C scores. Default is "1".
#'   * \code{"1"}: Normalize by rank scores.
#'   * \code{"2"}: Normalize by -log10 transformation.
#'   * \code{"3"}: Binarize the scores (values > 0 become 1; values ≤ 0 become 0).
#'   * \code{"4"}: Normalize by min-max scaling.
#' @param element Vector of unique row or column elements present in the Hi-C matrix.
#'
#' @importFrom stats xtabs ecdf
#' @return A normalized Hi-C contact matrix.
#' @export
#'
#' @examples
#' hicweight <- data.frame(
#'   element1 = c(
#'     "CHN1", "CHN1", "CHN1", "CHN1", "CHN1", "chr2_174901074_174901475", "chr2_174901074_174901475",
#'     "chr2_174901074_174901475", "chr2_174901074_174901475", "chr2_174990102_174990619", "chr2_174990102_174990619",
#'     "chr2_174990102_174990619", "chr2_175048293_175048704", "chr2_175048293_175048704", "chr2_175167038_175169195"
#'   ),
#'   element2 = c(
#'     "chr2_174901074_174901475", "chr2_174990102_174990619", "chr2_175048293_175048704", "chr2_175167038_175169195",
#'     "chr2_175187467_175187717", "chr2_174990102_174990619", "chr2_175048293_175048704", "chr2_175167038_175169195",
#'     "chr2_175187467_175187717", "chr2_175048293_175048704", "chr2_175167038_175169195", "chr2_175187467_175187717",
#'     "chr2_175167038_175169195", "chr2_175187467_175187717", "chr2_175187467_175187717"
#'   ),
#'   score = c(
#'     0.0016062358, 0.0076253167, 0.0031021570, 0.0012148262, 0.0010788392, 0.0019968580, 0.0012957845,
#'     0.0009198567, 0.0011448870, 0.0023991885, 0.0010785481, 0.0009664508, 0.0012865316, 0.0012122647, 0.0057781939
#'   )
#' )
#' element <- unique(c(hicweight$element1, hicweight$element2))
#' weight <- normalizeHiCWeights(hicweight, "1", element)
#' head(weight)
normalizeHiCWeights <- function(hicweight, normalizeMethod = "1", element) {
  # Validate the normalization method input
  if (!normalizeMethod %in% c("1", "2", "3", "4")) {
    stop("You must choose the normalized method provided by SCEG-HiC\n")
  }

  # Method 1: Rank score normalization
  if (normalizeMethod == "1") {
    contacts <- hicweight

    # Rank scores in descending order, then normalize to 0–1 range
    contacts$score <- (rank(-contacts$score, ties.method = "min") - 1) / (nrow(contacts) - 1)

    # Convert long-format data into a matrix form
    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = contacts))

    # Initialize a zero matrix with all elements
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    # Fill in non-zero scores
    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix

    # Make the matrix symmetric
    chrcontact1 <- chrcontact + t(chrcontact)

    # Set the diagonal to 1 (self-contact)
    diag(chrcontact1) <- 1

    # Method 2: Negative log transform followed by min-max normalization
  } else if (normalizeMethod == "2") {
    contacts <- hicweight
    contacts$score <- -log(contacts$score)
    contacts$score <- (contacts$score - min(contacts$score)) / (max(contacts$score) - min(contacts$score))

    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = contacts))
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix
    chrcontact1 <- chrcontact + t(chrcontact)
    diag(chrcontact1) <- 1

    # Method 3: Binarization (presence/absence of contact), then invert
  } else if (normalizeMethod == "3") {
    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = hicweight))
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix
    chrcontact1 <- chrcontact + t(chrcontact)

    # Set all non-zero contacts to 1 (binarization)
    chrcontact1[chrcontact1 > 0] <- 1

    # Invert: 1 becomes 0 (contact), 0 becomes 1 (no contact)
    chrcontact1 <- 1 - chrcontact1

    # Method 4: Reversed min-max normalization
  } else {
    contacts <- hicweight
    contacts$score <- (max(contacts$score) - contacts$score) / (max(contacts$score) - min(contacts$score))

    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = contacts))
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix
    chrcontact1 <- chrcontact + t(chrcontact)
    diag(chrcontact1) <- 1
  }

  return(chrcontact1)
}

#' Gene network reconstruction using prior biological knowledge
#'
#' This function reconstructs a gene (or regulatory element) network by integrating prior biological knowledge into the estimation of the precision matrix.
#' Instead of applying a uniform penalty to all element pairs, this method assigns different penalty values based on prior knowledge regarding whether two elements should be connected.
#' Such prior information can be derived or obtained from biological data, such as a normalized bulk average Hi-C matrix.
#'
#' @param S A symmetric p-by-p covariance matrix representing co-variation between elements.
#' @param sample An integer indicating the number of samples used to estimate the covariance matrix.
#' @param normhicMatrix A normalized average Hi-C contact matrix used as prior biological knowledge for penalization.
#' @param element A character vector containing the row/column names (element IDs) of the Hi-C matrix.
#'
#' @return A symmetric matrix of partial correlation coefficients between elements.
#' @import glasso
#' @import dplyr
#' @export
#'
wglasso <- function(S, sample, normhicMatrix, element) {
  # Add small constant to the diagonal for numerical stability
  diag(S) <- diag(S) + 1e-4

  num <- length(element)
  rholist <- seq(from = 0, to = 1, by = 0.01)
  BICs <- NULL
  ENs <- NULL

  # Select optimal rho by minimizing Bayesian Information Criterion (BIC)
  for (i in 1:(length(rholist))) {
    g <- glasso(S, rholist[i] * normhicMatrix, nobs = sample)

    # Extract estimated inverse covariance (precision) matrix
    mat <- g$wi
    diag(mat) <- 0

    # Binary adjacency: 1 if edge exists, 0 otherwise
    mat[which(mat != 0)] <- 1
    edge.num <- sum(mat) / 2
    degree <- colSums(mat)
    ENs <- c(ENs, edge.num)

    # Compute modified log-likelihood and BIC
    lamda <- 1 - log(sample) / (4 * log(num))
    ln <- sample / 2 * (2 / sample * g$loglik + sum(abs(g$wi)) * rholist[i])
    BIC <- -2 * ln + edge.num * log(sample)
    BICs <- c(BICs, BIC)
  }

  # Choose rho with lowest BIC
  rho <- rholist[which.min(BICs)]
  g <- glasso(S, rho * normhicMatrix)
  mat <- g$wi
  colnames(mat) <- element
  rownames(mat) <- element

  # Convert precision matrix to partial correlation matrix
  diag_sqrt <- sqrt(diag(mat))
  mat1 <- -mat / (diag_sqrt %*% t(diag_sqrt))
  return(mat1)
}

#' Run the SCEG-HiC Model
#'
#' This function runs the SCEG-HiC model to infer gene-enhancer (gene-peak) links using preprocessed
#' single-cell multi-omics data and bulk average Hi-C data as a prior. It calculates a partial correlation
#' matrix via the weighted graphical lasso (wglasso) approach and selects gene-peak pairs accordingly.
#'
#' @param SCEGdata A list containing data preprocessed by \code{SCEG-HiC}.
#' @param HiCWeights  A data frame of raw bulk average Hi-C interaction weights.
#' Output of \code{calculateHiCWeights()}.
#' @param method The method used for gene network reconstruction. Default is \code{"wglasso"}.
#' @param focus_gene A character vector of gene symbols to focus on.
#' @param normalizeMethod Method used to normalize the average Hi-C scores. Default is "1".
#'   * \code{"1"}: Normalize by rank scores.
#'   * \code{"2"}: Normalize by -log10 transformation.
#'   * \code{"3"}: Binarize the scores (values > 0 become 1; values ≤ 0 become 0).
#'   * \code{"4"}: Normalize by min-max scaling.
#' @param cutoff Threshold for selecting gene-peak pairs. Default is \code{NULL}.
#' If \code{aggregate = TRUE}, we recommend setting \code{cutoff = 0.01}.
#' If \code{aggregate = FALSE}, we recommend \code{cutoff = 0.001}.
#' @param verbose Logical. Should progress messages and warnings be printed?
#'
#' @importFrom stats sd cor
#' @return A \code{data.frame} containing correlation values and metadata for each gene-peak pair.
#' @export
#'
#' @examples
#' data(multiomic_small)
#' SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#' fpath <- system.file("extdata", package = "SCEGHiC")
#' gene <- c("TRABD2A", "GNLY", "MFSD6", "CTLA4", "LCLAT1", "NCK2", "GALM", "TMSB10", "ID2", "CXCR4")
#' weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
#' results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
#' results_SCEGHiC[results_SCEGHiC$gene == "CTLA4", ]
#'
Run_SCEG_HiC <- function(SCEGdata, HiCWeights, method = "wglasso", focus_gene, normalizeMethod = "1", cutoff = NULL, verbose = TRUE) {
  # Extract ATAC data and ensure unique row names
  atac <- SCEGdata[["atac"]]
  rownames(atac) <- gsub("-", "_", rownames(atac))
  peak <- rownames(atac)

  # If RNA data available, extract and check dimension consistency
  if ("rna" %in% names(SCEGdata)) {
    rna <- SCEGdata[["rna"]]
    if (ncol(rna) != ncol(atac)) {
      stop("RNA and ATAC matrices must have the same number of columns (cells).")
    }

    # Filter HiCWeights to genes present in RNA data
    HiCWeights <- HiCWeights[which(names(HiCWeights) %in% rownames(rna))]
  }

  # Filter focus genes to those with HiCWeights available
  focus_gene <- intersect(names(HiCWeights), focus_gene)
  if (verbose) {
    message(paste0("Total predicted genes: ", length(focus_gene)))
  }

  results_list <- list()
  for (gene in focus_gene) {
    if (verbose) message(paste0("Running model for gene: ", gene))

    geneinfo <- HiCWeights[[gene]]

    # Construct response variable Y
    if ("rna" %in% names(SCEGdata)) {
      Y <- as.matrix(rna[which(rownames(rna) == gene), ])
    } else {
      promoters <- geneinfo[["promoters"]]
      matched_promoters <- intersect(promoters, rownames(atac))
      if (length(matched_promoters) == 0) next

      if (length(matched_promoters) == 1) {
        Y <- as.numeric(atac[matched_promoters, ])
      } else {
        Y <- colSums(atac[matched_promoters, ])
      }
    }

    # Construct predictor matrix X1 for enhancers
    if (sd(Y) > 0) {
      enhancers <- geneinfo[["enhancers"]]
      matched_enhancers <- intersect(enhancers, rownames(atac))
      if (length(matched_enhancers) == 0) next

      X1 <- t(as.matrix(atac[matched_enhancers, , drop = FALSE]))
      colnames(X1) <- matched_enhancers

      # Extract Hi-C contact matrix subset for gene and enhancers
      element <- c(gene, matched_enhancers)
      hicweight <- geneinfo[["contact"]]
      hicweight <- hicweight[hicweight$element1 %in% element & hicweight$element2 %in% element, ]

      if (nrow(hicweight) < 2) {
        if (verbose) message(paste0("Gene ", gene, " has too few Hi-C contact entries (<2). Skipping."))
        next
      }

      # Normalize Hi-C weights
      genecontact <- normalizeHiCWeights(hicweight, normalizeMethod, element)


      # Combine Y and X1 into a single matrix X2
      X2 <- cbind(Y, X1)
      colnames(X2)[1] <- gene

      # Remove variables with zero standard deviation
      std_devs <- apply(X2, 2, sd)
      zero_sd_vars <- std_devs == 0
      if (sum(zero_sd_vars) > 0) {
        # If almost all variables have zero variance, skip this gene
        if (sum(zero_sd_vars) >= ncol(X2) - 1) next

        X2 <- X2[, !zero_sd_vars, drop = FALSE]
        element <- colnames(X2)
        id <- match(element, colnames(genecontact))
        id <- id[!is.na(id)]
        genecontact <- genecontact[id, id, drop = FALSE]
      }
      # Calculate correlation matrix
      s <- cor(X2)

      # Run model based on method
      if (method == "wglasso") {
        result_matrix <- wglasso(s, sample = ncol(atac), normhicMatrix = genecontact, element = colnames(X2))
      } else {
        stop("Unsupported method. Currently only 'wglasso' is supported.")
      }

      # Format results into a data.frame
      if ("rna" %in% names(SCEGdata)) {
        result_df <- data.frame(
          gene = gene,
          peak = colnames(X2)[-1],
          score = result_matrix[-1, 1]
        )
      } else {
        promoters <- geneinfo[["promoters"]]
        result_df <- data.frame(
          gene = gene,
          peak1 = promoters[1],
          peak2 = colnames(X2)[-1],
          score = result_matrix[-1, 1]
        )
      }

      results_list[[length(results_list) + 1]] <- result_df
    }
  }

  # Combine all gene results
  results_list <- results_list[!sapply(results_list, is.null)]
  if (length(results_list) == 0) {
    return(NULL)
  }
  final_result <- do.call(rbind, results_list)

  # Apply score cutoff if provided
  if (!is.null(cutoff)) {
    final_result <- final_result[final_result$score > cutoff, ]
  }

  return(final_result)
}
