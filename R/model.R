#' The bulk average Hi-C matrix normalization
#'
#' @param hicweight Raw the bulk average Hi-C matrix.
#' @param normalizeMethod The method of normalizing the average Hi-C.Default is "1".
#' * \code{"1"}, rank score to normalize the average Hi-C.
#' * \code{"2"}, min-max to normalize the average Hi-C.
#' * \code{"3"}, 0-1 to normalize the average Hi-C.
#' * \code{"4"}, -log10 to normalize the average Hi-C.
#' @param element The row or col elements of the average Hi-C matrix.
#'
#' @importFrom stats xtabs
#' @return A matrix with normalizing average Hi-C.
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
  if (!normalizeMethod %in% c("1", "2", "3", "4")) {
    stop("You must choose the normalized method provided by SCEG-HiC\n")
  }

  if (normalizeMethod == "1") {
    contacts <- hicweight[order(hicweight$score, decreasing = T), ]
    dd <- seq(0, 1, length.out = length(unique(hicweight$score)))
    dd1 <- as.matrix(table(contacts$score))
    contacts$score <- rep(dd, times = rev(dd1))
    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = contacts))
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix
    chrcontact1 <- chrcontact + t(chrcontact)
    diag(chrcontact1) <- 1
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
  } else if (normalizeMethod == "3") {
    xtabs_matrix <- as.matrix(xtabs(score ~ element1 + element2, data = hicweight))
    chrcontact <- matrix(0,
      nrow = length(element), ncol = length(element),
      dimnames = list(element, element)
    )

    chrcontact[rownames(xtabs_matrix), colnames(xtabs_matrix)] <- xtabs_matrix
    chrcontact1 <- chrcontact + t(chrcontact)
    chrcontact1[chrcontact1 > 0] <- 1
    chrcontact1 <- 1 - chrcontact1
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

#' Gene network reconstruction by integration of prior biological knowledge
#'
#' The method, instead of using a single penalty, different penalty values were applied for element pairs
#' based on a priori knowledge as to whether the two element should be connected. The a priori information
#' can be calculated or retrieved from other biological data, e.g., the bulk average Hi-C matrix.
#'
#' @param S Covariance matrix:p by p matrix (symmetric)
#' @param sample The number of samples.
#' @param normhicMatrix A matrix with normalizing average Hi-C.
#' @param element The row or col elements of the average Hi-C matrix.
#'
#' @return A matrix of partial correlation coefficients between elements.
#' @import glasso
#' @import dplyr
#' @export
#'
wglasso <- function(S, sample, normhicMatrix, element) {
  diag(S) <- diag(S) + 1e-4
  num <- length(element)
  rholist <- seq(from = 0.0, to = 1, by = 0.01)
  BICs <- NULL
  ENs <- NULL
  # wglasso, select rho based on BIC(Bayesian information criterion)
  for (i in 1:(length(rholist))) {
    g <- glasso(S, rholist[i] * normhicMatrix, nobs = num)
    mat <- g$wi
    diag(mat) <- 0
    mat[which(mat != 0)] <- 1
    edge.num <- sum(mat) / 2
    degree <- colSums(mat)
    ENs <- c(ENs, edge.num)
    lamda <- 1 - log(sample) / (4 * log(num))
    ln <- sample / 2 * (2 / num * g$loglik + sum(abs(g$wi)) * rholist[i])
    BIC <- -2 * ln + edge.num * log(sample)
    BICs <- c(BICs, BIC)
  }
  rho <- rholist[which.min(BICs)]
  g <- glasso(S, rho * normhicMatrix)
  mat <- g$wi
  colnames(mat) <- element
  rownames(mat) <- element
  diag_sqrt <- sqrt(diag(mat))
  mat1 <- -mat / (diag_sqrt %*% t(diag_sqrt))
  return(mat1)
}

#' Run SCEG-HiC's model
#'
#' The model uses data preprocessed by SCEG-HIC, the bulk average Hi-C as the penalty matrix,
#' calculates the partial correlation coefficient matrix between elements using the wglasso model,
#' and selects the gene-peak.
#'
#' @param SCEGdata Preprocessed data for SCEG-HiC.
#' @param HiCWeights Raw the bulk average Hi-C matrix.Run the result of code{'calculateHiCWeights'}
#' @param method The method of gene network reconstruction.Default is wglasso.
#' @param focus_gene The focused genes.
#' @param normalizeMethod The method of normalizing the average Hi-C.Default is "1".
#' * \code{'1'},-rank score to normalize the average Hi-C.
#' * \code{'2'},-min-max to normalize the average Hi-C.
#' * \code{'3'},-0-1 to normalize the average Hi-C.
#' * \code{'4'},--log10 to normalize the average Hi-C.
#' @param cutoff The threshold of gene-peak.Default is NULL.If \code{"aggregate"} is TRUE ,we recommend the cutoff of 0.01.
#' If \code{"aggregate"} is FALSE ,we recommend the cutoff of 0.001.
#' @param verbose Logical, should warning and info messages be printed?
#'
#' @importFrom stats sd cor
#' @return A data.frame of correlations for each gene-peak pair.
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
  if ("rna" %in% names(SCEGdata)) {
    rna <- SCEGdata[["rna"]]
    HiCWeights <- HiCWeights[which(names(HiCWeights) %in% rownames(rna))]
  }
  atac <- SCEGdata[["atac"]]
  rownames(atac) <- gsub("-", "_", rownames(atac))
  peak <- rownames(atac)

  focus_gene <- intersect(names(HiCWeights), focus_gene)
  if (verbose) {
    message(paste0("The predicted genes are ", length(focus_gene), " in total."))
  }
  SCEG_HiC_result <- list()
  for (n in 1:length(focus_gene)) {
    if (verbose) {
      message(paste0("Run model of  ", focus_gene[n]))
    }
    gene <- focus_gene[n]
    geneinfo <- HiCWeights[[gene]]
    if ("rna" %in% names(SCEGdata)) {
      Y <- as.matrix(rna[which(rownames(rna) == gene), ])
    } else {
      promoters <- geneinfo[["promoters"]]
      if (length(promoters) > 1) {
        Y <- colSums(atac[which(rownames(atac) %in% promoters), ])
      } else if (length(promoters) == 1) {
        Y <- atac[which(rownames(atac) %in% promoters), ]
      } else {
        next
      }
    }

    if (sd(Y) > 0) {
      enhancers <- geneinfo[["enhancers"]]
      if (length(enhancers) == 1) {
        X1 <- t(as.matrix(atac[which(rownames(atac) %in% enhancers), ]))
        if (dim(X1)[2] == 0) {
          next
        }
      } else {
        X1 <- as.matrix(atac[which(rownames(atac) %in% enhancers), ])
      }

      hicweight <- geneinfo[["contact"]]
      element <- c(gene, enhancers)
      if ((dim(X1)[1] == dim(atac)[2]) || (dim(X1)[1] == 0)) {
        next
      }
      if (dim(X1)[1] != length(enhancers)) {
        enhancers <- rownames(X1)
        element <- c(gene, enhancers)
        hicweight <- hicweight[which(hicweight$element1 %in% element), ]
        hicweight <- hicweight[which(hicweight$element2 %in% element), ]
      }

      genecontact <- normalizeHiCWeights(hicweight, normalizeMethod, element)

      X2 <- cbind(Y, t(X1))
      colnames(X2) <- element
      std_devs <- apply(X2, 2, sd)
      if (sum(std_devs == 0) > 0) {
        if (sum(std_devs == 0) == length(element) - 1) {
          next
        }
        X2 <- X2[, std_devs != 0]
        id <- which(element %in% colnames(X2))
        genecontact <- genecontact[id, id]
        element <- colnames(X2)
      }
      s <- cor(X2)
      result_matrix <- wglasso(s, dim(atac)[2], genecontact, element)

      if ("rna" %in% names(SCEGdata)) {
        SCEG_HiC_result[[n]] <- data.frame(gene = gene, peak = element[-1], score = result_matrix[, 1][-1])
      } else {
        SCEG_HiC_result[[n]] <- data.frame(gene = gene, peak1 = promoters[1], peak2 = element[-1], score = result_matrix[, 1][-1])
      }
    }
  }
  SCEG_HiC_result <- do.call(rbind, SCEG_HiC_result)
  if (!is.null(cutoff)) {
    SCEG_HiC_result <- SCEG_HiC_result[SCEG_HiC_result$score > cutoff, ]
  }

  return(SCEG_HiC_result)
}
