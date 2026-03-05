#' Annotate transcription start sites (TSS)
#'
#' This function annotates transcription start sites based on the specified species and genome assembly version.
#'
#' @param species Character string specifying the species name. Supported values are "Homo sapiens" or "Mus musculus".
#' @param genome Character string specifying the genome assembly. Supported values are "hg38", "hg19", "mm10", or "mm9".
#'
#' @return A data.frame containing TSS annotations filtered by species and genome, including chromosome,
#' gene name, and TSS genomic loci information.
#' @export
#'
#' @examples
#' # TSS annotations for the hg38
#' tssdata <- annotateTSS("Homo sapiens", "hg38")
#' head(tssdata)
#'
#' # TSS annotations for the mm10
#' tssdata <- annotateTSS("Mus musculus", "mm10")
#' head(tssdata)
annotateTSS <- function(species, genome) {
  # Validate input parameters
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    stop("Species must be 'Homo sapiens' or 'Mus musculus'")
  }
  if (!genome %in% c("hg19", "hg38", "mm10", "mm9")) {
    stop("Genome must be one of 'hg19', 'hg38', 'mm10', or 'mm9'")
  }

  # Load TSS data from package internal datasets
  switch(genome,
    hg19 = {
      tssdata <- SCEGHiC::hg19TSSdata
    },
    hg38 = {
      tssdata <- SCEGHiC::hg38TSSdata
    },
    mm10 = {
      tssdata <- SCEGHiC::mm10TSSdata
    },
    mm9 = {
      tssdata <- SCEGHiC::mm9TSSdata
    }
  )

  genes <- lapply(tssdata$TargetGene, function(x) strsplit(x, "[|]")[[1]][1])
  genes <- lapply(genes, function(x) strsplit(x, "[.]")[[1]][1])
  genes <- unlist(genes)
  tssdata$TargetGene <- genes

  # Remove duplicated genes, keep first occurrence
  tssdata <- tssdata[!duplicated(tssdata$TargetGene), ]

  # Filter chromosomes based on species:
  # Human: chr1-22 and chrX
  # Mouse: chr1-19 and chrX
  valid_chromosomes <- if (species == "Homo sapiens") {
    paste0("chr", c(1:22, "X"))
  } else {
    paste0("chr", c(1:19, "X"))
  }

  tssdata <- tssdata[tssdata$chr %in% valid_chromosomes, ]

  return(tssdata)
}

#' Convert genomic coordinates of peaks between genome assemblies
#'
#' This function converts genomic coordinates for a set of peaks to a specified genome assembly.
#' The \code{SCEGHiC} package provides bulk average Hi-C data for the hg38 and mm10 assemblies,
#' so coordinate conversion is required for other assemblies such as hg19 and mm9.
#'
#' @param split_peaks A matrix or data.frame with separated peak coordinates, typically obtained by splitting peak strings (e.g., "chr2-973888-974755") into columns for chromosome, start, and end.
#' @param genome Character string specifying the genome assembly. Supported values are "hg19" or "mm9".
#'
#' @importFrom  IRanges IRanges
#' @importFrom  GenomicRanges GRanges
#' @importFrom  S4Vectors Rle
#' @import rtracklayer
#' @return A data.frame containing peaks after genome coordinate conversion, with updated chromosome, start, and end positions.
#' @export
#'
#' @examples
#' peaks <- c(
#'   "chr2-973888-974755", "chr2-3195868-3196264", "chr2-3242398-3243508", "chr2-6863856-6866834", "chr2-6913141-6914394",
#'   "chr2-6924663-6924972", "chr2-6999095-6999554", "chr2-7007713-7008819", "chr2-7083945-7084032", "chr2-7097148-7098029"
#' )
#' split_peaks <- do.call("rbind", strsplit(as.character(peaks), "-", fixed = TRUE))
#' peakinfo <- liftOverPeaks(split_peaks, "hg19")
#' head(peakinfo)
#'
liftOverPeaks <- function(split_peaks, genome) {
  # Check input
  if (ncol(split_peaks) != 3) {
    stop("Input split_peaks must have  3 columns: chromosome, start, end")
  }

  chep <- data.frame(R1.chrom = split_peaks[, 1], R1.start = as.numeric(split_peaks[, 2]), R1.end = as.numeric(split_peaks[, 3]), stringsAsFactors = FALSE)

  # Determine chain file path based on genome assembly
  if (grepl("hg", genome)) {
    f <- system.file("extdata", "hg19ToHg38.over.chain", package = "SCEGHiC")
  } else if (genome == "mm9") {
    f <- system.file("extdata", "mm9ToMm10.over.chain", package = "SCEGHiC")
  } else {
    stop("Unsupported genome version. Supported: 'hg19', 'mm9'")
  }

  # Import liftover chain
  liftover.chain <- rtracklayer::import.chain(f)

  # Create GRanges object from peaks
  chep$line <- 1:nrow(chep)
  chepr1 <- with(chep, GRanges(seqnames = Rle(R1.chrom), ranges = IRanges(start = R1.start, end = R1.end), id = line))

  # Perform liftover
  chepr1.lift <- unlist(liftOver(chepr1, liftover.chain))

  # Construct data.frame from liftover results
  df <- data.frame(iranges = chepr1.lift)
  colnames(chep) <- c("R1.chrom", "R1.start", "R1.end", "iranges.id")

  # Merge lifted data with original peaks on peak_id
  datas <- merge(df, chep, by = "iranges.id", all.x = T)

  # Filter to keep only those where lifted chromosome matches original chromosome (optional)
  datas <- datas[datas$iranges.seqnames == datas$R1.chrom, ]

  # Remove duplicated peaks (keep first occurrence)
  datas1 <- datas[which(!duplicated(datas$iranges.id)), ]

  # Create output with new peak coordinates and original peak string
  liftov <- data.frame(chr = datas1$iranges.seqnames, start = datas1$iranges.start, end = datas1$iranges.end, peak = paste(datas1$R1.chrom, datas1$R1.start, datas1$R1.end, sep = "_"))
  liftov$midpoint <- (as.numeric(liftov$start) + as.numeric(liftov$end)) / 2

  return(liftov)
}

# Declaration a global variable
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}

#' Calculate promoters, enhancers, and Hi-C matrix for genes
#'
#' This function obtains enhancers and promoters of each gene, and calculates the bulk average Hi-C matrix of contacts between the gene and its enhancers.
#'
#' @param SCEGdata Preprocessed data for SCEG-HiC.
#' @param species Character string specifying the species name. Supported values are "Homo sapiens" or "Mus musculus".
#' @param genome Character string specifying the genome assembly. Supported values are "hg38", "hg19", "mm10", or "mm9".
#' @param focus_gene Character vector of gene names to focus on.
#' @param averHicPath Path to the bulk average Hi-C data.
#' @param TSSwindow Numeric specifying the number of base pairs to extend upstream and downstream around each TSS to define promoters. Default is 1000 bp (total window size of 2000 bp).
#' @param upstream Numeric specifying the number of base pairs upstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param downstream Numeric specifying the number of base pairs downstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param verbose Logical. Should progress messages and warnings be printed?
#'
#' @importFrom cicero find_overlapping_coordinates
#' @importFrom utils combn
#' @importFrom  data.table fread
#' @import reader
#' @import R.utils
#' @return A named list where each element corresponds to a focus gene and contains:
#' \itemize{
#'   \item \code{promoters}: data.frame of promoter regions.
#'   \item \code{enhancers}: data.frame of enhancer regions.
#'   \item \code{contact}: bulk average Hi-C contact matrix between the gene and enhancers.
#' }
#' @export
#'
#' @examples
#' data(multiomic_small)
#' SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#' fpath <- system.file("extdata", package = "SCEGHiC")
#' gene <- c("TRABD2A", "GNLY", "MFSD6", "CTLA4", "LCLAT1", "NCK2", "GALM", "TMSB10", "ID2", "CXCR4")
#' weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
#' weight[["CTLA4"]]$promoters
#' weight[["CTLA4"]]$enhancers
#' head(weight[["CTLA4"]]$contact)
#'
calculateHiCWeights <- function(SCEGdata, species, genome, focus_gene, averHicPath, TSSwindow = 1000, upstream = 250000, downstream = 250000, verbose = TRUE) {
  # Validate input parameters
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    stop("Species must be 'Homo sapiens' or 'Mus musculus'")
  }
  if (!genome %in% c("hg19", "hg38", "mm10", "mm9")) {
    stop("Genome must be one of 'hg19', 'hg38', 'mm10', or 'mm9'")
  }

  if (is.null(focus_gene) || length(focus_gene) == 0) {
    stop("Please provide at least one focus gene.")
  }

  all.peaks <- rownames(SCEGdata$atac)
  split_peaks <- do.call("rbind", strsplit(as.character(all.peaks), "-", fixed = TRUE))

  # For hg19 and mm9, perform liftover; else just parse coordinates
  if (genome %in% c("hg19", "mm9")) {
    peakinfo <- liftOverPeaks(split_peaks, genome)
  } else {
    peakinfo <- data.frame(chr = split_peaks[, 1], start = as.numeric(split_peaks[, 2]), end = as.numeric(split_peaks[, 3]), stringsAsFactors = FALSE)
    peakinfo$peak <- paste(peakinfo$chr, peakinfo$start, peakinfo$end, sep = "_")
    peakinfo$midpoint <- (as.numeric(peakinfo$start) + as.numeric(peakinfo$end)) / 2
  }
  all.peaks <- peakinfo$peak

  # Get TSS annotations for the focus genes
  tssdata <- annotateTSS(species, genome)
  genome.info <- tssdata[which(tssdata$TargetGene %in% focus_gene), ]

  if (nrow(genome.info) == 0) {
    stop("Failed to retrieve TSS info for focus genes. Check gene names and genome version.")
  }

  chromosome <- as.character(unique(genome.info$chr))


  HiCWeights <- list()

  # Loop through chromosomes containing focus genes
  for (i in chromosome) {
    if (verbose) message(paste0("Processing chromosome ", i, "..."))

    # Construct file path to average Hi-C file and check existence
    file_path <- paste0(averHicPath, "/", i, "/", i, ".bed.gz") # ifelse(species=="Homo sapiens",".bed.gz",".avg.gz")
    if (!file.exists(file_path)) {
      stop(paste0("Hi-C file not found: ", file_path))
    }

    averHiC <- fread(file_path)
    colnames(averHiC) <- c("bin1", "bin2", "score")

    genome.info.used <- genome.info[genome.info$chr == i, ]
    focus_markers <- genome.info.used$TargetGene
    Starts <- genome.info.used$TargetGeneTSS

    # Use bulk Hi-C TSS annotations for weighting (hg38/mm10 depending on species)
    if (species == "Homo sapiens") {
      weight_tssdata <- annotateTSS(species, "hg38")
    } else {
      weight_tssdata <- annotateTSS(species, "mm10")
    }
    weight_Starts <- weight_tssdata[match(focus_markers, weight_tssdata$TargetGene), ]$TargetGeneTSS

    if (length(focus_markers) == 0) {
      warning(paste0("No focus gene TSS found on ", i, ". Skipping."))
      next
    }

    if (verbose) message(paste0("Found ", length(focus_markers), " TSS loci on ", i, "."))

    for (n in 1:length(focus_markers)) {
      if (verbose) message(paste0("Calculating Hi-C weights for gene ", focus_markers[n], "..."))

      # Define promoter and enhancer search regions around TSS
      p1 <- paste(i, ":", Starts[n] - TSSwindow, "-", Starts[n] + TSSwindow, sep = "")
      p2 <- paste(i, ":", Starts[n] - upstream, "-", Starts[n] + downstream, sep = "")
      promoters <- find_overlapping_coordinates(all.peaks, p1)
      enhancers <- find_overlapping_coordinates(all.peaks, p2)
      enhancers <- setdiff(enhancers, promoters)

      # Return empty results if no enhancers found
      if (length(enhancers) > 0) {
        peakinfo1 <- peakinfo[which(peakinfo$peak %in% enhancers), ]

        # Calculate bins at 5kb resolution for gene and enhancers
        element <- c(focus_markers[n], enhancers)
        pairs <- as.data.frame((t(combn(element, 2))))
        colnames(pairs) <- c("element1", "element2")
        bin1 <- floor(rep(c(weight_Starts[n], peakinfo1$midpoint[-length(enhancers)]), length(enhancers):1) / 5000) * 5000
        indices <- unlist(lapply(1:length(enhancers), function(x) x:length(enhancers)))
        bin2 <- floor(c(peakinfo1$midpoint[indices]) / 5000) * 5000
        pairs$bin1 <- pmin(bin1, bin2)
        pairs$bin2 <- pmax(bin1, bin2)
        hicweight <- averHiC[pairs, on = .(bin1, bin2), nomatch = 0]
        hicweight$score[is.na(hicweight$score)] <- 0
        hicweight1 <- hicweight[, c(4, 5, 3)]
        result <- list()
        result[["promoters"]] <- promoters
        result[["enhancers"]] <- enhancers
        result[["contact"]] <- hicweight1
        HiCWeights[[focus_markers[n]]] <- result
      }
    }
  }

  if (verbose) message("Finished calculating Hi-C weights for all genes.")

  return(HiCWeights)
}
