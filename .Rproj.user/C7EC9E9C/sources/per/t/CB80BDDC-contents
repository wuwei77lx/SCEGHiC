#' Annotate the provided Transcription Start Sites
#'
#' This function annotates transcription start sites based on species and genome versions
#' @param species String indicating the name of the species.It can be either "Homo sapiens","Mus musculus".
#' @param genome String indicating the name of the genome assembly.It can be either "hg38",
#' "hg19", "mm10", "mm9".
#'
#' @return A data.frame with each gene TSS locus
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
  unik <- !duplicated(genes)
  tssdata <- tssdata[unik, ]
  if (species == "Homo sapiens") {
    CHR <- c()
    for (i in c(1:22, "X")) {
      CHR <- c(CHR, paste0("chr", i))
    }
    tssdata <- tssdata[which(tssdata$chr %in% CHR), ]
  } else {
    CHR <- c()
    for (i in c(1:19, "X")) {
      CHR <- c(CHR, paste0("chr", i))
    }
    tssdata <- tssdata[which(tssdata$chr %in% CHR), ]
  }
  return(tssdata)
}

#' Peaks corresponding genome version conversion
#'
#' The function convert genomic coordinates based on provided peaks.\code{SCEGHIC} provides the bulk
#' average Hi-C genome versions of hg19 and mm10, so   require genome coordinate conversion for hg38 and mm10.
#' @param split_peaks Peaks of genome coordinate conversion.
#' @param genome String indicating the name of the genome assembly.It can be either "hg38","mm9".
#'
#' @importFrom  IRanges IRanges
#' @importFrom  GenomicRanges GRanges
#' @importFrom  S4Vectors Rle
#' @import rtracklayer
#' @return A data.frame with peaks after genome coordinate conversion.
#' @export
#'
#' @examples
#' peaks <- c(
#'   "chr2-973888-974755", "chr2-3195868-3196264", "chr2-3242398-3243508", "chr2-6863856-6866834", "chr2-6913141-6914394",
#'   "chr2-6924663-6924972", "chr2-6999095-6999554", "chr2-7007713-7008819", "chr2-7083945-7084032", "chr2-7097148-7098029"
#' )
#' split_peaks <- do.call("rbind", strsplit(as.character(peaks), "-", fixed = TRUE))
#' peakinfo <- liftOverPeaks(split_peaks, "hg38")
#' head(peakinfo)
liftOverPeaks <- function(split_peaks, genome) {
  chep <- data.frame(R1.chrom = split_peaks[, 1], R1.start = as.numeric(split_peaks[, 2]), R1.end = as.numeric(split_peaks[, 3]))

  if (grepl("hg", genome)) {
    f <- system.file("extdata", "hg38ToHg19.over.chain", package = "SCEGHiC")
  } else {
    f <- system.file("extdata", "mm9ToMm10.over.chain", package = "SCEGHiC")
  }
  f <- sub(".gz", "", f)

  liftover.chain <- import.chain(f)
  chep$line <- 1:nrow(chep)
  chepr1 <- with(chep, GRanges(seqnames = Rle(R1.chrom), ranges = IRanges(start = R1.start, end = R1.end), id = line))
  chepr1.lift <- unlist(liftOver(chepr1, liftover.chain))
  df <- data.frame(iranges = chepr1.lift)
  colnames(chep) <- c("R1.chrom", "R1.start", "R1.end", "iranges.id")
  datas <- merge(df, chep, by = "iranges.id", all.x = T)
  datas <- datas[datas$iranges.seqnames == datas$R1.chrom, ]
  datas1 <- datas[which(!duplicated(datas$iranges.id)), ]
  liftov <- data.frame(chr = datas1$iranges.seqnames, start = datas1$iranges.start, end = datas1$iranges.end, peak = paste(datas1$R1.chrom, datas1$R1.start, datas1$R1.end, sep = "_"))
  liftov$midpoint <- (as.numeric(liftov$start) + as.numeric(liftov$end)) / 2

  return(liftov)
}

# Declaration a global variable
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("."))
}
#' Calculate  promoters, enhancers, and Hi-C matrix of genes
#'
#' The function obtains enhancers and promoters of each gene, and calculates the bulk average Hi-C matrix of the contact between the gene and enhancers.
#'
#' @param SCEGdata Preprocessed data for SCEG-HiC
#' @param species String indicating the name of the species.It can be either "Homo sapiens","Mus musculus".
#' @param genome String indicating the name of the genome assembly.It can be either "hg38",
#' "hg19", "mm10", "mm9".
#' @param focus_gene The focused genes.
#' @param averHicPath Path to the bulk average Hi-C.
#' @param TSSwindow Numeric specifying the window size (in base pairs) to pad around either side of each TSS, to fetch gene's promoters.
#' Default is 1000 bp (so 2000 bp window is drawn arund each TSS).
#' @param upstream Numeric specifying the window size (in base pairs) to pad around upstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param downstream Numeric specifying the window size (in base pairs) to pad around downstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param verbose Logical, should warning and info messages be printed?
#'
#' @importFrom cicero find_overlapping_coordinates
#' @importFrom utils combn
#' @importFrom  data.table fread
#' @import reader
#' @import R.utils
#' @return A list with enhancers, promoters, and Hi-C matrices for each \code{focus_gene}
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
calculateHiCWeights <- function(SCEGdata, species, genome, focus_gene, averHicPath, TSSwindow = 1000, upstream = 250000, downstream = 250000, verbose = TRUE) {
  if (!genome %in% c("hg19", "hg38", "mm10", "mm9")) {
    stop("You must specify one of hg19, hg38 , mm10 or mm9 as a genome build for currently supported TSS annotations..\n")
  }
  if (!species %in% c("Homo sapiens", "Mus musculus")) {
    stop("SCEG-HiC only provides predictions for human and mouse species\n")
  }

  all.peaks <- rownames(SCEGdata$atac)
  split_peaks <- do.call("rbind", strsplit(as.character(all.peaks), "-", fixed = TRUE))
  if (genome %in% c("hg38", "mm9")) {
    peakinfo <- liftOverPeaks(split_peaks, genome)
  } else {
    peakinfo <- data.frame(chr = split_peaks[, 1], start = as.numeric(split_peaks[, 2]), end = as.numeric(split_peaks[, 3]))
    peakinfo$peak <- paste(peakinfo$chr, peakinfo$start, peakinfo$end, sep = "_")
    peakinfo$midpoint <- (as.numeric(peakinfo$start) + as.numeric(peakinfo$end)) / 2
  }

  all.peaks <- peakinfo$peak
  tssdata <- annotateTSS(species, genome)
  genome.info <- tssdata[which(tssdata$TargetGene %in% focus_gene), ]
  chromosome <- as.character(unique(genome.info$chr))

  if (length(chromosome) == 0) {
    stop("Failed to obtain any TSS locus of gene..\n")
  }

  HiCWeights <- list()
  for (i in chromosome) {
    file_path <- paste0(averHicPath, "/", i, "/", i, ".bed.gz") # ifelse(species=="Homo sapiens",".bed.gz",".avg.gz")
    averHiC <- fread(file_path)
    colnames(averHiC) <- c("bin1", "bin2", "score")
    genome.info.used <- genome.info[genome.info$chr == i, ]
    focus_markers <- genome.info.used$TargetGene
    Starts <- genome.info.used$TargetGeneTSS

    # bulk average Hi-C genome versions of hg19 and mm10
    if (species == "Homo sapiens") {
      weight_tssdata <- annotateTSS(species, "hg19")
      weight_Starts <- weight_tssdata[match(focus_markers, weight_tssdata$TargetGene), ]$TargetGeneTSS
    } else {
      weight_tssdata <- annotateTSS(species, "mm10")
      weight_Starts <- weight_tssdata[match(focus_markers, weight_tssdata$TargetGene), ]$TargetGeneTSS
    }

    if (length(focus_markers) == 0) {
      stop(paste0("Failed to obtain any TSS locus of genes from chromosome ", sub("chr", "", i), "..\n"))
    } else {
      message(paste0("Successfully obtained ", length(focus_markers), " TSS loci of genes from chromosome ", sub("chr", "", i), "."))
    }

    for (n in 1:length(focus_markers)) {
      if (verbose) {
        message(paste0("Calculate weight of  ", focus_markers[n]))
      }
      p1 <- paste(i, ":", Starts[n] - TSSwindow, "-", Starts[n] + TSSwindow, sep = "")
      p2 <- paste(i, ":", Starts[n] - upstream, "-", Starts[n] + downstream, sep = "")
      promoters <- find_overlapping_coordinates(all.peaks, p1)
      enhancers <- find_overlapping_coordinates(all.peaks, p2)
      enhancers <- setdiff(enhancers, promoters)
      if (length(enhancers) > 0) {
        peakinfo1 <- peakinfo[which(peakinfo$peak %in% enhancers), ]
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
  return(HiCWeights)
}
