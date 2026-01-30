#' Plot SCEG-HiC predicted enhancer-gene links
#'
#' Visualizes the predicted enhancer-gene links from SCEG-HiC for a given gene by plotting them in a genome region view.
#'
#' @param SCEG_HiC_Result A data.frame containing the output from \code{Run_SCEG_HiC()}.
#' @param species Character string specifying the species name. Supported values are "Homo sapiens" or "Mus musculus".
#' @param genome Character string specifying the genome assembly. Supported values are "hg38", "hg19", "mm10", or "mm9".
#' @param focus_gene A character vector of gene symbols to focus on.
#' @param cutoff Threshold for selecting gene-peak pairs. Default is \code{NULL}.
#' If \code{aggregate = TRUE}, we recommend setting \code{cutoff = 0.01}.
#' If \code{aggregate = FALSE}, we recommend \code{cutoff = 0.001}.
#' @param upstream Numeric specifying the number of base pairs upstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param downstream Numeric specifying the number of base pairs downstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param gene_anno A \code{data.frame} or \code{NULL}. If provided, should be compatible with the \code{\link[Gviz]{GeneRegionTrack-class}}.
#'   Column names must not contain \code{NA}.
#'
#' @importFrom cicero plot_connections
#' @import Gviz
#' @return A genome browser-style plot displaying the focused gene and predicted enhancer-gene connections.
#' @export
#'
#' @examples
#' #' data(multiomic_small)
#' SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#' fpath <- system.file("extdata", package = "SCEGHiC")
#' gene <- c("TRABD2A", "GNLY", "MFSD6", "CTLA4", "LCLAT1", "NCK2", "GALM", "TMSB10", "ID2", "CXCR4")
#' weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
#' results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
#' connections_Plot(results_SCEGHiC, species = "Homo sapiens", genome = "hg38", focus_gene = "CTLA4", cutoff = 0.01, gene_anno = NULL)
connections_Plot <- function(SCEG_HiC_Result, species, genome, focus_gene,
                             cutoff = NULL, upstream = 250000, downstream = 250000,
                             gene_anno = NULL) {
  # Check that only one gene is selected
  if (length(focus_gene) != 1) {
    stop("Please specify only one gene to visualize.")
  }

  # Filter results for the specified gene
  gene_hits <- which(SCEG_HiC_Result$gene %in% focus_gene)
  if (length(gene_hits) == 0) {
    message("No results found for the specified gene. Please try another gene.")
    return(NULL)
  }

  # Subset and format results
  result_gene <- SCEG_HiC_Result[gene_hits, , drop = FALSE]
  if ("peak2" %in% colnames(result_gene)) {
    result_gene <- result_gene[, c("gene", "peak2", "score")]
    colnames(result_gene) <- c("gene", "peak", "score")
  }

  # Annotate TSS
  tssdata <- SCEGHiC::annotateTSS(species, genome)
  gene_tss <- tssdata[tssdata$TargetGene == focus_gene, ]

  # Construct connections data frame
  conns <- data.frame(
    Peak1 = paste(gene_tss$chr, gene_tss$TargetGeneTSS, gene_tss$TargetGeneTSS + 1, sep = "_"),
    Peak2 = result_gene$peak,
    coaccess = result_gene$score
  )
  rownames(conns) <- NULL

  # Default cutoff
  if (is.null(cutoff)) {
    cutoff <- 0.01
  }

  # Generate plot
  plot_connections(
    conns,
    chr = gene_tss$chr,
    minbp = as.numeric(gene_tss$TargetGeneTSS) - upstream,
    maxbp = as.numeric(gene_tss$TargetGeneTSS) + downstream,
    gene_model = gene_anno,
    coaccess_cutoff = cutoff,
    connection_width = 0.5,
    connection_color = "orange",
    peak_color = "grey",
    collapseTranscripts = "longest"
  )
}

# Convert region argument to genomic coordinates
# Region can be a string, name of a gene, or GRanges object
#' @importFrom methods is
FindRegion <- function(object, region, sep = c("-", "-"), assay = NULL, extend.upstream = 0, extend.downstream = 0) {
  if (!is(object = region, class2 = "GRanges")) {
    # first try to convert to coordinates, if not lookup gene
    region <- tryCatch(
      expr = suppressWarnings(
        expr = StringToGRanges(regions = region, sep = sep)
      ),
      error = function(x) {
        region <- LookupGeneCoords(
          object = object,
          assay = assay,
          gene = region
        )
        return(region)
      }
    )
    if (is.null(x = region)) {
      stop("Gene not found")
    }
  }
  region <- suppressWarnings(expr = Extend(
    x = region,
    upstream = extend.upstream,
    downstream = extend.downstream
  ))
  return(region)
}

#' Plot linked genomic elements
#'
#' Visualize predicted links between genomic elements (e.g., genes and regulatory peaks) within a specified genomic region.
#'
#' @param links  A data.frame containing enhancer-gene links.
#' @param gene_anno A \code{data.frame} or \code{NULL}. If provided, should be compatible with the \code{\link[Gviz]{GeneRegionTrack-class}}.
#'   Column names must not contain \code{NA}.
#' @param gene A character string specifying the target gene to highlight.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object specifying the genomic region to plot.
#' @param lowcolor Color for the lowest scores in the gradient.
#' @param highcolor Color for the highest scores in the gradient.
#' @param titlename A character string specifying the title of the plot.
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot geom_hline aes theme_classic xlim xlab
#' ylab theme element_blank scale_color_gradient2 aes_string
#' @import ggforce
#' @return A \code{\link[ggplot2]{ggplot}} object visualizing enhancer-gene links.
#' @export
#'
#' @examples
#' links <- data.frame(
#'   gene = "CTLA4", peak = c("chr2_203623664_203623982", "chr2_203730093_203731243", "chr2_203992800_203993979"),
#'   score = c(0.03623216, 0.14814205, 0.43240254)
#' )
#' library(GenomicRanges)
#' region <- GRanges(seqnames = Rle("chr2"), ranges = IRanges(start = 203617771, end = 204117772))
#' geneanno <- annotateTSS("Homo sapiens", "hg38")
#' colnames(geneanno) <- c("chr", "gene", "tss")
#' links <- LinksPlot(links, geneanno, gene = "CTLA4", region, lowcolor = "blue", highcolor = "orange", titlename = "SCEG-HiC")
#' links
LinksPlot <- function(links, geneanno, gene, region, lowcolor, highcolor, titlename) {
  if (length(gene) != 1) {
    stop("You only choose one gene to display\n")
  }

  colnames(links) <- c("gene", "peak", "score")
  links.df <- merge(links, geneanno, by = "gene")
  links.df$peak <- gsub("_", "-", links.df$peak)
  links.df$midpoint <- apply(links.df, 1, function(x) {
    round((as.numeric(unlist(strsplit(x[2], "-"))[2]) + as.numeric(unlist(strsplit(x[2], "-"))[3])) / 2)
  })
  links.df$start <- pmin(links.df$tss, links.df$midpoint)
  links.df$end <- pmax(links.df$tss, links.df$midpoint)
  links.use <- links.df[links.df$gene == gene, ]
  links.use <- links.use[links.use$start >= start(x = region) & links.use$end <= end(x = region), ]
  links.use$group <- seq_len(length.out = nrow(x = links.use))
  df <- data.frame(
    x = c(
      links.use$start,
      (links.use$start + links.use$end) / 2,
      links.use$end
    ),
    y = c(
      rep(x = 0, nrow(x = links.use)),
      rep(x = -1, nrow(x = links.use)),
      rep(x = 0, nrow(x = links.use))
    ),
    group = rep(x = links.use$group, 3),
    score = rep(links.use$score, 3)
  )
  min.color <- min(0, min(df$score))
  p <- ggplot(data = df) +
    ggforce::geom_bezier(
      mapping = aes(x = .data$x, y = .data$y, group = .data$group, color = .data$score)
    )
  p <- p +
    geom_hline(yintercept = 0, color = "grey") +
    scale_color_gradient2(
      low = lowcolor, mid = "grey", high = highcolor,
      limits = c(min.color, max(df$score)),
      n.breaks = 3
    )
  p <- p +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    ylab(titlename) +
    xlab(label = paste0(seqnames(x = region), " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
}

#' Plot SCEG-HiC Links and Tn5 Insertion Frequency Over a Genomic Region
#'
#' This function visualizes Tn5 insertion signal across a given genomic window centered on a focus gene,
#' and overlays predicted and validated gene-enhancer (gene-peak) links.
#'
#' The output includes:
#' \itemize{
#'   \item{Tn5 insertion signal across cell groups}
#'   \item{Links predicted by SCEG-HiC}
#'   \item{Hi-C validated links}
#'   \item{Correlation-based links}
#'   \item{eQTL variants (optional)}
#'   \item{Peak positions}
#' }
#'
#' @param object A Seurat object.
#' @param focus_gene A character vector of gene symbols to focus on.
#' @param species Character string specifying the species name. Supported values are "Homo sapiens" or "Mus musculus".
#' @param genome Character string specifying the genome assembly. Supported values are "hg38", "hg19", "mm10", or "mm9".
#' @param assay Character or vector. Assay(s) to use. The first assay determines gene annotations and link metadata.
#' @param upstream Numeric specifying the number of base pairs upstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param downstream Numeric specifying the number of base pairs downstream of each TSS to define enhancers. Default is 250,000 bp (250 kb).
#' @param HIC_Result Data.frame. Hi-C gene-peak interaction result with columns `gene`, `peak`, `score`.
#' @param HIC_cutoff Numeric. Score threshold for Hi-C links. Default: 0.
#' @param SCEG_HiC_Result A data.frame containing the output from \code{Run_SCEG_HiC()}.
#' @param SCEG_HiC_cutoff Numeric. Score threshold for SCEG-HiC links. Recommend 0.01 for aggregated, 0.001 for single-cell.
#' @param correlation  Data.frame. Correlation-based gene-peak links.
#' @param cells Character vector. Subset of cells to include.
#' @param cellnames Character vector. Name(s) of one or more metadata columns used to group the cells.
#' Default is the current cell identities.
#' @param eqtl.positions Numeric vector. Genomic positions of eQTL variants (optional).
#'
#' @importFrom patchwork wrap_plots
#' @importFrom Signac PeakPlot CoveragePlot
#' @import Seurat
#' @import grid
#' @import ggplot2
#' @return A \code{\link[patchwork]{patchwork}} plot.
#' @export
#' @examples
#' data(multiomic_small)
#' SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#' fpath <- system.file("extdata", package = "SCEGHiC")
#' gene <- c("TRABD2A", "GNLY", "MFSD6", "CTLA4", "LCLAT1", "NCK2", "GALM", "TMSB10", "ID2", "CXCR4")
#' weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
#' results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
#' fpath <- system.file("extdata", "multiomic_small_atac_fragments.tsv.gz", package = "SCEGHiC")
#' library(Signac)
#' frags <- CreateFragmentObject(path = fpath, cells = colnames(multiomic_small))
#' Fragments(multiomic_small) <- frags
#' coverPlot(multiomic_small, focus_gene = "CTLA4", species = "Homo sapiens", genome = "hg38", assay = "peaks", SCEG_HiC_Result = results_SCEGHiC, SCEG_HiC_cutoff = 0.01)
coverPlot <- function(object, focus_gene, species, genome, assay = NULL, upstream = 250000, downstream = 250000, HIC_Result = NULL, HIC_cutoff = 0, SCEG_HiC_Result = NULL, SCEG_HiC_cutoff = 0, correlation = NULL, cells = NULL, cellnames = NULL, eqtl.positions = NULL) {
  # Fetch TSS data and locate gene
  tssdata <- SCEGHiC::annotateTSS(species, genome)
  colnames(tssdata) <- c("chr", "gene", "tss")
  gene.chr <- as.character(tssdata[tssdata$gene == focus_gene, "chr"])
  gene.tss <- as.numeric(tssdata[tssdata$gene == focus_gene, "tss"])

  # Define genomic region
  region <- FindRegion(
    object = object,
    region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1),
    extend.upstream = upstream,
    extend.downstream = downstream
  )

  # Hi-C links
  HIC_link <- if (!is.null(HIC_Result)) {
    HIC_Result <- HIC_Result[HIC_Result$score > HIC_cutoff, ]
    LinksPlot(HIC_Result, tssdata, focus_gene, region, "blue", "black", "Hi-C")
  } else {
    NULL
  }

  # SCEG-HiC links
  SCEG_HiC_link <- if (!is.null(SCEG_HiC_Result)) {
    if ("peak2" %in% colnames(SCEG_HiC_Result)) {
      SCEG_HiC_Result <- SCEG_HiC_Result[, c("gene", "peak2", "score")]
      colnames(SCEG_HiC_Result) <- c("gene", "peak", "score")
    }
    SCEG_HiC_Result <- SCEG_HiC_Result[SCEG_HiC_Result$score > SCEG_HiC_cutoff, ]
    LinksPlot(SCEG_HiC_Result, tssdata, focus_gene, region, "blue", "orange", "SCEG-HiC")
  } else {
    NULL
  }

  # Correlation links
  correlation_link <- if (!is.null(correlation)) {
    LinksPlot(correlation, tssdata, focus_gene, region, "blue", "red", "correlation")
  } else {
    NULL
  }

  # Peak positions
  p1 <- PeakPlot(
    object,
    region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1),
    assay = assay,
    color = "red",
    extend.upstream = upstream,
    extend.downstream = downstream
  )

  # Tn5 coverage signal
  p2 <- CoveragePlot(
    object = object,
    assay = assay,
    region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1),
    extend.upstream = upstream,
    extend.downstream = downstream,
    peaks = FALSE,
    links = FALSE,
    idents = cells,
    group.by = cellnames
  )

  obj.groups <- if (is.null(cellnames)) {
    as.character(Idents(object))
  } else {
    as.character(object@meta.data[[cellnames]])
  }

  bulk.height <- if (is.null(cells)) {
    nident <- length(unique(obj.groups))
    (1 / nident) * 10
  } else {
    min(2 * length(cells), 10)
  }

  # eQTL layer
  p3 <- if (!is.null(eqtl.positions)) {
    eqtl.df <- data.frame(pos = eqtl.positions)
    start.pos <- start(region)
    end.pos <- end(region)
    eqtl.df$pos <- pmax(pmin(eqtl.df$pos, end.pos), start.pos)

    ggplot(eqtl.df) +
      geom_point(aes(x = pos, y = 0), color = "black", shape = 17, size = 3) +
      theme_classic() +
      xlab(paste0(gene.chr, " position (bp)")) +
      xlim(c(start(region), end(region))) +
      ylab("eQTL") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    NULL
  }

  # Combine plots
  heights <- c(1.5, bulk.height, 1.5, 1, 1, 1, 1.5)
  plot.list <- list(HIC_link +
    theme(legend.position = "none"), p2[[1]][[1]], correlation_link, p2[[1]][[2]], p3, p1, SCEG_HiC_link)

  final_plot <- CombineTracks(
    plotlist = plot.list,
    heights = heights
  ) & theme(
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

  return(final_plot)
}
