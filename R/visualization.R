#' Plot links of SCEG-HiC on focused genes
#'
#' @param SCEG_HiC_Result SCEG_HiC result
#' @param species String indicating the name of the species.It can be either "Homo sapiens","Mus musculus".
#' @param genome String indicating the name of the genome assembly.It can be either "hg38",
#' "hg19", "mm10", "mm9".
#' @param focus_gene The focused a gene.
#' @param cutoff The threshold of gene-peak.Default is NULL.If \code{"aggregate"} is TRUE ,we recommend the cutoff of 0.01.
#' If \code{"aggregate"} is FALSE ,we recommend the cutoff of 0.001.
#' @param upstream Numeric specifying the window size (in base pairs) to pad around upstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param downstream Numeric specifying the window size (in base pairs) to pad around downstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param gene_anno Either \code{NULL} or a data.frame. The data.frame should be in a form compatible with the Gviz function.
#'   \code{\link[Gviz]{GeneRegionTrack-class}} (cannot have NA as column names).
#'
#' @importFrom cicero plot_connections
#' @import Gviz
#' @return A gene region and gene-peak pairs plot
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
connections_Plot <- function(SCEG_HiC_Result, species, genome, focus_gene, cutoff = NULL, upstream = 250000, downstream = 250000, gene_anno = NULL) {
  if (length(focus_gene) != 1) {
    stop("You only choose one gene to display\n")
  }

  # take out the result of marker
  if (length(which(SCEG_HiC_Result$gene %in% focus_gene)) > 0) {
    SCEG_HiC_Result <- SCEG_HiC_Result[which(SCEG_HiC_Result$gene %in% focus_gene), , drop = FALSE]
    if ("peak2" %in% colnames(SCEG_HiC_Result)) {
      SCEG_HiC_Result <- SCEG_HiC_Result[, c("gene", "peak2", "score")]
      colnames(SCEG_HiC_Result) <- c("gene", "peak", "score")
    }
    tssdata <- SCEGHiC::annotateTSS(species, genome)
    conns <- data.frame(
      Peak1 = paste(tssdata[tssdata$TargetGene == focus_gene, ]$chr, tssdata[tssdata$TargetGene == focus_gene, ]$TargetGeneTSS, tssdata[tssdata$TargetGene == focus_gene, ]$TargetGeneTSS + 1, sep = "_"),
      Peak2 = SCEG_HiC_Result$peak, coaccess = SCEG_HiC_Result$score
    )
    rownames(conns) <- NULL
    if (is.null(cutoff)) {
      cutoff <- 0.01
    }
    plot_connections(conns, tssdata[tssdata$TargetGene == focus_gene, ]$chr, as.numeric(tssdata[tssdata$TargetGene == focus_gene, ]$TargetGeneTSS) - upstream, as.numeric(tssdata[tssdata$TargetGene == focus_gene, ]$TargetGeneTSS) + downstream,
      gene_model = gene_anno,
      coaccess_cutoff = cutoff,
      connection_width = .5,
      connection_color = "orange",
      peak_color = "grey",
      collapseTranscripts = "longest"
    )
  } else {
    message("Please try other markers!")
  }
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
#' Display links between pairs of genomic elements within a given region of the genome.
#'
#' @param links  The result of gene-peak pairs.
#' @param geneanno Either \code{NULL} or a data.frame. The data.frame should be in a form compatible with the Gviz function.
#'   \code{\link[Gviz]{GeneRegionTrack-class}} (cannot have NA as column names).
#' @param gene The focused genes.
#' @param region A genomic region to plot.
#' @param lowcolor Colours for low ends of the gradient.
#' @param highcolor Colours for high ends of the gradient.
#' @param titlename The title of linking genomic elements.
#'
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges start end
#' @importFrom GenomeInfoDb seqnames
#' @importFrom ggplot2 ggplot geom_hline aes theme_classic xlim xlab
#' ylab theme element_blank scale_color_gradient2 aes_string
#' @import ggforce
#' @return Returns a \code{\link[ggplot2]{ggplot}} object
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
      mapping = aes_string(x = "x", y = "y", group = "group", color = "score")
    ) +
    theme(legend.position = "none")
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
      axis.text.y = element_blank(), legend.position = "none"
    ) +
    ylab(titlename) +
    xlab(label = paste0(seqnames(x = region), " position (bp)")) +
    xlim(c(start(x = region), end(x = region)))
}

#' Plot links of SCEG-HiC on focused genes and Tn5 insertion frequency over a region
#'
#' Plot frequency of Tn5 insertion events for different groups of cells within
#' given regions of the genome.  And plot links for gene-peak  predicted by
#' SCEG-HiC within given regions of the genome.
#'
#' Additional information can be layered on the coverage plot by setting several
#' different options in the CoveragePlot function. This includes showing:
#' \itemize{
#' \item{gene annotations}
#' \item{peak positions}
#' \item{gene-peak results predicted by correlation}
#' \item{gene-peak results validated by Hi-C}
#' }
#'
#' @param object A Seurat object.
#' @param focus_gene The focused a gene.
#' @param species String indicating the name of the species.It can be either "Homo sapiens","Mus musculus".
#' @param genome String indicating the name of the genome assembly.It can be either "hg38",
#' "hg19", "mm10", "mm9".
#' @param assay Name of the assay to plot. If a list of assays is provided,  data from each assay will be
#' shown overlaid on each track. The first assay in the list will define the assay used for gene annotations,
#' links, and peaks (if shown). The order of assays given defines the plotting order.
#' @param upstream Numeric specifying the window size (in base pairs) to pad around upstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param downstream Numeric specifying the window size (in base pairs) to pad around downstream of each TSS, to fetch gene's enhancers.
#' Default is 250 kb (so 500 kb window is drawn arund each TSS).
#' @param HIC_Result Hi-C result.
#' @param HIC_cutoff The threshold of gene-peak validated by Hi-C.Default is 0.
#' @param SCEG_HiC_Result SCEG_HiC result.
#' @param SCEG_HiC_cutoff The threshold of gene-peak predicted by SCEG-HiC.Default is 0.If \code{"aggregate"} is TRUE ,we recommend the cutoff of 0.01.
#' If \code{"aggregate"} is FALSE ,we recommend the cutoff of 0.001.
#' @param correlation  Correlation between elements result.
#' @param cells Which cells to plot. Default all cells.
#' @param cellnames Name of one or more metadata columns to group (color) the cells by. Default is the current cell identities.
#'
#' @importFrom patchwork wrap_plots
#' @importFrom Signac PeakPlot CoveragePlot
#' @import Seurat
#' @import grid
#' @import ggplot2
#' @return Returns a \code{\link[patchwork]{patchwork}} object
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
coverPlot <- function(object, focus_gene, species, genome, assay = NULL, upstream = 250000, downstream = 250000, HIC_Result = NULL, HIC_cutoff = 0, SCEG_HiC_Result = NULL, SCEG_HiC_cutoff = 0, correlation = NULL, cells = NULL, cellnames = NULL) {
  tssdata <- SCEGHiC::annotateTSS(species, genome)
  colnames(tssdata) <- c("chr", "gene", "tss")
  gene.chr <- as.character(tssdata[tssdata$gene == focus_gene, ]$chr)
  gene.tss <- as.numeric(tssdata[tssdata$gene == focus_gene, ]$tss)
  region <- FindRegion(
    object = object,
    region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1),
    extend.upstream = upstream,
    extend.downstream = downstream
  )
  if (!(is.null(HIC_Result))) {
    HIC_Result <- HIC_Result[HIC_Result$score > HIC_cutoff, ]
    HIC_link <- LinksPlot(HIC_Result, tssdata, focus_gene, region, "blue", "black", "Hi-C")
  } else {
    HIC_link <- NULL
  }
  if (!(is.null(SCEG_HiC_Result) & is.null(SCEG_HiC_cutoff))) {
    if ("peak2" %in% colnames(SCEG_HiC_Result)) {
      SCEG_HiC_Result <- SCEG_HiC_Result[, c("gene", "peak2", "score")]
      colnames(SCEG_HiC_Result) <- c("gene", "peak", "score")
    }
    SCEG_HiC_Result <- SCEG_HiC_Result[SCEG_HiC_Result$score > SCEG_HiC_cutoff, ]
    SCEG_HiC_link <- LinksPlot(SCEG_HiC_Result, tssdata, focus_gene, region, "blue", "orange", "SCEG-HiC")
  } else {
    SCEG_HiC_link <- NULL
  }
  if (!is.null(correlation)) {
    correlation_link <- LinksPlot(correlation, tssdata, focus_gene, region, "blue", "red", "correlation")
  } else {
    correlation_link <- NULL
  }
  p1 <- PeakPlot(object, region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1), assay = assay, color = "red", extend.upstream = upstream, extend.downstream = downstream)
  if (is.null(cellnames)) {
    p2 <- CoveragePlot(
      object = object,
      assay = assay,
      region = paste0(gene.chr, "-", gene.tss, "-", gene.tss + 1),
      extend.upstream = upstream,
      extend.downstream = downstream,
      peaks = FALSE,
      links = FALSE,
      idents = cells
    )
    obj.groups <- as.character(Idents(object))
  } else {
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
    obj.groups <- as.character(object@meta.data[[cellnames]])
  }

  if (is.null(cells)) {
    nident <- length(x = unique(x = obj.groups))
    bulk.height <- (1 / nident) * 10
  } else {
    bulk.height <- min(2 * length(cells), 10)
  }

  heights <- c(1.5, bulk.height, 1.5, 1, 1, 1.5)
  p <- CombineTracks(
    plotlist = list(HIC_link, p2[[1]][[1]], correlation_link, p2[[1]][[2]], p1, SCEG_HiC_link),
    heights = heights
  ) & theme(
    legend.key.size = unit(x = 1 / 2, units = "lines"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8)
  )
  return(p)
}
