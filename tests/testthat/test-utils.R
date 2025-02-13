test_that("liftOverPeaks works", {
  peaks <- c(
    "chr2-973888-974755", "chr2-3195868-3196264", "chr2-3242398-3243508", "chr2-6863856-6866834", "chr2-6913141-6914394",
    "chr2-6924663-6924972", "chr2-6999095-6999554", "chr2-7007713-7008819", "chr2-7083945-7084032", "chr2-7097148-7098029"
  )
  split_peaks <- do.call("rbind", strsplit(as.character(peaks), "-", fixed = TRUE))
  peakinfo <- liftOverPeaks(split_peaks, "hg38")
  liftover_peaks <- do.call("rbind", strsplit(peakinfo$peak, "_", fixed = TRUE))
  expect_equal(sum(peakinfo$start == liftover_peaks[, 2]), 0)
  expect_equal(sum(peakinfo$end == liftover_peaks[, 3]), 0)
  expect_equal(class(peakinfo$midpoint), "numeric")
})

test_that("annotateTSS  works", {
  hicweight <- data.frame(
    element1 = c(
      "CHN1", "CHN1", "CHN1", "CHN1", "CHN1", "chr2_174901074_174901475",
      "chr2_174901074_174901475", "chr2_174901074_174901475", "chr2_174901074_174901475",
      "chr2_174990102_174990619", "chr2_174990102_174990619", "chr2_174990102_174990619",
      "chr2_175048293_175048704", "chr2_175048293_175048704", "chr2_175167038_175169195"
    ),
    element2 = c(
      "chr2_174901074_174901475", "chr2_174990102_174990619", "chr2_175048293_175048704",
      "chr2_175167038_175169195", "chr2_175187467_175187717", "chr2_174990102_174990619",
      "chr2_175048293_175048704", "chr2_175167038_175169195", "chr2_175187467_175187717",
      "chr2_175048293_175048704", "chr2_175167038_175169195", "chr2_175187467_175187717",
      "chr2_175167038_175169195", "chr2_175187467_175187717", "chr2_175187467_175187717"
    ),
    score = c(
      0.0016062358, 0.0076253167, 0.0031021570, 0.0012148262, 0.0010788392, 0.0019968580,
      0.0012957845, 0.0009198567, 0.0011448870, 0.0023991885, 0.0010785481, 0.0009664508,
      0.0012865316, 0.0012122647, 0.0057781939
    )
  )
  element <- unique(c(hicweight$element1, hicweight$element2))
  weight <- normalizeHiCWeights(hicweight, "1", element)
  expect_equal(dim(weight)[1], dim(weight)[2])
  expect_equal(sum(diag(weight)), 6)
})

test_that("LinksPlot  works", {
  links <- data.frame(
    gene = "CTLA4", peak = c("chr2_203623664_203623982", "chr2_203730093_203731243", "chr2_203992800_203993979"),
    score = c(0.03623216, 0.14814205, 0.43240254)
  )
  library(GenomicRanges)
  region <- GRanges(seqnames = Rle("chr2"), ranges = IRanges(start = 203617771, end = 204117772))
  geneanno <- annotateTSS("Homo sapiens", "hg38")
  colnames(geneanno) <- c("chr", "gene", "tss")
  links <- LinksPlot(links, geneanno, gene = "CTLA4", region, lowcolor = "blue", highcolor = "orange", titlename = "SCEG-HiC")
  vdiffr::expect_doppelganger(
    "LinksPlot_SCEGHiC",
    links
  )
})
