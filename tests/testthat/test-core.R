test_that("process_data works", {
  SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
  expect_equal(class(SCEGdata), "list")
  expect_true("rna" %in% names(SCEGdata))
  expect_true("atac" %in% names(SCEGdata))
  expect_error(process_data(multiomic_small, celltype = "0"))
  expect_error(process_data(multiomic_small, aggregate = F))
})

test_that("calculateHiCWeights works", {
  SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
  gene <- multiomic_small@assays[["RNA"]]@data@Dimnames[[1]]
  fpath <- system.file("extdata", package = "SCEGHiC")
  weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
  expect_equal(length(weight), 85)
  expect_error(calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = "GLB1", averHicPath = fpath))
})

test_that("Run_SCEG_HiC works", {
  SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
  gene <- multiomic_small@assays[["RNA"]]@data@Dimnames[[1]]
  fpath <- system.file("extdata", package = "SCEGHiC")
  weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
  results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
  expect_equal(class(results_SCEGHiC), "data.frame")
})


test_that("connections_Plot works", {
  SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
  gene <- multiomic_small@assays[["RNA"]]@data@Dimnames[[1]]
  fpath <- system.file("extdata", package = "SCEGHiC")
  weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
  results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
  vdiffr::expect_doppelganger(
    "connections_Plot_CTLA4",
    connections_Plot(results_SCEGHiC,
      species = "Homo sapiens", genome = "hg38",
      focus_gene = "CTLA4", cutoff = 0.01, gene_anno = NULL
    )
  )
})

test_that("coverPlot works", {
  library(Signac)
  SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
  gene <- multiomic_small@assays[["RNA"]]@data@Dimnames[[1]]
  fpath <- system.file("extdata", package = "SCEGHiC")
  weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
  results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
  expect_error(coverPlot(multiomic_small,
    focus_gene = "CTLA4", species = "Homo sapiens", genome = "hg38",
    assay = "peaks", SCEG_HiC_Result = results_SCEGHiC, SCEG_HiC_cutoff = 0.01
  ))
  fpath <- system.file("extdata", "multiomic_small_atac_fragments.tsv.gz", package = "SCEGHiC")
  frags <- CreateFragmentObject(path = fpath, cells = colnames(multiomic_small))
  Fragments(multiomic_small) <- frags
  vdiffr::expect_doppelganger(
    "coverPlot_CTLA4",
    coverPlot(multiomic_small,
      focus_gene = "CTLA4", species = "Homo sapiens", genome = "hg38",
      assay = "peaks", SCEG_HiC_Result = results_SCEGHiC, SCEG_HiC_cutoff = 0.01
    )
  )
})
