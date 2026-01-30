
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCEG-HiC (Predicting enhancer-gene links from single-cell multi-omics data by integrating prior Hi-C information)

<!-- badges: start -->

<!-- badges: end -->

## Overview

SCEG-HiC predicts enhancer–gene links by integrating multi-omics
single-cell data (either paired scATAC-seq/RNA-seq or scATAC-seq alone)
with three-dimensional chromatin conformation information derived from
bulk average Hi-C data. It employs the weighted graphical lasso
(wglasso) model to incorporate average bulk Hi-C data, effectively
regularizing the correlation matrix with the prior Hi-C contact matrix
as a penalty term.

## Installation

### Required software

SCEG-HiC runs in the [R statistical computing
environment](https://www.r-project.org/). It requires R version 4.1.0 or
higher, [Bioconductor](https://bioconductor.org/) version 3.14 or
higher, and Seurat 4.0 or higher to access the latest features.

To install Bioconductor, open an R session and run:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
```

Next, install a few Bioconductor packages that are not installed
automatically:

``` r
BiocManager::install(c(
  'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
  'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
  'SummarizedExperiment', 'batchelor', 'HDF5Array',
  'terra', 'ggrastr', 'Gviz', 'rtracklayer', 'GenomeInfoDb', 'GenomicRanges'
))
```

Installation of other dependencies

- Install the Signac pacakge:
  `devtools::install_github("timoast/signac", ref = "develop")`. If you
  encounter any issues, please check the [Signac
  documentation](https://stuartlab.org/signac/).
- Install the Cicero package:
  `devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")`.
  If you encounter any issues, please check the [Cicero installation
  guide](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/#installing-cicero).

Now, you can install the development version of SCEG-HiC from
[GitHub](https://github.com/) with:

``` r
# If you haven't installed devtools yet, uncomment and run:
# install.packages("devtools")

# Install the development version of SCEG-HiC from GitHub
devtools::install_github("wuwei77lx/SCEGHiC")
```

If you prefer the stable release version from CRAN, run:

``` r
# Install the released version from CRAN
install.packages("SCEGHiC")
```

### Testing the installation

To verify that SCEG-HiC installed correctly, start a new R session and
run:

``` r
library(SCEGHiC)
```

If no errors appear, the package is successfully loaded and ready to
use.

## Quickstart

This basic example demonstrates how to analyze paired scATAC-seq/RNA-seq
data using SCEG-HiC:

``` r
library(SCEGHiC)
library(Signac)

# Load example multi-omics dataset
data(multiomic_small)

# Preprocess the data (aggregation)
SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)
#> Generating aggregated data
#> Aggregating cluster 0
#> Sample cells randomly.
#> There are 11 samples
#> Aggregating cluster 1
#> Sample cells randomly.
#> There are 11 samples

# Define genes of interest
gene <- c("TRABD2A", "GNLY", "MFSD6", "CTLA4", "LCLAT1", "NCK2", "GALM", "TMSB10", "ID2", "CXCR4")

# Get path to example average Hi-C data
fpath <- system.file("extdata", package = "SCEGHiC")

# Calculate Hi-C based weights for enhancer-gene pairs
weight <- calculateHiCWeights(SCEGdata, species = "Homo sapiens", genome = "hg38", focus_gene = gene, averHicPath = fpath)
#> Processing chromosome chr2...
#> Found 10 TSS loci on chr2.
#> Calculating Hi-C weights for gene TRABD2A...
#> Calculating Hi-C weights for gene GNLY...
#> Calculating Hi-C weights for gene MFSD6...
#> Calculating Hi-C weights for gene CXCR4...
#> Calculating Hi-C weights for gene CTLA4...
#> Calculating Hi-C weights for gene LCLAT1...
#> Calculating Hi-C weights for gene NCK2...
#> Calculating Hi-C weights for gene ID2...
#> Calculating Hi-C weights for gene GALM...
#> Calculating Hi-C weights for gene TMSB10...
#> Finished calculating Hi-C weights for all genes.

# Run the SCEG-HiC model
results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene)
#> Total predicted genes: 10
#> Running model for gene: TRABD2A
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.43"
#> Running model for gene: GNLY
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.19"
#> Running model for gene: MFSD6
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.22"
#> Running model for gene: CXCR4
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.14"
#> Running model for gene: CTLA4
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.17"
#> Running model for gene: LCLAT1
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.41"
#> Running model for gene: NCK2
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.25"
#> Running model for gene: ID2
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.13"
#> Running model for gene: GALM
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.11"
#> Running model for gene: TMSB10
#> [1] "The optimal penalty parameter (rho) selected by BIC is: 0.44"

# Arc plot visualization predicted enhancer-gene links for CTLA4
connections_Plot(results_SCEGHiC, species = "Homo sapiens", genome = "hg38", focus_gene = "CTLA4", cutoff = 0.01, gene_anno = NULL)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" style="display: block; margin: auto;" />

``` r
# Load fragment data for coverage plotting
frag_path <- system.file("extdata", "multiomic_small_atac_fragments.tsv.gz", package = "SCEGHiC")
frags <- CreateFragmentObject(path = frag_path, cells = colnames(multiomic_small))
#> Computing hash
Fragments(multiomic_small) <- frags

# Coverage plot and visualize the links of CTLA4
coverPlot(multiomic_small, focus_gene = "CTLA4", species = "Homo sapiens", genome = "hg38",
          assay = "peaks", SCEG_HiC_Result = results_SCEGHiC, SCEG_HiC_cutoff = 0.01)
#> Warning in .merge_two_Seqinfo_objects(x, y): The 2 combined objects have no sequence levels in common. (Use
#>   suppressWarnings() to suppress this warning.)
```

<img src="man/figures/README-unnamed-chunk-3-2.png" width="100%" style="display: block; margin: auto;" />

<details>

<summary>

**Session Info**
</summary>

``` r
sessionInfo()
#> R version 4.4.2 (2024-10-31)
#> Platform: x86_64-conda-linux-gnu
#> Running under: Rocky Linux 9.6 (Blue Onyx)
#> 
#> Matrix products: default
#> BLAS/LAPACK: /home/liangxuan/conda/envs/test/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Asia/Shanghai
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] harmony_1.2.3                     Rcpp_1.0.14                      
#>  [3] Matrix_1.6-5                      dplyr_1.1.4                      
#>  [5] EnsDb.Mmusculus.v79_2.99.0        SCEGHiC_0.0.0.9000               
#>  [7] SeuratDisk_0.0.0.9021             hdf5r_1.3.12                     
#>  [9] ggplot2_3.5.1                     BSgenome.Hsapiens.UCSC.hg38_1.4.5
#> [11] BSgenome_1.74.0                   rtracklayer_1.66.0               
#> [13] BiocIO_1.16.0                     Biostrings_2.74.1                
#> [15] XVector_0.46.0                    EnsDb.Hsapiens.v86_2.99.0        
#> [17] ensembldb_2.30.0                  AnnotationFilter_1.30.0          
#> [19] GenomicFeatures_1.58.0            AnnotationDbi_1.68.0             
#> [21] Biobase_2.66.0                    GenomicRanges_1.58.0             
#> [23] GenomeInfoDb_1.42.1               IRanges_2.40.1                   
#> [25] S4Vectors_0.44.0                  BiocGenerics_0.52.0              
#> [27] Seurat_5.2.0                      SeuratObject_5.0.2               
#> [29] sp_2.1-4                          Signac_1.14.9001                 
#> 
#> loaded via a namespace (and not attached):
#>   [1] R.methodsS3_1.8.2           dichromat_2.0-0.1          
#>   [3] progress_1.2.3              urlchecker_1.0.1           
#>   [5] nnet_7.3-20                 goftest_1.2-3              
#>   [7] vctrs_0.6.5                 spatstat.random_3.3-2      
#>   [9] digest_0.6.37               png_0.1-8                  
#>  [11] ggrepel_0.9.6               deldir_2.0-4               
#>  [13] parallelly_1.41.0           pkgdown_2.1.1              
#>  [15] MASS_7.3-64                 reshape2_1.4.4             
#>  [17] httpuv_1.6.15               withr_3.0.2                
#>  [19] ggrastr_1.0.2               xfun_0.50                  
#>  [21] ellipsis_0.3.2              survival_3.8-3             
#>  [23] commonmark_1.9.2            memoise_2.0.1              
#>  [25] rcmdcheck_1.4.0             ggbeeswarm_0.7.2           
#>  [27] profvis_0.4.0               systemfonts_1.1.0          
#>  [29] ragg_1.3.3                  zoo_1.8-12                 
#>  [31] pbapply_1.7-2               R.oo_1.27.0                
#>  [33] Formula_1.2-5               prettyunits_1.2.0          
#>  [35] KEGGREST_1.46.0             promises_1.3.2             
#>  [37] httr_1.4.7                  restfulr_0.0.15            
#>  [39] globals_0.16.3              fitdistrplus_1.2-2         
#>  [41] ps_1.8.1                    rstudioapi_0.17.1          
#>  [43] UCSC.utils_1.2.0            miniUI_0.1.1.1             
#>  [45] generics_0.1.3              processx_3.8.5             
#>  [47] base64enc_0.1-3             curl_6.0.1                 
#>  [49] zlibbioc_1.52.0             polyclip_1.10-7            
#>  [51] xopen_1.0.1                 GenomeInfoDbData_1.2.13    
#>  [53] SparseArray_1.6.0           xtable_1.8-4               
#>  [55] stringr_1.5.1               desc_1.4.3                 
#>  [57] evaluate_1.0.3              S4Arrays_1.6.0             
#>  [59] BiocFileCache_2.14.0        hms_1.1.3                  
#>  [61] irlba_2.3.5.1               colorspace_2.1-1           
#>  [63] filelock_1.0.3              ROCR_1.0-11                
#>  [65] reticulate_1.40.0           spatstat.data_3.1-4        
#>  [67] magrittr_2.0.3              lmtest_0.9-40              
#>  [69] later_1.4.1                 lattice_0.22-6             
#>  [71] spatstat.geom_3.3-4         future.apply_1.11.3        
#>  [73] scattermore_1.2             XML_3.99-0.17              
#>  [75] cowplot_1.1.3               matrixStats_1.5.0          
#>  [77] RcppAnnoy_0.0.22            Hmisc_5.2-2                
#>  [79] pillar_1.10.1               nlme_3.1-166               
#>  [81] cicero_1.3.9                compiler_4.4.2             
#>  [83] RSpectra_0.16-2             stringi_1.8.7              
#>  [85] tensor_1.5                  minqa_1.2.8                
#>  [87] SummarizedExperiment_1.36.0 devtools_2.4.5             
#>  [89] GenomicAlignments_1.42.0    plyr_1.8.9                 
#>  [91] crayon_1.5.3                abind_1.4-8                
#>  [93] bit_4.5.0.1                 fastmatch_1.1-6            
#>  [95] NCmisc_1.2.0                codetools_0.2-20           
#>  [97] textshaping_1.0.1           monocle3_1.3.7             
#>  [99] bslib_0.8.0                 slam_0.1-55                
#> [101] biovizBase_1.54.0           plotly_4.10.4              
#> [103] mime_0.12                   splines_4.4.2              
#> [105] fastDummies_1.7.4           dbplyr_2.5.0               
#> [107] interp_1.1-6                knitr_1.49                 
#> [109] blob_1.2.4                  lme4_1.1-36                
#> [111] fs_1.6.5                    listenv_0.9.1              
#> [113] checkmate_2.3.2             Rdpack_2.6.2               
#> [115] pkgbuild_1.4.5              Gviz_1.50.0                
#> [117] tibble_3.2.1                callr_3.7.6                
#> [119] tweenr_2.0.3                pkgconfig_2.0.3            
#> [121] tools_4.4.2                 cachem_1.1.0               
#> [123] RhpcBLASctl_0.23-42         rbibutils_2.3              
#> [125] RSQLite_2.3.9               viridisLite_0.4.2          
#> [127] DBI_1.2.3                   fastmap_1.2.0              
#> [129] rmarkdown_2.29              scales_1.3.0               
#> [131] grid_4.4.2                  usethis_3.1.0              
#> [133] ica_1.0-3                   Rsamtools_2.22.0           
#> [135] sass_0.4.9                  patchwork_1.3.0            
#> [137] FNN_1.1.4.1                 dotCall64_1.2              
#> [139] VariantAnnotation_1.52.0    RANN_2.6.2                 
#> [141] rpart_4.1.24                farver_2.1.2               
#> [143] reformulas_0.4.0            yaml_2.3.10                
#> [145] VGAM_1.1-12                 roxygen2_7.3.3             
#> [147] latticeExtra_0.6-30         MatrixGenerics_1.18.1      
#> [149] foreign_0.8-88              cli_3.6.3                  
#> [151] purrr_1.0.2                 lifecycle_1.0.4            
#> [153] uwot_0.2.2                  sessioninfo_1.2.2          
#> [155] backports_1.5.0             BiocParallel_1.40.0        
#> [157] gtable_0.3.6                rjson_0.2.23               
#> [159] ggridges_0.5.6              progressr_0.15.1           
#> [161] testthat_3.2.3              parallel_4.4.2             
#> [163] jsonlite_1.8.9              RcppHNSW_0.6.0             
#> [165] bitops_1.0-9                bit64_4.5.2                
#> [167] assertthat_0.2.1            brio_1.1.5                 
#> [169] Rtsne_0.17                  glasso_1.11                
#> [171] spatstat.utils_3.1-2        jquerylib_0.1.4            
#> [173] spatstat.univar_3.1-1       R.utils_2.12.3             
#> [175] lazyeval_0.2.2              shiny_1.10.0               
#> [177] htmltools_0.5.8.1           sctransform_0.4.1          
#> [179] rappdirs_0.3.3              glue_1.8.0                 
#> [181] spam_2.11-0                 httr2_1.0.7                
#> [183] RCurl_1.98-1.16             rprojroot_2.0.4            
#> [185] jpeg_0.1-10                 gridExtra_2.3              
#> [187] boot_1.3-31                 igraph_2.0.3               
#> [189] R6_2.5.1                    tidyr_1.3.1                
#> [191] SingleCellExperiment_1.28.1 RcppRoll_0.3.1             
#> [193] labeling_0.4.3              cluster_2.1.8              
#> [195] pkgload_1.4.0               nloptr_2.1.1               
#> [197] DelayedArray_0.32.0         tidyselect_1.2.1           
#> [199] vipor_0.4.7                 ProtGenerics_1.38.0        
#> [201] htmlTable_2.4.3             ggforce_0.4.2              
#> [203] xml2_1.5.0                  vdiffr_1.0.8               
#> [205] future_1.34.0               munsell_0.5.1              
#> [207] KernSmooth_2.23-26          data.table_1.16.4          
#> [209] htmlwidgets_1.6.4           RColorBrewer_1.1-3         
#> [211] biomaRt_2.62.0              rlang_1.1.4                
#> [213] spatstat.sparse_3.1-0       spatstat.explore_3.3-4     
#> [215] remotes_2.5.0               reader_1.0.6               
#> [217] beeswarm_0.4.0
```

</details>

See the [documentation website](https://wuwei77lx.github.io/SCEGHiC/)
for more information!

## The bulk average Hi-C data

The human cell types used for averaging are based on [34 Hi-C
datasets](https://www.encodeproject.org/annotations/ENCSR382HAW/) from
the ENCODE project.

The mouse cell types used for averaging are: two embryonic stem cell
types (mESC1, mESC2), CH12LX, CH12F3, fiber, epithelium, and B cells.

The bulk average Hi-C data can be generated using the [Activity by
Contact (ABC)](https://www.nature.com/articles/s41588-019-0538-0)
model’s `makeAverageHiC.py` script.

### Download Links for Average Hi-C Data

- **Human bulk average Hi-C (~54.1 GB)**

  Download from:
  [ENCFF134PUN.bed.gz](https://www.encodeproject.org/files/ENCFF134PUN/@@download/ENCFF134PUN.bed.gz).

After downloading, extract the human bulk average Hi-C using the
[Activity by Contact
(ABC)](https://www.nature.com/articles/s41588-019-0538-0) model’s
`extract_avg_hic.py` script:

``` sh
python code/Hi_C/extract_avg_hic.py --avg_hic_bed_file ../ENCFF134PUN.bed.gz --output_dir ../
```

- **Mouse bulk average Hi-C (~13.9 GB)**

  Download from:
  [10.5281/zenodo.14849886](https://zenodo.org/record/14849886).

For more details about bulk average Hi-C data, please visit:
<https://github.com/wuwei77lx/compare_model>.

## Example

In SCEG-HiC, you can choose either aggregation or single-cell retention
approach:

- **Aggregation approach**: Aggregates binarized scATAC-seq data across
  cell types using k-nearest neighbor smoothing to reduce data sparsity.
  This approach captures a broader spectrum of enhancer-gene links
  across cell types, with slightly reduced prediction accuracy.

- **Single-cell retention approach**: Normalizes scATAC-seq data within
  individual cell types to individual cell signals. This method achieves
  higher precision and accuracy, albeit identifying fewer enhancer-gene
  links

**Recommendation**: To balance accuracy and coverage, we implemented
both preprocessing strategies in SCEG-HiC, with aggregation designated
as the default. The single-cell retention approach can be optionally
used when higher precision within specific cell types is desired.

For more details and real data examples, please visit:

- [SCEG-HiC on paired scATAC-seq/RNA-seq data of PBMC
  (aggregation)](https://wuwei77lx.github.io/SCEGHiC/articles/PBMC_multiomic_aggregation.html)

- [SCEG-HiC on paired scATAC-seq/RNA-seq data of PBMC
  CD4T](https://wuwei77lx.github.io/SCEGHiC/articles/PBMC_multiomic_CD4T.html)

- [SCEG-HiC on paired scATAC-seq/RNA-seq data of mouse skin
  (aggregation)](https://wuwei77lx.github.io/SCEGHiC/articles/mouse_skin_multiomic_aggregation.html)

- [SCEG-HiC on scATAC-seq data alone from human COVID-19
  monocytes](https://wuwei77lx.github.io/SCEGHiC/articles/human_covid_19_scATAC_seq_monocytes.html)

Alternatively, when applying SCEG-HiC to a new tissue, the penalty
strength can be adjusted using the alpha parameter (`alpha * rho`):

- **Increase alpha** to prioritize high-confidence enhancer-gene links,
  e.g., when the new tissue is developmentally similar to the datasets
  used to construct the bulk averaged Hi-C prior or when the single-cell
  data are noisy.

- **Decrease alpha** to explore a broader set of potential links, e.g.,
  for rare, unique, or disease-specific tissues.

``` r
# Example: run SCEG-HiC with scaled penalty
results_SCEGHiC <- Run_SCEG_HiC(SCEGdata, weight, focus_gene = gene, alpha=1.5)
```

**Notes:**

- `alpha = 1` uses the original penalty (default).

- `alpha > 1` strengthens the penalty.

- `alpha < 1`weakens the penalty.

## Help

If you have any questions, comments, or suggestions, please contact Xuan
Liang at <liangxuan2022@sinh.ac.cn>.
