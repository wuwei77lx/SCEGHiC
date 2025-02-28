% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{process_data}
\alias{process_data}
\title{Create SCEG-HiC preprocessed input data}
\usage{
process_data(
  object,
  aggregate = TRUE,
  celltype = NULL,
  rna_assay = "RNA",
  atac_assay = "peaks",
  cellnames = NULL,
  atacbinary = TRUE,
  k_neigh = 50,
  max_overlap = 0.8,
  size_factor_normalize = TRUE,
  reduction.name = NULL,
  seed = 123,
  verbose = TRUE
)
}
\arguments{
\item{object}{Seurat object.}

\item{aggregate}{Logically,  should the data be clustered ?If FALSE, \code{celltype} must be specified as one cell type.}

\item{celltype}{Which cells to generate data.}

\item{rna_assay}{Name of assay containing gene expression information.If \code{aggregate} is TRUE, name of assay containing raw gene expression
information.If \code{aggregate} is FALSE, name of assay containing normalize gene expression information.}

\item{atac_assay}{Name of assay containing peak information.If \code{aggregate} is TRUE, name of assay containing raw peak information.
If \code{aggregate} is FALSE, name of assay containing normalize peak information.}

\item{cellnames}{Name of one or more metadata columns to group the cells by. Default is the current cell identities.}

\item{atacbinary}{Logical, whether the aggregated scATAC-seq data need binary?}

\item{k_neigh}{Number of cells to be aggregated per group.}

\item{max_overlap}{The maximum overlapping ratio of two groups.}

\item{size_factor_normalize}{Logical, should accessibility values be normalized by size factor?}

\item{reduction.name}{The reduction name of extracting the cell coordinates used for aggregating.}

\item{seed}{Random seed.}

\item{verbose}{Logical, should warning and info messages be printed?}
}
\value{
Preprocessed data for SCEG-HiC
}
\description{
Function to create preprocessed input data for SCEG-HiC.\code{calculateHiCWeights} and
\code{Run_SCEG_HiC}takes as input the generated data.This function can generate clustered data
or non clustered data for a specific cell type.
}
\examples{
data(multiomic_small)
# aggregated with paired scRNA-seq and scATAC-seq data
SCEGdata <- process_data(multiomic_small, k_neigh = 5, max_overlap = 0.5)

# select "1" cells with paired scRNA-seq and scATAC-seq data
SCEGdata <- process_data(multiomic_small, aggregate = FALSE, celltype = "1")

}
