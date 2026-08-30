#' Run CIDR Clustering
#'
#' Applies the CIDR clustering procedure to the data assay of a
#' `SingleCellExperiment` object for one or more requested numbers of clusters.
#' Cluster assignments are added to the object's `colData`.
#'
#' For a cluster size `k`, the resulting column is named `"CIDR_cs{k}"`.
#'
#' @param sce A `SingleCellExperiment` object containing a `"data"` assay.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   to fit.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return The input `SingleCellExperiment` object with additional CIDR
#'   cluster-assignment columns in `colData`.
#'
#' @references
#' Lin P, Troup M, Ho JWK. CIDR: Ultrafast and accurate clustering through
#' imputation for single-cell RNA-seq data.
#' Genome Biology. 2017;18:59.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce <- run_cidr(sce, cluster_sizes = 2:4)
#' SummarizedExperiment::colData(sce)
#'
#' @export
run_cidr <- function(sce, cluster_sizes, verbose = TRUE) {

  d <- SummarizedExperiment::assay(sce, "data")
  scd <- cidr::scDataConstructor(d)
  scd <- cidr::determineDropoutCandidates(scd)
  scd <- cidr::wThreshold(scd)
  scd <- cidr::scDissim(scd, threads = 1)
  scd <- cidr::scPCA(scd, plotPC = FALSE)
  scd <- cidr::nPC(scd)

  for (k in cluster_sizes) {

    if (verbose) {
      message(paste0("CIDR: clustering with k = ", k))
    }

    fit <- cidr::scCluster(
      scd,
      nCluster = k,
      nPC = scd@nPC,
      cMethod = "ward.D2"
    )

    sce[[paste0("CIDR_cs", k)]] <- fit@clusters
  }

  return(sce)
}