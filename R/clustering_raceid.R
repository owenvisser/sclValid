#' Run RaceID Clustering
#'
#' Applies RaceID clustering to the data assay of a `SingleCellExperiment`
#' object for one or more requested numbers of clusters. Cluster assignments
#' are added to the object's `colData`.
#'
#' If feature or observation names are absent from the input matrix, temporary
#' names are generated internally because RaceID requires them.
#'
#' For a cluster size `k`, the resulting column is named `"RaceID_cs{k}"`.
#'
#' @param sce A `SingleCellExperiment` object containing a `"data"` assay.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   to fit.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return The input `SingleCellExperiment` object with additional RaceID
#'   cluster-assignment columns in `colData`.
#'
#' @references
#' Grün D, Lyubimova A, Kester L, et al.
#' Single-cell messenger RNA sequencing reveals rare intestinal cell types.
#' Nature. 2015;525:251–255.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce <- run_raceid(sce, cluster_sizes = 2:4)
#' SummarizedExperiment::colData(sce)
#'
#' @export
run_raceid <- function(sce, cluster_sizes, verbose = TRUE) {

  d <- SummarizedExperiment::assay(sce, "data")

  if (is.null(rownames(d))) {
    rownames(d) <- paste0("gene_", seq_len(nrow(d)))
  }

  if (is.null(colnames(d))) {
    colnames(d) <- paste0("cell_", seq_len(ncol(d)))
  }

  scm <- RaceID::SCseq(d)

  scm <- RaceID::filterdata(
    scm,
    mintotal = 1,
    minexpr = 1,
    minnumber = 1,
    bmode = "RaceID",
    verbose = FALSE
  )

  scm <- RaceID::compdist(
    scm,
    metric = "pearson",
    FSelect = FALSE
  )

  for (k in cluster_sizes) {

    if (verbose) {
      message("RaceID: clustering with k = ", k)
    }
    
    fit <- RaceID::clustexp(
      scm,
      sat = FALSE,
      cln = k,
      clustnr = 30,
      bootnr = 50,
      rseed = 17000,
      FUNcluster = "kmedoids",
      verbose = FALSE
    )

    sce[[paste0("RaceID_cs", k)]] <- fit@cluster$kpart
  }

  return(sce)
}