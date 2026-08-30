#' Run a User-Defined Clustering Method
#'
#' Applies a user-supplied clustering function to the `"data"` assay of a
#' `SingleCellExperiment` object. The clustering function must return one
#' cluster assignment for each observation in the dataset.
#'
#' The clustering assignments are stored in `colData` using the name supplied
#' by `name`. Information required to rerun the clustering method is retained
#' internally so that stability-based validation measures can reproduce the
#' clustering on perturbed datasets.
#'
#' @param sce A `SingleCellExperiment` object containing a `"data"` assay.
#' @param fun A clustering function. The first argument supplied to `fun` is
#'   the numeric data matrix with features in rows and observations in columns.
#'   The function must return one cluster assignment per observation.
#' @param name Character string specifying the name of the clustering solution.
#' @param ... Additional arguments passed to `fun`.
#'
#' @return The input `SingleCellExperiment` object with the custom clustering
#'   assignments added to `colData`.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#'
#' my_kmeans <- function(data, centers) {
#'   stats::kmeans(t(data), centers = centers)$cluster
#' }
#'
#' sce <- run_custom_clustering(
#'   sce,
#'   fun = my_kmeans,
#'   name = "MyKmeans_cs3",
#'   centers = 3
#' )
#'
#' SummarizedExperiment::colData(sce)
#'
#' @export
run_custom_clustering <- function(
    sce,
    fun,
    name,
    ...
) {

  data <- SummarizedExperiment::assay(
    sce,
    "data"
  )

  assignments <- fun(
    data,
    ...
  )

  if (length(assignments) != ncol(sce)) {
    stop(
      "The custom clustering function must return one cluster assignment ",
      "per observation."
    )
  }

  sce[[name]] <- assignments

  metadata_sce <- S4Vectors::metadata(sce)

  if (is.null(metadata_sce$custom_clusterings)) {
    metadata_sce$custom_clusterings <- list()
  }

  metadata_sce$custom_clusterings[[name]] <- list(
    fun = fun,
    args = list(...)
  )

  S4Vectors::metadata(sce) <- metadata_sce

  sce
}