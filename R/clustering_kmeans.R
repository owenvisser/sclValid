#' Run K-Means Clustering on Reduced Coordinates
#'
#' Internal helper that applies K-means clustering to a reduced-dimensional
#' representation for one or more requested numbers of clusters.
#'
#' @param x A numeric matrix of reduced-dimensional coordinates, with
#'   observations in rows.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   to fit.
#' @param prefix Character string used to identify the dimensionality
#'   reduction in the returned result names, such as `"PCA"` or `"TSNE"`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return A named list of cluster-assignment vectors. For a cluster size `k`,
#'   result names follow the form `"{prefix}_Kmeans_cs{k}"`.
#'
#' @keywords internal
.run_kmeans <- function(x, cluster_sizes, prefix, verbose = TRUE) {

  results <- list()

  for (k in cluster_sizes) {

    if (verbose) {
      message(paste0(prefix, " K-means: clustering with k = ", k))
    }

    fit <- stats::kmeans(
      x,
      centers = k,
      iter.max = 100,
      nstart = 10,
      trace = FALSE
    )

    results[[paste0(prefix, "_Kmeans_cs", k)]] <- fit$cluster
  }

  return(results)
}