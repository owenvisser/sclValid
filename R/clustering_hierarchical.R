#' Run Hierarchical Clustering on Reduced Distances
#'
#' Internal helper that performs Ward hierarchical clustering using a
#' precomputed distance object and returns cluster assignments for one or more
#' requested numbers of clusters.
#'
#' @param d A distance object containing pairwise distances between
#'   observations.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   to obtain from the hierarchical clustering tree.
#' @param prefix Character string used to identify the dimensionality
#'   reduction in the returned result names, such as `"PCA"` or `"TSNE"`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return A named list of cluster-assignment vectors. For a cluster size `k`,
#'   result names follow the form `"{prefix}_HC_cs{k}"`.
#'
#' @keywords internal
.run_hierarchical <- function(d, cluster_sizes, prefix, verbose = TRUE) {

  results <- list()

  tree <- stats::hclust(
    d = d,
    method = "ward.D"
  )

  for (k in cluster_sizes) {

    if (verbose) {
      message(paste0(prefix, " hierarchical: clustering with k = ", k))
    }

    results[[paste0(prefix, "_HC_cs", k)]] <-
      stats::cutree(
        tree = tree,
        k = k
      )
  }

  return(results)
}