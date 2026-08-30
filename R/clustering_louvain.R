#' Run Louvain Clustering on a Nearest-Neighbor Graph
#'
#' Internal helper that constructs mutual nearest-neighbor graphs from a
#' precomputed distance matrix and applies Louvain community detection across
#' combinations of neighborhood size and resolution.
#'
#' @param dmat A numeric matrix containing pairwise distances between
#'   observations.
#' @param knn An integer vector specifying the numbers of nearest neighbors
#'   used to construct the mutual nearest-neighbor graphs.
#' @param resolution A numeric vector specifying Louvain resolution values.
#' @param prefix Character string used to identify the dimensionality
#'   reduction in the returned result names, such as `"PCA"` or `"TSNE"`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return A named list of cluster-assignment vectors. Result names follow
#'   the form `"{prefix}_Louvain_knn{knn}_res{resolution}"`.
#'
#' @keywords internal
.run_louvain <- function(dmat, knn, resolution, prefix, verbose = TRUE) {

  results <- list()

  for (k_nn in knn) {

    graph <- cccd::nng(
      dx = dmat,
      mutual = TRUE,
      k = k_nn
    )

    for (res in resolution) {

      if (verbose) {
        message(
          paste0(
            prefix,
            " Louvain: clustering with knn = ",
            k_nn,
            ", resolution = ",
            res
          )
        )
      }

      fit <- igraph::cluster_louvain(
        graph,
        resolution = res
      )

      results[[
        paste0(
          prefix,
          "_Louvain_knn",
          k_nn,
          "_res",
          res
        )
      ]] <- fit$membership
    }
  }

  return(results)
}