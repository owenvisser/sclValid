#' Calculate In-Group Proportion
#'
#' Calculates the In-Group Proportion (IGP), which summarizes the proportion
#' of observations whose nearest neighbor belongs to the same cluster.
#'
#' The implementation follows the overall IGP construction used by Pihur
#' and colleagues.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param nn A list containing nearest-neighbor indices for each observation.
#'
#' @return A numeric in-group proportion, with larger values indicating
#'   stronger agreement between clusters and local neighborhood structure.
#'
#' @references
#' Kapp AV, Tibshirani R. Are clusters found in one dataset present in another
#' dataset? Biostatistics. 2007;8(1):9-31.
#'
#' Pihur V, Brock GN, Datta S.
#' Cluster validation for microarray data: an appraisal.
#' In: Advances in Multivariate Statistical Methods.
#' World Scientific; 2009:79-94.
#'
#' @examples
#' clusters <- c(1, 1, 2, 2)
#' x <- matrix(
#'   c(
#'     0, 0, 3, 3,
#'     0, 1, 3, 4
#'   ),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#' nn <- list(c(2, 3),  c(1, 3),  c(4, 2),  c(3, 2))
#' validation_igp(clusters, nn)
#'
#' @export
validation_igp <- function(clusters, nn) {

  unique_clusters <- unique(clusters)

  score <- 0

  for (cluster in unique_clusters) {

    cluster_index <- which(clusters == cluster)

    nearest_in_cluster <- sum(
      vapply(
        cluster_index,
        function(i) nn[[i]][1] %in% cluster_index,
        logical(1)
      )
    )

    score <- score + nearest_in_cluster / length(cluster_index)
  }

  score / length(unique_clusters)
}