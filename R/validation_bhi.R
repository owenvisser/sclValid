#' Calculate Biological Homogeneity Index
#'
#' Calculates the Biological Homogeneity Index (BHI), which measures
#' homogeneity of known biological labels within clusters.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param labels A vector containing the known label for each observation.
#'
#' @return A numeric biological homogeneity index between zero and one, with
#'   larger values indicating greater biological homogeneity.
#'
#' @references
#' Pihur V, Brock GN, Datta S.
#' Cluster validation for microarray data: an appraisal.
#' In: Advances in Multivariate Statistical Methods.
#' World Scientific; 2009:79-94.
#'
#' @examples
#' clusters <- c(1, 1, 2, 2)
#' labels <- c("A", "A", "B", "B")
#' validation_bhi(clusters, labels)
#'
#' @export
validation_bhi <- function(clusters, labels) {

  unique_clusters <- unique(clusters)
  score <- 0

  for (cluster in unique_clusters) {

    cluster_index <- which(clusters == cluster)
    n_cluster <- length(cluster_index)

    matches <- 0

    if (n_cluster > 1) {

      for (j in 1:(n_cluster - 1)) {
        for (h in (j + 1):n_cluster) {

          if (labels[cluster_index[j]] == labels[cluster_index[h]]) {
            matches <- matches + 1
          }
        }
      }

      score <- score + matches / (n_cluster * (n_cluster - 1) / 2)
    }
  }

  score / length(unique_clusters)
}