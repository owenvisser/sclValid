#' Calculate Adjusted Rand Index
#'
#' Calculates the adjusted Rand index between a clustering solution and known
#' observation labels. The adjusted Rand index measures agreement between two
#' partitions while correcting for agreement expected by chance.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param labels A vector containing the known label for each observation.
#'
#' @return A numeric adjusted Rand index.
#'
#' @references
#' Hubert L, Arabie P. Comparing partitions.
#' Journal of Classification. 1985;2:193-218.
#' 
#' @examples
#' clusters <- c(1, 1, 2, 2)
#' labels <- c("A", "A", "B", "B")
#' validation_ari(clusters, labels)
#'
#' @export
validation_ari <- function(clusters, labels) {

  mclust::adjustedRandIndex(
    clusters,
    labels
  )
}