#' Calculate the Dunn Index
#'
#' Calculates the Dunn Index (DI), defined as the minimum inter-cluster
#' distance divided by the maximum intra-cluster distance.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param distmat A numeric matrix containing pairwise distances between
#'   observations.
#'
#' @return A numeric Dunn index. Larger values indicate greater separation
#'   between clusters relative to within-cluster dispersion.
#'
#' @references
#' Dunn JC. Well-separated clusters and optimal fuzzy partitions.
#' Journal of Cybernetics. 1974;4(1):95-104.
#'
#' Pihur V, Datta S, Datta S. Weighted rank aggregation of cluster validation
#' measures: a Monte Carlo cross-entropy approach.
#' Bioinformatics. 2007;23(13):1607-1615.
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
#' distmat <- as.matrix(stats::dist(t(x)))
#' validation_dunn(clusters, distmat)
#'
#' @export
validation_dunn <- function(clusters, distmat) {

  DunnIndexCpp(
    ocp = clusters,
    distmat = distmat
  )
}