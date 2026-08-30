#' Calculate Average Silhouette Width
#'
#' Calculates the average Silhouette Width (SW) for a clustering solution.
#' For each observation, the silhouette value compares the average distance
#' to observations in its own cluster with the average distance to the nearest
#' alternative cluster.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param distmat A numeric matrix containing pairwise distances between
#'   observations.
#'
#' @return A numeric average silhouette width. Values closer to one indicate
#'   better-separated and more compact clusters.
#'
#' @references
#' Rousseeuw PJ. Silhouettes: a graphical aid to the interpretation and
#' validation of cluster analysis.
#' Journal of Computational and Applied Mathematics. 1987;20:53-65.
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
#' validation_silhouette(clusters, distmat)
#'
#' @export
validation_silhouette <- function(clusters, distmat) {

  SilhouetteDistancecpp(
    ocp = clusters,
    distmat = distmat
  )
}