#' Calculate Connectivity
#'
#' Calculates the Connectivity (CN) measure for a clustering solution using a
#' nearest-neighbor representation of the data. The measure penalizes nearby
#' observations that are assigned to different clusters.
#'
#' @param clusters A vector containing one cluster assignment per observation.
#' @param nn A list containing nearest-neighbor indices for each observation.
#' @param h Number of nearest neighbors to consider. Defaults to `5`.
#'
#' @return A numeric connectivity score. Lower values indicate better
#'   agreement between cluster assignments and local neighborhood structure.
#'
#' @references
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
#' validation_connectivity(clusters, nn, h = 2)
#'
#' @export
validation_connectivity <- function(clusters, nn, h = 5) {

  score <- 0

  for (i in seq_along(clusters)) {

    check <- clusters[i] != clusters[nn[[i]][1:h]]

    if (any(check)) {
      score <- score + sum(1 / which(check))
    }
  }

  return(score)
}