#' Calculate Average Distance Between Means
#'
#' Calculates an adapted Average Distance Between Means (ADM) stability
#' measure based on distances between cluster centroids in full and perturbed
#' clustering solutions.
#'
#' The implementation is adapted from the clustering stability measures
#' described by Datta and Datta.
#'
#' @param clusters A vector containing cluster assignments obtained from the
#'   full dataset.
#' @param reduced_clusters A list of clustering-assignment vectors obtained
#'   from the reduced datasets.
#' @param data A numeric matrix with features in rows and observations in
#'   columns.
#' @param index A list specifying the observations removed for each reduced
#'   clustering solution.
#'
#' @return A numeric average-distance-between-means stability score. Lower
#'   values indicate greater stability in cluster locations.
#'
#' @references
#' Adapted from:
#' Datta S, Datta S. Comparisons and validation of statistical clustering
#' techniques for microarray gene expression data.
#' Bioinformatics. 2003;19(4):459-466.
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
#'
#' reduced_clusters <- list(
#'   c(1, 2, 2),
#'   c(1, 2, 2),
#'   c(1, 1, 2),
#'   c(1, 1, 2)
#' )
#'
#' index <- as.list(1:4)
#'
#' validation_average_distance_means(
#'   clusters,
#'   reduced_clusters,
#'   x,
#'   index
#' )
#'
#' @export
validation_average_distance_means <- function(
    clusters,
    reduced_clusters,
    data,
    index
) {

  score <- 0
  n <- length(index)
  full_index <- seq_along(clusters)

  for (i in seq_len(n)) {

    reduced_i <- reduced_clusters[[i]]
    full_i <- clusters[-index[[i]]]
    n_i <- length(full_i)

    retained_index <- full_index[-index[[i]]]

    tab <- as.data.frame(
      table(reduced_i, full_i)
    )

    tab <- tab[tab[, 3] != 0, ]

    values <- apply(
      tab,
      1,
      function(x) {

        x <- as.numeric(x)

        reduced_members <- retained_index[
          which(reduced_i == x[1])
        ]

        full_members <- retained_index[
          which(full_i == x[2])
        ]

        reduced_centroid <- rowMeansC(
          as.matrix(data[, reduced_members, drop = FALSE])
        )

        full_centroid <- rowMeansC(
          as.matrix(data[, full_members, drop = FALSE])
        )

        x[3] * euclideanDistance(
          reduced_centroid,
          full_centroid
        )
      }
    )

    score <- score + sum(values) / n_i
  }

  score / n
}