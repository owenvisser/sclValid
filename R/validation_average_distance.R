#' Calculate Average Distance
#'
#' Calculates an adapted Average Distance (AD) stability measure by comparing
#' clustering solutions obtained from the full dataset with corresponding
#' solutions obtained after perturbing the dataset.
#'
#' The implementation is adapted from the clustering stability measures
#' described by Datta and Datta.
#'
#' @param clusters A vector containing cluster assignments obtained from the
#'   full dataset.
#' @param reduced_clusters A list of clustering-assignment vectors obtained
#'   from the reduced datasets.
#' @param distmat A numeric matrix containing pairwise distances between
#'   observations in the full dataset.
#' @param index A list specifying the observations removed for each reduced
#'   clustering solution.
#'
#' @return A numeric average-distance stability score. Lower values indicate
#'   greater agreement between full and perturbed clustering solutions.
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
#' distmat <- as.matrix(stats::dist(t(x)))
#'
#' validation_average_distance(
#'   clusters,
#'   reduced_clusters,
#'   distmat,
#'   index
#' )
#'
#' @export
validation_average_distance <- function(
    clusters,
    reduced_clusters,
    distmat,
    index
) {

  score <- 0
  n <- length(index)

  for (i in seq_len(n)) {

    reduced_i <- reduced_clusters[[i]]
    full_i <- clusters[-index[[i]]]
    n_i <- length(full_i)

    tab <- as.data.frame(
      table(reduced_i, full_i)
    )

    tab <- tab[tab[, 3] != 0, ]

    values <- apply(
      tab,
      1,
      function(x) {

        x <- as.numeric(x)

        reduced_index <- which(reduced_i == x[1])
        full_index <- which(full_i == x[2])

        x[3] * mean(
          distmat[reduced_index, full_index]
        )
      }
    )

    score <- score + sum(values) / n_i
  }

  score / n
}