#' Calculate Average Proportion of Non-Overlap
#'
#' Calculates an adapted Average Proportion of Non-Overlap (APN) stability
#' measure by comparing cluster memberships between full and perturbed
#' clustering solutions.
#'
#' The implementation is adapted from the clustering stability measures
#' described by Datta and Datta.
#'
#' @param clusters A vector containing cluster assignments obtained from the
#'   full dataset.
#' @param reduced_clusters A list of clustering-assignment vectors obtained
#'   from the reduced datasets.
#' @param index A list specifying the observations removed for each reduced
#'   clustering solution.
#'
#' @return A numeric average proportion of non-overlap. Lower values indicate
#'   more stable clustering.
#'
#' @references
#' Adapted from:
#' Datta S, Datta S. Comparisons and validation of statistical clustering
#' techniques for microarray gene expression data.
#' Bioinformatics. 2003;19(4):459-466.
#' 
#' @examples
#' clusters <- c(1, 1, 2, 2)
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
#' validation_apn(
#'   clusters,
#'   reduced_clusters,
#'   index
#' )
#'
#' @export
validation_apn <- function(
    clusters,
    reduced_clusters,
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

        numerator <- length(
          intersect(
            which(reduced_i == x[1]),
            which(full_i == x[2])
          )
        )

        denominator <- length(
          which(full_i == x[2])
        )

        x[3] * (1 - numerator / denominator)
      }
    )

    score <- score + sum(values) / n_i
  }

  score / n
}