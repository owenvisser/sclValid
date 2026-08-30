#' Calculate Biological Stability Index
#'
#' Calculates an adapted Biological Stability Index (BSI), which evaluates
#' the stability of clustering assignments within known biological groups
#' across dataset perturbations.
#'
#' The implementation is adapted from the biological stability measure
#' described by Datta and Datta.
#'
#' @param clusters A vector containing cluster assignments obtained from the
#'   full dataset.
#' @param labels A vector containing the known label for each observation.
#' @param reduced_clusters A list of clustering-assignment vectors obtained
#'   after removing observations and reclustering the reduced datasets.
#'
#' @return A numeric biological stability index between zero and one, with
#'   larger values indicating greater stability with respect to known labels.
#'
#' @references
#' Adapted from:
#' Datta S, Datta S. Methods for evaluating clustering algorithms for gene
#' expression data using a reference set of functional classes.
#' BMC Bioinformatics. 2006;7:1-9.
#'
#' @examples
#' clusters <- c(1, 1, 2, 2)
#' labels <- c("A", "A", "B", "B")
#'
#' reduced_clusters <- list(
#'   c(1, 2, 2),
#'   c(1, 2, 2),
#'   c(1, 1, 2),
#'   c(1, 1, 2)
#' )
#'
#' validation_bsi(
#'   clusters = clusters,
#'   labels = labels,
#'   reduced_clusters = reduced_clusters
#' )
#' 
#' @export
validation_bsi <- function(
    clusters,
    labels,
    reduced_clusters
) {

  labels_numeric <- as.numeric(
    as.factor(labels)
  )

  BioSIndexCpp(
    oc = clusters,
    bd = labels_numeric,
    rc = reduced_clusters
  )
}