.scale_validation <- function(x, lower = -1e6, upper = 1e6) {

  x <- pmax(
    lower,
    pmin(x, upper)
  )

  range_x <- max(x) - min(x)

  if (range_x == 0) {
    return(rep(1, length(x)))
  }

  (x - min(x)) / range_x
}

#' Scale Clustering Validation Measures
#'
#' Scales clustering validation measures to the interval from zero to one so
#' that larger values consistently represent better clustering performance.
#'
#' Validation measures for which smaller values indicate better performance
#' are reversed after scaling. Before scaling, numeric values are bounded to
#' the interval from `-1e6` to `1e6`.
#'
#' Raw validation scores stored in `metadata(sce)$validation` are not modified.
#' The scaled results are stored separately in
#' `metadata(sce)$validation_scaled`.
#'
#' @param sce A `SingleCellExperiment` containing validation results produced
#'   by `run_validation()`.
#'
#' @return The input `SingleCellExperiment` with scaled validation results
#'   added to `metadata(sce)$validation_scaled`.
#'
#' @details
#' The measures DI, IGP, SW, ARI, BSI, and BHI are treated as measures for
#' which larger values indicate better performance.
#'
#' AD, ADM, APN, and CN are treated as measures for which smaller values
#' indicate better performance and are reversed after scaling.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce <- run_raceid(sce, cluster_sizes = 2:3)
#' sce <- run_validation(sce, measures = c("ARI", "DI", "IGP", "SW"))
#' sce <- scale_validation(sce)
#'
#' S4Vectors::metadata(sce)$validation_scaled
#'
#' @export
scale_validation <- function(sce) {

  validation <- S4Vectors::metadata(sce)$validation

  if (is.null(validation)) {
    stop(
      "No validation results were found. ",
      "Run `run_validation()` first."
    )
  }

  higher_better <- c(
    "DI",
    "IGP",
    "SW",
    "ARI",
    "BSI",
    "BHI"
  )

  lower_better <- c(
    "AD",
    "ADM",
    "APN",
    "CN"
  )

  present_higher <- intersect(
    higher_better,
    colnames(validation)
  )

  present_lower <- intersect(
    lower_better,
    colnames(validation)
  )

  scaled <- validation

  for (measure in present_higher) {

    scaled[[measure]] <- .scale_validation(
      as.numeric(validation[[measure]])
    )
  }

  for (measure in present_lower) {

    scaled[[measure]] <- 1 - .scale_validation(
      as.numeric(validation[[measure]])
    )
  }

  metadata_sce <- S4Vectors::metadata(sce)
  metadata_sce$validation_scaled <- scaled
  S4Vectors::metadata(sce) <- metadata_sce

  return(sce)
}