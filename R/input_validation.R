#' Check Input Data
#'
#' Checks whether an input matrix and corresponding label vector are suitable
#' for use with sclValid. The matrix must be numeric, finite, and contain no
#' missing values. Rows represent features and columns represent observations.
#' The label vector must contain exactly one label for each column of the matrix.
#'
#' @param x A numeric matrix with features in rows and observations in columns.
#' @param labels A vector or factor containing one label for each column of `x`.
#'
#' @return Invisibly returns `TRUE` if all input checks are passed.
#'
#' @examples
#' check_input(example_data, example_labels)
#'
#' @export
check_input <- function(x, labels) {

  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }

  if (!is.numeric(x)) {
    stop("`x` must be a numeric matrix.")
  }

  if (nrow(x) == 0 || ncol(x) == 0) {
    stop("`x` must contain at least one row and one column.")
  }

  if (anyNA(x)) {
    stop("`x` contains missing values.")
  }

  if (any(!is.finite(x))) {
    stop("`x` contains non-finite values.")
  }

  if (!is.atomic(labels) && !is.factor(labels)) {
    stop("`labels` must be a vector or factor.")
  }

  if (length(labels) != ncol(x)) {
    stop(
      "`labels` must contain one label for each column of `x`."
    )
  }

  if (anyNA(labels)) {
    stop("`labels` contains missing values.")
  }

  invisible(TRUE)
}

#' Create a SingleCellExperiment Object
#'
#' Creates a `SingleCellExperiment` object from a numeric feature matrix and
#' a corresponding vector of observation labels. The input matrix is stored
#' as the `"data"` assay and the supplied labels are stored in `colData`.
#'
#' @param x A numeric matrix with features in rows and observations in columns.
#' @param labels A vector or factor containing one label for each column of `x`.
#'
#' @return A `SingleCellExperiment` object containing the input matrix and
#'   observation labels.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce
#' 
#' @export
make_sce <- function(x, labels) {

  check_input(x, labels)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(data = x),
    colData = data.frame(label = labels)
  )

  return(sce)
}