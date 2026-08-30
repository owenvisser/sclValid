.aggregate_ranks <- function(scores,
                             importance = NULL,
                             k = NULL,
                             seed = 0) {

  scores <- as.matrix(scores)

  if (is.null(importance)) {
    importance <- rep(1, nrow(scores))
  }

  if (length(importance) != nrow(scores)) {
    stop(
      "`importance` must have one value for each validation measure."
    )
  }

  if (anyNA(scores)) {
    stop(
      "Rank aggregation cannot be performed when validation scores contain NA values."
    )
  }

  if (ncol(scores) == 1) {
    return(colnames(scores))
  }

  rank_lists <- t(
    apply(
      scores,
      1,
      function(x) {
        colnames(scores)[order(x, decreasing = TRUE)]
      }
    )
  )

  rank_weights <- t(
    apply(
      scores,
      1,
      function(x) {
        x[order(x, decreasing = TRUE)]
      }
    )
  )

  if (is.null(k)) {
    k <- ncol(scores)
  }

  fit <- RankAggreg::RankAggreg(
    x = rank_lists,
    seed = seed,
    weights = rank_weights,
    k = min(k, ncol(scores)),
    importance = importance,
    standardizeWeights = FALSE,
    N = max(100, 10 * k^2),
    verbose = FALSE
  )

  fit$top.list
}

#' Rank Clustering Solutions
#'
#' Aggregates scaled clustering validation measures to produce an overall
#' ranking of clustering solutions.
#'
#' Rank aggregation is performed in two stages. First, clustering solutions
#' are compared separately within each cluster size, and the highest-ranked
#' solution is retained for each size. The retained solutions are then
#' aggregated across cluster sizes to produce the final ranking.
#'
#' Validation scores are obtained from
#' `S4Vectors::metadata(sce)$validation_scaled`. All scores are assumed to
#' have been scaled so that larger values indicate better clustering
#' performance.
#'
#' @param sce A `SingleCellExperiment` containing scaled validation results
#'   produced by `scale_validation()`.
#' @param measures Optional character vector specifying validation measures
#'   to include in the aggregation. If `NULL`, all validation measures
#'   available in the scaled validation table are used.
#' @param importance Optional numeric vector specifying the relative
#'   importance of the selected validation measures. If `NULL`, all measures
#'   receive equal importance.
#' @param seed Integer random seed supplied to `RankAggreg`.
#'
#' @return The input `SingleCellExperiment` with the final aggregated ranking
#'   stored in `metadata(sce)$ranking`.
#'
#' @references
#' Pihur V, Datta S, Datta S. Weighted rank aggregation of cluster validation
#' measures: a Monte Carlo cross-entropy approach.
#' Bioinformatics. 2007;23(13):1607-1615.
#'
#' Pihur V, Datta S, Datta S. RankAggreg, an R package for weighted rank
#' aggregation.
#' BMC Bioinformatics. 2009;10:1-10.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#'
#' sce <- run_cidr(
#'   sce,
#'   cluster_sizes = 2:3
#' )
#'
#' sce <- run_validation(
#'   sce,
#'   measures = c("ARI", "BHI", "DI", "IGP", "SW"),
#'   verbose = FALSE,
#'   progress = FALSE
#' )
#'
#' sce <- scale_validation(sce)
#'
#' sce <- run_rank_aggregation(sce)
#'
#' S4Vectors::metadata(sce)$ranking
#' 
#' @export
run_rank_aggregation <- function(
    sce,
    measures = NULL,
    importance = NULL,
    seed = 0
) {

  validation <- S4Vectors::metadata(sce)$validation_scaled

  if (is.null(validation)) {
    stop(
      "No scaled validation results were found. Run `scale_validation()` first."
    )
  }

  measure_columns <- setdiff(
    colnames(validation),
    c("clustering", "cluster_size")
  )

  if (is.null(measures)) {
    measures <- measure_columns
  } else {
    missing_measures <- setdiff(measures, measure_columns)

    if (length(missing_measures) > 0) {
      stop(
        "The following validation measures were not found: ",
        paste(missing_measures, collapse = ", ")
      )
    }
  }

  if (length(measures) == 0) {
    stop("At least one validation measure must be selected.")
  }

  if (is.null(importance)) {
    importance <- rep(1, length(measures))
  }

  if (length(importance) != length(measures)) {
    stop(
      "`importance` must have one value for each selected validation measure."
    )
  }

  sizes <- unique(validation$cluster_size)

  size_winners <- character(length(sizes))

  for (i in seq_along(sizes)) {

    size_i <- sizes[i]

    dat_i <- validation[
      validation$cluster_size == size_i,
      ,
      drop = FALSE
    ]

    scores_i <- t(
      as.matrix(
        dat_i[, measures, drop = FALSE]
      )
    )

    colnames(scores_i) <- dat_i$clustering

    ranking_i <- .aggregate_ranks(
      scores = scores_i,
      importance = importance,
      k = min(ncol(scores_i), 10),
      seed = seed
    )

    size_winners[i] <- ranking_i[1]
  }

  finalists <- validation[
    validation$clustering %in% size_winners,
    ,
    drop = FALSE
  ]

  final_scores <- t(
    as.matrix(
      finalists[, measures, drop = FALSE]
    )
  )

  colnames(final_scores) <- finalists$clustering

  final_ranking <- .aggregate_ranks(
    scores = final_scores,
    importance = importance,
    k = min(ncol(final_scores), 20),
    seed = seed
  )

  ranking <- finalists[
    match(final_ranking, finalists$clustering),
    c("clustering", "cluster_size", measures),
    drop = FALSE
  ]

  ranking <- data.frame(
    rank = seq_len(nrow(ranking)),
    ranking,
    row.names = NULL,
    check.names = FALSE
  )

  metadata_sce <- S4Vectors::metadata(sce)
  metadata_sce$ranking <- ranking
  S4Vectors::metadata(sce) <- metadata_sce

  sce
}