#' Run Clustering on PCA and t-SNE Representations
#'
#' Performs dimensionality reduction followed by one or more clustering
#' procedures on the `"data"` assay of a `SingleCellExperiment` object.
#' Principal component analysis (PCA) and t-distributed stochastic neighbor
#' embedding (t-SNE) representations are computed only when requested and are
#' reused across clustering procedures to avoid repeated dimensionality
#' reduction and distance calculations.
#'
#' Supported clustering methods include K-means, hierarchical clustering,
#' Louvain community detection, and Leiden community detection. Cluster
#' assignments are added to the object's `colData`.
#'
#' @param sce A `SingleCellExperiment` object containing a `"data"` assay.
#' @param reductions Character vector specifying the dimensionality reductions
#'   to use. Supported values are `"PCA"` and `"TSNE"`.
#' @param methods Character vector specifying the clustering methods to apply.
#'   Supported values are `"kmeans"`, `"hierarchical"`, `"louvain"`, and
#'   `"leiden"`.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   used for K-means and hierarchical clustering.
#' @param knn An integer vector specifying the numbers of nearest neighbors
#'   used to construct graphs for Louvain and Leiden clustering.
#' @param resolution A numeric vector specifying resolution values for Louvain
#'   and Leiden clustering.
#' @param pca_dims Maximum number of principal components to retain. The actual
#'   number is limited by the dimensions of the input data. Defaults to `50`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return The input `SingleCellExperiment` object with additional
#'   cluster-assignment columns in `colData`.
#'
#' @details
#' K-means and hierarchical clustering results are named according to the
#' dimensionality reduction and requested number of clusters. For example,
#' `"PCA_Kmeans_cs3"` contains PCA-based K-means assignments for three
#' clusters.
#'
#' Louvain and Leiden results include both the nearest-neighbor and resolution
#' parameters in their column names, for example
#' `"PCA_Louvain_knn10_res0.5"`.
#'
#' When t-SNE is requested, PCA distances are first calculated and used as the
#' input distance representation for t-SNE.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#'
#' sce <- run_reduction_clustering(
#'   sce,
#'   reductions = c("PCA", "TSNE"),
#'   methods = c("kmeans", "hierarchical"),
#'   cluster_sizes = 2:4
#' )
#'
#' SummarizedExperiment::colData(sce)
#'
#' @export
run_reduction_clustering <- function(
    sce,
    reductions = c("PCA", "TSNE"),
    methods = c("kmeans", "hierarchical", "louvain", "leiden"),
    cluster_sizes = 2:4,
    knn = c(10, 25, 40),
    resolution = seq(0.1, 0.9, 0.1),
    pca_dims = 50,
    verbose = TRUE
) {

  # Extract the data matrix
  d <- SummarizedExperiment::assay(sce, "data")

  # Standardize user input
  reductions <- toupper(reductions)
  methods <- tolower(methods)

  # Check requested reductions
  valid_reductions <- c("PCA", "TSNE")

  if (!all(reductions %in% valid_reductions)) {
    stop(
      "`reductions` must contain only: ",
      paste(valid_reductions, collapse = ", ")
    )
  }

  # Check requested clustering methods
  valid_methods <- c(
    "kmeans",
    "hierarchical",
    "louvain",
    "leiden"
  )

  if (!all(methods %in% valid_methods)) {
    stop(
      "`methods` must contain only: ",
      paste(valid_methods, collapse = ", ")
    )
  }

  # Determine whether distances are needed
  need_dist <- any(
    methods %in% c(
      "hierarchical",
      "louvain",
      "leiden"
    )
  )

  # Determine whether a full distance matrix is needed
  need_dmat <- any(
    methods %in% c(
      "louvain",
      "leiden"
    )
  )

  # Helper for adding returned clustering vectors to colData
  add_results <- function(sce, results) {

    for (nm in names(results)) {
      sce[[nm]] <- results[[nm]]
    }

    sce
  }


  # ============================================================
  # PCA
  # ============================================================

  # PCA is always required when PCA is requested.
  # It is also required internally to calculate t-SNE.
  if ("PCA" %in% reductions || "TSNE" %in% reductions) {

    if (verbose) {
      message("Calculating PCA...")
    }

    max_pca_dims <- min(
      pca_dims,
      ncol(d) - 1,
      nrow(d)
    )

    pca <- stats::prcomp(
      x = t(d),
      rank. = max_pca_dims
    )$x
  }


  # ============================================================
  # PCA clustering
  # ============================================================

  if ("PCA" %in% reductions) {

    if ("kmeans" %in% methods) {

      results <- .run_kmeans(
        x = pca,
        cluster_sizes = cluster_sizes,
        prefix = "PCA",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if (need_dist) {

      if (verbose) {
        message("Calculating PCA distances...")
      }

      pca_dist <- stats::dist(
        pca,
        method = "euclidean"
      )
    }


    if ("hierarchical" %in% methods) {

      results <- .run_hierarchical(
        d = pca_dist,
        cluster_sizes = cluster_sizes,
        prefix = "PCA",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if (need_dmat) {
      pca_dmat <- as.matrix(pca_dist)
    }


    if ("louvain" %in% methods) {

      results <- .run_louvain(
        dmat = pca_dmat,
        knn = knn,
        resolution = resolution,
        prefix = "PCA",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if ("leiden" %in% methods) {

      results <- .run_leiden(
        dmat = pca_dmat,
        knn = knn,
        resolution = resolution,
        prefix = "PCA",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }
  }


  # ============================================================
  # t-SNE
  # ============================================================

  if ("TSNE" %in% reductions) {

    # Your original t-SNE procedure uses PCA distances
    if (!exists("pca_dist", inherits = FALSE)) {

      if (verbose) {
        message("Calculating PCA distances for t-SNE...")
      }

      pca_dist <- stats::dist(
        pca,
        method = "euclidean"
      )
    }

    pca_dmat <- as.matrix(pca_dist)

    if (verbose) {
      message("Calculating t-SNE...")
    }

    perplexity <- (nrow(pca_dmat) - 1) / 3

    tsne <- Rtsne::Rtsne(
      pca_dmat,
      dims = 3,
      perplexity = perplexity,
      theta = 0.5,
      check_duplicates = FALSE,
      pca = FALSE,
      is_distance = TRUE,
      Y_init = NULL,
      normalize = FALSE
    )$Y


    # ==========================================================
    # t-SNE clustering
    # ==========================================================

    if ("kmeans" %in% methods) {

      results <- .run_kmeans(
        x = tsne,
        cluster_sizes = cluster_sizes,
        prefix = "TSNE",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if (need_dist) {

      if (verbose) {
        message("Calculating t-SNE distances...")
      }

      tsne_dist <- stats::dist(
        tsne,
        method = "euclidean"
      )
    }


    if ("hierarchical" %in% methods) {

      results <- .run_hierarchical(
        d = tsne_dist,
        cluster_sizes = cluster_sizes,
        prefix = "TSNE",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if (need_dmat) {
      tsne_dmat <- as.matrix(tsne_dist)
    }


    if ("louvain" %in% methods) {

      results <- .run_louvain(
        dmat = tsne_dmat,
        knn = knn,
        resolution = resolution,
        prefix = "TSNE",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }


    if ("leiden" %in% methods) {

      results <- .run_leiden(
        dmat = tsne_dmat,
        knn = knn,
        resolution = resolution,
        prefix = "TSNE",
        verbose = verbose
      )

      sce <- add_results(sce, results)
    }
  }


  return(sce)
}