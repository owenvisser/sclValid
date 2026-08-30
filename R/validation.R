#' Calculate Clustering Validation Measures
#'
#' Calculates internal, external, and stability-based validation measures for
#' clustering solutions stored in the `colData` of a
#' `SingleCellExperiment`.
#'
#' Clustering solutions are detected automatically from columns produced by the
#' built-in clustering functions in `sclValid`, as well as from user-defined
#' clustering methods registered with `run_custom_clustering()`. This allows
#' clustering methods not implemented directly in `sclValid` to be evaluated
#' using the same validation framework. A subset of clustering solutions may
#' instead be selected with `clusterings`.
#'
#' Shared quantities such as pairwise distances and nearest-neighbor lists are
#' calculated once and reused across validation measures. When stability
#' measures are requested, reduced datasets are reclustered internally and the
#' resulting reduced clustering solutions are used for comparison with the
#' full-data solutions.
#'
#' Validation results are stored in
#' `S4Vectors::metadata(sce)$validation`, with one row per clustering solution
#' and one column per requested validation measure.
#'
#' @param sce A `SingleCellExperiment` containing a `"data"` assay and one or
#'   more clustering-assignment columns in `colData`.
#' @param measures Character vector specifying validation measures to
#'   calculate. Supported values are `"AD"`, `"ADM"`, `"APN"`, `"ARI"`,
#'   `"BHI"`, `"BSI"`, `"CN"`, `"DI"`, `"IGP"`, and `"SW"`.
#' @param clusterings Optional character vector specifying particular
#'   clustering columns in `colData` to validate. If `NULL`, all recognized
#'   clustering solutions are used.
#' @param index Optional list specifying observations removed for stability
#'   calculations. If `NULL`, leave-one-observation-out removal is used.
#' @param nn_k Number of nearest neighbors to calculate when a requested
#'   validation measure requires nearest-neighbor information. Defaults to
#'   `50`.
#' @param connectivity_h Number of nearest neighbors used by the connectivity
#'   measure. Defaults to `5`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during validation. Defaults to `FALSE`.
#' @param progress Logical indicating whether a progress bar should be
#'   displayed while reduced clustering solutions are generated for
#'   stability-based validation measures. Defaults to `TRUE`.
#'
#' @return The input `SingleCellExperiment` with validation results added to
#'   `metadata(sce)$validation`.
#'
#' @details
#' ARI, BHI, and BSI require known observation labels stored in the `"label"`
#' column of `colData`. If complete labels are not available, these measures
#' are skipped while validation measures that do not require known labels may
#' still be calculated.
#'
#' AD, ADM, APN, and BSI require reduced clustering solutions and therefore
#' trigger reclustering of reduced versions of the dataset.
#'
#' User-defined clustering methods may be registered with
#' `run_custom_clustering()`. The supplied clustering function should accept
#' the expression matrix as its first argument, with features in rows and
#' observations in columns, and must return one cluster assignment per
#' observation. Additional arguments may be supplied through
#' `run_custom_clustering()`.
#'
#' A general custom clustering workflow is:
#'
#' \preformatted{
#' my_clustering <- function(data, ...) {
#'   clusters <- clustering_method(data, ...)
#'   clusters
#' }
#'
#' sce <- run_custom_clustering(
#'   sce,
#'   fun = my_clustering,
#'   name = "MyClustering"
#' )
#'
#' sce <- run_validation(sce)
#' }
#'
#' The custom function and its supplied arguments are retained internally so
#' that stability-based validation measures can rerun the same clustering
#' procedure on reduced datasets.
#'
#' By default, `run_validation()` displays a progress bar during reduced
#' reclustering while suppressing detailed clustering messages. Set
#' `progress = FALSE` to disable the progress bar, or `verbose = TRUE` to
#' display additional status messages. The progress bar is only displayed
#' for stability-based validation measures that require reduced clustering
#' solutions.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce <- run_raceid(sce, cluster_sizes = 2:3)
#'
#' sce <- run_validation(
#'   sce,
#'   measures = c("ARI", "BHI", "DI", "IGP", "SW")
#' )
#'
#' S4Vectors::metadata(sce)$validation
#'
#' @export
run_validation <- function(
    sce,
    measures = c(
      "AD", "ADM", "APN",
      "ARI", "BHI", "BSI",
      "CN", "DI", "IGP", "SW"
    ),
    clusterings = NULL,
    index = NULL,
    nn_k = 50,
    connectivity_h = 5,
    verbose = FALSE,
    progress = TRUE
) {

  measures <- toupper(measures)

  valid_measures <- c(
    "AD", "ADM", "APN",
    "ARI", "BHI", "BSI",
    "CN", "DI", "IGP", "SW"
  )

  if (!all(measures %in% valid_measures)) {
    stop(
      "`measures` must contain only: ",
      paste(valid_measures, collapse = ", ")
    )
  }

  data <- SummarizedExperiment::assay(sce, "data")
  cd <- SummarizedExperiment::colData(sce)

  # ------------------------------------------------------------
  # Find clustering solutions
  # ------------------------------------------------------------

  cluster_names <- grep(
    paste0(
      "^(",
      "pcaReduceM_cs|",
      "pcaReduceS_cs|",
      "RaceID_cs|",
      "PCA_|",
      "TSNE_",
      ")"
    ),
    colnames(cd),
    value = TRUE
  )

  custom_clusterings <-
    S4Vectors::metadata(sce)$custom_clusterings

  if (!is.null(custom_clusterings)) {
    cluster_names <- unique(
      c(
        cluster_names,
        names(custom_clusterings)
      )
    )
  }

  if (length(cluster_names) == 0) {
    stop("No clustering results were found in `colData(sce)`.")
  }

  # If specific clustering solutions were requested, keep only those
  if (!is.null(clusterings)) {

    missing_clusterings <- setdiff(
      clusterings,
      cluster_names
    )

    if (length(missing_clusterings) > 0) {
      stop(
        "The following clustering results were not found in `colData(sce)`: ",
        paste(missing_clusterings, collapse = ", ")
      )
    }

    cluster_names <- clusterings
  }

# ------------------------------------------------------------
# Check whether known labels are available
# ------------------------------------------------------------

  has_labels <-
    "label" %in% colnames(cd) &&
    length(sce$label) == ncol(sce) &&
    !any(is.na(sce$label)) &&
    !any(as.character(sce$label) == "")

  known_measures <- c("ARI", "BHI", "BSI")

  if (!has_labels && any(measures %in% known_measures)) {

    warning(
      "No complete `label` column was found. ",
      "ARI, BHI, and BSI will be skipped."
    )

    measures <- setdiff(
      measures,
      known_measures
    )
  }

  # ------------------------------------------------------------
  # Determine which shared calculations are needed
  # ------------------------------------------------------------

  need_dist <- any(
    measures %in% c("AD", "DI", "SW")
  )

  need_nn <- any(
    measures %in% c("CN", "IGP")
  )

  need_stability <- any(
    measures %in% c("AD", "ADM", "APN", "BSI")
  )

  if (need_dist) {

    if (verbose) {
      message("Calculating expression distance matrix...")
    }

    distmat <- dist_mat(data)
  }

  if (need_nn) {

    if (verbose) {
      message("Calculating nearest neighbors...")
    }

    nn_k_use <- min(
      nn_k,
      ncol(sce) - 1
    )

    nn <- list_nn(
      data,
      k = nn_k_use
    )
  }

  # ------------------------------------------------------------
  # Stability removal index
  # ------------------------------------------------------------

  if (need_stability) {

    if (is.null(index)) {

      # Default: leave one observation out at a time
      index <- as.list(
        seq_len(ncol(sce))
      )
    }

    if (verbose) {
      message(
        "Generating reduced clustering solutions for ",
        length(index),
        " removal sets..."
      )
    }

    reduced <- .validation_recluster(
      sce = sce,
      cluster_names = cluster_names,
      index = index,
      verbose = verbose,
      progress = progress
    )
  }

  # ------------------------------------------------------------
  # Calculate validation measures
  # ------------------------------------------------------------

  validation <- vector(
    "list",
    length(cluster_names)
  )

  names(validation) <- cluster_names

  if (verbose) {
    message("Calculating validation measures...")
  }

  for (cluster_name in cluster_names) {

    clusters <- cd[[cluster_name]]

    result <- list(
      clustering = cluster_name,
      cluster_size = length(unique(clusters))
    )

    if ("AD" %in% measures) {

      result$AD <- validation_average_distance(
        clusters = clusters,
        reduced_clusters = reduced[[cluster_name]],
        distmat = distmat,
        index = index
      )
    }

    if ("ADM" %in% measures) {

      result$ADM <- validation_average_distance_means(
        clusters = clusters,
        reduced_clusters = reduced[[cluster_name]],
        data = data,
        index = index
      )
    }

    if ("APN" %in% measures) {

      result$APN <- validation_apn(
        clusters = clusters,
        reduced_clusters = reduced[[cluster_name]],
        index = index
      )
    }

    if ("ARI" %in% measures) {

      result$ARI <- validation_ari(
        clusters = clusters,
        labels = sce$label
      )
    }

    if ("BHI" %in% measures) {

      result$BHI <- validation_bhi(
        clusters = clusters,
        labels = sce$label
      )
    }

    if ("BSI" %in% measures) {

      result$BSI <- validation_bsi(
        clusters = clusters,
        labels = sce$label,
        reduced_clusters = reduced[[cluster_name]]
      )
    }

    if ("CN" %in% measures) {

      result$CN <- validation_connectivity(
        clusters = clusters,
        nn = nn,
        h = connectivity_h
      )
    }

    if ("DI" %in% measures) {

      result$DI <- validation_dunn(
        clusters = clusters,
        distmat = distmat
      )
    }

    if ("IGP" %in% measures) {

      result$IGP <- validation_igp(
        clusters = clusters,
        nn = nn
      )
    }

    if ("SW" %in% measures) {

      result$SW <- validation_silhouette(
        clusters = clusters,
        distmat = distmat
      )
    }

    validation[[cluster_name]] <- result
  }

  validation <- do.call(
    rbind,
    lapply(
      validation,
      function(x) {
        as.data.frame(
          x,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(validation) <- NULL

  metadata_sce <- S4Vectors::metadata(sce)
  metadata_sce$validation <- validation
  S4Vectors::metadata(sce) <- metadata_sce

  return(sce)
}


.validation_recluster <- function(
    sce,
    cluster_names,
    index,
    verbose = TRUE,
    progress = TRUE
) {

  data <- SummarizedExperiment::assay(
    sce,
    "data"
  )

  has_labels <-
    "label" %in%
    colnames(
      SummarizedExperiment::colData(sce)
    )

  reduced <- stats::setNames(
    vector(
      "list",
      length(cluster_names)
    ),
    cluster_names
  )

  if (progress) {
    cli::cli_progress_bar(
      "Creating reduced clustering solutions",
      total = length(index)
    )
  }

  for (i in seq_along(index)) {

    remove <- index[[i]]

    keep <- setdiff(
      seq_len(ncol(sce)),
      remove
    )

    reduced_sce <-
      SingleCellExperiment::SingleCellExperiment(
        assays = list(
          data = data[, keep, drop = FALSE]
        )
      )

    if (has_labels) {
      reduced_sce$label <- sce$label[keep]
    }

    # ----------------------------------------------------------
    # pcaReduce M
    # ----------------------------------------------------------

    prm_names <- grep(
      "^pcaReduceM_cs",
      cluster_names,
      value = TRUE
    )

    if (length(prm_names) > 0) {

      sizes <- as.integer(
        sub(
          "^pcaReduceM_cs",
          "",
          prm_names
        )
      )

      reduced_sce <- run_pcareduce(
        reduced_sce,
        cluster_sizes = sizes,
        method = "M",
        verbose = FALSE
      )
    }

    # ----------------------------------------------------------
    # pcaReduce S
    # ----------------------------------------------------------

    prs_names <- grep(
      "^pcaReduceS_cs",
      cluster_names,
      value = TRUE
    )

    if (length(prs_names) > 0) {

      sizes <- as.integer(
        sub(
          "^pcaReduceS_cs",
          "",
          prs_names
        )
      )

      reduced_sce <- run_pcareduce(
        reduced_sce,
        cluster_sizes = sizes,
        method = "S",
        verbose = FALSE
      )
    }

    # ----------------------------------------------------------
    # RaceID
    # ----------------------------------------------------------

    raceid_names <- grep(
      "^RaceID_cs",
      cluster_names,
      value = TRUE
    )

    if (length(raceid_names) > 0) {

      sizes <- as.integer(
        sub(
          "^RaceID_cs",
          "",
          raceid_names
        )
      )

      reduced_sce <- run_raceid(
        reduced_sce,
        cluster_sizes = sizes,
        verbose = FALSE
      )
    }

    # ----------------------------------------------------------
    # PCA-based methods
    # ----------------------------------------------------------

    pca_names <- grep(
      "^PCA_",
      cluster_names,
      value = TRUE
    )

    if (length(pca_names) > 0) {

      reduced_sce <- .validation_run_reduction_set(
        sce = reduced_sce,
        names = pca_names,
        reduction = "PCA",
        verbose = FALSE
      )
    }

    # ----------------------------------------------------------
    # t-SNE-based methods
    # ----------------------------------------------------------

    tsne_names <- grep(
      "^TSNE_",
      cluster_names,
      value = TRUE
    )

    if (length(tsne_names) > 0) {

      reduced_sce <- .validation_run_reduction_set(
        sce = reduced_sce,
        names = tsne_names,
        reduction = "TSNE",
        verbose = FALSE
      )
    }

    # ----------------------------------------------------------
    # Custom clustering methods
    # ----------------------------------------------------------

    custom_clusterings <-
      S4Vectors::metadata(sce)$custom_clusterings

    if (!is.null(custom_clusterings)) {

      custom_names <- intersect(
        cluster_names,
        names(custom_clusterings)
      )

      if (length(custom_names) > 0) {

        reduced_data <-
          SummarizedExperiment::assay(
            reduced_sce,
            "data"
          )

        for (nm in custom_names) {

          custom <- custom_clusterings[[nm]]

          assignments <- do.call(
            custom$fun,
            c(list(reduced_data),custom$args)
          )
          if (length(assignments) != ncol(reduced_sce)) {
            stop(
              "Custom clustering function `",
              nm,
              "` did not return one cluster assignment per observation."
            )
          }

          reduced_sce[[nm]] <- assignments
        }
      }
    }

    # ----------------------------------------------------------
    # Save only the requested reduced solutions
    # ----------------------------------------------------------

    reduced_cd <-
      SummarizedExperiment::colData(
        reduced_sce
      )

    for (nm in cluster_names) {

      reduced[[nm]][[i]] <-
        reduced_cd[[nm]]
    }
    if (progress) {
      cli::cli_progress_update()
    }
  }

  if (progress) {
    cli::cli_progress_done()
  }

  reduced
}


.validation_run_reduction_set <- function(
    sce,
    names,
    reduction,
    verbose = FALSE
) {

  methods <- character()
  cluster_sizes <- integer()
  knn <- numeric()
  resolution <- numeric()

  # K-means
  km_names <- grep(
    paste0("^", reduction, "_Kmeans_cs"),
    names,
    value = TRUE
  )

  if (length(km_names) > 0) {

    methods <- c(methods, "kmeans")

    cluster_sizes <- c(
      cluster_sizes,
      as.integer(
        sub(
          paste0("^", reduction, "_Kmeans_cs"),
          "",
          km_names
        )
      )
    )
  }

  # Hierarchical
  hc_names <- grep(
    paste0("^", reduction, "_HC_cs"),
    names,
    value = TRUE
  )

  if (length(hc_names) > 0) {

    methods <- c(
      methods,
      "hierarchical"
    )

    cluster_sizes <- c(
      cluster_sizes,
      as.integer(
        sub(
          paste0("^", reduction, "_HC_cs"),
          "",
          hc_names
        )
      )
    )
  }

  # Louvain / Leiden
  graph_names <- grep(
    paste0(
      "^",
      reduction,
      "_(Louvain|Leiden)_knn"
    ),
    names,
    value = TRUE
  )

  if (length(graph_names) > 0) {

    if (
      any(
        grepl(
          paste0("^", reduction, "_Louvain_"),
          graph_names
        )
      )
    ) {
      methods <- c(methods, "louvain")
    }

    if (
      any(
        grepl(
          paste0("^", reduction, "_Leiden_"),
          graph_names
        )
      )
    ) {
      methods <- c(methods, "leiden")
    }

    knn <- as.numeric(
      sub(
        ".*_knn([0-9]+)_res.*",
        "\\1",
        graph_names
      )
    )

    resolution <- as.numeric(
      sub(
        ".*_res",
        "",
        graph_names
      )
    )
  }

  methods <- unique(methods)
  cluster_sizes <- unique(cluster_sizes)
  knn <- unique(knn)
  resolution <- unique(resolution)

  # Values are irrelevant if the corresponding
  # method is not requested.
  if (length(cluster_sizes) == 0) {
    cluster_sizes <- 2L
  }

  if (length(knn) == 0) {
    knn <- 10
  }

  if (length(resolution) == 0) {
    resolution <- 1
  }

  run_reduction_clustering(
    sce,
    reductions = reduction,
    methods = methods,
    cluster_sizes = cluster_sizes,
    knn = knn,
    resolution = resolution,
    verbose = verbose
  )
}