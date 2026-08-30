#' Run pcaReduce Clustering
#'
#' Applies the pcaReduce clustering procedure to the data assay of a
#' `SingleCellExperiment` object for one or more requested numbers of clusters.
#' Cluster assignments are added to the object's `colData`.
#'
#' The pcaReduce method may be run using either the `"M"` or `"S"` merging
#' procedure. For a cluster size `k`, resulting columns are named
#' `"pcaReduceM_cs{k}"` or `"pcaReduceS_cs{k}"`, depending on the selected
#' method.
#'
#' @param sce A `SingleCellExperiment` object containing a `"data"` assay.
#' @param cluster_sizes An integer vector specifying the numbers of clusters
#'   to fit.
#' @param method Character string specifying the pcaReduce merging method.
#'   Must be either `"M"` or `"S"`.
#' @param nbt Number of pcaReduce iterations to perform. Defaults to `100`.
#' @param verbose Logical indicating whether detailed progress messages should
#'   be displayed during clustering`.
#'
#' @return The input `SingleCellExperiment` object with additional pcaReduce
#'   cluster-assignment columns in `colData`.
#'
#' @references
#' Žurauskienė J, Yau C.
#' pcaReduce: Hierarchical clustering of single cell transcriptional profiles.
#' BMC Bioinformatics. 2016;17:140.
#'
#' @examples
#' sce <- make_sce(example_data, example_labels)
#' sce <- run_pcareduce(
#'   sce,
#'   cluster_sizes = 2:4,
#'   method = "M"
#' )
#' SummarizedExperiment::colData(sce)
#'
#' @export
run_pcareduce <- function(sce,
                          cluster_sizes,
                          method = "M",
                          nbt = 100,
                          verbose = TRUE) {

  d <- SummarizedExperiment::assay(sce, "data")

  top <- max(max(cluster_sizes) - 1, 2)
  bot <- max(min(cluster_sizes), 2)

  if (top == 2 && bot == 2) {
    index <- 2
    size_index <- (top + 1):bot
  } else {
    index <- rev(seq_along(cluster_sizes))
    size_index <- (top + 1):bot
  }

  results <- .pcareduce(
    D_t = t(d),
    q = top,
    nbt = nbt,
    method = method
  )

  for (i in index) {

    k <- size_index[i]

    if (verbose) {
      message(
        "pcaReduce ", method,
        ": clustering with k = ", k
      )
    }

    partitions <- lapply(
      results,
      function(x) clue::as.cl_partition(x[, i])
    )

    combined <- clue::cl_consensus(
      partitions,
      method = "HE",
      control = list(
        nruns = 50,
        k = k
      )
    )

    assignments <- apply(
      combined$.Data,
      1,
      function(row) which(row == max(row))
    )

    sce[[paste0("pcaReduce", method, "_cs", k)]] <- assignments
  }

  return(sce)
}

# Internal compatibility fix for the old pcaReduce package.
.pcareduce_mergeM <- function(dat, cl_id, cent, K) {

  P <- matrix(0, K, K)
  a <- utils::combn(seq_len(K), 2)
  v <- c()

  for (l in seq_len(0.5 * (K * K - K))) {

    ind1 <- a[1, l]
    ind2 <- a[2, l]

    x1 <- dat[which(cl_id == ind1), ]
    x2 <- dat[which(cl_id == ind2), ]

    if (is.numeric(x1) && is.null(dim(x1))) {
      C1 <- diag(K - 1)
      n1 <- 1
      mean1 <- x1
    } else {
      C1 <- stats::cov(x1)
      n1 <- nrow(x1)
      mean1 <- colMeans(x1)
    }

    if (is.numeric(x2) && is.null(dim(x2))) {
      C2 <- diag(K - 1)
      n2 <- 1
      mean2 <- x2
    } else {
      C2 <- stats::cov(x2)
      n2 <- nrow(x2)
      mean2 <- colMeans(x2)
    }

    pi1 <- n1 / (n1 + n2)
    pi2 <- n2 / (n1 + n2)

    p2 <- .pcareduce_DMnorm_max(
      rbind(x1, x2),
      pi1 * mean1 + pi2 * mean2,
      pi1 * C1 + pi2 * C2,
      K
    )

    P[ind1, ind2] <- p2
    P[ind2, ind1] <- NA
    v <- c(v, p2)
  }

  diag(P) <- NA

  max_value <- which(P == max(v), arr.ind = TRUE)

  if (nrow(max_value) > 1) {
    to_merge <- max_value[sample(seq_len(nrow(max_value)), 1), ]
  } else {
    to_merge <- max_value
  }

  c(to_merge, pi1, pi2)
}


.pcareduce <- function(D_t, nbt, q, method) {

  Y <- pcaMethods::prep(
    D_t,
    scale = "none",
    center = TRUE
  )

  pca_out <- pcaMethods::pca(
    Y,
    method = "svd",
    center = FALSE,
    nPcs = q
  )

  x <- pca_out@scores

  if (method == "S") {

    Cl_history <- list()

    for (t in seq_len(nbt)) {

      dat <- x
      q <- ncol(dat)
      K <- q + 1

      KM <- stats::kmeans(dat, K)

      cent <- KM$centers
      cl_id <- KM$cluster
      Cl_mat <- c(cl_id)

      Q <- q - 1

      for (i in seq_len(Q)) {

        mrg <- .pcareduce_mergeS(
            dat,
            cl_id,
            cent,
            K
            )

        cl_id[which(cl_id == mrg[[1]][2])] <- mrg[[1]][1]

        cent[mrg[[1]][1], ] <-
          mrg[[2]][1] * cent[mrg[[1]][1], ] +
          mrg[[2]][2] * cent[mrg[[1]][2], ]

        cent <- cent[-mrg[[1]][2], -ncol(cent)]

        a <- unique(cl_id)
        Omega <- seq(1, max(a), 1)
        b <- setdiff(Omega, a)
        N <- length(b)

        if (N > 0) {
          for (ii in seq_len(N)) {
            cl_id[cl_id > b[N + 1 - ii]] <-
              cl_id[cl_id > b[N + 1 - ii]] - 1
          }
        }

        K <- length(unique(cl_id))
        dat <- dat[, -ncol(dat), drop = FALSE]
        Cl_mat <- cbind(Cl_mat, cl_id)
      }

      Cl_history[[t]] <- Cl_mat
    }

    return(Cl_history)
  }

  if (method == "M") {

    Cl_history <- list()

    for (t in seq_len(nbt)) {

      dat <- x
      q <- ncol(dat)
      K <- q + 1

      KM <- stats::kmeans(dat, K)

      cent <- KM$centers
      cl_id <- KM$cluster
      Cl_mat <- c(cl_id)

      Q <- q - 1

      for (i in seq_len(Q)) {

        mrg_ind <- .pcareduce_mergeM(
          dat,
          cl_id,
          cent,
          K
        )

        cl_id[which(cl_id == mrg_ind[2])] <- mrg_ind[1]

        cent[mrg_ind[1], ] <-
          mrg_ind[3] * cent[mrg_ind[1], ] +
          mrg_ind[4] * cent[mrg_ind[2], ]

        cent <- cent[-mrg_ind[2], -ncol(cent)]

        a <- unique(cl_id)
        Omega <- seq(1, max(a), 1)
        b <- setdiff(Omega, a)
        N <- length(b)

        if (N > 0) {
          for (ii in seq_len(N)) {
            cl_id[cl_id > b[N + 1 - ii]] <-
              cl_id[cl_id > b[N + 1 - ii]] - 1
          }
        }

        K <- length(unique(cl_id))
        dat <- dat[, -ncol(dat), drop = FALSE]
        Cl_mat <- cbind(Cl_mat, cl_id)
      }

      Cl_history[[t]] <- Cl_mat
    }

    return(Cl_history)
  }

  stop("`method` must be either 'M' or 'S'.")
}

.pcareduce_DMnorm_max <- function(x, m, C, K) {

  # Force initial symmetry
  C <- (C + t(C)) / 2

  eig <- eigen(C, symmetric = TRUE)

  # Ensure strictly positive eigenvalues
  eig$values[eig$values < 1e-8] <- 1e-8

  C <- eig$vectors %*%
    diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)

  # Reconstruction can introduce tiny numerical asymmetry
  C <- (C + t(C)) / 2

  # Retain original pcaReduce regularization
  C <- C + 0.05 * diag(nrow(C))

  # Be absolutely explicit for mnormt
  C <- (C + t(C)) / 2

  sum(
    mnormt::dmnorm(
      x,
      m,
      C,
      log = TRUE
    )
  )
}

.pcareduce_DMnorm <- function(x, m, C, K) {

  # Force initial symmetry
  C <- (C + t(C)) / 2

  eig <- eigen(C, symmetric = TRUE)

  # Ensure strictly positive eigenvalues
  eig$values[eig$values < 1e-8] <- 1e-8

  C <- eig$vectors %*%
    diag(eig$values, nrow = length(eig$values)) %*%
    t(eig$vectors)

  # Reconstruction can introduce tiny numerical asymmetry
  C <- (C + t(C)) / 2

  # Retain original pcaReduce regularization
  C <- C + 0.05 * diag(nrow(C))

  # Be absolutely explicit for mnormt
  C <- (C + t(C)) / 2

  sum(
    mnormt::dmnorm(
      x,
      m,
      C,
      log = TRUE
    )
  )
}

.pcareduce_mergeS <- function(dat, cl_id, cent, K) {

  P <- matrix(0, K, K)
  a <- utils::combn(seq_len(K), 2)
  v <- c()

  for (l in seq_len(0.5 * (K * K - K))) {

    ind1 <- a[1, l]
    ind2 <- a[2, l]

    x1 <- dat[which(cl_id == ind1), ]
    x2 <- dat[which(cl_id == ind2), ]

    if (is.numeric(x1) && is.null(dim(x1))) {
      C1 <- diag(K - 1)
      n1 <- 1
      mean1 <- x1
    } else {
      C1 <- stats::cov(x1)
      n1 <- nrow(x1)
      mean1 <- colMeans(x1)
    }

    if (is.numeric(x2) && is.null(dim(x2))) {
      C2 <- diag(K - 1)
      n2 <- 1
      mean2 <- x2
    } else {
      C2 <- stats::cov(x2)
      n2 <- nrow(x2)
      mean2 <- colMeans(x2)
    }

    pi1 <- n1 / (n1 + n2)
    pi2 <- n2 / (n1 + n2)

    p2 <- .pcareduce_DMnorm(
      rbind(x1, x2),
      pi1 * mean1 + pi2 * mean2,
      pi1 * C1 + pi2 * C2,
      K
    )

    P[ind1, ind2] <- p2
    P[ind2, ind1] <- NA
    v <- c(v, p2)
  }

  diag(P) <- NA

  a <- t(a)
  p <- P[a]

  p <- exp(p - max(p))

  ii <- sample(
    x = seq_len(nrow(a)),
    size = 1,
    prob = p
  )

  to_merge <- a[ii, ]

  list(
    to_merge,
    c(pi1, pi2),
    p,
    ii
  )
}