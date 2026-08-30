#' Example Single-Cell Expression Data
#'
#' A reduced example expression matrix derived from the single-cell RNA-seq
#' dataset reported by Biase et al. The original dataset contains cell
#' identities established in that study and is available from the NCBI Gene
#' Expression Omnibus under accession GSE57249.
#'
#' For use in sclValid examples, 1000 rows were randomly sampled from the
#' original expression matrix while retaining all 49 observations. Rows
#' represent features and columns represent cells.
#'
#' @format A numeric matrix with 1000 rows and 49 columns.
#'
#' @source
#' Biase et al. single-cell RNA-seq dataset, NCBI Gene Expression Omnibus,
#' accession GSE57249.
"example_data"


#' Example Cell Labels
#'
#' Cell identity labels corresponding to the 49 columns of
#' \code{example_data}. The labels are derived from the cell identities
#' established in the Biase et al. study associated with GEO accession
#' GSE57249.
#'
#' For use in sclValid, the original cell groups were relabeled as `"A"`,
#' `"B"`, and `"C"`.
#'
#' @format A character vector of length 49.
#'
#' @source
#' Biase et al. single-cell RNA-seq dataset, NCBI Gene Expression Omnibus,
#' accession GSE57249.
"example_labels"