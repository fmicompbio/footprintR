#' Get names of assays containing read-level data
#'
#' @keywords internal
#' @noRd
#'
#' @param se A \code{SummarizedExperiment} object.
#'
#' @author Charlotte Soneson
#'
#' @return A (possibly empty) character vector with the names of the assays of
#' se containing read-level data.
#'
#' @importFrom SummarizedExperiment assayNames assay
.getReadLevelAssayNames <- function(se) {
    isReadLevel <- vapply(SummarizedExperiment::assayNames(se), function(nm) {
        !is.null(dim(SummarizedExperiment::assay(se, nm)[1, 1]))
    }, FALSE)
    SummarizedExperiment::assayNames(se)[isReadLevel]
}

#' Check internal consistency of SummarizedExperiment object
#'
#' All assays with read-level data must have the same number and order of
#' the reads, which must also agree with the order in \code{se$QC} if that
#' exists. All assays must have the same column names, which must also
#' agree with the column names of the object, and the \code{sample} column
#' in the \code{colData}.
#'
#' @keywords internal
#' @noRd
#'
#' @param se A \code{SummarizedExperiment object}.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Charlotte Soneson
#'
#' @return Silently returns \code{NULL}. If the object is not valid, an error
#' will be raised.
#'
#' @importFrom SummarizedExperiment colData assayNames assay
.checkSEValidity <- function(se, verbose = FALSE) {
    if (verbose) {
        message("Checking consistency of sample names")
    }
    stopifnot("sample" %in% colnames(SummarizedExperiment::colData(se)))

    for (an in SummarizedExperiment::assayNames(se)) {
        stopifnot(colnames(SummarizedExperiment::assay(
            se, an, withDimnames = FALSE)) == se$sample)
    }
    stopifnot(rownames(SummarizedExperiment::colData(se)) == se$sample)
    if ("QC" %in% colnames(SummarizedExperiment::colData(se))) {
        if (!is.null(dim(se$QC[[1]]))) {
            stopifnot(names(se$QC) == se$sample)
        }
    }

    rlAssays <- .getReadLevelAssayNames(se)
    if (length(rlAssays) > 0) {
        if (verbose) {
            message("Read-level assay found")
        }
        ## Choose one assay as the reference to compare to
        refAssay <- rlAssays[1]
        refReads <- lapply(SummarizedExperiment::assay(se, refAssay),
                           colnames)
        for (an in setdiff(rlAssays, refAssay)) {
            if (verbose) {
                message("Comparing ", refAssay, " and ", an)
            }
            for (sn in se$sample) {
                if (!all(colnames(SummarizedExperiment::assay(se, an)[[sn]]) ==
                         refReads[[sn]])) {
                    stop("Mismatching reads for assays ", refAssay, " and ",
                         an, ", sample ", sn)
                }
            }
        }
        if ("QC" %in% colnames(SummarizedExperiment::colData(se))) {
            if (!is.null(dim(se$QC[[1]]))) {
                if (verbose) {
                    message("QC information found, checking consistency")
                }
                for (sn in se$sample) {
                    if (!all(rownames(se$QC[[sn]]) == refReads[[sn]])) {
                        stop("Mismatching reads for assay ", refAssay, " and ",
                             "sample QC data, sample ", sn)
                    }
                }
            }
        }
    }
}
