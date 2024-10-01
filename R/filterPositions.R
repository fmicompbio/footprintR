#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowRanges assayNames assay
#'
.pruneAmbiguousStrandPositions <- function(se, assay.type = "Nvalid",
                                           verbose = FALSE) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = rowRanges(se), type = "GPos")
    .assertScalar(x = assay.type, type = "character",
                  validValues = assayNames(se))
    .assertVector(x = rownames(se), type = "character")
    .assertScalar(x = verbose, type = "logical")

    # Group positions by chromosome and position
    pGroup <- split(rownames(se), f = paste0(seqnames(rowRanges(se)), ":",
                                             pos(rowRanges(se))))
    pGroup <- pGroup[lengths(pGroup) > 1]

    # For all groups of >1 row, find the one with lowest total count and
    # record the row name for later removal
    tmpmat <- as.matrix(assay(se, assay.type)[unlist(pGroup, use.names = FALSE), ])
    if (assay.type %in% .getReadLevelAssayNames(se)) {
        rs <- rowSums(tmpmat >= 0, na.rm = TRUE)
    } else {
        rs <- rowSums(tmpmat, na.rm = TRUE)
    }
    posToRemove <- unlist(lapply(pGroup, function(pg) {
        pg[-which.max(rs[pg])]
    }))

    # Remove the recorded positions
    if (length(posToRemove) > 0) {
        se <- se[!rownames(se) %in% posToRemove, ]
        if (verbose) {
            message(length(posToRemove), " rows removed to ensure that each ",
                    "genomic position is represented by at most one row")
        }
    } else {
        if (verbose) {
            message("No genomic positions represented by multiple rows found")
        }
    }

    # Check that removal worked as expected, i.e. that all remaining
    # positions are only present once
    # Leave it as a warning so that the user has a chance to inspect the
    # output
    if (length(unique(paste0(seqnames(rowRanges(se)), ":",
                             pos(rowRanges(se))))) != nrow(se)) {
        # nocov start
        warning("Pruning of ambiguous positions failed - the object still ",
                "contains positions represented by multiple rows")
        # nocov end
    }

    se
}
