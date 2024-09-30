#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowRanges assayNames assay
#'
.pruneAmbiguousStrandPositions <- function(se, assay.type = "Nvalid") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = rowRanges(se), type = "GPos")
    .assertScalar(x = assay.type, type = "character",
                  validValues = assayNames(se))
    .assertVector(x = rownames(se), type = "character")

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
    }

    se
}
