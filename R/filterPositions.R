#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment assay assayNames
#'
.filterPositionsByCoverage <- function(se, assayName = "Nvalid", minCov = 1,
                                       minNbrSamples = NULL) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = minCov, type = "numeric")
    .assertScalar(x = minNbrSamples, type = "numeric", allowNULL = TRUE)

    # If assayName is a read-level assay, first calculate the number of
    # non-NA values in each row
    if (assayName %in% .getReadLevelAssayNames(se)) {
        mat <- assay(flattenReadLevelAssay(se, assayName = assayName,
                                           statistics = "Nvalid", keepReads = FALSE,
                                           verbose = FALSE),
                     "Nvalid")
    } else {
        mat <- assay(se, assayName)
    }

    if (is.null(minNbrSamples)) {
        ## use the total coverage (sum across all samples)
        keep <- which(rowSums(mat) >= minCov)
    } else {
        keep <- which(rowSums(mat >= minCov) >= minNbrSamples)
    }

    se[keep, ]
}

#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowData
#' @importFrom Biostrings vcountPattern
#' @importFrom cli cli_abort
#'
.keepPositionsBySequenceContext <- function(se, sequenceContext = NULL) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = sequenceContext, type = "character", allowNULL = TRUE)

    if (!is.null(sequenceContext)) {
        if (is.null(rowData(se)$sequenceContext)) {
            cli_abort("No sequence context found in {.code rowData(se)$sequenceContext}")
        }
        .assertVector(x = rowData(se)$sequenceContext,
                      type = "DNAStringSet")
        nmatch <- Reduce("+", lapply(sequenceContext, function(pat) {
            vcountPattern(pat, rowData(se)$sequenceContext, fixed = "subject")
        }), init = rep(0, nrow(se)))
        se <- se[nmatch > 0, ]
    }
    se
}

#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment assay
#' @importFrom SparseArray rowSums is_nonna
#'
.removeAllNAPositions <- function(se, assayName = "mod_prob") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))

    # Get requested assay and convert to a single NaMatrix
    mat <- as.matrix(assay(se, assayName))

    # Find positions to keep and subset se
    keep <- which(rowSums(is_nonna(mat)) > 0)
    se[keep, ]
}

#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowRanges assayNames assay
#' @importFrom BiocGenerics pos
#' @importFrom GenomeInfoDb seqnames
#' @importFrom cli cli_warn
#'
.pruneAmbiguousStrandPositions <- function(se, assayName = "Nvalid",
                                           verbose = FALSE) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = rowRanges(se), type = "GPos")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertVector(x = rownames(se), type = "character")
    .assertScalar(x = verbose, type = "logical")

    # Group positions by chromosome and position
    pGroup <- split(x = rownames(se),
                    f = paste0(seqnames(rowRanges(se)),
                               ":", pos(rowRanges(se))))
    pGroup <- pGroup[lengths(pGroup) > 1]

    # For all groups of >1 row, find the one with lowest total count and
    # record the row name for later removal
    tmpmat <- as.matrix(assay(se, assayName)[
        unlist(pGroup, use.names = FALSE), ])
    if (assayName %in% .getReadLevelAssayNames(se)) {
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
        .message(
            paste0("{length(posToRemove)} row{?s} removed to ensure that each ",
                   "genomic position is represented by at most one row"))
    } else {
        .message("No genomic positions represented by multiple rows found")
    }

    # Check that removal worked as expected, i.e. that all remaining
    # positions are only present once
    # Leave it as a warning so that the user has a chance to inspect the
    # output
    if (length(unique(paste0(seqnames(rowRanges(se)), ":",
                             pos(rowRanges(se))))) != nrow(se)) {
        # nocov start
        cli_warn(paste0(
            "Pruning of ambiguous positions failed - the object still ",
            "contains positions represented by multiple rows"))
        # nocov end
    }

    se
}

#' Filter positions
#'
#' Filter positions based on any combination of sequence context,
#' coverage, repetition of the same position, and presence of non-NA values.
#' Filters are applied in the order they are specified to the \code{filters}
#' argument. Any filter type can be repeated an arbitrary number of times.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param filters A character vector. All values must be one of
#'     \code{"sequenceContext"}, \code{"coverage"}, \code{"repeated.positions"}
#'     and \code{"all.na"}. Filters are applied in the order specified by
#'     this vector.
#' @param sequenceContext A character vector with sequence contexts to
#'     retain. To apply this filter, the \code{"sequenceContext"} column must
#'     be present in \code{rowData(se)} (see \code{addSeqContext}).
#' @param assayNameCov A character scalar indicating the assay to use to
#'     define the coverage. If this is a read-level assay, coverage is first
#'     calculated using \code{flattenReadLevelAssay(..., statistics = "Nvalid")}.
#' @param minCov A numeric scalar indicating the lowest acceptable
#'     coverage in order to keep a position.
#' @param minNbrSamples A numeric scalar, or \code{NULL}. If \code{NULL}
#'     (default), the row sum of \code{assayNameCov} (i.e., the total
#'     coverage across all samples) is used for the coverage filtering. If
#'     not \code{NULL}, a position is required to have at least \code{minCov}
#'     coverage in at least \code{minNbrSamples} to be retained.
#' @param assayNameAmbig A character scalar indicating the assay to use to
#'     decide which row to retain if multiple rows represent the same
#'     genomic position (on different strands). The row with the largest row
#'     sum in this assay is retained.
#' @param assayNameNA A character scalar indicating the assay to use as the
#'     basis for filtering out positions with NA values across all reads.
#'     This should be a read level assay.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A filtered \code{SummarizedExperiment}.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- flattenReadLevelAssay(se)
#' se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = reffile)
#' sefilt <- filterPositions(se, c("sequenceContext", "coverage", "all.na"),
#'                           minCov = 5, sequenceContext = "TAG")
#'
#' @importFrom SparseArray colSums is_nonna
#' @importFrom SummarizedExperiment assay
filterPositions <- function(se,
                            filters = c("sequenceContext", "coverage",
                                        "all.na"),
                            sequenceContext = NULL,
                            assayNameCov = "Nvalid",
                            minCov = 1,
                            minNbrSamples = NULL,
                            assayNameAmbig = "Nvalid",
                            assayNameNA = "mod_prob") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = filters, type = "character",
                  validValues = c("sequenceContext", "coverage",
                                  "repeated.positions", "all.na"))
    .assertScalar(x = assayNameNA, type = "character",
                  validValues = .getReadLevelAssayNames(se),
                  allowNULL = TRUE)

    for (f in filters) {
        if (f == "sequenceContext") {
            se <- .keepPositionsBySequenceContext(
                se, sequenceContext = sequenceContext
            )
        } else if (f == "coverage") {
            se <- .filterPositionsByCoverage(
                se, assayName = assayNameCov, minCov = minCov,
                minNbrSamples = minNbrSamples
            )
        } else if (f == "repeated.positions") {
            se <- .pruneAmbiguousStrandPositions(
                se, assayName = assayNameAmbig
            )
        } else if (f == "all.na") {
            se <- .removeAllNAPositions(
                se, assayName = assayNameNA
            )
        }
    }

    ## Remove reads that are NA in all retained positions
    if (!is.null(assayNameNA)) {
        readsToKeep <- lapply(assay(se, assayNameNA),
                              function(x) {
                                  which(colSums(is_nonna(x)) > 0)
                              })
        se <- subsetReads(se, readsToKeep)
    }

    se
}
