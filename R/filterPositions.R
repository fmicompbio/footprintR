#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment assay assayNames
#' 
.filterPositionsByCoverage <- function(se, assay.type = "Nvalid", min.cov = 1) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assay.type, type = "character", 
                  validValues = assayNames(se))
    .assertScalar(x = min.cov, type = "numeric")
    
    # If assay.type is a read-level assay, first calculate the number of 
    # non-NA values in each row
    if (assay.type %in% .getReadLevelAssayNames(se)) {
        mat <- assay(addReadsSummary(se, assay.type = assay.type, 
                                     statistics = "Nvalid", keep.reads = FALSE,
                                     verbose = FALSE), "Nvalid")
    } else {
        mat <- assay(se, assay.type)
    }
    
    keep <- rownames(mat)[which(rowSums(mat) >= min.cov)]
    se[keep, ]
}

#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowData
#' @importFrom Biostrings vcountPattern
#' 
.keepPositionsBySequenceContext <- function(se, sequence.context = NULL) {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = sequence.context, type = "character", allowNULL = TRUE)
    
    if (!is.null(sequence.context)) {
        if (is.null(rowData(se)$sequence.context)) {
            stop("No sequence context found in `rowData(se)$sequence.context`")
        }
        nmatch <- Reduce("+", lapply(sequence.context, function(pat) {
            vcountPattern(pat,
                          rowData(se)$sequence.context,
                          fixed = FALSE)
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
.removeAllNAPositions <- function(se, assay.type = "mod_prob") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assay.type, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    
    # Get requested assay and convert to a single NaMatrix
    mat <- as.matrix(assay(se, assay.type))
    
    # Find positions to keep and subset se
    keep <- which(rowSums(is_nonna(mat)) > 0)
    se[keep, ]
} 

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

#' Filter positions
#' 
#' Filter positions based on any combination of sequence context, 
#' coverage, repetition of the same position, and presence of non-NA values. 
#' Filters are applied in the order they are specified to the \code{filters}
#' argument. Any filter type can be repeated an arbitrary number of times.
#' 
#' @param se A \code{SummarizedExperiment} object.
#' @param filters A character vector. All values must be one of 
#'     \code{"sequence.context"}, \code{"coverage"}, \code{"repeated.positions"}
#'     and \code{"all.na"}. Filters are applied in the order specified by 
#'     this vector. 
#' @param sequence.context A character vector with sequence contexts to 
#'     retain. To apply this filter, the \code{"sequence.context"} column must 
#'     be present in \code{rowData(se)} (see \code{addSeqContext}).
#' @param assay.type.cov A character scalar indicating the assay to use to 
#'     define the coverage. If this is a read-level assay, coverage is first 
#'     calculated using \code{addReadsSummary(..., statistics = "Nvalid")}. 
#' @param min.cov A numeric scalar indicating the lowest acceptable 
#'     coverage in order to keep a position.
#' @param assay.type.ambig A character scalar indicating the assay to use to 
#'     decide which row to retain if multiple rows represent the same 
#'     genomic position (on different strands). The row with the largest row 
#'     sum in this assay is retained. 
#' @param assay.type.na A character scalar indicating the assay to use as the 
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
#'                  modbase = "a", verbose = FALSE)
#' se <- addReadsSummary(se)
#' se <- addSeqContext(se, sequence.context.width = 3, sequence.reference = reffile)
#' sefilt <- filterPositions(se, c("sequence.context", "coverage", "all.na"),
#'                           min.cov = 5, sequence.context = "TAG")
filterPositions <- function(se, 
                            filters = c("sequence.context", "coverage",
                                        "all.na"),
                            sequence.context = NULL,
                            assay.type.cov = "Nvalid",
                            min.cov = 1,
                            assay.type.ambig = "Nvalid",
                            assay.type.na = "mod_prob") {
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = filters, type = "character", 
                  validValues = c("sequence.context", "coverage",
                                  "repeated.positions", "all.na"))

    for (f in filters) {
        if (f == "sequence.context") {
            se <- .keepPositionsBySequenceContext(
                se, sequence.context = sequence.context
            )
        } else if (f == "coverage") {
            se <- .filterPositionsByCoverage(
                se, assay.type = assay.type.cov, min.cov = min.cov
            )
        } else if (f == "repeated.positions") {
            se <- .pruneAmbiguousStrandPositions(
                se, assay.type = assay.type.ambig
            )
        } else if (f == "all.na") {
            se <- .removeAllNAPositions(
                se, assay.type = assay.type.na
            )
        }
    }
    se
}
