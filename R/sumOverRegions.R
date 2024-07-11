#' Sum number of modified and total reads for each region.
#'
#' @description
#' This function takes a \code{RangedSummarizedExperiment} with modification
#' counts (collapsed over reads, such as created by \code{\link{readBedMethyl}})
#' and sums the counts for all features that overlap provided regions.
#'
#' @param se A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     object with assays \code{"Nmod"} and \code{"Nvalid"}, typically
#'     returned by \code{\link{readBedMethyl}}.
#' @param regions A \code{\link[GenomicRegions]{GRanges}} object with target
#'     regions.
#' @param keepZero Logical scalar. If \code{FALSE} (the default), only elements
#'     from \code{regions} that contain at least one feature from \code{se} will
#'     be returned. If \code{TRUE}, all elements of \code{regions} will be
#'     returned, potentially with zero counts.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#' @param ... Additional parameters for \code{\link[scuttle]{aggregateAcrossFeatures}},
#'     such as \code{BPPARAM} to run the summing in parallel.
#'
#' @return A \code{\link[SummarizedExepriment]{RangedSummarizedExperiment}}
#'     with up to \code{length(regions)} rows (exactly \code{length(regions)}
#'     rows if \code{keepZero = TRUE}) and \code{ncol(se)} columns.
#'     \code{colData(se)}, but not \code{rowData(se)} will be preserved.
#'
#' @author Michael Stadler
#'
#' @examples
#' # example bedMethyl file
#' bmfile <- system.file("extdata", "modkit_pileup_1.bed.gz", package = "footprintR")
#'
#' # read into a RangedSummarizedExperiment
#' se <- readBedMethyl(bmfile)
#'
#' # collaps it to a single region of interest
#' regions <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(start = 6940000, end = 7000000, names = "a"))
#' sumOverRegions(se, regions)
#'
#' @seealso \code{\link[scuttle]{aggregateAcrossFeatures}} that is used to
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays colData
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom scuttle aggregateAcrossFeatures
#'
#' @export
sumOverRegions <- function(se, regions, keepZero = FALSE, verbose = FALSE, ...) {
    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertVector(x = regions, type = "GRanges")
    .assertScalar(x = keepZero, type = "logical")
    .assertScalar(x = verbose, type = "logical")

    # find overlaps of features with regions
    if (verbose) {
        message("finding overlaps between features and regions")
    }
    ov <- findOverlaps(query = se, subject = regions)
    if ((n_drop <- sum(!seq.int(nrow(se)) %in% queryHits(ov))) > 0) {
        warning("dropping ", n_drop, " of ", nrow(se), " positions (",
                round(n_drop * 100 / nrow(se), 3), "%) that do not overlap any tile")
    }

    # aggregate counts
    if (verbose) {
        message("aggregating counts")
    }
    u <- unique(subjectHits(ov))
    tmp <- scuttle::aggregateAcrossFeatures(x = se[queryHits(ov), ],
                                            ids = factor(subjectHits(ov), levels = u),
                                            use.assay.type = c("Nmod", "Nvalid"),
                                            ...)
    rownames(tmp) <- NULL

    # keep zeros if requested
    if (keepZero) {
        if (verbose) {
            "filling in zero count regions"
        }
        assayList <- lapply(assays(tmp), function(x) {
            xx <- matrix(data = 0, nrow = length(regions), ncol = ncol(se),
                         dimnames = list(NULL, colnames(se)))
            xx[u, ] <- x
            xx
        })
        rrng <- regions
    } else {
        assayList <- assays(tmp)
        rrng <- regions[u]
    }

    # generate new RangedSummarizedExperiment
    if (verbose) {
        message("generating new RangedSummarizedExperiment")
    }
    seC <- SummarizedExperiment(assays = assayList,
                                rowRanges = rrng,
                                colData = colData(se))
    if (!is.null(names(regions))) {
        rownames(seC) <- names(regions)[u]
    }

    # return results
    return(seC)
}
