#' Remove all-NA reads
#'
#' Remove reads (columns) that contain only NA values from a read-level
#' assay (\code{DataFrame} of \code{NaArray}s).
#'
#' @param x A \code{\link[S4Vectors]{DataFrame}} with
#'     \code{\link[SparseArray]{NaArray}} objects in its columns (typically
#'     a read-level assay as returned by \code{\link{readModBam}}).
#' @param prune A logical scalar. If \code{TRUE} (the default), samples
#'     (columns of the \code{DataFrame}) for which the NA-read filtering retains none of
#'     the reads will be completely removed. If \code{FALSE},
#'     such samples are retained as a zero-column \code{NAMatrix}).
#'
#' @returns A \code{DataFrame} of \code{NaArray}s with all columns containing at
#'     least one non-\code{NA} value.
#'
#' @importFrom S4Vectors endoapply
#' @importFrom SparseArray colSums is_nonna
#' @importFrom BiocGenerics rownames rownames<-
#'
#' @noRd
#' @keywords internal
.removeAllNAReads <- function(x, prune = TRUE) {
    .assertScalar(x = prune, type = "logical")

    rnms <- rownames(x)
    x <- endoapply(x, function(y) {
        if (!is.null(dim(y))) {
            y <- y[, colSums(is_nonna(y)) > 0, drop = FALSE]
        }
        y
    })
    rownames(x) <- rnms

    if (prune) {
        x <- x[, unlist(lapply(x, ncol), use.names = FALSE) > 0, drop = FALSE]
    }

    return(x)
}

#' Filter reads
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param assayName A character scalar providing the name of a read-level
#'     assay in \code{se}. This assay will be used to extract read names, as
#'     well as to filter out any read that is not overlapping any of the
#'     positions in the object.
#' @param readInfoCol A character scalar providing the name of the column in
#'     \code{colData} that contains read info. Can be \code{NULL} if no such
#'     column exists.
#' @param qcCol A character scalar providing the name of the column in
#'     \code{colData} that contains quality metrics (calculated by
#'     \code{calcReadStats}). Can be \code{NULL} if no such column exists.
#' @param minQscore A numeric scalar representing the smallest acceptable
#'     read-level Qscore. Reads with Qscore below this value will be filtered
#'     out.
#' @param maxEntropy A numeric scalar representing the largest acceptable
#'     read-level entropy. Reads with entropy above this value will be filtered
#'     out.
#' @param maxFracLowConf A numeric scalar representing the maximally acceptable
#'     fraction of low-confidence modified base calls in a read. Reads with
#'     a fraction of low confidence calls greater than this value will be
#'     filtered out.
#' @param minReadLength A numeric scalar representing the smallest acceptable
#'     read length. Reads that are shorter than this value will be filtered
#'     out.
#' @param minAlignedLength A numeric scalar representing the smallest acceptable
#'     aligned length. Reads with aligned length shorter than this value will
#'     be filtered out.
#' @param minAlignedFraction A numeric scalar representing the smallest
#'     acceptable aligned fraction of a read. Reads where the aligned fraction
#'     is smaller than this value will be filtered out.
#' @param minCoveredFraction A numeric scalar giving the minimal fraction of
#'     \code{region} that a read alignment needs to cover to be retained.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region to be used for the \code{minCoveredFraction} filter.
#'     Alternatively, the region can be specified as a character scalar (e.g.
#'     "chr1:1200-1300") that can be coerced into a \code{GRanges} object. If
#'     \code{NULL} (the default), all reads are retained.
#' @param prune A logical scalar. If \code{TRUE} (the default), samples for
#'     which the filtering retains none of the reads will be completely removed
#'     from the returned \code{SummarizedExperiment} (also from \code{colData}
#'     and from assays that do not store read-level data). If \code{FALSE},
#'     such samples are retained (in the assays with read-level data as a
#'     zero-column \code{SparseMatrix}).
#' @param onlyStats A logical scalar. If \code{FALSE} (the default), the
#'     \code{SummarizedExperiment} object will be filtered according to the
#'     provided thresholds. If \code{TRUE}, the filter statistics are calculated
#'     and returned, but the object is not subset.
#' @param removeAllNApos A logical scalar. If \code{TRUE}, remove all positions
#'     for which all values in \code{assayName} are \code{NA} after the read
#'     filtering.
#'
#' @author Charlotte Soneson, Michael Stadler
#' @export
#'
#' @returns If \code{onlyStats} is \code{FALSE}, a filtered
#' \code{SummarizedExperiment} object. The metadata of this
#' object contains a slot named \code{filteredOutReads}, which tabulate all
#' reads that are filtered out, together with the reason(s) for exclusion.
#' If \code{onlyStats} is \code{TRUE}, only this table is returned.
#'
#' @examples
#' library(SummarizedExperiment)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfile = modbamfile, regions = "chr1:6920000-6995000",
#'                  modbase = "a", verbose = TRUE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- addReadStats(se, name = "QC",
#'                    BPPARAM = BiocParallel::SerialParam())
#'
#' ## Filter se
#' ## ... by quality score and aligned length
#' sefilt <- filterReads(se, minQscore = 14, minAlignedLength = 10000)
#' ## ... by covered fraction of the specified region
#' sefilt <- filterReads(se, minCoveredFraction = 1, region = "chr1:6930000-6935000")
#'
#' ## Only calculate filter stats
#' filtstats <- filterReads(se, minQscore = 14, minAlignedLength = 10000,
#'                          onlyStats = TRUE)
#' filtstats
#'
#' ## Visualize filter stats in UpSet plot, e.g. with ComplexUpset
#' if (require(ComplexUpset)) {
#'     ComplexUpset::upset(as.data.frame(filtstats$s1),
#'                         intersect = colnames(filtstats$s1))
#' }
#'
#' @importFrom SparseArray SVT_SparseArray rowSums colSums is_nonna nnawhich
#' @importFrom SummarizedExperiment colData assay rowRanges
#' @importFrom BiocGenerics pos
#' @importFrom GenomicRanges GRanges pintersect
#' @importFrom IRanges IRanges
#'
filterReads <- function(se, assayName = "mod_prob",
                        readInfoCol = "readInfo", qcCol = "QC",
                        minQscore = 0, maxEntropy = Inf,
                        maxFracLowConf = 1, minReadLength = 0,
                        minAlignedLength = 0, minAlignedFraction = 0,
                        minCoveredFraction = 0, region = NULL,
                        prune = TRUE, onlyStats = FALSE,
                        removeAllNApos = TRUE) {
    ## Input checks
    .assertVector(x = se, type = "SummarizedExperiment")
    .checkSEValidity(se)
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = readInfoCol, type = "character", allowNULL = TRUE,
                  validValues = colnames(colData(se)))
    .assertScalar(x = qcCol, type = "character", allowNULL = TRUE,
                  validValues = colnames(colData(se)))
    .assertScalar(x = minQscore, type = "numeric")
    .assertScalar(x = maxEntropy, type = "numeric")
    .assertScalar(x = maxFracLowConf, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minReadLength, type = "numeric")
    .assertScalar(x = minAlignedLength, type = "numeric")
    .assertScalar(x = minAlignedFraction, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minCoveredFraction, type = "numeric", rngIncl = c(0, 1))
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    }
    .assertScalar(x = region, type = "GRanges", allowNULL = TRUE)
    .assertScalar(x = prune, type = "logical")
    .assertScalar(x = onlyStats, type = "logical")
    .assertScalar(x = removeAllNApos, type = "logical")

    ## Initialize sparse logical array for each sample, which will be TRUE
    ## for reads that are filtered out with respect to the different criteria
    ## Remark: Could move this to a global constant
    filterNames <- c("Qscore", "Entropy", "FracLowConf", "ReadLength",
                     "AlignedLength", "AlignedFraction", "CoveredFraction",
                     "AllNA")
    readsToRemove <- lapply(
        structure(colnames(se), names = colnames(se)),
        function(nm) {
            SVT_SparseArray(
                dim = c(ncol(assay(se, assayName)[[nm]]),
                        length(filterNames)),
                dimnames = list(colnames(assay(se, assayName)[[nm]]),
                                filterNames),
                type = "logical"
        )}
    )

    for (nm in colnames(se)) {
        ## Extract read info and QC
        if (!is.null(readInfoCol) && !is.null(se[[readInfoCol]]) &&
            (nm %in% names(se[[readInfoCol]]))) {
            ri <- se[[readInfoCol]][[nm]]
        } else {
            ri <- NULL
        }
        if (!is.null(qcCol) && !is.null(se[[qcCol]]) &&
            (nm %in% names(se[[qcCol]]))) {
            qc <- se[[qcCol]][[nm]]
        } else {
            qc <- NULL
        }

        ## Quality score
        if (!is.null(ri) && "qscore" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$qscore < minQscore),
                                "Qscore"] <- TRUE
        }

        ## KS entropy
        if (!is.null(qc) && "SEntrModProb" %in% colnames(qc)) {
            readsToRemove[[nm]][which(qc$SEntrModProb > maxEntropy),
                                "Entropy"] <- TRUE
        }

        ## Fraction of low-confidence modification calls
        if (!is.null(qc) && "FracLowConf" %in% colnames(qc)) {
            readsToRemove[[nm]][which(qc$FracLowConf > maxFracLowConf),
                                "FracLowConf"] <- TRUE
        }

        ## Read length
        if (!is.null(ri) && "read_length" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$read_length < minReadLength),
                                "ReadLength"] <- TRUE
        }

        ## Aligned length
        if (!is.null(ri) && "aligned_length" %in% colnames(ri)) {
            readsToRemove[[nm]][ which(ri$aligned_length < minAlignedLength),
                                 "AlignedLength"] <- TRUE
        }

        ## Aligned fraction
        if (!is.null(ri) && "aligned_fraction" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$aligned_fraction < minAlignedFraction),
                                "AlignedFraction"] <- TRUE
        }

        ## Region coverage fraction
        if (!is.null(region) && minCoveredFraction > 0) {
            tmp <- nnawhich(assay(se, assayName)[[nm]], arr.ind = TRUE)
            tmpp <- split(pos(rowRanges(se))[tmp[, 1]], tmp[, 2])
            tmpgr <- GRanges(
                seqnames = seqnames(rowRanges(se))[1],
                ranges = IRanges(start = unlist(lapply(tmpp, min),
                                                use.names = FALSE),
                                 end = unlist(lapply(tmpp, max),
                                              use.names = FALSE)))
            cvgFrac <- width(
                pintersect(tmpgr, region, ignore.strand = TRUE,
                           drop.nohit.ranges = FALSE)) / width(region)
            readsToRemove[[nm]][
                as.numeric(names(tmpp))[which(cvgFrac < minCoveredFraction)],
                "CoveredFraction"] <- TRUE
        }

        ## NA in all positions
        readsToRemove[[nm]][colnames(
            assay(se, assayName)[[nm]])[colSums(
                is_nonna(assay(se, assayName)[[nm]])) == 0],
            "AllNA"] <- TRUE
    }

    ## Subset
    readsToRemove <- lapply(readsToRemove, function(rr) {
        rr[rowSums(rr, na.rm = TRUE) > 0, ]
    })
    if (onlyStats) {
        return(readsToRemove)
    } else {
        sesub <- subsetReads(se = se, reads = lapply(readsToRemove, rownames),
                             prune = prune, invert = TRUE,
                             removeAllNApos = removeAllNApos,
                             assayNameNA = assayName)
        metadata(sesub)$filteredOutReads <- readsToRemove

        return(sesub)
    }
}
