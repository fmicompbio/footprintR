# global vector with default read stats (that will be calculated if
# stats = NULL in calcReadStats)
# exclude "SEntrModProb"
defaultReadStats <- c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "Lag1DModProb", "ACModProb", "PACModProb")

#' Calculate or add summary statistics for read-level base modification data
#'
#' @description
#' \code{calcReadStats} calculates various per-read summary statistics on
#' modification probabilities or calls from a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with genomic
#' positions in rows and samples in columns. \code{addReadStats} adds them to
#' the \code{colData} under \code{name}. See details for more information on the
#' statistics that are calculated.
#'
#' @param se A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     object with assay \code{assayName} typically returned by
#'     \code{\link{readModkitExtract}} or \code{\link{readModBam}}.
#' @param assayName A character scalar specifying the assay of \code{se}
#'     containing the read-level data to be summarized. Typically, this assay
#'     contains modification probabilities.
#' @param stats Character vector specifying which statistics to calculate.
#'     When set to \code{NULL} all available statistics except the sample
#'     entropy are calculated. See details for available read statistics.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object limiting the
#'     positions included in the calculations to the ones overlapping the
#'     corresponding genomic regions. Alternatively, regions can be
#'     specified as a character vector (e.g. "chr1:1200-1300") that can be
#'     coerced into a \code{GRanges} object.
#' @param sequenceContext A character vector with sequence context(s)
#'     to include in the calculations. Only positions that match one of the
#'     provided sequence contexts will be included. Sequence contexts can be
#'     provided using IUPAC redundancy codes. The sequence contexts of modified
#'     bases are obtained from \code{rowData(se)$sequenceContext} and thus
#'     requires that \code{se} contains the appropriate information, for example
#'     by setting the \code{sequenceContextWidth} and \code{sequenceReference}
#'     arguments of \code{\link{readModkitExtract}} when it was generated,
#'     or by adding it using \code{\link{addSeqContext}}.
#' @param minNobsPpos A numeric scalar value >=1 indicating the minimum
#'     coverage on individual positions for them to be included in the
#'     calculations. In high coverage data this is an effective filter for
#'     removing spurious modbases, typically the result of erroneous
#'     basecalling. The default \code{NULL} sets its value to Q3-0.5*IQR, where
#'     Q3 and IQR are the third quartile and interquartile range of the coverage
#'     distribution estimated from the data in \code{se}.
#' @param minNobsPread A numeric scalar with the minimum number of observed
#'     modifiable bases per read for it to be included in the calculations.
#'     \code{NA} values are returned for the reads that do not pass this
#'     threshold.
#' @param LowConf A numeric scalar with the minimum call confidence below which
#'     calls are considered "low confidence".
#' @param LagRange A numeric vector of two values (minimum and maxium) defining
#'     the range of lags for the calculation of autocorrelation and partial
#'     autocorrelation (see details section).
#' @param name For \code{addReadStats} only: A character scalar specifying the
#'     name to be used to store the result in the
#'     \code{\link[SummarizedExperiment]{colData}} of the output.
#' @param ... For \code{addReadStats} only: Additional arguments passed on to
#'     \code{calcReadStats}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for calculating
#'     read statistics. The default value (\code{\link[BiocParallel]{bpparam}})
#'     will select an appropriate value for the current environment, or the
#'     default parallel backend registered using
#'     \code{\link[BiocParallel]{register}}.
#' @param verbose If \code{TRUE}, report on progress.
#'
#' @details
#' \code{calcReadStats} calculates a collection of location/scatter statistics
#' and information theoretic/signal-processing metrics for the modification
#' probability, confidence or modification call value vectors across individual
#' reads (data in assay \code{assayName}). Only bases matching the criteria
#' given by\code{regions}, \code{sequenceContext}, \code{minNobsPpos} and
#' \code{minNobsPread} are included in the calculations. The values of these
#' filtering parameters are stored in the attribute of the output.
#'
#' \code{stats} selects the summaries to be calculated. Currently available
#' values are:
#' \describe{
#'     \item{MeanModProb}{: Mean modification probability across the read.}
#'     \item{FracMod}{: Fraction of confidently called bases (either modified
#'         or unmodified), that are modified. By default, all bases are
#'         considered confidently called.}
#'     \item{MeanConf}{:  Mean call confidence across the read.}
#'     \item{MeanConfUnm}{: Mean call confidence confined to unmodified bases
#'         (modifiable bases with modification probability < 0.5).}
#'     \item{MeanConfMod}{: Mean call confidence confined to modified bases
#'         (modifiable bases with modification probability >= 0.5).}
#'     \item{FracLowConf}{: Fraction of modifiable bases called with low
#'         confidence (call confidence < \code{LowConf}).}
#'     \item{IQRModProb}{: Interquartile range of modification probabilities
#'         across the read.}
#'     \item{sdModProb}{: Standard deviation of modification probabilities
#'         across the read.}
#'     \item{SEntrModProb}{: Sample entropy of the modification probability
#'         signal. Sample entropy is a metric assessing the complexity of
#'         one-dimensional physiological signals. The higher the sample entropy,
#'         the more irregular, unpredictable and therefore complex the signal.
#'         See [wikipedia:Sample_entropy](https://en.wikipedia.org/wiki/Sample_entropy)
#'         for more details.}
#'     \item{Lag1DModProb}{: Mean Lag1 differences of modification calls,
#'         defined as: \code{mean(Mod[i]-Mod[i-1])}, where \code{Mod} is a
#'         \code{{0,1}} modification call.}
#'     \item{ACModProb}{: Autocorrelation of the modification probability values
#'         for lags in the range \code{LagRange}. This range typically covers
#'         the signal of nucleosome periodicity.}
#'     \item{PACModProb}{: Partial autocorrelation of the modification
#'         probability values for lags in the range \code{LagRange}. This range
#'         typically covers the signal of nucleosome periodicity.}
#'  }
#'
#' @return
#' For \code{calcReadStats}, a \code{SimpleList} object with summary statistics
#' for the samples (columns) in \code{se}.
#'
#' For \code{addReadStats}, a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with summary statistics added to the \code{name} column of
#' \code{\link[SummarizedExperiment]{colData}}.
#'
#' @author Panagiotis Papapasaikas, Charlotte Soneson, Michael Stadler
#' @name calcReadStats
#'
#' @examples
#' # load example data
#' library(SummarizedExperiment)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfile = modbamfile, regions = "chr1:6940000-6955000",
#'            modbase = "a", verbose = TRUE)
#'
#' readStats <- calcReadStats(se)
#' readStats$s1
#'
#' se_withReadStats <- addReadStats(se, name = "QC")
#' se_withReadStats$QC$s1
#' metadata(se_withReadStats$QC$s1)
#'
#' @importFrom S4Vectors metadata make_zero_col_DFrame SimpleList
#' @importFrom SummarizedExperiment assay
#' @importFrom SparseArray rowSums nnawhich nnavals
#' @importFrom stats sd IQR acf pacf na.pass
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocGenerics colnames
#' @importFrom BiocParallel bplapply bpparam
#'
#' @export
calcReadStats <- function(se,
                          assayName = "mod_prob",
                          stats = NULL,
                          regions = NULL,
                          sequenceContext = NULL,
                          minNobsPpos = 0,
                          minNobsPread = 0,
                          LowConf = 0.7,
                          LagRange = c(12, 64),
                          BPPARAM = bpparam(),
                          verbose = FALSE) {
    # define functions to calculate summary statistics
    statFunctions <- list(
        MeanModProb = mean,
        FracMod = function(x, c = 0.5) {
            sum(x >= (0.5 + (c - 0.5))) / sum(abs(0.5 - x) >= (c - 0.5))
        },
        MeanConf = function(x) {
            mean(pmax(x, 1 - x))
        },
        MeanConfUnm = function(x) {
            mean((1 - x)[x < 0.5])
        },
        MeanConfMod = function(x) {
            mean(x[x >= 0.5])
        },
        FracLowConf = function(x, c = LowConf) {
            sum(abs(0.5 - x) < (c - 0.5)) / length(x)
        },
        IQRModProb = IQR,
        sdModProb = sd,
        SEntrModProb = function(x) {
            if (length(x) > 64) {
                sampleEntropy(x, 2L, 0.2)
            } else {
                NA
            }
        },
        Lag1DModProb = function(x) {
            xC <- as.numeric(x >= 0.5)
            mean(abs(diff(xC, lag = 1)))
        },
        ACModProb = function(x, lag.max = max(LagRange),
                             xrange = LagRangeValues) {
            if (length(x) > lag.max) {
                acf(x, na.action = na.pass, lag.max = lag.max,
                           plot = FALSE)$acf[xrange]
            } else {
                rep(0, length(xrange))
            }
        },
        PACModProb = function(x, lag.max = max(LagRange),
                              xrange = LagRangeValues) {
            if (length(x) > lag.max) {
                pacf(x, na.action = na.pass, lag.max = lag.max,
                            plot = FALSE)$acf[xrange]
            } else {
                rep(0, length(xrange))
            }
        }
    )

    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertVector(x = stats, type = "character", allowNULL = TRUE,
                  validValues = names(statFunctions))
    if (is.character(regions)) {
        regions <- as(regions, "GRanges")
    }
    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)
    .assertVector(x = sequenceContext, type = "character", allowNULL = TRUE)
    .assertScalar(x = minNobsPpos, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minNobsPread, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = LowConf, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = LagRange, type = "vector", rngIncl = c(1, 256), len = 2)
    LagRangeValues <- seq(LagRange[1], LagRange[2])
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")

    # Subset se by region
    if (!is.null(regions)) {
        se <- subsetByOverlaps(x = se, ranges = regions)
    }

    # Subset by sequenceContext
    se <- .keepPositionsBySequenceContext(se, sequenceContext = sequenceContext)

    # Calculate statistics for each sample
    out <- SimpleList(lapply(
        structure(colnames(se), names = colnames(se)), function(nm) {
            sesub <- .filterPositionsByCoverage(
                se[, nm], assayName = assayName, minCov = minNobsPpos,
                minNbrSamples = NULL)

            mat <- assay(sesub, assayName)[[nm]]

            # Non-NA indices
            NNAind <- nnawhich(mat, arr.ind = TRUE)

            # Coverage per row (i.e per position)
            Nobs <- rep(0, nrow(mat))
            TBL <- table(NNAind[, 1])
            Nobs[as.numeric(names(TBL))] <- unclass(TBL)

            # Create list of non-NA row indices per column (i.e per read)
            NNAind_byCol <- split(NNAind[, 1], NNAind[, 2])
            names(NNAind_byCol) <- colnames(mat)[as.numeric(names(NNAind_byCol))]

            # List of non-zero observations by column (i.e by read):
            NNAvals <- nnavals(mat)
            NNAvals_byCol <- split(NNAvals, NNAind[, 2])
            names(NNAvals_byCol) <- colnames(mat)[as.numeric(names(NNAvals_byCol))]

            # Number of (valid) observations per read:
            NobsReads <- lengths(NNAind_byCol)

            # Collapsed mod probs per position:
            MeanModProb <- rowSums(mat) / Nobs

            # Include in calculations only reads with sufficient Number of observations:
            use.reads <- colnames(mat)[NobsReads >= minNobsPread]

            if (!is.null(stats)) {
                param_names <- stats
            } else {
                param_names <- defaultReadStats
            }

            # Iterate over param_names and add columns to stats_res
            do.call(cbind, bplapply(param_names, function(param, mycolnames = colnames(mat),
                                                          myuse.reads = use.reads,
                                                          mystatFunctions = statFunctions,
                                                          myNNAvals_byCol = NNAvals_byCol,
                                                          myLagRangeValues = LagRangeValues) {
                stats_res <- make_zero_col_DFrame(nrow = length(mycolnames))
                row.names(stats_res) <- mycolnames
                if (param %in% c("ACModProb", "PACModProb")) {
                    stats_res[[param]] <- lapply(
                        structure(mycolnames, names = mycolnames), function(r) {
                            if (r %in% myuse.reads) {
                                mystatFunctions[[param]](myNNAvals_byCol[[r]])
                            } else {
                                rep(NA, length(myLagRangeValues))
                            }
                        })
                } else {
                    stats_res[[param]] <- rep(NA, length(mycolnames))
                    stats_res[myuse.reads, param] <- vapply(myuse.reads, function(r) {
                        statFunctions[[param]](myNNAvals_byCol[[r]])
                    }, numeric(1))
                }
                stats_res
            }, BPPARAM = BPPARAM))
        })
    )

    # add filtering parameters to `out`
    metadata(out) <- list(regions = regions,
                          sequenceContext = sequenceContext,
                          minNobsPpos = minNobsPpos,
                          minNobsPread = minNobsPread,
                          Lags = LagRangeValues)

    return(out)
}

#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#'
#' @export
#' @rdname calcReadStats
addReadStats <- function(se, ..., name = "QC") {

    .assertScalar(x = name, type = "character")

    colData(se)[[name]] <- calcReadStats(se = se, ...)
    metadata(se)$readLevelData$colDataColumns <- c(
        metadata(se)$readLevelData$colDataColumns, name
    )
    return(se)
}

