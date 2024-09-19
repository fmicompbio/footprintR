#' Calculate a variety of read-level modified basecalling summary statistics
#'
#' @description
#' This function calculates various per-read summary statistics on modification
#' probabilities or calls from a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with genomic positions in rows and reads in columns. See details
#' for more information on the statistics that are calculated.
#'
#' @param se A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     object with assay \code{"mod_prob"} typically returned by
#'     \code{\link{readModkitExtract}} or \code{\link{readModBam}}.
#' @param stats Character vector specifying which statistics to calculate.
#'     When set to \code{NULL} all available statistics are calculated. See
#'     details for a complete list of available read statistics.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object limiting the
#'     positions included in the calculations to the ones overlapping the
#'     corresponding genomic regions. Alternatively, regions can be
#'     specified as a character vector (e.g. "chr1:1200-1300") that can be
#'     coerced into a \code{GRanges} object.
#' @param sequence.context A character vector with sequence context(s)
#'     to include in the calculations. Only positions that match one of the
#'     provided sequence contexts will be included. Sequence contexts can be
#'     provided using IUPAC redundancy codes. The sequence contexts of modified
#'     bases are obtained from \code{rowData(se)$sequence.context} and thus
#'     requires that \code{se} contains the appropriate information, for example
#'     by setting the \code{sequence.context} and \code{sequence.reference}
#'     arguments of \code{\link{readModkitExtract}} when it was generated,
#'     or by adding it using \code{\link{addSeqContext}}.
#' @param min.Nobs.ppos A numeric scalar value >=1 indicating the minimum
#'     coverage on individual positions for them to be included in the
#'     calculations. In high coverage data this is an effective filter for
#'     removing spurious modbases, typically the result of erroneous
#'     basecalling. The default \code{NULL} sets its value to Q3-0.5*IQR, where
#'     Q3 and IQR are the third quartile and interquartile range of the coverage
#'     distribution estimated from the data in \code{se}.
#' @param min.Nobs.pread A numeric scalar with the minimum number of observed
#'     modifiable bases per read for it to be included in the calculations.
#'     \code{NA} values are returned for the reads that do not pass this
#'     threshold.
#' @param LowConf A numeric scalar with the minimum call confidence below which
#'     calls are considered "low confidence".
#' @param LagRange A numeric vector of two values (minimum and maxium) defining
#'     the range of lags for the calculation of autocorrelation and partial
#'     autocorrelation (see details section).
#' @param verbose If \code{TRUE}, report on progress.
#'
#' @details
#' Calculates a collection of location/scatter statistics and information
#' theoretic/signal-processing metrics for the modification probability,
#' confidence or modification call value vectors across individual reads. When
#' \code{sequence.context}, \code{min.coverage} or \code{min.Nobs.ppos} filters
#' are enforced, only modifiable bases passing the filters are included in the
#' calculations. The following statistics are available:
#' \describe{
#'     \item{MeanModProb}{: Mean modification probability across the read.}
#'     \item{FracMod}{: Fraction of confidently modified bases, defined as the
#'         ratio of modifiable bases with modification probability >= 0.5 over
#'         all modifiable bases.}
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
#' @return A \code{SimpleList} object with QC statistics for the samples in
#' \code{se}.
#'
#' @author Panagiotis Papapasaikas, Charlotte Soneson
#'
#' @examples
#' library(SummarizedExperiment)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfile = modbamfile, regions = "chr1:6940000-6955000",
#'            modbase = "a", verbose = TRUE)
#' ReadStats <- calcReadStats(se)
#' ReadStats[["s1"]]
#'
#' @importFrom BiocGenerics start
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom S4Vectors metadata make_zero_col_DFrame SimpleList
#' @importFrom SparseArray rowSums nnawhich nnavals
#' @importFrom SummarizedExperiment assay rowData assayNames
#' @importFrom stats sd IQR acf pacf na.pass
#' @importFrom Biostrings vcountPattern
#' @importFrom IRanges subsetByOverlaps
#' @importFrom rlang .data
#'
#' @export
calcReadStats <- function(se,
                          regions = NULL,
                          sequence.context = NULL,
                          stats = NULL,
                          min.Nobs.ppos = NULL,
                          min.Nobs.pread = 0,
                          LowConf = 0.7,
                          LagRange = c(12, 64),
                          verbose = FALSE) {

    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (!all(c("mod_prob") %in% SummarizedExperiment::assayNames(se))) {
        stop("`se` needs to have a 'mod_prob' assay")
    }
    if (is.character(regions)) {
        regions <- as(regions, "GRanges")
    }
    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)
    .assertVector(x = sequence.context, type = "character", allowNULL = TRUE)
    .assertScalar(x = LowConf, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = LagRange, type = "vector", rngIncl = c(1, 256), len = 2)
    LagRangeValues <- seq(LagRange[1], LagRange[2])
    statFunctions <- list(
        MeanModProb = mean,
        FracMod = function(x, c = 0.5) {
            sum(x >= (0.5 + (c - 0.5))) / sum(abs(0.5 - x) > (c - 0.5))
        },
        MeanConf = function(x) {
            mean(pmax(x, 1 - x))
        },
        MeanConfUnm = function(x) {
            mean((1 - x)[x < 0.5])
        },
        MeanConfMod = function(x) {
            mean((x)[x >= 0.5])
        },
        FracLowConf = function(x, c = LowConf) {
            sum(abs(0.5 - x) < (c - 0.5)) / length(x)
        },
        IQRModProb = function(x) {
            stats::IQR(x)
        },
        sdModProb = function(x) {
            stats::sd(x)
        },
        SEntrModProb = function(x) {
             if (length(x) > 64) {
                 sampleEntropy(x, 2L, 0.2)
             } else {
                 NA
             }
         },
        Lag1DModProb = function(x) {
            xC <- 1 * (x > 0.5)
            mean(abs(diff(xC, lag = 1)))
        },
        ACModProb = function(x, lag.max = max(LagRange),
                             xrange = LagRangeValues) {
            if (length(x) > lag.max) {
                stats::acf(x, na.action = stats::na.pass, lag.max = lag.max,
                           plot = FALSE)$acf[xrange]
            } else {
                rep(0, length(xrange))
            }
        },
        PACModProb = function(x, lag.max = max(LagRange),
                              xrange = LagRangeValues) {
            if (length(x) > lag.max) {
                stats::pacf(x, na.action = stats::na.pass, lag.max = lag.max,
                            plot = FALSE)$acf[xrange]
            } else {
                rep(0, length(xrange))
            }
        }
    )
    .assertVector(x = stats, type = "character", allowNULL = TRUE,
                  validValues = names(statFunctions))
    .assertScalar(x = min.Nobs.ppos, type = "numeric", allowNULL = TRUE,
                  rngIncl = c(1, Inf))
    .assertScalar(x = min.Nobs.pread, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")

    # Subset se by region
    if (!is.null(regions)) {
        se <- subsetByOverlaps(x = se, ranges = regions)
        # % removed
    }

    # Subset by sequence.context
    if (!is.null(sequence.context)) {
        if (is.null(SummarizedExperiment::rowData(se)$sequence.context)) {
            stop("No sequence context found in `rowData(se)$sequence.context`")
        }
        nmatch <- Reduce("+", lapply(sequence.context, function(pat) {
            vcountPattern(pat,
                          SummarizedExperiment::rowData(se)$sequence.context,
                          fixed = FALSE)
        }))
        se <- se[nmatch > 0, ]
        # % removed
    }

    # Calculate statistics for each sample
    out <- SimpleList(lapply(colnames(se), function(nm) {
        mat <- SummarizedExperiment::assay(se, "mod_prob")[[nm]]

        # Non-NA indices:
        NNAind <- SparseArray::nnawhich(mat, arr.ind = TRUE)

        # Coverage per row (i.e per position)
        Nobs <- rep(0, nrow(mat))
        TBL <- table(NNAind[, 1])
        Nobs[as.numeric(names(TBL))] <- unclass(TBL)

        # Subset positions by coverage
        if (is.null(min.Nobs.ppos)) {
            min.cov <- stats::quantile(Nobs, 0.75) -
                0.5 * stats::IQR(Nobs)
        } else{
            min.cov <- min.Nobs.ppos
        }
        min.cov <- max(floor(min.cov), 1)
        idx <- which(Nobs >= min.cov)
        mat <- mat[idx, ]
        Nobs <- Nobs[idx]
        if (verbose) {
            message(
                "(", nm, ") Applied coverage filter\nPositions with coverage < ",
                min.cov, " removed.")
        }

        # Create list of non-zero row indices per column (i.e per read)
        NNAind <- SparseArray::nnawhich(mat, arr.ind = TRUE)
        NNAind_byCol <- split(NNAind[, 1], NNAind[, 2])
        names(NNAind_byCol) <- colnames(mat)[as.numeric(names(NNAind_byCol))]

        # List of non-zero observations by column (i.e by read):
        NNAvals <- SparseArray::nnavals(mat)
        NNAvals_byCol <- split(NNAvals, NNAind[, 2])
        names(NNAvals_byCol) <- colnames(mat)[as.numeric(names(NNAvals_byCol))]

        # Number of (valid) observations per read:
        NobsReads <- lengths(NNAind_byCol)

        ## TODO:
        # Add stats on removed positions / reads for each filter

        # Collapsed mod probs per position:
        MeanModProb <- SparseArray::rowSums(mat) / Nobs

        # "collapsed methylation" over the same positions as the read-level observations
        # READSTATS_6mA$meanMeth_CL <- sapply(1:ncol(Probs_6mA), function(x) {
        #     obs <- NNAindL[[x]]
        #     mean(collapsed_6mA_f[obs],na.rm=TRUE)
        # })

        # Include in calculations only reads with sufficient Number of observations:
        if (min.Nobs.pread > 0) {
            use.reads <- colnames(mat)[NobsReads > min.Nobs.pread]
        } else {
            use.reads <- colnames(mat)
        }

        if (!is.null(stats)) {
            param_names <- stats
        } else {
            param_names <- names(statFunctions)
        }

        # Iterate through the logicals and add columns to stats_res
        # if the parameter is TRUE
        stats_res <- S4Vectors::make_zero_col_DFrame(nrow = ncol(mat))
        row.names(stats_res) <- colnames(mat)
        metadata(stats_res) <- list(min.Nobs.ppos = min.cov,
                                    Lags = LagRangeValues,
                                    stats = param_names)
        for (param in param_names) {
            if (param %in% c("ACModProb", "PACModProb")) {
                stats_res[[param]] <- rep(list(rep(NA, length(LagRangeValues))),
                                          ncol(mat))
                stats_res[use.reads, param] <- I(lapply(use.reads, function(r) {
                    v <- NNAvals_byCol[[r]]
                    statFunctions[[param]](v)
                }))
            } else {
                stats_res[[param]] <- rep(NA, ncol(mat))
                stats_res[use.reads, param] <- vapply(use.reads, function(r) {
                    v <- NNAvals_byCol[[r]]
                    statFunctions[[param]](v)
                }, numeric(1))
            }
        }
        stats_res$sample <- nm
        stats_res
    }))
    names(out) <- colnames(se)

    return(out)
}

