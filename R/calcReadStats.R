# -- Global variables ----------------------------------------------------------
# global vector with default read stats (that will be calculated if
# stats = NULL in calcReadStats)
# exclude "SEntrModProb"
defaultReadStats <- c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "ACModProb", "PACModProb", "SNR", "SignalVar", "NoiseVar")
# global vector with all available read stats functions
allReadStats <- c(defaultReadStats, "SEntrModProb")
# global vector with available signal-to-noise read stats functions
snrStats <- c("SNR", "SignalVar", "NoiseVar")

# -- Helper functions to Individual read statistics ----------------------------
# the helper functions
# - are currently not exported
# - have the same name as the statistic in `defaultReadStats`
# - have at least two arguments:
#      probList: a list of numeric mod_prob vectors for each read
#      useReads: a vector defining the elements of probList for which to
#                calculate the statistic (results for the remaining elements
#                are NA)
#      ...     : optional arguments for specific statistics
# - return a vector or list with length(probList) elements
#' @noRd
#' @keywords internal
MeanModProb <- function(probList, useReads, ...) {
    statsRes <- vapply(probList, mean, numeric(1L), USE.NAMES = FALSE)
    statsRes[setdiff(seq_along(statsRes), useReads)] <- NA
    statsRes
}

#' @noRd
#' @keywords internal
FracMod <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        mean(probList[[r]] >= 0.5)
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
MeanConf <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        mean(pmax(probList[[r]], 1 - probList[[r]]))
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
MeanConfUnm <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        mean((1 - probList[[r]])[probList[[r]] < 0.5])
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
MeanConfMod <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        mean(probList[[r]][probList[[r]] >= 0.5])
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
FracLowConf <- function(probList, useReads, lowConf = 0.7, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        sum(abs(0.5 - probList[[r]]) < (lowConf - 0.5)) / length(probList[[r]])
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
#' @importFrom stats IQR
IQRModProb <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        IQR(probList[[r]])
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
#' @importFrom stats var
# remark: this is one of the few cases where working on the NAmatrix directly
#         would spead up things (about 2-fold), thanks to the SparseArray::colSds
#         but for consistency we keep the list-of-mod_prob version
sdModProb <- function(probList, useReads, ...) {
    statsRes <- rep(NA, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        sqrt(var(probList[[r]]))
    }, numeric(1))
    statsRes
}

#' @noRd
#' @keywords internal
SEntrModProb <- function(probList, useReads,
                         sampen_m=2, sampen_r=0.2, sampen_maxStarts=1000, sampen_nThreads=1, ...) {
    statsRes <- rep(NA_real_, length(probList))
    statsRes[useReads] <- vapply(useReads, function(r) {
        if (length(probList[[r]]) > 64) {
            sampleEntropy(probList[[r]], sampen_m, sampen_r, sampen_maxStarts, sampen_nThreads)
        } else {
            NA_real_
        }
    }, numeric(1))
    statsRes
}



#' @noRd
#' @keywords internal
#' @importFrom stats acf na.pass
ACModProb <- function(probList, useReads, xrange = 12:64, ...) {
    lagMax <- max(xrange)
    statsRes <- lapply(lengths(probList), function(i) rep(NA, length(xrange)))
    statsRes[useReads] <- lapply(useReads, function(r) {
        if (length(probList[[r]]) > lagMax) {
            acf(probList[[r]], na.action = na.pass, lag.max = lagMax,
                plot = FALSE)$acf[xrange]
        } else {
            rep(0, length(xrange))
        }
    })
    statsRes
}

#' @noRd
#' @keywords internal
#' @importFrom stats pacf na.pass
PACModProb <- function(probList, useReads, xrange = 12:64, ...) {
    lagMax <- max(xrange)
    statsRes <- lapply(lengths(probList), function(i) rep(NA, length(xrange)))
    statsRes[useReads] <- lapply(useReads, function(r) {
        if (length(probList[[r]]) > lagMax) {
            pacf(probList[[r]], na.action = na.pass, lag.max = lagMax,
                 plot = FALSE)$acf[xrange]
        } else {
            rep(0, length(xrange))
        }
    })
    statsRes
}


#' Internal: estimate per-read SNR, signal and noise  (NA-gap aware)
#' 1. **Total variance** = `var(x, na.rm = TRUE)`
#' 2. **Noise variance** ≈ `0.5 * Var(Δx)` where Δx are lag-1 differences that may skip
#'    up to *k* missing values. This follows from error propagation and the assumption
#'    of low varying x in adjacent measurements: Var(Δx)≈2*Var(x)
#' 3. A **noise floor** is imposed by fitting a background noise model:
#'    `noise = pmax(noise, b0 + b1 * mean(x))`,
#'    where *b0*, *b1* are obtained from a robust linear fit
#'    (`quantile 0.1 – 0.9`) of noise ~ mean(x).
#' 4. **Signal variance** = `pmax(total - noise, eps)` with a small floor *eps*.
#' 5. **SNR** = `log2(signal / noise)`.
#'
#' @importFrom stats var quantile lm coef
#' @importFrom dplyr between
#' @importFrom cli cli_warn
#' @importFrom stats setNames quantile lm coef
#'
#' @noRd
#' @keywords internal
.estimateSNRprobList <- function(probList, idxList, k = 2L, min_diffs = -1L,
                                 floor_pars = NULL, eps = 1e-3, ...) {
    nReads <- length(probList)

    # Noise estimation using cpp function:
    comps <- lapply(seq_len(nReads), function(i) {
        estimateNoise(
            probs = probList[[i]],
            read_pos = idxList[[i]],
            k = as.integer(k),
            min_diffs = as.integer(min_diffs)  # -1 triggers C++ auto per-read
        )
    })
    m_means <- vapply(comps, `[[`, numeric(1), "mean")
    noise_r <- vapply(comps, `[[`, numeric(1), "noise_raw")
    ndiffs <- vapply(comps, `[[`, numeric(1), "ndiffs")

    # Robust noise floor fit (if not provided)
    if (is.null(floor_pars)) {
        keep <- between(
            m_means,
            quantile(m_means, 0.10, na.rm = TRUE),
            quantile(m_means, 0.90, na.rm = TRUE)
        ) & !is.na(noise_r)

        if (sum(keep, na.rm = TRUE) < 16L) {
            cli_warn("Too few points to estimate noise floor ({sum(keep)}); raw noise variances are used.")
            floor_pars <- c(NA_real_, NA_real_)  # bookkeeping
            use_mode <- "raw"
        } else {
            # noise_raw ~ mean
            fit <- lm(noise_r[keep] ~ m_means[keep])
            floor_pars <- coef(fit) # c(b0, b1)
            use_mode <- "floor"
        }
    } else {
        if (all(is.finite(floor_pars))) {
            use_mode <- "floor"
        } else {
            use_mode <- "raw"
        }
    }

    # Per-read noise/signal/SNR in C++
    # Build betas/features for the chosen model (intercept + mean)
    if (all(!is.na(floor_pars))) {
        betas <- floor_pars
    } else {
        betas <- numeric(0L)  # raw mode
    }

    # run estimateSNR per read
    snr_sig_noise <- lapply(seq_len(nReads), function(i) {
        ci <- comps[[i]]
        if (length(ci) < 3L || !is.finite(ci[1]) || !is.finite(ci[2]) ||
            !is.finite(ci[3]) || ndiffs[i] < 2) {
            return(c(snr = NA_real_, signal = NA_real_, noise = NA_real_))
        }
        if (length(betas)) {
            feats <- c(1, m_means[i])
        } else {
            feats <- numeric(0L)
        }
        res <- estimateSNR(
            totalVar = ci[["total"]],
            noiseRaw = ci[["noise_raw"]],
            eps = eps,
            betas = betas,
            features = feats,
            noise_mode = use_mode # "raw" (fallback) or "floor"
        )
        c(snr = res[["snr"]], signal = res[["signal"]], noise = res[["noise"]])
    })

    snr_v <- vapply(snr_sig_noise, `[[`, numeric(1), "snr")
    signal_v <- vapply(snr_sig_noise, `[[`, numeric(1), "signal")
    noise_v <- vapply(snr_sig_noise, `[[`, numeric(1), "noise")

    list(
        snr = snr_v,
        signal = signal_v,
        noise = noise_v,
        floor_pars = setNames(floor_pars, c("intercept", "slope"))
    )
}

#' @noRd
#' @keywords internal
SNR <- function(probList, idxList, useReads, ...) {
    snr_res <- .estimateSNRprobList(probList = probList,
                                    idxList = idxList, ...)
    out <- rep(NA_real_, length(probList))
    out[useReads] <- snr_res$snr[useReads]
    out
}

#' @noRd
#' @keywords internal
SignalVar <- function(probList, idxList, useReads, ...) {
    snr_res <- .estimateSNRprobList(probList = probList,
                                    idxList  = idxList, ...)
    out <- rep(NA_real_, length(probList))
    out[useReads] <- snr_res$signal[useReads]
    out
}

#' @noRd
#' @keywords internal
NoiseVar <- function(probList, idxList, useReads, ...) {
    snr_res <- .estimateSNRprobList(probList = probList,
                                    idxList  = idxList, ...)
    out <- rep(NA_real_, length(probList))
    out[useReads] <- snr_res$noise[useReads]
    out
}

# -- Main function to calculate read statistics or add them to an SE -----------

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
#'     basecalling.
#' @param minNobsPread A numeric scalar with the minimum number of observed
#'     modifiable bases per read for it to be included in the calculations.
#'     \code{NA} values are returned for the reads that do not pass this
#'     threshold.
#' @param LowConf A numeric scalar with the minimum call confidence below which
#'     calls are considered "low confidence".
#' @param LagRange A numeric vector of two values (minimum and maximum) defining
#'     the range of lags for the calculation of autocorrelation and partial
#'     autocorrelation (see details section).
#' @param EntrControl Optional named list with elements
#'   \code{m}, \code{r}, \code{maxStarts}, \code{nThreads} to control
#'   Sample Entropy calculation. Missing elements use defaults: m=2L, r=0.2,
#'   maxStarts=1000, nThreads=1L. See also \code{\link{sampleEntropy}}
#' @param name For \code{addReadStats} only: A character scalar specifying the
#'     name to be used to store the result in the
#'     \code{\link[SummarizedExperiment]{colData}} of the output.
#' @param ... For \code{addReadStats} only: Additional arguments passed on to
#'     \code{calcReadStats}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for calculating
#'     read statistics. The default value is
#'     (\code{\link[BiocParallel]{MulticoreParam}(4L, RNGseed = 42L)}).
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
#'     \item{ACModProb}{: Autocorrelation of the modification probability values
#'         for lags in the range \code{LagRange}. This range typically covers
#'         the signal of nucleosome periodicity.}
#'     \item{PACModProb}{: Partial autocorrelation of the modification
#'         probability values for lags in the range \code{LagRange}. This range
#'         typically covers the signal of nucleosome periodicity.}
#'     \item{NoiseVar}{: raw **Read Noise variance** estimated as
#'         \eqn{0.5\,\mathrm{Var}(\Delta x)}, where \eqn{\Delta x} are
#'         lag-1 differences that may skip up to \eqn{k} consecutive NAs.
#'         A floor is applied: \eqn{\mathrm{noise} =
#'         \max(\mathrm{rawNoise},\, b_0 + b_1 \bar{x})}.}
#'     \item{SignalVar}{: **Read Signal variance**
#'         \eqn{\max(\mathrm{totalVar}-\mathrm{NoiseVar},\,\varepsilon)}.}
#'     \item{SNR}{: **Read Signal-to-Noise Ratio**
#'         \eqn{\log_2(\mathrm{SignalVar}/\mathrm{NoiseVar})}.}
#'
#'  }
#' When SNR-related statistics are requested,
#' a robust noise floor `b0 + b1 * mean(x)` is fitted per sample (see *NoiseVar*).
#' The fitted coefficients are stored in the result metadata (see below). If
#' SNR-related statistics are not computed, `NA`/`NA` placeholders are stored
#' for uniformity.
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
#'            modbase = "a", verbose = TRUE,
#'            BPPARAM = BiocParallel::SerialParam())
#'
#' readStats <- calcReadStats(se, BPPARAM = BiocParallel::SerialParam())
#' readStats$s1
#'
#' se_withReadStats <- addReadStats(se, name = "QC",
#'                                  BPPARAM = BiocParallel::SerialParam())
#' se_withReadStats$QC$s1
#' metadata(se_withReadStats$QC$s1)
#'
#' @importFrom S4Vectors metadata make_zero_col_DFrame SimpleList
#' @importFrom SummarizedExperiment assay
#' @importFrom SparseArray rowSums nnawhich nnavals
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocGenerics colnames pos
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom utils modifyList
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
                          EntrControl = NULL,
                          BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                          verbose = FALSE) {
    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertVector(x = stats, type = "character", allowNULL = TRUE,
                  validValues = allReadStats)
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

    if (is.null(EntrControl)) {
        EntrControl <- list()
    }

    SEctrl <- modifyList(
        list(m = 2L, r = 0.2, maxStarts = 1000, nThreads = 1L),
        EntrControl
    )


    # Assert SE ctrl arguments
    .assertScalar(SEctrl$m, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(SEctrl$r, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(SEctrl$maxStarts, type = "numeric",  rngIncl = c(-1, Inf))
    .assertScalar(SEctrl$nThreads, type = "numeric", rngIncl = c(1, Inf))


    # Subset se by region
    if (!is.null(regions)) {
        se <- subsetByOverlaps(x = se, ranges = regions)
    }

    # Subset by sequenceContext
    se <- .keepPositionsBySequenceContext(se, sequenceContext = sequenceContext)

    # Calculate statistics for each sample
    sample_out <- lapply(
        structure(colnames(se), names = colnames(se)), function(nm) {
            sesub <- .filterPositionsByCoverage(
                se[, nm], assayName = assayName, minCov = minNobsPpos,
                minNbrSamples = NULL)

            mat <- assay(sesub, assayName)[[nm]]
            POS <- pos(rowRanges(sesub))

            # Non-NA indices
            NNAind <- nnawhich(mat, arr.ind = TRUE)

            # Positions of observed measurements
            idxPos_byCol <- split(POS[NNAind[, 1]], NNAind[, 2])
            names(idxPos_byCol) <- colnames(mat)[as.numeric(names(idxPos_byCol))]

            # Create list of non-NA row indices per column (i.e per read)
            NNAind_byCol <- split(NNAind[, 1], NNAind[, 2])
            names(NNAind_byCol) <- colnames(mat)[as.numeric(names(NNAind_byCol))]

            # List of non-zero observations by column (i.e by read):
            NNAvals <- nnavals(mat)
            NNAvals_byCol <- split(NNAvals, NNAind[, 2])
            names(NNAvals_byCol) <- colnames(mat)[as.numeric(names(NNAvals_byCol))]

            # Number of (valid) observations per read:
            NobsReads <- lengths(NNAind_byCol)

            # Include in calculations only reads with sufficient Number of observations
            useReads <- which(NobsReads >= minNobsPread)

            if (!is.null(stats)) {
                param_names <- stats
            } else {
                param_names <- defaultReadStats
            }

            # get snr_res / floor_pars
            if (any(param_names %in% snrStats)) {
                snr_res <- .estimateSNRprobList(
                    probList = NNAvals_byCol,
                    idxList = idxPos_byCol
                )
                # always store a 2-length named vector
                floor_pars <- snr_res$floor_pars  # c(intercept=..., slope=...)
            } else {
                snr_res <- NULL
                # uniform shape even when SNR stats were not computed
                floor_pars <- c(intercept = NA_real_, slope = NA_real_)
            }

            # Iterate over param_names and add columns to statsRes
            stats_df <- do.call(
                cbind,
                bplapply(
                    param_names,
                    function(param) {

                        statsRes <- make_zero_col_DFrame(nrow = length(colnames(mat)))
                        row.names(statsRes) <- colnames(mat)
                        vec <- rep(NA_real_, length(colnames(mat)))
                        names(vec) <- colnames(mat)

                        if (param %in% snrStats) {
                            src <- switch(param,
                                          SNR = snr_res$snr,
                                          SignalVar = snr_res$signal,
                                          NoiseVar = snr_res$noise)
                            vec[useReads] <- src[useReads]
                        } else {
                            helper_args <- list(probList = NNAvals_byCol,
                                                idxList = idxPos_byCol,
                                                useReads = useReads,
                                                lowConf = LowConf,
                                                xrange = LagRangeValues,
                                                # SampEn controls
                                                sampen_m = SEctrl$m,
                                                sampen_r = SEctrl$r,
                                                sampen_maxStarts = SEctrl$maxStarts,
                                                sampen_nThreads = SEctrl$nThreads)
                            vec[names(NNAvals_byCol)] <- do.call(param, helper_args)
                        }

                        statsRes[[param]] <- vec
                        statsRes
                    },
                    BPPARAM = BPPARAM))

            list(stats = stats_df, floor_pars = floor_pars)
        })

    # Assemble stats output from the sample_out lapply returns:
    out <- SimpleList(lapply(sample_out, `[[`, "stats"))

    # Assemble floor_pars output from the sample_out lapply returns:
    noiseCoefs_by_sample <- lapply(sample_out, `[[`, "floor_pars")

    # add filtering parameters to `out`
    metadata(out) <- list(
        regions = regions,
        sequenceContext = sequenceContext,
        minNobsPpos = minNobsPpos,
        minNobsPread = minNobsPread,
        Lags = LagRangeValues,
        snr_noise_coef = noiseCoefs_by_sample, # named list: one c(intercept, slope) per sample
        snr_config = list(k = 2L, min_diffs = NULL, eps = 1e-3)
    )

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

