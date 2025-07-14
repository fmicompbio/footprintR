# -- Global variables ----------------------------------------------------------
# global vector with default read stats (that will be calculated if
# stats = NULL in calcReadStats)
# exclude "SEntrModProb"
defaultReadStats <- c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "ACModProb", "PACModProb","SNR", "SignalVar", "NoiseVar")
# global vector with all available read stats functions
allReadStats <- c(defaultReadStats, "SEntrModProb")
# global vector with available signal-to-noise read stats functions
SNRstats <- c("SNR", "SignalVar", "NoiseVar")

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
    stats_res <- vapply(probList, mean, numeric(1), USE.NAMES = FALSE)
    stats_res[setdiff(seq_along(stats_res), useReads)] <- NA
    stats_res
}

#' @noRd
#' @keywords internal
FracMod <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        mean(probList[[r]] >= 0.5)
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
MeanConf <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        mean(pmax(probList[[r]], 1 - probList[[r]]))
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
MeanConfUnm <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        mean((1 - probList[[r]])[probList[[r]] < 0.5])
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
MeanConfMod <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        mean(probList[[r]][probList[[r]] >= 0.5])
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
FracLowConf <- function(probList, useReads, lowConf = 0.7, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        sum(abs(0.5 - probList[[r]]) < (lowConf - 0.5)) / length(probList[[r]])
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
#' @importFrom stats IQR
IQRModProb <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        IQR(probList[[r]])
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
#' @importFrom stats var
# remark: this is one of the few cases where working on the NAmatrix directly
#         would spead up things (about 2-fold), thanks to the SparseArray::colSds
#         but for consistency we keep the list-of-mod_prob version
sdModProb <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        sqrt(var(probList[[r]]))
    }, numeric(1))
    stats_res
}

#' @noRd
#' @keywords internal
SEntrModProb <- function(probList, useReads, ...) {
    stats_res <- rep(NA, length(probList))
    stats_res[useReads] <- vapply(useReads, function(r) {
        if (length(probList[[r]]) > 64) {
            sampleEntropy(probList[[r]], 2L, 0.2)
        } else {
            NA
        }
    }, numeric(1))
    stats_res
}



#' @noRd
#' @keywords internal
#' @importFrom stats acf na.pass
ACModProb <- function(probList, useReads, xrange = 12:64, ...) {
    lagMax <- max(xrange)
    stats_res <- lapply(lengths(probList), function(i) rep(NA, length(xrange)))
    stats_res[useReads] <- lapply(useReads, function(r) {
        if (length(probList[[r]]) > lagMax) {
            acf(probList[[r]], na.action = na.pass, lag.max = lagMax,
                plot = FALSE)$acf[xrange]
        } else {
            rep(0, length(xrange))
        }
    })
    stats_res
}

#' @noRd
#' @keywords internal
#' @importFrom stats pacf na.pass
PACModProb <- function(probList, useReads, xrange = 12:64, ...) {
    lagMax <- max(xrange)
    stats_res <- lapply(lengths(probList), function(i) rep(NA, length(xrange)))
    stats_res[useReads] <- lapply(useReads, function(r) {
        if (length(probList[[r]]) > lagMax) {
            pacf(probList[[r]], na.action = na.pass, lag.max = lagMax,
                 plot = FALSE)$acf[xrange]
        } else {
            rep(0, length(xrange))
        }
    })
    stats_res
}


#' Internal: estimate per-read SNR, signal and noise  (NA-gap aware)
#' 1. **Total variance**  = `var(x, na.rm = TRUE)`  
#' 2. **Noise variance**  ≈ `0.5 * Var(Δx)` where Δx are lag-1 differences that may skip
#'    up to *k* missing values. This follows from error propagation and the assumption
#'    of low varying x in nearby measurements: Var(Δx)≈2*Var(x)  
#' 3. A **noise floor** is imposed:  
#'    `noise = pmax(noise, b0 + b1 * mean(x))`,  
#'    where *b0*, *b1* are obtained from a robust linear fit
#'    (`quantile 0.1 – 0.9`) of noise ~ mean(x).  
#' 4. **Signal variance** = `pmax(total - noise, eps)` with a small floor *eps*.  
#' 5. **SNR**            =  `log2(signal / noise)`.
#'
#' @importFrom stats var quantile lm coef
#' @importFrom dplyr between
#' @importFrom cli cli_warn
#'
#' @noRd
#' @keywords internal
.estimate_snr_probList <- function(probList, idxList,
                                   k          = 2L,
                                   min_diffs  = NULL,
                                   floor_pars = NULL,
                                   eps        = 1e-3,
                                   ...)
{
    nReads <- length(probList)
    if (is.null(min_diffs))
        min_diffs <- max(16L, floor(0.05 * median(lengths(idxList))))
    
    ## means + total variances ------------------------------------------------
    m_means <- vapply(probList, mean, numeric(1), na.rm = TRUE)
    total_v <- vapply(probList, var , numeric(1), na.rm = TRUE)
    
    ## noise variance: lag-1 diffs that jump ≤ k NAs --------------------------
    noise_v <- vapply(seq_len(nReads), function(i) {
        x   <- probList[[i]]
        idx <- idxList [[i]]
        if (length(idx) < 2L) return(NA_real_)
        gaps  <- diff(idx) - 1L
        d     <- diff(x)[gaps <= k]
        if (length(d) >= min_diffs) var(d) / 2 else NA_real_
    }, numeric(1))
    
    ## robust noise floor -----------------------------------------------------
    if (is.null(floor_pars)) {
        keep <- dplyr::between(
            m_means,
            quantile(m_means, .10, na.rm = TRUE),
            quantile(m_means, .90, na.rm = TRUE)) &
            !is.na(noise_v)
        
        if (sum(keep) < 16) {
            ## Fallback: keep raw noise variances, no floor – but warn the user
            cli_warn("Too few points to estimate noise floor ({sum(keep)}); raw noise variances are used.")
            floor_pars   <- c(NA_real_, NA_real_)        # returned for bookkeeping
            fitted_floor <- rep(-Inf, length(noise_v))   # pmax() leaves noise_v untouched
        } else {
            floor_pars   <- coef(lm(noise_v[keep] ~ m_means[keep]))
            fitted_floor <- floor_pars[1] + floor_pars[2] * m_means
        }
    } else {
        fitted_floor <- floor_pars[1] + floor_pars[2] * m_means
    }
    
    noise_v <- pmax(noise_v, fitted_floor, na.rm = TRUE)
    
    ## signal + SNR -----------------------------------------------------------
    signal_v <- pmax(total_v - noise_v, eps)
    snr_v    <- log2(signal_v / noise_v)
    
    list(snr = snr_v,
         signal = signal_v,
         noise  = noise_v,
         floor_pars = setNames(floor_pars, c("intercept", "slope")))
}



#' @noRd
#' @keywords internal
SNR <- function(probList, idxList, useReads, ...) {
    snr_res  <- .estimate_snr_probList(probList = probList,
                                       idxList  = idxList, ...)
    out      <- rep(NA_real_, length(probList))
    out[useReads] <- snr_res$snr[useReads]
    out
}

#' @noRd
#' @keywords internal
SignalVar <- function(probList, idxList, useReads, ...) {
    snr_res  <- .estimate_snr_probList(probList = probList,
                                       idxList  = idxList, ...)
    out      <- rep(NA_real_, length(probList))
    out[useReads] <- snr_res$signal[useReads]
    out
}

#' @noRd
#' @keywords internal
NoiseVar <- function(probList, idxList, useReads, ...) {
    snr_res  <- .estimate_snr_probList(probList = probList,
                                       idxList  = idxList, ...)
    out      <- rep(NA_real_, length(probList))
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
#' @importFrom BiocGenerics colnames
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom BiocGenerics pos
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
            POS <- BiocGenerics::pos(rowRanges(sesub))
            
            # Non-NA indices
            NNAind <- nnawhich(mat, arr.ind = TRUE)
            
            # Positions of observed measurements
            idxPos_byCol <- split(POS[NNAind[,1]], NNAind[,2])             
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
            
            if (any(param_names %in% SNRstats)) {
                snr_res <- .estimate_snr_probList(
                    probList = NNAvals_byCol,
                    idxList  = idxPos_byCol
                )
            }
            
            # Iterate over param_names and add columns to stats_res
            do.call(
                cbind,
                bplapply(
                    param_names,
                    function(param) {
                        
                        stats_res <- make_zero_col_DFrame(nrow = length(colnames(mat)))
                        row.names(stats_res) <- colnames(mat)
                        vec <- rep(NA_real_, length(colnames(mat))); names(vec) <- colnames(mat)
                        
                        if (param %in% SNRstats) {
                            src <- switch(param,
                                          SNR       = snr_res$snr,
                                          SignalVar = snr_res$signal,
                                          NoiseVar  = snr_res$noise)
                            vec[useReads] <- src[useReads]
                        } else {
                            helper_args <- list(probList = NNAvals_byCol,
                                                idxList  = idxPos_byCol,
                                                useReads = useReads,
                                                lowConf  = LowConf,
                                                xrange   = LagRangeValues)
                            vec[names(NNAvals_byCol)] <- do.call(param, helper_args)
                        }
                        
                        stats_res[[param]] <- vec
                        stats_res
                    },
                    BPPARAM = BPPARAM))
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

