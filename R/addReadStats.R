#' Add a variety of read-level modified basecalling summary statistics to a SummarizedExperiment
#'
#' @description
#' This function calculates various per-read summary statistics on modification
#' probabilities or calls from a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with genomic positions in rows and reads in columns and adds them to the
#' \code{colData}.
#' See details for more information on the statistics that are calculated.
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
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with \code{colData} filtered for positions according to \code{regions},
#'     \code{sequence.context} and \code{min.Nobs.ppos} arguments and extended
#'     to include the read statistics in its row- and column-data.
#'
#' @author Panagiotis Papapasaikas
#'
#' @examples
#' library(SummarizedExperiment)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfile = modbamfile, regions = "chr1:6940000-6955000",
#'            modbase = "a", verbose = TRUE)
#' se_withReadStats <- addReadStats(se)
#' rowData(se_withReadStats)
#' colData(se_withReadStats)
#'
#' @importFrom BiocGenerics start
#' @import ggplot2
#' @importFrom tidyr gather
#' @importFrom S4Vectors metadata make_zero_col_DFrame
#' @importFrom SparseArray rowSums
#' @importFrom SummarizedExperiment assay rowData assayNames
#' @importFrom stats sd IQR acf pacf na.pass
#' @importFrom Biostrings vcountPattern
#' @importFrom rlang .data
#'
#' @export
addReadStats <- function(se,
                         regions = NULL,
                         sequence.context = NULL,
                         stats = NULL,
                         min.Nobs.ppos = NULL,
                         min.Nobs.pread = 0,
                         LowConf = 0.7,
                         LagRange = c(12, 64),
                         verbose = FALSE) {

    se$QC <- calcReadStats(
        se = se,
        regions = regions,
        sequence.context = sequence.context,
        stats = stats,
        min.Nobs.ppos = min.Nobs.ppos,
        min.Nobs.pread = min.Nobs.pread,
        LowConf = LowConf,
        LagRange = LagRange,
        verbose = verbose
    )[colnames(se)]

    return(se)
}

