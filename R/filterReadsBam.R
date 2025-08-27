#' Write bam records from \code{infile} to \code{outfile} if they pass filter
#' criteria.
#'
#' For each bam file in \code{infiles}, parse alignments, calculate read
#' statistics and write the alignment to the corresponding output file from
#' \code{outfiles} if the read passes all criteria defined by the filtering
#' arguments. Filters are processed hierarchically: If a read does not pass a
#' given filter, the remaining filters will not be examined and the processing
#' continues with the next read.
#' The filters are examined in this order: \code{keepUnmapped},
#' \code{keepSecondary}, \code{keepSupplementary}, \code{minReadLength},
#' \code{minAlignedLength}, \code{minAlignedFraction}, \code{minQscore},
#' \code{minSNR}, \code{maxFracLowConf}, \code{maxEntropy}.
#'
#' @inheritParams filterReads
#'
#' @param infiles Character vector with name(s) of the input bam file(s).
#' @param outfiles Character vector with name(s) of the output bam file(s).
#'     Needs to have the same length as \code{infiles}.
#' @param modbase Character scalar defining the modified base to analyze
#'     (used by \code{maxEntropy} and \code{maxFracLowConf}).
#' @param indexOutfiles Logical scalar. If \code{TRUE} (the default) create
#'     a bam index file (\code{.bai} file) for each of the generated
#'     \code{outfiles}.
#' @param overwriteOutfiles Logical scalar. If \code{FALSE} (the default),
#'     existing \code{outfiles} will not be overwritten and the function will
#'     abort with an error message.
#' @param keepUnmapped,keepSecondary,keepSupplementary Logical scalars
#'     indicating whether to keep unmapped, secondary or supplementary
#'     alignments.
#' @param maxEntropy A numeric scalar representing the largest acceptable
#'     read-level entropy. Reads without modified-base calls or with entropy
#'     above this value will be filtered out. A value of \code{Inf} deactivates
#'     the entropy filter.
#' @param minSNR Numeric scalar. Minimum acceptable read Signal-to-Noise Ratio (SNR).
#'   Reads with SNR below this value are filtered out. Set to \code{-Inf}
#'   to disable SNR filtering. The SNR is computed per read as
#'   \eqn{\log_2(\mathrm{SignalVar}/\mathrm{NoiseVar})} where
#'   \eqn{\mathrm{NoiseVar} \approx 0.5\,\mathrm{Var}(\Delta x)} using
#'   adjacent methylation differences that can jump up to \eqn{k=2} gaps, and
#'   \eqn{\mathrm{SignalVar} = \max(\mathrm{Var}(x) - \mathrm{NoiseVar}, \varepsilon)} with \eqn{\varepsilon=10^{-3}}.
#' @param noiseCoef Numeric vector of length 2 giving the background noise model
#'   coefficients \code{(b0, b1)}. If provided, these are used to impose a floor
#'   on the estimated per-read noise variance (see also \code{\link{calcReadStats}}).
#' @param LowConf A numeric scalar with the minimum call confidence below which
#'     calls are considered "low confidence".
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for some of the steps
#'     in \code{filterReadsBam()}. The default value is
#'     (\code{\link[BiocParallel]{MulticoreParam}(4L, RNGseed = 42L)}).
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return \code{filterReadsBam} is called for its side effect of generating
#'     new bam files containing the subset of bam records from input bam files
#'     that pass all filtering criteria. In addition, it returns a
#'     \code{data.frame} with one row per \code{infiles} giving the numbers of
#'     bam records that were read in \code{total}, that were \code{retained} in
#'     the \code{outfiles} and that were filtered-out by reason of exclusion.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))
#' res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
#'                       modbase = "a", indexOutfiles = FALSE, minReadLength = 6746,
#'                       minAlignedLength = 6896, minAlignedFraction = 0.56, minSNR=-0.768,
#'                       minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
#'                       BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
#' res
#' unlink(filtbamfiles)
#'
#' @author Michael Stadler
#'
#' @importFrom BiocParallel MulticoreParam bplapply bpnworkers bpworkers<-
#' @importFrom cli cli_abort cli_alert_info
#'
#' @export
filterReadsBam <- function(infiles,
                           outfiles,
                           modbase,
                           indexOutfiles = TRUE,
                           overwriteOutfiles = FALSE,
                           keepUnmapped = TRUE,
                           keepSecondary = TRUE,
                           keepSupplementary = TRUE,
                           minReadLength = 0,
                           minAlignedLength = 0,
                           minAlignedFraction = 0,
                           minQscore = 0.0,
                           minSNR = -Inf,
                           maxFracLowConf = 1.0,
                           maxEntropy = Inf,
                           noiseCoef = c(NA_real_, NA_real_),
                           LowConf = 0.7,
                           BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                           verbose = FALSE) {
    # validate arguments
    .assertVector(x = infiles, type = "character")
    if (any(i <- !file.exists(infiles))) {
        cli_abort("not all {.arg infiles} exist: {infiles[i]}")
    }
    .assertVector(x = outfiles, type = "character", len = length(infiles))
    if (!overwriteOutfiles && any(i <- file.exists(outfiles))) {
        cli_abort("existing {.arg outfiles} would be overwritten: {outfiles[i]}")
    }
    .assertScalar(x = keepUnmapped, type = "logical")
    .assertScalar(x = keepSecondary, type = "logical")
    .assertScalar(x = keepSupplementary, type = "logical")
    .assertScalar(x = modbase, type = "character")
    .assertScalar(x = indexOutfiles, type = "logical")
    .assertScalar(x = minReadLength, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedLength, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedFraction, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minQscore, type = "numeric")
    .assertScalar(x = minSNR, type = "numeric")
    .assertScalar(x = maxEntropy, type = "numeric")
    .assertScalar(x = LowConf, type = "numeric", rngIncl = c(0.5, 1))
    .assertScalar(x = maxFracLowConf, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = noiseCoef, len=2, type = "numeric")
    

    # determine the number of parallel threads to be used for
    # bam files (chromosomes) and decompression of bam records
    ncpuTotal <- bpnworkers(BPPARAM)
    ncpuDecompression <- 2L # accept some level of over-subscription

    # iterate over bam files
    res <- data.frame(
        sample = if (is.null(names(infiles))) paste0("s", seq_along(infiles)) else names(infiles),
        infile = infiles,
        outfile = outfiles,
        do.call(
            rbind,
            lapply(seq_along(infiles), function(i) {
                # get chromosome names
                if (ncpuTotal > 1) {
                    chrs <- paste0(getChromosomeNamesFromBam(infiles[i]), ":")
                    if (keepUnmapped) {
                        chrs <- c(chrs, "*")
                    }
                } else {
                    chrs <- "."
                }

                # create temporary output file names
                tmpsamfiles <- tempfile(pattern = sprintf("file%05d_", seq_along(chrs)),
                                        fileext = ".sam")

                # filter in parallel
                if (verbose) {
                    cli_alert_info(
                        paste0("start filtering of {.file {infiles[i]}} ",
                               "using {ncpuTotal} thread{?s}"))
                }
                res1PerChr <- bplapply(
                    seq_along(chrs)[order(rep_len(seq.int(ncpuTotal), length(chrs)))],
                    function(j,
                             myinfile = infiles[i],
                             myoutfile = tmpsamfiles[j],
                             mymodbase = modbase,
                             myregion = chrs[j],
                             myKeepUnmapped = keepUnmapped,
                             myKeepSecondary = keepSecondary,
                             myKeepSupplementary = keepSupplementary,
                             myMinReadLength = as.integer(minReadLength),
                             myMinAlignedLength = as.integer(minAlignedLength),
                             myMinAlignedFraction = minAlignedFraction,
                             myMinQscore = minQscore,
                             myMinSNR = minSNR,
                             myMaxFracLowConf = maxFracLowConf,
                             myMaxEntropy = ifelse(is.finite(maxEntropy), maxEntropy, -1.0),
                             myLowConf = LowConf,
                             myNThreads = ncpuDecompression,
                             myverbose = FALSE) {
                        filter_modbam_cpp(infile = myinfile,
                                          outfile = myoutfile,
                                          modbase = mymodbase,
                                          region = myregion,
                                          includeHeader = TRUE,
                                          keepUnmapped = myKeepUnmapped,
                                          keepSecondary = myKeepSecondary,
                                          keepSupplementary = myKeepSupplementary,
                                          minReadLength = myMinReadLength,
                                          minAlignedLength = myMinAlignedLength,
                                          minAlignedFraction = myMinAlignedFraction,
                                          minQscore = myMinQscore,
                                          minSNR = myMinSNR,
                                          maxFracLowConf = myMaxFracLowConf,
                                          maxEntropy = myMaxEntropy,
                                          LowConf = myLowConf,
                                          noiseCoefB0 = noiseCoef[1],
                                          noiseCoefB1 = noiseCoef[2],
                                          nThreads = myNThreads,
                                          verbose = myverbose)
                    }, BPPARAM = BPPARAM)

                # merge partial outputs
                if (verbose) {
                    cli_alert_info("merging {length(tmpsamfiles)} filtered chunks")
                }
                concatenate_hts_files(input_files = tmpsamfiles,
                                      output_file = outfiles[i],
                                      ncpu = ncpuTotal)
                unlink(tmpsamfiles)

                # sum and return filter statistics
                res1 <- Reduce(f = "+", x = res1PerChr)
                if (verbose) {
                    cli_alert_info(
                        paste0("done filtering: retained {res1['retained']} ",
                               "of {res1['total']} records ({round(res1['retained']/res1['total']*100, 1)}%)"))
                }
                return(res1)
            })))

    if (indexOutfiles) {
        .message("indexing {length(outfiles)} output file{?s}")
        idxfiles <- bplapply(outfiles, function(fn) {
            index_bam_cpp(infile = fn) # nocov
        }, BPPARAM = BPPARAM)
    }

    return(res)
}

