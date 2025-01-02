#' Write bam records from \code{infile} to \code{outfile} if they pass filter
#' criteria.
#'
#' For each bam file in \code{infiles}, parse alignments, calculate read
#' statistics and write the alignment to the corresponding output file from
#' \code{outfiles} if the read passes all criteria defined by the filtering
#' arguments. Filters are processed hierarchically: If a read does not pass a
#' given filter, the remaining filters will not be examined and the processing
#' continues with the next read.
#' The filters are examined in this order: \code{minReadLength},
#' \code{minAlignedLength}, \code{minAlignedFraction}, \code{minQscore},
#' \code{maxFracLowConf}, \code{maxEntropy}.
#'
#' @param infiles Character vector with name(s) of the input bam file(s).
#' @param outfiles Character vector with name(s) of the output bam file(s).
#'     Needs to have the same length as \code{infiles}.
#' @param modbase Character scalar defining the modified base to analyze
#'     (used by \code{maxEntropy} and \code{maxFracLowConf}).
#' @param indexOutfiles Logical scalar. If \code{TRUE} (the default) create
#'     a bam index file (\code{.bai} file) for each of the generated
#'     \code{outfiles}.
#' @param minReadLength A numeric scalar representing the smallest acceptable
#'     read length. Reads that are shorter than this value will be filtered
#'     out.
#' @param minAlignedLength A numeric scalar representing the smallest acceptable
#'     aligned length. Reads with aligned length shorter than this value will
#'     be filtered out.
#' @param minAlignedFraction A numeric scalar representing the smallest
#'     acceptable aligned fraction of a read. Reads where the aligned fraction
#'     is smaller than this value will be filtered out.
#' @param minQscore A numeric scalar representing the smallest acceptable
#'     read-level Qscore. Reads with Qscore below this value will be filtered
#'     out.
#' @param maxFracLowConf A numeric scalar representing the maximally acceptable
#'     fraction of low-confidence modified base calls in a read. Reads with
#'     a fraction of low confidence calls greater than this value will be
#'     filtered out.
#' @param maxEntropy A numeric scalar representing the largest acceptable
#'     read-level entropy. Reads with entropy above this value will be filtered
#'     out. A value of \code{Inf} deactivates the entropy filter.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for some of the steps
#'     in \code{filterReadsBam()}. The default value is
#'     (\code{\link[BiocParallel]{MulticoreParam}(4L, RNGseed = 42L)}).
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A \code{data.frame} with one row per \code{infile} giving the numbers
#'     of filtered out records per reason for exclusion.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))
#' res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
#'                       modbase = "a", indexOutfiles = FALSE, minReadLength = 6746,
#'                       minAlignedLength = 6896, minAlignedFraction = 0.56,
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
                           minReadLength = 0,
                           minAlignedLength = 0,
                           minAlignedFraction = 0,
                           minQscore = 0.0,
                           maxFracLowConf = 1.0,
                           maxEntropy = Inf,
                           BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                           verbose = FALSE) {
    # validate arguments
    .assertVector(x = infiles, type = "character")
    if (any(i <- !file.exists(infiles))) {
        cli_abort(paste0("not all `infiles` exist: ",
                         paste(infiles[i], collapse = ", ")))
    }
    .assertVector(x = outfiles, type = "character", len = length(infiles))
    if (any(i <- file.exists(outfiles))) {
        cli_abort(paste0("existing `outfiles` would be overwritten: ",
                         paste(outfiles[i], collapse = ", ")))
    }
    .assertScalar(x = modbase, type = "character")
    .assertScalar(x = indexOutfiles, type = "logical")
    .assertScalar(x = minReadLength, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedLength, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedFraction, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minQscore, type = "numeric")
    .assertScalar(x = maxEntropy, type = "numeric")
    .assertScalar(x = maxFracLowConf, type = "numeric", rngIncl = c(0, 1))

    # determine the number of parallel threads to be used for
    # bam files (preferred) and decompression of bam records (if available)
    # (accept some level of over-subscription)
    ncpuTotal <- bpnworkers(BPPARAM)
    if (is(BPPARAM, "MulticoreParam") || is(BPPARAM, "SnowParam")) {
        ncpuFiles <- min(ncpuTotal, length(infiles))
        oversubscriptionRate <- 2.0
        ncpuDecompression <- min(8L, max(1L, as.integer(
            floor(oversubscriptionRate * ncpuTotal / ncpuFiles))))
        bpworkers(BPPARAM) <- ncpuFiles
        on.exit(bpworkers(BPPARAM) <- ncpuTotal)
    } else {
        ncpuDecompression <- 1L
    }

    # iterate over bam files
    res <- data.frame(
        sample = if (is.null(names(infiles))) paste0("s", seq_along(infiles)) else names(infiles),
        infile = infiles,
        outfile = outfiles,
        do.call(
            rbind,
            bplapply(seq_along(infiles),
                     function(i,
                              myinfile = infiles[i],
                              myoutfile = outfiles[i],
                              mymodbase = modbase,
                              myMinReadLength = as.integer(minReadLength),
                              myMinAlignedLength = as.integer(minAlignedLength),
                              myMinAlignedFraction = minAlignedFraction,
                              myMinQscore = minQscore,
                              myMaxEntropy = ifelse(is.finite(maxEntropy), maxEntropy, -1.0),
                              myMaxFracLowConf = maxFracLowConf,
                              myNThreads = ncpuDecompression,
                              myverbose = verbose) {
                         if (myverbose) {
                             cli_alert_info(
                                 paste0("opening input file {.file {myinfile}} ",
                                        "using {myNThreads} thread{?s}"))
                         }
                         res1 <- filter_modbam_cpp(infile = myinfile,
                                                   outfile = myoutfile,
                                                   modbase = mymodbase,
                                                   minReadLength = myMinReadLength,
                                                   minAlignedLength = myMinAlignedLength,
                                                   minAlignedFraction = myMinAlignedFraction,
                                                   minQscore = myMinQscore,
                                                   maxEntropy = myMaxEntropy,
                                                   maxFracLowConf = myMaxFracLowConf,
                                                   nThreads = myNThreads,
                                                   verbose = myverbose)
                         if (myverbose) {
                             cli_alert_info(
                                 paste0("done filtering: retained {res1['retained']} ",
                                        "of {res1['total']} records ({round(res1['retained']/res1['total']*100, 1)}%)"))
                         }
                         return(res1)
                     }, BPPARAM = BPPARAM)))

    if (indexOutfiles) {
        .message("indexing {length(outfiles)} output file{?s}")
        idxfiles <- bplapply(outfiles, function(fn) {
            index_bam_cpp(infile = fn) # nocov
        }, BPPARAM = BPPARAM)
    }

    return(res)
}
