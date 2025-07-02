#' Count pairs of modified bases by distance and modification state
#'
#' For all pairs of bases with modification calls in a read, tabulate the
#' number of pairs at a given distance with a given modification state.
#'
#' @param bamfile Character scalar with the path to a \code{modBAM}
#'     file, containing information about base modifications in \code{MM} and
#'     \code{ML} tags. The \code{bamfile} must have an index.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic regions to extract the reads from. Alternatively, regions can be
#'     specified as a character vector (e.g. "chr1:1200-1300", "chr2:-6000",
#'     "chr1:10-", "chrM" or ".").
#' @param modbase Character scalar defining the modified base.
#' @param threshUnmod,threshMod Numeric scalars defining how to convert
#'     modification probabilities \code{p} to modification states. Bases with
#'     \code{p < threshUnmod} will be considered unmodified, and bases with
#'     \code{p >= threshMod} modified. All remaining bases will be ignored.
#' @param windowSize Numeric scalar giving the maximum window size
#'     covering pairs of modified bases to consider in pair-counting mode.
#'     A window size of 1 corresponds to a single base, a size of 2 to
#'     directly adjacent bases, etc.
#' @param minMapQ Numeric scalar giving the minimal mapping quality to include
#'     alignments in pair-counting mode.
#' @param minAlignedLength Numeric scalar giving the minimal alignment length
#'     to include alignments in pair-counting mode.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for decompressing
#'     bam records.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Charlotte Soneson, Michael Stadler
#'
#' @return A \code{DataFrame} with \code{windowSize} rows and five columns.
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' res <- countStatePairs(bamfile = modbamfile,
#'                        regions = "chr1",
#'                        modbase = "a", windowSize = 300,
#'                        BPPARAM = BiocParallel::SerialParam())
#' res
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bpnworkers
#' @importFrom cli cli_abort
#'
#' @export
countStatePairs <- function(bamfile,
                            regions = ".",
                            modbase,
                            threshUnmod = 0.5,
                            threshMod = 0.5,
                            windowSize = 200,
                            minMapQ = 0,
                            minAlignedLength = 0,
                            BPPARAM = BiocParallel::MulticoreParam(4L),
                            verbose = FALSE) {
    .assertScalar(x = bamfile, type = "character")
    if (!file.exists(bamfile)) {
        cli_abort("{.arg bamfile} does not exist: {.file {bamfile}}")
    }
    if (is(regions, "GRanges")) {
        regions <- as.character(regions, ignore.strand = TRUE)
    }
    .assertVector(x = regions, type = "character")
    .assertScalar(x = modbase, type = "character")
    .assertValidModbase(modbase)
    .assertScalar(x = threshUnmod, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = threshMod, type = "numeric", rngIncl = c(0, 1))
    if (threshUnmod > threshMod) {
        cli_abort(paste0(
            "{.arg threshUnmod} ({threshUnmod}) must be less than or equal ",
            "to {.arg threshMod} ({threshMod})"))
    }
    .assertScalar(x = windowSize, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = minMapQ, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedLength, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")

    ncpuDecompression <- bpnworkers(BPPARAM)

    resL <- read_modbam_cpp(inname_str = bamfile,
                            regions = regions,
                            modbase = modbase,
                            n_alns_to_sample = 0L,
                            tnames_for_sampling = character(0),
                            variantRefNames = character(0),
                            variantRefPositions = integer(0),
                            threshUnmod = threshUnmod,
                            threshMod = threshMod,
                            windowSize = as.integer(windowSize),
                            minMapQ = as.integer(minMapQ),
                            minAlignedLength = as.integer(minAlignedLength),
                            n_threads = ncpuDecompression,
                            verbose = verbose)

    res <- DataFrame(S = seq.int(windowSize),
                     unmod_unmod = resL$pair_counts[, 1],
                     unmod_mod = resL$pair_counts[, 2],
                     mod_unmod = resL$pair_counts[, 3],
                     mod_mod = resL$pair_counts[, 4])
    return(res)
}
