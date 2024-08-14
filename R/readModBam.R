#' Read base modifications from a bam file.
#'
#' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
#' section 1.7) and return a list of information on modified bases.
#'
#' @param bamfile Character scalar specifying the path to a \code{modBAM}
#'     file. The BAM file must have an index.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic regions to extract the reads from. Alternatively, regions can be
#'     specified as a character scalar (e.g. "chr1:1200-1300") that can be
#'     coerced into a \code{GRanges} object. Note that the reads are not
#'     trimmed to the boundaries of the specified ranges. As a result, returned
#'     positions will typically extend out of the specified regions.
#' @param modbase Character scalar defining the modified base to extract.
#' @param seqinfo \code{NULL} or a \code{\link[GenomeInfoDb]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Useful to set the sorting order of sequence names.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with genomic positions in rows and reads in columns.
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' readModBam(bamfile = modbamfile, regions = "chr1:6940000-6955000",
#'            modbase = "a", verbose = TRUE)
#'
#' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
#'     SAM ML and MM tags for base modifications.
#'
#' @author Michael Stadler
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#' @importFrom SparseArray SparseArray
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges GPos sort match
#' @importFrom S4Vectors DataFrame
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics pos strand
#'
#' @export
readModBam <- function(bamfile,
                       regions,
                       modbase,
                       seqinfo = NULL,
                       verbose = FALSE) {
    # digest arguments
    .assertScalar(x = bamfile, type = "character")
    if (!file.exists(bamfile)) {
        stop("bamfile (", bamfile, ") does not exist.")
    }
    if (is.character(regions)) {
        regions <- as(regions, "GRanges")
    }
    .assertScalar(x = regions, type = "GRanges")
    # for valid values of `modbase`, see https://samtools.github.io/hts-specs/SAMtags.pdf (section 1.7)
    .assertScalar(x = modbase, type = "character",
                  validValues = c("m","h","f","c","C","g","e","b","T",
                                  "U","a","A","o","G","n","N"))
    if (!is.null(seqinfo)) {
        if (!is(seqinfo, "Seqinfo") &&
            (!is.numeric(seqinfo) || is.null(names(seqinfo)))) {
            stop("`seqinfo` must be `NULL`, a `Seqinfo` object or a named",
                 " numeric vector with genomic sequence lengths.")
        }
    }
    .assertScalar(x = verbose, type = "logical")

    # extract modification probabilities from `bamfile`
    resL <- read_modbam(inname_str = bamfile,
                        regions = as.character(regions, ignore.strand = TRUE),
                        modbase = modbase,
                        verbose = verbose)

    # remove unaligned and convert coordinates 0-based to 1-based
    keep <- resL$ref_position != -1
    if (verbose) {
        message("filtering out ", sum(!keep), " of ", length(keep),
                " unaligned modification events (e.g. soft-masked)")
    }
    resL <- lapply(resL, "[", keep)
    resL$ref_position <- resL$ref_position + 1L

    # convert inferred `call_prob` to our minimal value as in readModkitExtract()
    # (inferred means that the modification was omitted from the BAM file, e.g.
    #  dorado omits base modification probabilities less than 0.05, and
    #  read_modbam returns a call_probability of -1 for these)
    resL$call_prob[resL$call_prob == -1] <- 0.02

    # convert to SummarizedExperiment
    gpos_all <- GenomicRanges::GPos(seqnames = resL$ref_name,
                                    pos = resL$ref_position,
                                    strand = resL$strand,
                                    seqinfo = seqinfo)
    gpos <- GenomicRanges::sort(unique(gpos_all))
    read_name_unique <- unique(resL$read_name)
    modmat <- SparseArray::SparseArray(Matrix::sparseMatrix(
            i = GenomicRanges::match(gpos_all, gpos),
            j = match(resL$read_name, read_name_unique),
            x = resL$call_prob,
            dims = c(length(gpos), length(read_name_unique)),
            dimnames = list(NULL, read_name_unique)))

    # create SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(mod_prob = modmat),
        rowRanges = gpos,
        colData = S4Vectors::DataFrame(
            row.names = read_name_unique
        ),
        metadata = list()
    )
    rownames(se) <- paste0(
        GenomeInfoDb::seqnames(SummarizedExperiment::rowRanges(se)),
        ":", BiocGenerics::pos(SummarizedExperiment::rowRanges(se)), ":",
        BiocGenerics::strand(SummarizedExperiment::rowRanges(se)))
    colnames(se) <- rownames(SummarizedExperiment::colData(se))

    se
}