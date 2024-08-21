#' Extract the sequence context around positions of interest.
#'
#' @description
#' This function will extract a sequence context of \code{sequence.context.width}
#' bases around the center of the regions defined in \code{x} from
#' \code{sequence.reference}.
#'
#' @param x A \code{\link[GenomicRanges]{GRanges}} object defining the regions
#'     of interest. The extracted sequences will correspond to the regions
#'     defined as \code{resize(x, width = sequence.context.width, fix = "center"}.
#' @param sequence.context.width A numeric scalar giving the width of the
#'     sequence context to be extracted from the reference
#'     (\code{sequence.reference} argument). This must be an odd number
#'     so that the sequence can be centered on the modified base.
#'     If \code{sequence.context.width = 0} (the default), no
#'     sequence context will be extracted.
#' @param sequence.reference A \code{\link[BSgenome]{BSgenome}} object, or a
#'     character scalar giving the path to a fasta formatted file with reference
#'     sequences, or a \code{\link[Biostrings]{DNAStringSet}} object.
#'     The sequence context (see \code{sequence.context.width} argument) will be
#'     extracted from these sequences.
#'
#' @return A \code{\link[Biostrings]{DNAStringSet}} object of the same length
#'     as \code{x} with extracted sequence context. All elements are guaranteed
#'     to have identical length (if the sequence context extends to before the
#'     start or beyond the end of a  reference sequence, it will be padded with
#'     'N' bases.
#'
#' @author Michael Stadler
#'
#' @examples
#' # file with sequence in fasta format of length 6957060
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' # define some regions at the end of the reference sequence
#' regions <- GenomicRanges::GRanges(
#'     "chr1", IRanges::IRanges(start = 6957060 - c(4, 2, 0),
#'     width = 1, names = c("a","b","c")))
#'
#' # extract sequence context (note the padding with N's)
#' seqContext(regions, 7, reffile)
#'
#' @seealso \code{\link[GenomicRanges]{resize}}, \code{\link[Biostrings]{DNAStringSet}}
#'
#' @importFrom GenomicRanges GRanges resize trim
#' @importFrom GenomeInfoDb seqlengths seqlengths<-
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom BSgenome getSeq
#' @importFrom methods as is
#'
#' @export
extractSeqContext <- function(x,
                       sequence.context.width,
                       sequence.reference) {
    # digest arguments
    .assertVector(x = x, type = "GRanges")
    .assertScalar(x = sequence.context.width, type = "numeric", rngIncl = c(1, 1000))
    if (sequence.context.width %% 2 == 0) {
        sequence.context.width <- sequence.context.width + 1
        warning("`sequence.context.width` was increased to ", sequence.context.width,
                " (must be an odd number)")
    }
    if (!is(sequence.reference, "BSgenome") &&
        !is(sequence.reference, "DNAStringSet") &&
        !(is.character(sequence.reference) && file.exists(sequence.reference))) {
        stop("`sequence.reference` must be either a BSgenome object, ",
             "a DNAStringSet object, or a path to a fasta file.")
    }

    # resize x
    xcontext <- GenomicRanges::resize(x, width = sequence.context.width, fix = "center")

    # obtain reference sequences
    if (is.character(sequence.reference)) {
        ref <- Biostrings::readDNAStringSet(sequence.reference)
        names(ref) <- sub(" .*$", "", names(ref))
    } else {
        ref <- sequence.reference
    }
    # seqlengths(xcontext) <- seqlengths(ref)

    # extract sequences
    Npre <- pmax(0L, 1L - start(xcontext))
    Npost <- pmax(0L, end(xcontext) - seqlengths(ref)[as.character(seqnames(xcontext))])
    if (any(Npre > 0) || any(Npost > 0)) {
        suppressWarnings(seqlengths(xcontext) <- seqlengths(ref))
        xcontext <- GenomicRanges::trim(xcontext)
        seqcontext <- DNAStringSet(
            x = paste0(strrep("N", Npre),
                       getSeq(ref, xcontext),
                       strrep("N", Npost)),
            use.names = FALSE)
    } else {
        seqcontext <- getSeq(ref, xcontext)
    }

    # return results
    if (!is.null(names(x))) {
        names(seqcontext) <- names(x)
    }
    return(seqcontext)
}
