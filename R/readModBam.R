#' Read base modifications from a bam file.
#'
#' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
#' section 1.7) and return a list of information on modified bases.
#'
#' @param bamfiles Character vector with one or several paths of \code{modBAM}
#'     files, containing information about base modifications in \code{MM} and
#'     \code{ML} tags. If \code{bamfiles} is a named vector, the names are used
#'     as sample names and prefixes for the column (read) names in the returned
#'     \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#'     Otherwise, the prefixes will be \code{s1}, ..., \code{sN}, where \code{N}
#'     is the length of \code{bamfiles}. All \code{bamfiles} must have an index.
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
#' @param ncpu A numeric scalar giving the number of parallel CPU threads to
#'     to use for some of the steps in \code{readModBam}.
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
#' @importFrom BiocGenerics do.call cbind pos strand
#' @importFrom parallel mclapply
#'
#' @export
readModBam <- function(bamfiles,
                       regions,
                       modbase,
                       seqinfo = NULL,
                       ncpu = 1L,
                       verbose = FALSE) {
    # digest arguments
    .assertVector(x = bamfiles, type = "character")
    if (any(i <- !file.exists(bamfiles))) {
        stop("not all `bamfiles` exist: ", paste(bamfiles[i], collapse = ", "))
    }
    if (is.character(regions)) {
        regions <- as(regions, "GRanges")
    }
    .assertScalar(x = regions, type = "GRanges")
    # for valid values of `modbase`, see
    # https://samtools.github.io/hts-specs/SAMtags.pdf (section 1.7)
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
    .assertScalar(x = ncpu, type = "numeric")
    .assertScalar(x = verbose, type = "logical")

    # make sure bamfiles has sample names
    if (is.null(names(bamfiles))) {
        names(bamfiles) <- paste0("s", seq_along(bamfiles))
    }

    # extract modification probabilities from `bamfiles`
    if (verbose) {
        message("extracting base modifications from modBAM files")
    }
    regions_str <- as.character(regions, ignore.strand = TRUE)
    resLL <- parallel::mclapply(bamfiles, function(bamfile) {
        # extract modifications (returned list is similar to modkit extract
        # output, see https://nanoporetech.github.io/modkit/intro_extract.html)
        resL <- read_modbam_cpp(inname_str = bamfile,
                                regions = regions_str,
                                modbase = modbase,
                                verbose = verbose)

        # convert 0-based ref_position to 1-based
        resL$ref_position <- resL$ref_position + 1L

        # convert inferred `mod_prob` to our minimal value as in
        # readModkitExtract(). Inferred means that the modification was omitted
        # from the BAM file, e.g. DORADO omits base modification probabilities
        # less than 0.05, and read_modbam_cpp returns a call_probability of -1 for
        # these.
        resL$mod_prob[resL$mod_prob == -1] <- 0.02
        resL
    }, mc.cores = ncpu)

    # create GPos objects for each input
    gposL <- parallel::mclapply(resLL, function(resL) {
        GenomicRanges::GPos(seqnames = resL$chrom, pos = resL$ref_position,
                            strand = resL$ref_strand, seqinfo = seqinfo)
    }, mc.cores = ncpu)

    # create combined GPos, reduce to unique positions
    if (verbose) {
        message("finding unique genomic positions...", appendLF = FALSE)
    }
    gpos <- GenomicRanges::sort(unique(do.call(c, unname(gposL))))
    if (verbose) {
        message("collapsed ", sum(lengths(gposL)), " positions to ",
                length(gpos), " unique ones")
    }

    # extract read names
    readL <- lapply(resLL, function(resL) unique(resL$read_id))

    # modified probability
    modmat <- S4Vectors::make_zero_col_DFrame(nrow = length(gpos))
    for (nm in names(bamfiles)) {
        x <- resLL[[nm]]
        # only record observed values
        modmat[[nm]] <- SparseArray::SparseArray(Matrix::sparseMatrix(
            i = GenomicRanges::match(gposL[[nm]], gpos),
            j = match(x$read_id, readL[[nm]]),
            x = x$mod_prob,
            dims = c(length(gpos), length(readL[[nm]])),
            dimnames = list(NULL, paste0(nm, "-", readL[[nm]]))
        ))
    }
    # reduce to a single sparse matrix
    readnames <- do.call(c, lapply(modmat, colnames))
    samplenames <- rep(names(modmat), vapply(modmat, ncol, 0))
    modmat <- BiocGenerics::do.call(BiocGenerics::cbind, modmat)

    # create SummarizedExperiment object
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(mod_prob = modmat),
        rowRanges = gpos,
        colData = S4Vectors::DataFrame(
            row.names = readnames,
            sample = samplenames
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
