#' Read base modifications from bam file(s).
#'
#' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
#' section 1.7) and return a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with information on modified bases.
#'
#' @param bamfiles Character vector with one or several paths of \code{modBAM}
#'     files, containing information about base modifications in \code{MM} and
#'     \code{ML} tags. If \code{bamfiles} is a named vector, the names are used
#'     as sample names and prefixes for read names. Otherwise, the prefixes will
#'     be \code{s1}, ..., \code{sN}, where \code{N} is the length of
#'     \code{bamfiles}. All \code{bamfiles} must have an index.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic regions to extract the reads from. Alternatively, regions can be
#'     specified as a character vector (e.g. "chr1:1200-1300") that can be
#'     coerced into a \code{GRanges} object. Note that the reads are not
#'     trimmed to the boundaries of the specified ranges. As a result, returned
#'     positions will typically extend out of the specified regions.
#'     If \code{nAlnsToSample} is set to a non-zero value, \code{regions} is
#'     ignored.
#' @param modbase Character vector defining the modified base for each sample.
#'     If \code{modbase} is a named vector, the names should correspond to
#'     the names of \code{bamfiles}. Otherwise, it will be assumed that the
#'     elements are in the same order as the files in \code{bamfiles}. If
#'     \code{modbase} has length 1, the same modified base will be used for
#'     all samples.
#' @param nAlnsToSample A numeric scalar. If non-zero, \code{regions} is ignored
#'     and approximately \code{nAlnsToSample} randomly selected alignments on
#'     \code{seqnamesToSampleFrom} are read from each of the \code{bamfiles}.
#' @param seqnamesToSampleFrom A character vector with one or several sequence
#'     names (chromosomes) from which to sample alignments from (only used if
#'     `nAlnsToSample` is greater than zero).
#' @param seqinfo \code{NULL} or a \code{\link[GenomeInfoDb]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Useful to set the sorting order of sequence names.
#' @param sequence.context.width,sequence.reference Define the sequence
#'     context to be extracted around modified bases. By default (
#'     \code{sequence.context.width = 0}), no sequence context will be
#'     extracted, otherwise it will be returned in \code{rowData(x)$sequence.context}.
#'     See \code{\link{addSeqContext}} for details.
#' @param variantPositions An optional \code{GPos} object with seqnames and
#'     coordinates of single nucleotide variant positions, to be used to
#'     construct read labels for allele-specific analysis. Ignored if \code{NULL}
#'     or \code{nAlnsToSample > 0} (sampling-mode).
#' @param ncpu A numeric scalar giving the number of parallel CPU threads to
#'     to use for some of the steps in \code{readModBam} (e.g. the number of
#'     bam files to process in parallel).
#' @param ncpuDecompression A numeric scalar giving the number of parallel CPU
#'     threads to use for decompressing bam records. When reading from multiple
#'     bam files in parallel (\code{ncpu > 1}), the total number of threads
#'     used may be as many as \code{ncpu * ncpuDecompression}.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with genomic positions in rows and samples in columns. The assay
#'     \code{"mod_prob"} contains per-read modification probabilities,
#'     with each column (sample) corresponding to a position-by-read
#'     \code{\link[SparseArray]{NaMatrix}}.
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'            modbase = "a", verbose = TRUE)
#'
#' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
#'     SAM ML and MM tags for base modifications.
#'
#' @author Michael Stadler
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos sort match subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics do.call cbind pos strand
#' @importFrom parallel mclapply
#'
#' @export
readModBam <- function(bamfiles,
                       regions = NULL,
                       modbase,
                       nAlnsToSample = 0,
                       seqnamesToSampleFrom = "chr19",
                       seqinfo = NULL,
                       sequence.context.width = 0,
                       sequence.reference = NULL,
                       variantPositions = NULL,
                       ncpu = 1L,
                       ncpuDecompression = 2L,
                       verbose = FALSE) {
    # digest arguments
    .assertVector(x = bamfiles, type = "character")
    if (any(i <- !file.exists(bamfiles))) {
        stop("not all `bamfiles` exist: ", paste(bamfiles[i], collapse = ", "))
    }
    if (is.null(names(bamfiles))) {
        names(bamfiles) <- paste0("s", seq_along(bamfiles))
    } else if (any(duplicated(names(bamfiles)))) {
        stop("`names(bamfiles)` are not unique")
    }
    if (is.character(regions)) {
        regions <- as(regions, "GRanges")
    }
    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)
    if (length(modbase) == 1) {
        modbase <- rep(modbase, length(bamfiles))
    }
    .assertVector(x = modbase, type = "character", len = length(bamfiles))
    if (is.null(names(modbase))) {
        names(modbase) <- names(bamfiles)
    } else {
        if (!all(names(modbase) %in% names(bamfiles))) {
            stop("names of `modbase` and `bamfiles` don't agree")
        }
    }
    # for valid values of `modbase`, see
    # https://samtools.github.io/hts-specs/SAMtags.pdf (section 1.7)
    if (any(i <- !modbase %in% c("m","h","f","c","C","g","e","b","T",
                                 "U","a","A","o","G","n","N"))) {
        stop("invalid `modbase` values: ",
             paste(unique(modbase[i]), collapse = ", "))
    }
    .assertScalar(x = nAlnsToSample, type = "numeric", rngIncl = c(0, Inf))
    if (nAlnsToSample > 0) {
        if (length(regions) > 0) {
            warning("Ignoring `regions` because `nAlnsToSample` is greater than zero")
        }
        regions <- GRanges()
        if (!is.null(variantPositions)) {
            warning("Ignoring `variantPositions` because `nAlnsToSample` is greater than zero")
        }
        variantPositions <- NULL
    } else {
        if (length(regions) == 0) {
            stop("`regions` must contain at least one genomic range if not in sampling mode")
        }
    }
    .assertVector(x = seqnamesToSampleFrom, type = "character")
    if (!is.null(seqinfo)) {
        if (!is(seqinfo, "Seqinfo") &&
            (!is.numeric(seqinfo) || is.null(names(seqinfo)))) {
            stop("`seqinfo` must be `NULL`, a `Seqinfo` object or a named",
                 " numeric vector with genomic sequence lengths.")
        }
    }
    .assertScalar(x = sequence.context.width, type = "numeric", rngIncl = c(0, 1000))
    .assertVector(x = variantPositions, type = "GPos", allowNULL = TRUE)
    .assertScalar(x = ncpu, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = ncpuDecompression, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = verbose, type = "logical")

    # sort and subset variantPositions
    if (length(variantPositions) > 0) {
        variantPositions <- sort(subsetByOverlaps(x = variantPositions,
                                                  ranges = regions,
                                                  ignore.strand = TRUE))
        variantRefNames <- as.character(seqnames(variantPositions))
        variantRefPositions <- pos(variantPositions)
    } else {
        variantRefNames <- character(0L)
        variantRefPositions <- integer(0L)
    }

    # extract modification probabilities from `bamfiles`
    .message("extracting base modifications from modBAM files")
    regions_str <- as.character(regions, ignore.strand = TRUE)
    resLL <- mclapply(structure(names(bamfiles), names = names(bamfiles)),
                                function(nm) {
        # extract modifications (returned list is similar to modkit extract
        # output, see https://nanoporetech.github.io/modkit/intro_extract.html)
        resL <- read_modbam_cpp(inname_str = bamfiles[nm],
                                regions = regions_str,
                                modbase = modbase[nm],
                                n_alns_to_sample = as.integer(nAlnsToSample),
                                tnames_for_sampling = seqnamesToSampleFrom,
                                variantRefNames = variantRefNames,
                                variantRefPositions = as.integer(variantRefPositions),
                                n_threads = as.integer(ncpuDecompression),
                                verbose = verbose)

        # convert 0-based ref_position to 1-based
        resL$ref_position <- resL$ref_position + 1L

        # convert inferred `mod_prob` to our minimal value as in
        # readModkitExtract(). Inferred means that the modification was omitted
        # from the BAM file, e.g. DORADO omits base modification probabilities
        # less than 0.05, and read_modbam_cpp returns a mod_prob of -1 for
        # these.
        resL$mod_prob[resL$mod_prob == -1] <- 0.02
        resL
    }, mc.cores = ncpu)

    # create GPos objects for each input
    gposL <- mclapply(resLL, function(resL) {
        GPos(seqnames = resL$chrom, pos = resL$ref_position,
             strand = resL$ref_mod_strand, seqinfo = seqinfo)
    }, mc.cores = ncpu)

    # create combined GPos, reduce to unique positions
    .message("finding unique genomic positions...")
    gpos <- sort(unique(do.call(c, unname(gposL))))
    .message("collapsed {sum(lengths(gposL))} positions to {length(gpos)} unique ones")

    # add sequence context
    if (sequence.context.width > 0) {
        .message("extracting sequence contexts")
        mcols(gpos)$sequence.context <- extractSeqContext(
            x = as(gpos, "GRanges"),
            sequence.context.width = sequence.context.width,
            sequence.reference = sequence.reference)
    }

    # extract unique read names
    readL <- lapply(resLL, function(resL) resL$read_df$read_id)

    # modified probability
    modmat <- make_zero_col_DFrame(nrow = length(gpos))
    readdfL <- SimpleList()
    for (nm in names(bamfiles)) {
        x <- resLL[[nm]]
        if (length(x$read_id) > 0) {
            namat <- NaArray(dim = c(length(gpos), length(readL[[nm]])),
                             dimnames = list(NULL, paste0(nm, "-", readL[[nm]])),
                             type = "double")
            i <- match(gposL[[nm]], gpos)
            j <- match(x$read_id, readL[[nm]])
            namat[cbind(i, j)] <- x$mod_prob
            modmat[[nm]] <- namat
            rownames(x$read_df) <- paste0(nm, "-", x$read_df$read_id)
            x$read_df$read_id <- NULL
            x$read_df$aligned_fraction <- x$read_df$aligned_length / x$read_df$read_length
            readdfL[[nm]] <- DataFrame(x$read_df)
        } else {
            modmat[[nm]] <- NaArray(dim = c(length(gpos), 0), type = "double")
            readdfL[[nm]] <- DataFrame(qscore = numeric(0),
                                       read_length = integer(0),
                                       aligned_length = integer(0),
                                       aligned_fraction = numeric(0),
                                       variant_label = character(0))
        }
    }

    # create SummarizedExperiment object
    se <- SummarizedExperiment(
        assays = list(mod_prob = modmat),
        rowRanges = gpos,
        colData = DataFrame(
            row.names = names(bamfiles),
            sample = names(bamfiles),
            modbase = modbase[names(bamfiles)],
            n_reads = unlist(lapply(readdfL, nrow), use.names = FALSE),
            read_info = readdfL
        ),
        metadata = list(readLevelData = list(assayNames = "mod_prob",
                                             colDataColumns = "read_info"),
                        variantPositions = variantPositions)
    )
    if (nrow(se) > 0) {
        rownames(se) <- paste0(
            seqnames(rowRanges(se)), ":", pos(rowRanges(se)), ":",
            strand(rowRanges(se)))
        colnames(se) <- rownames(colData(se))
    }

    se
}
