#' Read base modifications from bam file(s)
#'
#' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
#' section 1.7) and return a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with information on modified bases. Implicitly called bases will get
#' a modification probability of zero. Secondary and supplementary alignments
#' are ignored.
#'
#' @param bamfiles Character vector with one or several paths of \code{modBAM}
#'     files, containing information about base modifications in \code{MM} and
#'     \code{ML} tags. If \code{bamfiles} is a named vector, the names are used
#'     as sample names and prefixes for read names. Otherwise, the prefixes will
#'     be \code{s1}, ..., \code{sN}, where \code{N} is the length of
#'     \code{bamfiles}. All \code{bamfiles} must have an index.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic regions to extract the reads from. Alternatively, regions can be
#'     specified as a character vector (e.g. "chr1:1200-1300", "chr2:-6000" or
#'     "chrM") that will be coerced into a \code{GRanges} object. If the end
#'     coordinate is not provided (for example in "chr1:10-", which means "to
#'     the end of chr1"), it will be obtained from \code{seqinfo} or default
#'     to a large value if \code{seqinfo} is not provided. Note that
#'     unless \code{trim=TRUE}, the reads are not trimmed to the boundaries of
#'     the specified ranges. As a result, returned positions will typically
#'     extend out of the specified regions. If \code{nAlnsToSample} is set to a
#'     non-zero value, \code{regions} is ignored.
#' @param modbase Character vector defining the modified base for each sample.
#'     If \code{modbase} is a named vector, the names should correspond to
#'     the names of \code{bamfiles}. Otherwise, it will be assumed that the
#'     elements are in the same order as the files in \code{bamfiles}. If
#'     \code{modbase} has length 1, the same modified base will be used for
#'     all samples.
#' @param level Character scalar specifying the level of returned modification
#'     data. Supported values are:
#'     \describe{
#'         \item{"read"}{: Extracts modification probabilities for individual
#'             reads into an assay called \code{"mod_prob"}. This is the
#'             default if \code{nAlnsToSample} is non-zero or
#'             \code{variantPositions} is not \code{NULL}.}
#'         \item{"quickread"}{: Like "read", but runs faster, does not support
#'             read sampling, and does not return all annotations currently
#'             provided by "read". This is the default if \code{nAlnsToSample}
#'             is zero and \code{variantPositions} is \code{NULL}.}
#'         \item{"summary"}{: Counts the total and modified bases for each
#'             position and strand and returns them in assays named
#'             \code{"Nvalid"} and \code{"Nmod"}, respectively.}
#'     }
#' @param sampleAnnot A \code{data.frame} (or \code{NULL}) providing annotations
#'     for the samples. It must contain at least one column, named
#'     \code{"sample"}, which must contain all the values of
#'     \code{names(bamfiles)}. The provided annotations will be propagated to
#'     the returned \code{SummarizedExperiment} object.
#' @param nAlnsToSample A numeric scalar. If non-zero, \code{regions} is ignored
#'     and approximately \code{nAlnsToSample} randomly selected alignments on
#'     \code{seqnamesToSampleFrom} are read from each of the \code{bamfiles}.
#'     In order to make the results reproducible, make sure to set the
#'     \code{RNGseed} argument in the provided \code{BPPARAM} object (see
#'     below). Please note that secondary and supplementary alignments in
#'     \code{bamfiles} contribute to the total number of alignments but will not
#'     be sampled, thus the number of returned alignments may be lower than
#'     \code{nAlnsToSample}.
#' @param seqnamesToSampleFrom A character vector with one or several sequence
#'     names (chromosomes) from which to sample alignments from (only used if
#'     \code{nAlnsToSample} is greater than zero).
#' @param seqinfo \code{NULL} or a \code{\link[GenomeInfoDb]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Useful to set the sorting order of sequence names.
#' @param sequenceContextWidth,sequenceReference Define the sequence
#'     context to be extracted around modified bases. By default (
#'     \code{sequenceContextWidth = 0}), no sequence context will be
#'     extracted, otherwise it will be returned in \code{rowData(x)$sequenceContext}.
#'     See \code{\link{addSeqContext}} for details.
#' @param variantPositions An optional \code{GPos} object with seqnames and
#'     coordinates of single nucleotide variant positions, to be used to
#'     construct read labels for allele-specific analysis. Ignored if \code{NULL}
#'     or \code{nAlnsToSample > 0} (sampling-mode).
#' @param modProbThreshold A numeric scalar, indicating the modification
#'     probability threshold to use to classify a base as 'modified' or
#'     'unmodified'. Only used if \code{level} is \code{"summary"}.
#' @param trim A logical scalar. If \code{TRUE}, the returned
#'     \code{SummarizedExperiment} object will only contain the positions
#'     overlapping the specified \code{regions}. If \code{FALSE} (default),
#'     the object will be extended to all positions covered by the reads
#'     overlapping \code{regions}. In both cases, only reads overlapping
#'     the specified \code{regions} are included.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use for some of the steps
#'     in \code{readModBam()}. The default value is
#'     (\code{\link[BiocParallel]{MulticoreParam}(4L, RNGseed = 42L)}).
#'     If randomly sampling reads (\code{nAlnsToSample > 0}), make sure to set
#'     the \code{RNGseed} argument when constructing the \code{BPPARAM} object
#'     for reproducible results (see also
#'     \code{vignette("Random_Numbers", package = "BiocParallel")}).
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
#'            modbase = "a", verbose = TRUE,
#'            BPPARAM = BiocParallel::SerialParam())
#'
#' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
#'     SAM ML and MM tags for base modifications.
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match seqnames start end
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics do.call cbind pos strand sort
#' @importFrom BiocParallel bplapply MulticoreParam bpnworkers bpworkers<-
#' @importFrom methods is
#' @importFrom cli cli_abort cli_warn
#'
#' @export
readModBam <- function(bamfiles,
                       regions = NULL,
                       modbase,
                       level = ifelse(nAlnsToSample == 0 & is.null(variantPositions),
                                      "quickread", "read"),
                       sampleAnnot = NULL,
                       nAlnsToSample = 0,
                       seqnamesToSampleFrom = "chr19",
                       seqinfo = NULL,
                       sequenceContextWidth = 0,
                       sequenceReference = NULL,
                       variantPositions = NULL,
                       modProbThreshold = 0.5,
                       trim = FALSE,
                       BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                       verbose = FALSE) {
    # digest arguments
    .assertVector(x = bamfiles, type = "character")
    if (length(bamfiles) == 0) {
        cli_abort("{.arg bamfiles} must not be an empty vector")
    }
    if (any(i <- !file.exists(bamfiles))) {
        cli_abort("not all {.arg bamfiles} exist: {.file {bamfiles[i]}}")
    }
    if (is.null(names(bamfiles))) {
        names(bamfiles) <- paste0("s", seq_along(bamfiles))
    } else if (any(duplicated(names(bamfiles)))) {
        cli_abort("{.code names(bamfiles)} are not unique")
    }
    .assertScalar(x = level, type = "character",
                  validValues = c("read", "summary", "quickread"))
    .assertVector(x = sampleAnnot, type = "data.frame", allowNULL = TRUE)
    if (!is.null(sampleAnnot)) {
        if (!("sample" %in% colnames(sampleAnnot))) {
            cli_abort("{.arg sampleAnnot} must have at least a column named 'sample'")
        }
        if (!all(names(bamfiles) %in% sampleAnnot$sample)) {
            cli_abort(paste0(
                "Annotation information missing for some samples: ",
                 "{setdiff(names(bamfiles), sampleAnnot$sample)}"))
        }
    }
    if (is.character(regions)) {
        regions <- .regionStringToGRanges(regions = regions,
                                          seqinfo = seqinfo)
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
            cli_abort("names of {.arg modbase} and {.arg bamfiles} don't agree")
        }
    }
    # for valid values of `modbase`, see
    # https://samtools.github.io/hts-specs/SAMtags.pdf (section 1.7)
    if (any(i <- !modbase %in% c("m","h","f","c","C","g","e","b","T",
                                 "U","a","A","o","G","n","N"))) {
        cli_abort("invalid {.arg modbase} values: {unique(modbase[i])}")
    }
    .assertScalar(x = nAlnsToSample, type = "numeric", rngIncl = c(0, Inf))
    if (nAlnsToSample > 0 && level %in% c("summary", "quickread")) {
        cli_abort(paste0("Read sampling is not supported if {.arg level} is set ",
                         "to 'summary' or 'quickread'"))
    }
    if (nAlnsToSample > 0) {
        if (length(regions) > 0) {
            cli_warn("Ignoring {.arg regions} because {.arg nAlnsToSample} is greater than zero")
        }
        regions <- GRanges()
        if (!is.null(variantPositions)) {
            cli_warn("Ignoring {.arg variantPositions} because {.arg nAlnsToSample} is greater than zero")
        }
        variantPositions <- NULL
    } else {
        if (length(regions) == 0) {
            cli_abort("{.arg regions} must contain at least one genomic range if not in sampling mode")
        }
    }
    .assertScalar(x = modProbThreshold, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = seqnamesToSampleFrom, type = "character")
    if (!is.null(seqinfo)) {
        if (!is(seqinfo, "Seqinfo") &&
            (!is.numeric(seqinfo) || is.null(names(seqinfo)))) {
            cli_abort(paste0(
                "{.arg seqinfo} must be {.code NULL}, a {.cls Seqinfo} object ",
                "or a named {.cls numeric} vector with genomic sequence lengths."))
        }
    }
    .assertScalar(x = sequenceContextWidth, type = "numeric", rngIncl = c(0, 1000))
    .assertVector(x = variantPositions, type = "GPos", allowNULL = TRUE)
    .assertScalar(x = trim, type = "logical")
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")

    # sort and subset variantPositions
    if (length(variantPositions) > 0) {
        variantPositions <- sort(
            sort(subsetByOverlaps(x = variantPositions,
                                  ranges = regions,
                                  ignore.strand = TRUE)),
            ignore.strand = TRUE)
        variantRefNames <- as.character(seqnames(variantPositions))
        # make coordinates zero-based
        variantRefPositions <- pos(variantPositions) - 1L
    } else {
        variantRefNames <- character(0L)
        variantRefPositions <- integer(0L)
    }

    # determine the number of parallel threads to be used for
    # bam files (preferred) and decompression of bam records (if available)
    # (accept some level of over-subscription)
    ncpuTotal <- bpnworkers(BPPARAM)
    if (is(BPPARAM, "MulticoreParam") || is(BPPARAM, "SnowParam")) {
        ncpuFiles <- min(ncpuTotal, length(bamfiles))
        oversubscriptionRate <- 2.0
        ncpuDecompression <- min(8L, max(1L, as.integer(
            floor(oversubscriptionRate * ncpuTotal / ncpuFiles))))
        bpworkers(BPPARAM) <- ncpuFiles
        on.exit(bpworkers(BPPARAM) <- ncpuTotal)
    } else {
        ncpuDecompression <- 1L
    }

    # extract modification probabilities from `bamfiles`
    .message("extracting base modifications from modBAM files", noTimer = TRUE)
    # remark: unsafe to use as.character(GRanges) here, as a length-1 range
    #         would become e.g. "chr1:35000" (no end coordinate), which is
    #         interpreted by htslib as: "read all alignments overlapping chr1:35000-END_OF_chr1"
    regions_str <- paste0(seqnames(regions), ":", start(regions), "-", end(regions))
    resLL <- bplapply(
        structure(names(bamfiles), names = names(bamfiles)),
        function(nm,
                 mylevel = level,
                 bamf = bamfiles[nm],
                 myregions_str = regions_str,
                 mymodbase = modbase[nm],
                 mymodProbThreshold = modProbThreshold,
                 mynAlnsToSample = nAlnsToSample,
                 myseqnamesToSampleFrom = seqnamesToSampleFrom,
                 myvariantRefNames = variantRefNames,
                 myvariantRefPositions = variantRefPositions,
                 myncpuDecompression = ncpuDecompression,
                 myverbose = verbose) {

            if (mylevel == "read") {
                # extract modifications (returned list is similar to modkit extract
                # output, see https://nanoporetech.github.io/modkit/intro_extract.html)
                resL <- read_modbam_cpp(inname_str = bamf,
                                        regions = myregions_str,
                                        modbase = mymodbase,
                                        n_alns_to_sample = as.integer(mynAlnsToSample),
                                        tnames_for_sampling = myseqnamesToSampleFrom,
                                        variantRefNames = myvariantRefNames,
                                        variantRefPositions = as.integer(myvariantRefPositions),
                                        n_threads = as.integer(myncpuDecompression),
                                        verbose = myverbose)

                # convert inferred `mod_prob` to zero. Inferred means that the
                # modification was omitted from the BAM file, e.g. DORADO omits base
                # modification probabilities less than 0.05, and read_modbam_cpp returns
                # a mod_prob of -1 for these.
                resL$mod_prob[resL$mod_prob == -1] <- 0
            } else if (mylevel == "summary") {
                resL <- pileup_modbam_cpp(inname_str = bamf,
                                          regions = myregions_str,
                                          modbase = mymodbase,
                                          level = "summary",
                                          mod_prob_thresh = mymodProbThreshold,
                                          n_threads = as.integer(myncpuDecompression),
                                          verbose = myverbose)
            } else if (mylevel == "quickread") {
                resL <- pileup_modbam_cpp(inname_str = bamf,
                                          regions = myregions_str,
                                          modbase = mymodbase,
                                          level = "read",
                                          mod_prob_thresh = mymodProbThreshold,
                                          n_threads = as.integer(myncpuDecompression),
                                          verbose = myverbose)
                resL$mod_prob[resL$mod_prob == -1] <- 0
                resL$read_df <- resL$read_df[resL$read_df$read_id %in% resL$read_id, ]
            }
            resL
    }, BPPARAM = BPPARAM)

    # create GPos objects for each input
    gposL <- bplapply(resLL, function(resL, myseqinfo = seqinfo) {
        GenomicRanges::GPos(seqnames = resL$chrom, pos = resL$ref_position,
                            strand = resL$ref_mod_strand, seqinfo = myseqinfo)
    }, BPPARAM = BPPARAM)

    # create combined GPos, reduce to unique positions
    .message("finding unique genomic positions...")
    gpos <- sort(sort(unique(do.call(c, unname(gposL)))), ignore.strand = TRUE)
    .message("collapsed {sum(lengths(gposL))} positions to {length(gpos)} unique ones")

    # if trim=TRUE, trim GPos to only the indicated region
    if (trim) {
        gpos <- subsetByOverlaps(gpos, regions)
    }

    # add sequence context
    if (sequenceContextWidth > 0) {
        .message("extracting sequence contexts")
        mcols(gpos)$sequenceContext <- extractSeqContext(
            x = as(gpos, "GRanges"),
            sequenceContextWidth = sequenceContextWidth,
            sequenceReference = sequenceReference)
    }

    if (level %in% c("read", "quickread")) {
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
                # if trim=TRUE, not all positions in gposL may be present in gpos
                found <- which(!is.na(i))
                j <- match(x$read_id, readL[[nm]])
                namat[cbind(i[found], j[found])] <- x$mod_prob[found]
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
                                           variant_label = character(0),
                                           aligned_fraction = numeric(0))
            }
        }
    } else {
        Nmod <- matrix(0, nrow = length(gpos), ncol = length(bamfiles))
        Nvalid <- matrix(0, nrow = length(gpos), ncol = length(bamfiles))
        for (i in seq_along(bamfiles)) {
            nm <- names(bamfiles)[i]
            # if trim=TRUE, not all positions in gposL may be present in gpos
            j <- match(gposL[[nm]], gpos)
            found <- which(!is.na(j))
            Nmod[j[found], i] <- resLL[[nm]]$Nmod[found]
            Nvalid[j[found], i] <- resLL[[nm]]$Nvalid[found]
            FracMod <- Nmod / Nvalid
        }
    }

    # create SummarizedExperiment object
    cdata <- DataFrame(
        row.names = names(bamfiles),
        sample = names(bamfiles),
        modbase = modbase[names(bamfiles)]
    )
    if (level %in% c("read", "quickread")) {
        cdata <- cbind(
            cdata,
            DataFrame(n_reads = unlist(lapply(readdfL, nrow), use.names = FALSE),
                      readInfo = readdfL)
        )
    }
    if (!is.null(sampleAnnot) && any(colnames(sampleAnnot) != "sample")) {
        sampleAnnot <- sampleAnnot[match(cdata$sample, sampleAnnot$sample),
                                   colnames(sampleAnnot) != "sample",
                                   drop = FALSE]
        cdata <- cbind(cdata, sampleAnnot)
    }
    if (level %in% c("read", "quickread")) {
        stopifnot(names(modmat) == cdata$sample)
        se <- SummarizedExperiment(
            assays = list(mod_prob = modmat),
            rowRanges = gpos,
            colData = cdata,
            metadata = list(readLevelData = list(assayNames = "mod_prob",
                                                 colDataColumns = "readInfo"),
                            variantPositions = variantPositions)
        )
    } else {
        stopifnot(colnames(Nmod) == cdata$sample,
                  colnames(Nvalid) == cdata$sample,
                  colnames(FracMod) == cdata$sample)
        se <- SummarizedExperiment(
            assays = list(Nmod = Nmod,
                          Nvalid = Nvalid,
                          FracMod = FracMod),
            rowRanges = gpos,
            colData = cdata,
            metadata = list(readLevelData = list(assayNames = character(0),
                                                 colDataColumns = character(0)))
        )
    }
    if (nrow(se) > 0) {
        rownames(se) <- paste0(
            seqnames(rowRanges(se)), ":", pos(rowRanges(se)), ":",
            strand(rowRanges(se)))
        colnames(se) <- rownames(colData(se))
    }

    # Remove reads with all NA values
    if (level %in% c("read", "quickread")) {
        se <- filterReads(se, readInfoCol = NULL, qcCol = NULL, prune = FALSE)
    }
    se
}
