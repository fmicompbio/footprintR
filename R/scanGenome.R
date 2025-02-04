#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @keywords internal
#' @noRd
.tileChromosome <- function(tileSize,
                            windowSize,
                            windowStep,
                            chromName,
                            chromLength) {
    # tile a chromosome
    nWindowsPerTile <- floor(tileSize / windowStep)
    tileSize <- ((nWindowsPerTile - 1) * windowStep) + windowSize
    nTiles <- ceiling((chromLength - (windowSize - windowStep)) / (tileSize - (windowSize - windowStep)))
    s <- 1 + (seq.int(nTiles) - 1) * (tileSize - windowSize + windowStep)
    regs <- GRanges(chromName, IRanges(start = s, width = tileSize))
    regs
}

#' Generate counts for sequential windows in a single region
#'
#' Read modification data from \code{bamfiles} for a chunk of the genome
#' defined by \code{region} (by calling \code{\link{readModBam}} and aggregate
#' modified and total counts for windows.
#'
#' @param bamfiles Character vector with one or several modBam file names.
#'     Note that no read-filtering will be performed on the data from
#'     these files, so possibly the files should contain only filtered
#'     alignments.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object specifying which
#'     genomic region to extract the reads from. Alternatively, the regions can
#'     be specified as a character scalar (e.g. "chr1:1200-1300") that can be
#'     coerced into a \code{GRanges} object.
#' @param modbase Character vector defining the modified base to extract.
#' @param modProbThreshold A numeric scalar, indicating the modification
#'     probability threshold to use to classify a base as 'modified' or
#'     'unmodified'.
#' @param sampleAnnot A \code{data.frame} (or \code{NULL}) providing annotations
#'     for the samples. It must contain at least one column, named
#'     \code{"sample"}, which must contain all the values of
#'     \code{names(bamfiles)}. The provided annotations will be propagated to
#'     the returned \code{SummarizedExperiment} object.
#' @param seqinfo \code{NULL} or a \code{\link[GenomeInfoDb]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Useful to set the sorting order of sequence names.
#' @param sequenceContextWidth,sequenceReference Define the sequence
#'     context to be extracted around modified bases. By default (
#'     \code{sequenceContextWidth = 0}), no sequence context will be
#'     extracted, otherwise it will be returned in \code{rowData(x)$sequenceContext}.
#'     See \code{\link{addSeqContext}} for details.
#' @param sequenceContext A character vector with sequence contexts to
#'     retain. To apply this filter, the arguments \code{"sequenceContextWidth"}
#'     and \code{"sequenceReference"} must be set.
#' @param windowMode Character scalar defining how windows in \code{region}
#'     are define. Currently supported are:
#'     \describe{
#'         \item{"fixed"}{: The window is of fixed size (in number of bases),
#'         corresponding to the \code{windowSize} argument value.}
#'     }
#' @param windowSize Numeric scalar defining the size of windows (see
#'     \code{windowMode} argument).
#' @param windowStep Numeric scalar defining the step (shift) between
#'     start positions of consecutive windows. For non-overlapping
#'     consecutive windows, \code{windowStep} should be equal to
#'     \code{windowSize}.
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
#' @author Panagiotis Papasaikas, Sebastien Smallwood, Charlotte Soneson, Michael Stadler
#'
#' @returns A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with columns corresponding to samples (the elements of \code{bamfiles})
#'     and rows corresponding to windows in \code{region}. The object contains
#'     the assays \code{"Nmod"} and \code{"Nvalid"} containing the number of
#'     modified and total (valid) bases in each window and sample, respectively.
#'
#' @examples
#' modbamfiles <- system.file("extdata",
#'                            c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' gnmfasta <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#' quantifyWindowsInRegion(bamfiles = modbamfiles,
#'                         region = "chr1:6940000-6955000",
#'                         modbase = "a",
#'                         sequenceContextWidth = 1,
#'                         sequenceReference = gnmfasta,
#'                         sequenceContext = "A",
#'                         BPPARAM = BiocParallel::SerialParam())
#' quantifyWindowsInRegion(bamfiles = modbamfiles,
#'                         region = "chr1:6940000-6955000",
#'                         modbase = "a",
#'                         BPPARAM = BiocParallel::SerialParam())
#'
#' @importFrom SummarizedExperiment rowRanges colData
#' @importFrom GenomicRanges GRanges GPos
#' @importFrom IRanges IRanges start end findOverlaps
#' @importFrom S4Vectors queryHits subjectHits metadata
#' @importFrom cli cli_abort
#'
#' @export
quantifyWindowsInRegion <- function(bamfiles,
                                    region,
                                    modbase,
                                    modProbThreshold = 0.5,
                                    sampleAnnot = NULL,
                                    seqinfo = NULL,
                                    sequenceContextWidth = 0,
                                    sequenceReference = NULL,
                                    sequenceContext = NULL,
                                    windowMode = "fixed",
                                    windowSize = 24,
                                    windowStep = round(windowSize / 2),
                                    BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                                    verbose = FALSE) {
    # check parameters
    .assertScalar(x = region)
    .assertScalar(x = modbase, type = "character")
    .assertScalar(x = windowMode, type = "character",
                  validValues = c("fixed"))
    .assertScalar(x = windowSize, type = "numeric", rngIncl = c(1L, Inf))

    # read summary-level data
    se <- readModBam(bamfiles = bamfiles, regions = region, modbase = modbase,
                     level = "summary", sampleAnnot = sampleAnnot,
                     seqinfo = seqinfo, sequenceContextWidth = sequenceContextWidth,
                     sequenceReference = sequenceReference,
                     modProbThreshold = modProbThreshold,
                     trim = TRUE, BPPARAM = BPPARAM,
                     verbose = verbose)

    # filter positions
    if (!is.null(sequenceContext)) {
        se <- filterPositions(se = se, filters = "sequenceContext",
                              sequenceContext = sequenceContext,
                              assayNameNA = NULL)
    }

    if (nrow(se) > 0) {
        # define windows for aggregation
        if (identical(windowMode, "fixed")) {
            rng <- range(rowRanges(se), ignore.strand = TRUE)
            s <- seq(start(rng), end(rng) - windowSize + 1, by = windowStep)
            windowgr <- GRanges(seqnames = seqnames(rng),
                                ranges = IRanges(start = s, width = windowSize))
        }
        ov <- findOverlaps(query = rowRanges(se), subject = windowgr,
                           ignore.strand = TRUE)

        # aggregate counts in windows
        ## TODO: this block needs to be modular and just produce an assayList
        ##       (possibly with the correct, empty elements, replacing the "else"
        ##        below, conditionally putting rowRanges = GPos() if there is no data)
        ## needed modules (window quantification functions):
        ## - sum Nmod and Nvalid per window (current block)
        ## - extract FracMod vector per window and perform stft (phasing score per window)
        ##   (this will need a new scoreFunction, which extract an assay from an SE and
        ##    adds columns to mcols(rowRanges(SE)). this scoreFunction may be combined
        ##    with all "per-sample" window quantification functions)
        mNmod <- rowsum(x = assay(se, "Nmod")[queryHits(ov), ],
                        group = subjectHits(ov), reorder = TRUE)
        mNvalid <- rowsum(x = assay(se, "Nvalid")[queryHits(ov), ],
                          group = subjectHits(ov), reorder = TRUE)
        rnms <- as.numeric(rownames(mNmod))
        stopifnot(exprs = {
            identical(rownames(mNmod), rownames(mNvalid))
            all(diff(rnms) > 0)
        })

        # construct SummarizedExperiment
        seNew <- SummarizedExperiment(assays = list(Nmod = mNmod,
                                                    Nvalid = mNvalid),
                                      rowRanges = windowgr[rnms],
                                      colData = colData(se),
                                      metadata = metadata(se))
    } else {
        seNew <- SummarizedExperiment(
            assays = list(Nmod = matrix(nrow = 0, ncol = ncol(se)),
                          Nvalid = matrix(nrow = 0, ncol = ncol(se))),
            rowRanges = GPos(),
            colData = colData(se),
            metadata = metadata(se))
    }
    return(seNew)
}


#' Analyze counts for sequential windows in a single region
#'
#' Given a \code{SummarizedExperiment} with modified and total
#' base counts in assays \code{"Nmod"} and \code{"Nvalid"},
#' perform a pairwise statistical test for differential modification.
#'
#' @param se \code{SummarizedExperiment}, for example returned by
#'     \code{quantifyWindowsInRegion}. It is expected to at least
#'     contain assays for modified and total counts (given by
#'     \code{assayNameMod} and \code{assayNameValid}) and
#'     a \code{colData} column that defines the groups (given
#'     by \code{groupCol}).
#' @param assayNameMod,assayNameValid Character scalars that give
#'     the assay names in \code{se} containing the modified and
#'     total counts, respectively.
#' @param groupCol Character scalar giving the column in \code{colData(se)}
#'     that defines the groups of samples to be compared.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Panagiotis Papasaikas, Sebastien Smallwood, Charlotte Soneson, Michael Stadler
#'
#' @returns The \code{\link[GenomicRanges]{GRanges}} object constructed from
#'     the \code{\link[edgeR]{topTags}} output obtained for the statistical
#'     analysis, with an additional column named "dirNegLog10PValue", calculated
#'     as the sign of the logFC multiplied with the -log10(PValue).
#'
#' @examples
#' modbamfiles <- system.file("extdata",
#'                            c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
#'                              "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' se <- quantifyWindowsInRegion(bamfiles = modbamfiles,
#'                               region = "chr1:6940000-6955000", modbase = "a",
#'                               BPPARAM = BiocParallel::SerialParam())
#' se$group <- c("group1", "group1", "group2", "group2")
#' gr <- getDifferentiallyModifiedWindows(se, groupCol = "group")
#' class(gr)
#' head(gr)
#'
#' @importFrom SummarizedExperiment assayNames colData assay ncol
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols<- DataFrame
#' @importFrom stats model.matrix
#' @importFrom methods as
#' @importFrom cli cli_abort
#'
#' @export
getDifferentiallyModifiedWindows <- function(se,
                                             assayNameMod = "Nmod",
                                             assayNameValid = "Nvalid",
                                             groupCol = "group",
                                             verbose = FALSE) {
    # check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assayNameMod, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = assayNameValid, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = groupCol, type = "character",
                  validValues = colnames(colData(se)))
    if (length(unique(colData(se)[[groupCol]])) != 2) {
        cli_abort(paste0(
            "The group column in {.code colData(se)} ({groupCol}) ",
            "needs to have exactly two unique values."))
    }
    .assertScalar(x = verbose, type = "logical")
    .assertPackagesAvailable(pkgs = "edgeR")

    if (nrow(se) > 0) {
        # calculate library size and unmodified counts
        .message("calculating library sizes and normalization factors")
        libsizes <- colSums(assay(se, assayNameValid))
        nfacts <- edgeR::normLibSizes(assay(se, assayNameValid))
        cnt <- cbind(assay(se, assayNameMod),
                     assay(se, assayNameValid) - assay(se, assayNameMod))
        cd2 <- cbind(rbind(colData(se)[, c("sample", groupCol)],
                           colData(se)[, c("sample", groupCol)]),
                     data.frame(type = rep(c("mod", "unmod"), each = ncol(se))))
        cd2$group <- factor(cd2$group)

        # create design matrix
        .message("creating design matrix")
        dsgn <- model.matrix(~ 0 + sample, data = cd2)
        for (grp in levels(cd2$group)) {
            dsgn <- cbind(dsgn, cd2$group == grp & cd2$type == "mod")
        }
        colnames(dsgn)[ncol(dsgn) - c(1, 0)] <- levels(cd2$group)
        rownames(cnt) <- NULL
        colnames(cnt) <- rownames(dsgn) <- paste0(rep(colnames(se), 2),
                                                  rep(c(".mod", ".unmod"),
                                                      each = ncol(se)))

        # test for differential modification
        .message("testing for differential modifications")
        dgeL <- edgeR::DGEList(counts = cnt, lib.size = rep(libsizes, 2),
                               norm.factors = rep(nfacts, 2),
                               genes = as.data.frame(unname(rowRanges(se))))
        dgeL <- edgeR::estimateDisp(y = dgeL, design = dsgn)
        fit <- edgeR::glmFit(y = dgeL, design = dsgn)
        tst <- edgeR::glmLRT(
            glmfit = fit,
            contrast = (colnames(dsgn) == levels(cd2$group)[2]) -
                (colnames(dsgn) == levels(cd2$group)[1]))

        # coerce topTags to GRanges
        tt <- edgeR::topTags(object = tst, n = Inf, sort.by = "none")
        tt$table$dirNegLog10PValue <- sign(tt$table$logFC) * -log10(tt$table$PValue)
        gr <- as(tt$table, "GRanges")
    } else {
        gr <- GRanges()
        mcols(gr) <- DataFrame(logFC = numeric(0),
                               logCPM = numeric(0),
                               LR = numeric(0),
                               PValue = numeric(0),
                               FDR = numeric(0),
                               dirNegLog10PValue = numeric(0))
    }
    return(gr)
}

#' Identify regions of interest genome-wide.
#'
#' Given scores or statistical estimates for windows, identify
#' regions of interest by fusing consistent neighboring windows
#' along the genome.
#'
#' @param x \code{GRanges} object with window-based scores or estimates.
#'     Ranges correspond to windows, and columns in \code{mcols(x)} to scores.
#' @param scoreCol Character scalar giving the column name in \code{mcols(x)}
#'     to use for the analysis.
#' @param thresh A numeric scalar giving the minimal absolute window score
#'     (after smoothing, see \code{minperiod} argument) defining a region of
#'     interest. Higher values make the region detection more stringent.
#' @param minperiod Numeric scalar that defines the low-pass
#'     filter parameter used to smooth the scores for segmentation.
#'     \code{minperiod} gives the minimal period (in number of windows) for the
#'     critical frequency in the score signal that should be retained. The
#'     default value is suitable to retain signals that occur in three
#'     neighboring windows. The filtering is performed using
#'     \code{\link[signal]{filtfilt}} and thus requires the \code{signal}
#'     package to be installed.
#' @param maxGap Numeric scalar giving the maximal gap between neighboring
#'     windows, in base pairs, from the end of the first to the start of the
#'     next, to be fused into a single region of interest.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Sebastien Smallwood, Charlotte Soneson, Michael Stadler
#'
#' @returns A \code{\link[GenomicRanges]{GRanges}} object with identified
#'     regions of interest.
#'
#' @examples
#' modbamfiles <- system.file("extdata",
#'                            c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
#'                              "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' se <- quantifyWindowsInRegion(bamfiles = modbamfiles,
#'                               region = "chr1:6940000-6955000", modbase = "a",
#'                               BPPARAM = BiocParallel::SerialParam())
#' se$group <- c("group1", "group1", "group2", "group2")
#' gr <- getDifferentiallyModifiedWindows(se, groupCol = "group")
#' gr
#'
#' grFused <- fuseWindows(x = gr, scoreCol = "logFC", thresh = 5.0)
#' grFused
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom BiocGenerics sort
#' @importFrom S4Vectors mcols mcols<- subjectHits queryHits DataFrame
#' @importFrom dplyr mutate filter group_by ungroup group_split
#' @importFrom IRanges reduce findOverlaps
#' @importFrom cli cli_abort
#' @importFrom rlang .data
#'
#' @export
fuseWindows <- function(x,
                        scoreCol = "dirNegLog10PValue",
                        thresh = 3,
                        minperiod = 3,
                        maxGap = 50,
                        verbose = FALSE) {
    # check argument values
    .assertVector(x = x, type = "GRanges")
    .assertScalar(x = scoreCol, type = "character", validValues = colnames(mcols(x)))
    .assertScalar(x = thresh, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = minperiod, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = maxGap, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")
    .assertPackagesAvailable(pkgs = "signal")

    # smooth scores
    .message("smoothing windows")
    xdf <- as.data.frame(x) |>
        mutate(chunkId = cumsum(c(1, (start[-1] - end[-length(x)] > maxGap |
                                          seqnames[-1] != seqnames[-length(x)])))) |>
        group_by(.data$chunkId) |>
        mutate(sscore = .filterScores(.data[[scoreCol]],
                                      minperiod = minperiod,
                                      maxperiod = Inf,
                                      type = "low")) |>
        ungroup()

    # threshold
    .message("thresholding smoothed scores")
    xdfSel <- xdf |>
        dplyr::filter(abs(.data$sscore) >= thresh) |>
        mutate(direction = factor(ifelse(sign(.data$sscore) == -1, "down", "up"),
                                  levels = c("down", "up")))

    # summarise
    .message("summarise {nrow(xdfSel)} window{?s} into regions of interest")
    grL <- xdfSel |>
        group_by(.data$direction) |>
        group_split() |>
        lapply(function(x) {
            gr1 <- as(x, "GRanges") |>
                reduce(min.gapwidth = maxGap)
            ov1 <- findOverlaps(query = as(x, "GRanges"),
                                subject = gr1, type = "within")
            mcols(gr1)[[paste0(scoreCol, "Thresh")]] <- as.vector(
                tapply(X = x$sscore[queryHits(ov1)],
                       INDEX = subjectHits(ov1),
                       FUN = mean))
            mcols(gr1)[["numWindowsThresh"]] <- tabulate(subjectHits(ov1))
            gr1$direction <- x$direction[1]
            gr1
        })
    if (sum(lengths(grL)) > 0) {
        gr <- sort(sort(do.call(c, grL)), ignore.strand = TRUE)
        ov <- findOverlaps(query = x, subject = gr, type = "within")
        mcols(gr)[[scoreCol]] <- as.vector(
            tapply(X = xdf$sscore[queryHits(ov)],
                   INDEX = subjectHits(ov),
                   FUN = mean))
        mcols(gr)[["numWindows"]] <- tabulate(subjectHits(ov))
    } else {
        gr <- GRanges()
        mcols(gr) <- DataFrame(thresh = numeric(0),
                               numWindowsThresh = integer(0),
                               direction = factor(character(0),
                                                  levels = c("down", "up")),
                               score = numeric(0),
                               numWindows = numeric(0))
        colnames(mcols(gr))[c(1,4)] <- paste0(scoreCol, c("Thresh", ""))
    }

    return(gr)
}

#' Scan one or more chromosomes for high-scoring regions
#'
#' Given one or several modBam files and a vector with chromosome lengths,
#' quantify the modification data in windows, calculate a score for each
#' window and fuse neighboring high-scoring windows. For example, the score
#' could be a directional P value measuring differential modification between
#' two groups of modBam files, and the identified windows would then be
#' differentially modified regions (DMRs).
#'
#' @inheritParams quantifyWindowsInRegion
#' @inheritParams fuseWindows
#'
#' @param chromosomeLengths A named character vector with lengths of
#'     chromosomes to be analyzed, for example the return value of
#'     \code{seqlengths(BSgenomeObject)}.
#' @param scoreFunction A character scalar giving the name of the scoring
#'     function. This function should take as a first argument a
#'     \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, as
#'     generated by \code{quantifyWindowsInRegion}. Additional arguments may
#'     be passed to \code{scoreFunction} using \code{...}. The function must
#'     return a \code{\link[GenomicRanges]{GRanges}} object (typically
#'     corresponding to the \code{rowRanges} of the input) with the score
#'     added to the \code{mcols()}.
#' @param ... Additional parameters sent to \code{scoreFunction}.
#' @param tileSize A numeric scalar giving the approximate size of a
#'     genomic region for which the data is loaded at once. Chromosomes that
#'     are larger than \code{tileSize} will be automatically split into
#'     multiple tiles. Reduce this size if you would like to reduce the memory
#'     usage of the function.
#'
#' @author Charlotte Soneson, Michael Stadler
#'
#' @returns A \code{\link[GenomicRanges]{GRanges}} object with identified
#'     regions of interest.
#'
#' @examples
#' modbamfiles <- system.file("extdata",
#'                            c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
#'                              "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' gr <- scanForHighScoringRegions(bamfiles = modbamfiles,
#'                                 sampleAnnot = data.frame(sample = c("s1","s2","s3","s4"),
#'                                                          group = c("A","A","B","B")),
#'                                 chromosomeLengths = c(chr1 = 6955000),
#'                                 modbase = "a", BPPARAM = BiocParallel::SerialParam())
#' gr
#'
#' @importFrom BiocParallel MulticoreParam
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom cli cli_abort
#'
#' @export
scanForHighScoringRegions <- function(bamfiles,
                                      sampleAnnot,
                                      chromosomeLengths,
                                      scoreFunction = "getDifferentiallyModifiedWindows",
                                      modbase,
                                      modProbThreshold = 0.5,
                                      tileSize = 1e6,
                                      seqinfo = NULL,
                                      sequenceContextWidth = 0,
                                      sequenceReference = NULL,
                                      sequenceContext = NULL,
                                      windowMode = "fixed",
                                      windowSize = 24,
                                      windowStep = round(windowSize/2),
                                      scoreCol = "dirNegLog10PValue",
                                      thresh = 3,
                                      minperiod = 3,
                                      maxGap = 50,
                                      BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                                      verbose = FALSE,
                                      ...) {
    # check parameters
    .assertVector(x = chromosomeLengths, type = "numeric")
    if (is.null(names(chromosomeLengths)) ||
        any(duplicated(names(chromosomeLengths))) ||
        any(names(chromosomeLengths) == "")) {
        cli_abort("{.arg chromosomeLengths} must be a named vector")
    }
    .assertScalar(x = scoreFunction, type = "character")
    if (!exists(scoreFunction)) {
        cli_abort("{.arg scoreFunction} must be the name of an existing function")
    } else {
        scoreFunction <- get(scoreFunction)
    }
    .assertScalar(x = tileSize, type = "numeric", rngExcl = c(0, Inf))

    # loop over chromosomes
    gr <- do.call(c, lapply(names(chromosomeLengths), function(chr) {
        regs <- .tileChromosome(tileSize = tileSize,
                                windowSize = windowSize,
                                windowStep = windowStep,
                                chromName = chr,
                                chromLength = chromosomeLengths[chr])

        # quantify windows for each tile and merge
        tileL <- lapply(seq_along(regs), function(i) {
            quantifyWindowsInRegion(bamfiles = bamfiles,
                                    region = regs[i],
                                    modbase = modbase,
                                    modProbThreshold = modProbThreshold,
                                    sampleAnnot = sampleAnnot,
                                    seqinfo = seqinfo,
                                    sequenceContextWidth = sequenceContextWidth,
                                    sequenceReference = sequenceReference,
                                    sequenceContext = sequenceContext,
                                    windowMode = windowMode,
                                    windowSize = windowSize,
                                    windowStep = windowStep,
                                    BPPARAM = BPPARAM,
                                    verbose = verbose)
        })
        se <- do.call(rbind, tileL[unlist(lapply(tileL, function(x) nrow(x) > 0))])

        # calculate window scores
        grScores <- scoreFunction(se, ...)

        # fuse windows
        grScoresFused <- fuseWindows(x = grScores,
                                     scoreCol = scoreCol,
                                     thresh = thresh,
                                     minperiod = minperiod,
                                     maxGap = maxGap,
                                     verbose = verbose)

        # return fused windows
        return(grScoresFused)
    }))

    return(gr)
}
