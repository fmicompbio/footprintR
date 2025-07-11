# global data.frame of plot types and characteristics
plotRegionPlotTypes <- data.frame(
    name = c("Point", "Smooth", "PointSmooth",
             "Lollipop", "Heatmap", "GenomicRegion",
             "GenomicRegions", "BigWig"),
    type = c("summary", "summary", "summary",
             "reads", "reads", "annotation", "annotation",
             "file")
)

defaultFootprintColors <- c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4",
                            "#FDDAEC", "#FED9A6", "#FFFFCC", "#E5D8BD")

.convertRegionToModBaseSpace <- function(regdf, datadf) {
    if (is.factor(datadf$position)) {
        lvs <- as.numeric(levels(datadf$position))
        regdf$start <- vapply(regdf$start, function(s) {
            if (any(lvs >= s)) {
                as.character(min(lvs[lvs >= s]))
            } else {
                NA_character_
            }
        }, NA_character_)
        regdf$end <- vapply(regdf$end, function(s) {
            if (any(lvs <= s)) {
                as.character(max(lvs[lvs <= s]))
            } else {
                NA_character_
            }
        }, NA_character_)
        regdf <- regdf[rowSums(is.na(regdf)) == 0, ]
        regdf$start <- factor(regdf$start, levels = levels(datadf$position))
        regdf$end <- factor(regdf$end, levels = levels(datadf$position))
    }
    regdf
}

#' Plot single-molecule footprinting data for a single genomic region
#'
#' @description
#' The \code{plotRegion} function visualizes read-level or collapsed
#' single-molecule footprinting data, such as data imported using
#' \code{\link{readModkitExtract}}, \code{\link{readModBam}} or
#' \code{\link{readBedMethyl}}. The \code{plotReadsLollipop},
#' \code{plotReadsHeatmap}, \code{plotSummaryPointSmooth} and
#' \code{plotGenomicRegions} functions are helper functions for creating
#' single plot tracks. These are invoked by \code{plotRegion},
#' and typically do not need to be directly called by the user.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level or collapsed single-molecule footprinting data (positions
#'     in rows and samples in columns).
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region. Only data from \code{se} overlapping this region will be plotted.
#'     Alternatively, the region can be specified as a character scalar (e.g.
#'     "chr1:1200-1300") that can be coerced into a \code{GRanges} object. If
#'     \code{NULL} (the default), all the data on the first sequence in
#'     \code{se} will be visualized.
#' @param tracks A list of named lists, representing the tracks to generate.
#'     Each element of the outer list defines one track, and has to contain
#'     at least list entries named 'trackData' (the name of a suitable assay
#'     in \code{se}, for data tracks, or a
#'     \code{\link[GenomicRanges]{GRangesList}} object for the annotation
#'     tracks) and 'trackType' (the type of plot), plus any additional
#'     arguments to the respective plot function. Currently supported plot
#'     types are
#'     \describe{
#'         \item{\code{"Point"}}{: A point plot displaying values in the assay.}
#'         \item{\code{"Smooth"}}{: A smoothed line plot displaying values in
#'             the assay.}
#'         \item{\code{"PointSmooth"}}{: A point and smoothed line plot
#'             displaying values in the assay.}
#'         \item{\code{"Lollipop"}}{: Lollipop plot (filled circles with the
#'             color representing the values in the assay).}
#'         \item{\code{"Heatmap"}}{: Heatmap plot (tiles with the color
#'             representing the values in the assay).}
#'         \item{\code{"GenomicRegion"} or \code{"GenomicRegions"}}{: Genomic
#'             annotations (e.g., transcripts, peaks, CpG islands).}
#'     }
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the x-axis will be
#'     shown in the space of modified bases and contain only the positions at
#'     which there are modified bases in the data without any gaps between them.
#'     If \code{FALSE}, the x-axis will show the genomic coordinate on which
#'     the modified bases are typically irregularly spaced.
#' @param sequenceContext A character vector with sequence context(s)
#'     to plot. Only positions that match one of the provided sequence
#'     contexts will be included in the plot. Sequence contexts can be provided
#'     using IUPAC redundancy codes. The sequence contexts of modified bases are
#'     obtained from \code{rowData(se)$sequenceContext} and thus requires that
#'     \code{se} contains the appropriate information, for example by setting
#'     the \code{sequenceContextWidth} and \code{sequenceReference} arguments of
#'     \code{\link{readBedMethyl}}, \code{\link{readModBam}} or
#'     \code{\link{readModkitExtract}} when reading data, or by adding it using
#'     \code{\link{addSeqContext}}.
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used. Ignored if \code{modbaseSpace} is
#'     \code{TRUE}.
#' @param labelAccuracy A numeric scalar indicating the precision of the
#'     positions along the genomic axis. Will be passed to
#'     \code{\link[scales]{label_number}}. If \code{NULL} (default), a suitable
#'     value will be derived from \code{region}.
#' @param suppressTickLabels Logical scalar. If \code{TRUE}, suppress x-axis
#'     tick labels for all but the last panel.
#' @param minCoveredFraction A numeric scalar giving the minimal fraction of
#'     \code{region} that a read needs to cover to be plotted.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object with tracks selected by
#'     \code{tracks}.
#'
#' @author Charlotte Soneson, Michael Stadler
#' @name plotRegion
#'
#' @examples
#' # summarized data (5mC)
#' bmfiles <- system.file("extdata",
#'                        c("modkit_pileup_1.bed.gz", "modkit_pileup_2.bed.gz"),
#'                        package = "footprintR")
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' seA <- readBedMethyl(bmfiles, modbase = "m",
#'                      sequenceContextWidth = 3, sequenceReference = reffile,
#'                      BPPARAM = BiocParallel::SerialParam())
#'
#' plotRegion(seA, region = "chr1:6940000-6955000", sequenceContext = "GCH")
#' plotRegion(seA, region = "chr1:6940000-6955000", sequenceContext = "HCG")
#'
#' plotRegion(seA, region = "chr1:6940000-6955000",
#'            tracks = list(list(trackData = "Nvalid", trackType = "Smooth")))
#'
#' # read-level data (6mA)
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "footprintR")
#' seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
#'                          BPPARAM = BiocParallel::SerialParam())
#'
#' # Lollipop plot
#' plotRegion(seB, region = "chr1:6935800-6935900", minCoveredFraction = 0.95,
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
#'                               orderReads = "regionAvg")))
#' # Heatmap plots (observed only or interpolated)
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Heatmap")))
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
#'                               interpolate = TRUE)))
#'
#' # multiple plots, in 'modbase' space
#' plotRegion(seB, region = "chr1:6935400-6935450",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
#'                               size = 4),
#'                          list(trackData = "mod_prob", trackType = "Heatmap")),
#'            modbaseSpace = TRUE)
#'
#' # combine read-level and summary tracks,
#' # set relative heights of tracks, don't facet by sample,
#' # change titles of legends
#' seB <- flattenReadLevelAssay(seB, assayName = "mod_prob")
#' plotRegion(seB, region = "chr1:6935400-6935450",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
#'                               size = 4, legendTitle = "6mA",
#'                               facetBy = NULL),
#'                          list(trackData = "FracMod", trackType = "Smooth")),
#'            modbaseSpace = TRUE) +
#'     patchwork::plot_layout(heights = c(3, 2))
#'
#' @seealso \code{\link{readModBam}}, \code{\link{readModkitExtract}} and
#'     \code{\link{readBedMethyl}} for reading read-level and summarized
#'     footprinting data.
#'
#' @importFrom SummarizedExperiment assay assayNames rowRanges assays
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocGenerics nrow ncol
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom cli cli_abort cli_warn
#' @importFrom methods as
#'
#' @export
plotRegion <- function(
        se,
        region = NULL,
        tracks = list(list(trackData = "FracMod", trackType = "Point")),
        modbaseSpace = FALSE,
        sequenceContext = NULL,
        referenceCoordinate = NULL,
        labelAccuracy = NULL,
        suppressTickLabels = FALSE,
        minCoveredFraction = 0) {

    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    } else if (is.null(region)) {
        # get the range of covered positions on the first seqname
        region <- range(
            rowRanges(se)[seqnames(rowRanges(se)) == seqlevels(se)[1]],
            ignore.strand = TRUE)
    }
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertVector(x = tracks, type = "list", rngLen = c(1, Inf))
    assaysInUse <- c()
    for (i in seq_along(tracks)) {
        if (!is.list(tracks[[i]]) || length(tracks[[i]]) < 2) {
            cli_abort("{.code tracks[[{i}]]} has to be a {.cls list} of length >=2.")
        }
        if (is.null(names(tracks[[i]])) || any(names(tracks[[i]]) == "") ||
            any(!c("trackData", "trackType") %in% names(tracks[[i]]))) {
            cli_abort(paste("{.code tracks[[{i}]]} must be a named {.cls list}, and ",
                            "contain at least entries named 'trackData' and 'trackType'"))
        }
        if (!tracks[[i]]$trackType %in% plotRegionPlotTypes$name) {
            cli_abort("{.code tracks[[{i}]]$trackType} must be one of: {plotRegionPlotTypes$name}")
        }
        if (is.character(tracks[[i]]$trackData) &&
            length(tracks[[i]]$trackData) == 1 &&
            tracks[[i]]$trackData == "FracMod" &&
            !"FracMod" %in% assayNames(se)) {
            if (all(c("Nmod", "Nvalid") %in% assayNames(se))) {
                assay(se, "FracMod") <- assay(se, "Nmod") / assay(se, "Nvalid")
            } else {
                cli_abort(paste("Cannot plot 'FracMod' - need either an assay",
                                "called 'FracMod' or both 'Nmod' and 'Nvalid'",
                                "assays"))
            }
        }
        type_i <- plotRegionPlotTypes$type[match(tracks[[i]]$trackType,
                                                 plotRegionPlotTypes$name)]
        if (type_i %in% c("reads", "summary") &&
            !(is.character(tracks[[i]]$trackData) &&
              length(tracks[[i]]$trackData) == 1 &&
              tracks[[i]]$trackData %in% assayNames(se))) {
            cli_abort(paste("{.code tracks[[{i}]]$trackData} must be a ",
                            "{.cls character} scalar corresponding to a name ",
                            "of an assay in {.arg se}"))
        }
        if (type_i == "reads" &&
            !tracks[[i]]$trackData %in% .getReadLevelAssayNames(se)) {
            cli_abort(paste("{.code tracks[[{i}]]$trackData} must be the name ",
                            "of a read-level assay in {.arg se}"))
        }
        if (type_i == "summary" &&
            tracks[[i]]$trackData %in% .getReadLevelAssayNames(se)) {
            cli_abort(paste("{.code tracks[[{i}]]$trackData} must be the name ",
                            "of a summary assay in {.arg se}"))
        }
        if (type_i %in% c("reads", "summary")) {
            # add assay to list of assays that are required for the plots
            assaysInUse <- union(assaysInUse, tracks[[i]]$trackData)
        }
        if (type_i == "annotation" &&
            !(is(tracks[[i]]$trackData, "GRangesList") &&
              !is.null(names(tracks[[i]]$trackData)))) {
            cli_abort(paste("{.code tracks[[{i}]]$trackData} must be a named",
                            "{.cls GRangesList} object"))
        }
        if (type_i == "annotation" &&
            length(unlist(lapply(tracks[[i]]$trackData,
                                 function(y) unique(as.character(strand(y)))))) !=
            length(tracks[[i]]$trackData)) {
            cli_abort(paste("There are entries in {.code tracks[[{i}]]$trackData}",
                            "with mixed strand annotations"))
        }
        if (type_i == "file" &&
            !(is.character(tracks[[i]]$trackData) &&
              !is.null(names(tracks[[i]]$trackData)))) {
            cli_abort(paste("{.code tracks[[{i}]]$trackData} must be a named",
                            "character vector"))
        }

        if (modbaseSpace &&
            "interpolate" %in% names(tracks[[i]]) &&
            is.logical(tracks[[i]]$interpolate) &&
            length(tracks[[i]]$interpolate) == 1 &&
            tracks[[i]]$interpolate) {
            cli_warn(paste("Plotting in `modbaseSpace` is not allowed if",
                           "interpolate = TRUE (seen in tracks[[{i}]]).",
                           "Setting modbaseSpace=FALSE"))
            modbaseSpace <- FALSE
        }
        if (modbaseSpace &&
            tracks[[i]]$trackType %in% c("GenomicRegion", "GenomicRegions")) {
            cli_warn(paste("Plotting in `modbaseSpace` is not allowed if",
                           "GenomicRegion tracks are included.",
                           "Setting modbaseSpace=FALSE"))
            modbaseSpace <- FALSE
        }
        if (modbaseSpace && "footprintColumns" %in% names(tracks[[i]]) &&
            is.character(tracks[[i]]$footprintColumns) &&
            length(tracks[[i]]$footprintColumns) > 0) {
            cli_warn(paste("Plotting in `modbaseSpace` is not allowed if",
                           "footprintColumns are provided (seen in
                           tracks[[{i}]]. Setting modbaseSpace=FALSE"))
            modbaseSpace <- FALSE
        }
        if (modbaseSpace &&
            tracks[[i]]$trackType == "BigWig") {
            cli_warn(paste("Plotting in `modbaseSpace` is not allowed if",
                           "BigWig tracks are included.",
                           "Setting modbaseSpace=FALSE"))
            modbaseSpace <- FALSE
        }
    }
    .assertVector(x = sequenceContext, type = "character", allowNULL = TRUE)
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minCoveredFraction, type = "numeric", rngIncl = c(0, 1))

    if (modbaseSpace) {
        # relative coordinates are not meaningful in modbase space (as there
        # are no actual positions indicated anyway)
        referenceCoordinate <- NULL
    }

    # remove assays that are not used in any plot, for faster subsetting of se
    suppressWarnings(
        # currently, assigning to assays triggers a deprecation warning
        # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
        # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
        assays(se) <- assays(se)[assaysInUse]
    )

    # subset positions in SE (respecting the strand of the region)
    se <- subsetByOverlaps(x = se, ranges = region)
    se <- .keepPositionsBySequenceContext(
        se = se, sequenceContext = sequenceContext)

    if (nrow(se) == 0) {
        cli_abort("No positions retained for plotting!")
    }

    # subset reads in SE
    if (minCoveredFraction > 0) {
        for (nm in intersect(assaysInUse, .getReadLevelAssayNames(se))) {
            se <- filterReads(se = se, assayName = nm, readInfoCol = NULL,
                              qcCol = NULL, minCoveredFraction = minCoveredFraction,
                              region = region, prune = TRUE)
        }
    }

    if (ncol(se) == 0) {
        cli_abort("No reads retained for plotting!")
    }

    ## create plots
    pL <- vector("list", length = length(tracks))
    for (i in seq_along(tracks)) {
        tr <- tracks[[i]]
        trt <- plotRegionPlotTypes$type[
            match(tr$trackType, plotRegionPlotTypes$name)]
        if (trt %in% c("summary", "reads")) {
            args <- c(
                list(se = quote(se), region = region, assayName = tr$trackData,
                     modbaseSpace = modbaseSpace,
                     referenceCoordinate = referenceCoordinate,
                     labelAccuracy = labelAccuracy),
                tr[!names(tr) %in% c("trackData", "trackType", "se", "region",
                                     "assayName", "modbaseSpace", "doSmooth",
                                     "doPoint", "referenceCoordinate",
                                     "labelAccuracy")]
            )
        } else if (trt == "annotation") {
            args <- c(
                list(grl = quote(subsetByOverlaps(tr$trackData, region)),
                     region = region, labelAccuracy = labelAccuracy,
                     referenceCoordinate = referenceCoordinate),
                     tr[!names(tr) %in% c("trackData", "trackType", "grl",
                                          "region", "referenceCoordinate",
                                          "labelAccuracy")]
            )
        } else if (trt == "file") {
            args <- c(
                list(bwFiles = tr$trackData,
                     region = region, labelAccuracy = labelAccuracy,
                     referenceCoordinate = referenceCoordinate),
                tr[!names(tr) %in% c("trackData", "trackType", "region",
                                     "referenceCoordinate", "labelAccuracy",
                                     "bwFiles")]
            )
        }
        pL[[i]] <- switch(
            tr$trackType,
            Point = do.call(plotSummaryPointSmooth,
                            c(args, list(doSmooth = FALSE))),
            Smooth = do.call(plotSummaryPointSmooth,
                             c(args, list(doPoint = FALSE))),
            PointSmooth = do.call(plotSummaryPointSmooth, args),
            Lollipop = do.call(plotReadsLollipop, args),
            Heatmap = do.call(plotReadsHeatmap, args),
            GenomicRegion = do.call(plotGenomicRegions, args),
            GenomicRegions = do.call(plotGenomicRegions, args),
            BigWig = do.call(plotBigWig, args)
        )
    }

    ## assemble composite plot
    if (length(pL) > 1L) { # suppress x-axis labels for all but last plot
        for (i in seq.int(length(pL) - 1L)) {
            pL[[i]] <- pL[[i]] + labs(x = NULL)
            if (suppressTickLabels) {
                pL[[i]] <- pL[[i]] + theme(axis.text.x = element_blank())
            }
        }
    }
    p <- wrap_plots(pL, ncol = 1)
    if (!is.null(sequenceContext)) {
        p <- p + labs(caption = paste0("Sequence contexts: ",
                                       paste(sequenceContext, collapse = ", ")))
    }

    # return
    return(p)
}


## plot* functions for plotRegion() -------------------------------------------

#' @param bwFiles A named character vector with paths to one or more bigWig
#'     files to plot.
#' @param yAxisLabel A character scalar providing the label to use for the
#'     y-axis.
#'
#' @importFrom cli cli_abort
#' @importFrom dplyr bind_rows mutate group_by group_modify ungroup select
#' @importFrom BiocGenerics setdiff
#'
#' @export
#' @rdname plotRegion
#'
plotBigWig <- function(bwFiles,
                       region,
                       trackTitle = NULL,
                       legendTitle = NULL,
                       yAxisLabel = "Score",
                       showLegend = TRUE,
                       highlightRegions = NULL,
                       colors = NULL,
                       referenceCoordinate = NULL,
                       labelAccuracy = NULL,
                       yAxisRange = NULL) {

    .assertVector(x = bwFiles, type = "character")
    .assertVector(x = names(bwFiles), type = "character")
    if (any(i <- !file.exists(bwFiles))) {
        cli_abort("Not all bigWig files exist: {bwFiles[i]}")
    }
    if (any(i <- duplicated(names(bwFiles)))) {
        cli_abort("Duplicated file names: {unique(names(bwFiles)[i])}")
    }
    .assertPackagesAvailable("BiocIO")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = yAxisLabel, type = "character", allowNULL = FALSE)
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- GenomicRanges::pintersect(highlightRegions, region,
                                                      ignore.strand = TRUE,
                                                      drop.nohit.ranges = TRUE)
        highlightRegions <- highlightRegions[width(highlightRegions) > 0]
        # highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
        #                                             ignore.strand = TRUE)
    }
    .assertVector(x = colors, type = "character", allowNULL = TRUE)
    if (!is.null(colors)) {
        .assertVector(x = names(colors), type = "character")
        if (!all(names(bwFiles) %in% names(colors))) {
            cli_abort("Missing color specification for some values")
        }
    }
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)
    .assertVector(x = yAxisRange, type = "numeric", len = 2,
                  allowNULL = TRUE)

    # Import data and convert to data frames
    df <- do.call(
        bind_rows,
        lapply(structure(names(bwFiles), names = names(bwFiles)),
               function(nm) {
                   x <- BiocIO::import(bwFiles[nm], which = region)
                   y <- BiocGenerics::setdiff(region, x)
                   mcols(y)$score <- rep(0, length(y))
                   x <- sort(c(x, y))
                   as.data.frame(x) |>
                       mutate(idx = seq_along(.data$seqnames)) |>
                       group_by(.data$idx) |>
                       group_modify(~ data.frame(
                           sample = nm,
                           chr = .x$seqnames,
                           position = seq(.x$start, .x$end),
                           strand = .x$strand,
                           value = .x$score)) |>
                       ungroup() |>
                       select(c("position", "sample", "value"))
               }))
    if (!is.null(referenceCoordinate)) {
        # shift all ranges
        df$position <- df$position - referenceCoordinate
        if (!is.null(highlightRegions)) {
            highlightRegions <- shift(highlightRegions, -referenceCoordinate)
        }
        region <- shift(region, -referenceCoordinate)
    }

    # create base plot
    if (is.null(yAxisRange)) {
        yAxisRange <- c(0, NA)
    } else {
        yAxisRange[1] <- 0
    }
    p <- .createBaseplotSummary(df = df,
                                region = region,
                                trackTitle = trackTitle,
                                legendTitle = legendTitle,
                                showLegend = showLegend,
                                highlightRegions = highlightRegions,
                                groupBy = "sample",
                                colorBy = "sample",
                                colors = colors,
                                referenceCoordinate = referenceCoordinate,
                                labelAccuracy = labelAccuracy,
                                yAxisLabel = yAxisLabel,
                                yAxisRange = yAxisRange)

    # add geom
    p <- p + geom_area()

    # facet (if there are more than one sample)
    if (length(unique(df$sample)) > 1) {
        p <- p + facet_wrap(~ sample, ncol = 1) +
            theme(strip.background.x = element_blank(),
                  strip.text.x = element_text(
                      hjust = 0, margin = margin(t = 0, r = 0, b = 2, l = 0)))
    }

    # Remove top/right axes/boundaries
    p + theme(panel.border = element_blank(),
              axis.line = element_line(lineend = "square"))
}

#' @param assayName A character or numerical scalar selecting the assay to plot.
#'     This should be an existing read-level or summary assay, as appropriate
#'     for the plot type. A special case is the track name \code{"FracMod"}:
#'     If \code{se} does not contain a summary assay of that name, but
#'     \code{"Nmod"} and \code{"Nvalid"} assays are available, \code{"FracMod"}
#'     will be calculated from \code{assay(se, "Nmod") / assay(se, "Nvalid")}.
#' @param size A numeric scalar giving the size of the points (\code{size}
#'     argument of \code{\link[ggplot2]{geom_point}}).
#' @param stroke A numeric scalar giving the stroke (line width) of the point
#'     outlines (\code{stroke} argument of \code{\link[ggplot2]{geom_point}}).
#' @param drawRead A logical scalar. If \code{TRUE}, draw a horizontal line
#'     segment for each read from its start to its end.
#' @param orderReads A character scalar, or \code{NULL}. If \code{"cluster"},
#'     the position of reads on the y-axis will be reordered using
#'     \code{hclust(as.dist(sqrt(2 - 2 * cor(X, method = "pearson",
#'     use = "pairwise.complete"))))$order}, where \code{X} is
#'     \code{assay(x, assayName)} with zero values set to \code{NA} and averaged
#'     over windows of 25 nucleotides. If set to \code{"squish"}, the display
#'     will be compacted by placing multiple reads in the same row when
#'     possible. If set to \code{"regionAvg"}, the reads are sorted by increasing
#'     average modification probability in the window given by \code{orderRegion}.
#'     If \code{NULL}, no reordering is done.
#' @param orderRegion Either \code{NULL} or a length-one \code{GRanges} object.
#'     If \code{orderReads = "regionAvg"}, the \code{GRanges} object defines the
#'     window in which average modification probability is calculated to order
#'     the reads in read-level plots. A \code{NULL} value indicates that the
#'     entire plotted region should be used as the window.
#' @param windowWidth A numeric scalar giving the window width for which
#'     read-level data will be averaged before clustering. This should help
#'     to reduce the noise and allows to compare reads without any common
#'     modification calls, such as plus- and minus-strand reads with 6mA calls.
#' @param clustDist A character scalar defining the distance measure to use
#'     for clustering. Should be one of \code{"pearson"} or
#'     \code{"euclidean"}.
#' @param trackTitle A character scalar or \code{NULL}, giving the title of
#'     the track.
#' @param legendTitle A character scalar or \code{NULL}. If not \code{NULL},
#'     this will be the title of the track fill/color legend. If \code{NULL},
#'     the name of the assay (\code{assayName}, for heatmaps and lollipop plots)
#'     \code{"Sample"} (for summary plots), or \code{"strand"} (for genomic
#'     region plots) will be used.
#' @param showLegend A logical scalar, indicating whether or not to display
#'     the legend for the track.
#' @param highlightRegions A \code{\link[GenomicRanges]{GRanges}} object
#'     containing regions to highlight with a grey shading.
#' @param footprintColumns A character vector with names of columns from
#'     \code{colData(se)} containing footprints to display. Typically these
#'     columns are generated using \code{addFootprints}. If \code{NULL}, no
#'     footprints are displayed.
#' @param arglistFootprints A named list with arguments to be passed to
#'     \code{\link[ggplot2]{geom_tile}}, in addition to the \code{data},
#'     \code{mapping}, and \code{inherit.aes} arguments, which are set
#'     automatically. \code{arglistFootprints} can be either a single
#'     list of such arguments (in which case the same arguments will be used
#'     for all plotted footprint columns), or a named list with one entry per
#'     value in \code{footprintColumns}.
#' @param facetBy A character scalar indicating the sample annotation column
#'     to facet the plot by (if \code{NULL}, no faceting is done). By default,
#'     the plot will be facetted by 'sample', corresponding to the columns of
#'     \code{se}.
#' @param adjustFacetHeight A logical scalar. If \code{TRUE}, adjust the
#'     height of the facets by the number of reads in each of them. If
#'     \code{FALSE}, all facets have the same height.
#' @param fillColors A character scalar defining the continuous color palette to
#'     represent modification probabilities. If \code{fillColors}
#'     has length one, it is assumed to be the a supported value to pass to the
#'     \code{option} argument of \code{\link[ggplot2]{scale_fill_viridis_c}},
#'     optionally prefixed with a minus sign to in addition set
#'     \code{direction = -1}). If \code{fillColors} has more than one element,
#'     it is assumed to be a vector of colors to pass to the \code{colors}
#'     argument of \code{\link[ggplot2]{scale_colour_gradientn}}.
#'
#' @export
#' @rdname plotRegion
#'
#' @examples
#' library(GenomicRanges)
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "footprintR")
#' seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
#'                          BPPARAM = BiocParallel::SerialParam())
#' plotReadsLollipop(seB, region = as("chr1:6935400-6935450", "GRanges"),
#'                   assayName = "mod_prob",
#'                   highlightRegion = GRanges("chr1", IRanges(6935420, 6935430)))
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom dplyr filter
#' @importFrom cli cli_warn
#'
plotReadsLollipop <- function(se,
                              region,
                              assayName,
                              size = 3.0,
                              stroke = 0.5,
                              drawRead = TRUE,
                              orderReads = "cluster",
                              orderRegion = NULL,
                              clustDist = "euclidean",
                              windowWidth = 25,
                              modbaseSpace = FALSE,
                              trackTitle = NULL,
                              legendTitle = NULL,
                              yAxisLabel = "Reads",
                              showLegend = TRUE,
                              highlightRegions = NULL,
                              footprintColumns = NULL,
                              arglistFootprints = list(),
                              facetBy = "sample",
                              adjustFacetHeight = TRUE,
                              referenceCoordinate = NULL,
                              labelAccuracy = NULL,
                              fillColors = "-cividis") {

    argL <- .checkArgsReadLevelPlots(
        se = se, region = region, assayName = assayName, drawRead = drawRead,
        orderReads = orderReads, orderRegion = orderRegion,
        clustDist = clustDist, windowWidth = windowWidth,
        modbaseSpace = modbaseSpace, trackTitle = trackTitle,
        legendTitle = legendTitle, yAxisLabel = yAxisLabel,
        showLegend = showLegend, highlightRegions = highlightRegions,
        footprintColumns = footprintColumns, arglistFootprints = arglistFootprints,
        facetBy = facetBy, adjustFacetHeight = adjustFacetHeight,
        referenceCoordinate = referenceCoordinate,
        labelAccuracy = labelAccuracy, size = size, stroke = stroke)

    # prepare plot data
    df <- .preparePlotdataReads(x = se, assayName = assayName,
                                modbaseSpace = modbaseSpace,
                                referenceCoordinate = argL$referenceCoordinate,
                                extraColAnnots = setdiff(facetBy, "sample"),
                                orderReads = orderReads,
                                orderRegion = argL$orderRegion,
                                clustDist = clustDist,
                                windowWidth = windowWidth,
                                facetBy = facetBy)

    # create base plot
    p <- .createBaseplotReads(df = df, region = argL$region,
                              trackTitle = trackTitle,
                              legendTitle = ifelse(!is.null(legendTitle),
                                                   legendTitle, assayName),
                              yAxisLabel = yAxisLabel,
                              showLegend = showLegend,
                              highlightRegions = argL$highlightRegions,
                              facetBy = facetBy,
                              adjustFacetHeight = adjustFacetHeight,
                              referenceCoordinate = argL$referenceCoordinate,
                              labelAccuracy = labelAccuracy,
                              fillColors = fillColors)

    # add segments (round 1) - need to keep this before the footprints to
    # make sure that the read ordering is respected
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(
            df, groupVars = union(facetBy, "sample"))
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["plotRow"]],
                                  xend = .data[["end"]]
                              ), color = "transparent")
    }

    # add footprints
    for (fpc in footprintColumns) {
        fp <- .prepareFootprintsForPlot(fp = argL$footprints[[fpc]],
                                        plotdf = df, se = se, facetBy = facetBy)
        if (nrow(fp) > 0) {
            aL <- c(
                list(data = fp |> dplyr::filter(!is.na(.data$plotRow)),
                     mapping = aes(x = (as.numeric(.data$start) +
                                            as.numeric(.data$end)) / 2,
                                   y = .data$plotRow,
                                   width = abs(as.numeric(.data$end) -
                                                   as.numeric(.data$start))
                     ),
                     inherit.aes = FALSE),
                argL$arglistFootprints[[fpc]]
            )
            if (is.null(aL$height)) {
                aL$height <- 1
            }
            if (is.null(aL$fill)) {
                aL$fill <- argL$footprintColors[fpc]
            }
            repArgs <- duplicated(names(aL))
            if (any(repArgs)) {
                cli_warn("Ignoring pre-defined arguments: {.arg {names(aL)[repArgs]}}.")
            }
            aL <- aL[!repArgs]
            p <- p + do.call(geom_tile, aL)
        }
    }

    # add segments (round 2) - same data as above, just make sure that it
    # ends up on top of the footprints
    if (drawRead) {
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["plotRow"]],
                                  xend = .data[["end"]]
                              ), color = "gray80")
    }

    # add lollipops
    p <- p + geom_point(shape = 21, size = size,
                        stroke = stroke, color = "black")

    # return plot
    return(p)
}

#' @param linewidthTiles A numeric scalar, the line width of the border drawn
#'     around each measured base.
#' @param interpolate A logical scalar. If \code{TRUE}, the gaps between
#'     observations are filled in by linear interpolation.
#'
#' @export
#' @rdname plotRegion
#'
#' @examples
#' library(GenomicRanges)
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "footprintR")
#' seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
#'                          BPPARAM = BiocParallel::SerialParam())
#' plotReadsHeatmap(seB, region = as("chr1:6935400-6935450", "GRanges"),
#'                  assayName = "mod_prob",
#'                  highlightRegion = GRanges("chr1", IRanges(6935420, 6935430)))
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom dplyr filter
#' @importFrom cli cli_warn
#'
plotReadsHeatmap <- function(se,
                             region,
                             assayName,
                             drawRead = TRUE,
                             linewidthTiles = 0,
                             orderReads = "cluster",
                             orderRegion = NULL,
                             clustDist = "euclidean",
                             windowWidth = 25,
                             modbaseSpace = FALSE,
                             interpolate = FALSE,
                             trackTitle = NULL,
                             legendTitle = NULL,
                             yAxisLabel = "Reads",
                             showLegend = TRUE,
                             highlightRegions = NULL,
                             footprintColumns = NULL,
                             arglistFootprints = list(),
                             facetBy = "sample",
                             adjustFacetHeight = TRUE,
                             referenceCoordinate = NULL,
                             labelAccuracy = NULL,
                             fillColors = "-cividis") {

    argL <- .checkArgsReadLevelPlots(
        se = se, region = region, assayName = assayName, drawRead = drawRead,
        orderReads = orderReads, orderRegion = orderRegion,
        clustDist = clustDist, windowWidth = windowWidth,
        modbaseSpace = modbaseSpace, trackTitle = trackTitle,
        legendTitle = legendTitle, yAxisLabel = yAxisLabel,
        showLegend = showLegend, highlightRegions = highlightRegions,
        footprintColumns = footprintColumns, arglistFootprints = arglistFootprints,
        facetBy = facetBy, adjustFacetHeight = adjustFacetHeight,
        referenceCoordinate = referenceCoordinate,
        labelAccuracy = labelAccuracy, linewidthTiles = linewidthTiles,
        interpolate = interpolate)

    # prepare plot data
    df <- .preparePlotdataReads(x = se, assayName = assayName,
                                modbaseSpace = modbaseSpace,
                                interpolate = interpolate,
                                referenceCoordinate = argL$referenceCoordinate,
                                extraColAnnots = setdiff(facetBy, "sample"),
                                orderReads = orderReads,
                                orderRegion = argL$orderRegion,
                                clustDist = clustDist,
                                windowWidth = windowWidth,
                                facetBy = facetBy)

    # create base plot
    p <- .createBaseplotReads(df = df, region = argL$region,
                              trackTitle = trackTitle,
                              legendTitle = ifelse(!is.null(legendTitle),
                                                   legendTitle, assayName),
                              yAxisLabel = yAxisLabel,
                              showLegend = showLegend,
                              highlightRegions = argL$highlightRegions,
                              facetBy = facetBy,
                              adjustFacetHeight  = adjustFacetHeight,
                              referenceCoordinate = argL$referenceCoordinate,
                              labelAccuracy = labelAccuracy,
                              fillColors = fillColors)

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(
            df, groupVars = union(facetBy, "sample"))
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["plotRow"]],
                                  xend = .data[["end"]]
                              ), color = "gray80")
    }

    # add tiles
    p <- p + geom_tile(color = "gray20", width = 1, height = 1,
                       linewidth = linewidthTiles)

    # add footprints
    for (fpc in footprintColumns) {
        fp <- .prepareFootprintsForPlot(fp = argL$footprints[[fpc]],
                                        plotdf = df, se = se, facetBy = facetBy)
        if (nrow(fp) > 0) {
            aL <- c(
                list(data = fp |> dplyr::filter(!is.na(.data$plotRow)),
                     mapping = aes(x = (as.numeric(.data$start) +
                                            as.numeric(.data$end)) / 2,
                                   y = .data$plotRow,
                                   width = abs(as.numeric(.data$end) -
                                                   as.numeric(.data$start))
                     ),
                     inherit.aes = FALSE),
                argL$arglistFootprints[[fpc]]
            )
            if (is.null(aL$height)) {
                aL$height <- 1
            }
            if (is.null(aL$fill)) {
                aL$fill <- "transparent"
            }
            if (is.null(aL$color) && is.null(aL$colour)) {
                aL$color <- argL$footprintColors[fpc]
            }
            if (is.null(aL$linewidth)) {
                aL$linewidth <- 1.5
            }
            repArgs <- duplicated(names(aL))
            if (any(repArgs)) {
                cli_warn("Ignoring pre-defined arguments: {.arg {names(aL)[repArgs]}}.")
            }
            aL <- aL[!repArgs]
            p <- p + do.call(geom_tile, aL)
        }
    }

    # return plot
    return(p)
}

#' @param doPoint A logical scalar. If \code{TRUE}, show points in the plot.
#' @param arglistPoint A list with arguments to be sent to
#'     \code{\link[ggplot2]{geom_point}}.
#' @param doSmooth A logical scalar. If \code{TRUE}, show a smooth line in the
#'     plot.
#' @param arglistSmooth A list with arguments to be sent to
#'     \code{\link[ggplot2]{geom_line}}.
#' @param smoothMethod A character scalar indicating the method to use for
#'     smoothing. Current options are \code{"smoothSpline"} and
#'     \code{"rollingMean"} (default, linear interpolation of values to the
#'     single-base pair level (unless \code{modbaseSpace} is \code{TRUE}),
#'     followed by rolling mean calculation).
#' @param spar A numeric scalar typically in (0,1] specifying the desired
#'     degree of smoothing if \code{smoothMethod} is \code{"smoothSpline"}
#'     (\code{spar} argument of \code{\link[stats]{smooth.spline}}).
#' @param windowSize A numeric scalar specifying the window size for smoothing
#'     if \code{smoothMethod} is \code{"rollingMean"}.
#' @param groupBy A character scalar indicating the sample annotation column
#'     to group the points by for creating smoothed lines. By default,
#'     the points will be grouped by 'sample', corresponding to the columns of
#'     \code{se}. Typically, the \code{groupBy} column should provide a
#'     finer (or identical) partition compared to the \code{colorBy} column -
#'     for example, grouping by sample and coloring by condition.
#' @param colorBy A character scalar indicating the sample annotation column
#'     to color the points and smoothed lines by. By default, the points and
#'     lines will be colored by 'sample', corresponding to the columns of
#'     \code{se}.
#' @param colors A named character vector of colors to use for the unique
#'     values in the \code{colorBy} annotation column. If \code{NULL}
#'     (default), the default \code{ggplot2} colors will be used.
#' @param yAxisRange Numeric vector of length 2 giving the range to zoom in
#'     to on the y-axis. If \code{NULL} (default), will be determined from the
#'     data.
#'
#' @export
#' @rdname plotRegion
#'
#' @examples
#' library(GenomicRanges)
#' bmfiles <- system.file("extdata",
#'                        c("modkit_pileup_1.bed.gz", "modkit_pileup_2.bed.gz"),
#'                        package = "footprintR")
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' seA <- readBedMethyl(bmfiles, modbase = "m",
#'                      sequenceContextWidth = 3, sequenceReference = reffile,
#'                      BPPARAM = BiocParallel::SerialParam())
#' plotSummaryPointSmooth(seA, region = as("chr1:6940000-6955000", "GRanges"),
#'                        assayName = "Nvalid", doPoint = FALSE)
#'
#' @import ggplot2
#' @importFrom dplyr group_by ungroup group_modify across all_of
#' @importFrom rlang .data
#' @importFrom stats smooth.spline
#' @importFrom GenomicRanges shift pintersect
#' @importFrom IRanges subsetByOverlaps
#' @importFrom SummarizedExperiment assayNames
#' @importFrom BiocGenerics intersect
#' @importFrom zoo na.approx rollmean
#' @importFrom cli cli_abort
#'
plotSummaryPointSmooth <- function(se,
                                   region,
                                   assayName,
                                   doPoint = TRUE,
                                   arglistPoint = list(),
                                   doSmooth = TRUE,
                                   arglistSmooth = list(),
                                   smoothMethod = "rollingMean",
                                   spar = 0.01,
                                   windowSize = 15,
                                   modbaseSpace = FALSE,
                                   trackTitle = NULL,
                                   legendTitle = NULL,
                                   yAxisLabel = assayName,
                                   showLegend = TRUE,
                                   highlightRegions = NULL,
                                   groupBy = "sample",
                                   colorBy = "sample",
                                   colors = NULL,
                                   referenceCoordinate = NULL,
                                   labelAccuracy = NULL,
                                   yAxisRange = NULL) {

    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = doPoint, type = "logical")
    .assertVector(x = arglistPoint, type = "list")
    .assertScalar(x = doSmooth, type = "logical")
    .assertVector(x = arglistSmooth, type = "list")
    .assertScalar(x = smoothMethod, type = "character",
                  validValues = c("smoothSpline", "rollingMean"))
    if (smoothMethod == "rollingMean") {
        .assertScalar(x = windowSize, type = "numeric")
        if (windowSize %% 2 != 1) {
            cli_abort("{.arg windowSize} ({windowSize}) must be an odd integer")
        }
    } else if (smoothMethod == "smoothSpline") {
        .assertScalar(x = spar, type = "numeric")
    }
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = yAxisLabel, type = "character")
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- GenomicRanges::pintersect(highlightRegions, region,
                                                      ignore.strand = TRUE,
                                                      drop.nohit.ranges = TRUE)
        highlightRegions <- highlightRegions[width(highlightRegions) > 0]
        # highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
        #                                             ignore.strand = TRUE)
    }
    .assertScalar(x = groupBy, type = "character", allowNULL = TRUE)
    .assertScalar(x = colorBy, type = "character", allowNULL = TRUE)
    .assertVector(x = colors, type = "character", allowNULL = TRUE)
    if (!is.null(colors)) {
        .assertVector(x = names(colors), type = "character")
        if (!is.null(colorBy) &&
                     !all(colData(se)[[colorBy]] %in% names(colors))) {
            cli_abort("Missing color specification for some values")
        }
    }
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)
    .assertVector(x = yAxisRange, type = "numeric", len = 2,
                  allowNULL = TRUE)
    if (modbaseSpace) {
        referenceCoordinate <- NULL
    }

    if (!is.null(referenceCoordinate)) {
        # shift all ranges
        if (!is.null(highlightRegions)) {
            highlightRegions <- shift(highlightRegions, -referenceCoordinate)
        }
        region <- shift(region, -referenceCoordinate)
    }

    # prepare plot data
    df <- .preparePlotdataSummary(x = se, assayName = assayName,
                                  modbaseSpace = modbaseSpace,
                                  referenceCoordinate = referenceCoordinate,
                                  extraColAnnots = setdiff(c(groupBy, colorBy),
                                                           "sample"))

    # create base plot
    p <- .createBaseplotSummary(df = df,
                                region = region,
                                trackTitle = trackTitle,
                                legendTitle = legendTitle,
                                showLegend = showLegend,
                                highlightRegions = highlightRegions,
                                groupBy = groupBy,
                                colorBy = colorBy,
                                colors = colors,
                                referenceCoordinate = referenceCoordinate,
                                labelAccuracy = labelAccuracy,
                                yAxisLabel = yAxisLabel,
                                yAxisRange = yAxisRange)

    # add points
    if (doPoint) {
        p <- p + do.call(geom_point, arglistPoint)
    }

    if (doSmooth) {
        # helper function to compute smooth spline for each sample
        compute_smooth <- function(data) {
            ok <- is.finite(data[["value"]])
            if (smoothMethod == "smoothSpline") {
                smooth <- smooth.spline(
                    x = data[["position"]][ok],
                    y = data[["value"]][ok],
                    keep.data = FALSE,
                    spar = spar)
                data.frame(position = smooth$x,
                           value_smooth = smooth$y)
            } else if (smoothMethod == "rollingMean") {
                ok <- is.finite(data[["value"]])
                # first interpolate linearly
                if (is.factor(data$position)) {
                    data$position <- as.numeric(data$position)
                }
                xint <- seq(from = min(data[["position"]]) - (windowSize - 1) / 2,
                            to = max(data[["position"]]) + (windowSize - 1) / 2,
                            by = 1)
                yint <- rep(NA, length(xint))
                idxout <- match(data[["position"]][ok], xint)
                yint[idxout] <- data[["value"]][ok]
                yint[1] <- yint[idxout[1]]
                yint[length(yint)] <- yint[idxout[length(idxout)]]
                yint <- na.approx(yint, maxgap = Inf)
                # then take rolling mean
                yint <- rollmean(yint, k = windowSize,
                                 fill = c(yint[1], NA, yint[length(yint)]))
                idxout <- match(data[["position"]], xint)
                data.frame(position = xint[idxout],
                           value_smooth = yint[idxout])
            }
        }

        # apply the function to each group/color combination
        smooth_data <- df |>
            group_by(across(all_of(c(groupBy, colorBy)))) |>
            group_modify(~ compute_smooth(.x)) |>
            ungroup()

        # add the smoothed line
        arglistSmooth <- arglistSmooth[!names(arglistSmooth) %in%
                                           c("data", "inherit.aes",
                                             "mapping")]
        p <- p + do.call(geom_line,
                         c(list(data = smooth_data, inherit.aes = FALSE,
                                mapping = aes(x = .data[["position"]],
                                              y = .data[["value_smooth"]],
                                              group = .data[[groupBy]],
                                              color = .data[[colorBy]])),
                           arglistSmooth))
    }

    # return the plot
    return(p)
}


#' @param grl A named \code{\link[GenomicRanges]{GRangesList}} object where each
#'     entry corresponds to a transcript or genomic feature.
#' @param colorByStrand A logical scalar indicating whether or not to color
#'     features by strand.
#' @param displayNames A logical scalar indicating whether or not to display
#'     the names of the features in the plot.
#' @param labelSize A numeric scalar representing the font size of the displayed
#'     label (if \code{displayNames} is \code{TRUE}).
#' @param labelPosition A character scalar, either \code{"above"},
#'     \code{"below"} or \code{"inside"}, indicating whether to place the
#'     feature labels above, below or inside the respective feature.
#'
#' @export
#' @rdname plotRegion
#'
#' @examples
#' library(GenomicRanges)
#' plotGenomicRegions(grl = GRangesList(
#'     g1 = GRanges("chr1", IRanges(c(10, 30), c(20, 35)), "+"),
#'     cgi1 = GRanges("chr1", IRanges(15, 25), "*"),
#'     g2 = GRanges("chr1", IRanges(c(15, 25), c(20, 40)), "-")),
#'     region = as("chr1:1-50", "GRanges"),
#'     labelPosition = "inside",
#'     labelSize = 5)
#'
#' @import ggplot2
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BiocGenerics unlist start end
#' @importFrom S4Vectors mcols
#' @importFrom cli cli_abort
#'
plotGenomicRegions <- function(grl,
                               region,
                               colorByStrand = TRUE,
                               displayNames = TRUE,
                               labelSize = 3,
                               labelPosition = "above",
                               trackTitle = NULL,
                               legendTitle = NULL,
                               showLegend = TRUE,
                               referenceCoordinate = NULL,
                               labelAccuracy = NULL) {

    # check input arguments
    .assertVector(x = grl, type = "GRangesList")
    .assertVector(x = names(grl), type = "character")
    if (any(is.na(names(grl)))) {
        cli_abort("{.code NA} values are not allowed in {.code names(grl)}")
    }
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = colorByStrand, type = "logical")
    .assertScalar(x = displayNames, type = "logical")
    .assertScalar(x = labelSize, type = "numeric")
    .assertScalar(x = labelPosition, type = "character",
                  validValues = c("above", "below", "inside"))
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = showLegend, type = "logical")
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)

    # subset GRangesList to elements overlapping the provided region
    grl <- subsetByOverlaps(grl, region)

    # create two flattened objects - one with the full range of each feature,
    # and one with the individual building blocks
    fname_full <- names(grl)
    fname_unique_full <- make.unique(names(grl))
    fname_parts <- rep(fname_full, lengths(grl))
    fname_unique_parts <- rep(fname_unique_full, lengths(grl))

    fullRange <- unlist(range(grl), use.names = FALSE)
    mcols(fullRange)$fpname <- fname_full
    mcols(fullRange)$fpname_unique <- fname_unique_full
    fullRange <- as.data.frame(fullRange)
    fullRange$fpname_unique <- factor(fullRange$fpname_unique,
                                      levels = fname_unique_full)

    rangeParts <- unlist(grl, use.names = FALSE)
    mcols(rangeParts)$fpname <- fname_parts
    mcols(rangeParts)$fpname_unique <- fname_unique_parts
    rangeParts <- as.data.frame(rangeParts)
    rangeParts$fpname_unique <- factor(rangeParts$fpname_unique,
                                       levels = fname_unique_full)

    rng <- c(start(region) - 0.5, end(region) + 0.5)

    if (!is.null(referenceCoordinate)) {
        fullRange$start <- fullRange$start - referenceCoordinate
        fullRange$end <- fullRange$end - referenceCoordinate
        rangeParts$start <- rangeParts$start - referenceCoordinate
        rangeParts$end <- rangeParts$end - referenceCoordinate
        rng <- rng - referenceCoordinate
    }

    # plot
    gg <- ggplot() +
        geom_segment(data = fullRange,
                     mapping = aes(
                         x = .data[["start"]],
                         y = .data[["fpname_unique"]],
                         xend = .data[["end"]]
                     ), color = "gray80")
    if (colorByStrand) {
        gg <- gg +
            geom_rect(data = rangeParts,
                      mapping = aes(
                          xmin = .data[["start"]],
                          xmax = .data[["end"]],
                          ymin = as.numeric(.data[["fpname_unique"]]) - 0.25,
                          ymax = as.numeric(.data[["fpname_unique"]]) + 0.25,
                          fill = .data[["strand"]]
                      ), color = "gray20") +
            scale_fill_manual(values = c("+" = "#82b579",
                                         "-" = "#c79e9d",
                                         "*" = "gray80"),
                              breaks = c("+", "-", "*"),
                              labels = c("+", "-", ""))
    } else {
        gg <- gg +
            geom_rect(data = rangeParts,
                      mapping = aes(
                          xmin = .data[["start"]],
                          xmax = .data[["end"]],
                          ymin = as.numeric(.data[["fpname_unique"]]) - 0.25,
                          ymax = as.numeric(.data[["fpname_unique"]]) + 0.25
                      ), color = "gray20", fill = "gray80")
    }
    if (displayNames) {
        offset <- ifelse(labelPosition == "above", 0.25,
                         ifelse(labelPosition == "below", -0.25, 0))
        vjust <- ifelse(labelPosition == "above", -0.5,
                        ifelse(labelPosition == "below", 1.5, 0.5))
        gg <- gg +
            geom_text(
                data = fullRange,
                mapping = aes(
                    x = ifelse(.data[["strand"]] == "+",
                               pmax(.data[["start"]], rng[1]),
                               ifelse(.data[["strand"]] == "-",
                                      pmin(.data[["end"]], rng[2]),
                                      0.5 * pmax(.data[["start"]], rng[1]) +
                                          0.5 * pmin(.data[["end"]], rng[2]))),
                    y = as.numeric(.data[["fpname_unique"]]) + offset,
                    label = .data[["fpname"]],
                    vjust = vjust,
                    hjust = ifelse(.data[["strand"]] == "+", -0.1,
                                   ifelse(.data[["strand"]] == "-", 1.1, 0.5))
                ),
                size = labelSize)
    }
    gg <- gg +
        labs(title = trackTitle,
             fill = ifelse(!is.null(legendTitle), legendTitle, "strand"),
             x = ifelse(is.null(referenceCoordinate),
                        paste0("Position on ",
                               as.character(seqnames(region))),
                        paste0("Position relative to ",
                               as.character(seqnames(region)), ":",
                               referenceCoordinate))) +
        theme_bw() +
        theme(legend.position = ifelse(showLegend, "right", "none"),
              legend.text = element_text(size = 16),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title.y = element_blank(),
              panel.border = element_blank(),
              axis.line.x = element_line(color = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    if (is.null(labelAccuracy)) {
        labelAccuracy <- 10^round(log10((rng[2] - rng[1]) / max(abs(rng))))
    }
    gg <- gg + coord_cartesian(xlim = rng) +
        scale_x_continuous(
            expand = c(0, 0),
            labels = label_number(
                accuracy = labelAccuracy,
                scale_cut = c(0, ` Kb` = 1000, ` Mb` = 1e+06, ` Gb` = 1e+9)))

    gg
}

## helper functions used above -------------------------------------------------

#' Create a ggplot2 fill scale based on \code{fillColors} argument
#'
#' @importFrom ggplot2 scale_fill_viridis_c scale_fill_gradientn
#' @keywords internal
#' @noRd
.createFillScale <- function(fillColors) {
    # check arguments
    .assertVector(x = fillColors, type = "character", rngLen = c(1, Inf))

    scl <- NULL
    if (identical(length(fillColors), 1L)) {
        # scale_fill_viridis_c
        # ... extract direction
        direction <- 1
        if (identical(substr(fillColors, 1, 1), "-")) {
            direction <- -1
            fillColors <- substr(fillColors, 2, nchar(fillColors))
        }

        # ... create scale
        scl <- scale_fill_viridis_c(begin = 0, end = 1,
                                    option = fillColors,
                                    direction = direction,
                                    na.value = "beige")
    } else {
        # scale_fill_gradientn
        scl <- scale_fill_gradientn(colors = fillColors)
    }

    return(scl)
}

#' Check arguments for read-level plot functions, and generate required objects
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom SummarizedExperiment assayNames
#' @importFrom IRanges ranges
#' @importFrom GenomicRanges shift pintersect
#' @importFrom BiocGenerics intersect
#' @importFrom S4Vectors endoapply
#' @importFrom cli cli_abort
.checkArgsReadLevelPlots <- function(se, region, assayName, drawRead,
                                     orderReads, orderRegion,
                                     clustDist, windowWidth, modbaseSpace,
                                     trackTitle, legendTitle, yAxisLabel,
                                     showLegend, highlightRegions,
                                     footprintColumns, arglistFootprints,
                                     facetBy, adjustFacetHeight,
                                     referenceCoordinate, labelAccuracy,
                                     size = 0, stroke = 0,
                                     linewidthTiles = 0, interpolate = FALSE) {
    # check arguments
    # ... shared arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = drawRead, type = "logical")
    .assertScalar(x = orderReads, type = "character", allowNULL = TRUE,
                  validValues = c("cluster", "squish", "regionAvg"))
    .assertScalar(x = orderRegion, type = "GRanges", allowNULL = TRUE)
    .assertScalar(x = clustDist, type = "character",
                  validValues = c("pearson", "euclidean"))
    .assertScalar(x = windowWidth, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = yAxisLabel, type = "character")
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges", allowNULL = TRUE)
    .assertVector(x = footprintColumns, type = "character", allowNULL = TRUE,
                  validValues = .getReadLevelColDataNames(se))
    .assertVector(x = arglistFootprints, type = "list", allowNULL = TRUE)
    if (length(arglistFootprints) > 0) {
        .assertVector(x = names(arglistFootprints), type = "character")
    }
    .assertScalar(x = facetBy, type = "character", allowNULL = TRUE)
    .assertScalar(x = adjustFacetHeight, type = "logical")
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)

    # ... lollipop-specific arguments
    .assertScalar(x = size, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = stroke, type = "numeric", rngIncl = c(0, Inf))

    # ... heatmap-specific arguments
    .assertScalar(x = linewidthTiles, type = "numeric")
    .assertScalar(x = interpolate, type = "logical")

    # adjust arguments if necessary
    if (is.null(orderRegion)) {
        orderRegion <- region
    }
    if (!is.null(highlightRegions)) {
        highlightRegions <- GenomicRanges::pintersect(highlightRegions, region,
                                                      ignore.strand = TRUE,
                                                      drop.nohit.ranges = TRUE)
        highlightRegions <- highlightRegions[width(highlightRegions) > 0]
        # highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
        #                                             ignore.strand = TRUE)
    }
    if (modbaseSpace) {
        referenceCoordinate <- NULL
    }
    if (!is.null(footprintColumns)) {
        footprintColors <- defaultFootprintColors[seq_along(footprintColumns)]
        names(footprintColors) <- footprintColumns

        footprints <- lapply(
            structure(footprintColumns,
                      names = footprintColumns),
            function(nm) lapply(se[[nm]], function(y) {
                endoapply(y, function(z) {
                    BiocGenerics::intersect(z, ranges(region))
                })
            }))
        if (!any(footprintColumns %in% names(arglistFootprints))) {
            # assume that arglistFootprints is a list of arguments for plotting
            arglistFootprints <- lapply(
                setNames(footprintColumns, footprintColumns),
                function(x) arglistFootprints)
        } else {
            if (any(!names(arglistFootprints) %in% footprintColumns)) {
                cli_abort(paste0("Can't unambiguously interpret ",
                                 "{.code names(arglistFootprints)} as either ",
                                 "function arguments or footprint column names."))
            }
            for (mf in setdiff(footprintColumns,
                               names(arglistFootprints))) {
                arglistFootprints[[mf]] <- list()
            }
            are_lists <- vapply(arglistFootprints, is.list, FALSE)
            if (any(!are_lists)) {
                cli_abort(paste0("{.arg arglistFootprints} entries must be ",
                                 "lists, the following are not: ",
                                 "{.code {names(arglistFootprints)[!are_lists]}}."))
            }
        }
    } else {
        footprintColors <- NULL
        footprints <- NULL
    }
    if (!is.null(referenceCoordinate)) {
        # shift all ranges
        if (!is.null(highlightRegions)) {
            highlightRegions <- shift(highlightRegions, -referenceCoordinate)
        }
        if (!is.null(footprintColumns)) {
            footprints <- lapply(
                footprints,
                function(x) lapply(x,
                                   function(y) shift(y, -referenceCoordinate)))
        }
        region <- shift(region, -referenceCoordinate)
    }

    # return possibly adjusted arguments
    return(list(highlightRegions = highlightRegions,
                footprintColors = footprintColors,
                arglistFootprints = arglistFootprints,
                referenceCoordinate = referenceCoordinate,
                footprints = footprints, region = region,
                orderRegion = orderRegion))

}

#' Create data.frame from SummarizedExperiment for summary-level data
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the "position"
#'     column in the return data frame is categorical, instead of giving
#'     the numeric position in the genome.
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used.
#' @param extraColAnnots A character vector (or \code{NULL}) with names of
#'     columns in \code{colData(x)} to add to the generated data frame.
#'
#' @importFrom BiocGenerics start colnames rownames
#' @importFrom SummarizedExperiment assay colData
#' @importFrom cli cli_abort
#' @importFrom dplyr left_join bind_cols
#'
#' @noRd
#' @keywords internal
.preparePlotdataSummary <- function(x,
                                    assayName,
                                    modbaseSpace = FALSE,
                                    referenceCoordinate = NULL,
                                    extraColAnnots = NULL) {
    assaydat <- assay(x, assayName)
    i <- which(is.finite(assaydat), arr.ind = TRUE)
    df <- data.frame(
        position = start(x)[i[,"row"]],
        sample = colnames(x)[i[,"col"]],
        value = assaydat[i])
    if (modbaseSpace) {
        df$position <- factor(df$position,
                              levels = unique(sort(df$position,
                                                   decreasing = FALSE)))
    } else if (!is.null(referenceCoordinate)) {
        df$position <- df$position - referenceCoordinate
    }
    if (length(extraColAnnots) != 0) {
        colExists <- extraColAnnots %in% colnames(colData(x))
        if (!all(colExists)) {
            cli_abort("Some requested columns are not present in {.code colData(x)}")
        }
        df <- df |>
            left_join(bind_cols(
                sample = rownames(colData(x)),
                as.data.frame(colData(x)[, extraColAnnots, drop = FALSE])
            ), by = "sample")
    }
    return(df)
}

#' Create data.frame from SummarizedExperiment for read-level data
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data (positions in rows and reads in
#'     columns).
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the "position"
#'     column in the return data frame is categorical, instead of giving
#'     the numeric position in the genome.
#' @param interpolate A logical scalar. If \code{TRUE}, the gaps between
#'     observations are filled in by linear interpolation.
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used.
#' @param extraColAnnots A character vector (or \code{NULL}) with names of
#'     columns in \code{colData(x)} to add to the generated data frame.
#'
#' @importFrom BiocGenerics start colnames rownames unlist sort
#' @importFrom IRanges IRanges disjointBins
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SparseArray nnawhich nnavals colSums is_nonna
#' @importFrom S4Vectors endoapply
#' @importFrom cli cli_abort
#' @importFrom dplyr left_join bind_cols group_by summarise group_split filter
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.preparePlotdataReads <- function(x,
                                  assayName,
                                  modbaseSpace = FALSE,
                                  interpolate = FALSE,
                                  referenceCoordinate = NULL,
                                  extraColAnnots = NULL,
                                  orderReads = "cluster",
                                  orderRegion = NULL,
                                  clustDist = "euclidean",
                                  windowWidth = 25,
                                  facetBy = NULL) {
    assaydat <- assay(x, assayName)
    assaydat <- .removeAllNAReads(assaydat, prune = TRUE)
    # `assayName` columns are grouped reads -> flatten
    sample_ids <- rep(colnames(assaydat), unlist(lapply(assaydat, ncol)))
    assaydat <- as.matrix(assaydat)
    if (interpolate) {
        assaydat <- .interpolateColumns(assaydat, start(x))
        df <- data.frame(
            position = attr(assaydat, "pos"),
            read = factor(rep(colnames(assaydat), each = nrow(assaydat)),
                          levels = colnames(assaydat)),
            sample = rep(sample_ids, each = nrow(assaydat)),
            value = as.vector(assaydat)) |>
            filter(!is.na(.data$value))
    } else {
        i <- nnawhich(assaydat, arr.ind = TRUE)
        df <- data.frame(
            position = start(x)[i[, 1]],
            read = factor(colnames(assaydat)[i[, 2]],
                          levels = colnames(assaydat)),
            sample = sample_ids[i[, 2]],
            value = nnavals(assaydat))
    }
    if (modbaseSpace) {
        df$position <- factor(df$position,
                              levels = unique(sort(df$position,
                                                   decreasing = FALSE)))
    } else if (!is.null(referenceCoordinate)) {
        df$position <- df$position - referenceCoordinate
    }
    if (length(extraColAnnots) != 0) {
        colExists <- extraColAnnots %in% colnames(colData(x))
        if (!all(colExists)) {
            cli_abort("Some requested columns are not present in {.code colData(x)}")
        }
        df <- df |>
            left_join(bind_cols(
                sample = rownames(colData(x)),
                as.data.frame(colData(x)[, extraColAnnots, drop = FALSE])
            ), by = "sample")
    }

    # order reads
    if (!is.null(orderReads) && orderReads %in% c("cluster", "regionAvg")) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(x = x, assayName = assayName,
                                               method = orderReads,
                                               windowWidth = windowWidth,
                                               orderRegion = orderRegion,
                                               clustDist = clustDist))
        df$plotRow <- df$read
    } else if (!is.null(orderReads) && orderReads == "squish") {
        if (!is.null(facetBy)) {
            tmp <- df |> group_by(.data[[facetBy]], .data$read) |>
                summarise(start = min(as.numeric(.data$position)),
                          end = max(as.numeric(.data$position)))
        } else {
            tmp <- df |> group_by(.data$read) |>
                summarise(start = min(as.numeric(.data$position)),
                          end = max(as.numeric(.data$position)))
        }
        tmp <- tmp |> group_split() |>
            lapply(function(y) sort(IRanges(start = y$start, end = y$end,
                                            names = as.character(y$read)))) |>
            as("IRangesList") |>
            disjointBins() |>
            unlist()
        df$plotRow <- factor(tmp[as.character(df$read)])
    } else if (is.null(orderReads)) {
        df$plotRow <- df$read
    }

    return(df)
}

#' Return a base ggplot2 plot (no geometries yet) for summary-level data
#'
#' @param df A \code{\link{data.frame}} with the plot data (typically
#'     created by \code{\link{.preparePlotdataSummary}}.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region.
#' @param trackTitle A character scalar or \code{NULL}, giving the title of
#'     the track.
#' @param legendTitle A character scalar (or \code{NULL}) giving the title for
#'     the color legend.
#' @param showLegend A logical scalar indicating whether or not to show the
#'     legend for the plot.
#' @param highlightRegions A \code{\link[GenomicRanges]{GRanges}} object
#'     containing regions to highlight with a grey shading.
#' @param groupBy A character scalar indicating the column
#'     to group the points by for creating smoothed lines.
#' @param colorBy A character scalar indicating the column
#'     to color the points and smoothed lines by.
#' @param colors A named character vector of colors to use for the unique
#'     values in the \code{colorBy} annotation column. If \code{NULL}
#'     (default), the default \code{ggplot2} colors will be used.
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used.
#' @param labelAccuracy A numeric scalar indicating the precision of the
#'     positions along the genomic axis. Will be passed to
#'     \code{\link[scales]{label_number}}. If \code{NULL} (default), a suitable
#'     value will be derived from \code{region}.
#' @param yAxisLabel A character scalar providing the label to use for the
#'     y-axis.
#' @param yAxisRange Numeric vector of length 2 giving the range to zoom in
#'     to on the y-axis. If \code{NULL} (default), will be determined from the
#'     data.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.createBaseplotSummary <- function(df,
                                   region,
                                   trackTitle,
                                   legendTitle,
                                   showLegend,
                                   highlightRegions,
                                   groupBy,
                                   colorBy,
                                   colors,
                                   referenceCoordinate,
                                   labelAccuracy,
                                   yAxisLabel,
                                   yAxisRange) {
    if (is.null(yAxisRange)) {
        # set the y-axis range manually to make sure it stays consistent
        # with/without points
        yAxisRange <- range(df$value) + c(-1, 1) * 0.05 * diff(range(df$value))
    }
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["value"]],
                      group = .data[[groupBy]],
                      color = .data[[colorBy]],
                      fill = .data[[colorBy]])) +
        labs(x = ifelse(is.numeric(df$position),
                        ifelse(is.null(referenceCoordinate),
                               paste0("Position on ",
                                      as.character(seqnames(region))),
                               paste0("Position relative to ",
                                      as.character(seqnames(region)), ":",
                                      referenceCoordinate)),
                        paste0("Modified positions in ",
                               as.character(seqnames(region)),
                               ":", levels(df$position)[1], "-",
                               levels(df$position)[nlevels(df$position)])),
             y = yAxisLabel,
             color = ifelse(!is.null(legendTitle), legendTitle, colorBy),
             fill = ifelse(!is.null(legendTitle), legendTitle, colorBy),
             title = trackTitle) +
        theme_bw() +
        theme(legend.position = ifelse(showLegend, "right", "none"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())

    if (is.factor(df$position)) {
        p0 <- p0 + theme(axis.text.x = element_blank()) +
            coord_cartesian(ylim = yAxisRange) +
            scale_x_discrete(expand = expansion(mult = 0, add = 0.5))
    } else {
        p0 <- .addCoordAxisFormat(p0 = p0, region = region,
                                  labelAccuracy = labelAccuracy,
                                  yAxisRange = yAxisRange)
    }

    if (!is.null(colors)) {
        p0 <- p0 + scale_color_manual(values = colors) +
            scale_fill_manual(values = colors)
    }

    if (!is.null(highlightRegions)) {
        dfhr <- data.frame(highlightRegions)
        dfhr <- .convertRegionToModBaseSpace(regdf = dfhr, datadf = df)
        if (nrow(dfhr) > 0) {
            p0 <- p0 +
                geom_rect(
                    data = dfhr,
                    mapping = aes(xmin = as.numeric(start) - 0.0,
                                  xmax = as.numeric(end) + 0.0,
                                  ymin = -Inf, ymax = Inf),
                    fill = "gray90",
                    inherit.aes = FALSE
                )
        }
    }

    return(p0)
}

#' Return a base ggplot2 plot (no geometries yet) for read-level data
#'
#' @param df A \code{\link{data.frame}} with the plot data (typically
#'     created by \code{\link{.preparePlotdataReads}}.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region.
#' @param trackTitle A character scalar or \code{NULL}, giving the title of
#'     the track.
#' @param legendTitle A character scalar (or \code{NULL}) giving the title for
#'     the fill legend.
#' @param showLegend A logical scalar indicating whether or not to show the
#'     legend for the plot.
#' @param highlightRegions A \code{\link[GenomicRanges]{GRanges}} object
#'     containing regions to highlight with a grey shading.
#' @param facetBy A character scalar indicating the sample annotation column
#'     to facet the plot by (if \code{NULL}, no faceting is done).
#' @param adjustFacetHeight A logical scalar. If \code{TRUE}, adjust the
#'     height of the facets by the number of reads in each of them. If
#'     \code{FALSE}, all facets have the same height.
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used.
#' @param labelAccuracy A numeric scalar indicating the precision of the
#'     positions along the genomic axis. Will be passed to
#'     \code{\link[scales]{label_number}}. If \code{NULL} (default), a suitable
#'     value will be derived from \code{region}.
#' @param fillColors A character scalar defining the continuous color palette to
#'     represent modification probabilities. If \code{fillColors}
#'     has length one, it is assumed to be the a supported value to pass to the
#'     \code{option} argument of \code{\link[ggplot2]{scale_fill_viridis_c}},
#'     optionally prefixed with a minus sign to in addition set
#'     \code{direction = -1}). If \code{fillColors} has more than one element,
#'     it is assumed to be a vector of colors to pass to the \code{colors}
#'     argument of \code{\link[ggplot2]{scale_fill_gradientn}}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom ggforce facet_col
#'
#' @noRd
#' @keywords internal
.createBaseplotReads <- function(df,
                                 region,
                                 trackTitle,
                                 legendTitle,
                                 yAxisLabel,
                                 showLegend,
                                 highlightRegions,
                                 facetBy,
                                 adjustFacetHeight,
                                 referenceCoordinate,
                                 labelAccuracy,
                                 fillColors = "-cividis") {
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["plotRow"]],
                      fill = .data[["value"]])) +
        .createFillScale(fillColors) +
        labs(x = ifelse(is.numeric(df$position),
                        ifelse(is.null(referenceCoordinate),
                               paste0("Position on ",
                                      as.character(seqnames(region))),
                               paste0("Position relative to ",
                                      as.character(seqnames(region)), ":",
                                      referenceCoordinate)),
                        paste0("Modified positions in ",
                               as.character(seqnames(region)),
                               ":", levels(df$position)[1], "-",
                               levels(df$position)[nlevels(df$position)])),
             y = yAxisLabel,
             fill = legendTitle,
             title = trackTitle) +
        theme_bw() +
        theme(legend.position = ifelse(showLegend, "right", "none"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background.x = element_blank(),
              strip.text.x = element_text(
                  hjust = 0, margin = margin(t = 0, r = 0, b = 2, l = 0)))

    if (!is.null(facetBy)) {
        if (adjustFacetHeight) {
            # adjust facet height to the number of reads
            p0 <- p0 +
                facet_col(~ .data[[facetBy]], scales = "free_y", space = "free")
        } else {
            p0 <- p0 +
                facet_wrap(~ .data[[facetBy]], ncol = 1, scales = "free_y")
        }
    }
    if (is.factor(df$position)) {
        p0 <- p0 + theme(axis.text.x = element_blank())
    } else {
        p0 <- .addCoordAxisFormat(p0 = p0, region = region,
                                  labelAccuracy = labelAccuracy)
    }

    if (!is.null(highlightRegions)) {
        dfhr <- data.frame(highlightRegions)
        dfhr <- .convertRegionToModBaseSpace(regdf = dfhr, datadf = df)
        if (nrow(dfhr) > 0) {
            p0 <- p0 +
                geom_rect(
                    data = dfhr,
                    mapping = aes(xmin = as.numeric(start) - 0.0,
                                  xmax = as.numeric(end) + 0.0,
                                  ymin = -Inf, ymax = Inf),
                    fill = "gray90",
                    inherit.aes = FALSE
                )
        }
    }

    return(p0)
}

#' Prepare footprint coordinates for plotting
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr left_join bind_cols rename mutate
#' @importFrom SummarizedExperiment colData
.prepareFootprintsForPlot <- function(fp, plotdf, se, facetBy) {
    # unpack list
    fp <- lapply(fp, function(x) {
        as.data.frame(unlist(x, use.names = TRUE))
    })
    fp <- cbind(do.call(
        rbind, fp),
        sample = rep(names(fp), vapply(fp, nrow, 0))) |>
        rename(read = "names") |>
        mutate(read = factor(.data$read, levels = levels(plotdf$read)))
    fp$plotRow <- plotdf$plotRow[match(fp$read, plotdf$read)]
    fp <- .convertRegionToModBaseSpace(
        regdf = fp, datadf = plotdf)
    # add facetting variable
    if (!is.null(facetBy) && facetBy != "sample") {
        fp <- fp |>
            left_join(bind_cols(
                sample = rownames(colData(se)),
                as.data.frame(colData(se)[, facetBy, drop = FALSE])
            ), by = "sample")
    }
    return(fp)
}

#' Per-read summarize a read-level data.frame
#'
#' @param dfReads A \code{data.frame} object to summarize, typically generated
#'     by \code{\link{.preparePlotdataReads}}.
#'
#' @importFrom dplyr group_by summarise across all_of
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.summarizePlotdataPerRead <- function(dfReads, groupVars = "sample") {
    dfReads |>
        group_by(.data[["read"]], .data[["plotRow"]]) |>
        summarise(
            start = ifelse(is.factor(.data[["position"]]),
                           levels(.data[["position"]])[1],
                           min(.data[["position"]])),
            end = ifelse(is.factor(.data[["position"]]),
                         levels(.data[["position"]])[nlevels(.data[["position"]])],
                         max(.data[["position"]])),
            across(all_of(groupVars), unique),
            .groups = "drop")
}

#' Return ordered read identifiers
#'
#' @description
#' Returns ordered read identifiers (\code{colnames(x)} according to
#' \code{method}.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param method A character scalar with the ordering method to be used.
#'     Supported values are:
#'     \describe{
#'         \item{\code{"cluster"}}{, which orders reads using hierarchical
#'         clustering based on Pearson correlation distance
#'         (\code{as.dist(sqrt(2 - 2 * cor(X)))}, if
#'         \code{clustDist = "pearson"} or Euclidean distance (if
#'         \code{clustDist = "euclidean"}). The input to the distance calculation
#'         is \code{assay(x, assayName)} with zero values set to \code{NA} and
#'         averaged over windows of \code{windowWidth} nucleotides. The
#'         distance calculations can be further limited to a specific region
#'         by specifying \code{orderRegion}.}
#'         \item{\code{"regionAvg"}}{, which orders reads increasingly by the
#'         average modification probability in \code{orderRegion}.}
#'     }
#' @param windowWidth A numeric scalar giving the window width for which
#'     read-level data will be averaged before clustering. This should help
#'     to reduce the noise and allows to compare reads without any common
#'     modification calls, such as plus- and minus-strand reads with 6mA calls.
#' @param orderRegion Either \code{NULL} or a length-one \code{GRanges} object.
#'     If \code{method = "regionAvg"}, the \code{GRanges} object defines the
#'     window in which average modification probability is calculated to order
#'     the reads in read-level plots. If \code{method = "cluster"}, the
#'     object defines the region within which the clustering is calculated.
#'     A \code{NULL} value indicates that the entire plotted region should be
#'     used as the window.
#' @param clustDist A character scalar defining the distance measure to use
#'     for clustering. Should be one of \code{"pearson"} or
#'     \code{"euclidean"}.
#'
#' @importFrom BiocGenerics colnames start ncol
#' @importFrom SummarizedExperiment assay rowRanges
#' @importFrom IRanges overlapsAny
#' @importFrom stats cor as.dist hclust
#' @importFrom SparseArray colMeans
#'
#' @noRd
#' @keywords internal
.orderReads <- function(x,
                        assayName,
                        method = c("cluster", "regionAvg"),
                        windowWidth = 25,
                        orderRegion = NULL,
                        clustDist = c("pearson", "euclidean")) {
    method <- match.arg(method)
    clustDist <- match.arg(clustDist)

    # extract and flatten assay matrix
    X <- as.matrix(assay(x, assayName))
    res <- colnames(X)

    # select rows that fall into orderRegion
    sel <- which(overlapsAny(query = rowRanges(x), subject = orderRegion,
                             ignore.strand = TRUE))

    if (identical(method, "cluster")) {
        if (ncol(X) > 1) {
            # group positions into bins of windowWidth
            bin <- findInterval(
                x = start(x)[sel],
                vec = seq(from = min(start(x)[sel]),
                          to = ceiling(max(end(x)[sel]) / windowWidth) * windowWidth + 1,
                          by = windowWidth),
                rightmost.closed = TRUE, left.open = FALSE)
            iByBin <- split(seq.int(length(sel)), bin)
            XX <- do.call(rbind, lapply(iByBin, function(i) {
                colMeans(X[sel[i], , drop = FALSE], na.rm = TRUE)
            }))
            # calculate distances between reads
            D <- .calcDist(XX, clustDist = clustDist)
            # cluster reads and return order
            cl <- hclust(D, method = "ward.D2")
            res <- colnames(X)[cl$order]
        }
    } else if (identical(method, "regionAvg")) {
        avg <- colMeans(X[sel, , drop = FALSE], na.rm = TRUE)
        res <- colnames(X)[order(avg, na.last = TRUE, decreasing = TRUE)]
    }
    return(res)
}

#' Calculate distances between columns of matrix
#'
#' @noRd
#' @keywords internal
.calcDist <- function(X, clustDist = "euclidean") {
    if (clustDist == "pearson") {
        D <- as.dist(sqrt(2 - 2 * cor(X, method = "pearson",
                                      use = "pairwise.complete")))
        D[is.na(D)] <- 2.0
    } else if (clustDist == "euclidean") {
        D <- dist(t(X), method = "euclidean")
        D[is.na(D)] <- sqrt(nrow(X))
    } else if (clustDist == "cosine") {
        D <- as.dist(1 - (t(X) %*% X) /
                         sqrt(cbind(colSums(X ^ 2)) %*% rbind(colSums(X ^ 2))))
        # in principle, the maximal distance is 1 since the values in X are
        # non-negative - in general though, cosine distances are in [0, 2]
        D[is.na(D)] <- 2.0
    }
    D
}

#' Add formatting for base-space x-axis to ggplot object
#'
#' @description
#' Add limits and axis label formatting for a base-space x-axis with
#' genomic coordinates to a ggplot object.
#'
#' @param p0 A \code{ggplot} object to which to add the axis formatting.
#'     x-axis values are assumed to be in \code{p0$data$position}.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region.
#' @param labelAccuracy The desired accuracy of the labels - if \code{NULL} it
#'     will be automatically determined.
#' @param yAxisRange Numeric vector of length 2 giving the range to zoom in
#'     to on the y-axis. If \code{NULL} (default), will be determined from the
#'     data.
#'
#' @import ggplot2
#' @importFrom scales label_number
#'
#' @noRd
#' @keywords internal
.addCoordAxisFormat <- function(p0, region, labelAccuracy, yAxisRange = NULL) {
    rng <- c(start(region) - 0.5, end(region) + 0.5)
    if (is.null(labelAccuracy)) {
        labelAccuracy <- 10^round(log10((rng[2] - rng[1]) / max(abs(rng))))
    }
    p0 <- p0 + coord_cartesian(xlim = rng, ylim = yAxisRange) +
        scale_x_continuous(
            expand = c(0, 0),
            labels = label_number(
                accuracy = labelAccuracy,
                scale_cut = c(0, ` Kb` = 1000, ` Mb` = 1e+06, ` Gb` = 1e+9)))
    return(p0)
}
