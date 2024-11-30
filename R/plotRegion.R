# global data.frame of plot types and characteristics
plotRegionPlotTypes <- data.frame(
    name = c("Point", "Smooth", "PointSmooth",
             "Lollipop", "Heatmap", "GenomicRegion"),
    type = c("summary", "summary", "summary",
             "reads", "reads", "annotation")
)


#' Plot single-molecule footprinting data for a single genomic region.
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
#'         \item{\code{"GenomicRegion"}}{: Genomic annotations (e.g.,
#'             transcripts, peaks, CpG islands).}
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
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop")))
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
        labelAccuracy = NULL) {

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
            cli_abort("tracks[[{i}]] has to be a list of length >=2.")
        }
        if (is.null(names(tracks[[i]])) || any(names(tracks[[i]]) == "") ||
            any(!c("trackData", "trackType") %in% names(tracks[[i]]))) {
            cli_abort(paste("tracks[[{i}]] must be a named list, and contain at",
                            "least entries named 'trackData' and 'trackType'"))
        }
        if (!tracks[[i]]$trackType %in% plotRegionPlotTypes$name) {
            cli_abort(paste("tracks[[{i}]]$trackType must be one of: ",
                            paste(plotRegionPlotTypes$name, collapse = ",")))
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
            cli_abort(paste("tracks[[{i}]]$trackData must be a character scalar",
                            "corresponding to a name of an assay in se"))
        }
        if (type_i == "reads" &&
            !tracks[[i]]$trackData %in% .getReadLevelAssayNames(se)) {
            cli_abort(paste("tracks[[{i}]]$trackData must be the name of a",
                            "read-level assay in se"))
        }
        if (type_i == "summary" &&
            tracks[[i]]$trackData %in% .getReadLevelAssayNames(se)) {
            cli_abort(paste("tracks[[{i}]]$trackData must be the name of a",
                            "summary assay in se"))
        }
        if (type_i %in% c("reads", "summary")) {
            # add assay to list of assays that are required for the plots
            assaysInUse <- union(assaysInUse, tracks[[i]]$trackData)
        }
        if (type_i == "annotation" &&
            !(is(tracks[[i]]$trackData, "GRangesList") &&
              !is.null(names(tracks[[i]]$trackData)))) {
            cli_abort(paste("tracks[[{i}]]$trackData must be a named",
                            "GRangesList object"))
        }
        if (type_i == "annotation" &&
            length(unlist(lapply(tracks[[i]]$trackData,
                                 function(y) unique(as.character(strand(y)))))) !=
            length(tracks[[i]]$trackData)) {
            cli_abort(paste("There are entries in tracks[[{i}]]$trackData",
                            "with mixed strand annotations"))
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
            tracks[[i]]$trackType == "GenomicRegion") {
            cli_warn(paste("Plotting in `modbaseSpace` is not allowed if",
                           "GenomicRegion tracks are included.",
                           "Setting modbaseSpace=FALSE"))
            modbaseSpace <- FALSE
        }
    }
    .assertVector(x = sequenceContext, type = "character", allowNULL = TRUE)
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)

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
    
    # subset se
    se <- subsetByOverlaps(x = se, ranges = region)
    se <- .keepPositionsBySequenceContext(
        se = se, sequenceContext = sequenceContext)

    if (nrow(se) == 0) {
        cli_abort("No positions retained for plotting!")
    }
    
    ## create plots
    pL <- vector("list", length = length(tracks))
    for (i in seq_along(tracks)) {
        tr <- tracks[[i]]
        trt <- plotRegionPlotTypes$type[
            match(tr$trackType, plotRegionPlotTypes$name)]
        if (trt %in% c("summary", "reads")) {
            args <- c(
                list(se = se, region = region, assayName = tr$trackData,
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
                list(grl = subsetByOverlaps(tr$trackData, region),
                     region = region, labelAccuracy = labelAccuracy,
                     referenceCoordinate = referenceCoordinate),
                     tr[!names(tr) %in% c("trackData", "trackType", "grl",
                                          "region", "referenceCoordinate",
                                          "labelAccuracy")]
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
            GenomicRegion = do.call(plotGenomicRegions, args)
        )
    }

    ## assemble composite plot
    if (length(pL) > 1L) { # suppress x-axis labels for all but last plot
        for (i in seq.int(length(pL) - 1L)) {
            pL[[i]] <- pL[[i]] + labs(x = element_blank())
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
#' @param orderReads A logical scalar. If \code{TRUE}, the position of reads
#'     on the y-axis will be reordered using \code{hclust(as.dist(1-cor(X)))$order},
#'     where \code{X} is \code{assay(x, assayName)} with zero values set to \code{NA}.
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
#' @param facetBy A character scalar indicating the sample annotation column 
#'     to facet the plot by (if \code{NULL}, no faceting is done). By default,
#'     the plot will be facetted by 'sample', corresponding to the columns of 
#'     \code{se}.
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
#' @importFrom SummarizedExperiment assayNames
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges shift
#' @importFrom BiocGenerics intersect
#'
plotReadsLollipop <- function(se,
                              region,
                              assayName,
                              size = 3.0,
                              stroke = 0.5,
                              drawRead = TRUE,
                              orderReads = TRUE,
                              modbaseSpace = FALSE,
                              trackTitle = NULL,
                              legendTitle = NULL,
                              showLegend = TRUE,
                              highlightRegions = NULL,
                              facetBy = "sample",
                              referenceCoordinate = NULL,
                              labelAccuracy = NULL) {
    
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = size, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = stroke, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = drawRead, type = "logical")
    .assertScalar(x = orderReads, type = "logical")
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
                                                    ignore.strand = TRUE)
    }
    .assertScalar(x = facetBy, type = "character", allowNULL = TRUE)
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)
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
    df <- .preparePlotdataReads(x = se, assayName = assayName,
                                modbaseSpace = modbaseSpace,
                                referenceCoordinate = referenceCoordinate, 
                                extraColAnnots = setdiff(facetBy, "sample"))

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(se, assayName))
    }

    # create base plot
    p <- .createBaseplotReads(df = df, region = region,
                              trackTitle = trackTitle,
                              legendTitle = ifelse(!is.null(legendTitle), legendTitle, assayName),
                              showLegend = showLegend,
                              highlightRegions = highlightRegions,
                              facetBy = facetBy,
                              referenceCoordinate = referenceCoordinate,
                              labelAccuracy = labelAccuracy)

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(
            df, groupVars = union(facetBy, "sample"))
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["read"]],
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
#' @importFrom SummarizedExperiment assayNames
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomicRanges shift
#' @importFrom BiocGenerics intersect
#'
plotReadsHeatmap <- function(se,
                             region,
                             assayName,
                             drawRead = TRUE,
                             linewidthTiles = 0,
                             orderReads = TRUE,
                             modbaseSpace = FALSE,
                             interpolate = FALSE,
                             trackTitle = NULL,
                             legendTitle = NULL,
                             showLegend = TRUE,
                             highlightRegions = NULL,
                             facetBy = "sample",
                             referenceCoordinate = NULL,
                             labelAccuracy = NULL) {
    
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = drawRead, type = "logical")
    .assertScalar(x = linewidthTiles, type = "numeric")
    .assertScalar(x = orderReads, type = "logical")
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertScalar(x = interpolate, type = "logical")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
                                                    ignore.strand = TRUE)
    }
    .assertScalar(x = facetBy, type = "character", allowNULL = TRUE)
    .assertScalar(x = referenceCoordinate, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = labelAccuracy, type = "numeric", allowNULL = TRUE)
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
    df <- .preparePlotdataReads(x = se, assayName = assayName,
                                modbaseSpace = modbaseSpace,
                                interpolate = interpolate,
                                referenceCoordinate = referenceCoordinate, 
                                extraColAnnots = setdiff(facetBy, "sample"))

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(se, assayName))
    }

    # create base plot
    p <- .createBaseplotReads(df = df, region = region,
                              trackTitle = trackTitle,
                              legendTitle = ifelse(!is.null(legendTitle), legendTitle, assayName),
                              showLegend = showLegend,
                              highlightRegions = highlightRegions,
                              facetBy = facetBy,
                              referenceCoordinate = referenceCoordinate,
                              labelAccuracy = labelAccuracy)

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(
            df, groupVars = union(facetBy, "sample"))
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["read"]],
                                  xend = .data[["end"]]
                              ), color = "gray80")
    }

    # add tiles
    p <- p + geom_tile(color = "gray20", width = 1, height = 1,
                       linewidth = linewidthTiles)

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
#'     \code{"rollingMean"} (linear interpolation of values to the single-base
#'     pair level (unless \code{modbaseSpace} is \code{TRUE}), followed by 
#'     rolling mean calculation).
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
#' @importFrom GenomicRanges shift
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
                                   smoothMethod = "smoothSpline",
                                   spar = 0.01,
                                   windowSize = 15,
                                   modbaseSpace = FALSE,
                                   trackTitle = NULL,
                                   legendTitle = NULL,
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
            cli_abort("windowSize must be an odd integer")
        }
    } else if (smoothMethod == "smoothSpline") {
        .assertScalar(x = spar, type = "numeric")
    }
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertScalar(x = trackTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = legendTitle, type = "character", allowNULL = TRUE)
    .assertScalar(x = showLegend, type = "logical")
    .assertVector(x = highlightRegions, type = "GRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- BiocGenerics::intersect(highlightRegions, region,
                                                    ignore.strand = TRUE)
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
                                yAxisLabel = assayName,
                                yAxisRange = yAxisRange)

    # add points
    if (!doPoint) {
        arglistPoint$color <- "transparent"
    }
    p <- p + do.call(geom_point, arglistPoint)
    
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
             fill = ifelse(!is.null(legendTitle), legendTitle, "strand")) +
        theme_bw() +
        theme(legend.position = ifelse(showLegend, "right", "none"),
              legend.text = element_text(size = 16),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.title = element_blank(),
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
            cli_abort("Some requested columns are not present in colData(x)")
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
#' @importFrom BiocGenerics start colnames rownames
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SparseArray nnawhich nnavals colSums is_nonna
#' @importFrom S4Vectors endoapply
#' @importFrom cli cli_abort
#' @importFrom dplyr left_join bind_cols
#'
#' @noRd
#' @keywords internal
.preparePlotdataReads <- function(x,
                                  assayName,
                                  modbaseSpace = FALSE,
                                  interpolate = FALSE,
                                  referenceCoordinate = NULL,
                                  extraColAnnots = NULL) {
    assaydat <- assay(x, assayName)
    assaydat <- .removeAllNAReads(assaydat, prune = TRUE)
    # `assayName` columns are grouped reads -> flatten
    sample_ids <- rep(colnames(x), unlist(lapply(assaydat, ncol)))
    assaydat <- as.matrix(assaydat)
    if (interpolate) {
        assaydat <- .interpolateColumns(assaydat, start(x))
        df <- data.frame(
            position = attr(assaydat, "pos"),
            read = factor(rep(colnames(assaydat), each = nrow(assaydat)),
                          levels = colnames(assaydat)),
            sample = rep(sample_ids, each = nrow(assaydat)),
            value = as.vector(assaydat))
    } else {
        i <- nnawhich(assaydat, arr.ind = TRUE)
        df <- data.frame(
            position = start(x)[i[,1]],
            read = factor(colnames(assaydat)[i[,2]], levels = colnames(assaydat)),
            sample = sample_ids[i[,2]],
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
            cli_abort("Some requested columns are not present in colData(x)")
        }
        df <- df |>
            left_join(bind_cols(
                sample = rownames(colData(x)),
                as.data.frame(colData(x)[, extraColAnnots, drop = FALSE])
            ), by = "sample")
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
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["value"]],
                      group = .data[[groupBy]],
                      color = .data[[colorBy]])) +
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
        p0 <- p0 + scale_color_manual(values = colors)
    }
    
    if (!is.null(highlightRegions)) {
        dfhr <- data.frame(highlightRegions)
        if (is.factor(df$position)) {
            lvs <- as.numeric(levels(df$position))
            dfhr$start <- vapply(dfhr$start, function(s) {
                if (any(lvs >= s)) {
                    as.character(min(lvs[lvs >= s]))
                } else {
                    NA_character_
                }
            }, NA_character_)
            dfhr$end <- vapply(dfhr$end, function(s) {
                if (any(lvs <= s)) {
                    as.character(max(lvs[lvs <= s]))
                } else {
                    NA_character_
                }
            }, NA_character_)
            dfhr <- dfhr[rowSums(is.na(dfhr)) == 0, ]
            dfhr$start <- factor(dfhr$start, levels = levels(df$position))
            dfhr$end <- factor(dfhr$end, levels = levels(df$position))
        }
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
#' @param referenceCoordinate A numeric scalar providing the coordinate position
#'     (on the reference sequence in \code{region}) used as an "anchor" to
#'     display relative positions. If \code{NULL} (the default), absolute
#'     genomic positions are used. 
#' @param labelAccuracy A numeric scalar indicating the precision of the 
#'     positions along the genomic axis. Will be passed to 
#'     \code{\link[scales]{label_number}}. If \code{NULL} (default), a suitable 
#'     value will be derived from \code{region}.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.createBaseplotReads <- function(df,
                                 region,
                                 trackTitle,
                                 legendTitle,
                                 showLegend,
                                 highlightRegions,
                                 facetBy, 
                                 referenceCoordinate,
                                 labelAccuracy) {
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["read"]],
                      fill = .data[["value"]])) +
        scale_fill_viridis_c(begin = 0, end = 1, option = "cividis",
                             direction = -1, na.value = "beige") +
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
             y = "Reads",
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
        p0 <- p0 +
            facet_wrap(~ .data[[facetBy]], ncol = 1, scales = "free_y")
    }
    if (is.factor(df$position)) {
        p0 <- p0 + theme(axis.text.x = element_blank())
    } else {
        p0 <- .addCoordAxisFormat(p0 = p0, region = region,
                                  labelAccuracy = labelAccuracy)
    }

    if (!is.null(highlightRegions)) {
        dfhr <- data.frame(highlightRegions)
        if (is.factor(df$position)) {
            lvs <- as.numeric(levels(df$position))
            dfhr$start <- vapply(dfhr$start, function(s) {
                if (any(lvs >= s)) {
                    as.character(min(lvs[lvs >= s]))
                } else {
                    NA_character_
                }
            }, NA_character_)
            dfhr$end <- vapply(dfhr$end, function(s) {
                if (any(lvs <= s)) {
                    as.character(max(lvs[lvs <= s]))
                } else {
                    NA_character_
                }
            }, NA_character_)
            dfhr <- dfhr[rowSums(is.na(dfhr)) == 0, ]
            dfhr$start <- factor(dfhr$start, levels = levels(df$position))
            dfhr$end <- factor(dfhr$end, levels = levels(df$position))
        }
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

#' Per-read summarize a read-level data.frame
#'
#' @param dfReads A \code{data.frame} object to summarize, typically generated
#'     by \code{\link{.preparePlotdataReads}}.
#'
#' @importFrom dplyr group_by summarise across
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.summarizePlotdataPerRead <- function(dfReads, groupVars = "sample") {
    dfReads |>
        group_by(.data[["read"]]) |>
        summarise(
            start = ifelse(is.factor(.data[["position"]]),
                           levels(.data[["position"]])[1],
                           min(.data[["position"]])),
            end = ifelse(is.factor(.data[["position"]]),
                         levels(.data[["position"]])[nlevels(.data[["position"]])],
                         max(.data[["position"]])),
            across(groupVars, unique),
            .groups = "drop")
}

#' Return ordered read identifiers
#'
#' @description
#' Returns ordered read identifiers (\code{colnames(x)} such that they follow
#' \code{hclust(as.dist(sqrt(2 - 2 * cor(X))))$order}, where \code{X} is
#' \code{assay(x, assayName)} with zero values set to \code{NA} and overaged over
#' windows of \code{windowWidth} nucleotides.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param windowWidth A numeric scalar giving the window width for which read-level
#'     data will be averaged. This should help to reduce the noise and
#'     allows to compare reads without any common modification calls, such
#'     as plus- and minus-strand reads with 6mA calls.
#'
#' @importFrom BiocGenerics colnames start
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor as.dist hclust
#' @importFrom SparseArray colMeans
#'
#' @noRd
#' @keywords internal
.orderReads <- function(x,
                        assayName,
                        windowWidth = 25) {
    # extract and flatten assay matrix
    X <- as.matrix(assay(x, assayName))
    # group positions into bins of windowWidth
    bin <- findInterval(
        x = start(x),
        vec = seq(from = min(start(x)),
                  to = ceiling(max(end(x)) / windowWidth) * windowWidth + 1,
                  by = windowWidth),
        rightmost.closed = TRUE, left.open = FALSE)
    iByBin <- split(seq.int(nrow(X)), bin)
    XX <- do.call(rbind, lapply(iByBin, function(i) {
        colMeans(X[i, , drop = FALSE], na.rm = TRUE)
    }))
    # calculate distances between reads
    D <- as.dist(sqrt(2 - 2 * cor(XX, method = "pearson",
                                  use = "pairwise.complete")))
    D[is.na(D)] <- 1.0
    # cluster reads and return order
    cl <- hclust(D, method = "ward.D2")
    return(colnames(X)[cl$order])
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
