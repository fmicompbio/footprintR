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
#' \code{plotGenomicRegions} functions are more low-level helper functions 
#' for creating single plot tracks. These are invoked by \code{plotRegion}, 
#' and typically do not need to be directly called by the user.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level or collapsed single-molecule footprinting data (positions
#'     in rows and reads or samples in columns). 
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region. Only data from \code{se} overlapping this region will be plotted.
#'     Alternatively, the region can be specified as a character scalar (e.g.
#'     "chr1:1200-1300") that can be coerced into a \code{GRanges} object. If
#'     \code{NULL} (the default), all the data on the first sequence in
#'     \code{se} will be visualized.
#' @param tracks A list of named lists, representing the tracks to generate.
#'     Each element of the outer list defines one track, and has to contain
#'     at least list entries named 'trackData' (the name of a suitable assay
#'     in \code{se}) and 'trackType' (the type of plot), plus any additional
#'     arguments to the respective plot function. Currently supported plot
#'     types are
#'     \describe{
#'         \item{\code{"Point"}}{: A point plot displaying values in the assay.}
#'         \item{\code{"Smooth"}}{: A smoothed line plot displaying values in the
#'             assay.}
#'         \item{\code{"PointSmooth"}}{: A point and smoothed line plot displaying
#'             values in the assay.}
#'         \item{\code{"Lollipop"}}{: Lollipop plot (filled circles with the
#'             color representing the values in the assay).}
#'         \item{\code{"Heatmap"}}{: Heatmap plot (tiles with the color
#'             representing the values in the assay).}
#'         \item{\code{"GenomicRegion"}}{: Genomic annotations (e.g., transcripts, 
#'             peaks, CpG islands).}
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
#'     the \code{sequenceContext} and \code{sequenceReference} arguments of
#'     \code{\link{readBedMethyl}} when it was generated, or by adding it using
#'     \code{\link{addSeqContext}}.
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
#'                      sequenceContextWidth = 3, sequenceReference = reffile)
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
#' seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit")
#'
#' # Lollipop plot
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop")))
#' # Heatmap plots (observed only or filled)
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Heatmap")))
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
#'                               interpolate = TRUE)))
#'
#' # multiple plots
#' plotRegion(seB, region = "chr1:6935400-6935450",
#'            tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
#'                               size = 4),
#'                          list(trackData = "mod_prob", trackType = "Heatmap")),
#'            modbaseSpace = TRUE)
#'
#' @seealso \code{\link{readModBam}}, \code{\link{readModkitExtract}} and
#'     \code{\link{readBedMethyl}} for reading read-level and summarized
#'     footprinting data.
#'
#' @importFrom BiocGenerics start intersect
#' @importFrom SummarizedExperiment assay assayNames rowData nrow rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom dplyr filter mutate arrange group_by ungroup
#' @importFrom Biostrings vcountPattern
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom rlang .data
#' @importFrom cli cli_abort cli_warn
#'
#' @export
plotRegion <- function(se,
                       region = NULL,
                       tracks = list(list(trackData = "FracMod", trackType = "Point")),
                       modbaseSpace = FALSE,
                       sequenceContext = NULL) {
    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    } else if (is.null(region)) {
        region <- GRanges(seqnames = seqlevels(se)[1],
                          ranges = IRanges(start = 1, end = .Machine$integer.max))
    }
    .assertScalar(x = region, type = "GRanges", allowNULL = TRUE)
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertVector(x = tracks, type = "list", rngLen = c(1, Inf))
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

    # subset se
    se <- subsetByOverlaps(x = se, ranges = region)
    se <- .keepPositionsBySequenceContext(se = se, sequenceContext = sequenceContext)
    
    ## create plots
    pL <- vector("list", length = length(tracks))
    for (i in seq_along(tracks)) {
        tr <- tracks[[i]]
        trt <- plotRegionPlotTypes$type[match(tr$trackType, plotRegionPlotTypes$name)]
        if (trt %in% c("summary", "reads")) {
            args <- c(
                list(se = se, region = region, assayName = tr$trackData, 
                     modbaseSpace = modbaseSpace),
                tr[!names(tr) %in% c("trackData", "trackType", "se", "region",
                                     "assayName", "modbaseSpace", "doSmooth", 
                                     "doPoint")]
            )
        } else if (trt == "annotation") {
            args <- c(
                list(grl = subsetByOverlaps(tr$trackData, region), 
                     region = region), 
                     tr[!names(tr) %in% c("trackData", "trackType", "grl", 
                                          "region")]
            )
        }
        pL[[i]] <- switch(
            tr$trackType,
            Point = do.call(plotSummaryPointSmooth, c(args, list(doSmooth = FALSE))),
            Smooth = do.call(plotSummaryPointSmooth, c(args, list(doPoint = FALSE))),
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
#' @param facetBySample A logical scalar indicating whether or not to facet
#'     the plot by sample (for read-level plots). 
#' 
#' @export
#' @rdname plotRegion
#' 
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise
#' @importFrom rlang .data
#' @importFrom BiocGenerics start nrow colnames
#' @importFrom SummarizedExperiment colData assay
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
                              facetBySample = TRUE) {
    .assertVector(x = highlightRegions, type = "GenomicRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        ## subset highlightRegions - mostly to ensure to only retain 
        ## regions on the right chromosome, so that we can only focus on 
        ## the positions below
        highlightRegions <- subsetByOverlaps(highlightRegions, 
                                             region,
                                             ignore.strand = TRUE)
    }
    
    # prepare plot data
    df <- .preparePlotdataReads(se, assayName, modbaseSpace)

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(se, assayName))
    }

    # create base plot
    p <- .createBaseplotReads(df = df, assayName = assayName, region = region, 
                              trackTitle = trackTitle,
                              legendTitle = legendTitle,
                              showLegend = showLegend,
                              highlightRegions = highlightRegions,
                              facetBySample = facetBySample)

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(df)
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["read"]],
                                  xend = .data[["end"]]
                              ), colour = "gray80")
    }

    # add lollipops
    p <- p + geom_point(shape = 21, size = size,
                        stroke = stroke, colour = "black")

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
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom BiocGenerics start nrow colnames
#' @importFrom SummarizedExperiment colData assay
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
                             facetBySample = TRUE) {
    .assertVector(x = highlightRegions, type = "GenomicRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- subsetByOverlaps(highlightRegions, 
                                             region,
                                             ignore.strand = TRUE)
    }
    
    # prepare plot data
    df <- .preparePlotdataReads(se, assayName, modbaseSpace, interpolate)

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(se, assayName))
    }

    # create base plot
    p <- .createBaseplotReads(df = df, assayName = assayName, region = region, 
                              trackTitle = trackTitle,
                              legendTitle = legendTitle,
                              showLegend = showLegend,
                              highlightRegions = highlightRegions,
                              facetBySample = facetBySample)

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(df)
        p <- p + geom_segment(data = dfRead, inherit.aes = FALSE,
                              mapping = aes(
                                  x = .data[["start"]],
                                  y = .data[["read"]],
                                  xend = .data[["end"]]
                              ), colour = "gray80")
    }

    # add tiles
    p <- p + geom_tile(colour = "gray20", width = 1, height = 1,
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
#' @param spar A numeric scalar typically in (0,1] specifying the desired
#'     degree of smoothing (\code{spar} argument of \code{\link[stats]{smooth.spline}}).
#'
#' @export
#' @rdname plotRegion
#' 
#' @import ggplot2
#' @importFrom BiocGenerics start nrow
#' @importFrom dplyr group_by arrange mutate ungroup group_modify
#' @importFrom rlang .data
#' @importFrom stats smooth.spline
#'
plotSummaryPointSmooth <- function(se,
                                   region,
                                   assayName,
                                   doPoint = TRUE,
                                   arglistPoint = list(),
                                   doSmooth = TRUE,
                                   arglistSmooth = list(),
                                   spar = 0.01,
                                   modbaseSpace = FALSE,
                                   trackTitle = NULL,
                                   legendTitle = NULL,
                                   showLegend = TRUE,
                                   highlightRegions = NULL) {
    
    .assertVector(x = highlightRegions, type = "GenomicRanges",
                  allowNULL = TRUE)
    if (!is.null(highlightRegions)) {
        highlightRegions <- subsetByOverlaps(highlightRegions, 
                                             region,
                                             ignore.strand = TRUE)
    }
    
    # prepare plot data
    df <- .preparePlotdataSummary(x = se, assayName = assayName,
                                  modbaseSpace = modbaseSpace)

    # create base plot
    p <- .createBaseplotSummary(df = df, assayName = assayName,
                                region = region, 
                                trackTitle = trackTitle,
                                legendTitle = legendTitle,
                                showLegend = showLegend,
                                highlightRegions = highlightRegions)

    # add points
    if (doPoint) {
        p <- p + do.call(geom_point, arglistPoint)
    }

    if (doSmooth) {
        # helper function to compute smooth spline for each sample
        compute_smooth <- function(data) {
            ok <- is.finite(data[["value"]])
            smooth <- smooth.spline(
                x = data[["position"]][ok],
                y = data[["value"]][ok],
                keep.data = FALSE,
                spar = spar)
            data.frame(position = smooth$x,
                       value_smooth = smooth$y)
        }

        # apply the function to each sample
        smooth_data <- df |>
            group_by(sample) |>
            group_modify(~ compute_smooth(.x)) |>
            ungroup()

        # add the smoothed line
        p <- p + geom_line(
            data = smooth_data, inherit.aes = FALSE,
            mapping = aes(x = .data[["position"]],
                          y = .data[["value_smooth"]],
                          colour = .data[["sample"]]))
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
                               showLegend = TRUE) {
    # check input arguments
    .assertVector(x = grl, type = "GRangesList")
    .assertVector(x = names(grl), type = "character")
    .assertScalar(x = region, type = "GRanges")
    .assertScalar(x = colorByStrand, type = "logical")
    .assertScalar(x = displayNames, type = "logical")
    .assertScalar(x = labelSize, type = "numeric")
    .assertScalar(x = labelPosition, type = "character", 
                  validValues = c("above", "below", "inside"))
    
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
    
    # plot
    gg <- ggplot() + 
        geom_segment(data = fullRange,
                     mapping = aes(
                         x = .data[["start"]],
                         y = .data[["fpname_unique"]],
                         xend = .data[["end"]]
                     ), colour = "gray80")
    if (colorByStrand) {
        gg <- gg + 
            geom_rect(data = rangeParts, 
                      mapping = aes(
                          xmin = .data[["start"]],
                          xmax = .data[["end"]],
                          ymin = as.numeric(.data[["fpname_unique"]]) - 0.25,
                          ymax = as.numeric(.data[["fpname_unique"]]) + 0.25,
                          fill = .data[["strand"]]
                      ), colour = "gray20") + 
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
                      ), colour = "gray20", fill = "gray80")
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
    
    acc <- 10^round(log10((rng[2] - rng[1]) / rng[2]))
    gg <- gg + coord_cartesian(xlim = rng) +
        scale_x_continuous(
            expand = c(0, 0), 
            labels = label_number(
                accuracy = acc,
                scale_cut = c(0, ` Kb` = 1000, ` Mb` = 1e+06, ` Bb` = 1e+12)))
    
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
#'
#' @importFrom BiocGenerics start colnames
#' @importFrom SummarizedExperiment assay
#'
#' @noRd
#' @keywords internal
.preparePlotdataSummary <- function(x,
                                    assayName,
                                    modbaseSpace = FALSE) {
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
#'
#' @importFrom BiocGenerics start colnames
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SparseArray nnawhich nnavals
#'
#' @noRd
#' @keywords internal
.preparePlotdataReads <- function(x,
                                  assayName,
                                  modbaseSpace = FALSE,
                                  interpolate = FALSE) {
    assaydat <- assay(x, assayName)
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
    }
    return(df)
}

#' Return a base ggplot2 plot (no geometries yet) for summary-level data
#'
#' @param df A \code{\link{data.frame}} with the plot data (typically
#'     created by \code{\link{.preparePlotdataSummary}}.
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param chr A character scaler with the sequence name that is being plotted
#'     (will be used to label the x-axis).
#' @param trackTitle A character scalar or \code{NULL}, giving the title of 
#'     the track.
#' @param legendTitle A character scalar (or \code{NULL}) giving the title for
#'     the color legend.
#' @param showLegend A logical scalar indicating whether or not to show the 
#'     legend for the plot.
#' @param highlightRegions A \code{\link[GenomicRanges]{GRanges}} object 
#'     containing regions to highlight with a grey shading.
#'
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.createBaseplotSummary <- function(df,
                                   assayName,
                                   region, 
                                   trackTitle,
                                   legendTitle,
                                   showLegend,
                                   highlightRegions) {
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["value"]],
                      colour = .data[["sample"]])) +
        labs(x = paste0("Position on ", as.character(seqnames(region))),
             y = assayName,
             colour = ifelse(!is.null(legendTitle), legendTitle, "Sample"),
             title = trackTitle) +
        theme_bw() +
        theme(legend.position = ifelse(showLegend, "right", "none"))

    if (!is.null(highlightRegions)) {
        p0 <- p0 + 
            geom_rect(
                data = data.frame(highlightRegions),
                mapping = aes(xmin = start, xmax = end,
                              ymin = -Inf, ymax = Inf),
                fill = "gray90",
                inherit.aes = FALSE
            )
    }
    
    if (is.numeric(df$position)) {
        p0 <- .addCoordAxisFormat(p0 = p0, region = region)
    }

    return(p0)
}

#' Return a base ggplot2 plot (no geometries yet) for read-level data
#'
#' @param df A \code{\link{data.frame}} with the plot data (typically
#'     created by \code{\link{.preparePlotdataReads}}.
#' @param assayName A character or numerical scalar selecting the assay to plot.
#' @param chr A character scalar with the sequence name that is being plotted
#'     (will be used to label the x-axis).
#' @param trackTitle A character scalar or \code{NULL}, giving the title of 
#'     the track.
#' @param legendTitle A character scalar (or \code{NULL}) giving the title for
#'     the fill legend.
#' @param showLegend A logical scalar indicating whether or not to show the 
#'     legend for the plot.
#' @param highlightRegions A \code{\link[GenomicRanges]{GRanges}} object 
#'     containing regions to highlight with a grey shading.
#' @param facetBySample A logical scalar indicating whether or not to facet
#'     the plot by sample. 
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom BiocGenerics colnames
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor as.dist
#'
#' @noRd
#' @keywords internal
.createBaseplotReads <- function(df,
                                 assayName,
                                 region,
                                 trackTitle,
                                 legendTitle,
                                 showLegend,
                                 highlightRegions,
                                 facetBySample) {
    p0 <- ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["read"]],
                      fill = .data[["value"]])) +
        scale_fill_viridis_c(begin = 0, end = 1, option = "cividis",
                             direction = -1, na.value = "beige") +
        labs(x = ifelse(is.numeric(df$position),
                        paste0("Position on ", as.character(seqnames(region))),
                        paste0("Modified positions in ", as.character(seqnames(region)),
                               ":", levels(df$position)[1], "-",
                               levels(df$position)[nlevels(df$position)])),
             y = "Reads",
             fill = ifelse(!is.null(legendTitle), legendTitle, assayName),
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

    if (facetBySample) {
        p0 <- p0 + 
            facet_wrap(~ .data[["sample"]], ncol = 1, scales = "free_y")
    }
    if (is.factor(df$position)) {
        p0 <- p0 + theme(axis.text.x = element_blank())

    } else {
        p0 <- .addCoordAxisFormat(p0 = p0, region = region)
    }
    
    if (!is.null(highlightRegions)) {
        p0 <- p0 + 
            geom_rect(
                data = data.frame(highlightRegions),
                mapping = aes(xmin = start, xmax = end,
                              ymin = -Inf, ymax = Inf),
                fill = "gray90",
                inherit.aes = FALSE
            )
    }
    
    return(p0)
}

#' Per-read summarize a read-level data.frame
#'
#' @param dfReads A \code{data.frame} object to summarize, typically generated
#'     by \code{\link{.preparePlotdataReads}}.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.summarizePlotdataPerRead <- function(dfReads) {
    dfReads |>
        group_by(.data[["read"]]) |>
        summarise(
            start = ifelse(is.factor(.data[["position"]]),
                           levels(.data[["position"]])[1],
                           min(.data[["position"]])),
            end = ifelse(is.factor(.data[["position"]]),
                         levels(.data[["position"]])[nlevels(.data[["position"]])],
                         max(.data[["position"]])),
            sample = unique(.data[["sample"]]),
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
    bin <- findInterval(x = start(x),
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
#'
#' @import ggplot2
#' @importFrom scales label_number
#'
#' @noRd
#' @keywords internal
.addCoordAxisFormat <- function(p0, region) {
    rng <- c(start(region) - 0.5, end(region) + 0.5)
    acc <- 10^round(log10((rng[2] - rng[1]) / rng[2]))
    p0 <- p0 + coord_cartesian(xlim = rng) +
        scale_x_continuous(
            expand = c(0, 0), 
            labels = label_number(
                accuracy = acc,
                scale_cut = c(0, ` Kb` = 1000, ` Mb` = 1e+06, ` Bb` = 1e+12)))
    return(p0)
}
