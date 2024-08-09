#' Plot single-molecule footprinting data for a single genomic region.
#'
#' @description
#' This function will visualize read-level or collapsed single-molecule
#' footprinting data, such as data imported using \code{\link{readModkitExtract}}
#' or \code{\link{readBedMethyl}}.
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
#' @param tracks.reads A named list where the names correspond to assay names
#'     of read-level assays in \code{se} and the values are character vectors
#'     with the plot types to make for each assay. Currently supported plot
#'     types are:
#'     \itemize{
#'         \item \code{"Lollipop"}: Lollipop plot (filled circles with the
#'             color representing the values in the assay).
#'         \item \code{"Heatmap"}: Heatmap plot (tiles with the color
#'             represeting the values in the assay).
#'     }
#'     If \code{NULL}, do not plot any read-level tracks.
#' @param tracks.summary A named list where the names correspond to assay names
#'     of summarized data in \code{se} and the values are character vectors with
#'     plot types to make for each assay. Currently supported plot types are:
#'     \itemize{
#'         \item \code{Point}: A point plot displaying values in the assay.
#'         \item \code{Smooth}: A smoothed line plot displaying values in the
#'             assay.
#'         \item \code{PointSmooth}: A point and smoothed line plot displaying
#'             values in the assay.
#'     }
#'     If \code{NULL}, do not plot any summary data tracks.
#'     A special case is the track name \code{"FracMod"}: If \code{se} does not
#'     contain an assay of that name, but \code{"Nmod"} and \code{"Nvalid"}
#'     assays are available, \code{"FracMod"} will be calculated from
#'     \code{assay(se, "Nmod") / assay(se, "Nvalid")}.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the x-axis will be
#'     shown in the space of modified bases and contain only the positions at
#'     which there are modified bases in the data without any gaps between them.
#'     If \code{FALSE}, the x-axis will show the genomic coordinate on which
#'     the modified bases are typically irregularly spaced.
#' @param sequence.context A character vector with sequence context(s)
#'     to plot. Only positions that match one of the provided sequence
#'     contexts will be included in the plot. Sequence contexts can be provided
#'     using IUPAC redundancy codes. The sequence contexts of modified bases are
#'     obtained from \code{rowData(se)$sequence.context} and thus requires that
#'     \code{se} contains the appropriate information, for example by setting
#'     the \code{sequence.context} and \code{sequence.reference} arguments of
#'     \code{\link{readBedMethyl}} when it was generated, or by adding it using
#'     \code{\link{seqContext}}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object with tracks selected by
#'     \code{tracks.reads} and \code{tracks.summary}.
#'
#' @author Charlotte Soneson, Michael Stadler
#'
#' @examples
#' # summarized data (5mC)
#' bmfiles <- system.file("extdata",
#'                        c("modkit_pileup_1.bed.gz", "modkit_pileup_2.bed.gz"),
#'                        package = "footprintR")
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' seA <- readBedMethyl(bmfiles, sequence.context = 3, sequence.reference = reffile)
#'
#' plotRegion(seA, region = "chr1:6940000-6955000", sequence.context = "GCH")
#' plotRegion(seA, region = "chr1:6940000-6955000", sequence.context = "HCG")
#'
#' plotRegion(seA, region = "chr1:6940000-6955000",
#'            tracks.summary = list(Nvalid = "Smooth"))
#'
#' # read-level data (6mA)
#' extractfiles <- system.file("extdata",
#'                             c("modkit_extract_rc_6mA_1.tsv.gz",
#'                               "modkit_extract_rc_6mA_2.tsv.gz"),
#'                             package = "footprintR")
#' seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit")
#'
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks.summary = NULL,
#'            tracks.reads = list(mod_prob = "Lollipop"))
#' plotRegion(seB, region = "chr1:6935800-6935900",
#'            tracks.summary = NULL,
#'            tracks.reads = list(mod_prob = "Heatmap"))
#'
#' plotRegion(seB, region = "chr1:6935400-6935450",
#'            tracks.summary = NULL,
#'            tracks.reads = list(mod_prob = c("Lollipop", "Heatmap")),
#'            modbaseSpace = TRUE)
#'
#' @seealso \code{\link{readModkitExtract}} and \code{\link{readBedMethyl}} for
#'     reading read-level and summarized footprinting data.
#'
#' @importFrom BiocGenerics start
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom dplyr filter mutate arrange group_by ungroup
#' @importFrom Biostrings vcountPattern
#' @import ggplot2
#' @importFrom patchwork wrap_plots
#' @importFrom rlang .data
#'
#' @export
plotRegion <- function(se,
                       region = NULL,
                       tracks.reads = NULL,
                       tracks.summary = list(FracMod = "Point"),
                       modbaseSpace = FALSE,
                       sequence.context = NULL) {
    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    } else if (is.null(region)) {
        region <- GRanges(seqnames = seqlevels(se)[1],
                          ranges = IRanges(start = 1, end = .Machine$integer.max))
    }
    .assertScalar(x = region, type = "GRanges", allowNULL = TRUE)
    .assertVector(x = tracks.reads, type = "list", allowNULL = TRUE)
    if (!is.null(tracks.reads)) {
        .assertVector(x = names(tracks.reads), type = "character",
                      allowNULL = FALSE, validValues = assayNames(se))
    }
    .assertVector(x = tracks.summary, type = "list", allowNULL = TRUE)
    if (!is.null(tracks.summary)) {
        .assertVector(x = names(tracks.summary), type = "character",
                      allowNULL = FALSE, validValues = c("FracMod", assayNames(se)))
    }
    .assertScalar(x = modbaseSpace, type = "logical")
    .assertVector(x = sequence.context, type = "character", allowNULL = TRUE)

    # subset se
    se <- subsetByOverlaps(x = se, ranges = region)
    if (!is.null(sequence.context)) {
        if (is.null(rowData(se)$sequence.context)) {
            stop("No sequence context found in `rowData(se)$sequence.context`")
        }
        nmatch <- Reduce("+", lapply(sequence.context, function(pat) {
            vcountPattern(pat, rowData(se)$sequence.context, fixed = FALSE)
        }))
        se <- se[nmatch > 0]
    }

    ## create plots
    pL <- list()
    ## ... summary tracks
    for (aname in names(tracks.summary)) {
        if (aname == "FracMod" && !"FracMod" %in% SummarizedExperiment::assayNames(se)) {
            if (all(c("Nmod", "Nvalid") %in% assayNames(se))) {
                assay(se, "FracMod") <- SummarizedExperiment::assay(se, "Nmod") / SummarizedExperiment::assay(se, "Nvalid")
            } else {
                stop("Cannot plot 'FracMod' - need either an assay called ",
                     "'FracMod' or both 'Nmod' and 'Nvalid' assays")
            }
        }
        for (ptype in tracks.summary[[aname]]) {
            pname <- paste0(aname, "_", ptype)
            pL[[pname]] <- switch(
                ptype,
                Point = .plotSummaryPointSmooth(x = se, aname = aname,
                                                doSmooth = FALSE,
                                                modbaseSpace = modbaseSpace),
                Smooth = .plotSummaryPointSmooth(x = se, aname = aname,
                                                 doPoint = FALSE,
                                                 modbaseSpace = modbaseSpace),
                PointSmooth = .plotSummaryPointSmooth(x = se, aname = aname,
                                                      modbaseSpace = modbaseSpace,
                                                      arglistPoint = list(alpha = 0.2))
            )
        }
    }
    ## ... read-level tracks
    for (aname in names(tracks.reads)) {
        for (ptype in tracks.reads[[aname]]) {
            pname <- paste0(aname, "_", ptype)
            pL[[pname]] <- switch(
                ptype,
                Lollipop = .plotReadsLollipop(x = se, aname = aname,
                                              modbaseSpace = modbaseSpace),
                Heatmap = .plotReadsHeatmap(x = se, aname = aname,
                                            modbaseSpace = modbaseSpace)
            )
        }
    }

    ## assemble composite plot
    if (length(pL) > 1L) { # suppress x-axis labels for all but last plot
        for (i in seq.int(length(pL) - 1L)) {
            pL[[i]] <- pL[[i]] + ggplot2::labs(x = ggplot2::element_blank())
        }
    }
    p <- patchwork::wrap_plots(pL, ncol = 1)
    if (!is.null(sequence.context)) {
        p <- p + ggplot2::labs(caption = paste0("Sequence contexts: ",
                               paste(sequence.context, collapse = ", ")))
    }

    # return
    return(p)
}


## .plot* functions for plotRegion() -------------------------------------------

#' Plot an individual track: read-level data lollipop plot.
#'
#' @description
#' This function creates a single plot track for an assay with read-level
#' data and is typically called by \code{\link{plotRegion}}.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data (positions in rows and reads in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param size A numeric scalar giving the size of the points (\code{size}
#'     argument of \code{\link[ggplot2]{geom_point}}).
#' @param drawRead A logical scalar. If \code{TRUE}, draw a horizontal line
#'     segment for each read from its start to its end.
#' @param orderReads A logical scalar. If \code{TRUE}, the position of reads
#'     on the y-axis will be reordered using \code{hclust(as.dist(1-cor(X)))$order},
#'     where \code{X} is \code{assay(x, aname)} with zero values set to \code{NA}.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the x-axis will
#'     only contain the positions of modified bases instead of all position in
#'     the genome. This can be useful to remove the gaps between modified
#'     bases for visualization.
#'
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise
#' @importFrom rlang .data
#' @importFrom BiocGenerics start nrow colnames
#' @importFrom SummarizedExperiment colData assay
#'
#' @noRd
#' @keywords internal
.plotReadsLollipop <- function(x,
                               aname,
                               size = 3.0,
                               drawRead = TRUE,
                               orderReads = TRUE,
                               modbaseSpace = FALSE) {
    # prepare plot data
    df <- .preparePlotdataReads(x, aname, modbaseSpace)

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(x, aname))
    }

    # create base plot
    p <- .createBaseplotReads(df, aname, unique(seqnames(x))[1])

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(df)
        p <- p + ggplot2::geom_segment(data = dfRead, inherit.aes = FALSE,
                                       mapping = ggplot2::aes(
                                           x = .data[["start"]],
                                           y = .data[["read"]],
                                           xend = .data[["end"]]
                                       ), colour = "gray80")
    }

    # add lollipops
    p <- p + ggplot2::geom_point(shape = 21, size = size, colour = "black")

    # return plot
    return(p)
}

#' Plot an individual track: read-level data heatmap plot.
#'
#' @description
#' This function creates a single plot track for an assay with read-level
#' data and is typically called by \code{\link{plotRegion}}.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data (positions in rows and reads in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param drawRead A logical scalar. If \code{TRUE}, draw a horizontal line
#'     segment for each read from its start to its end.
#' @param orderReads A logical scalar. If \code{TRUE}, the position of reads
#'     on the y-axis will be reordered using \code{hclust(as.dist(1-cor(X)))$order},
#'     where \code{X} is \code{assay(x, aname)} with zero values set to \code{NA}.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the x-axis will
#'     only contain the positions of modified bases instead of all position in
#'     the genome. This can be useful to remove the gaps between modified
#'     bases for visualization.
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom BiocGenerics start nrow colnames
#' @importFrom SummarizedExperiment colData assay
#'
#' @noRd
#' @keywords internal
.plotReadsHeatmap <- function(x,
                              aname,
                              drawRead = TRUE,
                              orderReads = TRUE,
                              modbaseSpace = FALSE) {
    # prepare plot data
    df <- .preparePlotdataReads(x, aname, modbaseSpace)

    # order reads
    if (orderReads) {
        df$read <- factor(as.character(df$read),
                          levels = .orderReads(x, aname))
    }

    # create base plot
    p <- .createBaseplotReads(df, aname, unique(seqnames(x))[1])

    # add segments
    if (drawRead) {
        dfRead <- .summarizePlotdataPerRead(df)
        p <- p + ggplot2::geom_segment(data = dfRead, inherit.aes = FALSE,
                                       mapping = ggplot2::aes(
                                           x = .data[["start"]],
                                           y = .data[["read"]],
                                           xend = .data[["end"]]
                                       ), colour = "gray80")
    }

    # add tiles
    p <- p + ggplot2::geom_tile(colour = "gray20", width = 1, height = 1)

    # return plot
    return(p)
}

#' Plot an individual track: summary data point and/or smooth line plot.
#'
#' @description
#' This function creates a single plot track for an assay with summary-level
#' data and is typically called by \code{\link{plotRegion}}.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param doPoint A logical scalar. If \code{TRUE}, show points in the plot.
#' @param arglistPoint A list with arguments to be sent to
#'     \code{\link[ggplot2]{geom_point}}.
#' @param doSmooth A logical scalar. If \code{TRUE}, show a smooth line in the
#'     plot.
#' @param arglistSmooth A list with arguments to be sent to
#'     \code{\link[ggplot2]{geom_line}}.
#' @param spar.smooth A numeric scalar typically in (0,1] specifying the desired
#'     degree of smoothing (\code{spar} argument of \code{\link[stats]{smooth.spline}}).
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the x-axis will
#'     only contain the positions of modified bases instead of all position in
#'     the genome. This can be useful to remove the gaps between modified
#'     bases for visualization.
#'
#' @import ggplot2
#' @importFrom BiocGenerics start nrow
#' @importFrom dplyr group_by arrange mutate ungroup bind_rows
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom stats smooth.spline
#'
#' @noRd
#' @keywords internal
.plotSummaryPointSmooth <- function(x,
                                    aname,
                                    doPoint = TRUE,
                                    arglistPoint = list(),
                                    doSmooth = TRUE,
                                    arglistSmooth = list(),
                                    spar.smooth = 0.01,
                                    modbaseSpace = FALSE) {
    # prepare plot data
    df <- .preparePlotdataSummary(x = x, aname = aname,
                                  modbaseSpace = modbaseSpace)

    # create base plot
    p <- ggplot2::ggplot(
        data = df,
        mapping = aes(x = .data[["position"]],
                      y = .data[["value"]],
                      colour = .data[["sample"]])) +
        ggplot2::labs(x = paste0("Position on ", unique(seqnames(x))[1]),
                      y = aname,
                      colour = "Sample") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "right")
    if (is.numeric(df$position)) {
        p <- p + ggplot2::coord_cartesian(xlim = range(df$position))
    }

    # add points
    if (doPoint) {
        p <- p + do.call(ggplot2::geom_point, arglistPoint)
    }

    if (doSmooth) {
        # helper function to compute smooth spline for each sample
        compute_smooth <- function (data) {
            ok <- is.finite(data[["value"]])
            smooth <- stats::smooth.spline(
                x= data[["position"]][ok],
                y = data[["value"]][ok],
                keep.data = FALSE,
                spar = spar.smooth)
            data.frame(position = smooth$x,
                       value_smooth = smooth$y,
                       sample = unique(data$sample))
        }

        # apply the function to each sample
        smooth_data <- df |>
            base::split(df[["sample"]]) |>
            purrr::map(compute_smooth) |>
            dplyr::bind_rows()

        # add the smoothed line
        p <- p + ggplot2::geom_line(
            data = smooth_data, inherit.aes = FALSE,
            mapping = aes(x = .data[["position"]],
                          y = .data[["value_smooth"]],
                          colour = .data[["sample"]]))
    }

    # return the plot
    return(p)
}


## helper functions used above -------------------------------------------------

#' Create data.frame from SummarizedExperiment for read-level data
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data (positions in rows and reads in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the "position"
#'     column in the return data frame is categorical, instead of giving
#'     the numeric position in the genome.
#'
#' @importFrom BiocGenerics start colnames
#' @importFrom SummarizedExperiment colData assay
#' @importFrom SparseArray nzwhich nzvals
#'
#' @noRd
#' @keywords internal
.preparePlotdataReads <- function(x, aname, modbaseSpace = FALSE) {
    assaydat <- SummarizedExperiment::assay(x, aname)
    if (!is.null(dim(assaydat[1,1]))) {
        # `aname` columns are grouped reads -> flatten
        sample_ids <- rep(colnames(x), unlist(lapply(assaydat, ncol)))
        assaydat <- as.matrix(assaydat)
    } else {
        sample_ids <- SummarizedExperiment::colData(x)[["sample"]]
    }
    i <- SparseArray::nzwhich(assaydat, arr.ind = TRUE)
    df <- data.frame(
        position = start(x)[i[,1]],
        read = factor(colnames(assaydat)[i[,2]], levels = colnames(assaydat)),
        sample = sample_ids[i[,2]],
        value = nzvals(assaydat))
    if (modbaseSpace) {
        df$position <- factor(df$position,
                              levels = unique(sort(df$position,
                                                   decreasing = FALSE)))
    }
    return(df)
}

#' Per-read summarize a read-level data.frame
#'
#' @param df A \code{data.frame} object to summarize, typically generated by
#'     \code{\link{.preparePlotdataReads}}.
#'
#' @importFrom dplyr group_by summarise
#' @importFrom rlang .data
#'
#' @noRd
#' @keywords internal
.summarizePlotdataPerRead <- function(df) {
    df |>
        dplyr::group_by(.data[["read"]]) |>
        dplyr::summarise(
            start = ifelse(is.factor(.data[["position"]]),
                           levels(.data[["position"]])[1],
                           min(.data[["position"]])),
            end = ifelse(is.factor(.data[["position"]]),
                         levels(.data[["position"]])[nlevels(.data[["position"]])],
                         max(.data[["position"]])),
            sample = unique(.data[["sample"]]),
            .groups = "drop")
}

#' Create data.frame from SummarizedExperiment for summary-level data
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param modbaseSpace A logical scalar. If \code{TRUE}, the "position"
#'     column in the return data frame is categorical, instead of giving
#'     the numeric position in the genome.
#'
#' @importFrom BiocGenerics start colnames
#' @importFrom SummarizedExperiment assay
#'
#' @noRd
#' @keywords internal
.preparePlotdataSummary <- function(x, aname, modbaseSpace = FALSE) {
    assaydat <- SummarizedExperiment::assay(x, aname)
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

#' Return ordered read identifiers
#'
#' @description
#' Returns ordered read identifiers (\code{colnames(x)} such that they follow
#' \code{hclust(as.dist(sqrt(2 - 2 * cor(X))))$order}, where \code{X} is
#' \code{assay(x, aname)} with zero values set to \code{NA} and overaged over
#' windows of \code{window_width} nucleotides.
#'
#' @param x A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with summary-level footprinting data (positions in rows and samples in
#'     columns).
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param window_width A numeric scalar giving the window width for which read-level
#'     data will be averaged. This should help to reduce the noise and
#'     allows to compare reads without any common modification calls, such
#'     as plus- and minus-strand reads with 6mA calls.
#'
#' @importFrom BiocGenerics colnames start
#' @importFrom SummarizedExperiment assay
#' @importFrom scuttle sumCountsAcrossFeatures
#' @importFrom stats cor as.dist
#'
#' @noRd
#' @keywords internal
.orderReads <- function(x, aname, window_width = 25) {
    # extract assay matrix and set zero to NA
    X <- SummarizedExperiment::assay(x, aname)
    if (!is.null(dim(X[1,1]))) {
        # `aname` columns are grouped reads -> flatten
        X <- as.matrix(X)
    }
    Y <- X != 0
    # group positions into bins of window_width
    bin <- findInterval(x = start(x),
                        vec = seq(from = min(start(x)),
                                  to = ceiling(max(end(x)) / window_width) * window_width + 1,
                                  by = window_width),
                        rightmost.closed = TRUE, left.open = FALSE)
    XX <- scuttle::sumCountsAcrossFeatures(x = X, ids = bin)
    YY <- scuttle::sumCountsAcrossFeatures(x = Y, ids = bin)
    XX <- XX / YY
    # calculate distances between reads
    suppressWarnings(
        D <- stats::as.dist(sqrt(2 - 2 * stats::cor(XX, method = "pearson",
                                                    use = "pairwise.complete")))
    )
    D[is.na(D)] <- 1.0
    # cluster reads and return order
    cl <- stats::hclust(D, method = "ward.D2")
    return(colnames(X)[cl$order])
}

#' Return a base ggplot2 plot (no geometries yet) for read-level data
#'
#' @param df A \code{\link{data.frame}} with the plot data (typically
#'     created by \code{\link{.preparePlotdataReads}}.
#' @param aname A character or numerical scalar selecting the assay to plot.
#' @param chr A character scaler with the sequence name that is being plotted
#'     (will be used to label the x-axis).
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom BiocGenerics colnames
#' @importFrom SummarizedExperiment assay
#' @importFrom stats cor as.dist
#'
#' @noRd
#' @keywords internal
.createBaseplotReads <- function(df, aname, chr) {
    p0 <- ggplot2::ggplot(
        data = df,
        mapping = ggplot2::aes(x = .data[["position"]],
                               y = .data[["read"]],
                               fill = .data[["value"]])) +
        ggplot2::scale_fill_gradient(low = "white", high = "black",
                                     na.value = "beige", limits = c(0, 1)) +
        ggplot2::facet_wrap(~ .data[["sample"]], ncol = 1, scales = "free_y") +
        ggplot2::labs(x = ifelse(is.numeric(df$position),
                                 paste0("Position on ", chr),
                                 paste0("Modified positions in ", chr,
                                        ":", levels(df$position)[1], "-",
                                        levels(df$position)[nlevels(df$position)])),
                      y = "Read",
                      fill = aname) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "right",
                       axis.text.y = element_blank(),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    if (is.factor(df$position)) {
        p0 <- p0 + ggplot2::theme(axis.text.x = element_blank())
    } else {
        p0 <- p0 + ggplot2::coord_cartesian(xlim = range(df$position))
    }

    return(p0)
}
