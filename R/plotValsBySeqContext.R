#' Plot assay values stratified by sequence context
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data. Rows should correspond to positions
#'     and columns to samples.
#' @param seqContextColumn Character scalar indicating which column of
#'     \code{rowData(se)} contains the sequence context.
#' @param assayName Character scalar indicating from which read-level assay the
#'     values to plot should be extracted.
#' @param aggregation Character scalar indicating if/how to aggregate values
#'     across (samples and) reads before plotting. Must be one of "none" (no
#'     aggregation, individual values from the assay are plotted) and "mean"
#'     (calculate average value for each position before plotting).
#' @param selectContextsBy Character scalar indicating how to select and order
#'     the contexts to show in the plot. Should be one of \code{"sample"} (in
#'     which case the selection of sequence contexts will be done independently
#'     for each sample, and only the ones selected for a given sample will
#'     be displayed in the plot panel for that sample; in this case
#'     \code{facetBy} must be \code{"sample"}), \code{"sample_union"} (
#'     in which case the selection of sequence contexts will be done for each
#'     sample, and the union of all selected contexts will be shown in the
#'     plot) or \code{"overall"} (in which case the selection of sequence
#'     contexts will be made without considering the sample information). Note
#'     that for \code{"sample_union"}, the plots may contain more than
#'     \code{topN + bottomN} contexts. For \code{"sample_union"} and
#'     \code{"overall"}, the order of the contexts in the plots will be
#'     determined by the average value across samples).
#' @param facetBy Character scalar indicating how to facet values in the plot.
#'     Must be either \code{NULL} (in which case no facetting is done, and
#'     values are aggregated across all samples in \code{se}) or \code{"sample"}
#'     (in which case values are aggregated within each sample, and plots are
#'     facetted accordingly).
#' @param plotType Character scalar indicating what type of plot to create.
#'     Should be one of \code{"violin"} (note that the violins will be plotted
#'     with \code{scale="width"}), \code{"bar"} and \code{"errorbar"}.
#' @param topN,bottomN Numeric scalars determining the number of sequence
#'     contexts to include in the plot. The \code{topN} contexts with the
#'     highest average values and the \code{bottomN} contexts with the lowest
#'     average values will be displayed. If \code{facetBy = "sample"}, the
#'     selection of contexts to show is performed separately within each sample.
#' @param flipCoord Logical scalar indicating whether the plot axes should be
#'     flipped (to display sequence contexts on the y-axis and values on the
#'     x-axis) or not.
#' @param yAxisLabel Character scalar providing the label to display on the
#'     value axis.
#'
#' @returns
#' A \code{ggplot} object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' modbamfile <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                           package = "footprintR")
#' ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
#'                                                 package = "footprintR"))
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1", modbase = "a",
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = ref)
#' se <- filterPositions(se, filters = "sequenceContext", sequenceContext = "NAN")
#' plotValsBySeqContext(se, assayName = "mod_prob", plotType = "violin",
#'                      topN = 3, bottomN = 3)
#' plotValsBySeqContext(se, assayName = "mod_prob", plotType = "bar",
#'                      topN = 3, bottomN = 3)
#' plotValsBySeqContext(se, assayName = "mod_prob", plotType = "bar",
#'                      topN = 3, bottomN = 3, facetBy = "sample")
#'
#' @importFrom dplyr group_by mutate ungroup arrange desc select distinct
#'     bind_rows slice_max slice_min filter summarize inner_join
#' @importFrom SummarizedExperiment rowData assay colnames assayNames
#' @importFrom ggplot2 ggplot aes geom_violin theme_bw geom_col geom_linerange
#'     geom_errorbar coord_flip
#' @importFrom SparseArray rowSums is_nonna
#' @importFrom tidytext reorder_within scale_x_reordered
#' @importFrom stats sd
#' @importFrom rlang .data
#'
plotValsBySeqContext <- function(se,
                                 seqContextColumn = "sequenceContext",
                                 assayName = "mod_prob",
                                 aggregation = "mean",
                                 selectContextsBy = "sample_union",
                                 facetBy = NULL,
                                 plotType = "violin",
                                 topN = 10,
                                 bottomN = 10,
                                 flipCoord = TRUE,
                                 yAxisLabel = "Modification probability") {
    # check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = seqContextColumn, type = "character",
                  validValues = colnames(rowData(se)))
    .assertScalar(x = aggregation, type = "character",
                  validValues = c("none", "mean"))
    .assertScalar(x = facetBy, type = "character", validValues = "sample",
                  allowNULL = TRUE)
    if (is.null(facetBy)) {
        .assertScalar(x = selectContextsBy, type = "character",
                      validValues = c("sample_union", "overall"))
    } else {
        .assertScalar(x = selectContextsBy, type = "character",
                      validValues = c("sample", "sample_union", "overall"))
    }
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = plotType, type = "character",
                  validValues = c("violin", "bar", "errorbar"))
    .assertScalar(x = topN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = bottomN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = flipCoord, type = "logical")
    .assertScalar(x = yAxisLabel, type = "character")

    # calculate plot/selection values
    # ... global ones (required if either facetBy = NULL or selectContextsBy = "overall")
    if (is.null(facetBy) || selectContextsBy == "overall") {
        if (aggregation == "mean") {
            a <- as.matrix(assay(se, assayName))
            dfGlobal <- data.frame(sample = "overall",
                                   seqContext = rowData(se)[[seqContextColumn]],
                                   vals = rowSums(a, na.rm = TRUE) / rowSums(is_nonna(a))) |>
                # remove positions with all NA values (could happen e.g. if the input
                # se is obtained by subsetting an se with more samples)
                dplyr::filter(!is.na(.data$vals))
        } else if (aggregation == "none") {
            a <- as.matrix(assay(se, assayName))
            nna <- nnawhich(a, arr.ind = TRUE)
            nnav <- nnavals(a)
            dfGlobal <- data.frame(sample = "overall",
                                   seqContext = rowData(se)[[seqContextColumn]][nna[, 1]],
                                   vals = nnav)
        }
    }

    # ... sample-wise ones (required if either facetBy = "sample" or
    #     selectContextsBy = "sample" or "sample_union")
    if ((!is.null(facetBy) && facetBy == "sample") ||
        selectContextsBy %in% c("sample", "sample_union")) {
        if (aggregation == "mean") {
            dfSample <- do.call(bind_rows, lapply(colnames(se), function(cn) {
                a <- assay(se, assayName)[[cn]]
                data.frame(sample = cn,
                           seqContext = rowData(se)[[seqContextColumn]],
                           vals = rowSums(a, na.rm = TRUE) / rowSums(is_nonna(a)))
            })) |>
                dplyr::filter(!is.na(.data$vals))
        } else if (aggregation == "none") {
            dfSample <- do.call(bind_rows, lapply(colnames(se), function(cn) {
                a <- assay(se, assayName)[[cn]]
                nna <- nnawhich(a, arr.ind = TRUE)
                nnav <- nnavals(a)
                data.frame(sample = cn,
                           seqContext = rowData(se)[[seqContextColumn]][nna[, 1]],
                           vals = nnav)
            }))
        }
    }

    # choose data frames to use for selection/plotting
    if (selectContextsBy %in% c("sample", "sample_union")) {
        dfSel <- dfSample
    } else {
        dfSel <- dfGlobal
    }

    if (is.null(facetBy)) {
        dfPlot <- dfGlobal
    } else {
        dfPlot <- dfSample
    }

    dfPlot$flipCoord <- flipCoord

    # calculate mean/sd for each context and select top/bottom ones to include
    dfSelSum <- dfSel |>
        group_by(.data$seqContext, .data$sample) |>
        summarize(valsMean = mean(.data$vals),
                  .groups = "drop")
    dfSelSum <- bind_rows(
        dfSelSum |> group_by(sample) |>
            slice_max(.data$valsMean, n = topN, with_ties = FALSE),
        dfSelSum |> group_by(sample) |>
            slice_min(.data$valsMean, n = bottomN, with_ties = FALSE)
    ) |>
        ungroup() |>
        distinct()

    # Subset dfPlot for plotting
    # If selectContextsBy = "sample", need to keep the link between sample and
    # context. Otherwise, the union of the selected contexts should all be
    # selected for all samples in dfPlot
    if (selectContextsBy == "sample") {
        # here we know that facetBy = "sample"
        dfPlot <- dfPlot |>
            inner_join(dfSelSum, by = c("sample", "seqContext")) |>
            mutate(orderCol = .data$sample)
    } else {
        contextsToKeep <- unique(dfSelSum$seqContext)
        dfPlot <- dfPlot |>
            dplyr::filter(seqContext %in% contextsToKeep) |>
            mutate(orderCol = "overall")
    }

    # plot
    if (plotType == "violin") {
        dfPlot <- dfPlot |>
            mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
                .data$flipCoord, .data$vals, -.data$vals),
                within = orderCol, fun = mean))
        gg <- ggplot(dfPlot, aes(x = .data$seqContext, y = .data$vals))  +
            geom_violin(scale = "width")
    } else if (plotType %in% c("bar", "errorbar")) {
        dfPlot <- dfPlot |>
            group_by(.data$seqContext, .data$sample, .data$orderCol,
                     .data$flipCoord) |>
            summarize(valsMean = mean(.data$vals),
                      valsSd = sd(.data$vals),
                      .groups = "drop") |>
            mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
                .data$flipCoord, .data$valsMean, -.data$valsMean),
                within = orderCol, fun = mean))
        gg <- ggplot(dfPlot,
                     aes(x = .data$seqContext, y = .data$valsMean)) +
            geom_col()
        if (plotType == "errorbar") {
            gg <- gg +
                geom_linerange(
                    aes(ymin = .data$valsMean,
                        ymax = .data$valsMean + .data$valsSd)
                ) +
                geom_errorbar(
                    aes(ymin = .data$valsMean + .data$valsSd,
                        ymax = .data$valsMean + .data$valsSd),
                    width = 0.2
                )
        }
    }
    gg <- gg +
        theme_bw() +
        labs(x = "Sequence context",
             y = paste0(yAxisLabel, ifelse(plotType == "violin", "",
                                           ifelse(plotType == "bar", " (mean)",
                                                  " (mean + sd)"))))
    if (!is.null(facetBy)) {
        gg <- gg + facet_wrap(~ sample,
                              scales = ifelse(flipCoord, "free_y", "free_x"))
    }

    if (flipCoord) {
        gg <- gg +
            coord_flip()
    } else {
        gg <- gg +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }

    gg + scale_x_reordered()
}
