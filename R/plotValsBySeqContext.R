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
#'     plot), \code{"overall"} (in which case the selection of sequence
#'     contexts will be made without considering the sample information) or
#'     \code{"sample_var"} (in which case the variance across samples will be
#'     used to select/order contexts). Note
#'     that for \code{"sample_union"}, the plots may contain more than
#'     \code{topN + bottomN} contexts. For \code{"sample_union"} and
#'     \code{"overall"}, the order of the contexts in the plots will be
#'     determined by the average value across samples).
#' @param facetBy Character scalar indicating how to facet values in the plot.
#'     Must be either \code{NULL} (in which case no facetting is done, and
#'     values are aggregated across all samples in \code{se}) or \code{"sample"}
#'     (in which case values are aggregated within each sample, and plots are
#'     facetted accordingly). Ignored if \code{plotType = "pairs"}.
#' @param fillBy Character scalar indicating how to fill the bars or violins.
#'     Must be either \code{NULL} (in which case a single value should be
#'     specified to \code{"fillColors"} and used for all bars/violins) or
#'     \code{"sample"} (in which case bars or violins will be split and filled
#'     by sample). Ignored if \code{plotType = "pairs"}.
#' @param fillColors Either a (preferably named) character vector defining the
#'     color to use for each sample if \code{fillBy = "sample"}, or a single
#'     color to use for all violins/bars.
#' @param plotType Character scalar indicating what type of plot to create.
#'     Should be one of \code{"violin"} (note that the violins will be plotted
#'     with \code{scale="width"}), \code{"bar"}, \code{"errorbar"} or
#'     code{"pairs"}.
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
#'     geom_errorbar coord_flip position_dodge scale_fill_manual
#' @importFrom SparseArray rowSums is_nonna
#' @importFrom tidytext reorder_within scale_x_reordered
#' @importFrom stats sd
#' @importFrom rlang .data
#' @importFrom cli cli_warn
#'
plotValsBySeqContext <- function(se,
                                 seqContextColumn = "sequenceContext",
                                 assayName = "mod_prob",
                                 aggregation = "mean",
                                 selectContextsBy = "sample_union",
                                 facetBy = NULL,
                                 fillBy = NULL,
                                 fillColors = "grey",
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
    .assertScalar(x = plotType, type = "character",
                  validValues = c("violin", "bar", "errorbar", "pairs"))
    if (plotType == "pairs") {
        .assertPackagesAvailable("GGally")
    }
    if (is.null(facetBy) || plotType == "pairs") {
        .assertScalar(x = selectContextsBy, type = "character",
                      validValues = c("sample_union", "overall", "sample_var"))
    } else {
        .assertScalar(x = selectContextsBy, type = "character",
                      validValues = c("sample", "sample_union", "overall",
                                      "sample_var"))
    }
    if (plotType == "pairs" && ncol(se) < 2) {
        cli_abort("{.arg se} must have at least two samples for a pairs plot")
    }
    .assertScalar(x = fillBy, type = "character", validValues = "sample",
                  allowNULL = TRUE)
    .assertVector(x = fillColors, type = "character", rngLen = c(1, Inf))
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    if (selectContextsBy %in% c("sample_var") && ncol(se) < 2) {
        cli_abort("{.arg se} must have at least two samples if {.arg selectContextsBy} is {.val {selectContextsBy}}")
    }
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
    #     selectContextsBy = "sample", "sample_union" or "sample_var")
    if ((!is.null(facetBy) && facetBy == "sample") ||
        (!is.null(fillBy) && fillBy == "sample") ||
        selectContextsBy %in% c("sample", "sample_union", "sample_var")) {
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
    if (selectContextsBy %in% c("sample", "sample_union", "sample_var")) {
        dfSel <- dfSample
    } else {
        dfSel <- dfGlobal
    }

    if (is.null(facetBy) && is.null(fillBy)) {
        dfPlot <- dfGlobal
    } else {
        dfPlot <- dfSample
    }

    dfPlot$flipCoord <- flipCoord
    dfPlot$selectContextsBy <- selectContextsBy

    # calculate mean/sd for each context and select top/bottom ones to include
    dfSelSum <- dfSel |>
        group_by(.data$seqContext, .data$sample) |>
        summarize(valsMean = mean(.data$vals),
                  .groups = "drop")
    if (selectContextsBy == "sample_var") {
        dfSelSum <- dfSelSum |>
            group_by(.data$seqContext) |>
            mutate(valsMeanVar = var(.data$valsMean))
        selCol <- "valsMeanVar"
    } else {
        selCol <- "valsMean"
    }
    dfSelSum <- bind_rows(
        dfSelSum |> group_by(sample) |>
            slice_max(.data[[selCol]], n = topN, with_ties = FALSE),
        dfSelSum |> group_by(sample) |>
            slice_min(.data[[selCol]], n = bottomN, with_ties = FALSE)
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
        if (selectContextsBy == "sample_var") {
            dfPlot <- dfPlot |>
                left_join(dfSelSum |>
                              dplyr::select(c("seqContext", "valsMeanVar")),
                          by = "seqContext")
        }
    }

    # plot
    if (plotType == "pairs") {
        dfPlot <- dfPlot |>
            group_by(.data$seqContext, .data$sample, .data$orderCol,
                     .data$flipCoord) |>
            summarize(valsMean = mean(.data$vals),
                      valsSd = sd(.data$vals),
                      .groups = "drop") |>
            dplyr::select(c("seqContext", "sample", "valsMean")) |>
            tidyr::pivot_wider(names_from = "sample", values_from = "valsMean")
        GGally::ggpairs(dfPlot, columns = setdiff(colnames(dfPlot), "seqContext")) +
            theme_bw()
    } else {
        if (plotType == "violin") {
            dfPlot <- dfPlot |>
                mutate(seqContext = reorder_within(
                    .data$seqContext,
                    by = ifelse(flipCoord, 1, -1) *
                        ifelse(.data$selectContextsBy == "sample_var",
                               .data$valsMeanVar, .data$vals),
                    within = orderCol,
                    fun = mean
                ))
            # dfPlot <- dfPlot |>
            #     mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
            #         .data$flipCoord, .data$vals, -.data$vals),
            #         within = orderCol, fun = mean))
            gg <- ggplot(dfPlot, aes(x = .data$seqContext, y = .data$vals))
            if (is.null(fillBy)) {
                gg <- gg +
                    geom_violin(scale = "width", fill = fillColors[1])
            } else {
                gg <- gg +
                    geom_violin(scale = "width", aes(fill = .data[[fillBy]]))
                if (length(fillColors) >= length(unique(dfPlot[[fillBy]]))) {
                    gg <- gg +
                        scale_fill_manual(values = fillColors)
                } else {
                    cli_warn("Not enough colors - using defaults")
                }
            }
        } else if (plotType %in% c("bar", "errorbar")) {
            dfPlot <- dfPlot |>
                group_by(.data$seqContext, .data$sample, .data$orderCol,
                         .data$flipCoord, .data$valsMeanVar) |>
                summarize(valsMean = mean(.data$vals),
                          valsSd = sd(.data$vals),
                          .groups = "drop") |>
                mutate(seqContext = reorder_within(
                    .data$seqContext,
                    by = ifelse(flipCoord, 1, -1) *
                        ifelse(.data$selectContextsBy == "sample_var",
                               .data$valsMeanVar, .data$valsMean),
                    within = orderCol,
                    fun = mean
                ))
                # mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
                #     .data$flipCoord, .data$valsMean, -.data$valsMean),
                #     within = orderCol, fun = mean))
            gg <- ggplot(dfPlot,
                         aes(x = .data$seqContext, y = .data$valsMean))
            if (is.null(fillBy)) {
                gg <- gg +
                    geom_col(fill = fillColors[1])
            } else {
                gg <- gg +
                    geom_col(aes(fill = .data[[fillBy]]),
                             position = position_dodge())
                if (length(fillColors) >= length(unique(dfPlot[[fillBy]]))) {
                    gg <- gg +
                        scale_fill_manual(values = fillColors)
                } else {
                    cli_warn("Not enough colors - using defaults")
                }
            }
            if (plotType == "errorbar") {
                if (is.null(fillBy)) {
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
                } else {
                    gg <- gg +
                        geom_linerange(
                            aes(group = .data[[fillBy]],
                                ymin = .data$valsMean,
                                ymax = .data$valsMean + .data$valsSd),
                            position = position_dodge(width = 0.9)
                        ) +
                        geom_errorbar(
                            aes(group = .data[[fillBy]],
                                ymin = .data$valsMean + .data$valsSd,
                                ymax = .data$valsMean + .data$valsSd),
                            width = 0.2,
                            position = position_dodge(width = 0.9)
                        )
                }
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
}
