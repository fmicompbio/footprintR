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
#' @param facetBy Character scalar indicating how to facet values in the plot.
#'     Must be either \code{NULL} (in which case no facetting is done, and
#'     values are aggregated across all samples in \code{se}) or \code{"sample"}
#'     (in which case values are aggregated within each sample, and plots are
#'     facetted accordingly).
#' @param plotType Character scalar indicating what type of plot to create.
#'     Should be one of \code{"violin"}, \code{"bar"} and \code{"errorbar"}.
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
#'     bind_rows slice_max slice_min filter summarize
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
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = facetBy, type = "character", validValues = "sample",
                  allowNULL = TRUE)
    .assertScalar(x = plotType, type = "character",
                  validValues = c("violin", "bar", "errorbar"))
    .assertScalar(x = topN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = bottomN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = flipCoord, type = "logical")
    .assertScalar(x = yAxisLabel, type = "character")

    # calculate plot values
    if (aggregation == "mean" & is.null(facetBy)) {
        a <- as.matrix(assay(se, assayName))
        df <- data.frame(sample = "1",
                         seqContext = rowData(se)[[seqContextColumn]],
                         vals = rowSums(a, na.rm = TRUE) / rowSums(is_nonna(a))) |>
            # remove positions with all NA values (could happen e.g. if the input
            # se is obtained by subsetting an se with more samples)
            dplyr::filter(!is.na(.data$vals))
    } else if (aggregation == "none" & is.null(facetBy)) {
        a <- as.matrix(assay(se, assayName))
        nna <- nnawhich(a, arr.ind = TRUE)
        nnav <- nnavals(a)
        df <- data.frame(sample = "1",
                         seqContext = rowData(se)[[seqContextColumn]][nna[, 1]],
                         vals = nnav)
    } else if (aggregation == "mean" & facetBy == "sample") {
        df <- do.call(bind_rows, lapply(colnames(se), function(cn) {
            a <- assay(se, assayName)[[cn]]
            data.frame(sample = cn,
                       seqContext = rowData(se)[[seqContextColumn]],
                       vals = rowSums(a, na.rm = TRUE) / rowSums(is_nonna(a)))
        })) |>
            dplyr::filter(!is.na(.data$vals))
    } else if (aggregation == "none" & facetBy == "sample") {
        df <- do.call(bind_rows, lapply(colnames(se), function(cn) {
            a <- assay(se, assayName)[[cn]]
            nna <- nnawhich(a, arr.ind = TRUE)
            nnav <- nnavals(a)
            data.frame(sample = cn,
                       seqContext = rowData(se)[[seqContextColumn]][nna[, 1]],
                       vals = nnav)
        }))
    }
    df$flipCoord <- flipCoord

    # calculate mean/sd for each context and select top/bottom ones to include
    dfsum <- df |>
        group_by(.data$seqContext, .data$sample, .data$flipCoord) |>
        summarize(valsMean = mean(.data$vals),
                  valsSd = sd(.data$vals),
                  .groups = "drop")
    dfsum <- bind_rows(
        dfsum |> group_by(sample) |>
            slice_max(.data$valsMean, n = topN, with_ties = FALSE),
        dfsum |> group_by(sample) |>
            slice_min(.data$valsMean, n = bottomN, with_ties = FALSE)
    ) |>
        ungroup() |>
        distinct()

    # plot
    if (plotType == "violin") {
        df <- df |>
            dplyr::filter(.data$seqContext %in% dfsum$seqContext) |>
            mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
                .data$flipCoord, .data$vals, -.data$vals),
                within = sample, fun = mean))
        gg <- ggplot(df, aes(x = .data$seqContext, y = .data$vals))  +
            geom_violin(scale = "width")
    } else if (plotType %in% c("bar", "errorbar")) {
        dfsum <- dfsum |>
            mutate(seqContext = reorder_within(.data$seqContext, by = ifelse(
                .data$flipCoord, .data$valsMean, -.data$valsMean),
                within = sample, fun = mean))
        gg <- ggplot(dfsum,
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
