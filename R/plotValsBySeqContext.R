#' Plot assay values stratified by sequence context
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with footprinting data. Rows should correspond to positions
#'     and columns to samples.
#' @param seqContextColumn Character scalar indicating which column of
#'     \code{rowData(se)} contains the sequence context.
#' @param assayName Character scalar indicating from which assay the values
#'     to plot will be extracted. The assay can be a read-level or summary-level
#'     assay. In both cases, the function will calculate the row means for
#'     each position across all (reads and) samples in \code{se}.
#' @param plotType Character scalar indicating what type of plot to create.
#'     Should be one of \code{"violin"}, \code{"bar"} and \code{"errorbar"}.
#' @param topN,bottomN Numeric scalars determining the number of sequence
#'     contexts to include in the plot. The \code{topN} contexts with the
#'     highest average values, and the \code{bottomN} contexts with the lowest
#'     average values will be displayed.
#' @param flipCoord Logical scalar indicating whether the plot axes should be
#'     flipped (to display sequence contexts on the y-axis and values on the
#'     x-axis) or not.
#' @param yAxisLabel Character scalar providing the label to display on the
#'     y-axis (corresponding to the values).
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
#'                                                 package = "footprintR"))
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1", modbase = "a")
#' se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = ref)
#' se <- filterPositions(se, filters = "sequenceContext", sequenceContext = "NAN")
#' plotValsBySeqContext(se, assayName = "mod_prob", plotType = "violin",
#'                      topN = 3, bottomN = 3)
#' plotValsBySeqContext(se, assayName = "mod_prob", plotType = "bar",
#'                      topN = 3, bottomN = 3)
#'
#' @importFrom dplyr group_by mutate ungroup arrange desc select distinct
#'     bind_rows slice_max slice_min filter summarize
#' @importFrom SummarizedExperiment rowData assay colnames assayNames
#' @importFrom ggplot2 ggplot aes geom_violin theme_bw geom_col geom_linerange
#'     geom_errorbar coord_flip
#' @importFrom SparseArray rowSums is_nonna
#'
plotValsBySeqContext <- function(se,
                                 seqContextColumn = "sequenceContext",
                                 assayName = "mod_prob",
                                 plotType = "violin",
                                 topN = 10,
                                 bottomN = 10,
                                 flipCoord = TRUE,
                                 yAxisLabel = "Modification probability") {
    # check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = seqContextColumn, type = "character",
                  validValues = colnames(rowData(se)))
    .assertScalar(x = assayName, type = "character",
                  validValues = assayNames(se))
    .assertScalar(x = plotType, type = "character",
                  validValues = c("violin", "bar", "errorbar"))
    .assertScalar(x = topN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = bottomN, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = flipCoord, type = "logical")
    .assertScalar(x = yAxisLabel, type = "character")

    # calculate row means
    a <- as.matrix(assay(se, assayName))
    df <- data.frame(seqContext = rowData(se)[[seqContextColumn]],
                     vals = rowSums(a, na.rm = TRUE) / rowSums(is_nonna(a)))

    # calculate mean/sd for each context and select top/bottom ones to include
    dfsum <- df |>
        group_by(seqContext) |>
        summarize(valsMean = mean(vals),
                  valsSd = sd(vals),
                  .groups = "drop")
    dfsum <- bind_rows(
        dfsum |> slice_max(valsMean, n = topN, with_ties = FALSE),
        dfsum |> slice_min(valsMean, n = bottomN, with_ties = FALSE)
    ) |>
        distinct() |>
        arrange(desc(valsMean))

    # determine order to plot contexts
    if (flipCoord) {
        levs <- rev(unique(dfsum$seqContext))
    } else {
        levs <- unique(dfsum$seqContext)
    }

    # plot
    if (plotType == "violin") {
        df <- df |>
            dplyr::filter(seqContext %in% dfsum$seqContext) |>
            mutate(seqContext = factor(seqContext, levels = levs))
        gg <- ggplot(df, aes(x = seqContext, y = vals))  +
            geom_violin()
    } else if (plotType == "bar") {
        dfsum <- dfsum |>
            mutate(seqContext = factor(seqContext, levels = levs))
        gg <- ggplot(dfsum,
               aes(x = seqContext, y = valsMean)) +
            geom_col()
    } else if (plotType == "errorbar") {
        dfsum <- dfsum |>
            mutate(seqContext = factor(seqContext, levels = levs))
        gg <- ggplot(dfsum,
               aes(x = seqContext, y = valsMean)) +
            geom_col() +
            geom_linerange(
                aes(ymin = valsMean, ymax = valsMean + valsSd)
            ) +
            geom_errorbar(
                aes(ymin = valsMean + valsSd, ymax = valsMean + valsSd),
                width = 0.2
            )
    }
    gg <- gg + theme_bw() +
        labs(x = "Sequence context",
             y = paste0(yAxisLabel, ifelse(plotType == "violin", "",
                                           ifelse(plotType == "bar", " (mean)",
                                                  " (mean + sd)"))))
    if (flipCoord) {
        gg <- gg + coord_flip()
    } else {
        gg <- gg +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    }
    gg
}
