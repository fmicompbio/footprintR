#' Plot single-molecule footprinting data for a single genomic region.
#'
#' @description
#' This function will visualize collapsed single-molecule footprinting data
#' (reads combined per genomic position), such as data imported using
#' \code{\link{readBedMethyl}}.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with collapsed single-molecule footprinting data (positions in rows and
#'     samples in columns, with assays \code{Nmod} and \code{Nvalid}).
#'     Typically, this is obtained by a call to \code{\link{readBedMethyl}}.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region that specifies the genomic region to visualize. Alternatively,
#'     the region can be specified as a character scalar (e.g. "chr1:12000-13000")
#'     that can be coerced into a \code{GRanges} object. If \code{NULL} (the
#'     default), all the data on the first sequence in \code{se} will be
#'     visualized.
#' @param min.coverage A numeric scalar giving the minimum coverage. sites that
#'     have \code{Nvalid} values below this will not be included in the
#'     visualization.
#' @param k.smooth A numeric scalar giving the number of neighboring positions
#'     over which to calculate a running median and show as smooth line(s).
#'     if \code{k = 0} (the default), no smoothing will be performed and the
#'     smoothed line is not shown.
#' @param sequence.context A character vector with sequence context(s)
#'     to plot. Only positions that match one of the provided sequence
#'     contexts will be included in the plot. Sequence contexts can be provided
#'     using IUPAC redundancy codes. The sequence contexts of modified bases are
#'     obtained from \code{rowData(se)$sequence.context} and thus requires that
#'     \code{se} contains the appropriate information, for example by setting
#'     the \code{sequence.context} and \code{sequence.reference} arguments of
#'     \code{\link{readBedMethyl}} when it was generated.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object
#'     with genomic positions in rows and samples (the unique names of
#'     \code{fnames}) in the columns.
#'
#' @author Michael Stadler
#'
#' @examples
#' bmfiles <- system.file("extdata",
#'                        c("modkit_pileup_1.bed.gz", "modkit_pileup_2.bed.gz"),
#'                        package = "footprintR")
#' reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
#'
#' se <- readBedMethyl(bmfiles, sequence.context = 3, sequence.reference = reffile)
#'
#' plotRegion(se, region = "chr1:6940000-6955000", sequence.context = "GCH")
#' plotRegion(se, region = "chr1:6940000-6955000", sequence.context = "HCG")
#' plotRegion(se, region = "chr1:6948000-6952000",
#'            min.coverage = 10, k.smooth = 7, sequence.context = "GCH")
#'
#' @seealso \code{\link{readBedMethyl}} for reading footprinting data.
#'
#' @importFrom BiocGenerics start
#' @importFrom SummarizedExperiment assay assayNames rowData
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom stats runmed
#' @importFrom dplyr filter mutate arrange group_by ungroup
#' @importFrom Biostrings vcountPattern
#' @import ggplot2
#' @importFrom rlang .data
#'
#' @export
plotRegion <- function(se, region = NULL, min.coverage = 0,
                       k.smooth = 0, sequence.context = NULL) {
    # digest arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (!all(c("Nmod", "Nvalid") %in% assayNames(se))) {
        stop("`se` needs to have assays 'Nmod' and 'Nvalid'")
    }
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    }
    .assertScalar(x = region, type = "GRanges", allowNULL = TRUE)
    if (is.null(region)) {
        region <- GRanges(seqnames = seqlevels(se)[1],
                          ranges = IRanges(start = 1, end = .Machine$integer.max))
    }
    .assertScalar(x = min.coverage, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = k.smooth, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = sequence.context, type = "character", allowNULL = TRUE)

    # subset `se`
    se <- subsetByOverlaps(x = se, ranges = region)
    if (!is.null(sequence.context)) {
        nmatch <- Reduce("+", lapply(sequence.context, function(pat) {
            vcountPattern(pat, rowData(se)$sequence.context, fixed = FALSE)
        }))
        se <- se[nmatch > 0]
    }

    # extract plot data
    df <- data.frame(
        position = start(se),
        sample = rep(colnames(se), each = nrow(se)),
        coverage = as.vector(assay(se, "Nvalid")),
        fraction_modified = as.vector(assay(se, "Nmod") /
                                          assay(se, "Nvalid")))
    if (min.coverage > 0) {
        df <- df |>
            dplyr::filter(.data[["coverage"]] >= min.coverage)
    }
    if (k.smooth > 0) {
        df <- df |>
            dplyr::group_by(sample) |>
            dplyr::arrange(.data[["position"]]) |>
            dplyr::mutate(fraction_modified_smooth = stats::runmed(
                x = .data[["fraction_modified"]],
                k = k.smooth, endrule = "constant")) |>
            dplyr::ungroup()
    }

    # create plot
    p <- ggplot(df, aes(.data[["position"]], .data[["fraction_modified"]],
                        colour = .data[["sample"]])) +
        geom_point(alpha = ifelse(k.smooth > 0, 0.2, 1)) +
        labs(x = paste0("Position on ", seqlevels(region)[1]),
             y = "Fraction modified",
             colour = "Sample") +
        theme_bw() +
        theme(legend.position = "bottom")
    if (k.smooth > 0) {
        p <- p + geom_line(
            inherit.aes = FALSE,
            mapping = aes(.data[["position"]],
                          .data[["fraction_modified_smooth"]],
                          colour = .data[["sample"]]))
    }
    if (!is.null(sequence.context)) {
        p <- p + labs(caption = paste0("Sequence contexts: ",
                                       paste(sequence.context, collapse = ", ")))
    }

    # return
    return(p)
}
