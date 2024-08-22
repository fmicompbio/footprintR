#' Summarize a read-level object to sample-level.
#'
#' @description
#' This function will take a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' object with read-level footprinting data (e.g. returned by
#' \code{\link{readModkitExtract}} or \code{\link{readModBam}}) and summarize
#' reads in each sample, for instance to generate modified and total counts at
#' each position for each sample.
#'
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with read-level footprinting data. Rows should correspond to positions
#'     and columns to samples.
#' @param assay.type A string or integer scalar specifying the assay of \code{se}
#'     containing the read-level data to be summarized. Typically, this assay
#'     contains modification probabilities.
#' @param statistics Character vector specifying the type of statistics to be
#'     computed. Currently supported values are "Nmod" (number of modification
#'     probabilities greater or equal to 0.5), "Nvalid" (number of overlapping
#'     reads), "Pmod" (average modification probability), "AvgConf" (average
#'     confidence of (non-)modification probabilities).
#' @param keep.reads A scalar logical. If \code{TRUE}, the read-level data
#'     from \code{assay.type} will be retained in an assay of the same name.
#' @param verbose If \code{TRUE}, report on progress.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with the same dimensions as \code{se} (positions in rows and samples in
#'     columns), and added assays corresponding to the requested statistics.
#'
#' @author Charlotte Soneson, Michael Stadler
#'
#' @examples
#' exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz", package = "footprintR")
#' se <- readModkitExtract(exfile, modbase = "a")
#' se
#'
#' se_summary <- addReadSummary(se, keep.reads = TRUE)
#' se_summary
#'
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}} for the
#'     returned object type, \code{\link{readModkitExtract}} for the function
#'     used to read the input files.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames
#'     assay assay<- rowRanges rowRanges<- rowData rowData<-
#' @importFrom S4Vectors endoapply
#' @importFrom SparseArray pmax nzvals nzvals<- rowSums
#'
#' @export
addReadSummary <- function(se,
                           assay.type = "mod_prob",
                           statistics = c("Nmod", "Nvalid", "FracMod"),
                           keep.reads = FALSE,
                           verbose = FALSE) {
    # digest arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    if (!"sample" %in% colnames(colData(se))) {
        stop("'se' needs to contain sample information in colData(se)$sample")
    }
    if (length(assay.type) != 1L ||
        (is.numeric(assay.type) &&
         (assay.type < 1 || assay.type > length(SummarizedExperiment::assays(se)))) ||
        (is.character(assay.type) && !assay.type %in% SummarizedExperiment::assayNames(se))) {
        stop("'assay.type' must be a string or integer scalar specifying the ",
             "assay of se containing the read-level data to be summarized.")
    }
    .assertVector(x = statistics, type = "character",
                  validValues = c("Nmod", "Nvalid", "FracMod",
                                  "Pmod", "AvgConf"))
    .assertScalar(x = keep.reads, type = "logical")
    .assertScalar(x = verbose, type = "logical")

    # add statistics that are indirectly required
    statistics_use <- union(
        statistics,
        unlist(list(Nmod = character(0),
                    Nvalid = character(0),
                    FracMod = c("Nmod", "Nvalid"),
                    Pmod = "Nvalid",
                    AvgConf = "Nvalid")[statistics],
               use.names = FALSE))

    # generate DataFrame with sample-grouped columns
    # if (verbose) {
    #     message("Grouping read-level data")
    # }
    # dfReads <- S4Vectors::make_zero_col_DFrame(nrow = nrow(se))
    # ids <- factor(colData(se)$sample)
    # iBySample <- split(seq.int(ncol(se)), ids)
    # for (s in seq_along(iBySample)) {
    #     dfReads[[levels(ids)[s]]] <- assay(se, assay.type)[, iBySample[[s]]]
    # }

    # summarize reads
    if (verbose) {
        message("Summarizing reads")
    }
    dfReads <- SummarizedExperiment::assay(se, assay.type)
    assL <- lapply(structure(statistics_use, names = statistics_use), function(statistic) {
        switch(statistic,
            Nmod = as.matrix(endoapply(dfReads, function(y) rowSums(y >= 0.5))),
            Nvalid = as.matrix(endoapply(dfReads, function(y) rowSums(y != 0))),
            # these ones will be calculated later:
            FracMod = NULL,
            Pmod = NULL,
            AvgConf = NULL
        )
    })

    # calculate statistics that require multiple inputs
    if ("FracMod" %in% statistics) {
        assL[["FracMod"]] <- assL[["Nmod"]] / assL[["Nvalid"]]
    }
    if ("Pmod" %in% statistics) {
        assL[["Pmod"]] <- as.matrix(endoapply(dfReads, rowSums)) / assL[["Nvalid"]]
    }
    if ("AvgConf" %in% statistics) {
        # confidence: max(mod_prob, 1 - mod_prob)
        assL[["AvgConf"]] <- as.matrix(
            endoapply(dfReads, function(y) {
                nzvals(y) <- pmax(nzvals(y), 1 - nzvals(y))
                rowSums(y)
            })
        ) / assL[["Nvalid"]]
    }

    # create summarized experiment
    if (verbose) {
        message("Creating SummarizedExperiment")
    }
    se_summary <- SummarizedExperiment::SummarizedExperiment(
        assays = assL[statistics],
        colData = SummarizedExperiment::colData(se)
    )
    if (!is.null(SummarizedExperiment::rowRanges(se))) {
        SummarizedExperiment::rowRanges(se_summary) <- SummarizedExperiment::rowRanges(se)
    } else {
        SummarizedExperiment::rowData(se_summary) <- SummarizedExperiment::rowData(se)
    }

    # keep read-level data
    if (keep.reads) {
        SummarizedExperiment::assay(se_summary, assay.type,
                                    withDimnames = FALSE) <- dfReads
    }

    # return
    return(se_summary)
}
