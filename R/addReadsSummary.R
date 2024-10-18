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
#' @param assay.type A character scalar specifying the assay of \code{se}
#'     containing the read-level data to be summarized. Typically, this assay
#'     contains modification probabilities.
#' @param statistics Character vector specifying the type of statistics to be
#'     computed. Currently supported values are "Nmod" (number of modification
#'     probabilities greater or equal to 0.5), "Nvalid" (number of overlapping
#'     reads), "Pmod" (average modification probability), "AvgConf" (average
#'     confidence of (non-)modification probabilities).
#' @param keep.reads A logical scalar. If \code{TRUE} (the default), the
#'     read-level data from \code{assay.type} will be retained in an assay of
#'     the same name.
#' @param replace.existing A logical scalar. If \code{TRUE} (the default), 
#'     any existing assays with the same name as the ones requested will be 
#'     overwritten. Otherwise, existing assays will be retained and the 
#'     corresponding summary statistic(s) will not be recalculated. 
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
#' se_summary <- addReadsSummary(se)
#' se_summary
#'
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}} for the
#'     returned object type, \code{\link{readModkitExtract}} for the function
#'     used to read the input files.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment assays assayNames
#'     assay assay<- rowRanges rowRanges<- rowData rowData<-
#' @importFrom S4Vectors endoapply metadata
#' @importFrom SparseArray pmax nnavals nnavals<- rowSums
#'
#' @export
addReadsSummary <- function(se,
                            assay.type = "mod_prob",
                            statistics = c("Nmod", "Nvalid", "FracMod"),
                            keep.reads = TRUE,
                            replace.existing = TRUE, 
                            verbose = FALSE) {
    # digest arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assay.type, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertVector(x = statistics, type = "character",
                  validValues = c("Nmod", "Nvalid", "FracMod",
                                  "Pmod", "AvgConf"))
    .assertScalar(x = keep.reads, type = "logical")
    .assertScalar(x = replace.existing, type = "logical")
    .assertScalar(x = verbose, type = "logical")
    
    # if replace.existing is FALSE, exclude all assays that already exist in se
    if (!replace.existing) {
        existing_assays <- intersect(statistics, SummarizedExperiment::assayNames(se))
        if (length(existing_assays) > 0) {
            warning("Assay(s) ", paste(existing_assays, ", "), 
                    " already exist and replace.existing is FALSE - will not ",
                    "recalculate these assays.")
            statistics <- setdiff(statistics, SummarizedExperiment::assayNames(se))
        }
    }

    # add statistics that are indirectly required
    statistics_use <- union(
        statistics,
        unlist(list(Nmod = character(0),
                    Nvalid = character(0),
                    FracMod = c("Nmod", "Nvalid"),
                    Pmod = "Nvalid",
                    AvgConf = "Nvalid")[statistics],
               use.names = FALSE))

    # summarize reads
    # Note: all summarized assays are dense matrices, because:
    # - Nmod, Nvalid and AvgConf are dense because rowSums is based on
    #   MatrixGenerics::rowSums, which returns a (dense) numeric vector
    # - Nvalid typically has few zeros, and the SparseArray version would be
    #   larger in memory than the dense one
    # - FracMod and Pmod have to be calculated using dense matrices and would
    #   have to be converted back to sparse objects  (no "/" method for
    #   SparseArray objects, as the result wouldn't be sparse)
    if (verbose) {
        message("Summarizing reads")
    }
    dfReads <- SummarizedExperiment::assay(se, assay.type)
    assL <- lapply(structure(statistics_use, names = statistics_use), function(statistic) {
        switch(statistic,
            Nmod = as.matrix(endoapply(dfReads, function(y) rowSums(y >= 0.5, na.rm = TRUE))),
            Nvalid = as.matrix(endoapply(dfReads, function(y) rowSums(y >= 0, na.rm = TRUE))),
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
        assL[["Pmod"]] <- as.matrix(endoapply(dfReads, rowSums, na.rm = TRUE)) / assL[["Nvalid"]]
    }
    if ("AvgConf" %in% statistics) {
        # confidence: max(mod_prob, 1 - mod_prob)
        assL[["AvgConf"]] <- as.matrix(
            endoapply(dfReads, function(y) {
                SparseArray::nnavals(y) <- pmax(SparseArray::nnavals(y),
                                                1 - SparseArray::nnavals(y))
                rowSums(y, na.rm = TRUE)
            })
        ) / assL[["Nvalid"]]
    }

    # create summarized experiment
    if (verbose) {
        message("Creating SummarizedExperiment")
    }
    tmpList <- as.list(SummarizedExperiment::assays(se))
    tmpList[statistics] <- lapply(assL[statistics], function(a) {
        rownames(a) <- rownames(se)
        a
    })
    SummarizedExperiment::assays(se) <- tmpList

    # keep read-level data
    if (!keep.reads) {
        rlAssays <- .getReadLevelAssayNames(se)
        suppressWarnings(
            # currently, assigning to assays triggers a deprecation warning
            # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
            # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
            SummarizedExperiment::assays(se) <- 
                SummarizedExperiment::assays(se)[setdiff(
                    SummarizedExperiment::assayNames(se), rlAssays)]
        )
        ## Remove read-level assays from the metadata
        S4Vectors::metadata(se)$readLevelData$assayNames <- character(0)
    }

    # return
    return(se)
}
