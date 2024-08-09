#' Get read-level data for a region from a modBAM file.
#'
#' @description
#' This function returns a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#' object with read-level footprinting data extracted for a specific \code{region}
#' from a \code{bamfile} in modBAM format.
#' \code{getReadDataForRegion()} is a convenience wrapper around three other
#' functions that could also be called independently from one another:
#' \describe{
#'     \item{\code{\link{modkitExtract}}}{runs \code{modkit} and extracts the
#'         read-level base modification information from a modBAM file.}
#'     \item{\code{\link{readModkitExtract}}}{reads the resulting output from
#'         \code{modkit} into a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'         object.}
#'     \item{\code{\link{filterReadData}}}{calculates quality measures
#'         and optionally filters out data from low quality reads.}
#' }
#'
#' @param bamfile Character scalar with the name of a modBAM file (a BAM
#'     file that stores modified-base data).
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region for which data should be extracted from \code{bamfile}.
#'     Alternatively, the region can be specified as a character scalar (e.g.
#'     "chr1:1200-1300") that can be coerced into a \code{GRanges} object.
#' @param arglist.modkitExtract A named list with further argument to pass to
#'     \code{\link{modkitExtract}} (must not contain \code{bamfile} or
#'     \code{regions}, which are set automatically).
#' @param arglist.readModkitExtract A named list with further argument to pass to
#'     \code{\link{readModkitExtract}} (must not contain \code{fnames}, which
#'     is set automatically).
#' @param arglist.filterReadData A named list with further argument to pass to
#'     \code{\link{filterReadData}}. Alternatively, it can be set to \code{NULL},
#'     in which case no filtering will be performed.
#' @param assay.type A string scalar specifying the assay name to store
#'     modification probabilities.
#' @param verbose If \code{TRUE}, report on progress.
#'
#' @return A \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object
#'     with positions of modified bases in rows and reads in the columns.
#'
#' @author Michael Stadler
#'
#' @examples
#' \dontrun{
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
#' se <- getReadDataForRegion(modbamfile,
#'                            region = "chr1:6940000-6955000",
#'                            arglist.readModkitExtract = list(modbase = "a"))
#' se
#' }
#'
#' @seealso \code{\link[SummarizedExperiment]{SummarizedExperiment}} for the
#'     returned object type, \code{\link{modkitExtract}} for calling
#'     \code{modkit}, \code{\link{readModkitExtract}} for the function
#'     used to read the \code{modkit} output and \code{\link{filterReadData}}
#'     for filtering read-level footprinting data.
#'
#'
#' @export
getReadDataForRegion <- function(bamfile,
                                 region,
                                 arglist.modkitExtract = list(modkit_bin = "modkit"),
                                 arglist.readModkitExtract = list(),
                                 arglist.filterReadData = NULL,
                                 assay.type = "mod_prob",
                                 verbose = FALSE) {
    # digest arguments
    .assertScalar(x = bamfile, type = "character")
    if (!file.exists(bamfile)) {
        stop("bamfile (", bamfile, ") does not exist.")
    }
    if (is.character(region) && length(region) == 1L) {
        region <- as(region, "GRanges")
    }
    .assertScalar(x = region, type = "GRanges")
    .assertVector(x = arglist.modkitExtract, type = "list")
    if (length(arglist.modkitExtract) > 0L && is.null(names(arglist.modkitExtract))) {
        stop("`arglist.modkitExtract` must be an empty or a named list.")
    }
    if (any(c("bamfile", "regions") %in% names(arglist.modkitExtract))) {
        stop("`arglist.modkitExtract` must not contain 'bamfile' or 'regions'",
             " - these are set automatically.")
    }
    .assertVector(x = arglist.readModkitExtract, type = "list")
    if (length(arglist.readModkitExtract) > 0L && is.null(names(arglist.readModkitExtract))) {
        stop("`arglist.readModkitExtract` must be an empty or a named list.")
    }
    if (any(c("fnames") %in% names(arglist.readModkitExtract))) {
        stop("`arglist.readModkitExtract` must not contain 'fnames'",
             " - this is set automatically.")
    }
    .assertVector(x = arglist.filterReadData, type = "list", allowNULL = TRUE)
    if (length(arglist.filterReadData) > 0L && is.null(names(arglist.filterReadData))) {
        stop("`arglist.filterReadData` must be NULL or a named list.")
    }
    .assertScalar(x = assay.type, type = "character")
    .assertScalar(x = verbose, type = "logical")

    # extract data from `bamfile`
    if (verbose) {
        message("extracting modification data in ", as.character(region))
    }
    if (!"out_read_calls" %in% names(arglist.modkitExtract)) {
        arglist.modkitExtract[["out_read_calls"]] <- tempfile(
            pattern = "modkit_extract_", fileext = ".tsv")
        on.exit(expr = unlink(arglist.modkitExtract[["out_read_calls"]]),
                add = FALSE)
    }
    if (!"verbose" %in% names(arglist.modkitExtract)) {
        arglist.modkitExtract[["verbose"]] <- verbose
    }
    arglist.modkitExtract[["bamfile"]] <- bamfile
    arglist.modkitExtract[["regions"]] <- region
    outs <- do.call(modkitExtract, arglist.modkitExtract)

    # read data into RangedSummarizedExperiment
    if (verbose) {
        message("importing data")
    }
    if (!"verbose" %in% names(arglist.readModkitExtract)) {
        arglist.readModkitExtract[["verbose"]] <- verbose
    }
    arglist.readModkitExtract[["fnames"]] <- outs[["read-calls"]]
    se <- do.call(readModkitExtract, arglist.readModkitExtract)

    # filter data
    if (!is.null(arglist.filterReadData)) {
        if (verbose) {
            message("filtering data")
        }
        # TODO
    }

    # return
    return(se)
}
