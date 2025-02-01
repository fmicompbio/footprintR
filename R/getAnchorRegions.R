#' Extract data for one or more anchor regions
#'
#' Extract assay values for each read (or sample) and each position in a set of
#' anchor regions of the same width. The anchor regions are defined by their
#' genomic midpoint coordinate and the width.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param assayName Character vector, the assay(s) from which to extract values.
#' @param regionMidpoints Either a \code{GPos} object or a character vector
#'     that can be coerced into a \code{GPos} object, representing the midpoints
#'     of the desired anchor regions.
#' @param regionWidth Integer scalar, the desired width of the anchor regions.
#'     Must be an odd value.
#' @param anchorName Character scalar that gives the "sequence name" of the
#'     aligned anchor regions and will be used in generating the
#'     \code{rowRanges} of the return value.
#' @param prune Logical scalar. If \code{TRUE} (the default), samples for
#'     which there are no reads overlapping any of the anchor regions in any
#'     of the read-level assays in \code{assayName}
#'     will be completely removed from the returned \code{SummarizedExperiment}
#'     (also from \code{colData}). If \code{FALSE},
#'     such samples are retained (in assays with read-level data as a
#'     zero-column \code{NAMatrix}, in other assays as a dense matrix with a
#'     single column of NA values).
#' @param ignore.strand Logical scalar, whether to ignore the strand information
#'     when matching anchor regions with observations. Will be passed on to
#'     \code{GenomicRanges::match()}. If \code{TRUE}, a pruning step will be
#'     applied to ensure that each genomic position is represented by at most
#'     one row in the object. If multiple rows are found corresponding to the
#'     same position (on different strands), the one with the highest number of
#'     supporting reads (determined by the \code{Nvalid}, \code{Nmod},
#'     or \code{mod_prob} assay, in this preference order, depending on what
#'     is present in \code{se}) will be retained.
#' @param reverseMinusStrandRegions Logical scalar. If \code{TRUE}, data
#'     extracted from regions on the negative strand will be reversed before
#'     they are concatenated with the rest of the data, so that negative
#'     relative positions correspond to regions upstream of the region
#'     midpoint.
#' @param verbose Logical scalar. If \code{TRUE}, report on progress.
#'
#' @return A \code{SummarizedExperiment} with rows representing relative
#' positions within an anchor region (the midpoint of the region corresponds
#' to a relative position of 0) and columns representing samples. Each
#' column of the assay is an \code{NaArray} (if \code{assayName} is a
#' read-level assay) or a dense matrix (otherwise), with columns representing
#' read-anchor region (or sample-anchor region) combinations. The \code{region}
#' column of the \code{colData} records which anchor region a given column
#' corresponds to. In addition, all assays will be designated as 'read-level'
#' assays (each column represents a sample, which is in turn represented by
#' multiple columns in the actual data, corresponding to the different regions).
#' This allows downstream analysis similar to that for read-level data,
#' including flattening and plotting.
#'
#' @examples
#' library(SummarizedExperiment)
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- flattenReadLevelAssay(se)
#' ar <- getAnchorRegions(se, assayName = c("mod_prob", "FracMod", "Nvalid"),
#'                        regionMidpoints = c("chr1:6929389:-", "chr1:6935630:-"),
#'                        regionWidth = 9)
#'
#' ## Modification probabilities
#' assay(ar)
#'
#' ## Region assignment
#' colData(ar)
#'
#' @author Charlotte Soneson
#'
#' @importFrom GenomicRanges GPos match strand
#' @importFrom SparseArray NaArray cbind colSums is_nonna
#' @importFrom S4Vectors split endoapply make_zero_col_DFrame DataFrame metadata
#' @importFrom IRanges DataFrameList
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData assay
#' @importFrom methods as
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics pos strand
#' @importFrom cli cli_abort cli_warn
#'
#' @export
getAnchorRegions <- function(se,
                             assayName = "mod_prob",
                             regionMidpoints,
                             regionWidth,
                             anchorName = "anchor",
                             prune = TRUE,
                             ignore.strand = FALSE,
                             reverseMinusStrandRegions = FALSE,
                             verbose = FALSE) {

    # Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertVector(x = assayName, type = "character",
                  validValues = assayNames(se))
    names(assayName) <- assayName
    .checkSEValidity(se, verbose = FALSE)
    if (is.character(regionMidpoints)) {
        # This checks already that each region has width 1
        regionMidpoints <- as(regionMidpoints, "GPos")
    }
    .assertVector(x = regionMidpoints, type = "GPos")
    .assertScalar(x = regionWidth, type = "numeric", rngIncl = c(1, Inf))
    regionWidth <- round(regionWidth)
    if (regionWidth %% 2 == 0) {
        cli_abort("{.arg regionWidth} ({regionWidth}) must be an odd integer")
    }
    .assertScalar(x = anchorName, type = "character")
    .assertScalar(x = prune, type = "logical")
    .assertScalar(x = ignore.strand, type = "logical")
    .assertScalar(x = reverseMinusStrandRegions, type = "logical")
    .assertScalar(x = verbose, type = "logical")

    if (any(strand(regionMidpoints) == "*") && !ignore.strand) {
        cli_warn(paste0("The strand of some region midpoints is undefined. ",
                        "Setting {.arg ignore.strand} to {.code TRUE}."))
        ignore.strand <- TRUE
    }

    # If ignore.strand = TRUE, first make sure that there is no
    # position that is represented by two or more rows in the SE
    if (ignore.strand) {
        .message("Checking for positions represented by multiple rows")
        covAssays <- c("Nvalid", "Nmod", "mod_prob")
        if (!any(covAssays %in% assayNames(se))) {
            cli_abort(paste0(
                "None of the preferred coverage assays ",
                 "(Nvalid, Nmod, mod_prob) are available."))
        }
        se <- .pruneAmbiguousStrandPositions(
            se, assayName = intersect(covAssays, assayNames(se))[1],
            verbose = verbose)
    }

    # Create list of GPos objects, one for each region
    .message("Creating list of GPos objects for the regions")
    regions <- lapply(split(regionMidpoints), function(gp) {
        GPos(seqnames = seqnames(gp),
             pos = seq(pos(gp) - floor((regionWidth - 1) / 2),
                       pos(gp) + ceiling((regionWidth - 1) / 2)),
             strand = strand(gp))
    })
    names(regions) <- vapply(regions, function(reg) {
        paste0(seqnames(reg)[1], ":", min(pos(reg)), "-", max(pos(reg)),
               ":", strand(reg)[1])
    }, "")

    # Create new assays
    .message("Subsetting assays to selected regions")
    assayL <- lapply(assayName, function(atp) {
        if (atp %in% .getReadLevelAssayNames(se)) {
            # read-level assays (DataFrames with NaArrays)
            endoapply(assay(se, atp), function(mat) {
                # extract a submatrix for each region (with all reads)
                # prefix read IDs with region ID
                mats <- lapply(names(regions), function(regnm) {
                    reg <- regions[[regnm]]
                    # initialize NaArray with one row per base in the region
                    # no row names, as they won't be interpretable across
                    # regions anyway
                    namat <- NaArray(
                        dim = c(length(reg), ncol(mat)),
                        dimnames = list(NULL, paste0(regnm, "-", colnames(mat))),
                        type = "double")
                    # populate NaArray with values from observed positions
                    m <- match(reg, rowRanges(se), ignore.strand = ignore.strand)
                    newrow <- which(!is.na(m))
                    oldrow <- m[newrow]
                    i <- rep(newrow, ncol(mat))
                    j <- rep(seq_len(ncol(mat)), each = length(newrow))
                    val <- as.vector(mat[oldrow, ])
                    idx <- which(!is.na(val))
                    namat[cbind(i[idx], j[idx])] <- val[idx]
                    if (reverseMinusStrandRegions && as.character(strand(reg)[1]) == "-") {
                        # reverse relative positions
                        namat <- namat[seq(nrow(namat), 1), , drop = FALSE]
                    }
                    namat
                })
                # cbind matrices from different regions
                cbmat <- do.call(cbind, mats)
                # only keep read-region pairs with at least one non-NA value
                keep_reads <- which(colSums(is_nonna(cbmat)) > 0)
                cbmat[, keep_reads, drop = FALSE]
            })
        } else {
            # summary assays (matrices)
            # generate a DataFrame assay with dense matrices as columns, where
            # in each of these, one column is a sample-region pair
            assayDF <- make_zero_col_DFrame(nrow = regionWidth)
            for (s in colnames(assay(se, atp))) {
                # create sample-region matrix
                srmat <- do.call(cbind, lapply(names(regions), function(regnm) {
                    reg <- regions[[regnm]]
                    newmat <- matrix(NA,
                                     nrow = length(reg),
                                     ncol = 1,
                                     dimnames = list(NULL, paste0(regnm, "-", s)))
                    m <- match(reg, rowRanges(se), ignore.strand = ignore.strand)
                    newrow <- which(!is.na(m))
                    oldrow <- m[newrow]
                    newmat[newrow, 1] <- assay(se, atp)[oldrow, s]
                    if (reverseMinusStrandRegions && as.character(strand(reg)[1]) == "-") {
                        # reverse relative positions
                        newmat <- newmat[seq(nrow(newmat), 1), , drop = FALSE]
                    }
                    newmat
                }))
                assayDF[[s]] <- srmat
            }
            assayDF
        }
    })

    # Record the region corresponding to each column and add to colData
    .message("Assembling SummarizedExperiment object")
    cold <- DataFrame(
        sample = colnames(se)
    )
    for (atp in names(assayL)) {
        regs <- SimpleList(lapply(assayL[[atp]], function(m) {
            DataFrame(
                id = colnames(m),
                region = sub("^([^:]+:[0-9]+-[0-9]+:([+]|[-]|[*])).*$",
                             "\\1", colnames(m))
            )
        }))
        cold[[paste0("region_", atp)]] <- regs
    }

    # Assemble SE
    suppressWarnings({
        # currently, assigning to assays triggers a deprecation warning
        # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
        # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
        seout <- SummarizedExperiment(
            assays = DataFrameList(assayL),
            rowRanges = GPos(seqnames = anchorName,
                             pos = seq_len(regionWidth) - floor((regionWidth + 1) / 2)),
            colData = cold,
            metadata = list(readLevelData = list(
                assayNames = names(assayL),
                colDataColumns = paste0("region_", names(assayL))
            ))
        )
    })

    # Drop samples without reads if prune=TRUE
    if (prune && length(intersect(assayName, .getReadLevelAssayNames(se))) > 0) {
        keepSamples <- c()
        for (atp in intersect(assayName, .getReadLevelAssayNames(se))) {
            # check only read-level assays
            keepSamples <- union(
                keepSamples,
                colnames(assay(seout, atp)[vapply(assay(seout, atp), ncol, 0) > 0]))
        }
        seout <- seout[, keepSamples]
        .message("Dropping {ncol(se) - ncol(seout)} sample{?s} without reads")
    }

    seout
}


