#' Extract data for one or more anchor regions
#'
#' Extract assay values for each read (or sample) and each position in a set of
#' anchor regions of the same width. The anchor regions are defined by their
#' genomic midpoint coordinate and the width.
#'
#' @param se \code{SummarizedExperiment} object.
#' @param assay.type Character vector, the assay(s) from which to extract values.
#' @param regionMidpoints Either a \code{GPos} object or a character vector
#'     that can be coerced into a \code{GPos} object, representing the midpoints
#'     of the desired anchor regions.
#' @param regionWidth Integer scalar, the desired width of the anchor regions.
#' @param prune Logical scalar. If \code{TRUE} (the default), samples for
#'     which there are no reads overlapping any of the anchor regions in any
#'     of the read-level assays in \code{assay.type}
#'     will be completely removed from the returned \code{SummarizedExperiment}
#'     (also from \code{colData}). If \code{FALSE},
#'     such samples are retained (in assays with read-level data as a
#'     zero-column \code{NAMatrix}, in other assays as a dense matrix with a
#'     single column of NA values).
#' @param ignore.strand Logical scalar, whether to ignore the strand information
#'     when matching anchor regions with observations. Will be passed on to
#'     \code{GenomicRanges::match()}.
#'
#' @return A \code{SummarizedExperiment} with rows representing relative
#' positions within an anchor region (the midpoint of the region corresponds
#' to a relative position of 0) and columns representing samples. Each
#' column of the assay is an \code{NaArray} (if \code{assay.type} is a
#' read-level assay) or a dense matrix (otherwise), with columns representing
#' read-anchor region (or sample-anchor region) combinations. The \code{region}
#' column of the \code{colData} records which anchor region a given column
#' corresponds to.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
#'                  modbase = "a", verbose = FALSE)
#' se <- addReadsSummary(se, keep.reads = TRUE)
#' ar <- getAnchorRegions(se, assay.type = c("mod_prob", "FracMod", "Nvalid"),
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
#' @importFrom GenomicRanges GPos match
#' @importFrom SparseArray NaArray cbind colSums
#' @importFrom S4Vectors split endoapply make_zero_col_DFrame DataFrame
#' @importFrom IRanges DataFrameList
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData assay
#' @importFrom methods as
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics pos strand
#' @importFrom stringr str_extract
#'
#' @export
getAnchorRegions <- function(se,
                             assay.type = "mod_prob",
                             regionMidpoints,
                             regionWidth,
                             prune = TRUE,
                             ignore.strand = FALSE) {

    # Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    if ((is.numeric(assay.type) &&
         (any(assay.type < 1) || any(assay.type > length(SummarizedExperiment::assays(se))))) ||
        (is.character(assay.type) && any(!assay.type %in% SummarizedExperiment::assayNames(se)))) {
        stop("'assay.type' must be a string or integer vector specifying the ",
             "assay of se containing the data to be extracted")
    }
    if (is.character(assay.type)) {
        names(assay.type) <- assay.type
    }
    .checkSEValidity(se, verbose = FALSE)
    if (is.character(regionMidpoints)) {
        # This checks already that each region has width 1
        regionMidpoints <- as(regionMidpoints, "GPos")
    }
    .assertVector(x = regionMidpoints, type = "GPos")
    .assertScalar(x = regionWidth, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = prune, type = "logical")
    .assertScalar(x = ignore.strand, type = "logical")

    if (any(GenomicRanges::strand(regionMidpoints) == "*") &&
        !ignore.strand) {
        warning("The strand of some region midpoints is undefined ",
                "Setting `ignore.strand` to TRUE.")
        ignore.strand <- TRUE
    }

    # Create list of GPos objects, one for each region
    regions <- lapply(S4Vectors::split(regionMidpoints), function(gp) {
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
    assayL <- lapply(assay.type, function(atp) {
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
                    namat <- SparseArray::NaArray(
                        dim = c(length(reg), ncol(mat)),
                        dimnames = list(NULL, paste0(regnm, "-", colnames(mat))),
                        type = "double")
                    # populate NaArray with values from observed positions
                    m <- GenomicRanges::match(reg, rowRanges(se),
                                              ignore.strand = ignore.strand)
                    newrow <- which(!is.na(m))
                    oldrow <- m[newrow]
                    i <- rep(newrow, ncol(mat))
                    j <- rep(seq_len(ncol(mat)), each = length(newrow))
                    namat[cbind(i, j)] <- c(as.matrix(mat[oldrow, ]))
                    namat
                })
                # cbind matrices from different regions
                cbmat <- do.call(SparseArray::cbind, mats)
                # only keep read-region pairs with at least one non-NA value
                keep_reads <- which(SparseArray::colSums(cbmat, na.rm = TRUE) > 0)
                cbmat[, keep_reads, drop = FALSE]
            })
        } else {
            # summary assays (matrices)
            # generate a DataFrame assay with dense matrices as columns, where
            # in each of these, one column is a sample-region pair
            assayDF <- S4Vectors::make_zero_col_DFrame(nrow = regionWidth)
            for (s in colnames(assay(se, atp))) {
                # create sample-region matrix
                srmat <- do.call(cbind, lapply(names(regions), function(regnm) {
                    reg <- regions[[regnm]]
                    newmat <- matrix(NA,
                                     nrow = length(reg),
                                     ncol = 1,
                                     dimnames = list(NULL, paste0(regnm, "-", s)))
                    m <- GenomicRanges::match(reg, rowRanges(se),
                                              ignore.strand = ignore.strand)
                    newrow <- which(!is.na(m))
                    oldrow <- m[newrow]
                    newmat[newrow, 1] <- assay(se, atp)[oldrow, s]
                    newmat
                }))
                assayDF[[s]] <- srmat
            }
            assayDF
        }
    })

    # Record the region corresponding to each column and add to colData
    cold <- DataFrame(
        sample = colnames(se)
    )
    for (atp in names(assayL)) {
        regs <- SimpleList(lapply(assayL[[atp]], function(m) {
            S4Vectors::DataFrame(
                id = colnames(m),
                region = stringr::str_extract(
                    colnames(m),
                    paste(sub("*", "\\*", sub("+", "\\+",
                                              paste0("^", names(regions)),
                                              fixed = TRUE), fixed = TRUE),
                          collapse = "|"))
            )
        }))
        cold[[paste0("region_", atp)]] <- regs
    }

    # Assemble SE
    suppressWarnings({
        # currently, assigning to assays triggers a depreceation warning
        # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
        # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
        seout <- SummarizedExperiment::SummarizedExperiment(
            assays = DataFrameList(assayL),
            rowData = S4Vectors::DataFrame(
                relpos = seq_len(regionWidth) - floor((regionWidth + 1) / 2)
            ),
            colData = cold,
            metadata = list()
        )
    })

    # Drop samples without reads if prune=TRUE
    if (prune) {
        keepSamples <- c()
        for (atp in intersect(.getReadLevelAssayNames(se),
                              assayNames(seout))) {
            # check only read-level assays
            keepSamples <- union(
                keepSamples,
                colnames(assay(seout, atp)[vapply(assay(seout, atp), ncol, 0) > 0]))
        }
        seout <- seout[, keepSamples]
    }

    seout
}


