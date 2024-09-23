#' Extract data for one or more anchor regions
#'
#' Extract modification probabilities (or other read-level values) for each
#' read and each position in a set of anchor regions of the same width. The
#' anchor regions are defined by their genomic midpoint coordinate as well as
#' the width.
#'
#' @param se \code{SummarizedExperiment} object with at least one read-level
#'     assay.
#' @param assay.type Character scalar, the assay from which to extract
#'     read-level information.
#' @param regionMidpoints Either a \code{GPos} object or a character vector
#'     that can be coerced into a \code{GPos} object, representing the midpoints
#'     of the desired anchor regions.
#' @param regionWidth Integer scalar, the desired width of the anchor regions.
#'
#' @return A \code{SummarizedExperiment} with rows representing relative
#' positions within an anchor region (the midpoint of the region corresponds
#' to a relative position of 0) and columns representing samples. Each
#' column of the assay is an \code{NaArray}, with columns representing
#' read-anchor region combinations. The \code{region} column of the
#' \code{colData} records which anchor region a given column corresponds to.
#'
#' @examples
#' modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
#'                            package = "footprintR")
#' se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
#'                  modbase = "a", verbose = FALSE)
#' ar <- getAnchorRegions(se, assay.type = "mod_prob",
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
#' @importFrom GenomicRanges GPos
#' @importFrom SparseArray NaArray
#' @importFrom S4Vectors endoapply
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowData
#'
#' @export
getAnchorRegions <- function(se, assay.type = "mod_prob", regionMidpoints,
                             regionWidth) {
    ## TODO: Allow this to be used also for non-read-level assays

    # Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    if (length(assay.type) != 1L ||
        (is.numeric(assay.type) &&
         (assay.type < 1 || assay.type > length(SummarizedExperiment::assays(se)))) ||
        (is.character(assay.type) && !assay.type %in% SummarizedExperiment::assayNames(se))) {
        stop("'assay.type' must be a string or integer scalar specifying the ",
             "assay of se containing the read-level data to be summarized.")
    }
    .assertScalar(x = regionWidth, type = "numeric", rngIncl = c(1, Inf))
    if (is.character(regionMidpoints)) {
        # This checks already that each region has width 1
        regionMidpoints <- as(regionMidpoints, "GPos")
    }
    .assertVector(x = regionMidpoints, type = "GPos")

    # Create list of GPos objects, one for each region
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

    # Create new NaArray for each sample
    newmat <- endoapply(assay(se, assay.type), function(mat) {
        mats <- lapply(names(regions), function(regnm) {
            reg <- regions[[regnm]]
            namat <- SparseArray::NaArray(
                dim = c(length(reg), ncol(mat)),
                dimnames = list(NULL, paste0(regnm, "-", colnames(mat))),
                type = "double")
            m <- GenomicRanges::match(reg, rowRanges(se))
            ## TODO: Make sure that this is robust e.g. if there are no overlaps
            newrow <- which(!is.na(m))
            oldrow <- m[newrow]
            i <- rep(newrow, ncol(mat))
            j <- rep(seq_len(ncol(mat)), each = length(newrow))
            namat[cbind(i, j)] <- c(as.matrix(mat[oldrow, ]))
            namat
            ## Q: Should we keep all reads, or only those overlapping a given region?
            ## TODO: Here we can also summarize by region
        })
        do.call(SparseArray::cbind, mats)
    })

    # Record the region corresponding to each column
    # Can be a regular vector if summarization is applied
    regs <- SimpleList(lapply(newmat, function(m) {
        data.frame(id = colnames(m),
                   region = stringr::str_extract(colnames(m),
                                                 paste(paste0("^", names(regions)),
                                                       collapse = "|"))
        )
    }))

    # Assemble and return SE
    al <- list(tmp = newmat)
    if (is.character(assay.type)) {
        names(al) <- assay.type
    } else {
        names(al) <- NULL
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = al,
        rowData = S4Vectors::DataFrame(
            relpos = seq_len(nrow(newmat)) - floor((regionWidth + 1) / 2)
        ),
        colData = S4Vectors::DataFrame(
            sample = colnames(se),
            region = regs
        ),
        metadata = list()
    )
}


