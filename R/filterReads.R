#' Filter reads
#' 
#' @param se A \code{SummarizedExperiment} object.
#' @param assay.type.read A character scalar providing the name of a read-level
#'     assay in \code{se}. This assay will be used to extract read names, as 
#'     well as to filter out any read that is not overlapping any of the 
#'     positions in the object.
#' @param readInfoCol A character scalar providing the name of the column in 
#'     \code{colData} that contains read info. Can be \code{NULL} if no such 
#'     column exists. 
#' @param qcCol A character scalar providing the name of the column in 
#'     \code{colData} that contains quality metrics (calculated by 
#'     \code{calcReadStats}). Can be \code{NULL} if no such column exists. 
#' @param minQscore A numeric scalar representing the smallest acceptable 
#'     read-level Qscore. Reads with Qscore below this value will be filtered
#'     out.
#' @param minEntropy A numeric scalar representing the smallest acceptable 
#'     read-level entropy. Reads with entropy below this value will be filtered
#'     out.
#' @param minReadLength A numeric scalar representing the smallest acceptable 
#'     read length. Reads that are shorter than this value will be filtered
#'     out.
#' @param minAlignedLength A numeric scalar representing the smallest acceptable 
#'     aligned length. Reads with aligned length shorter than this value will 
#'     be filtered out.
#' @param minAlignedFraction A numeric scalar representing the smallest 
#'     acceptable aligned fraction of a read. Reads where the aligned fraction 
#'     is smaller than this value will be filtered out.
#' @param prune A logical scalar. If \code{TRUE} (the default), samples for
#'     which the filtering retains none of the reads will be completely removed
#'     from the returned \code{SummarizedExperiment} (also from \code{colData}
#'     and from assays that do not store read-level data). If \code{FALSE},
#'     such samples are retained (in the assays with read-level data as a
#'     zero-column \code{SparseMatrix}).
#' 
#' @author Charlotte Soneson
#' @export
#' 
#' @returns A filtered \code{SummarizedExperiment} object. The metadata of this
#' object contains a slot named \code{filteredOutReads}, which tabulate all
#' reads that are filtered out, together with the reason(s) for exclusion. 
#' 
#' @examples 
#' library(SummarizedExperiment)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfile = modbamfile, regions = "chr1:6920000-6995000",
#'            modbase = "a", verbose = TRUE)
#' se <- addReadStats(se, name = "QC")
#' sefilt <- filterReads(se, minQscore = 14, minAlignedLength = 10000)
#' 
#' 
#' @importFrom SparseArray SVT_SparseArray rowSums
#' 
filterReads <- function(se, assay.type.read = "mod_prob", 
                        readInfoCol = "read_info", qcCol = "QC",
                        minQscore = 0, minEntropy = 0, 
                        minReadLength = 0, minAlignedLength = 0, 
                        minAlignedFraction = 0, prune = TRUE) {
    ## Input checks
    .assertVector(x = se, type = "SummarizedExperiment")
    .checkSEValidity(se)
    .assertScalar(x = assay.type.read, type = "character", 
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = readInfoCol, type = "character", allowNULL = TRUE)
    .assertScalar(x = qcCol, type = "character", allowNULL = TRUE)
    .assertScalar(x = minQscore, type = "numeric")
    .assertScalar(x = minEntropy, type = "numeric")
    .assertScalar(x = minReadLength, type = "numeric")
    .assertScalar(x = minAlignedLength, type = "numeric")
    .assertScalar(x = minAlignedFraction, type = "numeric")
    .assertScalar(x = prune, type = "logical")
    
    ## Initialize sparse logical array for each sample, which will be TRUE 
    ## for reads that are filtered out with respect to the different criteria
    ## Remark: Could move this to a global constant
    filterNames <- c("Qscore", "Entropy", "ReadLength", "AlignedLength",
                     "AlignedFraction", "AllNA")
    readsToRemove <- lapply(
        structure(colnames(se), names = colnames(se)),
        function(nm) {
            SparseArray::SVT_SparseArray(
                dim = c(ncol(assay(se, assay.type.read)[[nm]]), 
                        length(filterNames)), 
                dimnames = list(colnames(assay(se, assay.type.read)[[nm]]),
                                filterNames),
                type = "logical"
        )}
    )
    
    for (nm in colnames(se)) {
        ## Extract read info and QC
        if (!is.null(readInfoCol) && !is.null(se[[readInfoCol]]) && 
            (nm %in% names(se[[readInfoCol]]))) {
            ri <- se[[readInfoCol]][[nm]]
        } else {
            ri <- NULL
        }
        if (!is.null(qcCol) && !is.null(se[[qcCol]]) && 
            (nm %in% names(se[[qcCol]]))) {
            qc <- se[[qcCol]][[nm]]
        } else {
            qc <- NULL
        }
        
        ## Quality score
        if (!is.null(ri) && "qscore" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$qscore < minQscore), 
                                "Qscore"] <- TRUE
        }
        
        ## KS entropy
        if (!is.null(qc) && "EntrModProb" %in% colnames(qc)) {
            readsToRemove[[nm]][which(qc$SEntrModProb < minEntropy),
                                "Entropy"] <- TRUE
        }
        
        ## Fail rate
        
        
        ## Read length
        if (!is.null(ri) && "read_length" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$read_length < minReadLength),
                                "ReadLength"] <- TRUE
        }
        
        ## Aligned length
        if (!is.null(ri) && "aligned_length" %in% colnames(ri)) {
            readsToRemove[[nm]][ which(ri$aligned_length < minAlignedLength),
                                 "AlignedLength"] <- TRUE
        }
        
        ## Aligned fraction
        if (!is.null(ri) && "aligned_fraction" %in% colnames(ri)) {
            readsToRemove[[nm]][which(ri$aligned_fraction < minAlignedFraction),
                                "AlignedFraction"] <- TRUE
        }

        ## NA in all positions
        readsToRemove[[nm]][colnames(
            assay(se, assay.type.read)[[nm]][, SparseArray::colSums(
                assay(se, assay.type.read)[[nm]], 
                na.rm = TRUE) == 0]), 
            "AllNA"] <- TRUE
    }
    
    ## Subset
    readsToRemove <- lapply(readsToRemove, function(rr) {
        rr[SparseArray::rowSums(rr, na.rm = TRUE) > 0, ]
    })
    sesub <- subsetReads(se = se, reads = lapply(readsToRemove, rownames),
                         prune = prune, invert = TRUE)
    metadata(sesub)$filteredOutReads <- readsToRemove
    
    sesub
}
