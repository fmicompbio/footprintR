#' Calculate footprinting scores for a defined footprint
#'
#' This function calculates the footprinting scores corresponding to a
#' footprint in the form of a weight vector \code{wgt} for all reads in a
#' given read-level assay with modification probabilities. The score is based
#' on a cross-correlation of the modification probabilities (centered by 
#' subtracting \code{0.5}) with \code{wgt}, weighted by the minimum of 
#' \code{minweight} and the modification probability.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param wgt Numeric vector with weights that define the footprint to score.
#'     Typically centered at zero.
#' @param assayName A character scalar providing the name of a read-level
#'     assay in \code{se} that contains modification probabilities.
#' @param minconf Numeric scalar giving the minimal confidence of a base to be
#'     included in the calculation.
#' @param minweight Numeric scalar giving the minimal weight for position.
#'
#' @examples
#' # load example data
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#'
#' # footprint weights for nucleosome with flanks
#' # - weights are centered around zero (starting with -0.5 and 0.5)
#' # - weights are adjusted for their relative frequency, to give similar
#' #   importance to flanks and protected region
#' # - weights are repeated for bases in the (un)modified parts (15 bases in
#' #   flanks and 140 bases protected by the nucleosome)
#' wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
#'
#' # calculate scores
#' scoresList <- calcFootprintScores(se, wgt)
#' lapply(scoresList, nrow)
#' lapply(scoresList, head)
#'
#' @returns
#' A \code{list} with one entry per sample. The elements are \code{data.frame}s
#' with columns, \code{readId}, \code{pos}, \code{pmod} and \code{score}, in
#' (continuous) base space from \code{min(pos(rowRanges(se)))} to
#' \code{max(pos(rowRanges(se)))} (only including positions covered by reads in
#' a given sample).
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom SparseArray nnawhich nnavals
#' @importFrom BiocGenerics start
#' @importFrom dplyr group_by group_modify ungroup
#'
#' @noRd
#' @keywords internal
calcFootprintScores <- function(se,
                                wgt,
                                assayName = "mod_prob",
                                minconf = 0.7,
                                minweight = 0.05) {

    ## Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .checkSEValidity(se)
    .assertVector(x = wgt, type = "numeric")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = minconf, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minweight, type = "numeric", rngIncl = c(0, Inf))

    ## Extract list of data.frames from assay
    adatList <- lapply(assay(se, assayName), function(adat) {
        idx <- nnawhich(adat, arr.ind = TRUE)
        data.frame(
            readId = colnames(adat)[idx[, 2]],
            pos = start(se)[idx[, 1]],
            pmod = nnavals(adat)
        )
    })

    ## Calculate footprint scores
    scoresL <- lapply(adatList, function(adf) {
        if (nrow(adf) > 0) {
            adf |>
                group_by(readId) |>
                group_modify(~ calcFootprintScoreForRead(.x$pos, .x$pmod,
                                                         wgt = wgt,
                                                         minconf = minconf,
                                                         minweight = minweight)) |>
                ungroup() |>
                as.data.frame()
        } else {
            data.frame(readId = character(0),
                       pos = integer(0),
                       pmod = numeric(0),
                       score = numeric(0))
        }
    })

    return(scoresL)
}

#' Smooth scores using band-pass filter
#'
#' @importFrom signal butter filtfilt
#' @noRd
#' @keywords internal
.filterScores <- function(score, minperiod, maxperiod) {
    Wn <- 1 / c(maxperiod, minperiod)
    testar <- butter(n = 3, W = Wn, type = "pass")

    nnaIndex <- which(!is.na(score))
    nnaScores <- score[nnaIndex]
    nnaSScores <- filtfilt(testar, c(
        rev(nnaScores), nnaScores, rev(nnaScores)))[
            (length(nnaScores) + 1):(2 * length(nnaScores))]
    sscore <- rep(NA, length(score))
    sscore[nnaIndex] <- nnaSScores

    return(sscore)
}

#' Segment footprinting scores
#'
#' Segment footprinting scores, typically calculated by
#' \code{calcFootprintScores}, into high-scoring segments of a constant width.
#'
#' @param scoresList A \code{list} of \code{data.frame}s with footprint score
#'     as returned by \code{calcFootprintScores}.
#' @param minperiod,maxperiod Numeric scalars that define the band-pass
#'     filter parameters used to smooth the footprint scores. \code{minperiod}
#'     and \code{maxperiod} give the minimal and maximal periods (in base pairs)
#'     for the critical frequencies in the score signal that should be retained.
#'     The default values are suitable to retain nucleosomal signals.
#' @param thresh A numeric scalar giving the minimal score of a footprint.
#'     Higher values make the footprint detection more stringent.
#' @param lenRange A numeric vector with two elements giving the minimal and
#'     maximal number of consecutive score values (base pairs) that must be
#'     greater than \code{thresh} for a region to be included in the
#'     returned footprints.
#' @param width A numeric scalar giving the width of returned footprints.
#'
#' @returns
#' A \code{list} of \code{\link[IRanges]{IRangesList}} elements. Each element
#' contains one \code{\link[IRanges]{IRanges}} object with the identified
#' footprints of a single read (all of \code{with} width).
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @examples
#' # load example data
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#'
#' # calculate footprint scores
#' wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
#' scoresList <- calcFootprintScores(se, wgt)
#'
#' # segment the scores
#' irl <- segmentFootprintScores(scoresList, thresh = 0.05)
#' irl$s1
#'
#' @importFrom IRanges IRanges IRangesList coerce Views viewApply viewWhichMaxs
#'     resize start width
#' @importFrom dplyr group_by mutate ungroup select
#'
#' @noRd
#' @keywords internal
segmentFootprintScores <- function(scoresList,
                                   minperiod = 50,
                                   maxperiod = 300,
                                   thresh = 0.15,
                                   lenRange = c(10, 80),
                                   width = 140) {
    ## Check arguments
    .assertVector(x = scoresList, type = "list")
    for (i in seq_along(scoresList)) {
        .assertVector(x = scoresList[[i]], type = "data.frame")
        .assertVector(x = colnames(scoresList[[i]]), type = "character",
                      validValues = c("readId", "pos", "pmod", "score"))
    }
    .assertScalar(x = minperiod, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = maxperiod, type = "numeric", rngIncl = c(minperiod, Inf))
    .assertScalar(x = thresh, type = "numeric")
    .assertVector(x = lenRange, type = "numeric", rngIncl = c(1, Inf), len = 2L)
    .assertScalar(x = width, type = "numeric", rngIncl = c(1, Inf))
    .assertPackagesAvailable(pkgs = "signal")

    ## Smooth scores
    if (require("signal", quietly = TRUE)) {
        irListList <- lapply(
            structure(names(scoresList), names = names(scoresList)),
            function(nm) {
                if (nrow(scoresList[[nm]]) > 0) {
                    # smooth scores using a band-pass filter
                    dat <- scoresList[[nm]] |>
                        select(readId, pos, score) |>
                        group_by(readId) |>
                        mutate(sscore = .filterScores(score,
                                                      minperiod = minperiod,
                                                      maxperiod = maxperiod)) |>
                        ungroup()
                    iByReadId <- split(seq.int(nrow(dat)), 
                                       dat$readId)[unique(dat$readId)]
                    midList <- lapply(iByReadId, function(i) {
                        
                        # segment smoothed scores
                        sscore <- dat$sscore[i]
                        irpos <- as(!is.na(sscore) & sscore > thresh, "IRanges")
                        irpos <- irpos[width(irpos) >= lenRange[1] &
                                           width(irpos) < lenRange[2]]
                        if (length(irpos)) {
                            xposmax <- viewWhichMaxs(Views(sscore, irpos))
                        } else {
                            xposmax <- numeric(0)
                        }
                        
                        # return midpoints (location of score maxima)
                        dat$pos[i][xposmax]
                    })
                    
                    # create IRangesList
                    return(do.call(IRangesList,
                                   lapply(midList, function(mid) {
                                       resize(IRanges(start = mid, width = 1L),
                                              width = width, fix = "center")
                                   })))
                } else {
                    return(IRangesList())
                }
            })
        return(irListList)
    }
}

#' Identify and add footprint regions to a SummarizedExperiment
#'
#' This function identifies regions in each read that score highly for
#' a footprint provided by a weight vector \code{wgt} and adds them to a
#' column in \code{colData(se)} named \code{name}.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param wgt Numeric vector with weights that define the footprint to score.
#'     Typically centered at zero.
#' @param assayName A character scalar providing the name of a read-level
#'     assay in \code{se} that contains modification probabilities.
#' @param minconf,minweight,minperiod,maxperiod,thresh,lenRange,width 
#'     Additional arguments passed to helper functions that calculate and 
#'     segment footprint scores.
#' @param name Character scalar giving the column name in \code{colData(se)} in
#'     which the high-scoring footprints are stored as a list (over samples) or
#'     \code{\link[IRanges]{IRangesList}}s (over reads).
#' @param verbose If \code{TRUE}, report on progress.
#'
#' @examples
#' # load example data
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#'
#' # add nucleosome footprints
#' wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
#' se <- addFootprints(se, wgt, thresh = 0.05, name = "nucl")
#' se$nucl
#'
#' @returns
#' A \code{SummarizedExperiment} object.
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom SparseArray nnawhich nnavals
#' @importFrom BiocGenerics start
#' @importFrom dplyr group_by group_modify ungroup
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export
addFootprints <- function(se,
                          wgt,
                          assayName = "mod_prob",
                          minconf = 0.7,
                          minweight = 0.05,
                          minperiod = 50,
                          maxperiod = 300,
                          thresh = 0.15,
                          lenRange = c(10, 80),
                          width = 140,
                          name = "footprints",
                          verbose = TRUE) {
    ## Check arguments
    .assertScalar(x = name, type = "character")
    .assertScalar(x = verbose, type = "logical")
    # arglist <- list(...)
    # remark: when working with ... and do.call below, we make ~5 copies of
    #         se and as a result the call will be ~3s instead of ~0.3s
    #         not sure yet what the most elegant path is here...

    ## Calculate scores
    .message("calculating scores")
    # argsScore <- c(list(se = se, wgt = wgt, assayName = assayName),
    #                arglist[intersect(names(arglist),
    #                                  c("minconf", "minweight"))])
    # scoresL <- do.call(calcFootprintScores, argsScore)
    scoresL <- calcFootprintScores(se = se, wgt = wgt, assayName = assayName,
                                   minconf = minconf, minweight = minweight)

    ## Segment scores
    .message("segmenting scores")
    # argsSeg <- c(list(scoresList = scoresL),
    #              arglist[intersect(names(arglist),
    #                                c("minperiod", "maxperiod", "thresh",
    #                                  "lenRange", "width"))])
    # irlL <- do.call(segmentFootprintScores, argsSeg)
    irlL <- segmentFootprintScores(scoresList = scoresL,
                                   minperiod = minperiod,
                                   maxperiod = maxperiod,
                                   thresh = thresh,
                                   lenRange = lenRange,
                                   width = width)

    ## add to se
    .message("adding segments")
    # colData(se)[[name]] <- irlL
    cd <- colData(se)
    cd[[name]] <- irlL
    colData(se) <- cd
    metadata(se)$readLevelData$colDataColumns <- c(
        metadata(se)$readLevelData$colDataColumns, name
    )

    return(se)
}
