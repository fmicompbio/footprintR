# -- main exported functions ---------------------------------------------------

#' Calculate and segment scores for footprints and add them to a SummarizedExeriment
#'
#' @description
#' Given a read-level footprinting data and a weight vector defining a specific
#' footprint, calculate footprinting scores, identify high-scoring segments
#' and add them to the input SummarizedExeriment object.
#'
#' Typicallly, this can be done using \code{addFootprints} which will:
#' \enumerate{
#'     \item \emph{Calculate footprint scores} for each individual read using
#'         \code{calcFootprintScores}, given a \code{SummarizedExperiment} with
#'         a read-level assay and a vector of weights that defines the
#'         footprint (see \code{wgt} parameter).
#'     \item \emph{Segment footprint scores} to identify high-scoring segments
#'         in the scores calculated at step 1 using
#'         \code{segmentFootprintScores}.
#'     \item \emph{Store high-scoring segments} in the
#'         \code{SummarizedExperiment} as a column in \code{colData} of the
#'         object.
#' }
#'
#' If required, the functions \code{calcFootprintScores} or
#' \code{segmentFootprintScores} can also be called directly to perform only
#' one of the steps above.
#'
#' The footprint is defined by a weight vector \code{wgt} which describes its
#' length in bases and the shape of the expected modification profile. The score
#' is based on a cross-correlation of the modification probabilities (centered
#' by subtracting \code{0.5}) with \code{wgt}, weighted by the minimum of
#' \code{minweight} and the modification probability.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param wgt Numeric vector with weights that define the footprint to score.
#'     Typically centered at zero.
#' @param assayName A character scalar providing the name of a read-level
#'     assay in \code{se} that contains modification probabilities.
#' @param minconf Numeric scalar giving the minimal confidence of a base to be
#'     included in the calculation of footprint scores.
#' @param minweight Numeric scalar giving the minimal weight for position in the
#'     calculation of footprint scores.
#' @param minperiod,maxperiod Numeric scalars that define the band-pass
#'     filter parameters used to smooth the footprint scores for segmentation.
#'     \code{minperiod} and \code{maxperiod} give the minimal and maximal
#'     periods (in base pairs) for the critical frequencies in the score signal
#'     that should be retained. The default values are suitable to retain
#'     nucleosomal signals. The filtering is performed using
#'     \code{\link[signal]{filtfilt}} and thus requires the \code{signal}
#'     package to be installed.
#' @param thresh A numeric scalar giving the minimal smoothed score of a
#'     footprint. Higher values make the footprint detection more stringent.
#' @param lenRange A numeric vector with two elements giving the minimal and
#'     maximal number of consecutive score values (base pairs) that must be
#'     greater than \code{thresh} for a region to be included in the
#'     returned footprints.
#' @param width A numeric scalar giving the width of returned footprints.
#'     Footprints are created by placing a region of length \code{width}
#'     centered on the maximal smoothed score value in each identified
#'     high-scoring segment.
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
#' # "wgt": footprint weights for nucleosome with flanks
#' # - weights are centered around zero (starting with -0.5 and 0.5)
#' # - weights are adjusted for their relative frequency, to give similar
#' #   importance to flanks and protected region
#' # - weights are repeated for bases in the (un)modified parts (15 bases in
#' #   flanks and 140 bases protected by the nucleosome)
#' wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
#' se <- addFootprints(se, wgt, thresh = 0.05, name = "nucl")
#' se$nucl
#'
#' # only calculate score
#' scoresList <- calcFootprintScores(se, wgt)
#' lapply(scoresList, nrow)
#' lapply(scoresList, head)
#'
#' # segment calculated scores
#' irl <- segmentFootprintScores(scoresList, thresh = 0.05)
#' irl$s1
#'
#' @returns
#' \describe{
#'     \item{For \code{addFootprints}}{, a \code{SummarizedExperiment} object
#'     with footprints in \code{colData(se)[[name]]}.}
#'     \item{For \code{calcFootprintScores}}{, a \code{list} with one entry per
#'     sample (column) in \code{se}. The list elements are \code{data.frame}s
#'     with columns, \code{readId}, \code{pos}, \code{pmod} and \code{score}, in
#'     (continuous) base space from \code{min(pos(rowRanges(se)))} to
#'     \code{max(pos(rowRanges(se)))} (only including positions covered by reads
#'     in a given sample).}
#'     \item{For segmentFootprintScores}{, a \code{list} (over samples) of
#'     \code{\link[IRanges]{IRangesList}} elements (over reads). Each element
#'     contains one \code{\link[IRanges]{IRanges}} object with the identified
#'     footprints of a single read (all of length \code{with} bases).}
#' }
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
#' @rdname footprintScoring
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

    ## Calculate scores
    scoresL <- calcFootprintScores(se = se, wgt = wgt, assayName = assayName,
                                   minconf = minconf, minweight = minweight,
                                   verbose = verbose)

    ## Segment scores
    irlL <- segmentFootprintScores(scoresList = scoresL,
                                   minperiod = minperiod,
                                   maxperiod = maxperiod,
                                   thresh = thresh,
                                   lenRange = lenRange,
                                   width = width,
                                   verbose = verbose)

    ## add to se
    cd <- colData(se)
    cd[[name]] <- irlL
    colData(se) <- cd
    metadata(se)$readLevelData$colDataColumns <- union(
        metadata(se)$readLevelData$colDataColumns, name
    )

    return(se)
}

#' @importFrom SummarizedExperiment assay
#' @importFrom SparseArray nnawhich nnavals
#' @importFrom BiocGenerics start
#' @importFrom dplyr mutate group_by group_modify ungroup
#' @importFrom rlang .data
#'
#' @export
#' @rdname footprintScoring
calcFootprintScores <- function(se,
                                wgt,
                                assayName = "mod_prob",
                                minconf = 0.7,
                                minweight = 0.05,
                                verbose = TRUE) {

    ## Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .checkSEValidity(se)
    .assertVector(x = wgt, type = "numeric")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = minconf, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minweight, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = verbose, type = "logical")

    .message("calculating footprint scores")

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
                mutate(readId = factor(.data$readId, levels = unique(.data$readId))) |>
                group_by(.data$readId) |>
                group_modify(~ calcFootprintScoreForRead(.x$pos, .x$pmod,
                                                         wgt = wgt,
                                                         minconf = minconf,
                                                         minweight = minweight)) |>
                ungroup() |>
                as.data.frame()
        } else {
            data.frame(readId = factor(),
                       pos = integer(0),
                       pmod = numeric(0),
                       score = numeric(0))
        }
    })

    return(scoresL)
}

#' @param scoresList A \code{list} of \code{data.frame}s with footprint scores
#'     as returned by \code{calcFootprintScores}.
#'
#' @importFrom IRanges IRanges IRangesList coerce Views viewApply viewWhichMaxs
#'     resize start width
#' @importFrom dplyr group_by mutate ungroup select
#' @importFrom rlang .data
#'
#' @export
#' @rdname footprintScoring
segmentFootprintScores <- function(scoresList,
                                   minperiod = 50,
                                   maxperiod = 300,
                                   thresh = 0.15,
                                   lenRange = c(10, 80),
                                   width = 140,
                                   verbose = TRUE) {
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
    .assertScalar(x = verbose, type = "logical")
    .assertPackagesAvailable(pkgs = "signal")

    .message("segmenting footprint scores")

    ## Smooth scores
    irListList <- lapply(
        structure(names(scoresList), names = names(scoresList)),
        function(nm) {
            if (nrow(scoresList[[nm]]) > 0) {
                # smooth scores using a band-pass filter
                dat <- scoresList[[nm]] |>
                    select("readId", "pos", "score") |>
                    group_by(.data$readId) |>
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

# -- helper functions ----------------------------------------------------------
#' Smooth scores using band-pass filter
#'
#' @noRd
#' @keywords internal
.filterScores <- function(score, minperiod, maxperiod) {
    .assertPackagesAvailable(pkgs = "signal")
    Wn <- 1 / c(maxperiod, minperiod)
    testar <- signal::butter(n = 3, W = Wn, type = "pass")

    nnaIndex <- which(!is.na(score))
    nnaScores <- score[nnaIndex]
    nnaSScores <- signal::filtfilt(testar, c(
        rev(nnaScores), nnaScores, rev(nnaScores)))[
            (length(nnaScores) + 1):(2 * length(nnaScores))]
    sscore <- rep(NA, length(score))
    sscore[nnaIndex] <- nnaSScores

    return(sscore)
}
