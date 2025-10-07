#' Estimate footprint lengths using nomeR
#'
#' Estimate footprint length spectrum by calling \code{countStatePairs}
#' followed by \code{nomeR::get_ftp_inference_summary}.
#'
#' @inheritParams countStatePairs
#' @param ftpLengths Integer vector of possible footprint lengths to score.
#' @param nomeRargs List providing additional arguments to
#'     \code{nomeR::get_ftp_inference_summary}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' fpl <- estimateFootprintLengths(bamfile = modbamfile,
#'                                 regions = "chr1",
#'                                 modbase = "a",
#'                                 verbose = FALSE,
#'                                 BPPARAM = BiocParallel::SerialParam())
#' fpl$PLOT
#' head(fpl$ESTIMATES$ftp_abundance_estimates)
#' fpl$ESTIMATES$ftp_protect_prob_estimate
#' fpl$ESTIMATES$bg_protect_prob_estimate
#'
#' @importFrom utils modifyList
#' @importFrom dplyr rename
#'
estimateFootprintLengths <- function(bamfile,
                                     regions = ".",
                                     modbase,
                                     threshUnmod = 0.5,
                                     threshMod = 0.5,
                                     windowSize = 200,
                                     minMapQ = 0,
                                     minAlignedLength = 0,
                                     ftpLengths = 20:200,
                                     nomeRargs = list(),
                                     BPPARAM = MulticoreParam(4L),
                                     verbose = FALSE) {
    # Check arguments (note that most arguments will be checked by
    # the call to countStatePairs)
    .assertVector(x = ftpLengths, type = "numeric", rngIncl = c(1, Inf))
    .assertVector(x = nomeRargs, type = "list")
    .assertPackagesAvailable("nomeR")

    # get state pair counts
    coocc <- countStatePairs(
        bamfile = bamfile,
        regions = regions,
        modbase = modbase,
        threshUnmod = threshUnmod,
        threshMod = threshMod,
        windowSize = windowSize,
        minMapQ = minMapQ,
        minAlignedLength = minAlignedLength,
        BPPARAM = BPPARAM,
        verbose = verbose
    )

    # rename columns for nomeR
    coocc <- as.data.frame(coocc) |>
        dplyr::rename(N00 = mod_mod,
                      N01 = mod_unmod,
                      N10 = unmod_mod,
                      N11 = unmod_unmod)

    # infer footprints with nomeR
    nomeRargs <- modifyList(
        x = nomeRargs,
        val = list(output_samples = 1000, iter = 10000, grad_samples = 1,
                   tol_rel_obj = 0.01, algorithm = "meanfield",
                   max_nruns = 3, refresh = ifelse(verbose, 100, 0))
    )
    nomeRargs$cooc_ctable <- coocc
    nomeRargs$ftp_length <- ftpLengths
    fpout <- do.call(nomeR::infer_footprints_vb, nomeRargs)
    fpout <- nomeR::get_ftp_inference_summary(
        infer_stanfit = fpout,
        plot = TRUE,
        show_plot = FALSE,
        suggest_ftps = FALSE
    )
    # Q: Should we return the suggested footprints here already? And use them in the addFootprintsNomeR function unless the user provides their own?
    fpout
}

#' Add footprints obtained using nomeR
#'
#' Given an estimated footprint length spectrum, estimate the positions of
#' footprints of different length in individual reads.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param assayName Character scalar, the name of the assay containing
#'     read-level data.
#' @param fpLengths The output of \code{estimateFootprintLengths}.
#' @param fpGroups A two-column matrix, where each row corresponds to one
#'     "group" of footprint lengths that should be considered equivalent,
#'     and the two columns correspond to the minimum and maximum footprint
#'     length in the group. For example, if all footprints between 120 and
#'     140 nt should be considered 'nucleosomes', there should be a row in
#'     \code{fpGroups} with rowname equal to 'nucleosomes', where the
#'     first element is 120 and the second element is 140.
#' @param threshUnmod,threshMod Numeric scalars used to classify observations
#'     as modified (modification probability >= threshMod), unmodified
#'     (modification probability < threshUnmod) or unknown (otherwise).
#' @param minProb Numeric scalar
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object that
#'     controls the number of parallel CPU threads to use.
#' @param verbose A logical scalar. If \code{TRUE}, report on progress.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "footprintR")
#' fpl <- estimateFootprintLengths(bamfile = modbamfile,
#'                                 regions = "chr1",
#'                                 modbase = "a",
#'                                 verbose = TRUE,
#'                                 BPPARAM = BiocParallel::SerialParam())
#' fpl$PLOT
#' fpGroups <- rbind(ftp30_39 = c(30, 39),
#'                   ftp145_155 = c(145, 155))
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6930000-6975000",
#'                  modbase = "a", verbose = TRUE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se <- addFootprintsNomeR(se = se, fpLengths = fpl, fpGroups = fpGroups,
#'                          minProb = 0.35)
#'
#' # plot identified footprints
#' plotRegion(se, region = "chr1:6930400-6931300",
#'            referenceCoordinate = 6935400,
#'            tracks = list(list(trackType = "Lollipop", trackData = "mod_prob",
#'                               footprintColumns = "ftp30_39",
#'                               arglistFootprints = list(fill = "pink",
#'                                                        height = 0.8))))
#'
#' @importFrom SummarizedExperiment assay colData colnames
#' @importFrom S4Vectors metadata
#' @importFrom IRanges IRanges IRangesList width
#' @importFrom BiocGenerics start
#' @importFrom dplyr mutate select all_of left_join group_by group_modify
#'     summarize filter
#' @importFrom SparseArray nnawhich
#' @importFrom purrr reduce
addFootprintsNomeR <- function(se,
                               assayName = "mod_prob",
                               fpLengths,
                               fpGroups,
                               threshUnmod = 0.5,
                               threshMod = 0.5,
                               minProb = 0.5,
                               BPPARAM = MulticoreParam(4L),
                               verbose = FALSE) {
    # TODO: verbose is currently not used
    # Check arguments
    .assertVector(x = se, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertVector(x = fpLengths, type = "list")
    .assertVector(x = fpGroups, type = "matrix")
    .assertScalar(x = threshUnmod, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = threshMod, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = minProb, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")
    .assertPackagesAvailable("nomeR")

    fpLengths$FTP_SUGGEST = fpGroups
    fpModels <- nomeR::get_ftp_models_for_prediction(fpLengths)

    # create a list of binary vectors (one vector per read), which will be
    # used as the input for predict_footprints
    # Q: Should these be in base space already (rather than observed base space)
    # Q: Does it matter if the reads are aligned to start in the same place or not (i.e., if there are leading NAs)
    binList <- lapply(setNames(colnames(se), colnames(se)), function(s) {
        tmp <- as.list(as.data.frame(as.matrix(assay(se, assayName)[[1]])))
        lapply(tmp, function(bl) {
            bl[bl >= threshMod] <- 1
            bl[bl < threshUnmod] <- 0
            bl[bl < threshMod & bl >= threshUnmod] <- NA
            1 - bl
        })
    })

    # unlist on the sample level to be able to run the footprint
    # prediction once
    binList <- unlist(binList, recursive = FALSE, use.names = TRUE)

    # predict footprints -> returns two data.frames (COVER_PROB and START_PROB)
    fpPred <- nomeR::predict_footprints(
        data = binList,
        footprint_models = fpModels$FTP_MODELS,
        bgprotectprob = fpModels$bgprotectprob,
        bgcoverprior = fpModels$bgcoverprior,
        ncpu = bpnworkers(BPPARAM)
    )

    # keep only COVER_PROB, add information about sample
    # Q: Should we use both START_PROB and COVER_PROB to determine the boundaries of the footprints? The approach below seems overly simplistic
    fpPred <- fpPred$COVER_PROB |>
        mutate(sample = sub("(.+)\\.(.+)", "\\1", seq)) |>
        mutate(seq = sub("(.+)\\.(.+)", "\\2", seq)) |>
        dplyr::select(-background)
    stopifnot(fpPred$sample %in% colnames(se))

    # merge cover probabilities for all lengths within each footprint group
    fpPred <- purrr::reduce(rownames(fpGroups), function(data, pat) {
        colnm <- grep(pat, names(data), value = TRUE)
        data |>
            mutate(!!pat := rowSums(across(all_of(colnm)))) |>
            select(-all_of(colnm))
    }, .init = fpPred)
    fpPred$pos <- start(se)[fpPred$pos]

    # interpolate/expand predicted probabilities to entire range
    .interpolate <- function(df, cols) {
        tmp <- as.matrix(df |> dplyr::select(all_of(cols)))
        pos <- df$pos
        .interpolateColumns(tmp, pos) |>
            as.data.frame() |>
            mutate(pos = seq(min(pos), max(pos))) |>
            left_join(df |> dplyr::select(-all_of(cols)),
                      by = "pos")
    }
    fpPred <- fpPred |>
        group_by(sample, seq) |>
        group_modify(~ .interpolate(.x, rownames(fpGroups)))

    # add information about the coverage range for each read, in order to
    # filter out regions not covered by a read before reporting footprints
    covRangeL <- lapply(assay(se, assayName), function(adat) {
        idx <- nnawhich(adat, arr.ind = TRUE)
        data.frame(
            seq = colnames(adat)[idx[, 2]],
            pos = start(se)[idx[, 1]]
        )
    })
    covRanges <- do.call(bind_rows, lapply(names(covRangeL), function(nm) {
        covRangeL[[nm]] |>
            group_by(seq) |>
            summarize(min = min(pos),
                      max = max(pos)) |>
            mutate(sample = nm)
    }))

    # for each footprint group, threshold probabilities and add to colData(se)
    cd <- colData(se)
    for (fpp in rownames(fpGroups)) {
        irlL <- lapply(setNames(colnames(se), colnames(se)), function(s) {
            # subset to positions with at least a certain probaility (minProb)
            # of being covered by the footprint group of interest
            aggrSub <- fpPred |>
                dplyr::select(all_of(c("sample", "seq", "pos", fpp))) |>
                dplyr::filter(.data[[fpp]] >= minProb,
                              sample == s) |>
                # restrict to range of read
                left_join(covRanges, by = c("seq", "sample")) |>
                dplyr::filter(pos >= min & pos <= max)

            # create IRangesList and subset to regions with length within
            # the indicated ones for the current footprint group
            do.call(
                IRangesList, lapply(
                    split(aggrSub, aggrSub$seq),
                    function(l) {
                        r <- IRanges::reduce(IRanges(start = l$pos,
                                                     end = l$pos))
                        r[width(r) >= fpGroups[fpp, 1] &
                              width(r) <= fpGroups[fpp, 2]]
                    }
                )
            )
        })
        ## add to se
        cd[[fpp]] <- irlL
        metadata(se)$readLevelData$colDataColumns <- union(
            metadata(se)$readLevelData$colDataColumns, fpp
        )
    }
    colData(se) <- cd

    se
}
