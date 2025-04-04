test_that(".filterScores works", {
    skip_if_not_installed("signal")

    xs <- pi * seq(0, 10, length.out = 100)
    ys <- sin(xs) + sin(xs / 3) + sin(xs * 10)

    # band-pass filter
    res1 <- .filterScores(score = ys, minperiod = 2, maxperiod = 10, type = "pass")
    expect_vector(res1, numeric())
    expect_type(res1, "double")
    expect_length(res1, length(xs))
    expect_true(cor(res1, sin(xs)) > cor(res1, ys))
    # plot(xs, ys, type = "l")
    # lines(xs, res1, col = "red")

    # low-pass filter
    res2 <- .filterScores(score = ys, minperiod = 2, type = "low")
    expect_vector(res2, numeric())
    expect_type(res2, "double")
    expect_length(res2, length(xs))
    expect_true(cor(res2, sin(xs) + sin(xs / 3)) > cor(res2, ys))
    # plot(xs, ys, type = "l")
    # lines(xs, res2, col = "red")

    # high-pass filter
    res3 <- .filterScores(score = ys, maxperiod = 10, type = "high")
    expect_vector(res3, numeric())
    expect_type(res3, "double")
    expect_length(res3, length(xs))
    expect_true(cor(res3, sin(xs * 10)) > cor(res2, ys))
    # plot(xs, ys, type = "l")
    # lines(xs, res3, col = "red")
})

test_that("addFootprint and helper functions work", {
    skip_if_not_installed("signal")

    # prepare data
    wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")
    se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
                     modbase = "a", verbose = FALSE,
                     BPPARAM = BiocParallel::SerialParam())
    adat <- assay(se, "mod_prob")
    SparseArray::nnavals(adat$s1) <- NA
    seEmpty <- SummarizedExperiment(assays = list(mod_prob = adat),
                                    metadata = metadata(se),
                                    colData = colData(se))

    # arguments
    expect_error(calcFootprintScores(se = "error"))
    expect_error(calcFootprintScores(se = se, wgt = "error"))
    expect_error(calcFootprintScores(se = se, wgt = wgt, assayName = "error"))
    expect_error(calcFootprintScores(se = se, wgt = wgt, minconf = "error"))
    expect_error(calcFootprintScores(se = se, wgt = wgt, minweight = "error"))
    expect_error(calcFootprintScores(se = se, wgt = wgt, verbose = "error"))
    expect_error(segmentFootprintScores(scoresList = "error"))
    expect_error(segmentFootprintScores(scoresList = list(data.frame(error = 1))))
    expect_error(segmentFootprintScores(scoresList = list(), minperiod = "error"))
    expect_error(segmentFootprintScores(scoresList = list(), maxperiod = "error"))
    expect_error(segmentFootprintScores(scoresList = list(), thresh = "error"))
    expect_error(segmentFootprintScores(scoresList = list(), lenRange = "error"))
    expect_error(segmentFootprintScores(scoresList = list(), lenRange = 1))
    expect_error(segmentFootprintScores(scoresList = list(), lenRange = c(-3, -1)))
    expect_error(segmentFootprintScores(scoresList = list(), width = "error"))
    expect_error(segmentFootprintScores(scoresList = list(), verbose = "error"))
    expect_error(addFootprints(se = se, name = 1:10))
    expect_error(addFootprints(se = se, verbose = "error"))

    # expected results
    # ... calcFootprintScores
    scoresListEmpty <- calcFootprintScores(seEmpty, wgt, verbose = FALSE)
    expect_vector(scoresListEmpty, list())
    expect_length(scoresListEmpty, ncol(se))
    expect_identical(nrow(scoresListEmpty$s1), 0L)
    suppressMessages(expect_message(
        scoresList <- calcFootprintScores(se, wgt)
    ))
    expect_vector(scoresList, list())
    expect_length(scoresList, ncol(se))
    for (i in seq_along(scoresList)) {
        expect_vector(scoresList[[i]], data.frame(readId = factor(levels = rownames(se$readInfo$s1)),
                                                  pos = integer(),
                                                  pmod = numeric(),
                                                  score = numeric()))
    }
    expect_identical(lapply(scoresList, nrow), list(s1 = 38549L))
    expect_identical(colnames(scoresList$s1), c("readId", "pos", "pmod", "score"))
    expect_identical(unique(as.character(scoresList$s1$readId)), rownames(se$readInfo$s1))
    expect_identical(range(scoresList$s1$pos), range(pos(rowRanges(se))))
    expect_equal(sum(scoresList$s1$score, na.rm = TRUE), 21.599059746360531165)

    # ... segmentFootprintScores
    irlEmpty <- segmentFootprintScores(scoresListEmpty, verbose = FALSE)
    expect_vector(irlEmpty, list())
    irlEmpty2 <- segmentFootprintScores(scoresList, thresh = 100, verbose = FALSE)
    expect_vector(irlEmpty2, list())
    expect_identical(sum(unlist(lapply(irlEmpty, lengths))),
                     sum(unlist(lapply(irlEmpty2, lengths))))
    suppressMessages(expect_message(
        irl <- segmentFootprintScores(scoresList, thresh = 0.05, width = 140)
    ))
    expect_vector(irl, list())
    expect_length(irl, ncol(se))
    expect_identical(names(irl), colnames(se))
    for (i in seq_along(irl)) {
        expect_s4_class(irl[[i]], "IRangesList")
        expect_identical(names(irl[[1]]), rownames(se$readInfo[[i]]))
    }
    expect_true(all(width(unlist(unname(irl$s1))) == 140L))
    expect_identical(lengths(irl$s1),
                     structure(c(15L, 1L, 2L),
                               names = rownames(se$readInfo$s1)))

    # ... addFootprints
    se <- addFootprints(se, wgt, thresh = 0.05, name = "nucl", verbose = FALSE)
    expect_true("nucl" %in% colnames(colData(se)))
    expect_identical(se$nucl, irl)
})
