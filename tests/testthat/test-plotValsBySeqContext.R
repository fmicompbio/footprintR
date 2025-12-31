test_that("plotValsBySeqContext works", {
    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")
    ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
                                                    package = "footprintR"))
    se <- readModBam(bamfiles = modbamfile, regions = "chr1", modbase = "a")
    se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = ref)
    se <- filterPositions(se, filters = "sequenceContext", sequenceContext = "NAN")
    se <- flattenReadLevelAssay(se, statistics = "Pmod")

    # fails with wrong arguments
    expect_error(plotValsBySeqContext(se = "error"),
                 ".se. must be of class .SummarizedExperiment.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = 1),
                 ".seqContextColumn. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "error"),
                 ".seqContextColumn. must be one of")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = 1),
                 ".assayName. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "error"),
                 ".assayName. must be one of")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = 1),
                 ".plotType. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "error"),
                 ".plotType. must be one of")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = "x"),
                 ".topN. must be of class .numeric.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = c(1, 2)),
                 ".topN. must have length 1")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = -1),
                 ".topN. must be between 0 and Inf")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = "x"),
                 ".bottomN. must be of class .numeric.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = c(1, 2)),
                 ".bottomN. must have length 1")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = -1),
                 ".bottomN. must be between 0 and Inf")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = 2, flipCoord = 1),
                 ".flipCoord. must be of class .logical.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = 2, flipCoord = c(TRUE, FALSE)),
                 ".flipCoord. must have length 1")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = 2, flipCoord = TRUE,
                                      yAxisLabel = 1),
                 ".yAxisLabel. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "violin",
                                      topN = 2, bottomN = 2, flipCoord = TRUE,
                                      yAxisLabel = c("one", "two")),
                 ".yAxisLabel. must have length 1")

    # works with correct arguments
    g1 <- plotValsBySeqContext(se = se, assayName = "mod_prob", plotType = "violin",
                               topN = 2, bottomN = 2, flipCoord = TRUE)
    g2 <- plotValsBySeqContext(se = se, assayName = "Pmod", plotType = "violin",
                               topN = 2, bottomN = 2, flipCoord = TRUE)
    expect_true(ggplot2::is_ggplot(g1))
    expect_true(ggplot2::is_ggplot(g2))
    expect_identical(g1@data, g2@data)

    g3 <- plotValsBySeqContext(se = se, assayName = "Pmod", plotType = "bar",
                               topN = 2, bottomN = 2, flipCoord = TRUE)
    expect_true(ggplot2::is_ggplot(g3))
    expect_s3_class(g3@data, "tbl_df")
    expect_equal(dim(g3@data), c(4L, 3L))

    g4 <- plotValsBySeqContext(se = se, assayName = "Pmod", plotType = "errorbar",
                               topN = 2, bottomN = 2, flipCoord = FALSE)
    expect_true(ggplot2::is_ggplot(g4))
    expect_s3_class(g4@data, "tbl_df")
    expect_equal(dim(g4@data), c(4L, 3L))
})
