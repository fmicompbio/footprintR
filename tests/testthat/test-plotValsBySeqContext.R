test_that("plotValsBySeqContext works", {
    # example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
                                                    package = "footprintR"))
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1", modbase = "a",
                     BPPARAM = BiocParallel::SerialParam())
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
                                      assayName = "mod_prob", aggregation = 1),
                 ".aggregation. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", aggregation = "error"),
                 ".aggregation. must be one of")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", facetBy = 1),
                 ".facetBy. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", facetBy = "error"),
                 ".facetBy. must be one of")
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
    g1 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "violin",
                               topN = Inf, bottomN = Inf, flipCoord = TRUE,
                               aggregation = "mean", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g1))
    expect_s3_class(g1@data, "data.frame")
    expect_identical(dim(g1@data), c(7975L + 6949L, 4L))

    g2 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "violin",
                               topN = Inf, bottomN = Inf, flipCoord = TRUE,
                               aggregation = "mean", facetBy = NULL)
    expect_true(ggplot2::is_ggplot(g2))
    expect_s3_class(g2@data, "data.frame")
    expect_identical(dim(g2@data), c(8112L, 4L))

    # check that subsetting to first sample gives identical values to
    # facetting within function
    g3a <- plotValsBySeqContext(se = se[, 1], assayName = "mod_prob",
                                plotType = "bar", aggregation = "mean",
                                facetBy = NULL, topN = 2, bottomN = 2,
                                flipCoord = TRUE)
    expect_true(ggplot2::is_ggplot(g3a))
    expect_s3_class(g3a@data, "tbl_df")
    expect_equal(dim(g3a@data), c(4L, 5L))
    expect_identical(as.character(g3a@data$seqContext), c("GAG___1", "CAG___1", "TAA___1", "TAT___1"))
    expect_equal(g3a@data$valsMean, c(0.1656683, 0.151384387, 0.051409203, 0.062818268))

    g3b <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                                plotType = "bar", aggregation = "mean",
                                facetBy = "sample", topN = 2, bottomN = 2,
                                flipCoord = TRUE)
    expect_true(ggplot2::is_ggplot(g3b))
    expect_s3_class(g3b@data, "tbl_df")
    expect_equal(dim(g3b@data), c(8L, 5L))
    tmp1 <- g3a@data
    levels(tmp1$seqContext) <- sub("1", "s1", levels(tmp1$seqContext))
    tmp <- g3b@data |> dplyr::filter(sample == "s1") |>
        dplyr::left_join(tmp1, by = "seqContext")
    expect_identical(tmp$valsMean.x, tmp$valsMean.y)
    expect_identical(tmp$valsSd.x, tmp$valsSd.y)
    expect_identical(g3b@data |> dplyr::filter(seqContext == "GAG___s1") |> dplyr::pull(valsMean), mean(SummarizedExperiment::assay(se, "Pmod")[SummarizedExperiment::rowData(se)$sequenceContext == "GAG", "s1"], na.rm = TRUE))
    expect_identical(g3b@data |> dplyr::filter(seqContext == "TAA___s2") |> dplyr::pull(valsSd), sd(SummarizedExperiment::assay(se, "Pmod")[SummarizedExperiment::rowData(se)$sequenceContext == "TAA", "s2"], na.rm = TRUE))

    g4 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "errorbar", aggregation = "mean",
                               topN = 2, bottomN = 2, flipCoord = FALSE)
    expect_true(ggplot2::is_ggplot(g4))
    expect_s3_class(g4@data, "tbl_df")
    expect_equal(dim(g4@data), c(4L, 5L))

    g5 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "errorbar", aggregation = "none",
                               topN = 2, bottomN = 2, flipCoord = FALSE)
    expect_true(ggplot2::is_ggplot(g5))
    expect_s3_class(g5@data, "tbl_df")
    expect_equal(dim(g5@data), c(4L, 5L))

    g6 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "errorbar", aggregation = "none",
                               topN = 2, bottomN = 2, flipCoord = FALSE,
                               facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g6))
    expect_s3_class(g6@data, "tbl_df")
    expect_equal(dim(g6@data), c(8L, 5L))

    g7 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "errorbar", aggregation = "mean",
                               topN = 2, bottomN = 2, flipCoord = FALSE,
                               facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g7))
    expect_s3_class(g7@data, "tbl_df")
    expect_equal(dim(g7@data), c(8L, 5L))
})
