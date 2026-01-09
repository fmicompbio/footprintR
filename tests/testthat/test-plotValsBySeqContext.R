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
                                      assayName = "mod_prob", selectContextsBy = 1),
                 ".selectContextsBy. must be of class .character.")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", selectContextsBy = "error"),
                 ".selectContextsBy. must be one of")
    expect_error(plotValsBySeqContext(se = se, seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", selectContextsBy = "sample",
                                      facetBy = NULL),
                 ".selectContextsBy. must be one of")
    expect_error(plotValsBySeqContext(se = se[, 1], seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", selectContextsBy = "sample_var"),
                 "must have at least two samples")
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
    expect_error(plotValsBySeqContext(se = se[, 1], seqContextColumn = "sequenceContext",
                                      assayName = "mod_prob", plotType = "pairs"),
                 "must have at least two samples for a pairs plot")
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
    # ... violin, all contexts, select by sample, mean aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "sample",
                              topN = Inf, bottomN = Inf, flipCoord = TRUE,
                              aggregation = "mean", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all rows with observed data should be part of the plot
    expect_identical(dim(g@data), c(7975L + 6949L, 7L))
    expect_identical(sum(g@data$seqContext == "GAT___s1"), 351L)
    expect_equal(unique(g@data$valsMean[g@data$seqContext == "GAT___s1"]), 0.14200813948)

    # ... violin, all contexts, select by sample, no aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "sample",
                              topN = Inf, bottomN = Inf, flipCoord = TRUE,
                              aggregation = "none", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all observed values should be part of the plot
    expect_identical(dim(g@data), c(29033L + 26095L, 7L))
    expect_equal(unique(g@data$valsMean[g@data$seqContext == "GAT___s1"]), 0.139679259)
    expect_equal(unique(g@data$valsMean[g@data$seqContext == "GAG___s1"]), 0.158524954)

    # ... violin, top/bottom context, select by sample, no aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "sample",
                              topN = 1, bottomN = 1, flipCoord = TRUE,
                              aggregation = "none", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    expect_equal(unique(g@data$valsMean[g@data$seqContext == "GAG___s1"]), 0.158524954)
    expect_equal(levels(g@data$seqContext), c("TAA___s1", "CAA___s2", "GAG___s1", "GAG___s2"))

    # ... violin, all contexts, select by overall, mean aggregation, no facet
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "overall",
                              topN = Inf, bottomN = Inf, flipCoord = TRUE,
                              aggregation = "mean", facetBy = NULL)
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all rows with observed data should be part of the plot
    expect_identical(dim(g@data), c(8112L, 6L))
    expect_identical(sum(g@data$seqContext == "GAT___overall"), 355L)
    expect_equal(mean(g@data$vals[g@data$seqContext == "GAT___overall"]), 0.160786660)

    # ... violin, top/bottom 2 contexts, select by sample_var, no aggregation, no facet
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "sample_var",
                              topN = 2, bottomN = 2, flipCoord = TRUE,
                              aggregation = "none", facetBy = NULL)
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all rows with observed data should be part of the plot
    expect_identical(dim(g@data), c(12884L, 7L))
    expect_identical(sum(g@data$seqContext == "AAG___overall"), 3871L)
    expect_equal(mean(g@data$vals[g@data$seqContext == "AAG___overall"]), 0.137676896)
    expect_identical(levels(g@data$seqContext),
                     c("TAG___overall", "CAA___overall", "GAG___overall", "AAG___overall"))

    # ... violin, all contexts, select by overall, mean aggregation, no facet,
    #     fill by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "violin", selectContextsBy = "overall",
                              topN = Inf, bottomN = Inf, flipCoord = TRUE,
                              aggregation = "mean", facetBy = NULL,
                              fillBy = "sample", fillColors = c("blue", "orange"))
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all rows with observed data should be part of the plot
    expect_identical(dim(g@data), c(7975L + 6949L, 6L))
    expect_identical(sum(g@data$seqContext == "GAT___overall" &
                             g@data$sample == "s1"), 351L)
    expect_equal(mean(g@data$vals[g@data$seqContext == "GAT___overall" &
                                      g@data$sample == "s1"]), 0.14200813948)
    # ... ... same but provide too few colors
    expect_warning(
        g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                                  plotType = "violin", selectContextsBy = "overall",
                                  topN = Inf, bottomN = Inf, flipCoord = TRUE,
                                  aggregation = "mean", facetBy = NULL,
                                  fillBy = "sample", fillColors = c("blue")),
        "Not enough colors")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "data.frame")
    # ... ... all rows with observed data should be part of the plot
    expect_identical(dim(g@data), c(7975L + 6949L, 6L))
    expect_identical(sum(g@data$seqContext == "GAT___overall" &
                             g@data$sample == "s1"), 351L)
    expect_equal(mean(g@data$vals[g@data$seqContext == "GAT___overall" &
                                      g@data$sample == "s1"]), 0.14200813948)

    # ... check that subsetting to first sample gives identical values to
    #     facetting within function
    g1 <- plotValsBySeqContext(se = se[, 1], assayName = "mod_prob",
                               plotType = "bar", aggregation = "mean",
                               facetBy = "sample", topN = 2, bottomN = 2,
                               flipCoord = TRUE, selectContextsBy = "sample")
    expect_true(ggplot2::is_ggplot(g1))
    expect_s3_class(g1@data, "tbl_df")
    expect_equal(dim(g1@data), c(4L, 7L))
    expect_identical(as.character(g1@data$seqContext), c("CAG___s1", "GAG___s1", "TAA___s1", "TAT___s1"))
    expect_equal(g1@data$valsMean, c(0.151384387, 0.1656683, 0.051409203, 0.062818268))

    g2 <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                               plotType = "bar", aggregation = "mean",
                               facetBy = "sample", topN = 2, bottomN = 2,
                               flipCoord = TRUE, selectContextsBy = "sample")
    expect_true(ggplot2::is_ggplot(g2))
    expect_s3_class(g2@data, "tbl_df")
    expect_equal(dim(g2@data), c(8L, 7L))
    tmp1 <- g1@data
    tmp <- g2@data |> dplyr::filter(sample == "s1") |>
        dplyr::left_join(tmp1, by = "seqContext")
    expect_identical(tmp$valsMean.x, tmp$valsMean.y)
    expect_identical(tmp$valsSd.x, tmp$valsSd.y)
    expect_identical(g2@data |> dplyr::filter(seqContext == "GAG___s1") |> dplyr::pull(valsMean), mean(SummarizedExperiment::assay(se, "Pmod")[SummarizedExperiment::rowData(se)$sequenceContext == "GAG", "s1"], na.rm = TRUE))
    expect_identical(g2@data |> dplyr::filter(seqContext == "TAA___s2") |> dplyr::pull(valsSd), sd(SummarizedExperiment::assay(se, "Pmod")[SummarizedExperiment::rowData(se)$sequenceContext == "TAA", "s2"], na.rm = TRUE))

    # ... errorbar, top/bottom 4 contexts, select by overall, mean aggregation, no facet
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "errorbar", aggregation = "mean",
                              topN = 4, bottomN = 4, flipCoord = FALSE,
                              selectContextsBy = "overall", facetBy = NULL)
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(8L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAT___overall"], 0.160786660)

    # ... errorbar, top/bottom 4 contexts, select by overall, mean aggregation,
    #     no facet, fill by sample
    expect_warning(
        g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                                  plotType = "errorbar", aggregation = "mean",
                                  topN = 4, bottomN = 4, flipCoord = FALSE,
                                  selectContextsBy = "overall", facetBy = NULL,
                                  fillBy = "sample"),
        "Not enough colors")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAT___overall" &
                                     g@data$sample == "s1"], 0.14200813948)
    # ... ... same but provide too few colors
    expect_warning(
        g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                                  plotType = "errorbar", aggregation = "mean",
                                  topN = 4, bottomN = 4, flipCoord = FALSE,
                                  selectContextsBy = "overall", facetBy = NULL,
                                  fillBy = "sample", fillColors = "red"),
        "Not enough colors")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAT___overall" &
                                     g@data$sample == "s1"], 0.14200813948)

    # ... errorbar, top/bottom 4 contexts, select by overall, no aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "errorbar", aggregation = "none",
                              topN = 4, bottomN = 4, flipCoord = FALSE,
                              selectContextsBy = "sample", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAG___s1"], 0.158524954)

    # ... bar, top/bottom 4 contexts, select by sample, no aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "bar", aggregation = "none",
                              topN = 4, bottomN = 4, flipCoord = FALSE,
                              selectContextsBy = "sample", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAG___s1"], 0.158524954)

    # ... bar, top/bottom 4 contexts, select by overall, no aggregation, facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "bar", aggregation = "none",
                              topN = 4, bottomN = 4, flipCoord = FALSE,
                              selectContextsBy = "overall", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAG___overall" & g@data$sample == "s1"], 0.158524954)

    # ... bar, top/bottom 4 contexts, select by overall, no aggregation, no facet
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "bar", aggregation = "none",
                              topN = 4, bottomN = 4, flipCoord = FALSE,
                              selectContextsBy = "overall", facetBy = NULL)
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(8L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAG___overall"], 0.180227711)

    # ... bar, top/bottom 4 contexts, select by overall, no aggregation, no facet,
    #     fill by sample
    expect_warning(
        g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                                  plotType = "bar", aggregation = "none",
                                  topN = 4, bottomN = 4, flipCoord = FALSE,
                                  selectContextsBy = "overall", facetBy = NULL,
                                  fillBy = "sample"),
        "Not enough colors")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(16L, 7L))
    expect_equal(g@data$valsMean[g@data$seqContext == "GAG___overall" &
                                     g@data$sample == "s1"], 0.158524954)

    # ... bar, top/bottom 2 contexts, select by sample_union, mean aggregation,
    #     facet by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              yAxisLabel = "New title",
                              plotType = "bar", aggregation = "mean",
                              topN = 2, bottomN = 2, flipCoord = FALSE,
                              selectContextsBy = "sample_union", facetBy = "sample")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(12L, 7L))
    expect_identical(levels(g@data$seqContext),
                     c("GAG___overall", "GAT___overall", "CAG___overall",
                       "TAT___overall", "CAA___overall", "TAA___overall"))

    # ... bar, top/bottom 2 contexts, select by sample_union, mean aggregation,
    #     facet by sample, change default color
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              yAxisLabel = "New title",
                              plotType = "bar", aggregation = "mean",
                              topN = 2, bottomN = 2, flipCoord = FALSE,
                              selectContextsBy = "sample_union",
                              facetBy = "sample", fillColors = "forestgreen")
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(12L, 7L))
    expect_identical(levels(g@data$seqContext),
                     c("GAG___overall", "GAT___overall", "CAG___overall",
                       "TAT___overall", "CAA___overall", "TAA___overall"))

    # ... bar, top/bottom 2 contexts, select by sample_union, mean aggregation,
    #     facet by sample, change default color by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              yAxisLabel = "New title",
                              plotType = "bar", aggregation = "mean",
                              topN = 2, bottomN = 2, flipCoord = FALSE,
                              selectContextsBy = "sample_union",
                              facetBy = "sample", fillBy = "sample",
                              fillColors = c("forestgreen", "orange"))
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(12L, 7L))
    expect_identical(levels(g@data$seqContext),
                     c("GAG___overall", "GAT___overall", "CAG___overall",
                       "TAT___overall", "CAA___overall", "TAA___overall"))

    # ... bar, top/bottom 2 contexts, select by sample_var, mean aggregation,
    #     fill by sample, change default color by sample
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              plotType = "bar", aggregation = "mean",
                              topN = 2, bottomN = 2, flipCoord = FALSE,
                              selectContextsBy = "sample_var",
                              facetBy = NULL, fillBy = "sample",
                              fillColors = c("forestgreen", "orange"))
    expect_true(ggplot2::is_ggplot(g))
    expect_s3_class(g@data, "tbl_df")
    expect_equal(dim(g@data), c(8L, 8L))
    expect_identical(levels(g@data$seqContext),
                     c("CAC___overall", "AAG___overall", "CAA___overall", "TAG___overall"))
    expect_equal(g@data$valsMeanVar[g@data$seqContext == "AAG___overall"][1], 0.001236437,
                 tolerance = 1e-6)
    expect_equal(g@data$valsMeanVar[g@data$seqContext == "TAG___overall"][1], 0.0001539375,
                 tolerance = 1e-6)

    # ... pairs, no selection
    g <- plotValsBySeqContext(se = se, assayName = "mod_prob",
                              yAxisLabel = "New title",
                              plotType = "pairs", aggregation = "mean",
                              topN = Inf, bottomN = 0, flipCoord = FALSE,
                              selectContextsBy = "sample_union",
                              facetBy = "sample", fillBy = "sample",
                              fillColors = c("forestgreen", "orange"))
    expect_s7_class(g, GGally::ggmatrix)
    expect_s3_class(g@data, "data.frame")
    expect_equal(dim(g@data), c(length(unique(rowData(se)$sequenceContext)), 3L))
    expect_equal(g@data[g@data$seqContext == "GAT", "s1"], 0.14200813948)

})
