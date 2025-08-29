test_that("plotReadStats works", {
    # example data
    bamf <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    se <- readModBam(bamf, regions = "chr1:6940000-6955000", modbase = "a",
                     BPPARAM = BiocParallel::SerialParam())
    expect_warning(
        se <- addReadStats(se, BPPARAM = BiocParallel::SerialParam()),
        "Too few points")

    setmp <- se
    names(setmp$QC) <- "wrong_name"
    expect_error(plotReadStats(setmp),
                 "names of .se\\$readInfo. and .se\\$QC. are not identical")
    rm(setmp)

    gg <- plotReadStats(se)
    expect_true(ggplot2::is_ggplot(gg))

    gg <- plotReadStats(se, readInfoCol = NULL)
    expect_true(ggplot2::is_ggplot(gg))

    gg <- plotReadStats(se, qcCol = NULL)
    expect_true(ggplot2::is_ggplot(gg))
    
    #SE is not normally part of the defaults:
    seSE <- addReadStats(se,stats="SEntrModProb", BPPARAM = BiocParallel::SerialParam())
    gg <- plotReadStats(se, qcCol = NULL)
    expect_true(ggplot2::is_ggplot(gg))
})
