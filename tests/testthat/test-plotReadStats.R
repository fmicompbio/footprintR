test_that("plotReadStats works", {
    # example data
    bamf <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    se <- readModBam(bamf, regions = "chr1:6940000-6955000", modbase = "a",
                     BPPARAM = BiocParallel::SerialParam())
    se <- addReadStats(se, BPPARAM = BiocParallel::SerialParam())

    setmp <- se
    names(setmp$QC) <- "wrong_name"
    expect_error(plotReadStats(setmp),
                 "names of .se\\$readInfo. and .se\\$QC. are not identical")
    rm(setmp)

    gg <- plotReadStats(se)
    expect_s3_class(gg, "ggplot")

    gg <- plotReadStats(se, readInfoCol = NULL)
    expect_s3_class(gg, "ggplot")

    gg <- plotReadStats(se, qcCol = NULL)
    expect_s3_class(gg, "ggplot")
})
