test_that("plotReadStats works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                          package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a")
    se <- addReadStats(se)

    gg <- plotReadStats(se)
    expect_s3_class(gg, "ggplot")
})
