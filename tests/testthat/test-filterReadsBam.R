test_that("filterReadsBam works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))

    expect_error(filterReadsBam(infiles = filtbamfiles, outfiles = filtbamfiles, modbase = "a"))
    tmp <- modbamfiles
    names(tmp) <- c("s1", "s1")
    expect_error(filterReadsBam(infiles = tmp, outfiles = filtbamfiles, modbase = "a"))
    rm(tmp)
    res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                          modbase = "a", minReadLength = 6746,
                          minAlignedLength = 6896, minAlignedFraction = 0.56,
                          minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                          BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
    expect_s3_class(res, "data.frame")
    expect_identical(dim(res), c(2L, 11L))
    expect_identical(colnames(res), c("sample", "infile", "outfile", "total",
                                      "retained", "filtered_minReadLength",
                                      "filtered_minAlignedLength",
                                      "filtered_minAlignedFraction",
                                      "filtered_minQscore",
                                      "filtered_maxFracLowConf",
                                      "filtered_maxEntropy"))
    expect_identical(res$sample, c("s1", "s2"))
    expect_identical(res$infile, modbamfiles)
    expect_identical(res$outfile, filtbamfiles)
    expect_identical(res$total, c(10, 10))
    expect_identical(res$retained, c(7, 7))
    expect_identical(res$filtered_minReadLength, c(0, 1))
    expect_identical(res$filtered_minAlignedLength, c(1, 0))
    expect_identical(res$filtered_minAlignedFraction, c(1, 0))
    expect_identical(res$filtered_minQscore, c(0, 1))
    expect_identical(res$filtered_maxFracLowConf, c(0, 1))
    expect_identical(res$filtered_maxEntropy, c(1, 0))
    expect_error(filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles, modbase = "a"))
    unlink(filtbamfiles)
})
