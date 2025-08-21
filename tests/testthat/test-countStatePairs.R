test_that("countStatePairs works", {
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")

    expect_error(countStatePairs(),
                 ".bamfile. is missing")
    expect_error(countStatePairs(bamfile = "error"),
                 ".bamfile. does not exist")
    expect_error(countStatePairs(bamfile = modbamfile),
                 ".modbase. is missing")
    expect_error(countStatePairs(bamfile = modbamfile, modbase = "X"),
                 "invalid .modbase.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = 1L,
                                 modbase = "a"),
                 "must be of class .character.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshUnmod = -1),
                 "must be between 0 and 1")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshMod = -1),
                 "must be between 0 and 1")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", threshUnmod = 0.8, threshMod = 0.3),
                 "must be less than or equal to")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", windowSize = "error"),
                 "must be of class .numeric.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", minMapQ = "error"),
                 "must be of class .numeric.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", minAlignedLength = -1),
                 "must be between 0 and Inf")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", BPPARAM = "error"),
                 "must be of class .BiocParallelParam.")
    expect_error(countStatePairs(bamfile = modbamfile, regions = ".",
                                 modbase = "a", verbose = "error"),
                 "must be of class .logical.")

    res1 <- countStatePairs(bamfile = modbamfile,
                            regions = GenomicRanges::GRanges("chr1:1-100000000"),
                            modbase = "a", windowSize = 300)
    res2 <- countStatePairs(bamfile = modbamfile, regions = ".",
                            modbase = "a", windowSize = 300)
    res3 <- countStatePairs(bamfile = modbamfile, regions = ".",
                            modbase = "h", windowSize = 300)

    resL <- list(res1, res2, res3)
    for (i in seq_along(resL)) {
        expect_s4_class(resL[[i]], "DataFrame")
        expect_identical(colnames(resL[[i]]), c("S", "unmod_unmod", "unmod_mod",
                                                "mod_unmod", "mod_mod"))
        expect_identical(dim(resL[[i]]), c(300L, 5L))
    }

    expect_identical(res1, res2)
    expect_identical(sum(as.matrix(res3[, -1])), 0)
})
