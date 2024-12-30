test_that("filterReadsBam works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))

    # non-existing input files
    expect_error(filterReadsBam(infiles = filtbamfiles, outfiles = filtbamfiles, modbase = "a"))

    # non-unique input file names
    tmp <- modbamfiles
    names(tmp) <- c("s1", "s1")
    expect_error(filterReadsBam(infiles = tmp, outfiles = filtbamfiles, modbase = "a"))
    rm(tmp)

    # expected results (filtering out exactly one read for each filter)
    suppressMessages(
        expect_message(
            res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                                  modbase = "a", indexOutfiles = TRUE, minReadLength = 6746,
                                  minAlignedLength = 6896, minAlignedFraction = 0.56,
                                  minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                                  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
        )
    )
    expect_true(all(file.exists(filtbamfiles)))
    expect_true(all(file.exists(paste0(filtbamfiles, ".bai"))))
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

    # pre-existing output files
    expect_error(filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles, modbase = "a"))
    unlink(filtbamfiles)
    unlink(paste0(filtbamfiles, ".bai"))

    # expected results (using default parameters that deactivates all filters)
    res2 <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                           modbase = "a", indexOutfiles = FALSE,
                           BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                           verbose = FALSE)
    expect_true(all(file.exists(filtbamfiles)))
    expect_true(all(!file.exists(paste0(filtbamfiles, ".bai"))))
    expect_s3_class(res2, "data.frame")
    expect_identical(dim(res2), c(2L, 11L))
    expect_identical(colnames(res2), colnames(res))
    expect_identical(res2$sample, c("s1", "s2"))
    expect_identical(res2$infile, modbamfiles)
    expect_identical(res2$outfile, filtbamfiles)
    expect_identical(res2$total, c(10, 10))
    expect_identical(res2$retained, c(10, 10))
    expect_identical(res2$filtered_minReadLength, c(0, 0))
    expect_identical(res2$filtered_minAlignedLength, c(0, 0))
    expect_identical(res2$filtered_minAlignedFraction, c(0, 0))
    expect_identical(res2$filtered_minQscore, c(0, 0))
    expect_identical(res2$filtered_maxFracLowConf, c(0, 0))
    expect_identical(res2$filtered_maxEntropy, c(0, 0))
    expect_identical(unname(tools::md5sum(modbamfiles)),
                     unname(tools::md5sum(filtbamfiles)))
    unlink(filtbamfiles)
})
