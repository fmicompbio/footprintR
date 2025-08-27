test_that("filterReadsBam works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))

    # non-existing input files
    expect_error(filterReadsBam(infiles = filtbamfiles,
                                outfiles = filtbamfiles,
                                modbase = "a"),
                 "not all .infiles. exist")

    # non-existing bam index
    tmpin <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = modbamfiles[1], to = tmpin))
    expect_error(filterReadsBam(infiles = tmpin, outfiles = filtbamfiles[1],
                                modbase = "a",
                                BPPARAM = BiocParallel::SerialParam()),
                 "Failed to load the index")
    unlink(c(tmpin, filtbamfiles[1]))

    # direct call to filter_modbam_cpp with verbose = TRUE
    suppressMessages({
        expect_length(filter_modbam_cpp(infile = modbamfiles[1],
                                        outfile = filtbamfiles[1],
                                        modbase = "a", region = ".",
                                        includeHeader = TRUE,
                                        verbose = TRUE),
                      12L)
    })
    unlink(filtbamfiles[1])

    # filter_modbam_cpp with unknown output extension
    tmpsam <- tempfile(fileext = ".error")
    expect_error(filter_modbam_cpp(infile = modbamfiles[1],
                                   outfile = tmpsam,
                                   modbase = "a", region = ".",
                                   includeHeader = TRUE,
                                   verbose = FALSE),
                 "Unknown .outfile. extension")

    # creating sam output from filter_modbam_cpp
    tmpsam <- tempfile(fileext = ".sam")
    res <- filter_modbam_cpp(infile = modbamfiles[1],
                             outfile = tmpsam,
                             modbase = "a", region = ".",
                             includeHeader = FALSE,
                             verbose = FALSE)
    expect_equal(res[["retained"]], length(readLines(tmpsam)))
    unlink(tmpsam)

    # miss-specified region
    expect_error(filter_modbam_cpp(infile = modbamfiles[1], outfile = filtbamfiles[1], modbase = "a", region = "ERROR"))
    unlink(filtbamfiles[1])

    # expected results (filtering out exactly one read for each filter)
    suppressMessages(
        expect_message(
            res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                                  modbase = "a", indexOutfiles = TRUE, minReadLength = 6746,
                                  minAlignedLength = 6896, minAlignedFraction = 0.56, minSNR=-0.768,
                                  minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                                  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
        )
    )
    expect_true(all(file.exists(filtbamfiles)))
    expect_true(all(file.exists(paste0(filtbamfiles, ".bai"))))
    expect_s3_class(res, "data.frame")
    expect_identical(dim(res), c(2L, 15L))
    expect_identical(colnames(res), c("sample", "infile", "outfile", "total",
                                      "retained", "filtered_unmapped",
                                      "filtered_secondary", "filtered_supplementary",
                                      "filtered_minReadLength",
                                      "filtered_minAlignedLength",
                                      "filtered_minAlignedFraction",
                                      "filtered_minQscore",
                                      "filtered_minSNR",
                                      "filtered_maxFracLowConf",
                                      "filtered_maxEntropy"))
    expect_identical(res$sample, c("s1", "s2"))
    expect_identical(res$infile, modbamfiles)
    expect_identical(res$outfile, filtbamfiles)
    expect_identical(res$total, c(10, 10))
    expect_identical(res$retained, c(6, 7))
    expect_identical(res$filtered_unmapped, c(0, 0))
    expect_identical(res$filtered_secondary, c(0, 0))
    expect_identical(res$filtered_supplementary, c(0, 0))
    expect_identical(res$filtered_minReadLength, c(0, 1))
    expect_identical(res$filtered_minAlignedLength, c(1, 0))
    expect_identical(res$filtered_minAlignedFraction, c(1, 0))
    expect_identical(res$filtered_minQscore, c(0, 1))
    expect_identical(res$filtered_minSNR, c(1, 0))
    expect_identical(res$filtered_maxFracLowConf, c(0, 1))
    expect_identical(res$filtered_maxEntropy, c(1, 0))

    # pre-existing output files
    expect_error(filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles, modbase = "a"))

    # pre-existing output files (overwriteOutfiles = TRUE)
    res0 <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                           modbase = "a", indexOutfiles = FALSE,
                           overwriteOutfiles = TRUE, minReadLength = 6746,
                           minAlignedLength = 6896, minAlignedFraction = 0.56, minSNR=-0.768,
                           minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                           BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
    expect_identical(res, res0)
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
    expect_identical(dim(res2), c(2L, 15L))
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
    expect_identical(res2$filtered_minSNR, c(0, 0))
    expect_identical(res2$filtered_maxFracLowConf, c(0, 0))
    expect_identical(res2$filtered_maxEntropy, c(0, 0))
    expect_identical(unname(tools::md5sum(modbamfiles)),
                     unname(tools::md5sum(filtbamfiles)))
    unlink(filtbamfiles)

    
    
    
    
    
    # non-primary alignments
    inbam <- system.file("extdata", "6mA_nonPrimary.bam", package = "footprintR")
    outbam <- tempfile(fileext = ".bam")
    # ... keeping them
    res3 <- filterReadsBam(infiles = inbam, outfiles = outbam,
                           modbase = "a", indexOutfiles = FALSE,
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
    expect_identical(res3[, c("total", "retained")], data.frame(total = 3, retained = 3))
    tmp3 <- Rsamtools::scanBam(file = outbam)
    expect_length(tmp3[[1]]$qname, 3L)
    unlink(outbam)
    # ... dropping them
    res4 <- filterReadsBam(infiles = inbam, outfiles = outbam,
                           modbase = "a", indexOutfiles = FALSE,
                           keepUnmapped = FALSE, keepSecondary = FALSE, keepSupplementary = FALSE,
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
    expect_identical(res4[, c("total", "retained", "filtered_unmapped",
                              "filtered_secondary", "filtered_supplementary")],
                     data.frame(total = 3, retained = 0, filtered_unmapped = 1,
                                filtered_secondary = 1, filtered_supplementary = 1))
    tmp4 <- Rsamtools::scanBam(file = outbam)
    expect_length(tmp4[[1]]$qname, 0L)
    unlink(outbam)
    
    
    # expected results (passing precalculated noiseCoefs that result in stricter filtering)
    suppressMessages(
        expect_message(
            res5 <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                                  modbase = "a", indexOutfiles = TRUE, minReadLength = 6746,
                                  minAlignedLength = 6896, minAlignedFraction = 0.56, minSNR=-0.768, noiseCoef=c(0,0.29),
                                  minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                                  BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
        )
    )
    expect_true(all(file.exists(filtbamfiles)))
    expect_true(all(file.exists(paste0(filtbamfiles, ".bai"))))
    expect_s3_class(res5, "data.frame")
    expect_identical(dim(res5), c(2L, 15L))
    expect_identical(res5$total, c(10, 10))
    expect_identical(res5$retained, c(5, 7))
    expect_identical(res5$filtered_minSNR, c(2, 1))
})
