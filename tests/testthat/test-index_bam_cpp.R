test_that("index_bam_cpp works", {
    # copy bam files (without index)
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    outbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))
    expect_true(all(file.copy(from = modbamfiles, to = outbamfiles)))
    testfile <- tempfile(fileext = ".bam")

    # fail when indexing a non-existing file
    expect_error(index_bam_cpp(infile = testfile))

    # fail when indexing an empty file
    writeLines(text = character(0), con = testfile)
    expect_error(index_bam_cpp(infile = testfile))

    # fail when indexing a text file
    writeLines(text = LETTERS, con = testfile, sep = "\n")
    expect_error(index_bam_cpp(infile = testfile))

    # succeed when indexing bam files
    res <- rep(NA, length(outbamfiles))
    for (i in seq_along(outbamfiles)) {
        res[i] <- index_bam_cpp(infile = outbamfiles[i])
    }
    expect_length(res, length(outbamfiles))
    expect_identical(res, paste0(outbamfiles, ".bai"))
    expect_true(all(file.exists(paste0(outbamfiles, ".bai"))))
    expect_s4_class(readModBam(outbamfiles, regions = "chr1:6930000-6931000",
                               modbase = "a", level = "summary",
                               BPPARAM = BiocParallel::SerialParam()),
                    "RangedSummarizedExperiment")

    # clean up
    unlink(testfile)
    unlink(outbamfiles)
    unlink(paste0(outbamfiles, ".bai"))
})
