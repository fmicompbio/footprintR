suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(GenomicRanges)
})

## -------------------------------------------------------------------------- ##
## Checks, readModBam
## -------------------------------------------------------------------------- ##
test_that("readModBam works", {
    # example data
    modbamfiles <- system.file("extdata",
                               c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    names(modbamfiles) <- c("sample1", "sample2")
    extractfiles <- system.file("extdata",
                                c("modkit_extract_rc_6mA_1.tsv.gz",
                                  "modkit_extract_rc_6mA_2.tsv.gz"),
                               package = "footprintR")
    names(extractfiles) <- names(modbamfiles)

    # invalid arguments
    expect_error(readModBam("error", "chr1:6940000-6955000", "a"),
                 "not all `bamfiles` exist")
    expect_error(readModBam(structure(unname(modbamfiles), names = c("s1", "s1")),
                            "chr1:6940000-6955000", "a"),
                 "are not unique")
    expect_error(readModBam(modbamfiles, "error", "a"),
                 "GRanges object must contain")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "Z"),
                 "invalid `modbase` values")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", c("a", "a", "a")),
                 "must have length")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000",
                            c(sample1 = "a", sample3 = "a")),
                 "names of `modbase` and `bamfiles` don't agree")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "a", "error"),
                 "`seqinfo` must be `NULL`, a `Seqinfo` object or")

    # expected results
    se0 <- readModkitExtract(fnames = extractfiles, modbase = "a")
    reg1 <- c("chr1:6940000-6955000", "chr1:6929000-6929500")
    reg2 <- GRanges("chr1", IRanges(start = 6940000, end = 6955000))
    reg3 <- rep("chr1:6940000-6955000", 3)
    reg4 <- "chr1:6940000-6955000"
    expect_message(expect_message(expect_message(expect_message(
        expect_message(expect_message(expect_message(expect_message(
            expect_message(
                se1 <- readModBam(bamfiles = modbamfiles,
                                  regions = reg1,
                                  modbase = "a", verbose = TRUE)
    )))))))))
    se2 <- readModBam(bamfiles = unname(modbamfiles),
                      regions = reg2,
                      modbase = "a",
                      verbose = FALSE)
    se3 <- readModBam(bamfiles = modbamfiles,
                      regions = reg3,
                      modbase = "a",
                      verbose = FALSE)
    se4 <- readModBam(bamfiles = modbamfiles,
                      regions = reg4,
                      modbase = c("a", "m"),
                      verbose = FALSE)

    # ... structure
    expect_s4_class(se1, "RangedSummarizedExperiment")
    expect_s4_class(se2, "RangedSummarizedExperiment")
    expect_s4_class(se3, "RangedSummarizedExperiment")
    expect_s4_class(se4, "RangedSummarizedExperiment")

    expect_s4_class(rowRanges(se1), "GPos")
    expect_s4_class(rowRanges(se2), "GPos")
    expect_s4_class(rowRanges(se3), "GPos")
    expect_s4_class(rowRanges(se4), "GPos")

    expected_coldata_names <- c("sample", "modbase", "n_reads", "qscore")
    expect_identical(colnames(colData(se1)), expected_coldata_names)
    expect_identical(colnames(colData(se2)), expected_coldata_names)
    expect_identical(colnames(colData(se3)), expected_coldata_names)
    expect_identical(colnames(colData(se4)), expected_coldata_names)

    expect_identical(se1$sample, names(modbamfiles))
    expect_identical(se2$sample, c("s1", "s2"))
    expect_identical(se3$sample, names(modbamfiles))
    expect_identical(se4$sample, names(modbamfiles))

    expect_identical(assayNames(se1), "mod_prob")
    expect_identical(assayNames(se2), "mod_prob")
    expect_identical(assayNames(se3), "mod_prob")
    expect_identical(assayNames(se4), "mod_prob")

    expect_s4_class(assay(se1, "mod_prob"), "DFrame")
    expect_s4_class(assay(se2, "mod_prob"), "DFrame")
    expect_s4_class(assay(se3, "mod_prob"), "DFrame")
    expect_s4_class(assay(se4, "mod_prob"), "DFrame")

    expect_s4_class(assay(se1, "mod_prob")[[1]], "SparseMatrix")
    expect_s4_class(assay(se2, "mod_prob")[[1]], "SparseMatrix")
    expect_s4_class(assay(se3, "mod_prob")[[1]], "SparseMatrix")
    expect_s4_class(assay(se4, "mod_prob")[[1]], "SparseMatrix")

    expect_equal(vapply(assay(se1, "mod_prob"), ncol, 0), se1$n_reads)
    expect_equal(vapply(assay(se2, "mod_prob"), ncol, 0), se2$n_reads)
    expect_equal(vapply(assay(se3, "mod_prob"), ncol, 0), se3$n_reads)
    expect_equal(vapply(assay(se4, "mod_prob"), ncol, 0), se4$n_reads)

    # ... content se1
    expect_identical(unname(se1$n_reads), c(4L, 6L))
    expect_identical(dim(se1), c(12571L, 2L))
    modprob0 <- as.matrix(assay(se0, "mod_prob"))
    modprob1 <- as.matrix(assay(se1, "mod_prob"))
    shared_rows <- intersect(rownames(modprob0), rownames(modprob1))
    shared_cols <- intersect(colnames(modprob0), colnames(modprob1))
    expect_length(shared_rows, 8534L)
    expect_length(shared_cols, 10L)
    nonzero <- as.vector(modprob1[shared_rows, shared_cols]) > 0 &
        as.vector(modprob0[shared_rows, shared_cols]) > 0
    # plot(as.vector(modprob1[shared_rows, shared_cols])[nonzero],
    #      as.vector(modprob0[shared_rows, shared_cols])[nonzero])
    expect_equal(as.vector(modprob1[shared_rows, shared_cols])[nonzero],
                 as.vector(modprob0[shared_rows, shared_cols])[nonzero],
                 tolerance = 1e-6)
    expect_identical(se1$sample, names(modbamfiles))
    expect_equal(se1$qscore,
                 S4Vectors::SimpleList(
                     sample1 = c(14.1428003311157, 16.0126991271973,
                                 21.1338005065918, 20.3082008361816),
                     sample2 = c(12.9041996002197, 9.67461013793945,
                                 15.0149002075195, 15.1365995407104,
                                 17.7175006866455, 13.6647996902466)))

    # ... content se2
    expect_identical(unname(se2$n_reads), c(3L, 2L))
    expect_identical(dim(se2), dim(se3))
    expect_identical(unname(as.matrix(assay(se2, "mod_prob"))),
                     unname(as.matrix(assay(se3, "mod_prob"))))
    expect_identical(sub("^s", "sample", colnames(se2)), colnames(se3))
    expect_identical(se2$sample, sub("sample", "s", se3$sample))
    expect_identical(unname(se2$qscore), unname(se3$qscore))

    # ... content se3
    expect_identical(unname(se3$n_reads), c(3L, 2L))
    expect_identical(dim(se3), c(9266L, 2L))
    modprob3 <- as.matrix(assay(se3, "mod_prob"))
    shared_rows <- intersect(rownames(modprob0), rownames(modprob3))
    shared_cols <- intersect(colnames(modprob0), colnames(modprob3))
    expect_length(shared_rows, 7809L)
    expect_length(shared_cols, 5L)
    nonzero <- as.vector(modprob3[shared_rows, shared_cols]) > 0 &
        as.vector(modprob0[shared_rows, shared_cols]) > 0
    # plot(as.vector(modprob3[shared_rows, shared_cols])[nonzero],
    #      as.vector(modprob0[shared_rows, shared_cols])[nonzero])
    expect_equal(as.vector(modprob3[shared_rows, shared_cols])[nonzero],
                 as.vector(modprob0[shared_rows, shared_cols])[nonzero],
                 tolerance = 1e-6)
    expect_equal(se3$qscore,
                 S4Vectors::SimpleList(
                     sample1 = c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                     sample2 = c(9.67461013793945, 13.6647996902466)))

    # ... content se4
    expect_identical(unname(se4$n_reads), c(3L, 0L))
    expect_identical(dim(se4), c(6057L, 2L))
    expect_identical(dim(as.matrix(assay(se4, "mod_prob"))), c(6057L, 3L))
    expect_equal(se4$qscore,
                 S4Vectors::SimpleList(sample1 = c(14.1428003311157,
                                                   16.0126991271973,
                                                   20.3082008361816),
                                       sample2 = numeric(0)))
})
