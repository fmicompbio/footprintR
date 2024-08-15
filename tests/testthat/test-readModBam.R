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
    expect_error(readModBam(modbamfiles, "error", "a"),
                 "GRanges object must contain")
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
                      modbase = "m",
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

    expect_identical(se1$sample, rep(names(modbamfiles), c(4, 6)))
    expect_identical(se2$sample, rep(c("s1", "s2"), c(3, 2)))
    expect_identical(se3$sample, rep(names(modbamfiles), c(3, 2)))
    expect_identical(se4$sample, character(0))

    expect_identical(assayNames(se1), "mod_prob")
    expect_identical(assayNames(se2), "mod_prob")
    expect_identical(assayNames(se3), "mod_prob")
    expect_identical(assayNames(se4), "mod_prob")

    # ... content se1
    expect_identical(dim(se1), c(12571L, 10L))
    ov_row <- findOverlaps(se0, se1)
    expect_length(ov_row, 8534L)
    ov_col <- match(colnames(se0), colnames(se1))
    expect_true(sum(!is.na(ov_col)) == 10L)
    # plot(as.vector(assay(se1, "mod_prob")[subjectHits(ov_row), na.omit(ov_col)]),
    #      as.vector(assay(se0, "mod_prob")[queryHits(ov_row), !is.na(ov_col)]))
    expect_true(cor(as.vector(assay(se1, "mod_prob")[subjectHits(ov_row), na.omit(ov_col)]),
                    as.vector(assay(se0, "mod_prob")[queryHits(ov_row), !is.na(ov_col)])) > 0.99)

    # ... content se2
    expect_identical(dim(se2), dim(se3))
    expect_identical(unname(assay(se2, "mod_prob")),
                     unname(assay(se3, "mod_prob")))

    # ... content se3
    expect_identical(dim(se3), c(9266L, 5L))
    ov_row <- findOverlaps(se0, se3)
    expect_length(ov_row, 7809L)
    ov_col <- match(colnames(se0), colnames(se3))
    expect_true(sum(!is.na(ov_col)) == 5L)
    # plot(as.vector(assay(se3, "mod_prob")[subjectHits(ov_row), na.omit(ov_col)]),
    #      as.vector(assay(se0, "mod_prob")[queryHits(ov_row), !is.na(ov_col)]))
    expect_true(cor(as.vector(assay(se3, "mod_prob")[subjectHits(ov_row), na.omit(ov_col)]),
                    as.vector(assay(se0, "mod_prob")[queryHits(ov_row), !is.na(ov_col)])) > 0.99)

    # ... content se4
    expect_identical(dim(se4), c(0L, 0L))
})
