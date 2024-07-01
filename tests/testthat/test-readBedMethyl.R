suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(GenomicRanges)
})

## -------------------------------------------------------------------------- ##
## Checks, readBedMethyl
## -------------------------------------------------------------------------- ##
test_that("readBedMethyl works", {
    # example data
    fname1 <- system.file("extdata", "modkit_pileup_1.bed.gz", package = "footprintR")
    fname2 <- system.file("extdata", "modkit_pileup_2.bed.gz", package = "footprintR")

    # invalid arguments
    expect_error(readBedMethyl("error"))
    expect_error(readBedMethyl(fname1, nrows = -1))
    expect_error(readBedMethyl(fname1, seqlens = "error"))
    expect_error(readBedMethyl(fname1, seqlens = c(100)))
    expect_error(readBedMethyl(fname1, verbose = "error"))

    # expected results
    se0 <- readBedMethyl(fnames = fname1, modbase = 'x')
    suppressMessages(
        expect_message(
            se1 <- readBedMethyl(fnames = fname1, modbase = 'm', verbose = TRUE)
        )
    )
    se2 <- readBedMethyl(fnames = c(s2 = fname2))
    se12 <- readBedMethyl(fnames = c(fname1, fname2))
    suppressMessages(
        expect_message(
            se11 <- readBedMethyl(fnames = c(s1 = fname1, s1 = fname2), verbose = TRUE)
        )
    )
    expect_s4_class(se0, "SummarizedExperiment")
    expect_s4_class(se1, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    expect_s4_class(se12, "SummarizedExperiment")
    expect_s4_class(se11, "SummarizedExperiment")
    expect_identical(dim(se0), c(0L, 1L))
    expect_identical(dim(se1), c(10000L, 1L))
    expect_identical(dim(se2), c(10000L, 1L))
    expect_identical(dim(se12), c(12020L, 2L))
    expect_identical(dim(se11), c(12020L, 1L))
    expect_identical(colnames(se1), "s1")
    expect_identical(colnames(se2), "s2")
    expect_identical(colnames(se12), c("s1","s2"))
    expect_identical(colnames(se11), "s1")
    expect_true(all(overlapsAny(rowRanges(se1), rowRanges(se12))))
    expect_true(all(overlapsAny(rowRanges(se2), rowRanges(se12))))
    expect_true(all(overlapsAny(rowRanges(se1), rowRanges(se11))))
    expect_true(all(overlapsAny(rowRanges(se2), rowRanges(se11))))
    expect_identical(assayNames(se1), c("Nmod", "Nvalid"))
    expect_identical(assayNames(se2), c("Nmod", "Nvalid"))
    expect_identical(assayNames(se12), c("Nmod", "Nvalid"))
    expect_identical(assayNames(se11), c("Nmod", "Nvalid"))
    i1to12 <- GenomicRanges::match(rowRanges(se1), rowRanges(se12))
    i2to12 <- GenomicRanges::match(rowRanges(se2), rowRanges(se12))
    expect_identical(assay(se1, "Nmod"),
                     assay(se12, "Nmod")[i1to12, 1, drop = FALSE])
    expect_identical(assay(se1, "Nvalid"),
                     assay(se12, "Nvalid")[i1to12, 1, drop = FALSE])
    expect_identical(assay(se2, "Nmod"),
                     assay(se12, "Nmod")[i2to12, 2, drop = FALSE])
    expect_identical(assay(se2, "Nvalid"),
                     assay(se12, "Nvalid")[i2to12, 2, drop = FALSE])
    expect_identical(sum(assay(se12, "Nmod")), sum(assay(se11, "Nmod")))
    expect_identical(sum(assay(se12, "Nvalid")), sum(assay(se11, "Nvalid")))
})
