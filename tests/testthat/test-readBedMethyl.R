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
    ref <- system.file("extdata", "reference.fa.gz", package = "footprintR")

    # invalid arguments
    expect_error(readBedMethyl("error", 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'x', 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = c(nonexistent = 'm'), 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, nrows = -1, 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', seqinfo = "error", 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', seqinfo = c(100), 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', 
                               seqinfo = c(chr2 = 1000), 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', 
                               sequenceContextWidth = -1, 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', 
                               sequenceContextWidth = 1, 
                               sequenceReference = NULL, 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', 
                               sequenceContextWidth = 1, 
                               sequenceReference = "error", 
                               BPPARAM = BiocParallel::SerialParam()))
    expect_error(readBedMethyl(fname1, modbase = 'm', BPPARAM = "error"))
    expect_error(readBedMethyl(fname1, modbase = 'm', 
                               BPPARAM = BiocParallel::SerialParam(),
                               verbose = "error"))
    expect_error(readBedMethyl(c(a = fname1, a = fname2), 
                               modbase = c(a = 'm', a = 'a'), 
                               BPPARAM = BiocParallel::SerialParam()))

    # expected results
    se0 <- readBedMethyl(fnames = fname1, modbase = 'a', 
                         BPPARAM = BiocParallel::SerialParam())
    suppressMessages(
        expect_message(
            se1 <- readBedMethyl(fnames = fname1, modbase = 'm',
                                 BPPARAM = BiocParallel::SerialParam(),
                                 sequenceReference = ref, verbose = TRUE)
        )
    )
    suppressMessages(
        expect_message(
            se2 <- readBedMethyl(fnames = c(s2 = fname2), modbase = 'm', 
                                 BPPARAM = BiocParallel::SerialParam(),
                                 sequenceContextWidth = 1,
                                 sequenceReference = ref, verbose = TRUE)
        )
    )

    se12 <- readBedMethyl(fnames = c(fname1, fname2), modbase = 'm',
                          sequenceContextWidth = 1, sequenceReference = ref, 
                          BPPARAM = BiocParallel::SerialParam())
    suppressMessages(
        expect_message(
            se11 <- readBedMethyl(fnames = c(s1 = fname1, s1 = fname2),
                                  modbase = 'm', verbose = TRUE, 
                                  BPPARAM = BiocParallel::SerialParam())
        )
    )
    expect_s4_class(se0, "SummarizedExperiment")
    expect_s4_class(se1, "SummarizedExperiment")
    expect_s4_class(se2, "SummarizedExperiment")
    expect_s4_class(se12, "SummarizedExperiment")
    expect_s4_class(se11, "SummarizedExperiment")
    expect_identical(colnames(colData(se0)), c("sample", "modbase"))
    expect_identical(colnames(colData(se1)), c("sample", "modbase"))
    expect_identical(colnames(colData(se2)), c("sample", "modbase"))
    expect_identical(colnames(colData(se12)), c("sample", "modbase"))
    expect_identical(colnames(colData(se11)), c("sample", "modbase"))
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
    expect_identical(rowSums(assay(se12, "Nmod")), rowSums(assay(se11, "Nmod")))
    expect_identical(rowSums(assay(se12, "Nvalid")), rowSums(assay(se11, "Nvalid")))
    expect_true("sequenceContext" %in% colnames(rowData(se2)))
    expect_equal(as.integer(table(as.character(rowData(se2)$sequenceContext))),
                 c(844L, 7535L, 801L, 820L))
})
