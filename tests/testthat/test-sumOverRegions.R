suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(BiocParallel)
})

## -------------------------------------------------------------------------- ##
## Checks, sumOverRegions
## -------------------------------------------------------------------------- ##
test_that("sumOverRegions works", {
    # example data
    bmfile <- system.file("extdata", "modkit_pileup_1.bed.gz", package = "footprintR")
    se <- readBedMethyl(bmfile)
    regions1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 6940000, end = 7000000, names = "a"))
    regions2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(start = seq(6930000, 7000000, by = 10000), width = 10000))

    # invalid arguments
    expect_error(sumOverRegions(se = "error"))
    expect_error(sumOverRegions(se = se, regions = "error"))
    expect_error(sumOverRegions(se = se, regions = regions1, verbose = "error"))

    # expected results
    expect_warning(s1 <- sumOverRegions(se = se, regions = regions1))
    expect_message(expect_message(expect_message(
        s2 <- sumOverRegions(se = se, regions = regions2, verbose = TRUE, BPPARAM = BiocParallel::MulticoreParam(1L))
    )))
    expect_s4_class(s1, "RangedSummarizedExperiment")
    expect_s4_class(s2, "RangedSummarizedExperiment")
    expect_identical(dim(s1), c(1L, 1L))
    expect_identical(dim(s2), c(sum(overlapsAny(regions2, se)), 1L))
    expect_identical(assayNames(s1), c("Nmod", "Nvalid"))
    expect_identical(assayNames(s2), c("Nmod", "Nvalid"))
    expect_identical(rownames(s1), "a")
    expect_null(rownames(s2))
    expect_identical(unname(assay(s1, "Nvalid")[1,1]), 128240)
    expect_identical(unname(assay(s2, "Nvalid")[, 1]), c(14972, 67103, 61137))
})
