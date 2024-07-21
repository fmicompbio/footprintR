suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
})

## -------------------------------------------------------------------------- ##
## Checks, getReadDataForRegion
## -------------------------------------------------------------------------- ##
test_that("getReadDataForRegion works", {
    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    reg <- as("chr1:6940000-6955000", "GRanges")

    # invalid arguments
    expect_error(getReadDataForRegion(bamfile = "error"), "does not exist")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = "error"),
                 ".chr:start-end.")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = 1L),
                 "must be of class 'GRanges'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.modkitExtract = FALSE),
                 "must be of class 'list'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.modkitExtract = list("nonames")),
                 "must be an empty or a named list")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.modkitExtract = list(bamfile = modbamfile)),
                 "must not contain 'bamfile' or 'regions'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.readModkitExtract = "error"),
                 "must be of class 'list'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.readModkitExtract = list("nonames")),
                 "must be an empty or a named list")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.readModkitExtract = list(fnames = "not_allowed")),
                 "must not contain 'fnames'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.filterReadData = "error"),
                 "must be of class 'list'")
    expect_error(getReadDataForRegion(bamfile = modbamfile, region = reg,
                                      arglist.filterReadData = list("nonames")),
                 "must be NULL or a named list")

    # expected results
    # these are tested conditionally if `modkit` is available in the PATH
    # together with modkitExtract() tests, see: test-modkitExtract.R
})
