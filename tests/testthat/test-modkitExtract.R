suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
})

## -------------------------------------------------------------------------- ##
## Checks, .modkitExtract
## -------------------------------------------------------------------------- ##
test_that(".modkitVersion works", {
    # invalid arguments
    expect_error(.modkitVersion(modkit_bin = 1L))

    # expected results
    expect_identical(.modkitVersion(modkit_bin = "error"), NA)
    rversion <- .modkitVersion(modkit_bin = file.path(R.home("bin"), "R"))
    expect_type(rversion, "character")
    expect_true(grepl("^R version ", rversion[1]))
})
