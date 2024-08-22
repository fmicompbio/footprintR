suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
})

## -------------------------------------------------------------------------- ##
## Checks, addReadsSummary
## -------------------------------------------------------------------------- ##
test_that("addReadsSummary works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz", package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a")
    se0 <- se
    colData(se0) <- NULL
    se1 <- se
    rowRanges(se1) <- NULL

    # invalid arguments
    expect_error(addReadsSummary(se = "error"))
    expect_error(addReadsSummary(se = se0))
    expect_error(addReadsSummary(se = se, assay.type = "error"))
    expect_error(addReadsSummary(se = se, statistics = "error"))
    expect_error(addReadsSummary(se = se, keep.reads = "error"))
    expect_error(addReadsSummary(se = se, verbose = "error"))

    # expected results
    expect_message(expect_message(
        s1 <- addReadsSummary(se = se, statistics = c("Nmod", "Nvalid", "FracMod", "Pmod", "AvgConf"), verbose = TRUE)
    ))
    s2 <- addReadsSummary(se = se, statistics = "FracMod", keep.reads = TRUE)
    s3 <- addReadsSummary(se = se1, statistics = "FracMod", keep.reads = TRUE)
    expect_s4_class(s1, "RangedSummarizedExperiment")
    expect_s4_class(s2, "RangedSummarizedExperiment")
    expect_s4_class(s3, "SummarizedExperiment")
    expect_identical(dim(s1), c(nrow(se), length(unique(se$sample))))
    expect_identical(dim(s2), c(nrow(se), length(unique(se$sample))))
    expect_identical(dim(s3), c(nrow(se), length(unique(se$sample))))
    expect_identical(assayNames(s1), c("Nmod", "Nvalid", "FracMod", "Pmod", "AvgConf"))
    expect_identical(assayNames(s2), c("FracMod", "mod_prob"))
    expect_identical(assayNames(s3), c("FracMod", "mod_prob"))
    expect_identical(rownames(s1), rownames(se))
    expect_identical(rownames(s2), rownames(se))
    expect_null(rownames(s3))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) >= 0.5),
                 sum(assay(s1, "Nmod")))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) != 0.0),
                 sum(assay(s1, "Nvalid")))
    expect_s4_class(assay(s2, "mod_prob"), "DataFrame")
    expect_s4_class(assay(s2, "mod_prob")[,1], "SparseMatrix")
    expect_identical(dim(assay(s2, "mod_prob")[,1]), dim(assay(se, "mod_prob")[,1]))
})
