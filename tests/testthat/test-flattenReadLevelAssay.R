suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
})

## -------------------------------------------------------------------------- ##
## Checks, flattenReadLevelAssay
## -------------------------------------------------------------------------- ##
test_that("flattenReadLevelAssay works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz", package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a", 
                            BPPARAM = BiocParallel::SerialParam())
    se0 <- se
    colData(se0) <- NULL

    # invalid arguments
    expect_error(flattenReadLevelAssay(se = "error"))
    expect_error(flattenReadLevelAssay(se = se, assayName = "error"))
    expect_error(flattenReadLevelAssay(se = se, statistics = "error"))
    expect_error(flattenReadLevelAssay(se = se, keepReads = "error"))
    expect_error(flattenReadLevelAssay(se = se, replaceExisting = "error"))
    expect_error(flattenReadLevelAssay(se = se, replaceExisting = c(TRUE, FALSE)))
    expect_error(flattenReadLevelAssay(se = se, verbose = "error"))

    # expected results
    suppressMessages({
        s1 <- flattenReadLevelAssay(se = se,
                                    statistics = c("Nmod", "Nvalid", "FracMod",
                                                   "Pmod", "AvgConf", "Mean",
                                                   "Sum"),
                                    modProbThreshold = 0.5,
                                    keepReads = FALSE, verbose = TRUE)
        })
    s2 <- flattenReadLevelAssay(se = se, statistics = "FracMod")
    s3 <- flattenReadLevelAssay(se = s2, statistics = "FracMod",
                                replaceExisting = TRUE)
    expect_identical(s2, s3)
    expect_warning(
        s3 <- flattenReadLevelAssay(se = s2, statistics = "FracMod",
                                    replaceExisting = FALSE))
    expect_identical(s2, s3)
    expect_s4_class(s1, "RangedSummarizedExperiment")
    expect_s4_class(s2, "RangedSummarizedExperiment")
    expect_identical(dim(s1), c(nrow(se), length(colnames(se))))
    expect_identical(dim(s2), c(nrow(se), length(colnames(se))))
    expect_identical(assayNames(s1), c("Nmod", "Nvalid", "FracMod", "Pmod", 
                                       "AvgConf", "Mean", "Sum"))
    expect_identical(assayNames(s2), c("mod_prob", "FracMod"))
    expect_identical(rownames(s1), rownames(se))
    expect_identical(rownames(s2), rownames(se))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) >= 0.5, na.rm = TRUE),
                 sum(assay(s1, "Nmod"), na.rm = TRUE))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) >= 0.0, na.rm = TRUE),
                 sum(assay(s1, "Nvalid"), na.rm = TRUE))
    expect_identical(assay(s1, "Pmod"), assay(s1, "Mean"))
    expect_s4_class(assay(s2, "mod_prob"), "DataFrame")
    expect_s4_class(assay(s2, "mod_prob")[,1], "NaMatrix")
    expect_identical(dim(assay(s2, "mod_prob")[,1]), dim(assay(se, "mod_prob")[,1]))
    expect_identical(assay(s1, "Sum") / assay(s1, "Nvalid"), assay(s1, "Mean"))
    
    # change modProbTreshold
    suppressMessages({
        s1 <- flattenReadLevelAssay(se = se,
                                    statistics = c("Nmod", "Nvalid", "FracMod",
                                                   "Pmod", "AvgConf", "Mean",
                                                   "Sum"),
                                    modProbThreshold = 0.75,
                                    keepReads = FALSE, verbose = TRUE)
    })
    expect_s4_class(s1, "RangedSummarizedExperiment")
    expect_identical(dim(s1), c(nrow(se), length(colnames(se))))
    expect_identical(assayNames(s1), c("Nmod", "Nvalid", "FracMod", "Pmod", 
                                       "AvgConf", "Mean", "Sum"))
    expect_identical(rownames(s1), rownames(se))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) >= 0.75, na.rm = TRUE),
                 sum(assay(s1, "Nmod"), na.rm = TRUE))
    expect_equal(sum(as.matrix(assay(se, "mod_prob")) >= 0.0, na.rm = TRUE),
                 sum(assay(s1, "Nvalid"), na.rm = TRUE))
    expect_identical(assay(s1, "Pmod"), assay(s1, "Mean"))
    expect_identical(assay(s1, "Sum") / assay(s1, "Nvalid"), assay(s1, "Mean"))
})
