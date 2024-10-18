suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(S4Vectors)
})

## -------------------------------------------------------------------------- ##
## Checks, subsetReads
## -------------------------------------------------------------------------- ##
test_that("subsetReads works", {
    # example data
    modbamfiles <- system.file("extdata",
                            c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                            package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles,
                     regions = "chr1:6940000-6955000", modbase = "a")
    se2 <- addReadsSummary(se, keep.reads = FALSE)
    se3 <- se
    assays(se3) <- SimpleList(mod_prob = assay(se, "mod_prob"),
                              mod_prob2 = assay(se, "mod_prob"))
    metadata(se3)$readLevelData$assayNames <- c("mod_prob", "mod_prob2")

    read_ids <- c("s1-233e48a7-f379-4dcf-9270-958231125563",
                  "s2-034b625e-6230-4f8d-a713-3a32cd96c298")

    # expected errors
    expect_error(subsetReads(se = "error", reads = read_ids),
                 "must be of class 'SummarizedExperiment'")
    expect_error(subsetReads(se = se, reads = c(read_ids, "error")),
                 "'reads' contains unknown identifiers")
    expect_error(subsetReads(se = se, reads = read_ids, prune = "error"),
                 "'prune' must be of class 'logical'")
    expect_error(subsetReads(se = se, reads = list(error = 1)),
                 "'reads' of type list must have names in")
    expect_error(subsetReads(se = se, reads = list(s1 = "error")),
                 "'reads' for sample 's1' contains unknown read names")
    expect_error(subsetReads(se, list(s2 = c(NA, "s2-034b625e-6230-4f8d-a713-3a32cd96c298"))),
                 "'reads' for sample 's2' contains unknown read names")
    expect_error(subsetReads(se = se, reads = list(s1 = 99)),
                 "'reads' for sample 's1' contains out-of-range indices")
    expect_error(subsetReads(se = se, reads = list(s1 = FALSE)),
                 "logical 'reads' for sample 's1' must be of length 3")
    expect_error(subsetReads(se = se, reads = TRUE),
                 "must be either a character vector")

    # expected results
    # ... no read-level assay
    expect_warning(seSub <- subsetReads(se = se2, reads = read_ids),
                   "'se' contains no read-level assays")
    expect_identical(seSub, se2)

    # ... several read-level assays
    seSub <- subsetReads(se = se3, reads = read_ids)
    expect_identical(assay(seSub, "mod_prob"), assay(seSub, "mod_prob2"))

    # ... prune
    seSub <- subsetReads(se = se, reads = list(s1 = 1), prune = FALSE)
    expect_identical(dim(seSub), c(7967L, 2L))
    expect_identical(nrow(colData(seSub)), 2L)
    expect_identical(ncols(assays(seSub)), c(mod_prob = 2L))
    expect_identical(colnames(seSub), c("s1", "s2"))
    seSub <- subsetReads(se = se, reads = list(s1 = 1), prune = TRUE)
    expect_identical(dim(seSub), c(7967L, 1L))
    expect_identical(nrow(colData(seSub)), 1L)
    expect_identical(ncols(assays(seSub)), c(mod_prob = 1L))
    expect_identical(colnames(seSub), c("s1"))

    # ... subset by character vector
    seSub <- subsetReads(se, c("s1-233e48a7-f379-4dcf-9270-958231125563",
                               "s2-034b625e-6230-4f8d-a713-3a32cd96c298"))
    expect_s4_class(seSub, "RangedSummarizedExperiment")

    # ... subset by list with numeric index
    seSub2 <- subsetReads(se, list(s1 = 1, s2 = 1))
    expect_identical(seSub2, seSub)

    # ... subset by list with logical index
    seSub2 <- subsetReads(se, list(s1 = c(TRUE, FALSE, FALSE), s2 = c(TRUE, FALSE)))
    expect_identical(seSub2, seSub)

    # ... subset by list with character index
    seSub <- subsetReads(se, list(s1 = "s1-233e48a7-f379-4dcf-9270-958231125563",
                                  s2 = "s2-034b625e-6230-4f8d-a713-3a32cd96c298"))
    expect_identical(seSub2, seSub)

    # ... invert
    seSub <- subsetReads(se, list(s1 = 2, s2 = 2))
    seSub2 <- subsetReads(se, list(s1 = c(1, 3), s2 = 1), invert = TRUE)
    expect_identical(seSub, seSub2)
})
