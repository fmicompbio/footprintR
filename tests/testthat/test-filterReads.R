test_that(".removeAllNAReads works", {
    naa1 <- naa2 <- NaArray(dim = c(10, 3))
    naa1[cbind(1:2, 1:2)] <- 1
    adat <- make_zero_col_DFrame(nrow = nrow(naa1))
    adat[["assay1"]] <- naa1
    adat[["assay2"]] <- naa2
    adat[["summary1"]] <- rowSums(naa1, na.rm = TRUE)
    adat[["summary2"]] <- rowSums(naa2, na.rm = FALSE)
    rownames(adat) <- paste0("r", seq.int(nrow(adat)))

    expect_error(.removeAllNAReads(x = adat, prune = "error"))

    res1 <- .removeAllNAReads(x = adat, prune = FALSE)
    res2 <- .removeAllNAReads(x = adat, prune = TRUE)

    expect_s4_class(res1, "DataFrame")
    expect_identical(dim(res1), c(10L, 4L))
    expect_identical(lapply(res1, ncol), list(assay1 = 2L, assay2 = 0L,
                                              summary1 = NULL, summary2 = NULL))
    expect_identical(rownames(res1), rownames(adat))

    expect_s4_class(res2, "DataFrame")
    expect_identical(dim(res2), c(10L, 2L))
    expect_identical(lapply(res2, ncol), list(assay1 = 2L, summary1 = NULL))
    expect_identical(rownames(res2), rownames(adat))
})

test_that("filterReads works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6990000",
                     modbase = "a", verbose = FALSE,
                     BPPARAM = BiocParallel::SerialParam())
    se <- flattenReadLevelAssay(se, keepReads = TRUE)
    se <- addReadStats(se, name = "qcc", stats = allReadStats,
                       BPPARAM = BiocParallel::SerialParam())

    expect_error(filterReads(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(filterReads(se = se, assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(filterReads(se = se, assayName = "error"),
                 "'assayName' must be one of")
    expect_error(filterReads(se = se, assayName = "Nvalid"),
                 "'assayName' must be one of")
    expect_error(filterReads(se = se, readInfoCol = 1),
                 "'readInfoCol' must be of class 'character'")
    expect_error(filterReads(se = se, readInfoCol = c("qcc", "readInfo")),
                 "'readInfoCol' must have length 1")
    expect_error(filterReads(se = se, readInfoCol = "missing"),
                 "'readInfoCol' must be one of")
    expect_error(filterReads(se = se, qcCol = 1),
                 "'qcCol' must be of class 'character'")
    expect_error(filterReads(se = se, qcCol = c("qcc", "readInfo")),
                 "'qcCol' must have length 1")
    expect_error(filterReads(se = se, qcCol = "missing"),
                 "'qcCol' must be one of")
    expect_error(filterReads(se = se, qcCol = "qcc", minQscore = "1"),
                 "'minQscore' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minQscore = c(1, 2)),
                 "'minQscore' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", maxEntropy = "1"),
                 "'maxEntropy' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", maxEntropy = c(1, 2)),
                 "'maxEntropy' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", maxFracLowConf = "1"),
                 "'maxFracLowConf' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", maxFracLowConf = c(0.5, 0.7)),
                 "'maxFracLowConf' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", minReadLength = "1"),
                 "'minReadLength' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minReadLength = c(1, 2)),
                 "'minReadLength' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedLength = "1"),
                 "'minAlignedLength' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedLength = c(1, 2)),
                 "'minAlignedLength' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedFraction = "1"),
                 "'minAlignedFraction' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedFraction = 2),
                 "'minAlignedFraction' must be within \\[0,1\\]")
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedFraction = c(0.5, 1)),
                 "'minAlignedFraction' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", prune = "1"),
                 "'prune' must be of class 'logical'")
    expect_error(filterReads(se = se, qcCol = "qcc", prune = c(TRUE, FALSE)),
                 "'prune' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", onlyStats = "1"),
                 "'onlyStats' must be of class 'logical'")
    expect_error(filterReads(se = se, qcCol = "qcc", onlyStats = c(TRUE, FALSE)),
                 "'onlyStats' must have length 1")

    ## Default arguments - no filtering
    out1 <- filterReads(se, qcCol = "qcc")
    setmp <- se
    metadata(setmp)$filteredOutReads = list(
        s1 = SVT_SparseArray(dim = c(0, 7),
                             dimnames = list(character(0),
                                             c("Qscore", "Entropy", "FracLowConf",
                                               "ReadLength", "AlignedLength",
                                               "AlignedFraction", "AllNA"))),
        s2 = SVT_SparseArray(dim = c(0, 7),
                             dimnames = list(character(0),
                                             c("Qscore", "Entropy", "FracLowConf",
                                               "ReadLength", "AlignedLength",
                                               "AlignedFraction", "AllNA"))))
    expect_identical(setmp, out1)

    ## Some filtering
    out1 <- filterReads(se, qcCol = "qcc", readInfoCol = "readInfo",
                        minQscore = 13, maxEntropy = 0.2,
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    expect_equal(dim(out1)[2], dim(se)[2])
    expect_equal(dim(out1), c(6534L, 2L))
    expect_equal(nrow(out1$qcc$s1), 2L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[c(2, 3)])
    expect_equal(nrow(out1$qcc$s2), 3L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[c(3, 5, 7)])
    expect_s4_class(metadata(out1)$filteredOutReads$s1, "SparseMatrix")
    expect_s4_class(metadata(out1)$filteredOutReads$s2, "SparseMatrix")
    expect_equal(dim(metadata(out1)$filteredOutReads$s1), c(8, 7))
    expect_equal(dim(metadata(out1)$filteredOutReads$s2), c(7, 7))
    expect_equal(colnames(metadata(out1)$filteredOutReads$s1),
                 c("Qscore", "Entropy", "FracLowConf", "ReadLength",
                   "AlignedLength", "AlignedFraction", "AllNA"))

    ## Only QC filtering
    out1 <- filterReads(se, qcCol = "qcc", readInfoCol = NULL,
                        minQscore = 13, maxEntropy = 0.2, maxFracLowConf = 0.1,
                        minReadLength = 8000, minAlignedLength = 5000,
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    ## calculate using length(which(rowSums(is_nonna(assay(se, "mod_prob")[["s1"]][, c(2, 3, 9)])) > 0 | rowSums(is_nonna(assay(se, "mod_prob")[["s2"]][, c(3, 5, 7, 10)])) > 0))
    expect_equal(dim(out1), c(6534L, 2L))
    ## calculate using se$qcc$s1$SEntrModProb < 0.2 & se$qcc$s1$FracLowConf < 0.1
    expect_equal(nrow(out1$qcc$s1), 2L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[c(2, 3)])
    expect_equal(nrow(out1$qcc$s2), 3L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[c(3, 5, 7)])

    ## Only read info filtering
    out1 <- filterReads(se, qcCol = NULL, readInfoCol = "readInfo",
                        minQscore = 13, maxEntropy = 0.2,
                        minReadLength = 8000, minAlignedLength = 5000,
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    expect_equal(dim(out1), c(7691L, 2L))
    expect_equal(nrow(out1$qcc$s1), 6L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[2:7])
    expect_equal(nrow(out1$qcc$s2), 5L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[3:7])

    ## Return filter stats only (compare to previous output)
    stats1 <- filterReads(se, qcCol = NULL, readInfoCol = "readInfo",
                          minQscore = 13, maxEntropy = 0.2,
                          minReadLength = 8000, minAlignedLength = 5000,
                          minAlignedFraction = 0.8, onlyStats = TRUE)
    expect_s4_class(stats1$s1, "SparseArray")
    expect_equal(dim(stats1$s1), c(4L, 7L))
    expect_equal(dim(stats1$s2), c(5L, 7L))
    expect_equal(stats1, metadata(out1)$filteredOutReads)
})

