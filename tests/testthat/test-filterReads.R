test_that("filterReads works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6990000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se, keep.reads = TRUE)
    se <- addReadStats(se, name = "qcc")
    
    expect_error(filterReads(se = "error"), 
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(filterReads(se = se, assay.type.read = 1),
                 "'assay.type.read' must be of class 'character'")
    expect_error(filterReads(se = se, assay.type.read = "error"),
                 "'assay.type.read' must be one of")
    expect_error(filterReads(se = se, assay.type.read = "Nvalid"),
                 "'assay.type.read' must be one of")
    expect_error(filterReads(se = se, readInfoCol = 1),
                 "'readInfoCol' must be of class 'character'")
    expect_error(filterReads(se = se, readInfoCol = c("qcc", "read_info")),
                 "'readInfoCol' must have length 1")
    expect_error(filterReads(se = se, readInfoCol = "missing"),
                 "'readInfoCol' must be one of")
    expect_error(filterReads(se = se, qcCol = 1),
                 "'qcCol' must be of class 'character'")
    expect_error(filterReads(se = se, qcCol = c("qcc", "read_info")),
                 "'qcCol' must have length 1")
    expect_error(filterReads(se = se, qcCol = "missing"),
                 "'qcCol' must be one of")
    expect_error(filterReads(se = se, qcCol = "qcc", minQscore = "1"), 
                 "'minQscore' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minQscore = c(1, 2)), 
                 "'minQscore' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", minEntropy = "1"), 
                 "'minEntropy' must be of class 'numeric'")
    expect_error(filterReads(se = se, qcCol = "qcc", minEntropy = c(1, 2)), 
                 "'minEntropy' must have length 1")
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
    expect_error(filterReads(se = se, qcCol = "qcc", minAlignedFraction = c(1, 2)), 
                 "'minAlignedFraction' must have length 1")
    expect_error(filterReads(se = se, qcCol = "qcc", prune = "1"), 
                 "'prune' must be of class 'logical'")
    expect_error(filterReads(se = se, qcCol = "qcc", prune = c(TRUE, FALSE)), 
                 "'prune' must have length 1")
    
    ## Default arguments - no filtering
    out1 <- filterReads(se, qcCol = "qcc")
    setmp <- se
    metadata(setmp)$filteredOutReads = list(
        s1 = SVT_SparseArray(dim = c(0, 6),
                             dimnames = list(character(0), 
                                             c("Qscore", "Entropy", "ReadLength", 
                                               "AlignedLength",
                                               "AlignedFraction", "AllNA"))),
        s2 = SVT_SparseArray(dim = c(0, 6),
                             dimnames = list(character(0), 
                                             c("Qscore", "Entropy", "ReadLength", 
                                               "AlignedLength",
                                               "AlignedFraction", "AllNA"))))
    expect_identical(setmp, out1)
    
    ## Some filtering
    out1 <- filterReads(se, qcCol = "qcc", readInfoCol = "read_info", 
                        minQscore = 13, minEntropy = 0.2, 
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    expect_equal(dim(out1)[2], dim(se)[2])
    expect_equal(dim(out1), c(7372L, 2L))
    expect_equal(nrow(out1$qcc$s1), 5L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[c(4, 5, 6, 7, 10)])
    expect_equal(nrow(out1$qcc$s2), 2L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[c(4, 6)])
    expect_s4_class(metadata(out1)$filteredOutReads$s1, "SparseMatrix")
    expect_s4_class(metadata(out1)$filteredOutReads$s2, "SparseMatrix")
    expect_equal(dim(metadata(out1)$filteredOutReads$s1), c(5, 6))
    expect_equal(dim(metadata(out1)$filteredOutReads$s2), c(8, 6))
    expect_equal(colnames(metadata(out1)$filteredOutReads$s1), 
                 c("Qscore", "Entropy", "ReadLength", "AlignedLength",
                   "AlignedFraction", "AllNA"))
    
    ## Only QC filtering
    out1 <- filterReads(se, qcCol = "qcc", readInfoCol = NULL, 
                        minQscore = 13, minEntropy = 0.2, 
                        minReadLength = 8000, minAlignedLength = 5000, 
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    expect_equal(dim(out1), c(8930L, 2L))
    expect_equal(nrow(out1$qcc$s1), 7L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[-c(2, 3, 9)])
    expect_equal(nrow(out1$qcc$s2), 6L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[c(1, 2, 4, 6, 8, 9)])
    
    ## Only read info filtering
    out1 <- filterReads(se, qcCol = NULL, readInfoCol = "read_info", 
                        minQscore = 13, minEntropy = 0.2, 
                        minReadLength = 8000, minAlignedLength = 5000, 
                        minAlignedFraction = 0.8)
    expect_s4_class(out1, "SummarizedExperiment")
    expect_equal(dim(out1), c(7691L, 2L))
    expect_equal(nrow(out1$qcc$s1), 6L)
    expect_equal(rownames(out1$qcc$s1), rownames(se$qcc$s1)[2:7])
    expect_equal(nrow(out1$qcc$s2), 5L)
    expect_equal(rownames(out1$qcc$s2), rownames(se$qcc$s2)[3:7])
})

    