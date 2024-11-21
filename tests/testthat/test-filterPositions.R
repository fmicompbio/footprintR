test_that(".filterPositionsByCoverage works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- flattenReadLevelAssay(se, keepReads = TRUE)

    ## Calculate coverage
    cov_total <- rowSums(as.matrix(as.matrix(assay(se, "mod_prob")) >= 0), na.rm = TRUE)
    cov_bysample <- as.matrix(endoapply(assay(se, "mod_prob"), function(x) {
        rowSums(x >= 0, na.rm = TRUE)
    }))
    expect_equal(rowSums(cov_bysample), cov_total, ignore_attr = TRUE)

    expect_error(.filterPositionsByCoverage(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.filterPositionsByCoverage(se = se, assayName = "missing"),
                 "must be one of")
    expect_error(.filterPositionsByCoverage(se = se, assayName = c("Nvalid", "Nmod")),
                 "'assayName' must have length 1")
    expect_error(.filterPositionsByCoverage(se = se, assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(.filterPositionsByCoverage(se = se, assayName = "Nvalid",
                                            minCov = "1"),
                 "'minCov' must be of class 'numeric'")
    expect_error(.filterPositionsByCoverage(se = se, assayName = "Nvalid",
                                            minCov = c(1, 2)),
                 "'minCov' must have length 1")
    expect_error(.filterPositionsByCoverage(se = se, assayName = "Nvalid",
                                            minCov = 1, minNbrSamples = "1"),
                 "'minNbrSamples' must be of class 'numeric'")
    expect_error(.filterPositionsByCoverage(se = se, assayName = "Nvalid",
                                            minCov = 1, minNbrSamples = c(1, 2)),
                 "'minNbrSamples' must have length 1")

    se1 <- .filterPositionsByCoverage(se, assayName = "Nvalid", minCov = 10,
                                      minNbrSamples = NULL)
    se2 <- .filterPositionsByCoverage(se, assayName = "mod_prob", minCov = 10,
                                      minNbrSamples = NULL)
    expect_identical(se1, se2)
    w <- which(cov_total >= 10)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), names(w))
    expect_equal(nrow(se1), 1783L)

    se1 <- .filterPositionsByCoverage(se, assayName = "Nvalid", minCov = 5,
                                      minNbrSamples = 1)
    se2 <- .filterPositionsByCoverage(se, assayName = "mod_prob", minCov = 5,
                                      minNbrSamples = 1)
    expect_identical(se1, se2)
    w <- which(cov_bysample[, 1] >= 5 | cov_bysample[, 2] >= 5)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(se)[w])
    expect_equal(nrow(se1), 4212L)

    se1 <- .filterPositionsByCoverage(se, assayName = "Nvalid", minCov = 5,
                                      minNbrSamples = 2)
    se2 <- .filterPositionsByCoverage(se, assayName = "mod_prob", minCov = 5,
                                      minNbrSamples = 2)
    expect_identical(se1, se2)
    w <- which(cov_bysample[, 1] >= 5 & cov_bysample[, 2] >= 5)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(se)[w])
    expect_equal(nrow(se1), 1606L)

    se1 <- .filterPositionsByCoverage(se, assayName = "Nvalid", minCov = 5,
                                      minNbrSamples = 3)
    se2 <- .filterPositionsByCoverage(se, assayName = "mod_prob", minCov = 5,
                                      minNbrSamples = 3)
    expect_identical(se1, se2)
    expect_equal(nrow(se1), 0L)
})

test_that(".keepPositionsBySequenceContext works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)

    expect_error(.keepPositionsBySequenceContext(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.keepPositionsBySequenceContext(se = se,
                                                 sequenceContext = "ACT"),
                 "No sequence context found in `rowData(se)$sequenceContext`",
                 fixed = TRUE)
    expect_identical(se, .keepPositionsBySequenceContext(se = se))

    gnm <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
                                                    package = "footprintR"))
    se <- addSeqContext(se, sequenceContextWidth = 3,
                        sequenceReference = gnm)

    setmp <- se
    rowData(setmp)$sequenceContext <- as.character(rowData(setmp)$sequenceContext)
    expect_error(.keepPositionsBySequenceContext(se = setmp,
                                                 sequenceContext = "ACT"),
                 "must be of class 'DNAStringSet'")

    se1 <- .keepPositionsBySequenceContext(se = se, sequenceContext = "TAG")
    w <- which(as.character(rowData(se)$sequenceContext) == "TAG")
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(se)[w])
    ## Applying the same filter again doesn't do anything
    se2 <- .keepPositionsBySequenceContext(se = se1, sequenceContext = "TAG")
    expect_identical(se1, se2)

    ## Try with IUPAC
    se1 <- .keepPositionsBySequenceContext(se = se, sequenceContext = "WAG")
    w <- which(as.character(rowData(se)$sequenceContext) %in% c("TAG", "AAG"))
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(se)[w])

    se1 <- .keepPositionsBySequenceContext(se = se, sequenceContext = "NNN")
    expect_identical(rowData(se), rowData(se1))

    ## Manually insert Ns - make sure that these are not retained
    secopy <- se
    rowData(secopy)$sequenceContext[1:5] <- rep("NNN", 5)
    se1 <- .keepPositionsBySequenceContext(se = secopy, sequenceContext = "WAG")
    w <- which(as.character(rowData(secopy)$sequenceContext) %in% c("TAG", "AAG"))
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(secopy)[w])

    ## Padding
    sec5 <- addSeqContext(se, sequenceContextWidth = 5,
                          sequenceReference = gnm)
    sec3 <- addSeqContext(se, sequenceContextWidth = 3,
                          sequenceReference = gnm)
    se1 <- .keepPositionsBySequenceContext(se = sec5, sequenceContext = "NTAGN")
    se2 <- .keepPositionsBySequenceContext(se = sec3, sequenceContext = "TAG")
    se3 <- .keepPositionsBySequenceContext(se = sec5, sequenceContext = "TAG")
    expect_equal(nrow(se1), nrow(se2))
    expect_equal(rownames(se1), rownames(se2))
    expect_false(nrow(se1) == nrow(se3))
})

test_that(".removeAllNAPositions works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- flattenReadLevelAssay(se)
    # Subset reads to make sure there are positions with all NAs
    se <- subsetReads(se, reads = list(s1 = c(1, 2), s2 = c(1, 3)))

    expect_error(.removeAllNAPositions(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.removeAllNAPositions(se = se,
                                       assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(.removeAllNAPositions(se = se,
                                       assayName = c("mod_prob", "mod_prob")),
                 "'assayName' must have length 1")
    expect_error(.removeAllNAPositions(se = se,
                                       assayName = "missing"),
                 "'assayName' must be one of")
    expect_error(.removeAllNAPositions(se = se,
                                       assayName = "Nvalid"),
                 "'assayName' must be one of")

    se1 <- .removeAllNAPositions(se, assayName = "mod_prob")
    namat <- as.matrix(assay(se1, "mod_prob"))
    expect_s4_class(namat, "NaArray")
    expect_equal(ncol(namat), 4L)
    w <- which(rowSums(!is.na(as.matrix(as.matrix(assay(se, "mod_prob"))))) > 0)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), names(w))
    ## Applying the same filter again doesn't do anything
    se2 <- .removeAllNAPositions(se1, assayName = "mod_prob")
    expect_identical(se1, se2)
})

test_that(".pruneAmbiguousStrandPositions works", {
    # generate some example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- flattenReadLevelAssay(se)

    expect_error(.pruneAmbiguousStrandPositions(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assayName = "missing"),
                 "'assayName' must be one of")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assayName = c("Nvalid", "Nmod")),
                 "'assayName' must have length 1")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assayName = "Nvalid",
                                                verbose = "TRUE"),
                 "'verbose' must be of class 'logical'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assayName = "Nvalid",
                                                verbose = c(TRUE, FALSE)),
                 "'verbose' must have length 1")

    expect_equal(nrow(se), 9127L)
    expect_equal(length(unique(paste0(seqnames(rowRanges(se)),
                                      pos(rowRanges(se))))), 8955L)
    sefilt <- .pruneAmbiguousStrandPositions(se, assayName = "Nvalid",
                                             verbose = FALSE)
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
    expect_message(expect_message(expect_message({
        sefilt <- .pruneAmbiguousStrandPositions(se, assayName = "mod_prob",
                                                 verbose = TRUE)},
        "172 rows removed to ensure"), "172 rows removed to ensure")
    )
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
    expect_identical(se[rownames(sefilt), ], sefilt)
    expect_message(expect_message(expect_message({
        sefilt <- .pruneAmbiguousStrandPositions(sefilt, assayName = "mod_prob",
                                                 verbose = TRUE)},
        "No genomic positions represented by multiple rows found"),
        "No genomic positions represented by multiple rows found")
    )
})

test_that("filterPositions works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- flattenReadLevelAssay(se)
    se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = reffile)

    expect_error(filterPositions(se = "missing"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(filterPositions(se = se, filters = 1),
                 "'filters' must be of class 'character'")
    expect_error(filterPositions(se = se, filters = "error"),
                 "must be one of")

    sefilt <- filterPositions(se, c("sequenceContext", "coverage",
                                    "repeated.positions", "all.na"),
                              minCov = 5, sequenceContext = "TAG")
    expect_gte(min(rowSums(assay(sefilt, "Nvalid"))), 5L)
    expect_equal(nrow(sefilt), 251L)
    expect_true(all(as.character(rowData(sefilt)$sequenceContext) %in% c("TAG")))
    expect_false(any(duplicated(paste0(seqnames(rowRanges(sefilt)), ":",
                                       pos(rowRanges(sefilt))))))

    ## Filter out reads that are NA in all positions
    sefilt <- filterPositions(se, c("sequenceContext"),
                              sequenceContext = "ATG")
    expect_equal(nrow(sefilt), 8L)
    expect_true(all(as.character(rowData(sefilt)$sequenceContext) %in% c("ATG")))
    expect_equal(ncol(assay(sefilt, "mod_prob")$s1), 7L)
    expect_equal(colnames(assay(sefilt, "mod_prob")$s1),
                 colnames(assay(se, "mod_prob")$s1)[c(1, 3, 5, 6, 7, 8, 9)])
    expect_equal(colSums(is_nonna(assay(sefilt, "mod_prob")$s1)), c(1, 1, 1, 1, 1, 1, 1),
                 ignore_attr = TRUE)
    expect_equal(ncol(assay(sefilt, "mod_prob")$s2), 5L)
    expect_equal(colnames(assay(sefilt, "mod_prob")$s2),
                 colnames(assay(se, "mod_prob")$s2)[c(2, 4, 6, 7, 8)])
    expect_equal(colSums(is_nonna(assay(sefilt, "mod_prob")$s2)), c(1, 1, 1, 2, 3),
                 ignore_attr = TRUE)
})

