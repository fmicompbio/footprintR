test_that(".filterPositionsByCoverage works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se, keep.reads = TRUE)
    cov <- rowSums(as.matrix(as.matrix(assay(se, "mod_prob")) > 0), na.rm = TRUE)
    
    expect_error(.filterPositionsByCoverage(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.filterPositionsByCoverage(se = se, assay.type = "missing"),
                 "must be one of")
    expect_error(.filterPositionsByCoverage(se = se, assay.type = c("Nvalid", "Nmod")),
                 "'assay.type' must have length 1")
    expect_error(.filterPositionsByCoverage(se = se, assay.type = 1),
                 "'assay.type' must be of class 'character'")
    expect_error(.filterPositionsByCoverage(se = se, assay.type = "Nvalid", 
                                            min.cov = "1"),
                 "'min.cov' must be of class 'numeric'")
    expect_error(.filterPositionsByCoverage(se = se, assay.type = "Nvalid", 
                                            min.cov = c(1, 2)),
                 "'min.cov' must have length 1")
    
    se1 <- .filterPositionsByCoverage(se, assay.type = "Nvalid", min.cov = 10)
    se2 <- .filterPositionsByCoverage(se, assay.type = "mod_prob", min.cov = 10)
    expect_identical(se1, se2)
    w <- which(cov >= 10)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), names(w))
    expect_equal(nrow(se1), 1783L)
})

test_that(".keepPositionsBySequenceContext works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    
    expect_error(.keepPositionsBySequenceContext(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.keepPositionsBySequenceContext(se = se, 
                                                 sequence.context = "ACT"), 
                 "No sequence context found in `rowData(se)$sequence.context`",
                 fixed = TRUE)
    expect_identical(se, .keepPositionsBySequenceContext(se = se))
    
    gnm <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz", 
                                                    package = "footprintR"))
    se <- addSeqContext(se, sequence.context.width = 3, 
                        sequence.reference = gnm)
    se1 <- .keepPositionsBySequenceContext(se = se, sequence.context = "TAG")
    w <- which(as.character(rowData(se)$sequence.context) %in% c("TAG", "NNN"))
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), rownames(se)[w])
    ## Applying the same filter again doesn't do anything
    se2 <- .keepPositionsBySequenceContext(se = se1, sequence.context = "TAG")
    expect_identical(se1, se2)
})

test_that(".removeAllNAPositions works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se)
    # Subset reads to make sure there are positions with all NAs
    se <- subsetReads(se, reads = list(s1 = c(1, 2), s2 = c(1, 3)))
    
    expect_error(.removeAllNAPositions(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.removeAllNAPositions(se = se,
                                       assay.type = 1),
                 "'assay.type' must be of class 'character'")
    expect_error(.removeAllNAPositions(se = se,
                                       assay.type = c("mod_prob", "mod_prob")),
                 "'assay.type' must have length 1")
    expect_error(.removeAllNAPositions(se = se,
                                       assay.type = "missing"),
                 "'assay.type' must be one of")
    expect_error(.removeAllNAPositions(se = se,
                                       assay.type = "Nvalid"),
                 "'assay.type' must be one of")
    
    se1 <- .removeAllNAPositions(se, assay.type = "mod_prob")
    namat <- as.matrix(assay(se1, "mod_prob"))
    expect_s4_class(namat, "NaArray")
    expect_equal(ncol(namat), 4L)
    w <- which(rowSums(!is.na(as.matrix(as.matrix(assay(se, "mod_prob"))))) > 0)
    expect_equal(nrow(se1), length(w))
    expect_equal(rownames(se1), names(w))
    ## Applying the same filter again doesn't do anything
    se2 <- .removeAllNAPositions(se1, assay.type = "mod_prob")
    expect_identical(se1, se2)
})

test_that(".pruneAmbiguousStrandPositions works", {
    # generate some example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se)

    expect_error(.pruneAmbiguousStrandPositions(se = "error"),
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assay.type = 1),
                 "'assay.type' must be of class 'character'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assay.type = "missing"),
                 "'assay.type' must be one of")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assay.type = c("Nvalid", "Nmod")),
                 "'assay.type' must have length 1")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assay.type = "Nvalid",
                                                verbose = "TRUE"),
                 "'verbose' must be of class 'logical'")
    expect_error(.pruneAmbiguousStrandPositions(se = se,
                                                assay.type = "Nvalid",
                                                verbose = c(TRUE, FALSE)),
                 "'verbose' must have length 1")

    expect_equal(nrow(se), 9127L)
    expect_equal(length(unique(paste0(seqnames(rowRanges(se)),
                                      pos(rowRanges(se))))), 8955L)
    sefilt <- .pruneAmbiguousStrandPositions(se, assay.type = "Nvalid",
                                             verbose = FALSE)
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
    expect_message({
        sefilt <- .pruneAmbiguousStrandPositions(se, assay.type = "mod_prob",
                                                 verbose = TRUE)},
        "172 rows removed to ensure"
    )
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
    expect_identical(se[rownames(sefilt), ], sefilt)
    expect_message({
        sefilt <- .pruneAmbiguousStrandPositions(sefilt, assay.type = "mod_prob",
                                                 verbose = TRUE)},
        "No genomic positions represented by multiple rows found"
    )
})

test_that("filterPositions works", {
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se)
    se <- addSeqContext(se, sequence.context.width = 3, sequence.reference = reffile)
    
    expect_error(filterPositions(se = "missing"), 
                 "'se' must be of class 'SummarizedExperiment'")
    expect_error(filterPositions(se = se, filters = 1),
                 "'filters' must be of class 'character'")
    expect_error(filterPositions(se = se, filters = "error"),
                 "must be one of")

    sefilt <- filterPositions(se, c("sequence.context", "coverage",
                                    "repeated.positions", "all.na"),
                              min.cov = 5, sequence.context = "TAG")
    expect_gte(min(rowSums(assay(sefilt, "Nvalid"))), 5L)
    expect_equal(nrow(sefilt), 4934L)
    expect_true(all(as.character(rowData(sefilt)$sequence.context) %in% c("NNN", "TAG")))
    expect_false(any(duplicated(paste0(seqnames(rowRanges(sefilt)), ":",
                                       pos(rowRanges(sefilt))))))
})
