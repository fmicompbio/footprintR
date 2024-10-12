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
