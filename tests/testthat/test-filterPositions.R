test_that(".pruneAmbiguousStrandPositions works", {
    # generate some example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se, keep.reads = TRUE)

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

    expect_equal(nrow(se), 9127L)
    expect_equal(length(unique(paste0(seqnames(rowRanges(se)),
                                      pos(rowRanges(se))))), 8955L)
    sefilt <- .pruneAmbiguousStrandPositions(se, assay.type = "Nvalid")
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
    sefilt <- .pruneAmbiguousStrandPositions(se, assay.type = "mod_prob")
    expect_equal(nrow(sefilt), length(unique(paste0(seqnames(rowRanges(se)),
                                                    pos(rowRanges(se))))))
})
