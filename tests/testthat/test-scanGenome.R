test_that("genome scanning works", {
    ## example data
    modbamfiles <- system.file("extdata",
                            c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
                              "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
                            package = "footprintR")
    gnmfasta <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se0 <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                      modbase = "a", level = "summary",
                      modProbThreshold = 0.5,
                      BPPARAM = BiocParallel::SerialParam(),
                      verbose = FALSE)


    ## quantifyWindowsInRegion
    expect_error(quantifyWindowsInRegion(bamfiles = "error",
                                         region = "chr1:6940000-6955000",
                                         modbase = "a"))
    expect_error(quantifyWindowsInRegion(bamfiles = modbamfiles,
                                         region = "error",
                                         modbase = "a"))

    suppressMessages(expect_message(
        se1 <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                                       region = "chr1:6940000-6955000", modbase = "a",
                                       modProbThreshold = 0.5,
                                       windowMode = "fixed", windowSize = 24L,
                                       BPPARAM = BiocParallel::SerialParam(),
                                       verbose = TRUE)
    ))
    se1$group <- c("group1", "group1", "group2", "group2")
    expect_s4_class(se1, "RangedSummarizedExperiment")
    expect_identical(dim(se1), c(1312L, 4L))
    expect_true(all(width(SummarizedExperiment::rowRanges(se1)) == 24L))
    expect_identical(SummarizedExperiment::assayNames(se1),
                     c("Nmod", "Nvalid"))
    expect_identical(colSums(SummarizedExperiment::assay(se1, "Nvalid")),
                     c(s1 = 22587, s2 = 22587, s3 = 12878, s4 = 12878))
    ov <- findOverlaps(query = SummarizedExperiment::rowRanges(se0),
                       subject = SummarizedExperiment::rowRanges(se1))
    Nvalid <- SummarizedExperiment::assay(se0, "Nvalid")
    manualWindows <- do.call(rbind, lapply(split(queryHits(ov), rownames(se1)[subjectHits(ov)])[as.character(rownames(se1))],
                                           function(i) {
                                               colSums(Nvalid[i, , drop = FALSE])
                                           }))
    expect_identical(manualWindows, SummarizedExperiment::assay(se1, "Nvalid"))
    Nmod <- SummarizedExperiment::assay(se0, "Nmod")
    manualWindows <- do.call(rbind, lapply(split(queryHits(ov), rownames(se1)[subjectHits(ov)])[as.character(rownames(se1))],
                                           function(i) {
                                               colSums(Nmod[i, , drop = FALSE])
                                           }))
    expect_identical(manualWindows, SummarizedExperiment::assay(se1, "Nmod"))

    se2 <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                                   region = "chr1:6940000-6955000", modbase = "a",
                                   sampleAnnot = data.frame(sample = paste0("s", 1:4),
                                                            group = se1$group),
                                   windowMode = "fixed", windowSize = 24L,
                                   sequenceContextWidth = 1,
                                   sequenceReference = gnmfasta,
                                   sequenceContext = "A",
                                   BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(se2, "RangedSummarizedExperiment")
    expect_identical(dim(se2), c(1311L, 4L))
    expect_true(all(width(SummarizedExperiment::rowRanges(se2)) == 24L))
    i <- GenomicRanges::match(SummarizedExperiment::rowRanges(se2),
                              SummarizedExperiment::rowRanges(se1))
    expect_true(!any(is.na(i)))
    expect_identical(SummarizedExperiment::assayNames(se2),
                     c("Nmod", "Nvalid"))
    expect_identical(colSums(SummarizedExperiment::assay(se2, "Nvalid")),
                     c(s1 = 22241, s2 = 22241, s3 = 12358, s4 = 12358))
    expect_true(all(SummarizedExperiment::assay(se1, "Nmod")[i,] >= SummarizedExperiment::assay(se2, "Nmod")))
    expect_true(all(SummarizedExperiment::assay(se1, "Nvalid")[i,] >= SummarizedExperiment::assay(se2, "Nvalid")))
    expect_identical(SummarizedExperiment::colData(se1),
                     SummarizedExperiment::colData(se2))


    ## getDifferentiallyModifiedWindows
    expect_error(getDifferentiallyModifiedWindows(se = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, assayNameMod = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, assayNameValid = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, groupCol = "error"))
    se1$group <- c("group1", "group2", "group3", "group4")
    expect_error(getDifferentiallyModifiedWindows(se1, groupCol = "group"))
    se1$group <- c("group1", "group1", "group2", "group2")

    suppressMessages(expect_message(
        tab1 <- getDifferentiallyModifiedWindows(se1, groupCol = "group", verbose = TRUE)
    ))
    expect_s4_class(tab1, "TopTags")
    expect_identical(dim(tab1), c(nrow(se1), 10L))
    expect_identical(colnames(tab1),
                     c("seqnames", "start", "end", "width", "strand", "logFC",
                       "logCPM", "LR", "PValue", "FDR"))
    expect_s4_class(gr1 <- as(tab1$table, "GRanges"), "GRanges")

    tab2 <- getDifferentiallyModifiedWindows(se2, groupCol = "group")
    expect_s4_class(tab2, "TopTags")
    expect_identical(dim(tab2), c(nrow(se2), 10L))
    expect_identical(colnames(tab2),
                     c("seqnames", "start", "end", "width", "strand", "logFC",
                       "logCPM", "LR", "PValue", "FDR"))
    expect_s4_class(gr2 <- as(tab2$table, "GRanges"), "GRanges")
    i <- GenomicRanges::match(gr2, gr1)
    expect_true(!any(is.na(i)))
    expect_true(cor(gr1$logFC[i], gr2$logFC) > 0.98)

    ## fuseWindows
    expect_error(fuseWindows(x = "error"))
    expect_error(fuseWindows(x = gr1))
    expect_error(fuseWindows(x = gr1, scoreCol = "error"))
    expect_error(fuseWindows(x = gr1, scoreCol = "logFC"))

    suppressMessages(expect_message(
        gr1Fused <- fuseWindows(x = gr1, scoreCol = "logFC", thresh = 5.0, verbose = TRUE)
    ))
    expect_s4_class(gr1Fused, "GRanges")
    expect_length(gr1Fused, 33L)
    expect_identical(colnames(GenomicRanges::mcols(gr1Fused)),
                     c("logFCThresh", "numWindowsThresh", "direction", "logFC", "numWindows"))
    expect_identical(sum(width(gr1Fused)), 3024L)

    gr2Fused <- fuseWindows(x = gr2, scoreCol = "logFC", thresh = 5.0)
    expect_s4_class(gr2Fused, "GRanges")
    expect_length(gr2Fused, 31L)
    expect_identical(colnames(GenomicRanges::mcols(gr2Fused)),
                     c("logFCThresh", "numWindowsThresh", "direction", "logFC", "numWindows"))
    expect_identical(sum(width(gr2Fused)), 2988L)
})
