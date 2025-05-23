test_that(".calcDirChiSqP works", {
    expect_equal(.calcDirChiSqP(0, 0, 0, 0), NaN)
    expect_equal(.calcDirChiSqP(15, 20, 125, 107), -0.6643047)
    expect_error(.calcDirChiSqP(15, 20, 7, 12), "must be nonnegative")
})

test_that("genome scanning works (helper functions)", {
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

    ## .tileChromosome
    rg <- .tileChromosome(tileSize = 40, windowSize = 12,
                          windowStep = 1, chromName = "chr1", chromLength = 100)
    expect_s4_class(rg, "GRanges")
    expect_equal(start(rg), c(1, 41, 81))
    expect_equal(width(rg), rep(51, 3))

    rg <- .tileChromosome(tileSize = 40, windowSize = 12,
                          windowStep = 3, chromName = "chr1", chromLength = 100)
    expect_s4_class(rg, "GRanges")
    expect_equal(start(rg), c(1, 40, 79))
    expect_equal(width(rg), rep(48, 3))

    rg <- .tileChromosome(tileSize = 40, windowSize = 12,
                          windowStep = 6, chromName = "chr1", chromLength = 100)
    expect_s4_class(rg, "GRanges")
    expect_equal(start(rg), c(1, 37, 73))
    expect_equal(width(rg), rep(42, 3))

    rg <- .tileChromosome(tileSize = 20, windowSize = 12,
                          windowStep = 6, chromName = "chr1", chromLength = 100)
    expect_s4_class(rg, "GRanges")
    expect_equal(start(rg), c(1, 19, 37, 55, 73, 91))
    expect_equal(width(rg), rep(24, 6))

    rg <- .tileChromosome(tileSize = 20, windowSize = 12,
                          windowStep = 12, chromName = "chr1", chromLength = 100)
    expect_s4_class(rg, "GRanges")
    expect_equal(start(rg), c(1, 13, 25, 37, 49, 61, 73, 85, 97))
    expect_equal(width(rg), rep(12, 9))

    ## strandDiffFracMod
    se1 <- se0
    SummarizedExperiment::assayNames(se1) <- c("a", "b", "c")
    expect_error(strandDiffFracMod(se1, rg), ".se. must contain assays")
    rg <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(
            start = seq(6930000, 6941000, by = 1000),
            width = 1000,
            names = letters[seq.int(12)]))
    res <- strandDiffFracMod(se0, rg)
    expect_s4_class(res, "RangedSummarizedExperiment")
    expect_identical(dim(res), c(length(rg), ncol(se0)))
    expect_identical(SummarizedExperiment::rowRanges(res), rg)
    expect_identical(SummarizedExperiment::assayNames(res),
                     c("Nmodpos", "Nmodneg", "Nvalidpos", "Nvalidneg",
                       "dirNegLog10PValue", "FracModDiff"))
    sel <- IRanges::overlapsAny(SummarizedExperiment::rowRanges(se0), rg)
    cnt0 <- vapply(SummarizedExperiment::assays(se0)[c("Nmod", "Nvalid")],
                   function(a) sum(a[sel,]), 0.0)
    cnt1 <- vapply(SummarizedExperiment::assays(res)[c("Nmodpos", "Nmodneg", "Nvalidpos", "Nvalidneg")],
                   function(a) sum(a), 0.0)
    expect_identical(SummarizedExperiment::assay(res, "Nmodpos")[, "s3"],
                     setNames(c(11, 37, 58, 35, 32, 14, 42, 37, 49, 32, 24, 0),
                              letters[1:12]))
    expect_identical(cnt0[["Nmod"]], cnt1[["Nmodpos"]] + cnt1[["Nmodneg"]])
    expect_identical(cnt0[["Nvalid"]], cnt1[["Nvalidpos"]] + cnt1[["Nvalidneg"]])
    res <- strandDiffFracMod(se0, unname(rg))
    expect_identical(rownames(res), as.character(seq_along(rg)))
    res <- strandDiffFracMod(se0[numeric(0), ], rg)
    expect_s4_class(res, "RangedSummarizedExperiment")
    expect_identical(dim(res), c(0L, ncol(se0)))

    ## sumNmodNvalid
    rg <- .tileChromosome(tileSize = 1e6, windowSize = 1e6,
                          windowStep = 1e6, chromName = "chr1",
                          chromLength = 70e6)
    se1 <- se0
    SummarizedExperiment::assayNames(se1) <- c("a", "b", "c")
    expect_error(sumNmodNvalid(se1, rg), ".se. must contain assays")
    res <- sumNmodNvalid(se0, rg, includeEmpty = FALSE)
    expect_identical(SummarizedExperiment::assay(res, "Nmod")[1, ],
                     colSums(SummarizedExperiment::assay(se0, "Nmod")))
    expect_identical(SummarizedExperiment::assay(res, "Nvalid")[1, ],
                     colSums(SummarizedExperiment::assay(se0, "Nvalid")))

    # ... includeEmpty = TRUE
    rg <- .tileChromosome(tileSize = 1e6, windowSize = 1e6,
                          windowStep = 1e6, chromName = "chr1",
                          chromLength = 70e6)
    se1 <- se0
    SummarizedExperiment::assayNames(se1) <- c("a", "b", "c")
    res <- sumNmodNvalid(se0, rg, includeEmpty = TRUE)
    expect_identical(SummarizedExperiment::assay(res, "Nmod")[7, ],
                     colSums(SummarizedExperiment::assay(se0, "Nmod")))
    expect_identical(SummarizedExperiment::assay(res, "Nvalid")[7, ],
                     colSums(SummarizedExperiment::assay(se0, "Nvalid")))

    ## phasingScoreFourier
    rg <- .tileChromosome(tileSize = 1e6, windowSize = 1e6,
                          windowStep = 1e6, chromName = "chr1",
                          chromLength = 70e6)
    se1 <- se0
    SummarizedExperiment::assayNames(se1) <- c("a", "b", "c")
    expect_error(phasingScoreFourier(se1, rg), ".se. must contain assays")
    expect_error(phasingScoreFourier(
        se0, GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 1:2, width = 3:4))),
        "have the same width")
    expect_error(phasingScoreFourier(
        se0, GenomicRanges::GRanges("chr1", IRanges::IRanges(start = 1:2, width = 12)),
        numCoef = 5),
        "is not divisible")
    expect_error(phasingScoreFourier(se0, rg[c(1, 3, 4)]),
                 "need to be regularly spaced")
    windowgr <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(start = seq(0, 11) * 190 + 6930001,
                                  width = 4 * 190))
    res0a <- phasingScoreFourier(se = se0[numeric(0), ], gr = windowgr, numCoef = 5)
    res0b <- phasingScoreFourier(se = se0, gr = GenomicRanges::GRanges(), numCoef = 5)
    res1 <- phasingScoreFourier(se = IRanges::subsetByOverlaps(se0, windowgr),
                                gr = windowgr, numCoef = 5)
    res1a <- phasingScoreFourier(se = IRanges::subsetByOverlaps(se0, windowgr[7]),
                                 gr = windowgr[7], numCoef = 5)
    expect_identical(dim(res0a), c(0L, BiocGenerics::ncol(se0)))
    expect_identical(dim(res0b), c(0L, BiocGenerics::ncol(se0)))
    expect_identical(dim(res1), c(length(windowgr), BiocGenerics::ncol(se0)))
    expect_equal(assay(res1, "phasingScoreAbs")[7, , drop = FALSE],
                 assay(res1a, "phasingScoreAbs"), tolerance = 1e-4)
    expect_equal(assay(res1, "phasingScoreRel")[7, , drop = FALSE],
                 assay(res1a, "phasingScoreRel"), tolerance = 1e-4)
    expect_identical(SummarizedExperiment::assayNames(res0a),
                     c("phasingScoreAbs", "phasingScoreRel"))
    expect_identical(SummarizedExperiment::assayNames(res0b),
                     c("phasingScoreAbs", "phasingScoreRel"))
    expect_identical(SummarizedExperiment::assayNames(res1),
                     c("phasingScoreAbs", "phasingScoreRel"))
    expect_identical(SummarizedExperiment::rowRanges(res1), windowgr)
    ass1 <- SummarizedExperiment::assay(res1, "phasingScoreAbs")
    ass2 <- SummarizedExperiment::assay(res1, "phasingScoreRel")
    expect_type(ass1, "double")
    expect_identical(ass1[, 1], ass1[, 2])
    expect_identical(ass1[, 3], ass1[, 4])
    expect_true(all(ass1 >= 0))
    expect_type(ass2, "double")
    expect_identical(ass2[, 1], ass2[, 2])
    expect_identical(ass2[, 3], ass2[, 4])
    expect_true(all(ass2 >= 0 & ass2 <= 1.0))
    corr1 <- cor(ass1[, c(1, 3)], ass2[, c(1, 3)])
    expect_equal(which.max(corr1[, 1]), 1L, ignore_attr = TRUE)
    expect_equal(which.max(corr1[, 2]), 2L, ignore_attr = TRUE)
    se4 <- res1

    ## estimateNRLwindows
    rng <- GenomicRanges::GRanges("chr1", IRanges::IRanges(6930000, 6940000))
    windowSize <- 2000
    windowStep <- 1000
    s <- seq(start(rng), end(rng) - windowSize + 1, by = windowStep)
    windowgr <- GenomicRanges::GRanges(
        seqnames = seqnames(rng),
        ranges = IRanges::IRanges(start = s, width = windowSize))
    se1 <- readModBam(bamfiles = modbamfiles, regions = rng, level = "quickread",
                      modbase = "a", trim = TRUE,
                      BPPARAM = BiocParallel::SerialParam())
    seEmpty1 <- estimateNRLwindows(se = se1, gr = GRanges())
    expect_s4_class(seEmpty1, "RangedSummarizedExperiment")
    expect_identical(dim(seEmpty1), c(0L, length(modbamfiles)))
    expect_identical(SummarizedExperiment::assayNames(seEmpty1),
                     c("NRL", "NRL.CI95low", "NRL.CI95high"))
    gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1:2, width = 2000))
    seEmpty2 <- estimateNRLwindows(se = se1, gr = gr2)
    expect_s4_class(seEmpty2, "RangedSummarizedExperiment")
    expect_identical(dim(seEmpty2), c(length(gr2), length(modbamfiles)))
    expect_identical(SummarizedExperiment::assayNames(seEmpty2),
                     c("NRL", "NRL.CI95low", "NRL.CI95high"))
    expect_true(all(is.na(SummarizedExperiment::assay(seEmpty2, "NRL"))))
    seNRL <- estimateNRLwindows(se = se1, gr = windowgr)
    expect_s4_class(seNRL, "RangedSummarizedExperiment")
    expect_equal(assay(seNRL, "NRL")[, c(1,3)], assay(seNRL, "NRL")[, c(2,4)],
                 ignore_attr = TRUE)
    expect_true(all(assay(seNRL, "NRL") > 100))
    expect_true(all(assay(seNRL, "NRL") < 300))
    expect_true(all(assay(seNRL, "NRL.CI95low") < assay(seNRL, "NRL.CI95high")))

    ## quantifyWindowsInRegion
    expect_error(quantifyWindowsInRegion(bamfiles = "error",
                                         region = "chr1:6940000-6955000",
                                         modbase = "a"))
    expect_error(quantifyWindowsInRegion(bamfiles = modbamfiles,
                                         region = "error",
                                         modbase = "a"))
    expect_error(quantifyWindowsInRegion(bamfiles = modbamfiles,
                                         region = "chr1:6940000-6955000",
                                         modbase = "a",
                                         quantFunction = "error"),
                 ".quantFunction. must be the name of an existing function")
    expect_identical(dim(quantifyWindowsInRegion(bamfiles = modbamfiles,
                                                 region = "chr1:1-1000",
                                                 modbase = "a",
                                                 BPPARAM = BiocParallel::SerialParam())),
                     c(0L, 4L))

    # sumNmodNvalid, includeEmpty = FALSE
    suppressMessages(expect_message(
        se1 <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                                       region = "chr1:6940000-6955000", modbase = "a",
                                       modProbThreshold = 0.5,
                                       windowMode = "fixed", windowSize = 24L,
                                       BPPARAM = BiocParallel::SerialParam(),
                                       verbose = TRUE,
                                       quantFunction = "sumNmodNvalid",
                                       quantFunctionArgs = list(includeEmpty = FALSE))
    ))
    se1$group <- c("group1", "group1", "group2", "group2")
    expect_s4_class(se1, "RangedSummarizedExperiment")
    expect_identical(dim(se1), c(136L, 4L))
    expect_true(all(width(SummarizedExperiment::rowRanges(se1)) == 24L))
    expect_identical(SummarizedExperiment::assayNames(se1),
                     c("Nmod", "Nvalid", "FracMod"))
    expect_identical(colSums(SummarizedExperiment::assay(se1, "Nvalid")),
                     c(s1 = 1610, s2 = 1610, s3 = 554, s4 = 554))
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

    # sumNmodNvalid, includeEmpty = TRUE
    suppressMessages(expect_message(
        se1b <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                                        region = "chr1:6940000-6955000", modbase = "a",
                                        modProbThreshold = 0.5,
                                        windowMode = "fixed", windowSize = 24L,
                                        BPPARAM = BiocParallel::SerialParam(),
                                        verbose = TRUE,
                                        quantFunction = "sumNmodNvalid",
                                        quantFunctionArgs = list(includeEmpty = TRUE))
    ))
    se1b$group <- c("group1", "group1", "group2", "group2")
    expect_s4_class(se1b, "RangedSummarizedExperiment")
    expect_identical(dim(se1b), c(1249L, 4L))
    expect_identical(sum(rowSums(assay(se1b, "Nvalid")) > 0), 136L)
    expect_true(all(width(SummarizedExperiment::rowRanges(se1b)) == 24L))
    expect_identical(SummarizedExperiment::assayNames(se1b),
                     c("Nmod", "Nvalid", "FracMod"))
    expect_identical(colSums(SummarizedExperiment::assay(se1b, "Nvalid")),
                     c(s1 = 1610, s2 = 1610, s3 = 554, s4 = 554))
    ov <- findOverlaps(query = SummarizedExperiment::rowRanges(se0),
                       subject = SummarizedExperiment::rowRanges(se1b))
    Nvalid <- SummarizedExperiment::assay(se0, "Nvalid")
    manualWindows <- do.call(rbind, lapply(split(queryHits(ov), rownames(se1b)[subjectHits(ov)])[as.character(rownames(se1b))],
                                           function(i) {
                                               colSums(Nvalid[i, , drop = FALSE])
                                           }))
    rownames(manualWindows) <- seq.int(nrow(manualWindows))
    expect_identical(manualWindows, SummarizedExperiment::assay(se1b, "Nvalid"))
    Nmod <- SummarizedExperiment::assay(se0, "Nmod")
    manualWindows <- do.call(rbind, lapply(split(queryHits(ov), rownames(se1b)[subjectHits(ov)])[as.character(rownames(se1b))],
                                           function(i) {
                                               colSums(Nmod[i, , drop = FALSE])
                                           }))
    rownames(manualWindows) <- seq.int(nrow(manualWindows))
    expect_identical(manualWindows, SummarizedExperiment::assay(se1b, "Nmod"))

    # filter by sequence context
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
    expect_identical(dim(se2), c(136L, 4L))
    expect_true(all(width(SummarizedExperiment::rowRanges(se2)) == 24L))
    i <- GenomicRanges::match(SummarizedExperiment::rowRanges(se2),
                              SummarizedExperiment::rowRanges(se1))
    expect_true(!any(is.na(i)))
    expect_identical(SummarizedExperiment::assayNames(se2),
                     c("Nmod", "Nvalid", "FracMod"))
    expect_identical(colSums(SummarizedExperiment::assay(se2, "Nvalid")),
                     c(s1 = 1590, s2 = 1590, s3 = 540, s4 = 540))
    expect_true(all(SummarizedExperiment::assay(se1, "Nmod")[i,] >= SummarizedExperiment::assay(se2, "Nmod")))
    expect_true(all(SummarizedExperiment::assay(se1, "Nvalid")[i,] >= SummarizedExperiment::assay(se2, "Nvalid")))
    expect_identical(SummarizedExperiment::colData(se1),
                     SummarizedExperiment::colData(se2))

    # predefined windows
    set.seed(1L)
    selwindows <- sample(nrow(se2), 40)
    se3 <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                                   region = "chr1:6940000-6955000", modbase = "a",
                                   sampleAnnot = data.frame(sample = paste0("s", 1:4),
                                                            group = se1$group),
                                   windowMode = "predefined",
                                   windows = SummarizedExperiment::rowRanges(se2)[selwindows],
                                   sequenceContextWidth = 1,
                                   sequenceReference = gnmfasta,
                                   sequenceContext = "A",
                                   BPPARAM = BiocParallel::SerialParam())
    expect_identical(se2[selwindows,], se3)


    ## getDifferentiallyModifiedWindows
    expect_error(getDifferentiallyModifiedWindows(se = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, assayNameMod = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, assayNameValid = "error"))
    expect_error(getDifferentiallyModifiedWindows(se = se1, groupCol = "error"))
    expect_length(getDifferentiallyModifiedWindows(se = se1[numeric(0), ]), 0L)
    se1$group <- c("group1", "group2", "group3", "group4")
    expect_error(getDifferentiallyModifiedWindows(se1, groupCol = "group"))
    se1$group <- c("group1", "group1", "group2", "group2")

    suppressMessages(expect_message(
        gr1 <- getDifferentiallyModifiedWindows(se1, groupCol = "group", verbose = TRUE)
    ))
    expect_s4_class(gr1, "GRanges")
    expect_length(gr1, 136L)
    expect_identical(ncol(GenomicRanges::mcols(gr1)), 9L)
    expect_identical(colnames(GenomicRanges::mcols(gr1)),
                     c("logFC", "logCPM", "LR", "PValue", "FDR", "dirNegLog10PValue",
                       "FracMod_group1", "FracMod_group2", "DeltaFracMod"))

    gr2 <- getDifferentiallyModifiedWindows(se2, groupCol = "group")
    expect_s4_class(gr2, "GRanges")
    expect_length(gr2, 136L)
    expect_identical(ncol(GenomicRanges::mcols(gr2)), 9L)
    expect_identical(colnames(GenomicRanges::mcols(gr2)),
                     c("logFC", "logCPM", "LR", "PValue", "FDR", "dirNegLog10PValue",
                       "FracMod_group1", "FracMod_group2", "DeltaFracMod"))
    i <- GenomicRanges::match(gr2, gr1)
    expect_true(!any(is.na(i)))
    expect_true(cor(gr1$logFC[i], gr2$logFC) > 0.98)

    ## getDifferentialWindows
    se4$group <- c("cond1", "cond1", "cond2", "cond2")
    dsgn <- stats::model.matrix(~ group, data = SummarizedExperiment::colData(se4))
    cntr <- c(0, 1)
    expect_error(getDifferentialWindows(se = "error"))
    expect_error(getDifferentialWindows(se = se4, assayName = "error"))
    expect_error(getDifferentialWindows(se = se4, designMatrix = "error"))
    expect_error(getDifferentialWindows(se = se4, designMatrix = dsgn,
                                        contrast = "error"))
    expect_error(getDifferentialWindows(se = se4, designMatrix = dsgn,
                                        contrast = cntr, method = "error"))
    expect_error(getDifferentialWindows(se = se4, designMatrix = dsgn,
                                        contrast = cntr, method = "limma",
                                        verbose = "error"))
    res0 <- getDifferentialWindows(se = se4[numeric(0), ],
                                   designMatrix =  dsgn, contrast = cntr)
    expect_s4_class(res0, "GRanges")
    expect_length(res0, 0L)
    expect_named(GenomicRanges::mcols(res0),
                 c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
                   "dirNegLog10PValue"))
    res1a <- getDifferentialWindows(se = se4, designMatrix =  dsgn,
                                    contrast = cntr, method = "limma")
    expect_warning(res1b <- getDifferentialWindows(se = se4, designMatrix =  dsgn,
                                                   contrast = cntr, method = "edgeR"))
    expect_s4_class(res1a, "GRanges")
    expect_s4_class(res1b, "GRanges")
    expect_length(res1a, nrow(se4))
    expect_identical(GenomicRanges::ranges(res1a),
                     GenomicRanges::ranges(SummarizedExperiment::rowRanges(se4)))
    expect_identical(GenomicRanges::ranges(res1a),
                     GenomicRanges::ranges(res1b))
    expect_identical(ncol(GenomicRanges::mcols(gr1)), 9L)
    expect_named(GenomicRanges::mcols(res1a),
                 c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
                   "dirNegLog10PValue"))
    expect_named(GenomicRanges::mcols(res1b),
                 c("logFC", "logCPM", "LR", "PValue", "FDR", "dirNegLog10PValue"))
    expect_true(cor(res1a$logFC, res1b$logFC) > 0.9)

    ## getRangesWithAssayValues
    resL <- list(getRangesWithAssayValues(se0),
                 getRangesWithAssayValues(se0, "Nmod"),
                 getRangesWithAssayValues(se0, "FracMod"))
    tst <- lapply(resL, \(x) expect_s4_class(x, class = "GRanges"))
    tst <- lapply(resL, \(x) expect_length(x, nrow(se0)))
    tst <- lapply(resL, \(x) {
        xx <- x
        GenomicRanges::mcols(xx) <- NULL
        expect_identical(xx, SummarizedExperiment::rowRanges(se0))
    })
    tmpNmod <- SummarizedExperiment::assay(se0, "Nmod")
    colnames(tmpNmod) <- paste0("Nmod.", colnames(tmpNmod))
    expect_identical(as.matrix(GenomicRanges::mcols(resL[[1]])),
                     tmpNmod)
    expect_identical(as.matrix(GenomicRanges::mcols(resL[[2]])),
                     tmpNmod)
    tmpFracMod <- SummarizedExperiment::assay(se0, "FracMod")
    colnames(tmpFracMod) <- paste0("FracMod.", colnames(tmpFracMod))
    expect_identical(as.matrix(GenomicRanges::mcols(resL[[3]])),
                     tmpFracMod)

    ## processWindowScores
    expect_error(processWindowScores(x = "error"))
    expect_error(processWindowScores(x = gr1, scoreCol = "error"))
    expect_length(processWindowScores(x = gr1[numeric(0)]), 0L)

    # ... pass
    expect_identical(processWindowScores(x = gr1, scoreCol = "logFC", scoreAction = "pass"), gr1)

    # ... select
    grs <- processWindowScores(x = gr1, scoreCol = "logFC", scoreAction = "select")
    expect_identical(gr1$logFC, grs$logFC)
    expect_identical(colnames(mcols(grs)), "logFC")

    # ... smooth
    gr1smooth <- processWindowScores(x = gr1, scoreCol = "logFC", scoreAction = "smooth")
    expect_identical(IRanges::ranges(gr1), IRanges::ranges(gr1smooth))
    expect_identical(GenomicRanges::mcols(gr1)[, -1], GenomicRanges::mcols(gr1smooth)[, -1])
    expect_true(cor(GenomicRanges::mcols(gr1)[, 1], GenomicRanges::mcols(gr1smooth)[, 1]) > 0.9)

    # ... smoothFuse
    suppressMessages(expect_message(
        gr1Fused <- processWindowScores(x = gr1, scoreCol = "logFC", thresh = 5.0, verbose = TRUE)
    ))
    expect_s4_class(gr1Fused, "GRanges")
    expect_length(gr1Fused, 5L)
    expect_identical(colnames(GenomicRanges::mcols(gr1Fused)),
                     c("logFCThresh", "numWindowsThresh", "direction", "logFC", "numWindows"))
    expect_identical(sum(width(gr1Fused)), 732L)

    gr2Fused <- processWindowScores(x = gr2, scoreCol = "logFC", thresh = 5.0)
    expect_s4_class(gr2Fused, "GRanges")
    expect_length(gr2Fused, 5L)
    expect_identical(colnames(GenomicRanges::mcols(gr2Fused)),
                     c("logFCThresh", "numWindowsThresh", "direction", "logFC", "numWindows"))
    expect_identical(sum(width(gr2Fused)), 732L)
})

test_that("genome scanning works (wrapper function)", {
    # example data
    modbamfiles <- system.file("extdata",
                               c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
                                 "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    annotdf <- data.frame(sample = c("s1","s2","s3","s4"),
                          group = c("A","A","B","B"))
    chrlen <- c(chr1 = 6955000)

    expect_error(scanForHighScoringRegions(bamfiles = modbamfiles,
                                           sampleAnnot = annotdf,
                                           chromosomeLengths = 1e6))
    expect_error(scanForHighScoringRegions(bamfiles = modbamfiles,
                                           sampleAnnot = annotdf,
                                           chromosomeLengths = chrlen,
                                           scoreFunction = "error"))

    gr <- scanForHighScoringRegions(bamfiles = modbamfiles,
                                    sampleAnnot = annotdf,
                                    chromosomeLengths = chrlen,
                                    modbase = "a", BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(gr, "GRanges")
    expect_length(gr, 38L)

    gr2 <- scanForHighScoringRegions(bamfiles = modbamfiles,
                                     sampleAnnot = annotdf,
                                     chromosomeLengths = chrlen,
                                     scoreCol = c("dirNegLog10PValue", "logFC"),
                                     modbase = "a", BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(gr2$dirNegLog10PValue, "GRanges")
    expect_length(gr2$dirNegLog10PValue, 38L)
    expect_identical(gr, gr2$dirNegLog10PValue)
})
