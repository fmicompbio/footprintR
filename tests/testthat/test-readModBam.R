suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(GenomicRanges)
    library(Rsamtools)
    library(Biostrings)
})

## -------------------------------------------------------------------------- ##
## Checks, readModBam
## -------------------------------------------------------------------------- ##
test_that("readModBam works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    modbamfiles <- system.file("extdata",
                               c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    names(modbamfiles) <- c("sample1", "sample2")
    extractfiles <- system.file("extdata",
                                c("modkit_extract_rc_6mA_1.tsv.gz",
                                  "modkit_extract_rc_6mA_2.tsv.gz"),
                               package = "footprintR")
    names(extractfiles) <- names(modbamfiles)
    sample_annot <- data.frame(sample = c("sample1", "sample2"),
                               group = c("group1", "group1"),
                               condition = c("cond2", "cond1"))

    # invalid arguments
    expect_error(readModBam(bamfiles = "error",
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "not all .bamfiles. exist")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = NULL, modbase = "a", nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 ".regions. must contain at least one genomic range if not in sampling mode")
    expect_error(readModBam(bamfiles = structure(unname(modbamfiles),
                                                 names = c("s1", "s1")),
                            regions = "chr1:6940000-6955000", modbase = "a",
                            nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "are not unique")
    expect_error(readModBam(bamfiles = modbamfiles, regions = "error",
                            modbase = "a", nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "Failed to get bam iterator")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "Z", nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "invalid .modbase. values")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = c("a", "a", "a"), nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "must have length")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = c(sample1 = "a", sample3 = "a"),
                            nAlnsToSample = 0,
                            BPPARAM = BiocParallel::SerialParam()),
                 "names of .modbase. and .bamfiles. don't agree")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = -1,
                            BPPARAM = BiocParallel::SerialParam()),
                 "must be between 0 and Inf")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = "error",
                            BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .numeric.")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = 0,
                            level = 1,
                            BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .character.")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = 0,
                            level = "error",
                            BPPARAM = BiocParallel::SerialParam()),
                 "must be one of")
    expect_error(
        expect_warning(
            expect_warning(readModBam(bamfiles = modbamfiles,
                                      regions = "chr1:6940000-6955000",
                                      modbase = "a", nAlnsToSample = 10,
                                      seqnamesToSampleFrom = "error",
                                      BPPARAM = BiocParallel::SerialParam()),
                           "Ignoring .regions."),
            "Ignoring unknown target name"),
        "Cannot sample 10 alignments from a total of 0")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = 0,
                            seqnamesToSampleFrom = "chr1", seqinfo = "error",
                            BPPARAM = BiocParallel::SerialParam()),
                 ".seqinfo. must be .NULL.")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000", modbase = "a",
                            BPPARAM = -1),
                 ".BPPARAM. must be of class .BiocParallelParam.")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000", modbase = "a",
                            BPPARAM = BiocParallel::SerialParam, trim = 1),
                 ".trim. must be of class .logical.")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000", modbase = "a",
                            sampleAnnot = sample_annot[1, ],
                            BPPARAM = BiocParallel::SerialParam()),
                 "Annotation information missing")
    expect_error(readModBam(
        bamfiles = modbamfiles,
        regions = "chr1:6940000-6955000", modbase = "a",
        sampleAnnot = sample_annot[, c("group", "condition")],
        BPPARAM = BiocParallel::SerialParam()),
        ".sampleAnnot. must have at least a column")
    expect_error(readModBam(bamfiles = modbamfiles,
                            regions = "chr1:6940000-6955000",
                            modbase = "a", nAlnsToSample = 5,
                            level = "summary",
                            BPPARAM = BiocParallel::SerialParam()),
                 "Read sampling is not supported")

    # expected results
    se0 <- readModkitExtract(fnames = extractfiles, modbase = "a",
                             BPPARAM = BiocParallel::SerialParam())
    reg1 <- c("chr1:6940000-6955000", "chr1:6929000-6929500")
    reg2 <- GRanges("chr1", IRanges(start = 6940000, end = 6955000))
    reg3 <- rep("chr1:6940000-6955000", 3)
    reg4 <- "chr1:6940000-6955000"
    reg5 <- c("chr1:6941000-6941001", "chr1:6928000-6928001")
    suppressMessages({
        expect_message(
            se1 <- readModBam(bamfiles = modbamfiles, regions = reg1,
                              modbase = "a", level = "read", nAlnsToSample = 0,
                              sequenceContextWidth = 1, sequenceReference = ref,
                              seqnamesToSampleFrom = "chr1", verbose = TRUE,
                              BPPARAM = BiocParallel::SerialParam())
        )
    })
    se1sum <- readModBam(bamfiles = modbamfiles, regions = reg1,
                         modbase = "a", level = "summary", nAlnsToSample = 0,
                         sequenceContextWidth = 1, sequenceReference = ref,
                         seqnamesToSampleFrom = "chr1", verbose = FALSE,
                         BPPARAM = BiocParallel::SerialParam())
    se1quick <- readModBam(bamfiles = modbamfiles, regions = reg1,
                           modbase = "a", level = "quickread", nAlnsToSample = 0,
                           sequenceContextWidth = 1, sequenceReference = ref,
                           seqnamesToSampleFrom = "chr1", verbose = FALSE,
                           BPPARAM = BiocParallel::SerialParam())
    se2 <- readModBam(bamfiles = unname(modbamfiles),
                      regions = reg2,
                      modbase = "a",
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                      verbose = FALSE)
    se2sum <- readModBam(bamfiles = unname(modbamfiles),
                         regions = reg2,
                         modbase = "a", level = "summary",
                         nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                         BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                         verbose = FALSE)
    se2quick <- readModBam(bamfiles = unname(modbamfiles),
                           regions = reg2,
                           modbase = "a", level = "quickread",
                           nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                           BPPARAM = BiocParallel::MulticoreParam(workers = 2L),
                           verbose = FALSE)
    se3 <- readModBam(bamfiles = modbamfiles,
                      regions = reg3,
                      modbase = "a",
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      BPPARAM = BiocParallel::SerialParam(),
                      verbose = FALSE)
    se3sum <- readModBam(bamfiles = modbamfiles,
                         regions = reg3,
                         modbase = "a", level = "summary",
                         nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                         BPPARAM = BiocParallel::SerialParam(),
                         verbose = FALSE)
    se3quick <- readModBam(bamfiles = modbamfiles,
                           regions = reg3,
                           modbase = "a", level = "quickread",
                           nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
    se4 <- readModBam(bamfiles = modbamfiles,
                      regions = reg4,
                      modbase = c("a", "m"),
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      BPPARAM = BiocParallel::SerialParam(),
                      verbose = FALSE)
    se4sum <- readModBam(bamfiles = modbamfiles,
                         regions = reg4, level = "summary",
                         modbase = c("a", "m"),
                         nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                         BPPARAM = BiocParallel::SerialParam(),
                         verbose = FALSE)
    se4quick <- readModBam(bamfiles = modbamfiles,
                           regions = reg4, level = "quickread",
                           modbase = c("a", "m"),
                           nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                           BPPARAM = BiocParallel::SerialParam(),
                           verbose = FALSE)
    se5a <- readModBam(bamfiles = modbamfiles[1],
                       regions = reg5[1],
                       modbase = "a",
                       nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                       BPPARAM = BiocParallel::SerialParam(),
                       verbose = FALSE)
    se5b <- readModBam(bamfiles = modbamfiles[1],
                       regions = reg5[1:2],
                       modbase = "a",
                       nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                       BPPARAM = BiocParallel::SerialParam(),
                       verbose = FALSE)
    aln5a <- Rsamtools::scanBam(file = modbamfiles[1],
                                param = Rsamtools::ScanBamParam(
        what = "qname",
        which = GRanges(reg5[1])
    ))
    aln5b <- Rsamtools::scanBam(file = modbamfiles[1],
                                param = Rsamtools::ScanBamParam(
        what = "qname",
        which = GRanges(reg5[1:2])
    ))
    expect_warning(
        suppressMessages({
            expect_message(
                se6a  <- readModBam(bamfiles = modbamfiles[1],
                                    regions = NULL,
                                    modbase = "a",
                                    nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                                    variantPositions = GPos("chr1", pos = 63000000),
                                    BPPARAM = BiocParallel::MulticoreParam(2L, RNGseed = 55L),
                                    verbose = TRUE)
            )
        })
    )
    se6b  <- readModBam(bamfiles = modbamfiles[1],
                        regions = NULL,
                        modbase = "a",
                        nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                        BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                        verbose = FALSE)
    se7  <- readModBam(bamfiles = modbamfiles,
                       regions = reg5[1:2],
                       modbase = "a",
                       sampleAnnot = sample_annot,
                       BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                       verbose = FALSE)
    se7sum  <- readModBam(bamfiles = modbamfiles,
                          regions = reg5[1:2],
                          modbase = "a", level = "summary",
                          sampleAnnot = sample_annot,
                          BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                          verbose = FALSE)
    se7quick  <- readModBam(bamfiles = modbamfiles,
                            regions = reg5[1:2],
                            modbase = "a", level = "quickread",
                            sampleAnnot = sample_annot,
                            BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                            verbose = FALSE)
    se8 <- readModBam(bamfiles = modbamfiles,
                      regions = reg1,
                      modbase = "a",
                      BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                      trim = TRUE, verbose = FALSE)
    se8sum <- readModBam(bamfiles = modbamfiles,
                         regions = reg1,
                         modbase = "a", level = "summary",
                         BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                         trim = TRUE, verbose = FALSE)
    se8quick <- readModBam(bamfiles = modbamfiles,
                           regions = reg1,
                           modbase = "a", level = "quickread",
                           BPPARAM = BiocParallel::SerialParam(RNGseed = 55L),
                           trim = TRUE, verbose = FALSE)

    seL <- list(se1, se2, se3, se4, se5a, se5b, se6a, se6b, se8)
    seLsum <- list(se1sum, se2sum, se3sum, se4sum, se8sum)
    seLquick <- list(se1quick, se2quick, se3quick, se4quick, se8quick)

    # ... structure
    expected_coldata_names <- c("sample", "modbase", "n_reads", "readInfo")
    expected_coldata_names_summary <- c("sample", "modbase")
    expected_read_info_names <- c("qscore", "read_length", "aligned_length",
                                  "variant_label", "aligned_fraction")
    for (se in c(seL, list(se7), seLquick, list(se7quick))) {
        expect_s4_class(se, "RangedSummarizedExperiment")
        expect_s4_class(rowRanges(se), "GPos")
        expect_s4_class(colData(se)$readInfo, "SimpleList")
        res_se <- lapply(colData(se)$readInfo, function(df) {
            expect_s4_class(df, "DataFrame")
            expect_named(df, expected_read_info_names)
        })
        expect_identical(assayNames(se), "mod_prob")
        expect_s4_class(assay(se, "mod_prob"), "DFrame")
        expect_s4_class(assay(se, "mod_prob")[[1]], "NaMatrix")
        expect_equal(vapply(assay(se, "mod_prob"), ncol, 0), se$n_reads,
                     ignore_attr = TRUE)
    }
    for (se in c(seL, seLquick)) {
        expect_identical(colnames(colData(se)), expected_coldata_names)
    }
    expect_identical(colnames(colData(se7)), c(expected_coldata_names,
                                               "group", "condition"))
    expect_identical(colnames(colData(se7quick)), c(expected_coldata_names,
                                                    "group", "condition"))
    ## ... ... summary-level
    for (se in c(seLsum, list(se7sum))) {
        expect_s4_class(se, "RangedSummarizedExperiment")
        expect_s4_class(rowRanges(se), "GPos")
        expect_identical(assayNames(se), c("Nmod", "Nvalid", "FracMod"))
        expect_type(assay(se, "Nmod"), "double")
    }
    for (se in seLsum) {
        expect_identical(colnames(colData(se)), expected_coldata_names_summary)
    }
    expect_identical(colnames(colData(se7sum)), c(expected_coldata_names_summary,
                                                  "group", "condition"))

    expect_identical(colnames(se1), names(modbamfiles))
    expect_identical(colnames(se1sum), names(modbamfiles))
    expect_identical(colnames(se1quick), names(modbamfiles))
    expect_identical(colnames(se2), c("s1", "s2"))
    expect_identical(colnames(se2sum), c("s1", "s2"))
    expect_identical(colnames(se2quick), c("s1", "s2"))
    expect_identical(colnames(se3), names(modbamfiles))
    expect_identical(colnames(se3sum), names(modbamfiles))
    expect_identical(colnames(se3quick), names(modbamfiles))
    expect_identical(colnames(se4), names(modbamfiles))
    expect_identical(colnames(se4sum), names(modbamfiles))
    expect_identical(colnames(se4quick), names(modbamfiles))
    expect_identical(colnames(se5a), names(modbamfiles)[1])
    expect_identical(colnames(se5b), names(modbamfiles)[1])
    expect_identical(colnames(se6a), names(modbamfiles)[1])
    expect_identical(colnames(se6b), names(modbamfiles)[1])
    expect_identical(colnames(se7), names(modbamfiles))
    expect_identical(colnames(se7sum), names(modbamfiles))
    expect_identical(colnames(se7quick), names(modbamfiles))
    expect_identical(colnames(se8), names(modbamfiles))
    expect_identical(colnames(se8sum), names(modbamfiles))
    expect_identical(colnames(se8quick), names(modbamfiles))

    # ... content se1
    expect_identical(unname(se1$n_reads), c(4L, 6L))
    expect_identical(dim(se1), c(8691L, 2L))
    modprob0 <- as.matrix(assay(se0, "mod_prob"))
    modprob1 <- as.matrix(assay(se1, "mod_prob"))
    shared_rows <- intersect(rownames(modprob0), rownames(modprob1))
    shared_cols <- intersect(colnames(modprob0), colnames(modprob1))
    expect_length(shared_rows, 8615L)
    expect_length(shared_cols, 10L)
    nonzero <- as.vector(modprob1[shared_rows, shared_cols]) > 0 &
        as.vector(modprob0[shared_rows, shared_cols]) > 0
    # plot(as.vector(modprob1[shared_rows, shared_cols])[nonzero],
    #      as.vector(modprob0[shared_rows, shared_cols])[nonzero])
    expect_equal(as.vector(modprob1[shared_rows, shared_cols])[nonzero],
                 as.vector(modprob0[shared_rows, shared_cols])[nonzero],
                 tolerance = 1e-6)
    expect_identical(colnames(se1), names(modbamfiles))
    expect_identical(lapply(se1$readInfo, rownames),
                     lapply(assay(se1, "mod_prob"), colnames))
    expect_equal(lapply(se1$readInfo, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973,
                                 21.1338005065918, 20.3082008361816),
                     sample2 = c(12.9041996002197, 9.67461013793945,
                                 15.0149002075195, 15.1365995407104,
                                 17.7175006866455, 13.6647996902466)))
    expect_identical(lapply(se1$readInfo, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 9246L, 12277L),
                         sample2 = c(13108L, 11834L, 9674L, 10047L, 8973L, 10057L)
                     ))
    expect_identical(lapply(se1$readInfo, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 9227L, 12227L),
                         sample2 = c(9656L, 11234L, 9579L, 9967L, 8915L, 9898L)
                     ))
    expect_identical(lapply(se1$readInfo, "[[", "variant_label"),
                     lapply(structure(se1$n_reads, names = colnames(se1)), function(n) rep(NA_character_, n)))
    expect_equal(unclass(table(as.character(SummarizedExperiment::rowData(se1)$sequenceContext))),
                 c(A = 8108L, C = 128L, G = 393L, T = 62L), ignore_attr = TRUE)
    # ... compare to se1sum
    expect_identical(rownames(se1), rownames(se1sum))
    se1tmp <- flattenReadLevelAssay(se1)
    expect_identical(assay(se1tmp, "Nmod"), assay(se1sum, "Nmod"))
    expect_identical(assay(se1tmp, "Nvalid"), assay(se1sum, "Nvalid"))
    expect_identical(assay(se1tmp, "FracMod"), assay(se1sum, "FracMod"))
    expect_identical(colData(se1)[, c("sample", "modbase")],
                     colData(se1sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se1), rowRanges(se1sum))
    # ... compare to se1quick
    expect_identical(se1, se1quick)
    expect_identical(rownames(se1), rownames(se1quick))
    expect_identical(rowRanges(se1), rowRanges(se1quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se1, "mod_prob"),
                     assay(se1quick, "mod_prob"))
    expect_identical(metadata(se1), metadata(se1quick))
    expect_identical(colData(se1)[, c("sample", "modbase", "n_reads")],
                     colData(se1quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se1)$readInfo$sample1),
                     rownames(colData(se1quick)$readInfo$sample1))

    # ... content se2
    expect_identical(unname(se2$n_reads), c(3L, 2L))
    expect_identical(dim(se2), dim(se3))
    expect_identical(unname(as.matrix(assay(se2, "mod_prob"))),
                     unname(as.matrix(assay(se3, "mod_prob"))))
    expect_identical(sub("^s", "sample", colnames(se2)), colnames(se3))
    expect_identical(colnames(se2), sub("sample", "s", colnames(se3)))
    for (nm in expected_read_info_names) {
        expect_equal(lapply(se2$readInfo, "[[", nm),
                     lapply(se3$readInfo, "[[", nm),
                     ignore_attr = TRUE)
    }
    # ... compare to se2sum
    expect_identical(rownames(se2), rownames(se2sum))
    se2tmp <- flattenReadLevelAssay(se2)
    expect_identical(assay(se2tmp, "Nmod"), assay(se2sum, "Nmod"))
    expect_identical(assay(se2tmp, "Nvalid"), assay(se2sum, "Nvalid"))
    expect_identical(assay(se2tmp, "FracMod"), assay(se2sum, "FracMod"))
    expect_identical(colData(se2)[, c("sample", "modbase")],
                     colData(se2sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se2), rowRanges(se2sum))
    # ... compare to se2quick
    expect_identical(se2, se2quick)
    expect_identical(rownames(se2), rownames(se2quick))
    expect_identical(rowRanges(se2), rowRanges(se2quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se2, "mod_prob"),
                     assay(se2quick, "mod_prob"))
    expect_identical(metadata(se2), metadata(se2quick))
    expect_identical(colData(se2)[, c("sample", "modbase", "n_reads")],
                     colData(se2quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se2)$readInfo$sample1),
                     rownames(colData(se2quick)$readInfo$sample1))

    # ... content se3
    expect_identical(unname(se3$n_reads), c(3L, 2L))
    expect_identical(dim(se3), c(7967L, 2L))
    modprob3 <- as.matrix(assay(se3, "mod_prob"))
    shared_rows <- intersect(rownames(modprob0), rownames(modprob3))
    shared_cols <- intersect(colnames(modprob0), colnames(modprob3))
    expect_length(shared_rows, 7924L)
    expect_length(shared_cols, 5L)
    nonzero <- as.vector(modprob3[shared_rows, shared_cols]) > 0 &
        as.vector(modprob0[shared_rows, shared_cols]) > 0
    # plot(as.vector(modprob3[shared_rows, shared_cols])[nonzero],
    #      as.vector(modprob0[shared_rows, shared_cols])[nonzero])
    expect_equal(as.vector(modprob3[shared_rows, shared_cols])[nonzero],
                 as.vector(modprob0[shared_rows, shared_cols])[nonzero],
                 tolerance = 1e-6)
    expect_identical(lapply(se3$readInfo, rownames),
                     lapply(assay(se3, "mod_prob"), colnames))
    expect_equal(lapply(se3$readInfo, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                     sample2 = c(9.67461013793945, 13.6647996902466)))
    expect_identical(lapply(se3$readInfo, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 12277L),
                         sample2 = c(11834L, 10057L)
                     ))
    expect_identical(lapply(se3$readInfo, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 12227L),
                         sample2 = c(11234L, 9898L)
                     ))
    # ... compare to se3sum
    expect_identical(rownames(se3), rownames(se3sum))
    se3tmp <- flattenReadLevelAssay(se3)
    expect_identical(assay(se3tmp, "Nmod"), assay(se3sum, "Nmod"))
    expect_identical(assay(se3tmp, "Nvalid"), assay(se3sum, "Nvalid"))
    expect_identical(assay(se3tmp, "FracMod"), assay(se3sum, "FracMod"))
    expect_identical(colData(se3)[, c("sample", "modbase")],
                     colData(se3sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se3), rowRanges(se3sum))
    # ... compare to se3quick
    expect_identical(se3, se3quick)
    expect_identical(rownames(se3), rownames(se3quick))
    expect_identical(rowRanges(se3), rowRanges(se3quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se3, "mod_prob"),
                     assay(se3quick, "mod_prob"))
    expect_identical(metadata(se3), metadata(se3quick))
    expect_identical(colData(se3)[, c("sample", "modbase", "n_reads")],
                     colData(se3quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se3)$readInfo$sample1),
                     rownames(colData(se3quick)$readInfo$sample1))

    # ... content se4
    expect_identical(unname(se4$n_reads), c(3L, 0L))
    expect_identical(dim(se4), c(4772L, 2L))
    expect_identical(dim(as.matrix(assay(se4, "mod_prob"))), c(4772L, 3L))
    expect_identical(unlist(lapply(se4$readInfo, rownames), use.names = FALSE),
                     unlist(lapply(assay(se4, "mod_prob"), colnames), use.names = FALSE))
    expect_equal(lapply(se4$readInfo, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                     sample2 = numeric(0)))
    expect_identical(lapply(se4$readInfo, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 12277L),
                         sample2 = integer(0)
                     ))
    expect_identical(lapply(se4$readInfo, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 12227L),
                         sample2 = integer(0)
                     ))
    # ... compare to se4sum
    expect_identical(rownames(se4), rownames(se4sum))
    se4tmp <- flattenReadLevelAssay(se4)
    expect_identical(assay(se4tmp, "Nmod"), assay(se4sum, "Nmod"))
    expect_identical(assay(se4tmp, "Nvalid"), assay(se4sum, "Nvalid"))
    expect_identical(assay(se4tmp, "FracMod"), assay(se4sum, "FracMod"))
    expect_identical(colData(se4)[, c("sample", "modbase")],
                     colData(se4sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se4), rowRanges(se4sum))
    # ... compare to se4quick
    expect_identical(se4, se4quick)
    expect_identical(rownames(se4), rownames(se4quick))
    expect_identical(rowRanges(se4), rowRanges(se4quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se4, "mod_prob"),
                     assay(se4quick, "mod_prob"))
    expect_identical(metadata(se4), metadata(se4quick))
    expect_identical(colData(se4)[, c("sample", "modbase", "n_reads")],
                     colData(se4quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se4)$readInfo$sample1),
                     rownames(colData(se4quick)$readInfo$sample1))

    # ... content of se5a and se5b (se5a should be a subset of se5b)
    # ... ... check ground truth
    expect_identical(names(aln5a), names(aln5b)[1])
    expect_identical(aln5a[[1]], aln5b[[1]])
    expect_length(aln5a[[1]]$qname, 1L)
    expect_length(intersect(aln5a[[1]]$qname, aln5b[[2]]$qname), 0L)
    # ... ... check return values
    mp5a <- assay(se5a, "mod_prob")
    mp5b <- assay(se5b, "mod_prob")
    expect_true(paste0("sample1-", aln5a[[1]]$qname) %in% colnames(mp5a$sample1))
    expect_true(paste0("sample1-", aln5b[[1]]$qname) %in% colnames(mp5b$sample1))
    idx <- intersect(rownames(se5a), rownames(se5b))
    expect_identical(idx, rownames(se5a))
    expect_identical(mp5a[idx, "sample1"][, paste0("sample1-", aln5a[[1]]$qname)],
                     mp5b[idx, "sample1"][, paste0("sample1-", aln5b[[1]]$qname)])

    # ... content of se6a and se6b
    expect_identical(se6a, se6b)
    expect_identical(dim(assay(se6a)$sample1), c(5996L, 5L))

    # ... content of se7
    mp7 <- assay(se7, "mod_prob")
    expect_true(paste0("sample1-", aln5b[[1]]$qname) %in% colnames(mp7$sample1))
    idx <- rownames(se7)
    expect_identical(mp7[idx, "sample1"][, paste0("sample1-", aln5b[[1]]$qname)],
                     mp5b[idx, "sample1"][, paste0("sample1-", aln5b[[1]]$qname)])
    # ... compare to se7sum
    expect_identical(rownames(se7), rownames(se7sum))
    se7tmp <- flattenReadLevelAssay(se7)
    expect_identical(assay(se7tmp, "Nmod"), assay(se7sum, "Nmod"))
    expect_identical(assay(se7tmp, "Nvalid"), assay(se7sum, "Nvalid"))
    expect_identical(assay(se7tmp, "FracMod"), assay(se7sum, "FracMod"))
    expect_identical(colData(se7)[, c("sample", "modbase")],
                     colData(se7sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se7), rowRanges(se7sum))
    # ... compare to se7quick
    expect_identical(se7, se7quick)
    expect_identical(rownames(se7), rownames(se7quick))
    expect_identical(rowRanges(se7), rowRanges(se7quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se7, "mod_prob"),
                     assay(se7quick, "mod_prob"))
    expect_identical(metadata(se7), metadata(se7quick))
    expect_identical(colData(se7)[, c("sample", "modbase", "n_reads")],
                     colData(se7quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se7)$readInfo$sample1),
                     rownames(colData(se7quick)$readInfo$sample1))

    # ... content of se8 (like se1, but trimmed)
    expect_identical(unname(se8$n_reads), c(4L, 6L))
    expect_identical(dim(se8), c(1009L, 2L))
    modprob1 <- as.matrix(assay(se1, "mod_prob"))
    modprob8 <- as.matrix(assay(se8, "mod_prob"))
    shared_rows <- intersect(rownames(modprob8), rownames(modprob1))
    shared_cols <- intersect(colnames(modprob8), colnames(modprob1))
    expect_length(shared_rows, 1009L)
    expect_length(shared_cols, 10L)
    expect_identical(modprob1[shared_rows, ], modprob8[shared_rows, ])
    expect_identical(colnames(se8), names(modbamfiles))
    expect_identical(lapply(se8$readInfo, rownames),
                     lapply(assay(se8, "mod_prob"), colnames))
    expect_identical(se1$readInfo, se8$readInfo)
    # ... compare to se8sum
    expect_identical(rownames(se8), rownames(se8sum))
    se8tmp <- flattenReadLevelAssay(se8)
    expect_identical(assay(se8tmp, "Nmod"), assay(se8sum, "Nmod"))
    expect_identical(assay(se8tmp, "Nvalid"), assay(se8sum, "Nvalid"))
    expect_identical(assay(se8tmp, "FracMod"), assay(se8sum, "FracMod"))
    expect_identical(colData(se8)[, c("sample", "modbase")],
                     colData(se8sum)[, c("sample", "modbase")])
    expect_identical(rowRanges(se8), rowRanges(se8sum))
    # ... compare to se8quick
    expect_identical(se8, se8quick)
    expect_identical(rownames(se8), rownames(se8quick))
    expect_identical(rowRanges(se8), rowRanges(se8quick))
    ## in principle, there is no guarantee that the reads have to be in the
    ## same order (but here they are)
    expect_identical(assay(se8, "mod_prob"),
                     assay(se8quick, "mod_prob"))
    expect_identical(metadata(se8), metadata(se8quick))
    expect_identical(colData(se8)[, c("sample", "modbase", "n_reads")],
                     colData(se8quick)[, c("sample", "modbase", "n_reads")])
    expect_identical(rownames(colData(se8)$readInfo$sample1),
                     rownames(colData(se8quick)$readInfo$sample1))
})

test_that("readModBam correctly labels reads", {
    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")

    # extract alignments
    alns <- Rsamtools::scanBam(
        file = modbamfile,
        param = ScanBamParam(which = GRanges("chr1", IRanges(1, 1e8)),
                             what = c("qname", "rname", "strand", "pos",
                                      "qwidth", "cigar", "seq")))[[1]]

    # extract first 2 and last 3 positions of each alignment
    rnames <- rep(as.character(alns$rname), each = 5L)
    # ... this was calculated using: GenomicAlignments::cigarWidthAlongReferenceSpace(alns$cigar)
    alnwidth <- c(14973L, 11303L, 9254L, 12288L, 10047L, 9066L, 9052L, 8346L,
                  7678L, 6985L)
    rpos <- unlist(lapply(seq_along(alns$qname), function(i) {
        alns$pos[i] + c(0:1, alnwidth[i] - 3:1)
    }))
    # ... add three positions that don't overlap any read (one before, two after)
    rnames <- c("chr1", rnames, "chr2", "chr1")
    rpos <- c(6925829L, rpos, 6941630L, 6941639L)
    # ... convert to GPos
    varpos <- GPos(seqnames = rnames, pos = rpos,
                   names = c("miss1", rep(alns$qname, each = 5), "miss2", "miss3"))

    # calculate expected labels
    softmaskStart <- suppressWarnings(
        ifelse(grepl("^[0-9]+S", alns$cigar),
               as.integer(sub("^([0-9]+)S.+$", "\\1", alns$cigar)),
               0L))
    softmaskEnd <- suppressWarnings(
        ifelse(grepl("[0-9]+S$", alns$cigar),
               as.integer(sub("^.+?([0-9]+)S$", "\\1", alns$cigar)),
               0L))
    readLabelParts <- paste0(subseq(x = alns$seq, start = softmaskStart + 1L, width = 2L),
                             subseq(x = alns$seq, end = width(alns$seq) - softmaskEnd, width = 3L))

    # run readModBam
    se <- readModBam(bamfiles = modbamfile, modbase = "a", regions = varpos,
                     variantPositions = varpos,
                     BPPARAM = BiocParallel::SerialParam())
    varposToSortedIdx <- match(varpos, metadata(se)$variantPositions)

    # compare to expected labels
    extractLabelParts <- unlist(lapply(seq_along(alns$qname), function(i) {
        rid <- alns$qname[i]
        idx <- varposToSortedIdx[mcols(varpos)$names == rid]
        paste(
            strsplit(se$readInfo$s1$variant_label[match(paste0("s1-", rid),
                                                         rownames(se$readInfo$s1))],
                     "")[[1]][idx],
            collapse = "")
    }))
    expect_identical(sort(sort(varpos), ignore.strand = TRUE), metadata(se)$variantPositions)
    expect_identical(extractLabelParts, readLabelParts)
    expect_true(all(grepl("^-.*--$", se$readInfo$s1$variant_label))) # missed positions

    # positions that are known to be variables across reads
    varpos2 <- GPos(seqnames = "chr1", pos = c(6937731, 6937788, 6937843, 6937857,
                                               6937873, 6937931, 6937932, 6938070,
                                               6938109))

    se2 <- readModBam(bamfiles = modbamfile, modbase = "a", regions = varpos2,
                      variantPositions = varpos2,
                      BPPARAM = BiocParallel::SerialParam())
    bases <- c("A", "C", "G", "T", "-")
    expCnt <- matrix(
        as.integer(c(0, 8, 0, 2, 0,
                     0, 8, 0, 2, 0,
                     8, 0, 2, 0, 0,
                     8, 0, 2, 0, 0,
                     2, 0, 7, 0, 1,
                     0, 2, 7, 0, 1,
                     2, 0, 7, 0, 1,
                     2, 0, 7, 0, 1,
                     0, 2, 7, 0, 1)),
        ncol = length(bases), byrow = TRUE, dimnames = list(NULL, bases))
    obsCnt <- do.call(rbind, lapply(seq.int(9), function(i) {
        f <- factor(
            unlist(lapply(se2$readInfo$s1$variant_label, substr, i, i)),
            levels = bases
        )
        unclass(table(f))
    }))
    expect_identical(obsCnt, expCnt)
})
