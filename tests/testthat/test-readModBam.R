suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(GenomicRanges)
    library(Rsamtools)
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

    # invalid arguments
    expect_error(readModBam("error", "chr1:6940000-6955000", "a", 0),
                 "not all `bamfiles` exist")
    expect_error(readModBam(structure(unname(modbamfiles), names = c("s1", "s1")),
                            "chr1:6940000-6955000", "a", 0),
                 "are not unique")
    expect_error(readModBam(modbamfiles, "error", "a", 0),
                 "GRanges object must contain")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "Z", 0),
                 "invalid `modbase` values")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", c("a", "a", "a"), 0),
                 "must have length")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000",
                            c(sample1 = "a", sample3 = "a"), 0),
                 "names of `modbase` and `bamfiles` don't agree")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "a", -1),
                 "must be within .0,Inf.")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "a", "error"),
                 "must be of class 'numeric'")
    expect_error(
        expect_warning(
            expect_warning(readModBam(modbamfiles, "chr1:6940000-6955000", "a", 10, "error"),
                           "Ignoring `regions`"),
            "Ignoring unknown target name"),
        "Cannot sample 10 alignments from a total of 0")
    expect_error(readModBam(modbamfiles, "chr1:6940000-6955000", "a", 0, "chr1", "error"),
                 "`seqinfo` must be `NULL`, a `Seqinfo` object or")

    # expected results
    se0 <- readModkitExtract(fnames = extractfiles, modbase = "a")
    reg1 <- c("chr1:6940000-6955000", "chr1:6929000-6929500")
    reg2 <- GRanges("chr1", IRanges(start = 6940000, end = 6955000))
    reg3 <- rep("chr1:6940000-6955000", 3)
    reg4 <- "chr1:6940000-6955000"
    reg5 <- c("chr1:6941000-6941001", "chr1:6928000-6928001")
    expect_message(expect_message(expect_message(expect_message(
        expect_message(expect_message(expect_message(expect_message(
            expect_message(expect_message(
                se1 <- readModBam(bamfiles = modbamfiles, regions = reg1,
                                  modbase = "a", nAlnsToSample = 0,
                                  sequence.context.width = 1, sequence.reference = ref,
                                  seqnamesToSampleFrom = "chr1", verbose = TRUE)
    ))))))))))
    se2 <- readModBam(bamfiles = unname(modbamfiles),
                      regions = reg2,
                      modbase = "a",
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      verbose = FALSE)
    se3 <- readModBam(bamfiles = modbamfiles,
                      regions = reg3,
                      modbase = "a",
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      verbose = FALSE)
    se4 <- readModBam(bamfiles = modbamfiles,
                      regions = reg4,
                      modbase = c("a", "m"),
                      nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                      verbose = FALSE)
    se5a <- readModBam(bamfiles = modbamfiles[1],
                       regions = reg5[1],
                       modbase = "a",
                       nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
                       verbose = FALSE)
    se5b <- readModBam(bamfiles = modbamfiles[1],
                       regions = reg5[1:2],
                       modbase = "a",
                       nAlnsToSample = 0, seqnamesToSampleFrom = "chr1",
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
    set.seed(55L)
    expect_message(
        expect_message(
            expect_message(
                expect_message(
                    expect_message(
                        expect_message(
                            expect_message(
                                se6a  <- readModBam(bamfiles = modbamfiles[1],
                                                    regions = NULL,
                                                    modbase = "a",
                                                    nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                                                    verbose = TRUE),
                                "extracting base modifications"),
                            "opening input file"),
                        "sampling alignments with probability 0.5"),
                    "reading alignments overlapping"),
                "removed 2006 unaligned"),
            "finding unique genomic"),
        "collapsed 16095 positions to 6852")
    set.seed(55L)
    se6b  <- readModBam(bamfiles = modbamfiles[1],
                        regions = NULL,
                        modbase = "a",
                        nAlnsToSample = 5, seqnamesToSampleFrom = "chr1",
                        verbose = FALSE)

    seL <- list(se1, se2, se3, se4, se5a, se5b, se6a, se6b)

    # ... structure
    expected_coldata_names <- c("sample", "modbase", "n_reads", "read_info")
    expected_read_info_names <- c("read_id", "qscore", "read_length",
                                  "aligned_length", "aligned_fraction")
    for (se in seL) {
        expect_s4_class(se, "RangedSummarizedExperiment")
        expect_s4_class(rowRanges(se), "GPos")
        expect_identical(colnames(colData(se)), expected_coldata_names)
        expect_s4_class(colData(se)$read_info, "SimpleList")
        lapply(colData(se)$read_info, function(df) {
            expect_s3_class(df, "data.frame")
            expect_named(df, expected_read_info_names)
        })
        expect_identical(assayNames(se), "mod_prob")
        expect_s4_class(assay(se, "mod_prob"), "DFrame")
        expect_s4_class(assay(se, "mod_prob")[[1]], "NaMatrix")
        expect_equal(vapply(assay(se, "mod_prob"), ncol, 0), se$n_reads,
                     ignore_attr = TRUE)
    }

    expect_identical(se1$sample, names(modbamfiles))
    expect_identical(se2$sample, c("s1", "s2"))
    expect_identical(se3$sample, names(modbamfiles))
    expect_identical(se4$sample, names(modbamfiles))
    expect_identical(se5a$sample, names(modbamfiles)[1])
    expect_identical(se5b$sample, names(modbamfiles)[1])
    expect_identical(se6a$sample, names(modbamfiles)[1])
    expect_identical(se6b$sample, names(modbamfiles)[1])

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
    expect_identical(se1$sample, names(modbamfiles))
    expect_identical(lapply(se1$read_info, "[[", "read_id"),
                     lapply(assay(se1, "mod_prob"), colnames))
    expect_equal(lapply(se1$read_info, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973,
                                 21.1338005065918, 20.3082008361816),
                     sample2 = c(12.9041996002197, 9.67461013793945,
                                 15.0149002075195, 15.1365995407104,
                                 17.7175006866455, 13.6647996902466)))
    expect_identical(lapply(se1$read_info, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 9246L, 12277L),
                         sample2 = c(13108L, 11834L, 9674L, 10047L, 8973L, 10057L)
                     ))
    expect_identical(lapply(se1$read_info, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 9227L, 12227L),
                         sample2 = c(9656L, 11234L, 9579L, 9967L, 8915L, 9898L)
                     ))
    expect_equal(unclass(table(as.character(SummarizedExperiment::rowData(se1)$sequence.context))),
                 c(A = 8108L, C = 128L, G = 393L, T = 62L), ignore_attr = TRUE)

    # ... content se2
    expect_identical(unname(se2$n_reads), c(3L, 2L))
    expect_identical(dim(se2), dim(se3))
    expect_identical(unname(as.matrix(assay(se2, "mod_prob"))),
                     unname(as.matrix(assay(se3, "mod_prob"))))
    expect_identical(sub("^s", "sample", colnames(se2)), colnames(se3))
    expect_identical(se2$sample, sub("sample", "s", se3$sample))
    for (nm in setdiff(expected_read_info_names, "read_id")) {
        expect_equal(lapply(se2$read_info, "[[", nm),
                     lapply(se3$read_info, "[[", nm),
                     ignore_attr = TRUE)
    }

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
    expect_identical(lapply(se3$read_info, "[[", "read_id"),
                     lapply(assay(se3, "mod_prob"), colnames))
    expect_equal(lapply(se3$read_info, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                     sample2 = c(9.67461013793945, 13.6647996902466)))
    expect_identical(lapply(se3$read_info, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 12277L),
                         sample2 = c(11834L, 10057L)
                     ))
    expect_identical(lapply(se3$read_info, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 12227L),
                         sample2 = c(11234L, 9898L)
                     ))

    # ... content se4
    expect_identical(unname(se4$n_reads), c(3L, 0L))
    expect_identical(dim(se4), c(4772L, 2L))
    expect_identical(dim(as.matrix(assay(se4, "mod_prob"))), c(4772L, 3L))
    expect_identical(unlist(lapply(se4$read_info, "[[", "read_id"), use.names = FALSE),
                     unlist(lapply(assay(se4, "mod_prob"), colnames), use.names = FALSE))
    expect_equal(lapply(se4$read_info, "[[", "qscore"),
                 list(
                     sample1 = c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                     sample2 = numeric(0)))
    expect_identical(lapply(se4$read_info, "[[", "read_length"),
                     list(
                         sample1 = c(20058L, 11305L, 12277L),
                         sample2 = integer(0)
                     ))
    expect_identical(lapply(se4$read_info, "[[", "aligned_length"),
                     list(
                         sample1 = c(14801L, 11214L, 12227L),
                         sample2 = integer(0)
                     ))

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
    expect_identical(dim(assay(se6a)$sample1), c(6852L, 6L))
})
