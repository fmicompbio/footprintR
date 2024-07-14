suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
})

## -------------------------------------------------------------------------- ##
## Checks, readModkitExtract
## -------------------------------------------------------------------------- ##
test_that("readModkitExtract works", {
    # example data
    fnames <- c(s1_5mC = system.file("extdata", "modkit_extract_rc_5mC_1.tsv.gz",
                                     package = "footprintR"),
                s2_5mC = system.file("extdata", "modkit_extract_rc_5mC_2.tsv.gz",
                                     package = "footprintR"),
                s1_6mA = system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                                     package = "footprintR"),
                s2_6mA = system.file("extdata", "modkit_extract_rc_6mA_2.tsv.gz",
                                     package = "footprintR")
    )

    # invalid arguments
    expect_error(readModkitExtract("error"), "not all `fnames` exist")
    expect_error(readModkitExtract(fnames, modbase = NULL),
                 "'modbase' must not be NULL")
    expect_error(readModkitExtract(fnames, modbase = 1),
                 "'modbase' must be of class 'character'")
    expect_error(readModkitExtract(fnames, modbase = c("m", "a")),
                 "'modbase' must have length 4")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   filter = "error"),
                 "All values in 'filter' must be one of")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   filter = c(0.1, 0.2)),
                 "`filter` must be a named vector")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   filter = c(c = 0.1)),
                 "a filter threshold needs to be supplied")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   nrows = -1),
                 "'nrows' must be within")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   nrows = "error"),
                 "'nrows' must be of class 'numeric'")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   seqinfo = c(100)),
                 "`seqinfo` must be ")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   seqinfo = c(chr2 = 1000)),
                 "'seqnames' contains sequence names with no entries")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   ncpu = "error"),
                 "'ncpu' must be of class 'numeric'")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   ncpu = c(1, 2)),
                 "'ncpu' must have length 1")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   verbose = "error"),
                 "'verbose' must be of class 'logical'")
    expect_error(readModkitExtract(fnames, modbase = c("m", "m", "a", "a"),
                                   verbose = c(TRUE, FALSE)),
                 "'verbose' must have length 1")

    # expected results
    # ... single file, no filtering
    rme <- readModkitExtract(fnames = fnames["s1_5mC"], modbase = "m",
                             filter = NULL, nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(6432, 10)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = NULL),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 2471.8587)
    expect_equal(sum(SummarizedExperiment::assay(rme) > 0), 18531) ## number of rows in the original file

    # ... single file, manual filtering
    rme <- readModkitExtract(fnames = fnames["s1_5mC"], modbase = "m",
                             filter = c(`m` = 0.6, `-` = 0.5),
                             nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(6415, 10)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.6, `-` = 0.5)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 2418.5403)
    expect_equal(sum(SummarizedExperiment::assay(rme) > 0), 18434) ## number of rows in the original file

    # ... single file, automatic filtering
    rme <- readModkitExtract(fnames = fnames["s1_5mC"], modbase = "m",
                             filter = "modkit",
                             nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(5893, 10)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 1824.64774)
    expect_equal(sum(SummarizedExperiment::assay(rme) > 0), 15325) ## number of rows in the original file

    # ... multiple files, no filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c("m", "m", "a"),
                             filter = NULL, nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(18655, 30)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = NULL, s2_5mC = NULL, s1_6mA = NULL),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 9054.297)
    expect_equal(sum(ceiling(SummarizedExperiment::assay(rme))), 71750) ## total number of rows in the original files

    # ... multiple files, manual filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c("m", "m", "a"),
                             filter = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                             nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(18615, 30)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                      s2_5mC = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                      s1_6mA = c(`m` = 0.6, `a` = 0.4, `-` = 0.3)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 8899.0412)
    expect_equal(sum(ceiling(SummarizedExperiment::assay(rme))), 71467) ## total number of rows in the original files

    # ... multiple files, automatic filtering
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC",
                                               "s1_6mA")],
                             modbase = c("m", "m", "a"),
                             filter = "modkit",
                             nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    expect_s4_class(rme, "RangedSummarizedExperiment")
    expect_equal(dim(rme), c(17459, 30)) ## number of unique positions
    expect_length(SummarizedExperiment::assays(rme), 1)
    expect_named(SummarizedExperiment::assays(rme), "mod_prob")
    expect_s4_class(SummarizedExperiment::assay(rme), "SVT_SparseMatrix")
    expect_false(is.null(colnames(rme)))
    expect_length(S4Vectors::metadata(rme), 2)
    expect_named(S4Vectors::metadata(rme), c("modkit_threshold",
                                             "filter_threshold"))
    expect_equal(S4Vectors::metadata(rme)$modkit_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(S4Vectors::metadata(rme)$filter_threshold,
                 list(s1_5mC = c(`m` = 0.7988281, `-` = 0.9082031),
                      s2_5mC = c(`m` = 0.7988281, `-` = 0.9003906),
                      s1_6mA = c(`a` = -Inf, `-` = 0.8964844)),
                 ignore_attr = TRUE)
    expect_equal(sum(SummarizedExperiment::assay(rme)), 6672.8205)
    expect_equal(sum(ceiling(SummarizedExperiment::assay(rme))), 61228) ## total number of rows in the original files
})
