test_that("validity checks work", {
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
    rme <- readModkitExtract(fnames = fnames[c("s1_5mC", "s2_5mC", "s1_6mA")],
                             modbase = c(s1_6mA = "a", s1_5mC = "m", s2_5mC = "m"),
                             filter = c(`m` = 0.6, `a` = 0.4, `-` = 0.3),
                             nrows = Inf, seqinfo = NULL,
                             ncpu = 1L, verbose = FALSE)
    rme <- addReadStats(rme)
    rme_withreads <- addReadsSummary(rme, keep.reads = TRUE)
    rme_withoutreads <- addReadsSummary(rme, keep.reads = FALSE)

    ## Test .getReadLevelAssayNames
    expect_equal(.getReadLevelAssayNames(rme_withreads), "mod_prob")
    expect_equal(.getReadLevelAssayNames(rme_withoutreads), character(0))

    ## Test .checkSEValidity
    expect_no_error(.checkSEValidity(rme_withreads))
    expect_no_error(.checkSEValidity(rme_withoutreads))

    expect_message(
        expect_message(
            expect_message(
                expect_message(
                    expect_message(
                        .checkSEValidity(rme_withreads, verbose = TRUE),
                        "Checking assay names"),
                    "Checking row names"),
                "Checking consistency of sample names"),
            "Read-level assay found"),
        "QC information found, checking consistency")
    expect_message(
        expect_message(
            expect_message(
                .checkSEValidity(rme_withoutreads, verbose = TRUE),
                "Checking assay names"),
            "Checking row names"),
        "Checking consistency of sample names")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    expect_message(
        expect_message(
            expect_message(
                expect_message(
                    expect_message(
                        expect_message(
                            .checkSEValidity(rme1, verbose = TRUE),
                            "Checking assay names"),
                        "Checking row names"),
                    "Checking consistency of sample names"),
                "Read-level assay found"),
            "Comparing mod_prob and test"),
        "QC information found, checking consistency")

    rme1 <- rme_withreads
    assayNames(rme1) <- c("", "", "", "")
    expect_error(.checkSEValidity(rme1),
                 '!is.null(assayNames(se)) && all(assayNames(se) != "") && !any(duplicated(assayNames(se))) is not TRUE', fixed = TRUE)

    rme1 <- rme_withreads
    expect_equal(length(assays(rme1)), 4)
    assays(rme1) <- list(assays(rme1)[[1]], assays(rme1)[[2]], assays(rme1)[[3]],
                         assays(rme1)[[4]])
    expect_null(assayNames(rme1))
    expect_error(.checkSEValidity(rme1),
                 '!is.null(assayNames(se)) && all(assayNames(se) != "") && !any(duplicated(assayNames(se))) is not TRUE', fixed = TRUE)

    rme1 <- rme_withreads
    rme1$QC <- rme1$QC[c(3, 1, 2)]
    expect_error(.checkSEValidity(rme1),
                 "names(se$QC) == se$sample are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    SummarizedExperiment::colData(rme1) <- SummarizedExperiment::colData(rme1)[c(3, 1, 2), ]
    expect_error(.checkSEValidity(rme1),
                 "colnames(SummarizedExperiment::assay(se, an, withDimnames = FALSE))", fixed = TRUE)

    rme1 <- rme_withreads
    colnames(rme1) <- colnames(rme1)[c(3, 1, 2)]
    expect_error(.checkSEValidity(rme1),
                 "rownames(SummarizedExperiment::colData(se)) == se$sample are not all TRUE", fixed = TRUE)

    rme1 <- rme_withreads
    colnames(rme1) <- colnames(rme1)[c(3, 1, 2)]
    rme1$sample <- colnames(rme1)
    expect_error(.checkSEValidity(rme1),
                 "colnames(SummarizedExperiment::assay(se, an, withDimnames = FALSE))", fixed = TRUE)

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "mod_prob", withDimnames = FALSE) <-
        SummarizedExperiment::assay(rme1, "mod_prob")[, c(3, 1, 2)]
    expect_error(.checkSEValidity(rme1),
                 "colnames(SummarizedExperiment::assay(se, an, withDimnames = FALSE))", fixed = TRUE)

    rme1 <- rme_withreads
    rme1$QC[[2]] <- rme1$QC[[2]][1:5, ]
    expect_error(.checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and sample QC data, sample s2_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, 1:5]
    expect_error(.checkSEValidity(rme1),
                 "Mismatching reads for assay mod_prob and sample QC data, sample s1_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, 1:5]
    expect_error(.checkSEValidity(rme1),
                 "Mismatching reads for assays mod_prob and test, sample s1_5mC")

    rme1 <- rme_withreads
    SummarizedExperiment::assay(rme1, "test") <- SummarizedExperiment::assay(rme1, "mod_prob")
    N <- ncol(SummarizedExperiment::assay(rme1, "mod_prob")[[1]])
    set.seed(123L)
    SummarizedExperiment::assay(rme1, "mod_prob")[[1]] <-
        SummarizedExperiment::assay(rme1, "mod_prob")[[1]][, sample(seq_len(N), N)]
    expect_error(.checkSEValidity(rme1),
                 "Mismatching reads for assays mod_prob and test, sample s1_5mC")
})
