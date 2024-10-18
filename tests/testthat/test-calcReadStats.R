test_that("calcReadStats works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                          package = "footprintR")
    reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a")

    ## Expected errors
    expect_error(calcReadStats(se, assay.type = "error"), "must be one of: mod_prob")

    ## No coverage requirement
    rs <- calcReadStats(se, min.Nobs.ppos = 1)
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    expect_length(S4Vectors::metadata(rs), 5L)
    expect_named(S4Vectors::metadata(rs),
                 c("regions", "sequence.context", "min.Nobs.ppos",
                   "min.Nobs.pread", "Lags"))
    expect_equal(S4Vectors::metadata(rs)$min.Nobs.ppos, 1L)
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 12L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb") %in%
                        colnames(qc)))
    expect_equal(qc$MeanModProb,
                 colSums(assay(se)$s1, na.rm = TRUE) /
                     colSums(assay(se)$s1 > 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 colSums(assay(se)$s1 > 0.5, na.rm = TRUE) /
                     colSums(assay(se)$s1 > 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_type(S4Vectors::metadata(qc), "list")

    ## Default coverage requirement
    Nobs <- rowSums(SummarizedExperiment::assay(se)[["s1"]] > 0, na.rm = TRUE)
    thr <- max(floor(stats::quantile(Nobs, 0.75) -
                         0.5 * stats::IQR(Nobs)), 1L)
    idx <- which(Nobs >= thr)
    rs <- calcReadStats(se, verbose = TRUE, min.Nobs.ppos = thr)
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    expect_type(S4Vectors::metadata(rs), "list")
    expect_length(S4Vectors::metadata(rs), 5L)
    expect_named(S4Vectors::metadata(rs),
                 c("regions", "sequence.context", "min.Nobs.ppos",
                   "min.Nobs.pread", "Lags"))
    expect_equal(S4Vectors::metadata(rs)$min.Nobs.ppos, thr)
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 12L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb") %in%
                        colnames(qc)))
    expect_equal(qc$MeanModProb,
                 colSums(assay(se)$s1[idx, ], na.rm = TRUE) /
                     colSums(assay(se)$s1[idx, ] > 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 colSums(assay(se)$s1[idx, ] > 0.5, na.rm = TRUE) /
                     colSums(assay(se)$s1[idx, ] > 0, na.rm = TRUE),
                 ignore_attr = TRUE)

    ## Using `regions` and large LagRange
    rs1 <- calcReadStats(se, regions = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(6935000, 6935100)), LagRange = c(200, 256),
        min.Nobs.ppos = 5)
    rs2 <- calcReadStats(se, regions = "chr1:6935000-6935100",
                         LagRange = c(200, 256), min.Nobs.ppos = 5)
    expect_identical(rs1, rs2)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 12L))
    expect_equal(sum(rs1$s1$MeanModProb), 1.52775200714286)
    expect_true(all(vapply(rs1$s1$ACModProb, function(x) all(x == 0), TRUE)))
    expect_true(all(vapply(rs1$s1$PACModProb, function(x) all(x == 0), TRUE)))

    ## Using `sequence.context`, `min.Nobs.pread` and `stats`
    expect_error(calcReadStats(se, regions = "chr1:6935000-6935100",
                               sequence.context = c("TAA", "AAA")),
                 "No sequence context found")
    se1 <- addSeqContext(se, sequence.context.width = 3,
                         sequence.reference = reffile)
    expect_identical(colnames(rowData(se1)), "sequence.context")
    rs1 <- calcReadStats(se1, regions = "chr1:6935000-6936000",
                         sequence.context = c("TAA", "AAA"), min.Nobs.ppos = 5,
                         min.Nobs.pread = 1, stats = "MeanModProb")
    rs2 <- calcReadStats(se1, regions = "chr1:6935000-6936000",
                         sequence.context = "WAA", min.Nobs.ppos = 5,
                         min.Nobs.pread = 1, stats = "MeanModProb")
    expect_named(rs1, "s1")
    expect_named(rs2, "s1")
    # ignore metadata()$sequence.context (expected to differ, explicit vs. IUPAC code)
    meta_names <- setdiff(names(metadata(rs1)), "sequence.context")
    expect_identical(metadata(rs1)[meta_names], metadata(rs2)[meta_names])
    expect_identical(rs1$s1, rs2$s1)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 1L))
    expect_equal(sum(rs1$s1$MeanModProb), 0.65811757757861633067)
})

test_that("addReadStats works", {
    # example data
    exfiles <- system.file("extdata", c("modkit_extract_rc_6mA_1.tsv.gz",
                                        "modkit_extract_rc_6mA_2.tsv.gz"),
                           package = "footprintR")
    se <- readModkitExtract(exfiles, modbase = "a")
    se2 <- addReadStats(se, name = "qc2")
    se3 <- addReadStats(se, min.Nobs.pread = 2600, name = "qc2")

    # expected errors
    expect_error(addReadStats(se, name = -1), "must be of class 'character'")
    expect_error(addReadStats(se, name = c("a", "b")), "must have length 1")

    # expected results
    expect_s4_class(se2, "SummarizedExperiment")
    expect_equal(dim(se), dim(se2))
    expect_equal(assay(se), assay(se2))
    expect_null(se2[["QC"]])
    expect_s4_class(se2$qc2, "SimpleList")
    expect_length(se2$qc2, 2L)
    expect_named(se2$qc2, c("s1", "s2"))
    qc <- se2$qc2[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 12L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb") %in%
                        colnames(qc)))
    expect_identical(metadata(se2$qc2)$min.Nobs.pread, 0)
    expect_identical(metadata(se3$qc2)$min.Nobs.pread, 2600)
    na_rows <- lapply(endoapply(assay(se), function(x) colSums(is_nonna(x))),
                      function(y) which(y < 2600))
    expect_equal(na_rows, list(c(7,8,9,10), c(5,7,8,9,10)), ignore_attr = TRUE)
    expect_identical(se2$qc2$s1[-na_rows$s1,], se3$qc2$s1[-na_rows$s1,])
    expect_identical(se2$qc2$s2[-na_rows$s2,], se3$qc2$s2[-na_rows$s2,])
})
