test_that("calcReadStats works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                          package = "footprintR")
    reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a")

    ## Expected errors
    se2 <- se
    SummarizedExperiment::assayNames(se2) <- "error"
    expect_error(calcReadStats(se2), "needs to have a 'mod_prob' assay")
    rm(se2)

    ## No coverage requirement
    rs <- calcReadStats(se, min.Nobs.ppos = 1)
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 13L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb", "sample") %in%
                        colnames(qc)))
    expect_equal(qc$MeanModProb,
                 Matrix::colSums(assay(se)$s1) / Matrix::colSums(assay(se)$s1 > 0),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 Matrix::colSums(assay(se)$s1 > 0.5) / Matrix::colSums(assay(se)$s1 > 0),
                 ignore_attr = TRUE)
    expect_equal(unique(qc$sample), "s1")
    expect_type(S4Vectors::metadata(qc), "list")
    expect_length(S4Vectors::metadata(qc), 3L)
    expect_named(S4Vectors::metadata(qc), c("min.Nobs.ppos", "Lags", "stats"))
    expect_equal(S4Vectors::metadata(qc)$min.Nobs.ppos, 1L)

    ## Default coverage requirement
    Nobs <- rowSums(SummarizedExperiment::assay(se)[["s1"]] > 0)
    thr <- max(floor(stats::quantile(Nobs, 0.75) -
                         0.5 * stats::IQR(Nobs)), 1L)
    idx <- which(Nobs >= thr)
    expect_message(
        rs <- calcReadStats(se, verbose = TRUE)
    )
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 13L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb", "sample") %in%
                        colnames(qc)))
    expect_equal(qc$MeanModProb,
                 Matrix::colSums(assay(se)$s1[idx, ]) / Matrix::colSums(assay(se)$s1[idx, ] > 0),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 Matrix::colSums(assay(se)$s1[idx, ] > 0.5) / Matrix::colSums(assay(se)$s1[idx, ] > 0),
                 ignore_attr = TRUE)
    expect_equal(unique(qc$sample), "s1")
    expect_type(S4Vectors::metadata(qc), "list")
    expect_length(S4Vectors::metadata(qc), 3L)
    expect_named(S4Vectors::metadata(qc), c("min.Nobs.ppos", "Lags", "stats"))
    expect_equal(S4Vectors::metadata(qc)$min.Nobs.ppos, thr)

    ## Using `regions` and large LagRange
    rs1 <- calcReadStats(se, regions = GenomicRanges::GRanges("chr1", IRanges::IRanges(6935000, 6935100)), LagRange = c(200, 256))
    rs2 <- calcReadStats(se, regions = "chr1:6935000-6935100", LagRange = c(200, 256))
    expect_identical(rs1, rs2)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 13L))
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
                         sequence.context = c("TAA", "AAA"),
                         min.Nobs.pread = 1, stats = "MeanModProb")
    rs2 <- calcReadStats(se1, regions = "chr1:6935000-6936000",
                         sequence.context = "WAA",
                         min.Nobs.pread = 1, stats = "MeanModProb")
    expect_identical(rs1, rs2)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 2L))
    expect_equal(sum(rs1$s1$MeanModProb), 0.95680377182716)

    ## Test that addReadStats also works
    exfiles <- system.file("extdata", c("modkit_extract_rc_6mA_1.tsv.gz",
                                        "modkit_extract_rc_6mA_2.tsv.gz"),
                           package = "footprintR")
    se <- readModkitExtract(exfiles, modbase = "a")
    se2 <- addReadStats(se)

    expect_s4_class(se2, "SummarizedExperiment")
    expect_equal(dim(se), dim(se2))
    expect_equal(assay(se), assay(se2))
    expect_s4_class(SummarizedExperiment::colData(se2)$QC, "SimpleList")
    expect_length(SummarizedExperiment::colData(se2)$QC, 2L)
    expect_named(SummarizedExperiment::colData(se2)$QC, c("s1", "s2"))
    qc <- SummarizedExperiment::colData(se2)$QC[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 13L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm", "MeanConfMod",
                      "FracLowConf", "IQRModProb", "sdModProb", "SEntrModProb", "Lag1DModProb",
                      "ACModProb", "PACModProb", "sample") %in%
                        colnames(qc)))

})
