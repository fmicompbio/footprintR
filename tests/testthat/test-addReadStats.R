test_that("calcReadStats works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                          package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a")

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
    rs <- calcReadStats(se)
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
