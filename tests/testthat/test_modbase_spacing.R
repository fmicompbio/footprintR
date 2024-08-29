suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
})

## -------------------------------------------------------------------------- ##
## Checks, calcModbaseSpacing, estimateNRL and calcAndCountDist
## -------------------------------------------------------------------------- ##
test_that("calcModbaseSpacing(), estimateNRL() and calcAndCountDist() work properly", {
    ## create modified-base distances using 6mA data (20 reads) from chr1:6940000-6955000
    bamf <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                        package = "footprintR")
    se <- readModBam(bamf, "chr1:6940000-6955000", "a")
    expect_error(calcModbaseSpacing("error"),
                 "must be of class 'RangedSummarizedExperiment'")
    expect_error(calcModbaseSpacing(se, "error"),
                 "'assay.type' must be a string or integer")
    expect_error(calcModbaseSpacing(se, pool_reads = "error"),
                 "must be of class 'logical'")
    pg1 <- calcModbaseSpacing(se)
    pg1comb <- Reduce("+", pg1)
    pg2 <- calcModbaseSpacing(se, pool_reads = FALSE)
    pg2comb <- Reduce("+", endoapply(pg2, rowSums))
    pg3 <- calcModbaseSpacing(cbind(se, se))
    pg3comb <- Reduce("+", pg3)

    ## check invalid arguments
    expect_warning(estimateNRL(rep(0, 3000L)), "no non-zero distances")
    expect_warning(estimateNRL(pg1comb, usePeaks = 1:100), "less peaks detected")
    expect_error(calcAndCountDist(1:3, 3:1, numeric(3)), "must be sorted ascendingly")

    ## check if the plotting function runs
    p1 <- plotModbaseSpacing(x = pg1comb, hide = FALSE, usePeaks = 2:4)
    p2 <- plotModbaseSpacing(x = pg1comb, detailedPlots = TRUE, usePeaks = 2:4)
    expect_s3_class(p1, "ggplot")
    expect_identical(dim(p1$data), c(2720L, 3L))
    expect_s3_class(p2, "ggplot")
    expect_identical(dim(p2$data), c(3L, 2L))
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)
    print(p1)
    print(p2)
    dev.off()
    unlink(tf)

    ## check expected results
    expect_type(pg1, "list")
    expect_length(pg1, ncol(se))
    expect_true(all(unlist(lapply(pg1, class)) == "numeric"))
    expect_type(pg1comb, "double")
    expect_length(pg1comb, 1000L)
    expect_true(all(unlist(lapply(pg2, function(x) inherits(x, "matrix")))))
    expect_identical(pg1, lapply(pg2, rowSums))
    expect_true(all(pg2comb >= pg1comb))
    expect_true(all(pg1comb < pg3comb))

    nrl <- estimateNRL(pg1comb, usePeaks = 1:5)
    expect_type(nrl, "list")
    expect_equal(nrl$nrl, 183.9)
    expect_length(nrl$nrl.CI95, 2L)

    expect_equal(calcAndCountDist(c(1,2,4),c(1,3,5),numeric(4)), c(2,1,1,1))
    expect_equal(calcAndCountDist(c(1,3,5),c(1,2,4),numeric(4)), c(2,0,1,0))
})
