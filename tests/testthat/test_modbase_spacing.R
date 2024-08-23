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
    expect_error(calcModbaseSpacing("error"), "must be of class")
    expect_error(calcModbaseSpacing(se, "error"), "'assay.type' must be a string or integer")
    pg1 <- calcModbaseSpacing(se)
    pg2 <- calcModbaseSpacing(se, rmdup = FALSE)
    pg3 <- calcModbaseSpacing(cbind(se, se))

    ## check invalid arguments
    expect_warning(estimateNRL(rep(0, 3000L)), "no non-zero distances")
    expect_warning(estimateNRL(pg1, usePeaks = 1:100), "less peaks detected")
    expect_error(calcAndCountDist(1:3, 3:1, numeric(3)), "must be sorted ascendingly")

    ## check if the plotting function runs
    p1 <- plotModbaseSpacing(x = pg1, hide = FALSE)
    p2 <- plotModbaseSpacing(x = pg1, detailedPlots = TRUE)
    expect_s3_class(p1, "ggplot")
    expect_identical(dim(p1$data), c(2720L, 3L))
    expect_s3_class(p2, "ggplot")
    expect_identical(dim(p2$data), c(5L, 2L))
    tf <- tempfile(fileext = ".pdf")
    pdf(file = tf)
    print(p1)
    print(p2)
    dev.off()
    unlink(tf)

    ## check expected results
    expect_type(pg1, "double")
    expect_length(pg1, 1000L)
    expect_true(all(pg2 >= pg1))
    expect_true(all(pg1 < pg3))

    nrl <- estimateNRL(pg1, usePeaks = 1:5)
    expect_type(nrl, "list")
    expect_equal(nrl$nrl, 183.9)
    expect_length(nrl$nrl.CI95, 2L)

    expect_equal(calcAndCountDist(c(1,2,4),c(1,3,5),numeric(4)), c(2,1,1,1))
    expect_equal(calcAndCountDist(c(1,3,5),c(1,2,4),numeric(4)), c(2,0,1,0))
})
