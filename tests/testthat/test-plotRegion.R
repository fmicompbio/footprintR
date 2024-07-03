suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(ggplot2)
})

## -------------------------------------------------------------------------- ##
## Checks, plotRegion
## -------------------------------------------------------------------------- ##
test_that("plotRegion works", {
    # example data
    fname1 <- system.file("extdata", "modkit_pileup_1.bed.gz", package = "footprintR")
    fname2 <- system.file("extdata", "modkit_pileup_2.bed.gz", package = "footprintR")
    ref <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readBedMethyl(fnames = c(fname1, fname2), sequence.context.width = 3,
                        sequence.reference = ref)
    se0 <- se
    assayNames(se0) <- c("assay1", "assay2")

    # invalid arguments
    expect_error(plotRegion(se = "error"))
    expect_error(plotRegion(se = se0))
    expect_error(plotRegion(se = se, region = -1))
    expect_error(plotRegion(se = se, region = "error"))
    expect_error(plotRegion(se = se, min.coverage = -1))
    expect_error(plotRegion(se = se, min.coverage = "error"))
    expect_error(plotRegion(se = se, k.smooth = -1))
    expect_error(plotRegion(se = se, k.smooth = "error"))
    expect_error(plotRegion(se = se, sequence.context = 1))

    # expected results
    p1 <- plotRegion(se = se, region = "chr1:6948000-6952000")
    p2 <- plotRegion(se = se, min.coverage = 10)
    p3 <- plotRegion(se = se, k.smooth = 7)
    p4 <- plotRegion(se = se, sequence.context = c("GCH"))
    p5 <- plotRegion(se = se, sequence.context = c("GCA","GCC","GCT"))
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
    expect_s3_class(p4, "ggplot")
    expect_s3_class(p5, "ggplot")
    expect_identical(nrow(p1$data), 4792L)
    expect_identical(nrow(p2$data), 14559L)
    expect_identical(nrow(p3$data), 24040L)
    expect_identical(nrow(p4$data), 2512L)
    expect_identical(nrow(p5$data), 2512L)
    expect_identical(p4$data, p5$data)

    # make sure the plotting works
    tmpplot <- tempfile(fileext = ".png")
    expect_warning(ggsave(filename = tmpplot, plot = p1, width = 6, height = 6))
    ggsave(filename = tmpplot, plot = p2, width = 6, height = 6)
    expect_warning(expect_warning(ggsave(filename = tmpplot, plot = p3, width = 6, height = 6)))
    expect_warning(ggsave(filename = tmpplot, plot = p4, width = 6, height = 6))
    expect_warning(ggsave(filename = tmpplot, plot = p5, width = 6, height = 6))
    unlink(tmpplot)
})
