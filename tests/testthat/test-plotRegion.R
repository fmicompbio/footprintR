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
    se <- readBedMethyl(fnames = c(fname1, fname2), modbase = "m",
                        sequence.context.width = 3,
                        sequence.reference = ref)
    se0 <- se
    assayNames(se0) <- c("assay1", "assay2")
    fname3 <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz", package = "footprintR")
    seR <- readModkitExtract(fnames = fname3, modbase = 'a')
    seR2 <- addReadsSummary(se = seR)

    # invalid arguments
    expect_error(plotRegion(se = "error"))
    expect_error(plotRegion(se = se0))
    expect_error(plotRegion(se = se, region = -1))
    expect_error(plotRegion(se = se, region = "error"))
    expect_error(plotRegion(se = seR, tracks.reads = "error"))
    expect_error(plotRegion(se = seR, tracks.reads = list(mod_prob = "error")))
    expect_error(plotRegion(se = se, tracks.summary = "error"))
    expect_error(plotRegion(se = se, tracks.summary = list(FracMod = "error")))
    expect_error(plotRegion(se = se, modbaseSpace = "error"))
    expect_error(plotRegion(se = se, sequence.context = 1))
    expect_error(plotRegion(se = seR, sequence.context = "C"))

    # expected results
    p1 <- plotRegion(se = se, region = "chr1:6948000-6952000")
    p2 <- plotRegion(se = se, tracks.summary = list(Nvalid = "Point"))
    p3 <- plotRegion(se = se, tracks.summary = list(FracMod = "Smooth"))
    p4 <- plotRegion(se = se, sequence.context = c("GCH"), modbaseSpace = TRUE)
    p5 <- plotRegion(se = se, sequence.context = c("GCA","GCC","GCT"), modbaseSpace = TRUE)
    p6 <- plotRegion(se = seR,
                     tracks.reads = list(mod_prob = c("Lollipop", "Heatmap")),
                     tracks.summary = NULL)
    p7 <- plotRegion(se = seR, modbaseSpace = TRUE,
                     tracks.reads = list(mod_prob = c("Heatmap")),
                     tracks.summary = NULL)
    expect_warning(
        p8 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                         modbaseSpace = TRUE,
                         tracks.summary = list(FracMod = "Smooth"),
                         tracks.reads = list(mod_prob = c("Lollipop", "Heatmap", "HeatmapFilled")))
    )
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
    expect_s3_class(p4, "ggplot")
    expect_s3_class(p5, "ggplot")
    expect_s3_class(p6, "ggplot")
    expect_s3_class(p7, "ggplot")
    expect_s3_class(p8, "ggplot")
    expect_identical(nrow(p1$data), 4006L)
    expect_identical(nrow(p2$data), 24040L)
    expect_identical(nrow(p3$data), 20000L)
    expect_identical(nrow(p4$data), 2459L)
    expect_identical(nrow(p5$data), 2459L)
    expect_identical(p4$data, p5$data)
    expect_identical(nrow(p6$data), 29104L)
    expect_identical(nrow(p7$data), 29104L)
    expect_identical(nrow(p8$data), 500L)

    # make sure the plotting works
    tmpplot <- tempfile(fileext = ".png")
    expect_identical(ggsave(filename = tmpplot, plot = p1, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p2, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p3, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p4, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p5, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p6, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p7, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p8, width = 6, height = 6), tmpplot)
    unlink(tmpplot)
})
