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
                        sequenceContextWidth = 3,
                        sequenceReference = ref, 
                        BPPARAM = BiocParallel::SerialParam())
    se0 <- se
    assayNames(se0) <- c("assay1", "assay2")
    fname3 <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz", package = "footprintR")
    seR <- readModkitExtract(fnames = fname3, modbase = 'a', 
                             BPPARAM = BiocParallel::SerialParam())
    seR2 <- flattenReadLevelAssay(se = seR)

    # invalid arguments
    expect_error(plotRegion(se = "error"))
    expect_error(plotRegion(se = se0))
    expect_error(plotRegion(se = se, region = -1))
    expect_error(plotRegion(se = se, region = "error"))
    expect_error(plotRegion(se = seR, tracks = "error"))
    expect_error(plotRegion(se = seR, tracks = list(list(trackData = "mod_prob"))))
    expect_error(plotRegion(se = seR, tracks = list(list("mod_prob", "Lollipop"))))
    expect_error(plotRegion(se = seR, tracks = list(list(trackData = "mod_prob",
                                                         trackType = "error"))))
    expect_error(plotRegion(se = seR, tracks = list(list(trackData = "error",
                                                         trackType = "Lollipop"))))
    expect_error(plotRegion(se = se, tracks = list(list(trackData = "Nvalid",
                                                        trackType = "Lollipop"))))
    expect_error(plotRegion(se = seR, tracks = list(list(trackData = "mod_prob",
                                                         trackType = "Smooth"))))
    expect_error(plotRegion(se = se, modbaseSpace = "error"))
    expect_error(plotRegion(se = se, sequenceContext = 1))
    expect_error(plotRegion(se = seR2, sequenceContext = "C"),
                 "No sequence context found")
    expect_error(plotRegion(se = seR2, referenceCoordinate = "1"),
                 "'referenceCoordinate' must be of class 'numeric'")
    expect_error(plotRegion(se = seR2, referenceCoordinate = c(1, 2)),
                 "'referenceCoordinate' must have length 1")
    expect_error(plotRegion(se = seR2, tracks = list(list(trackType = "GenomicRegion",
                                                          trackData = 1))),
                 "must be a named GRangesList object")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "GenomicRegion",
        trackData = GenomicRanges::GRangesList(GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(c(1, 4), c(3, 7)), c("+", "-")
        ))))),
        "must be a named GRangesList object")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "GenomicRegion",
        trackData = GenomicRanges::GRangesList(x = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(c(1, 4), c(3, 7)), c("+", "-")
        ))))),
        "There are entries in")

    # expected results
    p1 <- plotRegion(se = se, region = "chr1:6948000-6952000")
    p2 <- plotRegion(se = se, tracks = list(list(trackData = "Nvalid",
                                                 trackType = "Point")))
    p3 <- plotRegion(se = se, tracks = list(list(trackData = "FracMod",
                                                 trackType = "Smooth")))
    p4 <- plotRegion(se = se, sequenceContext = c("GCH"), modbaseSpace = TRUE)
    p5 <- plotRegion(se = se, sequenceContext = c("GCA","GCC","GCT"), modbaseSpace = TRUE)
    p6 <- plotRegion(se = seR,
                     tracks = list(list(trackData = "mod_prob",
                                        trackType = "Lollipop"),
                                   list(trackData = "mod_prob",
                                        trackType = "Heatmap",
                                        highlightRegions = GenomicRanges::GRanges(
                                            "chr1", IRanges::IRanges(
                                                6926200, 6935400
                                            )
                                        ))))
    p7 <- plotRegion(se = seR, modbaseSpace = TRUE,
                     tracks = list(list(trackData = "mod_prob",
                                        trackType = "Heatmap")))
    expect_warning(
        p8 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                         modbaseSpace = TRUE,
                         tracks = list(list(trackData = "FracMod",
                                            trackType = "Smooth"),
                                       list(trackData = "mod_prob",
                                            trackType = "Lollipop"),
                                       list(trackData = "mod_prob",
                                            trackType = "Heatmap"),
                                       list(trackData = "mod_prob",
                                            trackType = "Heatmap",
                                            interpolate = TRUE)))
    )
    expect_warning(
        p9 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                         modbaseSpace = TRUE, 
                         tracks = list(list(trackType = "GenomicRegion",
                                            trackData = GenomicRanges::GRangesList(
                                                a = GenomicRanges::GRanges(
                                                    "chr1", IRanges::IRanges(
                                                        6935420, 6935440
                                                    ), "+"
                                                )
                                            )))),
        "is not allowed if"
    )
    p10 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                      tracks = list(list(trackType = "Lollipop",
                                         trackData = "mod_prob",
                                         highlightRegions = GenomicRanges::GRanges(
                                             "chr1", IRanges::IRanges(
                                                 6935410, 6935430
                                             )
                                         )),
                                    list(trackType = "Heatmap",
                                         trackData = "mod_prob",
                                         highlightRegions = GenomicRanges::GRanges(
                                             "chr1", IRanges::IRanges(
                                                 6935410, 6935430
                                             )
                                         )),
                                    list(trackType = "Smooth",
                                         trackData = "FracMod",
                                         highlightRegions = GenomicRanges::GRanges(
                                             "chr1", IRanges::IRanges(
                                                 6935410, 6935430
                                             )
                                         )),
                                    list(trackType = "GenomicRegion",
                                         trackData = GenomicRanges::GRangesList(
                                             a = GenomicRanges::GRanges(
                                                 "chr1", IRanges::IRanges(
                                                     6935420, 6935440
                                                 ), "+"
                                             )
                                         ),
                                         colorByStrand = FALSE)), 
                      referenceCoordinate = 6935400)
    
    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
    expect_s3_class(p3, "ggplot")
    expect_s3_class(p4, "ggplot")
    expect_s3_class(p5, "ggplot")
    expect_s3_class(p6, "ggplot")
    expect_s3_class(p7, "ggplot")
    expect_s3_class(p8, "ggplot")
    expect_s3_class(p9, "ggplot")
    expect_s3_class(p10, "ggplot")
    expect_identical(nrow(p1$data), 4006L)
    expect_identical(nrow(p2$data), 24040L)
    expect_identical(nrow(p3$data), 20000L)
    expect_identical(nrow(p4$data), 2459L)
    expect_identical(nrow(p5$data), 2459L)
    expect_identical(p4$data, p5$data)
    expect_identical(nrow(p6$data), 29104L)
    expect_identical(nrow(p7$data), 29104L)
    expect_identical(nrow(p8$data), 500L)
    expect_length(p9$data, 0L)

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
    expect_identical(ggsave(filename = tmpplot, plot = p9, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p10, width = 6, height = 6), tmpplot)
    unlink(tmpplot)
})
