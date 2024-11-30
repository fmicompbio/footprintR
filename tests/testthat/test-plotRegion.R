suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(ggplot2)
    library(GenomicRanges)
    library(patchwork)
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
    
    tmpgrl <- GenomicRanges::GRangesList(x = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(c(6929104, 6929106), 
                                 c(6941530, 6941530)), c("+", "+")
    ))
    names(tmpgrl) <- NA_character_
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "GenomicRegion", trackData = tmpgrl))),
        "NA values are not allowed")
    
    tmpgrl2 <- GenomicRanges::GRangesList(
        x = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(c(6929104, 6929106), 
                                     c(6941530, 6941530)), c("+", "+")),
        y = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(c(6929104, 6929106), 
                                     c(6941530, 6941530)), c("+", "+")
        ))
    names(tmpgrl2) <- c("x", NA_character_)
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "GenomicRegion", trackData = tmpgrl2))),
        "NA values are not allowed")

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
    p11 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                      modbaseSpace = TRUE, referenceCoordinate = 6935400,
                      tracks = list(list(trackType = "Lollipop", 
                                         trackData = "mod_prob")))
    p12 <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                      modbaseSpace = TRUE, 
                      tracks = list(list(trackType = "Lollipop", 
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
                                         ))))
    expect_error(
        plotRegion(se = seR2, region = "chr1:693-1093",
                   modbaseSpace = FALSE,
                   tracks = list(list(trackType = "Lollipop",
                                      trackData = "mod_prob",
                                      highlightRegions = GenomicRanges::GRanges(
                                          "chr1", IRanges::IRanges(
                                              6935410, 6935430
                                          )
                                      )),
                                 list(trackType = "Smooth",
                                      trackData = "FracMod"))),
        "No positions retained for plotting")
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450", 
                   tracks = list(list(trackType = "Smooth",
                                      trackData = "Nvalid",
                                      colors = c(X = "red")))),
        "Missing color specification"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450", 
                   tracks = list(list(trackType = "Smooth",
                                      trackData = "Nvalid",
                                      colorBy = "missing"))),
        "Some requested columns are not present"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450", 
                   tracks = list(list(trackType = "Heatmap",
                                      trackData = "mod_prob",
                                      facetBy = "missing"))),
        "Some requested columns are not present"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450", 
                   tracks = list(list(trackType = "Smooth",
                                      trackData = "Nvalid",
                                      smoothMethod = "rollingMean",
                                      windowSize = 4))),
        "windowSize must be an odd integer"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450", 
                   tracks = list(list(trackType = "Smooth",
                                      trackData = "Nvalid",
                                      smoothMethod = "missing"))),
        "All values in 'smoothMethod' must be one of"
    )
    
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
    expect_s3_class(p11, "ggplot")
    expect_s3_class(p12, "ggplot")
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
    expect_identical(ggsave(filename = tmpplot, plot = p11, width = 6, height = 6), tmpplot)
    expect_identical(ggsave(filename = tmpplot, plot = p12, width = 6, height = 6), tmpplot)
    unlink(tmpplot)
})

## The examples below are provided with the intention that they can be run 
## manually from time to time to also inspect the output (that is hard to 
## capture with automated unit tests). 
test_that("plotRegion works - manual inspection", {
    extractfiles <- system.file("extdata",
                                c("modkit_extract_rc_6mA_1.tsv.gz",
                                  "modkit_extract_rc_6mA_2.tsv.gz"),
                                package = "footprintR")
    seB <- readModkitExtract(extractfiles, modbase = "a", filter = "modkit",
                             BPPARAM = BiocParallel::SerialParam())
    seB <- flattenReadLevelAssay(seB)
    
    ## Annotation GRangesList
    grl <- GRangesList(
        CGI1 = GRanges(seqnames = "chr1", 
                       ranges = IRanges(start = 6935820, end = 6935850), 
                       strand = "*"),
        g1 = GRanges(seqnames = "chr1", 
                     ranges = IRanges(start = c(6934800, 6935870), end = c(6935820, 6935900)), 
                     strand = "+"),
        g1 = GRanges(seqnames = "chr1",
                     ranges = IRanges(start = c(6934800, 6935840), end = c(6935810, 6935950)),
                     strand = "-"),
        out = GRanges(seqnames = "chr2", 
                      ranges = IRanges(start = 6935820, end = 6935850), 
                      strand = "*")
    )
    
    ## Regions to highlight
    grh <- GRanges(seqnames = c("chr1", "chr1", "chr2"), 
                   ranges = IRanges(start = c(6935830, 6935700, 6935870),
                                    end = c(6935850, 6935820, 6935890)))
    
    ## Absolute base space, no reference coordinate
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        labelAccuracy = 1e-6,
        tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
                           size = 2, stroke = 0.25, legendTitle = "6mA", 
                           highlightRegions = grh),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = TRUE, labelSize = 3, 
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE))) + 
        plot_layout(heights = c(3, 1, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## ... don't color by strand, move labels
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA",
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = FALSE, labelSize = 3, 
                           labelPosition = "above", legendTitle = NULL),
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE, 
                           highlightRegions = grh))) + 
        plot_layout(heights = c(3, 1, 2))
    expect_s3_class(p, "ggplot")
    
    ## ... interpolate
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA",
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = TRUE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = FALSE, labelSize = 3, 
                           labelPosition = "above", legendTitle = NULL),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", 
                           highlightRegions = grh))) + 
        plot_layout(heights = c(3, 1, 2))
    expect_s3_class(p, "ggplot")
    
    ## referenceCoordinate = left border of plot
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = TRUE, labelSize = 2, 
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", 
                           highlightRegions = grh))) + 
            plot_layout(heights = c(3, 1, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## referenceCoordinate outside plot (left) - breaks in scale_cut!
    ## fixed in https://github.com/r-lib/scales/commit/6f2f979a81678c7cd5597b1d18cac78e9cf473c6, 
    ## but not yet on CRAN
    # expect_warning(p <- plotRegion(
    #     seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
    #     referenceCoordinate = 6934800,
    #     tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
    #                        legendTitle = "6mA", highlightRegions = grh,
    #                        orderReads = FALSE, trackTitle = "Heatmap",
    #                        facetBy = NULL, interpolate = FALSE,
    #                        linewidthTiles = 0.25),
    #                   list(trackData = grl, trackType = "GenomicRegion",
    #                        colorByStrand = TRUE, labelSize = 2,
    #                        labelPosition = "inside", legendTitle = NULL),
    #                   list(trackData = "mod_prob", trackType = "Lollipop",
    #                        legendTitle = "6mA", highlightRegions = grh,
    #                        orderReads = TRUE, facetBy = "sample",
    #                        size = 2, stroke = 0.5),
    #                   list(trackData = "Nvalid", trackType = "PointSmooth",
    #                        showLegend = FALSE, spar = 0.5,
    #                        trackTitle = "Smooth",
    #                        highlightRegions = grh))) +
    #         plot_layout(heights = c(3, 1, 3, 2)),
    #     "the standard deviation is zero")
    # expect_s3_class(p, "ggplot")

    ## referenceCoordinate outside plot (right)
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        referenceCoordinate = 6936800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = "modbase", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = TRUE, labelSize = 2, 
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", 
                           highlightRegions = grh))) + 
            plot_layout(heights = c(3, 1, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## referenceCoordinate in the middle of plot region
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE, 
        referenceCoordinate = 6935850,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion", 
                           colorByStrand = TRUE, labelSize = 2, 
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", groupBy = "modbase",
                           colorBy = "modbase",
                           highlightRegions = grh))) + 
            plot_layout(heights = c(3, 1, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## modbaseSpace = TRUE, change colors
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE, 
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", 
                           colors = c(s1 = "forestgreen", s2 = "firebrick1"),
                           highlightRegions = grh))) + 
            plot_layout(heights = c(3, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## ... with only smooth
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE, 
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", colorBy = "modbase",
                           highlightRegions = grh,
                           arglistSmooth = list(linewidth = 2)))) + 
            plot_layout(heights = c(3, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## ... change y-axis range
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE, 
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = FALSE, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = TRUE, facetBy = "sample", 
                           size = 2, stroke = 0.5), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth", colorBy = "modbase",
                           highlightRegions = grh, yAxisRange = c(3, 9)))) + 
            plot_layout(heights = c(3, 3, 2)),
        "the standard deviation is zero")
    expect_s3_class(p, "ggplot")
    
    ## ... compare smoothing methods
    p <- plotRegion(
        seB, region = "chr1:6935700-6936000", modbaseSpace = FALSE, 
        tracks = list(list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.1, 
                           trackTitle = "Smooth spline, spar = 0.1", 
                           smoothMethod = "smoothSpline"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth spline, spar = 0.5",
                           smoothMethod = "smoothSpline"), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, windowSize = 3, 
                           trackTitle = "Rolling mean, windowSize = 3",
                           smoothMethod = "rollingMean"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, windowSize = 15, 
                           trackTitle = "Rolling mean, windowSize = 15",
                           smoothMethod = "rollingMean")
        ))
    expect_s3_class(p, "ggplot")
    
    ## ... compare smoothing methods, in modbaseSpace
    p <- plotRegion(
        seB, region = "chr1:6935700-6936000", modbaseSpace = TRUE, 
        tracks = list(list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.1, 
                           trackTitle = "Smooth spline, spar = 0.1", 
                           smoothMethod = "smoothSpline"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5, 
                           trackTitle = "Smooth spline, spar = 0.5",
                           smoothMethod = "smoothSpline"), 
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, windowSize = 3, 
                           trackTitle = "Rolling mean, windowSize = 3",
                           smoothMethod = "rollingMean"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, windowSize = 15, 
                           trackTitle = "Rolling mean, windowSize = 15",
                           smoothMethod = "rollingMean")
        ))
    expect_s3_class(p, "ggplot")
})
