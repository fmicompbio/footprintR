suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(ggplot2)
    library(GenomicRanges)
    library(patchwork)
})

## Helper functions
test_that(".calcDist works", {
    set.seed(123L)
    X <- matrix(rnorm(15), ncol = 3)
    expect_equal(as.matrix(.calcDist(X, clustDist = "euclidean"))[1, 3],
                 sqrt(sum((X[, 1] - X[, 3]) ^ 2)))
    expect_equal(as.matrix(.calcDist(X, clustDist = "cosine"))[1, 3],
                 1 - sum(X[, 1] * X[, 3]) / sqrt(sum(X[, 1] ^ 2) * sum(X[, 3] ^ 2)))
    expect_equal(as.matrix(.calcDist(X, clustDist = "pearson"))[1, 3],
                 sqrt(2 - 2 * cor(X[, 1], X[, 3], method = "pearson")))

    X[1, 3] <- NA
    expect_equal(as.matrix(.calcDist(X, clustDist = "euclidean"))[1, 3],
                 sqrt(sum((X[2:5, 1] - X[2:5, 3]) ^ 2 * 5 / 4)))
    expect_equal(as.matrix(.calcDist(X, clustDist = "pearson"))[1, 3],
                 sqrt(2 - 2 * cor(X[, 1], X[, 3], method = "pearson",
                                  use = "pairwise.complete")))
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
    seR2 <- addFootprints(seR2, wgt = c(0.5, -0.5), verbose = FALSE)

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
                 ".referenceCoordinate. must be of class .numeric.")
    expect_error(plotRegion(se = seR2, referenceCoordinate = c(1, 2)),
                 ".referenceCoordinate. must have length 1")
    expect_error(plotRegion(se = seR2, tracks = list(list(trackType = "GenomicRegion",
                                                          trackData = 1))),
                 "must be a named .GRangesList. object")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "GenomicRegion",
        trackData = GenomicRanges::GRangesList(GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(c(1, 4), c(3, 7)), c("+", "-")
        ))))),
        "must be a named .GRangesList. object")
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
        ".NA. values are not allowed")

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
        ".NA. values are not allowed")

    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "Lollipop", trackData = "mod_prob",
        footprintColumns = "missing"
    ))), "All values in .footprintColumns. must be one of")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "Lollipop", trackData = "mod_prob",
        footprintColumns = "footprints",
        arglistFootprints = 1
    ))), ".arglistFootprints. must be of class .list.")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "Heatmap", trackData = "mod_prob",
        footprintColumns = "footprints",
        arglistFootprints = list(footprints = list(), nucl = list())
    ))), "Can't unambiguously interpret")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "Heatmap", trackData = "mod_prob",
        footprintColumns = "footprints",
        arglistFootprints = list(footprints = 1)
    ))), ".arglistFootprints. entries must be lists")
    expect_error(plotRegion(se = seR2, tracks = list(list(
        trackType = "Heatmap", trackData = "mod_prob",
        footprintColumns = "footprints",
        arglistFootprints = list(list())
    ))), ".names\\(arglistFootprints\\). must not be .NULL.")

    # expected results
    p1 <- plotRegion(se = se, region = "chr1:6948000-6952000")
    p2 <- plotRegion(se = se, tracks = list(list(trackData = "Nvalid",
                                                 trackType = "Point")))
    p3 <- plotRegion(se = se, tracks = list(list(trackData = "FracMod",
                                                 trackType = "Smooth",
                                                 smoothMethod = "smoothSpline")))
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
    expect_warning(
        p9b <- plotRegion(se = seR2, region = "chr1:6935400-6935450",
                          modbaseSpace = TRUE,
                          tracks = list(list(trackType = "GenomicRegions",
                                             trackData = GenomicRanges::GRangesList(
                                                 a = GenomicRanges::GRanges(
                                                     "chr1", IRanges::IRanges(
                                                         6935420, 6935440
                                                     ), "+"
                                                 )
                                             )))),
        "is not allowed if"
    )
    expect_identical(p9, p9b)
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
    p13 <- plotRegion(se = seR2, region = "chr1:6940000-6942000", minCoveredFraction = 0.5,
                      tracks = list(list(trackType = "Lollipop",
                                         trackData = "mod_prob", orderReads = NULL)))
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
        ".windowSize. \\(4\\) must be an odd integer"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6935400-6935450",
                   tracks = list(list(trackType = "Smooth",
                                      trackData = "Nvalid",
                                      smoothMethod = "missing"))),
        "All values in .smoothMethod. must be one of"
    )
    expect_error(
        plotRegion(se = seR2, region = "chr1:6940000-6942000", minCoveredFraction = 1.0,
                   tracks = list(list(trackType = "Lollipop",
                                      trackData = "mod_prob", orderReads = NULL))),
        "No reads retained for plotting"
    )

    expect_true(ggplot2::is_ggplot(p1))
    expect_true(ggplot2::is_ggplot(p2))
    expect_true(ggplot2::is_ggplot(p3))
    expect_true(ggplot2::is_ggplot(p4))
    expect_true(ggplot2::is_ggplot(p5))
    expect_true(ggplot2::is_ggplot(p6))
    expect_true(ggplot2::is_ggplot(p7))
    expect_true(ggplot2::is_ggplot(p8))
    expect_true(ggplot2::is_ggplot(p9))
    expect_true(ggplot2::is_ggplot(p10))
    expect_true(ggplot2::is_ggplot(p11))
    expect_true(ggplot2::is_ggplot(p12))
    expect_true(ggplot2::is_ggplot(p13))
    expect_identical(nrow(p1$data), 4006L)
    expect_identical(nrow(p2$data), 24040L)
    expect_identical(nrow(p3$data), 20000L)
    expect_identical(nrow(p4$data), 2459L)
    expect_identical(nrow(p5$data), 2459L)
    expect_identical(p4$data, p5$data)
    expect_identical(nrow(p6$data), 29104L)
    expect_identical(nrow(p7$data), 29104L)
    expect_identical(nrow(p8$data), 430L)
    expect_length(p9$data, 0L)
    expect_identical(nrow(p13$data), 476L)

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
    expect_identical(ggsave(filename = tmpplot, plot = p13, width = 6, height = 6), tmpplot)
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
    seB <- addFootprints(
        seB,
        wgt = rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170),
                  c(15, 140, 15)),
        name = "nucleosome", verbose = FALSE, thresh = 0.01)
    seB$footprint2 <- seB$nucleosome
    metadata(seB)$readLevelData$colDataColumns <- union(
        metadata(seB)$readLevelData$colDataColumns, "footprint2"
    )

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
    grlNAname <- grl
    names(grlNAname)[2] <- NA

    ## Regions to highlight
    grh <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
                   ranges = IRanges(start = c(6935830, 6935700, 6935870),
                                    end = c(6935850, 6935820, 6935890)))

    ## bigWig files
    bwfiles <- c(s1 = system.file("extdata", "ctcf_chip.bw", package = "footprintR"),
                 s2 = system.file("extdata", "ctcf_chip.bw", package = "footprintR"))

    ## NA in grl name
    expect_error(plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        tracks = list(list(trackData = grlNAname, trackType = "GenomicRegion",
                           colorByStrand = TRUE))))

    ## Repeated names in bwfiles
    expect_error(plotRegion(seB, region = "chr1:6935800-6935900",
                            tracks = list(list(trackData = bwfiles |> stats::setNames(c("s1", "s1")),
                                               trackType = "BigWig"))),
                 "Duplicated file names")
    ## No names for bwfiles
    expect_error(plotRegion(seB, region = "chr1:6935800-6935900",
                            tracks = list(list(trackData = bwfiles |> stats::setNames(NULL),
                                               trackType = "BigWig"))),
                 "must be a named character vector")
    ## Non-existing bwfiles
    expect_error(plotRegion(seB, region = "chr1:6935800-6935900",
                            tracks = list(list(trackData = c(s1 = "missing"),
                                               trackType = "BigWig"))),
                 "Not all bigWig files exist")
    ## Missing colors
    expect_error(plotRegion(seB, region = "chr1:6935800-6935900",
                            tracks = list(list(trackData = bwfiles,
                                               trackType = "BigWig",
                                               colors = c(s1 = "blue")))),
                 "Missing color specification")

    ## Absolute base space, no reference coordinate
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        labelAccuracy = 1e-6,
        tracks = list(list(trackData = "mod_prob", trackType = "Lollipop",
                           size = 2, stroke = 0.25, legendTitle = "6mA",
                           highlightRegions = grh, clustDist = "pearson"),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 3,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE))) +
        plot_layout(heights = c(3, 1, 2)),
        "the standard deviation is zero")
    expect_true(ggplot2::is_ggplot(p))

    ## ... don't color by strand, move labels, add footprints
    expect_warning({
        expect_warning({
            p <- plotRegion(
                seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
                tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                                   legendTitle = "6mA", orderReads = NULL,
                                   orderReverse = TRUE, trackTitle = "Heatmap",
                                   facetBy = NULL, footprintColumns = "nucleosome",
                                   arglistFootprints = list(nucleosome = list(inherit.aes = FALSE))),
                              list(trackData = "mod_prob", trackType = "Lollipop",
                                   legendTitle = "6mA",
                                   orderReads = NULL, trackTitle = "Heatmap",
                                   facetBy = NULL, footprintColumns = c("nucleosome", "footprint2"),
                                   arglistFootprints = list(nucleosome =
                                                                list(fill = "cyan",
                                                                     inherit.aes = FALSE))),
                              list(trackData = grl, trackType = "GenomicRegion",
                                   colorByStrand = FALSE, labelSize = 3,
                                   labelPosition = "above", legendTitle = NULL),
                              list(trackData = "Nvalid", trackType = "Smooth",
                                   showLegend = FALSE,
                                   highlightRegions = grh))) +
                plot_layout(heights = c(3, 3, 1, 2))
        })
    })
    expect_true(ggplot2::is_ggplot(p))

    ## ... define the region including strand
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900:-", modbaseSpace = FALSE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA",
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, footprintColumns = "nucleosome"),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA",
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, footprintColumns = "nucleosome",
                           arglistFootprints = list(fill = "cyan")),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = FALSE, labelSize = 3,
                           labelPosition = "above", legendTitle = NULL),
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE,
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 3, 1, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## ... interpolate, add bigwig tracks
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA",
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = TRUE,
                           linewidthTiles = 0.25,
                           footprintColumns = c("nucleosome", "footprint2")),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = FALSE, labelSize = 3,
                           labelPosition = "above", legendTitle = NULL),
                      list(trackData = bwfiles[1], trackType = "BigWig",
                           highlightRegions = grh, colors = c(s1 = "green")),
                      list(trackData = bwfiles, trackType = "BigWig",
                           yAxisRange = c(0, 10), yAxisLabel = "Score2"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 2, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## referenceCoordinate = left border of plot, squish+interpolate heatmap
    p <- plotRegion(
        seB, region = "chr1:6929237-6929337", modbaseSpace = FALSE,
        referenceCoordinate = 6929237,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", orderReverse = TRUE,
                           trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = TRUE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = bwfiles[1], trackType = "BigWig",
                           highlightRegions = grh),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", facetBy = "sample",
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 2, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## squish + no facet
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = NULL,
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## squish + facet
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", trackTitle = "Heatmap",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = "sample",
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## orderReads = "region" + no facet
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "regionAvg", orderRegion = grh[2],
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "regionAvg", orderRegion = grh[1],
                           facetBy = NULL, size = 2, stroke = 0.5)))
    expect_true(ggplot2::is_ggplot(p))

    ## orderReads = "region" + facet
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "regionAvg", orderRegion = grh[2],
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "regionAvg", orderRegion = grh[1],
                           facetBy = "sample", size = 2, stroke = 0.5)))
    expect_true(ggplot2::is_ggplot(p))

    ## cluster reads using different distance metrics
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", orderRegion = grh[2],
                           windowWidth = 15, clustDist = "euclidean",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", orderRegion = grh[1],
                           windowWidth = 15, clustDist = "euclidean",
                           facetBy = "sample", size = 2, stroke = 0.5)))
    expect_true(ggplot2::is_ggplot(p))

    ## facet, different number of reads per facet - adjust height
    setmp <- subsetReads(seB, c("s1-233e48a7-f379-4dcf-9270-958231125563",
                                "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9",
                                "s1-fc4646ce-66f9-401f-b968-e9b0cda14d61",
                                "s2-274d50aa-f060-4bcf-901e-4cab771295f6"))
    p <- plotRegion(
        setmp, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", trackTitle = "Heatmap",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = "sample",
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## facet, different number of reads per facet - don't adjust height
    ## put GenomicRegion track last - make sure axis text is still shown
    setmp <- subsetReads(seB, c("s1-233e48a7-f379-4dcf-9270-958231125563",
                                "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9",
                                "s1-fc4646ce-66f9-401f-b968-e9b0cda14d61",
                                "s2-274d50aa-f060-4bcf-901e-4cab771295f6"))
    p <- plotRegion(
        setmp, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", trackTitle = "Heatmap",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25, adjustFacetHeight = FALSE),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = "sample",
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome",
                           adjustFacetHeight = FALSE),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL))) +
        plot_layout(heights = c(3, 3, 2, 1))
    expect_true(ggplot2::is_ggplot(p))

    ## cluster, with single read
    setmp <- subsetReads(seB, c("s1-233e48a7-f379-4dcf-9270-958231125563"),
                         prune = TRUE)
    p <- plotRegion(
        setmp, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", trackTitle = "Heatmap",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = "sample",
                           size = 2, stroke = 0.5,
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
            plot_layout(heights = c(3, 1, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## referenceCoordinate outside plot (left) - breaks in scale_cut!
    ## fixed in https://github.com/r-lib/scales/commit/6f2f979a81678c7cd5597b1d18cac78e9cf473c6,
    ## but not yet on CRAN
    # expect_warning(p <- plotRegion(
    #     seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
    #     referenceCoordinate = 6934800,
    #     tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
    #                        legendTitle = "6mA", highlightRegions = grh,
    #                        orderReads = NULL, trackTitle = "Heatmap",
    #                        facetBy = NULL, interpolate = FALSE,
    #                        linewidthTiles = 0.25),
    #                   list(trackData = grl, trackType = "GenomicRegion",
    #                        colorByStrand = TRUE, labelSize = 2,
    #                        labelPosition = "inside", legendTitle = NULL),
    #                   list(trackData = "mod_prob", trackType = "Lollipop",
    #                        legendTitle = "6mA", highlightRegions = grh,
    #                        orderReads = "cluster", facetBy = "sample",
    #                        size = 2, stroke = 0.5),
    #                   list(trackData = "Nvalid", trackType = "PointSmooth",
    #                        showLegend = FALSE, spar = 0.5,
    #                        trackTitle = "Smooth",
    #                        highlightRegions = grh))) +
    #         plot_layout(heights = c(3, 1, 3, 2)),
    #     "the standard deviation is zero")
    # expect_true(ggplot2::is_ggplot(p))

    ## referenceCoordinate outside plot (right)
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6936800,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = "modbase", interpolate = FALSE,
                           linewidthTiles = 0.25,
                           footprintColumns = "nucleosome"),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "squish", facetBy = "sample",
                           size = 2, stroke = 0.5),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 1, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## referenceCoordinate in the middle of plot region
    expect_warning(expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        referenceCoordinate = 6935850,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", trackTitle = "Heatmap",
                           clustDist = "pearson",
                           facetBy = "sample", interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = grl, trackType = "GenomicRegion",
                           colorByStrand = TRUE, labelSize = 2,
                           labelPosition = "inside", legendTitle = NULL),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", facetBy = "sample",
                           clustDist = "pearson",
                           size = 2, stroke = 0.5),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth", groupBy = "modbase",
                           colorBy = "modbase",
                           highlightRegions = grh))) +
            plot_layout(heights = c(3, 1, 3, 2)),
        "the standard deviation is zero"), "the standard deviation is zero")
    expect_true(ggplot2::is_ggplot(p))

    ## modbaseSpace = TRUE, change colors
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", facetBy = "sample",
                           size = 2, stroke = 0.5),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           colors = c(s1 = "forestgreen", s2 = "firebrick1"),
                           highlightRegions = grh))) +
        plot_layout(heights = c(3, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## modbaseSpace = TRUE, footprints -> set modbaseSpace to FALSE
    expect_warning(expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", facetBy = "sample",
                           size = 2, stroke = 0.5, clustDist = "pearson",
                           footprintColumns = "nucleosome"),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth",
                           colors = c(s1 = "forestgreen", s2 = "firebrick1"),
                           highlightRegions = grh))) +
            plot_layout(heights = c(3, 3, 2)),
        "Plotting in `modbaseSpace` is not allowed"), "the standard deviation is zero")
    expect_true(ggplot2::is_ggplot(p))

    ## modbaseSpace = TRUE, bigwig -> set modbaseSpace to FALSE
    expect_warning(p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = bwfiles, trackType = "BigWig"))) +
            plot_layout(heights = c(3, 3)),
        "Plotting in `modbaseSpace` is not allowed if BigWig")
    expect_true(ggplot2::is_ggplot(p))

    ## ... with only smooth, suppressTickLabels = TRUE
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = FALSE,
        suppressTickLabels = TRUE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25, yAxisLabel = "MyAxis"),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, facetBy = "modbase",
                           size = 2, stroke = 0.5, yAxisLabel = "Axis2"),
                      list(trackData = "Nvalid", trackType = "Smooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth", colorBy = "modbase",
                           highlightRegions = grh, yAxisLabel = "Smooth",
                           arglistSmooth = list(linewidth = 2)))) +
        plot_layout(heights = c(3, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

    ## ... change y-axis range
    p <- plotRegion(
        seB, region = "chr1:6935800-6935900", modbaseSpace = TRUE,
        tracks = list(list(trackData = "mod_prob", trackType = "Heatmap",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = NULL, trackTitle = "Heatmap",
                           facetBy = NULL, interpolate = FALSE,
                           linewidthTiles = 0.25),
                      list(trackData = "mod_prob", trackType = "Lollipop",
                           legendTitle = "6mA", highlightRegions = grh,
                           orderReads = "cluster", facetBy = "sample",
                           size = 2, stroke = 0.5),
                      list(trackData = "Nvalid", trackType = "PointSmooth",
                           showLegend = FALSE, spar = 0.5,
                           trackTitle = "Smooth", colorBy = "modbase",
                           highlightRegions = grh, yAxisRange = c(3, 9)))) +
        plot_layout(heights = c(3, 3, 2))
    expect_true(ggplot2::is_ggplot(p))

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
    expect_true(ggplot2::is_ggplot(p))

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
    expect_true(ggplot2::is_ggplot(p))
})

## -------------------------------------------------------------------------- ##
## Checks, helper functions
## -------------------------------------------------------------------------- ##
test_that(".createBaseplotReads works", {
    expect_error(.createFillScale(NULL), "must not be .NULL.")
    expect_error(.createFillScale(1L), "must be of class .character.")

    expect_warning(scl1 <- .createFillScale("error"), "Option .error. does not exist")
    expect_s3_class(scl1, "ScaleContinuous")

    scl2 <- .createFillScale("-cividis")
    scl3 <- .createFillScale("cividis")
    scl4 <- .createFillScale("F")
    scl5 <- .createFillScale(c("red", "yellow", "blue"))

    expect_s3_class(scl2, "ScaleContinuous")
    expect_s3_class(scl3, "ScaleContinuous")
    expect_s3_class(scl4, "ScaleContinuous")
    expect_s3_class(scl5, "ScaleContinuous")

    expect_identical(scl2$palette(seq(0, 1, length.out = 10)),
                     scl3$palette(seq(1, 0, length.out = 10)))

    expect_identical(scl4$palette(c(0, 0.5, 1)),
                     c("#03051A", "#C52D4E", "#FAEBDD"))

    expect_identical(scl5$palette(c(0, 0.5, 1)),
                     unname(apply(grDevices::col2rgb(c("red", "yellow", "blue")),
                                  2, \(x) grDevices::rgb(x[1], x[2], x[3],
                                                         maxColorValue = 255))))
})