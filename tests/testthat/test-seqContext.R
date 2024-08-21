suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(Biostrings)
    library(BSgenome)
    library(withr)
})

## -------------------------------------------------------------------------- ##
## Checks, extractSeqContext
## -------------------------------------------------------------------------- ##
test_that("extractSeqContext works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    regions <- GenomicRanges::GRanges(seqnames = "chr1",
                                      ranges = IRanges::IRanges(start = 6957060 - c(4, 2, 0),
                                                                width = 1, names = c("x","y","z")))
    regions2 <- GenomicRanges::GRanges(seqnames = "chr1",
                                       ranges = IRanges::IRanges(start = 1 + c(0, 2, 4),
                                                                 width = 3, names = c("a","b","c")))
    se <- SummarizedExperiment(assays = matrix(1:3, ncol = 1), rowRanges = regions)

    # temporarily install custom BSgenome package
    bsgnmfile <- system.file("extdata", "BSgenome.Mmusculus.footprintR.reference_0.1.0.tar.gz", package = "footprintR")
    rlibdir <- tempfile(pattern = "Rlib")
    dir.create(rlibdir)
    local_libpaths(new = rlibdir, action = "prefix")
    install.packages(bsgnmfile, lib.loc = rlibdir, repos = NULL,
                     quiet = TRUE, verbose = FALSE)
    suppressPackageStartupMessages(suppressWarnings(
        library(BSgenome.Mmusculus.footprintR.reference, lib.loc = rlibdir, quietly = TRUE)
    ))
    gnm <- Biostrings::readDNAStringSet(ref)

    # invalid arguments
    expect_error(extractSeqContext(x = "error"))
    expect_error(extractSeqContext(x = regions, sequence.context.width = -1))
    expect_error(extractSeqContext(x = regions, sequence.context.width = "error"))
    expect_error(extractSeqContext(x = regions, sequence.context.width = 7))
    expect_error(extractSeqContext(x = regions, sequence.context.width = 7, sequence.reference = "error"))

    # expected results
    expect_warning(s1 <- extractSeqContext(x = regions, sequence.context.width = 6, sequence.reference = ref))
    s2 <- extractSeqContext(x = regions, sequence.context.width = 7, sequence.reference = gnm)
    s3 <- extractSeqContext(x = regions, sequence.context.width = 7, sequence.reference = BSgenome.Mmusculus.footprintR.reference)
    s4 <- extractSeqContext(x = unname(regions), sequence.context.width = 7, sequence.reference = gnm)
    s5 <- extractSeqContext(x = regions2, sequence.context.width = 7, sequence.reference = gnm)
    s6 <- extractSeqContext(x = resize(regions2, width = 1L, fix = "center"), sequence.context.width = 7, sequence.reference = gnm)
    s7 <- extractSeqContext(x = se, sequence.context.width = 7, sequence.reference = gnm)
    expect_s4_class(s1, "DNAStringSet")
    expect_s4_class(s2, "DNAStringSet")
    expect_s4_class(s3, "DNAStringSet")
    expect_s4_class(s4, "DNAStringSet")
    expect_s4_class(s5, "DNAStringSet")
    expect_s4_class(s6, "DNAStringSet")
    expect_s4_class(s7, "DNAStringSet")
    expect_identical(as.character(s1), c(x="AAAGGGG", y="AGGGGAN", z="GGGANNN"))
    expect_identical(s1, s2)
    expect_identical(s1, s3)
    expect_identical(unname(s1), s4)
    expect_identical(as.character(s5), c(a="NNNNNNN", b="NNNNNNN", c="NNNNNNN"))
    expect_identical(s5, s6)
    expect_identical(s7, s1)

    # clean up
    detach("package:BSgenome.Mmusculus.footprintR.reference", unload = TRUE,
           character.only = TRUE)
})

## -------------------------------------------------------------------------- ##
## Checks, addSeqContext
## -------------------------------------------------------------------------- ##
test_that("addSeqContext works", {
    # example data
    ref <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    regions <- GenomicRanges::GRanges(seqnames = "chr1",
                                      ranges = IRanges::IRanges(start = 6957060 - c(4, 2, 0),
                                                                width = 1, names = c("x","y","z")))
    se <- SummarizedExperiment(assays = matrix(1:3, ncol = 1), rowRanges = regions)

    # invalid arguments
    expect_error(addSeqContext(x = "error"))
    expect_error(addSeqContext(x = regions, sequence.context.width = -1))
    expect_error(addSeqContext(x = regions, sequence.context.width = 7, sequence.reference = "error"))

    # expected results
    expect_warning(se1 <- addSeqContext(x = se, sequence.context.width = 6, sequence.reference = ref))
    expect_s4_class(se1, "RangedSummarizedExperiment")
    expect_identical(dim(se), dim(se1))
    expect_true("sequence.context" %in% colnames(rowData(se1)))
    expect_identical(as.character(rowData(se1)$sequence.context),
                     c(x = "AAAGGGG", y = "AGGGGAN", z = "GGGANNN"))
})
