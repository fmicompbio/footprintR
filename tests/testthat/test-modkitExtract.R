suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(SparseArray)
})

## -------------------------------------------------------------------------- ##
## Checks, .modkitExtract
## -------------------------------------------------------------------------- ##
test_that(".modkitVersion works", {
    # invalid arguments
    expect_error(.modkitVersion(modkit_bin = 1L))

    # expected results
    expect_identical(.modkitVersion(modkit_bin = "error"), NA)
    rversion <- .modkitVersion(modkit_bin = file.path(R.home("bin"), "R"))
    expect_type(rversion, "character")
    expect_true(grepl("^R version ", rversion[1]))
})

test_that("modkitExtract works", {
    # check if modkit executable is available
    modkit_version <- .modkitVersion(modkit_bin = NULL)
    skip_if(is.na(modkit_version), "skipping because modkit not found")

    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    modbamfile2 <- system.file("extdata", "5mC_1_10reads.bam", package = "footprintR")
    reg <- as("chr1:6940000-6955000", "GRanges")
    regs <- as(c("chr1:6940000-6942400", "chr1:6942600-6950000"), "GRanges")

    # invalid arguments
    expect_error(modkitExtract(modkit_bin = "error"), "was not found")
    expect_error(modkitExtract(bamfile = "error", verbose = FALSE),
                 "BAM file not found")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = "error", verbose = FALSE),
                 "must be of class 'GRanges'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, num_reads = "error",
                               verbose = FALSE),
                 "must be of class 'numeric'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, out_extract_table = 1L,
                               verbose = FALSE),
                 "must be of class 'character'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, out_read_calls = 2L,
                               verbose = FALSE),
                 "must be of class 'character'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, out_log_file = 3L,
                               verbose = FALSE),
                 "must be of class 'character'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, modkit_args = 4L,
                               verbose = FALSE),
                 "must be of class 'character'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, tempdir_base = 5L,
                               verbose = FALSE),
                 "must be of class 'character'")
    expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                               regions = reg, verbose = "error"),
                 "must be of class 'logical'")

    # expected results
    tmp_tab <- tempfile()
    tmp_calls <- tempfile()
    tmp_log <- tempfile()

    # remark: some modkit console output cannot be suppressed using sink() or capture.output()
    # ... one bamfile, one region
    writeLines("TEST", tmp_log) # write into log file to trigger warning
    suppressMessages({
        expect_message(
            expect_warning(
                res10 <- modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                       regions = reg, num_reads = 10,
                                       out_extract_table = tmp_tab,
                                       out_read_calls = tmp_calls,
                                       out_log_file = tmp_log,
                                       verbose = TRUE)
            )
        )
    })
    expect_identical(res10, c(`extract-table` = normalizePath(tmp_tab),
                              `read-calls` = normalizePath(tmp_calls),
                              `run-log` = normalizePath(tmp_log)))
    lns <- readLines(tmp_tab)
    expect_length(lns, 13025L)
    expect_true(grepl("^read_id\tforward_read_position\tref_position\tchrom\tmod_strand\tref_strand", lns[1]))
    lns <- readLines(tmp_calls)
    expect_length(lns, 13025L)
    expect_true(all(c("chrom", "call_code", "read_id", "ref_strand",
                      "ref_position", "fail", "call_prob") %in%
                        strsplit(x = lns[1], split = "\t")[[1]]))
    lns <- readLines(tmp_log)
    expect_identical(lns[1], "TEST")
    expect_true(grepl("INFO.+processed", lns[length(lns)]))

    se1 <- readModkitExtract(tmp_calls, modbase = "a")
    suppressMessages(
        expect_message(
            se2 <- readModBam(bamfiles = modbamfile, regions = reg,
                              modbase = "a", verbose = TRUE)
        )
    )
    i1 <- overlapsAny(se1, se2)
    i2 <- match(se1[i1], se2)
    # workaround (missing NaArray methods)
    # expect_true(all(start(se1[!i1]) == 0 |
    #                     rowMaxs(assay(se1, "mod_prob")[!i1, ]) <= 0.02))
    expect_true(all(start(se1[!i1]) == 0 |
                        apply(assay(se1, "mod_prob")[!i1, ], 1, max, na.rm = TRUE) <= 0.02))
    expect_identical(rownames(colData(se1)), rownames(colData(se2)))
    expect_identical(rowData(se1)[i1,], rowData(se2)[i2,])
    # workaround (missing NaArray methods)
    # inz <- nzwhich(assay(se1, "mod_prob")[i1, ] > 0 &
    #                    assay(se2, "mod_prob")[i2, ] > 0, arr.ind = TRUE)
    # expect_equal(assay(se1, "mod_prob")[i1, ][inz],
    #              assay(se2, "mod_prob")[i2, ][inz], tolerance = 1e-6)
    inz1 <- nnawhich(assay(se1, "mod_prob")[i1, ] > 0, arr.ind = TRUE)
    inz2 <- nnawhich(assay(se2, "mod_prob")[i2, ] > 0, arr.ind = TRUE)
    inz <- inz1[paste0(inz1[,1], "_", inz1[,2]) %in% paste0(inz2[,1], "_", inz2[,2]), ]
    expect_equal(as.matrix(assay(se1, "mod_prob")[i1, ])[inz],
                 as.matrix(assay(se2, "mod_prob")[i2, ])[inz], tolerance = 1e-6)
    unlink(c(tmp_tab, tmp_calls, tmp_log))

    # ... one bamfile, no regions
    res2 <- modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                          num_reads = 3,
                          regions = NULL,
                          out_extract_table = NULL,
                          out_read_calls = NULL,
                          out_log_file = tmp_log,
                          verbose = FALSE)
    expect_identical(res2, c(`extract-table` = NA,
                             `read-calls` = NA,
                             `run-log` = normalizePath(tmp_log)))
    lns <- readLines(tmp_log)
    expect_true(grepl("INFO.+processed", lns[length(lns)]))
    unlink(c(tmp_tab, tmp_log))

    # ... one bamfile, two regions
    res3 <- modkitExtract(modkit_bin = NULL, bamfile = modbamfile2,
                          regions = regs,
                          out_extract_table = tmp_tab,
                          out_read_calls = tmp_calls,
                          out_log_file = NULL,
                          verbose = FALSE)
    expect_identical(res3, c(`extract-table` = normalizePath(tmp_tab),
                             `read-calls` = normalizePath(tmp_calls),
                             `run-log` = NA))
    lns <- readLines(tmp_tab)
    expect_length(lns, 9816L)
    expect_true(grepl("^read_id\tforward_read_position\tref_position\tchrom\tmod_strand\tref_strand", lns[1]))
    lns <- readLines(tmp_calls)
    expect_length(lns, 9929L)
    expect_true(grepl("^read_id\tforward_read_position\tref_position\tchrom\tmod_strand\tref_strand", lns[1]))
    suppressMessages(
        expect_message(
            se3 <- readModkitExtract(tmp_calls, modbase = "a", verbose = TRUE)
        )
    )
    expect_s4_class(se3, "RangedSummarizedExperiment")
    expect_identical(dim(se3), c(5496L, 1L))
    expect_identical(dim(as.matrix(assay(se3, "mod_prob"))), c(5496L, 4L))
    expect_true(all(overlapsAny(regs, se3)))
    unlink(c(tmp_calls))
})
