suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
    library(GenomicRanges)
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
    v <- .modkitVersion(modkit_bin = NULL)

    if (!is.na(v)) {
        # example data
        modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
        modbamfile2 <- system.file("extdata", "5mC_1_10reads.bam", package = "footprintR")
        reg <- as("chr1:6940000-6955000", "GRanges")
        regs <- as(c("chr1:6940000-6942400", "chr1:6942600-6950000"), "GRanges")

        # invalid arguments
        expect_error(modkitExtract(modkit_bin = "error"))
        expect_error(modkitExtract(bamfile = "error", verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = "error", verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, num_reads = "error",
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, out_extract_table = 1L,
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, out_read_calls = 2L,
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, out_log_file = 3L,
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, modkit_args = 4L,
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, tempdir_base = 5L,
                                   verbose = FALSE))
        expect_error(modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                                   regions = reg, verbose = "error"))

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
        se2 <- getReadDataForRegion(bamfile = modbamfile, region = reg,
                                    arglist.readModkitExtract = list(modbase = "a"))
        expect_identical(se1, se2)
        unlink(c(tmp_tab, tmp_calls, tmp_log))

        # ... one bamfile, no regions
        res2 <- modkitExtract(modkit_bin = NULL, bamfile = modbamfile,
                              regions = NULL,
                              out_extract_table = tmp_tab,
                              out_read_calls = NULL,
                              out_log_file = tmp_log,
                              verbose = FALSE)
        expect_identical(res2, c(`extract-table` = normalizePath(tmp_tab),
                                 `read-calls` = NA,
                                 `run-log` = normalizePath(tmp_log)))
        lns <- readLines(tmp_tab)
        expect_length(lns, 33300L)
        expect_true(grepl("^read_id\tforward_read_position\tref_position\tchrom\tmod_strand\tref_strand", lns[1]))
        lns <- readLines(tmp_log)
        expect_identical(lns[1], "TEST")
        expect_true(grepl("INFO.+processed", lns[length(lns)]))
        unlink(c(tmp_tab, tmp_log))

        # ... one bamfile, two regions
        res3 <- modkitExtract(modkit_bin = NULL, bamfile = modbamfile2,
                              regions = regs,
                              out_extract_table = NULL,
                              out_read_calls = tmp_calls,
                              out_log_file = NULL,
                              verbose = FALSE)
        expect_identical(res3, c(`extract-table` = NA,
                                 `read-calls` = normalizePath(tmp_calls),
                                 `run-log` = NA))
        lns <- readLines(tmp_calls)
        expect_length(lns, 9929L)
        expect_true(grepl("^read_id\tforward_read_position\tref_position\tchrom\tmod_strand\tref_strand", lns[1]))
        suppressMessages(
            expect_message(
                se3 <- readModkitExtract(tmp_calls, modbase = "a", verbose = TRUE)
            )
        )
        expect_s4_class(se3, "RangedSummarizedExperiment")
        expect_identical(dim(se3), c(5496L, 4L))
        expect_true(all(overlapsAny(regs, se3)))
        unlink(c(tmp_calls))
    }
})
