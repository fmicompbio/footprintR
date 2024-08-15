suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
})

## -------------------------------------------------------------------------- ##
## Checks, helper functions
## -------------------------------------------------------------------------- ##
test_that("complement works", {
    expect_error(complement(1))
    expect_error(complement(TRUE))
    expect_identical(complement("A"), "T")
    expect_identical(complement("a"), "T")
    expect_identical(complement("C"), "G")
    expect_identical(complement("c"), "G")
    expect_identical(complement("G"), "C")
    expect_identical(complement("g"), "C")
    expect_identical(complement("T"), "A")
    expect_identical(complement("t"), "A")
    expect_identical(complement("N"), "N")
    expect_identical(complement("n"), "N")
    expect_identical(complement("X"), "N")
})

test_that("get_unmodified_base works", {
    # see https://samtools.github.io/hts-specs/SAMtags.pdf (section 1.7)
    expect_error(get_unmodified_base(1))
    expect_error(get_unmodified_base(TRUE))
    ins <- c("m", "h", "f", "c", "C",
             "g", "e", "b", "T",
             "U",
             "a", "A",
             "o", "G",
             "n", "N", "X"
             )
    out <- rep(c("C", "T", "U", "A", "G", "N"),
               c(5, 4, 1, 2, 2, 3))
    for (i in seq_along(ins)) {
        expect_identical(get_unmodified_base(ins[i]), out[i])
    }
})

test_that("read_modbam works", {
    # example data
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")
    extractfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv",
                               package = "footprintR")

    # invalid arguments
    # ... non-existing bam file
    expect_error(read_modbam("error", "chr1", "a", FALSE))
    # ... no bam index
    tmpbam <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = modbamfile, to = tmpbam))
    expect_error(read_modbam(tmpbam, "chr1", "a", FALSE))
    unlink(tmpbam)

    # expected results
    df <- read.delim(extractfile)
    expect_message(expect_message(expect_message(
        res1 <- read_modbam(modbamfile, "chr1:6940000-6955000", "a", TRUE)
    )))
    res2 <- read_modbam(modbamfile, "chr1:", "a", FALSE)
    res3 <- read_modbam(modbamfile, "chr1:", "m", FALSE)

    # ... structure
    expect_type(res1, "list")
    expect_type(res2, "list")
    expect_type(res3, "list")

    expected_names <- c("read_id", "forward_read_position", "ref_position",
                        "chrom", "ref_strand", "call_code", "canonical_base",
                        "mod_prob")
    expect_named(res1, expected_names)
    expect_named(res2, expected_names)
    expect_named(res3, expected_names)

    expected_types <- c("character", "integer", "integer", "character",
                        "character", "character", "character", "double")
    for (i in seq_along(expected_names)) {
        expect_type(res1[[expected_names[i]]], expected_types[i])
        expect_type(res2[[expected_names[i]]], expected_types[i])
        expect_type(res3[[expected_names[i]]], expected_types[i])
    }

    # ... content res1
    expect_true(all(nchar(res1$call_code) == 1L))
    expect_true(all(res1$canonical_base == "A"))
    expect_true(all(res1$mod_prob == -1 | (res1$mod_prob >= 0 & res1$mod_prob <= 1.0)))
    for (nm in expected_names) {
        expect_length(res1[[nm]], 8350L)
    }
    expect_length(unique(res1$read_id), 3L)
    i1 <- match(paste0(res1$read_id, ":", res1$forward_read_position),
                paste0(df$read_id, ":", df$forward_read_position))
    expect_identical(sum(!is.na(i1)), 6616L)
    expect_identical(res1$ref_position[!is.na(i1)],
                     df$ref_position[i1[!is.na(i1)]])
    expect_true(all(
        res1$call_code[!is.na(i1)] == df$call_code[i1[!is.na(i1)]] |
            res1$mod_prob[!is.na(i1)] < 0.5))
    i1a <- ifelse(!is.na(i1) & res1$mod_prob > 0.5, i1, NA)
    expect_equal(
        ifelse(res1$call_code[!is.na(i1a)] == '-',
               1 - res1$mod_prob[!is.na(i1a)],
               res1$mod_prob[!is.na(i1a)]),
        df$call_prob[i1a[!is.na(i1a)]], tolerance = 1e-6)

    # ... content res2
    expect_true(all(res2$canonical_base == "A"))
    expect_true(
        all(paste0(res1$chrom, ":", res1$ref_position, ":", res1$ref_strand) %in%
            paste0(res2$chrom, ":", res2$ref_position, ":", res2$ref_strand)))
    for (nm in expected_names) {
        expect_length(res2[[nm]], 24005L)
    }
    expect_length(unique(res2$read_id), 10L)
    i2 <- match(paste0(res2$read_id, ":", res2$forward_read_position),
                paste0(df$read_id, ":", df$forward_read_position))
    expect_identical(sum(!is.na(i2)), 21561L)
    expect_identical(res2$ref_position[!is.na(i2)],
                     df$ref_position[i2[!is.na(i2)]])
    expect_true(all(
        res2$call_code[!is.na(i2)] == df$call_code[i2[!is.na(i2)]] |
            res2$mod_prob[!is.na(i2)] < 0.5))
    i2a <- ifelse(!is.na(i2) & res2$mod_prob > 0.5, i2, NA)
    expect_equal(
        ifelse(res2$call_code[!is.na(i2a)] == '-',
               1 - res2$mod_prob[!is.na(i2a)],
               res2$mod_prob[!is.na(i2a)]),
        df$call_prob[i2a[!is.na(i2a)]], tolerance = 1e-6)

    # ... content res3
    for (nm in expected_names) {
        expect_length(res3[[nm]], 0L)
    }

})
