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


## -------------------------------------------------------------------------- ##
## Checks, read_modbam_cpp
## -------------------------------------------------------------------------- ##
test_that("read_modbam_cpp works", {
    ## example data
    ## -------------------------------------------------------------------------
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")
    extractfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                               package = "footprintR")

    bam4 <- system.file("extdata", "6mA_simple.bam", package = "footprintR")
    bam5 <- system.file("extdata", "6mA_unaligned.bam", package = "footprintR")
    bam7 <- system.file("extdata", "6mA_mod-issue.bam", package = "footprintR")
    bam8 <- system.file("extdata", "6mA_too-many-mods.bam", package = "footprintR")

    ## invalid arguments
    ## -------------------------------------------------------------------------
    # ... non-existing bam file
    expect_error(read_modbam_cpp("error", "chr1", "a", FALSE))

    # ... no bam index
    tmpbam <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = modbamfile, to = tmpbam))
    expect_error(read_modbam_cpp(tmpbam, "chr1", "a", FALSE))
    unlink(tmpbam)

    # ... corrupted bam file
    tmpbam <- tempfile(fileext = ".bam")
    tmpbai <- paste0(tmpbam, ".bai")
    # ... ... copy only part of `modbamfile`
    con_in <- file(modbamfile, "rb")
    data <- readBin(con_in, what = "raw", n = 1e6)
    close(con_in)
    con_out <- file(tmpbam, "wb")
    writeBin(data[seq.int(length(data) - 77)], con_out)
    close(con_out)
    expect_true(file.copy(from = paste0(modbamfile, ".bai"), to = tmpbai))
    expect_error(read_modbam_cpp(tmpbam, "chr1", "a", FALSE))
    unlink(c(tmpbam, tmpbai))

    # ... requesting a region that is not contained in the bam header
    expect_error(read_modbam_cpp(bam4, "chr2", "a"))

    # ... MM/ML tags referring to position beyond read length
    expect_error(read_modbam_cpp(bam7, "chr1", "a", FALSE))

    # ... too many modifications on a single base
    expect_error(read_modbam_cpp(bam8, "chr1", "a", FALSE))

    ## expected results
    ## -------------------------------------------------------------------------
    # ... run read_modbam_cpp
    df <- read.delim(extractfile)
    expect_message(expect_message(expect_message(
        res1 <- read_modbam_cpp(modbamfile, "chr1:6940000-6955000", "a", TRUE)
    )))
    res2 <- read_modbam_cpp(modbamfile, "chr1:", "a", FALSE)
    res3 <- read_modbam_cpp(modbamfile, c("chr1", "chr2"), "m", FALSE)
    res4 <- read_modbam_cpp(bam4, "chr1", "a", FALSE)
    res5 <- read_modbam_cpp(bam5, "chr1", "a", FALSE)

    # ... results structure
    expect_type(res1, "list")
    expect_type(res2, "list")
    expect_type(res3, "list")
    expect_type(res4, "list")
    expect_type(res5, "list")

    expected_names <- c("read_id", "qscore", "forward_read_position",
                        "ref_position", "chrom", "ref_mod_strand", "call_code",
                        "canonical_base", "mod_prob")
    expect_named(res1, expected_names)
    expect_named(res2, expected_names)
    expect_named(res3, expected_names)
    expect_named(res4, expected_names)
    expect_named(res5, expected_names)

    expected_types <- c("character", "double", "integer", "integer",
                        "character", "character", "character", "character",
                        "double")
    for (i in seq_along(expected_names)) {
        expect_type(res1[[expected_names[i]]], expected_types[i])
        expect_type(res2[[expected_names[i]]], expected_types[i])
        expect_type(res3[[expected_names[i]]], expected_types[i])
        expect_type(res4[[expected_names[i]]], expected_types[i])
        expect_type(res5[[expected_names[i]]], expected_types[i])
    }

    # ... content res1
    expect_equal(res1$qscore,
                     rep(c(14.1428003311157, 16.0126991271973, 20.3082008361816),
                         c(4363L, 1750L,  2237L)))
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
    expect_equal(res2$qscore,
                 rep(c(14.1428003311157, 16.0126991271973, 21.1338005065918,
                       20.3082008361816, 16.0568008422852, 13.3486995697021,
                       13.7178001403809, 12.6245002746582, 16.3353996276855,
                       13.055100440979),
                     c(4363L, 1750L, 2925L, 2237L, 3078L,
                       2720L, 1085L, 2539L, 2412L, 896L)))
    expect_true(
        all(paste0(res1$chrom, ":", res1$ref_position, ":", res1$ref_mod_strand) %in%
            paste0(res2$chrom, ":", res2$ref_position, ":", res2$ref_mod_strand)))
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

    # ... content of res4
    expect_equal(res4, list(
        read_id = rep(c("artificial-read-1", "artificial-read-2"), c(5, 3)),
        qscore = rep(c(13.4761904761905, 13.24), c(5, 3)),
        forward_read_position = c(0L, 7L, 10L, 15L, 19L, 21L, 15L, 7L),
        ref_position = c(6940000L, 6940007L, 6940011L, 6940014L, 6940018L,
                         6940003L, 6940009L, 6940016L),
        chrom = rep("chr1", 8),
        ref_mod_strand = rep(c("+", "-"), c(5, 3)),
        call_code = c("a", "a", "-", "a", "a", "a", "-", "-"),
        canonical_base = rep("A", 8),
        mod_prob = c(0.134765625, 0.380859375, -1, 0.724609375, 0.998046875,
                     0.318359375, -1, -1)))

    # ... content of res5
    expect_identical(res5, list(
        read_id = character(0), qscore = numeric(0),
        forward_read_position = integer(0),
        ref_position = integer(0), chrom = character(0),
        ref_mod_strand = character(0), call_code = character(0),
        canonical_base = character(0), mod_prob = numeric(0)))
})
