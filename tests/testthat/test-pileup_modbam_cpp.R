test_that("pileup_modbam_cpp works", {
    ## example data ------------------------------------------------------------
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                              package = "footprintR")
    extractfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                               package = "footprintR")
    
    bam4 <- system.file("extdata", "6mA_simple.bam", package = "footprintR")
    bam5 <- system.file("extdata", "6mA_nonPrimary.bam", package = "footprintR")
    bam7 <- system.file("extdata", "6mA_mod-issue.bam", package = "footprintR")
    bam8 <- system.file("extdata", "6mA_too-many-mods.bam", package = "footprintR")
    
    ## invalid arguments -------------------------------------------------------
    # ... non-existing bam file
    expect_error(pileup_modbam_cpp(inname_str = "error", regions = "chr1",
                                   modbase = "a", n_threads = 2, level = "read",
                                   mod_prob_thresh = 0.5, verbose = FALSE),
                 "Could not open input file")
    
    # ... no bam index
    tmpbam <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = modbamfile, to = tmpbam))
    expect_error(pileup_modbam_cpp(inname_str = tmpbam, regions = "chr1",
                                   modbase = "a", n_threads = 2,
                                   verbose = FALSE),
                 "Failed to load the index")
    unlink(tmpbam)

    # ... requesting a region that is not contained in the bam header
    expect_error(pileup_modbam_cpp(inname_str = bam4, regions = "chr2",
                                   modbase = "a"),
                 "Failed to get bam iterator")
    
    # ... too many modifications on a single base
    expect_error(pileup_modbam_cpp(inname_str = bam8, regions = "chr1",
                                   modbase = "a", verbose = FALSE),
                 "More modifications than footprintR")
    
    ## expected results --------------------------------------------------------
    # ... run pileup_modbam_cpp
    df <- read.delim(extractfile)
    df$ref_position <- df$ref_position + 1L # modkit-extract has 0-based coordinates
    suppressMessages({
        res1 <- pileup_modbam_cpp(inname_str = modbamfile,
                                  regions = "chr1:6940000-6955000",
                                  modbase = "a", level = "read",
                                  n_threads = 2,
                                  verbose = TRUE)
    })
    res2 <- pileup_modbam_cpp(inname_str = modbamfile, 
                              regions = "chr1:", 
                              modbase = "a", level = "read",
                              n_threads = 1, 
                              verbose = FALSE)
    res3 <- pileup_modbam_cpp(inname_str = modbamfile, 
                              regions = c("chr1", "chr2"), 
                              modbase = "m", level = "read",
                              n_threads = 1, 
                              verbose = FALSE)
    res4 <- pileup_modbam_cpp(inname_str = bam4,
                              regions = "chr1",
                              modbase = "a", level = "read",
                              n_threads = 1,
                              verbose = FALSE)
    res5 <- pileup_modbam_cpp(inname_str = bam5, 
                              regions = "chr1", 
                              modbase = "a", level = "read",
                              n_threads = 1, 
                              verbose = FALSE)
    res6a <- pileup_modbam_cpp(inname_str = modbamfile, 
                               regions = "chr1:6941000-6941001", 
                               modbase = "a", level = "read",
                               n_threads = 1, 
                               verbose = FALSE)
    res6b <- pileup_modbam_cpp(inname_str = modbamfile, 
                               regions = c("chr1:6941000-6941001", "chr1:6928000-6928001"), 
                               modbase = "a", level = "read",
                               n_threads = 1, 
                               verbose = FALSE)
    aln6a <- Rsamtools::scanBam(file = modbamfile,
                                param = Rsamtools::ScanBamParam(
                                    what = "qname",
                                    which = GRanges("chr1:6941000-6941001")
                                ))
    aln6b <- Rsamtools::scanBam(file = modbamfile,
                                param = Rsamtools::ScanBamParam(
                                    what = "qname",
                                    which = GRanges(c("chr1:6941000-6941001", "chr1:6928000-6928001"))
                                ))
    
    # ... results structure
    expect_type(res1, "list")
    expect_type(res2, "list")
    expect_type(res3, "list")
    expect_type(res4, "list")
    expect_type(res5, "list")
    expect_type(res6a, "list")
    expect_type(res6b, "list")
    
    expected_names <- c(
        "chrom", "ref_position", "ref_mod_strand", "mod_prob", "read_id", "read_df")
    expect_named(res1, expected_names)
    expect_named(res2, expected_names)
    expect_named(res3, expected_names)
    expect_named(res4, expected_names)
    expect_named(res5, expected_names)
    expect_named(res6a, expected_names)
    expect_named(res6b, expected_names)
    
    expected_types <- c(
        "character", "integer", "character", "double",
        "character", "list")
    for (i in seq_along(expected_names)) {
        expect_type(res1[[expected_names[i]]], expected_types[i])
        expect_type(res2[[expected_names[i]]], expected_types[i])
        expect_type(res3[[expected_names[i]]], expected_types[i])
        expect_type(res4[[expected_names[i]]], expected_types[i])
        expect_type(res5[[expected_names[i]]], expected_types[i])
        expect_type(res6a[[expected_names[i]]], expected_types[i])
        expect_type(res6b[[expected_names[i]]], expected_types[i])
    }
    
    expect_s3_class(res1[["read_df"]], "data.frame")
    expect_s3_class(res2[["read_df"]], "data.frame")
    expect_s3_class(res3[["read_df"]], "data.frame")
    expect_s3_class(res4[["read_df"]], "data.frame")
    expect_s3_class(res5[["read_df"]], "data.frame")
    expect_s3_class(res6a[["read_df"]], "data.frame")
    expect_s3_class(res6b[["read_df"]], "data.frame")
    
    expected_df_colnames <- c("read_id", "qscore", "read_length", "aligned_length", "variant_label")
    expect_named(res1$read_df, expected_df_colnames)
    expect_named(res2$read_df, expected_df_colnames)
    expect_named(res3$read_df, expected_df_colnames)
    expect_named(res4$read_df, expected_df_colnames)
    expect_named(res5$read_df, expected_df_colnames)
    expect_named(res6a$read_df, expected_df_colnames)
    expect_named(res6b$read_df, expected_df_colnames)
    
    # ... content res1
    expect_identical(res1$read_df$read_id,
                     c("233e48a7-f379-4dcf-9270-958231125563",
                       "d52a5f6a-a60a-4f85-913e-eada84bfbfb9",
                       "92e906ae-cddb-4347-a114-bf9137761a8d"))
    expect_equal(res1$read_df$qscore,
                 c(14.1428003311157, 16.0126991271973, 20.3082008361816))
    expect_identical(res1$read_df$read_length, c(20058L, 11305L, 12277L))
    expect_identical(res1$read_df$aligned_length, c(14801L, 11214L, 12227L))
    expect_identical(res1$read_df$variant_label, rep(NA_character_, 3L))
    expect_identical(unname(sort(as.vector(table(res1$read_id)))),
                     c(3340L,  3597L, 4363L))
    expect_true(all(res1$mod_prob == -1 | (res1$mod_prob >= 0 & res1$mod_prob <= 1.0)))
    for (nm in setdiff(expected_names, "read_df")) {
        expect_length(res1[[nm]], 11300L)
    }
    expect_length(unique(res1$read_id), 3L)
    i1 <- match(paste0(res1$read_id, ":", res1$ref_position),
                paste0(df$read_id, ":", df$ref_position))
    expect_identical(sum(!is.na(i1)), 11183L)
    expect_identical(res1$ref_position[!is.na(i1)],
                     df$ref_position[i1[!is.na(i1)]])
    
    # ... content res2
    expect_identical(res2$read_df$read_id,
                     c("233e48a7-f379-4dcf-9270-958231125563", "d52a5f6a-a60a-4f85-913e-eada84bfbfb9",
                       "fc4646ce-66f9-401f-b968-e9b0cda14d61", "92e906ae-cddb-4347-a114-bf9137761a8d",
                       "6cf74134-e550-4c02-bd2b-91385422ee25", "5d45d8d2-d5f5-47ff-a9fa-f3fd6b7bd3c7",
                       "b6fea9db-c92d-4152-9d29-4d021bbc45e8", "49c1e21e-8cb0-415a-aba9-92912219c4bb",
                       "b0b20f04-931f-4f60-b3e4-0ee1f5666a61", "41ca0e97-11b3-454b-9741-bc373e29ef37"))
    expect_equal(res2$read_df$qscore,
                 c(14.1428003311157, 16.0126991271973, 21.1338005065918,
                   20.3082008361816, 16.0568008422852, 13.3486995697021,
                   13.7178001403809, 12.6245002746582, 16.3353996276855,
                   13.055100440979))
    expect_identical(res2$read_df$read_length,
                     c(20058L, 11305L, 9246L, 12277L, 10041L, 9044L, 9010L, 14736L, 7725L, 7013L))
    expect_identical(res2$read_df$aligned_length,
                     c(14801L, 11214L, 9227L, 12227L, 9968L, 8895L, 8891L, 8174L, 7637L, 6895L))
    expect_identical(res2$read_df$variant_label, rep(NA_character_, 10L))
    expect_identical(unname(as.vector(table(res2$read_id)[res2$read_df$read_id])),
                     c(4363L, 3340L, 2925L, 3597L, 3078L,
                       2720L, 2568L, 2539L, 2412L, 2003L))
    expect_true(
        all(paste0(res1$chrom, ":", res1$ref_position, ":", res1$ref_mod_strand) %in%
                paste0(res2$chrom, ":", res2$ref_position, ":", res2$ref_mod_strand)))
    for (nm in setdiff(expected_names, "read_df")) {
        expect_length(res2[[nm]], 29545L)
    }
    expect_length(unique(res2$read_id), 10L)
    i2 <- match(paste0(res2$read_id, ":", res2$ref_position),
                paste0(df$read_id, ":", df$ref_position))
    expect_identical(sum(!is.na(i2)), 29104L)
    expect_identical(res2$ref_position[!is.na(i2)],
                     df$ref_position[i2[!is.na(i2)]])
    expect_true(all(
        res2$call_code[!is.na(i2)] == df$call_code[i2[!is.na(i2)]] |
            res2$mod_prob[!is.na(i2)] < 0.5))
    
    # ... content res3
    for (nm in setdiff(expected_names, "read_df")) {
        expect_length(res3[[nm]], 0L)
    }
    ## This will still contain all reads - will be filtered out in the 
    ## readModBam R wrapper
    expect_identical(nrow(res3$read_df), 10L)
    
    # ... content of res4
    expect_equal(res4, list(
        chrom = rep("chr1", 8),
        ref_position = c(6940001L, 6940004L, 6940008L, 6940010L, 6940012L, 
                         6940015L, 6940017L, 6940019L),
        ref_mod_strand = c("+", "-", "+", "-", "+", "+", "-", "+"),
        mod_prob = c(0.134765625, 0.318359375, 0.380859375, -1, -1, 
                     0.724609375, -1, 0.998046875),
        read_id = c("artificial-read-1", "artificial-read-2", "artificial-read-1", 
                    "artificial-read-2", "artificial-read-1", "artificial-read-1", 
                    "artificial-read-2", "artificial-read-1"),
        read_df = data.frame(read_id = c("artificial-read-1", "artificial-read-2"),
                             qscore = c(13.4761904761905, 13.24),
                             read_length = c(21L, 25L),
                             aligned_length = c(19L, 23L),
                             variant_label = rep(NA_character_, 2L))))
    
    # ... content of res5
    expect_identical(res5, list(
        chrom = character(0), ref_position = integer(0), 
        ref_mod_strand = character(0), mod_prob = numeric(0),
        read_id = character(0), 
        read_df = data.frame(read_id = character(0), qscore = numeric(0),
                             read_length = integer(0), aligned_length = integer(0),
                             variant_label = character(0))))
    
    # ... content of res6a and res6b (res6a should be a subset of res6b)
    # ... ... check ground truth
    expect_identical(names(aln6a), names(aln6b)[1])
    expect_identical(aln6a[[1]], aln6b[[1]])
    expect_length(aln6a[[1]]$qname, 1L)
    expect_length(intersect(aln6a[[1]]$qname, aln6b[[2]]$qname), 0L)
    # ... ... check return values
    expect_true(all(aln6a[[1]]$qname == res6a$read_id))
    expect_true(aln6b[[1]]$qname %in% res6b$read_id)
    expect_identical(length(res6a$read_id), sum(aln6b[[1]]$qname == res6b$read_id))
    expect_identical(res6a$ref_position[res6a$read_id == aln6a[[1]]$qname],
                     res6b$ref_position[res6b$read_id == aln6b[[1]]$qname])
    idx <- match(paste(res6a$chrom, res6a$ref_position, res6a$ref_mod_strand, res6a$read_id),
                 paste(res6b$chrom, res6b$ref_position, res6b$ref_mod_strand, res6b$read_id))
    expect_true(all(!is.na(idx)))
    expect_identical(res6a$mod_prob, res6b$mod_prob[idx])
    expect_equal(res6a$read_df, res6b$read_df[2, , drop = FALSE],
                 ignore_attr = TRUE)
})

test_that("pileup_modbam_cpp works by comparing to read_modbam_cpp", {
    # helper function
    get_expected_result <- function(fname, reg, mb = "a", mod_prob_thresh = 0.5,
                                    level = "summary") {
        tmp0 <- read_modbam_cpp(inname_str = fname, regions = reg,
                                modbase = mb, n_alns_to_sample = 0,
                                tnames_for_sampling = "chr1",
                                variantRefNames = character(0),
                                variantRefPositions = integer(0),
                                n_threads = 1, verbose = FALSE)
        tmp <- tmp0[c("chrom", "ref_position", "ref_mod_strand", "mod_prob", 
                      "read_id")]
        if (level == "summary") {
            tmp <- tmp |> 
                as.data.frame() |>
                dplyr::group_by(chrom, ref_position, ref_mod_strand) |>
                dplyr::summarise(Nmod = sum(mod_prob >= mod_prob_thresh),
                                 Nvalid = dplyr::n(),
                                 .groups = "drop") |>
                dplyr::arrange(ref_position, ref_mod_strand)
        } else {
            tmp <- tmp |> 
                as.data.frame() |> 
                dplyr::arrange(ref_position, read_id, dplyr::desc(ref_mod_strand))
        }
        list(lst = as.list(tmp), df = tmp0$read_df)
    }

    # reading all alignments in a bam file (summary)
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    reg <- "."
    thresh <- 0.7
    res0 <- get_expected_result(modbamfile, reg, "a", thresh, level = "summary")

    # ... compare to pileup_modbam_cpp return value
    suppressMessages({
        res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                                 modbase = "a", mod_prob_thresh = thresh,
                                 n_threads = 1, verbose = TRUE, level = "summary")
    })
    expect_type(res, "list")
    expect_length(res, 5L)
    expect_named(res, c("chrom", "ref_position", "ref_mod_strand",
                        "Nmod", "Nvalid"))
    expect_identical(res0$lst, res)
    
    # reading all alignments in a bam file (read)
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    reg <- "."
    thresh <- 0.7
    res0 <- get_expected_result(modbamfile, reg, "a", thresh, level = "read")
    
    # ... compare to pileup_modbam_cpp return value
    suppressMessages({
        res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                                 modbase = "a", mod_prob_thresh = thresh,
                                 n_threads = 1, verbose = TRUE, level = "read")
    })
    expect_type(res, "list")
    expect_length(res, 6L)
    expect_named(res, c("chrom", "ref_position", "ref_mod_strand",
                        "mod_prob", "read_id", "read_df"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[1:5]))
    expect_identical(lengths(res0$lst), lengths(res[1:5]))
    res <- res[1:5] |> 
        as.data.frame() |> 
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_mod_strand)) |> 
        as.list()
    expect_identical(res0$lst, res)

    # reading alignments overlapping a region (summary)
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_2_10reads.bam", package = "footprintR")
    reg <- "chr1:6935830-6935830"
    thresh <- 0.3
    res0 <- get_expected_result(modbamfile, reg, "a", thresh, "summary")

    # ... compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                             modbase = "a", mod_prob_thresh = thresh,
                             n_threads = 1, verbose = FALSE, level = "summary")
    expect_type(res, "list")
    expect_length(res, 5L)
    expect_named(res, c("chrom", "ref_position", "ref_mod_strand",
                        "Nmod", "Nvalid"))
    expect_identical(res0$lst, res)
    
    # reading alignments overlapping a region (read)
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_2_10reads.bam", package = "footprintR")
    reg <- "chr1:6935830-6935830"
    thresh <- 0.3
    res0 <- get_expected_result(modbamfile, reg, "a", thresh, "read")
    
    # ... compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                             modbase = "a", mod_prob_thresh = thresh,
                             n_threads = 1, verbose = FALSE, level = "read")
    expect_type(res, "list")
    expect_length(res, 6L)
    expect_named(res, c("chrom", "ref_position", "ref_mod_strand",
                        "mod_prob", "read_id", "read_df"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    res <- res[1:5] |> 
        as.data.frame() |> 
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_mod_strand)) |> 
        as.list()
    expect_identical(res0$lst, res)
})
