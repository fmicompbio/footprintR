test_that("pileup_modbam_cpp works", {
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
