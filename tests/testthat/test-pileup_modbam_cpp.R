test_that("pileup_modbam_cpp works", {
    # helper function
    get_expected_result <- function(fname, reg, mb = "a", mod_prob_thresh = 0.5) {
        read_modbam_cpp(inname_str = fname, regions = reg,
                        modbase = mb, n_alns_to_sample = 0,
                        tnames_for_sampling = "chr1",
                        variantRefNames = character(0),
                        variantRefPositions = integer(0),
                        n_threads = 1, verbose = FALSE)[c("chrom", "ref_position", "mod_prob")] |>
            as.data.frame() |>
            dplyr::group_by(chrom, ref_position) |>
            dplyr::summarise(Nmod = sum(mod_prob >= mod_prob_thresh),
                             Nvalid = dplyr::n(),
                             .groups = "drop") |>
            dplyr::arrange(ref_position) |>
            as.list()
    }

    # reading all alignments in a bam file
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    reg <- "."
    thresh <- 0.7
    res0 <- get_expected_result(modbamfile, reg, "a", thresh)

    # ... compare to pileup_modbam_cpp return value
    suppressMessages({
        res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                                 modbase = "a", mod_prob_thresh = thresh,
                                 n_threads = 1, verbose = TRUE)
    })
    expect_type(res, "list")
    expect_length(res, 4L)
    expect_named(res, c("chrom", "ref_position", "Nmod", "Nvalid"))
    expect_identical(res0, res)

    # reading alignments overlapping a region
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_2_10reads.bam", package = "footprintR")
    reg <- "chr1:6935830-6935830"
    thresh <- 0.3
    res0 <- get_expected_result(modbamfile, reg, "a", thresh)

    # ... compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile, regions = reg,
                             modbase = "a", mod_prob_thresh = thresh,
                             n_threads = 1, verbose = FALSE)
    expect_type(res, "list")
    expect_length(res, 4L)
    expect_named(res, c("chrom", "ref_position", "Nmod", "Nvalid"))
    expect_identical(res0, res)
})
