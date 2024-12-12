test_that("pileup_modbam_cpp works", {
    # expected results
    modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
    res0 <- read_modbam_cpp(inname_str = modbamfile, regions = "chr1",
                            modbase = "a", n_alns_to_sample = 0,
                            tnames_for_sampling = "chr1",
                            variantRefNames = character(0),
                            variantRefPositions = integer(0),
                            n_threads = 1, verbose = FALSE)
    res0 <- as.data.frame(res0[c("chrom", "ref_position", "mod_prob")]) |>
        dplyr::group_by(chrom, ref_position) |>
        dplyr::summarise(Nmod = sum(mod_prob >= 0.7),
                         Nvalid = dplyr::n(),
                         .groups = "drop") |>
        dplyr::mutate(ref_position = ref_position + 1) |>
        dplyr::arrange(ref_position) |>
        dplyr::rename(ref_name = chrom, ref_pos = ref_position) |>
        dplyr::mutate(Nmod = as.numeric(Nmod),
                      Nvalid = as.numeric(Nvalid)) |>
        as.list()

    # compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile, modbase = "a",
                             mod_prob_thresh = 0.7, n_threads = 1,
                             verbose = TRUE)
    expect_type(res, "list")
    expect_length(res, 4L)
    expect_named(res, c("ref_name", "ref_pos", "Nmod", "Nvalid"))
    expect_identical(res0, res)
})
