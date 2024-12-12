test_that("pileup_modbam_cpp works", {
    # reading all alignments in a bam file
    # ... expected results
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
        dplyr::mutate(ref_position = ref_position + 1L) |>
        dplyr::arrange(ref_position) |>
        as.list()

    # ... compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile,
                             regions = ".", modbase = "a",
                             mod_prob_thresh = 0.7, n_threads = 1,
                             verbose = TRUE)
    expect_type(res, "list")
    expect_length(res, 4L)
    expect_named(res, c("chrom", "ref_position", "Nmod", "Nvalid"))
    expect_identical(res0, res)

    # reading alignments overlapping a region
    # ... expected results
    modbamfile <- system.file("extdata", "6mA_2_10reads.bam", package = "footprintR")
    res0 <- read_modbam_cpp(inname_str = modbamfile, regions = "chr1:6935830-6935830",
                            modbase = "a", n_alns_to_sample = 0,
                            tnames_for_sampling = "chr1",
                            variantRefNames = character(0),
                            variantRefPositions = integer(0),
                            n_threads = 1, verbose = FALSE)
    res0 <- as.data.frame(res0[c("chrom", "ref_position", "mod_prob")]) |>
        dplyr::group_by(chrom, ref_position) |>
        dplyr::summarise(Nmod = sum(mod_prob >= 0.3),
                         Nvalid = dplyr::n(),
                         .groups = "drop") |>
        dplyr::mutate(ref_position = ref_position + 1L) |>
        dplyr::arrange(ref_position) |>
        as.list()

    # ... compare to pileup_modbam_cpp return value
    res <- pileup_modbam_cpp(inname_str = modbamfile,
                             regions = "chr1:6935830-6935830", modbase = "a",
                             mod_prob_thresh = 0.3, n_threads = 1,
                             verbose = FALSE)
    expect_type(res, "list")
    expect_length(res, 4L)
    expect_named(res, c("chrom", "ref_position", "Nmod", "Nvalid"))
    expect_identical(res0, res)
})
