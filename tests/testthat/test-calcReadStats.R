test_that("read statistic functions work", {
    # example data
    # ... mod_prob matrix
    npos <- 1000L
    nreads <- 30L
    lagvals <- 12:64
    mat <- NaArray(dim = c(npos, nreads),
                   dimnames = list(paste0("pos", seq.int(npos)),
                                   paste0("r", seq.int(nreads))),
                   type = "double")
    set.seed(123L)
    rlens <- structure(sample(30:200, size = nreads),
                       names = colnames(mat))
    rstarts <- sort(sample(1:800, size = nreads))
    mat[cbind(unlist(lapply(seq.int(nreads), function(i) {
        rstarts[i] + seq.int(rlens[i]) - 1
    })), rep(1:30, rlens))] <- runif(sum(rlens), min = 0, max = 1)
    expect_equal(colSums(is_nonna(mat)), rlens)
    # ... list of non-NA values per read
    ind <- nnawhich(mat, arr.ind = TRUE)
    probList <- split(nnavals(mat), colnames(mat)[ind[, 2]])[colnames(mat)]
    idxList <- split(ind[, 1], colnames(mat)[ind[, 2]])[colnames(mat)]

    expect_identical(lengths(probList), rlens)
    # ... reads to include
    useReads <- sort(sample(nreads, size = nreads - 3L))
    # ... subsets (modified, unmodified)
    matMod <- matUnmod <- mat
    nnavals(matMod)[nnavals(mat) < 0.5] <- NA
    nnavals(matUnmod)[nnavals(mat) >= 0.5] <- NA
    # ... helper function to get call confidence
    .callConf <- function(x) {
        pmax(as.matrix(x), 1 - as.matrix(x))
    }

    # pre-compute SNR, Signal, Noise
    snr_res <- .estimate_snr_probList(probList, idxList)

    # expected values
    argL <- list(probList = probList, idxList = idxList, useReads = useReads,
                 lowConf = 0.7, xrange = lagvals)

    resL <- lapply(allReadStats, function(param) {
        res <- do.call(param, argL)
        expect_length(res, nreads)
        if (param %in% c("ACModProb","PACModProb")) {
            expect_type(res, "list")
            expect_identical(unname(lengths(res)), rep(length(lagvals), nreads))
            expect_equal(
                unname(res[useReads]),
                switch(param,
                       "ACModProb" = lapply(useReads, \(i) {
                           if (length(probList[[i]]) > 64) {
                               acf(probList[[i]], na.action = na.pass,
                                   lag.max = 64, plot = FALSE)$acf[lagvals]
                           } else {
                               rep(0, length(lagvals))
                           }
                       }),
                       "PACModProb" = lapply(useReads, \(i) {
                           if (length(probList[[i]]) > 64) {
                                pacf(probList[[i]], na.action = na.pass,
                                     lag.max = 64, plot = FALSE)$acf[lagvals]
                           } else {
                               rep(0, length(lagvals))
                           }
                       })
                ))
        } else {
            expect_type(res, "double")
            expect_true(all(is.na(res[-useReads])))
            expect_equal(
                res[useReads],
                switch(param,
                    "MeanModProb" = unname(
                        colMeans(mat[, useReads], na.rm = TRUE)),
                    "FracMod" = unname(
                        colMeans(mat[, useReads] > 0.5, na.rm = TRUE)),
                    "MeanConf" = unname(colMeans(.callConf(mat[, useReads]),
                                                 na.rm = TRUE)),
                    "MeanConfUnm" = unname(colMeans(.callConf(matUnmod[, useReads]),
                                                    na.rm = TRUE)),
                    "MeanConfMod" = unname(colMeans(.callConf(matMod[, useReads]),
                                                    na.rm = TRUE)),
                    "FracLowConf" = unname(colMeans(abs(0.5 - mat[, useReads]) < 0.2,
                                                    na.rm = TRUE)),
                    "SEntrModProb" = unlist(lapply(
                        useReads, \(i) {
                            if (length(probList[[i]]) > 64) {
                                sampleEntropy(probList[[i]], 2L, 0.2)
                            } else {
                                NA
                            }})),
                    "IQRModProb" = unname(
                        MatrixGenerics::colIQRs(as.matrix(mat[, useReads]),
                                                na.rm = TRUE)),
                    "sdModProb" = unname(SparseArray::colSds(mat[, useReads], na.rm = TRUE)),
                    "SignalVar" = unname(snr_res$signal[useReads]),
                    "NoiseVar" = unname(snr_res$noise[useReads]),
                    "SNR" = unname(snr_res$snr[useReads])
                ))
        }
    })
})

test_that("calcReadStats works", {
    # example data
    exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                          package = "footprintR")
    reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")
    se <- readModkitExtract(exfile, modbase = "a",
                            BPPARAM = BiocParallel::SerialParam())
    # ... with all-NA reads
    seNA <- se
    mp <- SummarizedExperiment::assay(seNA, "mod_prob")
    mp$s1 <- as.matrix(mp$s1)
    mp$s1[, 2:3] <- NA
    mp$s1 <- SparseArray::NaArray(mp$s1)
    SummarizedExperiment::assays(seNA, withDimnames = FALSE) <- list(mod_prob = mp)

    ## Expected errors
    expect_error(calcReadStats(se, assayName = "error",
                               BPPARAM = BiocParallel::SerialParam()),
                 "must be one of: mod_prob")

    ## No coverage requirement
    rs <- calcReadStats(se, minNobsPpos = 1,
                        stats = allReadStats,
                        BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    expect_length(S4Vectors::metadata(rs), 5L)
    expect_named(S4Vectors::metadata(rs),
                 c("regions", "sequenceContext", "minNobsPpos",
                   "minNobsPread", "Lags"))
    expect_equal(S4Vectors::metadata(rs)$minNobsPpos, 1L)
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 14L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "SEntrModProb", "ACModProb", "PACModProb",
                      "SignalVar", "NoiseVar", "SNR") %in% colnames(qc)))

    expect_equal(qc$MeanModProb,
                 colSums(assay(se)$s1, na.rm = TRUE) /
                     colSums(assay(se)$s1 >= 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 colSums(assay(se)$s1 >= 0.5, na.rm = TRUE) /
                     colSums(assay(se)$s1 >= 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_type(S4Vectors::metadata(qc), "list")

    ## Default coverage requirement
    Nobs <- rowSums(SummarizedExperiment::assay(se)[["s1"]] >= 0, na.rm = TRUE)
    thr <- max(floor(stats::quantile(Nobs, 0.75) -
                         0.5 * stats::IQR(Nobs)), 1L)
    idx <- which(Nobs >= thr)
    rs <- calcReadStats(se, verbose = TRUE, minNobsPpos = thr,
                        stats = allReadStats,
                        BPPARAM = BiocParallel::SerialParam())
    expect_s4_class(rs, "SimpleList")
    expect_length(rs, 1)
    expect_named(rs, "s1")
    expect_type(S4Vectors::metadata(rs), "list")
    expect_length(S4Vectors::metadata(rs), 5L)
    expect_named(S4Vectors::metadata(rs),
                 c("regions", "sequenceContext", "minNobsPpos",
                   "minNobsPread", "Lags"))
    expect_equal(S4Vectors::metadata(rs)$minNobsPpos, thr)
    qc <- rs[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 14L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "SEntrModProb", "ACModProb", "PACModProb",
                      "SignalVar", "NoiseVar", "SNR") %in% colnames(qc)))

    expect_equal(qc$MeanModProb,
                 colSums(assay(se)$s1[idx, ], na.rm = TRUE) /
                     colSums(assay(se)$s1[idx, ] >= 0, na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(qc$FracMod,
                 colSums(assay(se)$s1[idx, ] >= 0.5, na.rm = TRUE) /
                     colSums(assay(se)$s1[idx, ] >= 0, na.rm = TRUE),
                 ignore_attr = TRUE)

    ## Using `regions` and large LagRange
    rs1 <- calcReadStats(se, regions = GenomicRanges::GRanges(
        "chr1", IRanges::IRanges(6935000, 6935100)), LagRange = c(200, 256),
        minNobsPpos = 5, stats = allReadStats,
        BPPARAM = BiocParallel::SerialParam())
    rs2 <- calcReadStats(se, regions = "chr1:6935000-6935100",
                         LagRange = c(200, 256), minNobsPpos = 5,
                         stats = allReadStats,
                         BPPARAM = BiocParallel::SerialParam())
    expect_identical(rs1, rs2)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 15L))
    expect_equal(sum(rs1$s1$MeanModProb), 1.400375383766)
    expect_true(all(vapply(rs1$s1$ACModProb, function(x) all(x == 0), TRUE)))
    expect_true(all(vapply(rs1$s1$PACModProb, function(x) all(x == 0), TRUE)))

    ## Using `sequenceContext`, `minNobsPread` and `stats`
    expect_error(calcReadStats(se, regions = "chr1:6935000-6935100",
                               sequenceContext = c("TAA", "AAA"),
                               BPPARAM = BiocParallel::SerialParam()),
                 "No sequence context found")
    se1 <- addSeqContext(se, sequenceContextWidth = 3,
                         sequenceReference = reffile)
    expect_identical(colnames(rowData(se1)), "sequenceContext")
    rs1 <- calcReadStats(se1, regions = "chr1:6935000-6936000",
                         sequenceContext = c("TAA", "AAA"), minNobsPpos = 5,
                         minNobsPread = 1, stats = "MeanModProb",
                         BPPARAM = BiocParallel::SerialParam())
    rs2 <- calcReadStats(se1, regions = "chr1:6935000-6936000",
                         sequenceContext = "WAA", minNobsPpos = 5,
                         minNobsPread = 1, stats = "MeanModProb",
                         BPPARAM = BiocParallel::SerialParam())
    expect_named(rs1, "s1")
    expect_named(rs2, "s1")
    # ignore metadata()$sequenceContext (expected to differ, explicit vs. IUPAC code)
    meta_names <- setdiff(names(metadata(rs1)), "sequenceContext")
    expect_identical(metadata(rs1)[meta_names], metadata(rs2)[meta_names])
    expect_identical(rs1$s1, rs2$s1)
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 1L))
    expect_equal(sum(rs1$s1$MeanModProb), 0.4934760681446542)

    ## Using input with all-NA reads
    rs1 <- calcReadStats(se = seNA, stats = c("MeanModProb", "ACModProb"),
                         BPPARAM = BiocParallel::SerialParam())
    expect_named(rs1, "s1")
    expect_s4_class(rs1$s1, "DFrame")
    expect_identical(dim(rs1$s1), c(10L, 2L))
    expect_true(all(is.na(rs1$s1$MeanModProb[2:3])))
    expect_equal(sum(rs1$s1$MeanModProb[-(2:3)]), 0.909023830485028)
    expect_true(is.list(rs1$s1$ACModProb))
    expect_identical(lengths(rs1$s1$ACModProb, use.names = FALSE),
                     rep(c(53L, 1L, 53L), c(1, 2, 7)))
})

test_that("addReadStats works", {
    # example data
    exfiles <- system.file("extdata", c("modkit_extract_rc_6mA_1.tsv.gz",
                                        "modkit_extract_rc_6mA_2.tsv.gz"),
                           package = "footprintR")
    se <- readModkitExtract(exfiles, modbase = "a",
                            BPPARAM = BiocParallel::SerialParam())
    se2 <- addReadStats(se, name = "qc2", stats = allReadStats,
                        BPPARAM = BiocParallel::SerialParam())
    se3 <- addReadStats(se, minNobsPread = 2600, name = "qc2",
                        stats = allReadStats,
                        BPPARAM = BiocParallel::SerialParam())

    # expected errors
    expect_error(addReadStats(se, name = -1,
                              BPPARAM = BiocParallel::SerialParam()),
                 "must be of class .character.")
    expect_error(addReadStats(se, name = c("a", "b"),
                              BPPARAM = BiocParallel::SerialParam()),
                 "must have length 1")

    # expected results
    expect_s4_class(se2, "SummarizedExperiment")
    expect_equal(dim(se), dim(se2))
    expect_equal(assay(se), assay(se2))
    expect_null(se2[["QC"]])
    expect_s4_class(se2$qc2, "SimpleList")
    expect_length(se2$qc2, 2L)
    expect_named(se2$qc2, c("s1", "s2"))
    qc <- se2$qc2[["s1"]]
    expect_s4_class(qc, "DFrame")
    expect_equal(nrow(qc), 10L)
    expect_equal(ncol(qc), 14L)
    expect_true(all(c("MeanModProb", "FracMod", "MeanConf", "MeanConfUnm",
                      "MeanConfMod", "FracLowConf", "IQRModProb", "sdModProb",
                      "SEntrModProb", "ACModProb", "PACModProb",
                      "SignalVar", "NoiseVar", "SNR") %in% colnames(qc)))

    expect_identical(metadata(se2$qc2)$minNobsPread, 0)
    expect_identical(metadata(se3$qc2)$minNobsPread, 2600)
    na_rows <- lapply(endoapply(assay(se), function(x) colSums(is_nonna(x))),
                      function(y) which(y < 2600))
    expect_equal(na_rows, list(c(7,8,9,10), c(5,7,8,9,10)), ignore_attr = TRUE)
    expect_identical(se2$qc2$s1[-na_rows$s1,], se3$qc2$s1[-na_rows$s1,])
    expect_identical(se2$qc2$s2[-na_rows$s2,], se3$qc2$s2[-na_rows$s2,])
})
