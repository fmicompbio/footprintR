test_that("getAnchorRegions works", {
    # get example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                     modbase = "a", verbose = FALSE)
    se <- addReadsSummary(se, keep.reads = TRUE)

    # check that the function fails with the wrong input
    expect_error(getAnchorRegions(se = "error"), "'se' must be of class")
    expect_error(getAnchorRegions(se = se, assay.type = "missing"),
                 "'assay.type' must be of class 'character'")
    expect_error(getAnchorRegions(se = se, assay.type = c("mod_prob", "Nvalid"),
                                  regionMidpoints = 1),
                 "'regionMidpoints' must be of class 'GPos'")
    expect_error(getAnchorRegions(se = se, assay.type = c("mod_prob", "Nvalid"),
                                  regionMidpoints = "chr1:6920000-6940000"),
                 "all the ranges in the object to coerce to UnstitchedGPos")
    expect_error(getAnchorRegions(se = se, assay.type = c("mod_prob", "Nvalid"),
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = "9"),
                 "'regionWidth' must be of class 'numeric'")
    expect_error(getAnchorRegions(se = se, assay.type = c("mod_prob", "Nvalid"),
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = c(5, 8)),
                 "'regionWidth' must have length 1")
    expect_error(getAnchorRegions(se = se, assay.type = "mod_prob",
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = 9, prune = "TRUE"),
                 "'prune' must be of class 'logical'")
    expect_error(getAnchorRegions(se = se, assay.type = "mod_prob",
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = 9, prune = c(TRUE, FALSE)),
                 "'prune' must have length 1")
    expect_error(getAnchorRegions(se = se, assay.type = "mod_prob",
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = 9, prune = TRUE,
                                  ignore.strand = "FALSE"),
                 "'ignore.strand' must be of class 'logical'")
    expect_error(getAnchorRegions(se = se, assay.type = "mod_prob",
                                  regionMidpoints = "chr1:6930000:-",
                                  regionWidth = 9, prune = TRUE,
                                  ignore.strand = c(TRUE, FALSE)),
                 "'ignore.strand' must have length 1")


    # check that getAnchorRegions works with correct input
    # ... ignore.strand = TRUE
    # ... ... no pruning necessary, prune=TRUE/FALSE should give the same output
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929389:-", "chr1:6935630:-"),
                            regionWidth = 5, prune = FALSE,
                            ignore.strand = TRUE)
    ar2 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929389:-", "chr1:6935630:-"),
                            regionWidth = 5, prune = TRUE,
                            ignore.strand = TRUE)
    expect_identical(ar1, ar2)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), c("s1", "s2"))
    expect_equal(dim(ar1), c(5, 2))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-2, 2))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:6929387-6929391:-"), ignore.strand = TRUE)
    se_r2 <- subsetByOverlaps(se, GRanges("chr1:6935628-6935632:-"), ignore.strand = TRUE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    s1_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1) + length(s1_r2))
    expect_equal(ncol(assay(ar1, "mod_prob")$s2), length(s2_r1) + length(s2_r2))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1),
                 c(paste0("chr1:6929387-6929391:--", colnames(assay(se_r1, "mod_prob")$s1)[s1_r1]),
                   paste0("chr1:6935628-6935632:--", colnames(assay(se_r2, "mod_prob")$s1)[s1_r2])))
    expect_equal(colnames(assay(ar1, "mod_prob")$s2),
                 c(paste0("chr1:6929387-6929391:--", colnames(assay(se_r1, "mod_prob")$s2)[s2_r1]),
                   paste0("chr1:6935628-6935632:--", colnames(assay(se_r2, "mod_prob")$s2)[s2_r2])))
    # check sum of values
    expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE) +
                     sum(assay(se_r2, "mod_prob")$s1[, s1_r2], na.rm = TRUE))
    expect_equal(sum(assay(ar1, "mod_prob")$s2, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s2[, s2_r1], na.rm = TRUE) +
                     sum(assay(se_r2, "mod_prob")$s2[, s2_r2], na.rm = TRUE))
    # check values for Nvalid assay
    expect_equal(assay(ar1, "Nvalid")$s1[, 1][!is.na(assay(ar1, "Nvalid")$s1[, 1])],
                 assay(se_r1, "Nvalid")[, 1], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s1[, 2][!is.na(assay(ar1, "Nvalid")$s1[, 2])],
                 assay(se_r2, "Nvalid")[, 1], ignore_attr = TRUE)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, c("s1", "s2"))
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region,
                 c(rep("chr1:6929387-6929391:-", length(s1_r1)),
                   rep("chr1:6935628-6935632:-", length(s1_r2))))
    expect_equal(ar1$region_mod_prob$s2$id, colnames(assay(ar1, "mod_prob")$s2))
    expect_equal(ar1$region_mod_prob$s2$region,
                 c(rep("chr1:6929387-6929391:-", length(s2_r1)),
                   rep("chr1:6935628-6935632:-", length(s2_r2))))
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, c("s1", "s2"))
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, c("chr1:6929387-6929391:--s1",
                                            "chr1:6935628-6935632:--s1"))
    expect_equal(ar1$region_Nvalid$s1$region,
                 c("chr1:6929387-6929391:-", "chr1:6935628-6935632:-"))
    expect_equal(ar1$region_Nvalid$s2$id, c("chr1:6929387-6929391:--s2",
                                            "chr1:6935628-6935632:--s2"))
    expect_equal(ar1$region_Nvalid$s2$region,
                 c("chr1:6929387-6929391:-", "chr1:6935628-6935632:-"))

    # ... ignore.strand = FALSE
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929389:-", "chr1:6935630:-"),
                            regionWidth = 5, prune = FALSE,
                            ignore.strand = FALSE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), c("s1", "s2"))
    expect_equal(dim(ar1), c(5, 2))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-2, 2))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:6929387-6929391:-"), ignore.strand = FALSE)
    se_r2 <- subsetByOverlaps(se, GRanges("chr1:6935628-6935632:-"), ignore.strand = FALSE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    s1_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1) + length(s1_r2))
    expect_equal(ncol(assay(ar1, "mod_prob")$s2), length(s2_r1) + length(s2_r2))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1),
                 c(paste0("chr1:6929387-6929391:--", colnames(assay(se_r1, "mod_prob")$s1)[s1_r1]),
                   paste0("chr1:6935628-6935632:--", colnames(assay(se_r2, "mod_prob")$s1)[s1_r2])))
    expect_equal(colnames(assay(ar1, "mod_prob")$s2),
                 c(paste0("chr1:6929387-6929391:--", colnames(assay(se_r1, "mod_prob")$s2)[s2_r1]),
                   paste0("chr1:6935628-6935632:--", colnames(assay(se_r2, "mod_prob")$s2)[s2_r2])))
    # check sum of values
    expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE) +
                     sum(assay(se_r2, "mod_prob")$s1[, s1_r2], na.rm = TRUE))
    expect_equal(sum(assay(ar1, "mod_prob")$s2, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s2[, s2_r1], na.rm = TRUE) +
                     sum(assay(se_r2, "mod_prob")$s2[, s2_r2], na.rm = TRUE))
    # check values for Nvalid assay
    expect_equal(assay(ar1, "Nvalid")$s1[, 1][!is.na(assay(ar1, "Nvalid")$s1[, 1])],
                 assay(se_r1, "Nvalid")[, 1], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s1[, 2][!is.na(assay(ar1, "Nvalid")$s1[, 2])],
                 assay(se_r2, "Nvalid")[, 1], ignore_attr = TRUE)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, c("s1", "s2"))
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region,
                 c(rep("chr1:6929387-6929391:-", length(s1_r1)),
                   rep("chr1:6935628-6935632:-", length(s1_r2))))
    expect_equal(ar1$region_mod_prob$s2$id, colnames(assay(ar1, "mod_prob")$s2))
    expect_equal(ar1$region_mod_prob$s2$region,
                 c(rep("chr1:6929387-6929391:-", length(s2_r1)),
                   rep("chr1:6935628-6935632:-", length(s2_r2))))

    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, c("s1", "s2"))
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, c("chr1:6929387-6929391:--s1",
                                            "chr1:6935628-6935632:--s1"))
    expect_equal(ar1$region_Nvalid$s1$region,
                 c("chr1:6929387-6929391:-", "chr1:6935628-6935632:-"))
    expect_equal(ar1$region_Nvalid$s2$id, c("chr1:6929387-6929391:--s2",
                                            "chr1:6935628-6935632:--s2"))
    expect_equal(ar1$region_Nvalid$s2$region,
                 c("chr1:6929387-6929391:-", "chr1:6935628-6935632:-"))

    # ... ignore.strand = TRUE, in a region where the same position is covered
    #     by reads on both strands
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929338:-", "chr1:6938997:+"),
                            regionWidth = 5, prune = FALSE,
                            ignore.strand = TRUE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), c("s1", "s2"))
    expect_equal(dim(ar1), c(5, 2))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-2, 2))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:6929336-6929340:-"), ignore.strand = TRUE)
    se_r2 <- subsetByOverlaps(se, GRanges("chr1:6938995-6938999:+"), ignore.strand = TRUE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    s1_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r2 <- which(SparseArray::colSums(assay(se_r2, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1) + length(s1_r2))
    expect_equal(ncol(assay(ar1, "mod_prob")$s2), length(s2_r1) + length(s2_r2))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1),
                 c(paste0("chr1:6929336-6929340:--", colnames(assay(se_r1, "mod_prob")$s1)[s1_r1]),
                   paste0("chr1:6938995-6938999:+-", colnames(assay(se_r2, "mod_prob")$s1)[s1_r2])))
    expect_equal(colnames(assay(ar1, "mod_prob")$s2),
                 c(paste0("chr1:6929336-6929340:--", colnames(assay(se_r1, "mod_prob")$s2)[s2_r1]),
                   paste0("chr1:6938995-6938999:+-", colnames(assay(se_r2, "mod_prob")$s2)[s2_r2])))
    # check sum of values
    # TODO: currently not working
    # expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
    #              sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE) +
    #                  sum(assay(se_r2, "mod_prob")$s1[, s1_r2], na.rm = TRUE))
    # expect_equal(sum(assay(ar1, "mod_prob")$s2, na.rm = TRUE),
    #              sum(assay(se_r1, "mod_prob")$s2[, s2_r1], na.rm = TRUE) +
    #                  sum(assay(se_r2, "mod_prob")$s2[, s2_r2], na.rm = TRUE))
    # # check sum of values for Nvalid assay
    # expect_equal(sum(assay(ar1, "Nvalid")$s1[, 1], na.rm = TRUE),
    #              sum(assay(se_r1, "Nvalid")[, 1]), ignore_attr = TRUE)
    # expect_equal(sum(assay(ar1, "Nvalid")$s1[, 2], na.rm = TRUE),
    #              sum(assay(se_r2, "Nvalid")[, 1]), ignore_attr = TRUE)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, c("s1", "s2"))
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region,
                 c(rep("chr1:6929336-6929340:-", length(s1_r1)),
                   rep("chr1:6938995-6938999:+", length(s1_r2))))
    expect_equal(ar1$region_mod_prob$s2$id, colnames(assay(ar1, "mod_prob")$s2))
    expect_equal(ar1$region_mod_prob$s2$region,
                 c(rep("chr1:6929336-6929340:-", length(s2_r1)),
                   rep("chr1:6938995-6938999:+", length(s2_r2))))
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, c("s1", "s2"))
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, c("chr1:6929336-6929340:--s1",
                                            "chr1:6938995-6938999:+-s1"))
    expect_equal(ar1$region_Nvalid$s1$region,
                 c("chr1:6929336-6929340:-", "chr1:6938995-6938999:+"))
    expect_equal(ar1$region_Nvalid$s2$id, c("chr1:6929336-6929340:--s2",
                                            "chr1:6938995-6938999:+-s2"))
    expect_equal(ar1$region_Nvalid$s2$region,
                 c("chr1:6929336-6929340:-", "chr1:6938995-6938999:+"))

    # ... different region, where only s1 has overlapping reads
    # ... ... without pruning
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929015:+"),
                            regionWidth = 8, prune = FALSE,
                            ignore.strand = FALSE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), c("s1", "s2"))
    expect_equal(dim(ar1), c(8, 2))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-3, 4))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:6929012-6929019:+"), ignore.strand = FALSE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1))
    expect_equal(ncol(assay(ar1, "mod_prob")$s2), length(s2_r1))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1),
                 paste0("chr1:6929012-6929019:+-", colnames(assay(se_r1, "mod_prob")$s1)[s1_r1]))
    expect_equal(colnames(assay(ar1, "mod_prob")$s2), character(0))
    # check sum of values
    expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE))
    expect_equal(sum(assay(ar1, "mod_prob")$s2, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s2[, s2_r1], na.rm = TRUE))
    # check values for Nvalid assay
    expect_equal(assay(ar1, "Nvalid")$s1[, 1][!is.na(assay(ar1, "Nvalid")$s1[, 1])],
                 assay(se_r1, "Nvalid")[, 1], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s2[, 1][!is.na(assay(ar1, "Nvalid")$s2[, 1])],
                 assay(se_r1, "Nvalid")[, 2], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s2[, 1][!is.na(assay(ar1, "Nvalid")$s2[, 1])],
                 c(0, 0, 0), ignore_attr = TRUE)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, c("s1", "s2"))
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region,
                 rep("chr1:6929012-6929019:+", length(s1_r1)))
    expect_equal(ar1$region_mod_prob$s2$id, colnames(assay(ar1, "mod_prob")$s2))
    expect_equal(ar1$region_mod_prob$s2$region,
                 rep("chr1:6929012-6929019:+", length(s2_r1)))
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, c("s1", "s2"))
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, "chr1:6929012-6929019:+-s1")
    expect_equal(ar1$region_Nvalid$s1$region, "chr1:6929012-6929019:+")
    expect_equal(ar1$region_Nvalid$s2$id, "chr1:6929012-6929019:+-s2")
    expect_equal(ar1$region_Nvalid$s2$region, "chr1:6929012-6929019:+")

    # ... ... with pruning
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:6929015:+"),
                            regionWidth = 8, prune = TRUE,
                            ignore.strand = FALSE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), "s1")
    expect_equal(dim(ar1), c(8, 1))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-3, 4))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:6929012-6929019:+"), ignore.strand = FALSE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1),
                 paste0("chr1:6929012-6929019:+-", colnames(assay(se_r1, "mod_prob")$s1)[s1_r1]))
    # check sum of values
    expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE))
    # check values for Nvalid assay
    expect_equal(assay(ar1, "Nvalid")$s1[, 1][!is.na(assay(ar1, "Nvalid")$s1[, 1])],
                 assay(se_r1, "Nvalid")[, 1], ignore_attr = TRUE)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, "s1")
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region,
                 rep("chr1:6929012-6929019:+", length(s1_r1)))
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, "s1")
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, "chr1:6929012-6929019:+-s1")
    expect_equal(ar1$region_Nvalid$s1$region, "chr1:6929012-6929019:+")

    # ... region where neither sample has any overlapping reads
    # ... ... with pruning
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:692915:+"),
                            regionWidth = 8, prune = TRUE,
                            ignore.strand = FALSE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), character(0))
    expect_equal(dim(ar1), c(8, 0))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-3, 4))
    # get reads in each sample that overlap each of the regions
    expect_null(assay(ar1, "mod_prob")$s1)
    expect_null(assay(ar1, "mod_prob")$s2)
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_length(ar1$region_mod_prob, 0)
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_length(ar1$region_Nvalid, 0)

    # ... ... without pruning
    ar1 <- getAnchorRegions(se, assay.type = c("mod_prob", "Nvalid"),
                            regionMidpoints = c("chr1:692915:+"),
                            regionWidth = 8, prune = FALSE,
                            ignore.strand = FALSE)
    expect_s4_class(ar1, "SummarizedExperiment")
    expect_equal(colnames(ar1), c("s1", "s2"))
    expect_equal(dim(ar1), c(8, 2))
    expect_equal(SummarizedExperiment::assayNames(ar1), c("mod_prob", "Nvalid"))
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "Nvalid"), "DataFrame")
    expect_s4_class(SummarizedExperiment::assay(ar1, "mod_prob")$s1, "NaArray")
    expect_true(is.matrix(SummarizedExperiment::assay(ar1, "Nvalid")$s1))
    expect_equal(SummarizedExperiment::rowData(ar1)$relpos, seq(-3, 4))
    # get reads in each sample that overlap each of the regions
    se_r1 <- subsetByOverlaps(se, GRanges("chr1:692912-692919:+"), ignore.strand = FALSE)
    s1_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s1 >= 0, na.rm = TRUE) > 0)
    s2_r1 <- which(SparseArray::colSums(assay(se_r1, "mod_prob")$s2 >= 0, na.rm = TRUE) > 0)
    expect_equal(ncol(assay(ar1, "mod_prob")$s1), length(s1_r1))
    expect_equal(ncol(assay(ar1, "mod_prob")$s2), length(s2_r1))
    expect_equal(colnames(assay(ar1, "mod_prob")$s1), character(0))
    expect_equal(colnames(assay(ar1, "mod_prob")$s2), character(0))
    # check sum of values
    expect_equal(sum(assay(ar1, "mod_prob")$s1, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s1[, s1_r1], na.rm = TRUE))
    expect_equal(sum(assay(ar1, "mod_prob")$s2, na.rm = TRUE),
                 sum(assay(se_r1, "mod_prob")$s2[, s2_r1], na.rm = TRUE))
    # check values for Nvalid assay
    expect_equal(assay(ar1, "Nvalid")$s1[, 1][!is.na(assay(ar1, "Nvalid")$s1[, 1])],
                 assay(se_r1, "Nvalid")[, 1], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s1[, 1], rep(NA_real_, 8))
    expect_equal(assay(ar1, "Nvalid")$s2[, 1][!is.na(assay(ar1, "Nvalid")$s2[, 1])],
                 assay(se_r1, "Nvalid")[, 2], ignore_attr = TRUE)
    expect_equal(assay(ar1, "Nvalid")$s2[, 1], rep(NA_real_, 8))
    # colData
    expect_named(colData(ar1), c("sample", "region_mod_prob", "region_Nvalid"))
    expect_s4_class(colData(ar1), "DataFrame")
    expect_s4_class(ar1$region_mod_prob, "List")
    expect_named(ar1$region_mod_prob, c("s1", "s2"))
    expect_s4_class(ar1$region_mod_prob$s1, "DataFrame")
    expect_equal(ar1$region_mod_prob$s1$id, colnames(assay(ar1, "mod_prob")$s1))
    expect_equal(ar1$region_mod_prob$s1$region, character(0))
    expect_equal(ar1$region_mod_prob$s2$id, colnames(assay(ar1, "mod_prob")$s2))
    expect_equal(ar1$region_mod_prob$s2$region, character(0))
    expect_s4_class(ar1$region_Nvalid, "List")
    expect_named(ar1$region_Nvalid, c("s1", "s2"))
    expect_s4_class(ar1$region_Nvalid$s1, "DataFrame")
    expect_equal(ar1$region_Nvalid$s1$id, "chr1:692912-692919:+-s1")
    expect_equal(ar1$region_Nvalid$s1$region, "chr1:692912-692919:+")
    expect_equal(ar1$region_Nvalid$s2$id, "chr1:692912-692919:+-s2")
    expect_equal(ar1$region_Nvalid$s2$region, "chr1:692912-692919:+")
})
