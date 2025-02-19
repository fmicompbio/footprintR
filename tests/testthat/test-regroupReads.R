test_that("read regrouping works", {
    # get example data
    modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                            "6mA_2_10reads.bam"),
                               package = "footprintR")
    se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                     modbase = "a", verbose = FALSE,
                     variantPositions = GPos(seqnames = "chr1",
                                             pos = c(6940000, 6940500)),
                     BPPARAM = BiocParallel::SerialParam())
    se <- addReadStats(se, name = "QC", BPPARAM = BiocParallel::SerialParam())
    wgt <- rep(c(0.5, -0.5, 0.5) * c(140/170, 30/170, 140/170), c(15, 140, 15))
    se <- addFootprints(se, wgt, thresh = 0.05, name = "nucl", verbose = FALSE)
    # define read groups
    groups <- list(g1 = c("s1-233e48a7-f379-4dcf-9270-958231125563",
                          "s2-d03efe3b-a45b-430b-9cb6-7e5882e4faf8"),
                   g2 = "s1-92e906ae-cddb-4347-a114-bf9137761a8d",
                   g3 = c("s2-034b625e-6230-4f8d-a713-3a32cd96c298",
                          "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9"))

    # test that functions fail with wrong input
    expect_error(regroupReads(se = "1", readGroups = groups),
                 "must be of class .RangedSummarizedExperiment.")
    expect_error(regroupReads(se = se, readGroups = 1),
                 "must be of class .list.")
    expect_error(regroupReads(se = se, readGroups = unname(groups)),
                 "must not be")
    expect_error(regroupReadsByColData(se = se, colNames = 1, withinSample = TRUE),
                 "must be of class .character.")
    expect_error(regroupReadsByColData(se = se, colNames = "missing", withinSample = TRUE),
                 "must be one of")
    expect_error(regroupReadsByColData(se = se, colNames = "variant_label",
                                       withinSample = 1),
                 "must be of class .logical.")
    expect_error(regroupReadsByColData(se = se, colNames = "variant_label",
                                       withinSample = c(TRUE, FALSE)),
                 "must have length 1")

    # regroup reads based on predefined grouping
    sere <- regroupReads(se, readGroups = groups)
    expect_equal(lapply(assay(sere, "mod_prob"), ncol),
                 list(g1 = 2, g2 = 1, g3 = 2))
    expect_equal(colnames(as.matrix(assay(sere, "mod_prob"))),
                 paste0(rep(names(groups), lengths(groups)), "-",
                        unlist(groups)))
    expect_equal(assayNames(sere), "mod_prob")
    expect_equal(unname(as.matrix(assay(sere, "mod_prob"))),
                 unname(as.matrix(assay(se, "mod_prob"))[, c(1, 5, 3, 4, 2)]))
    tmp <- do.call(rbind, sere$readInfo)
    rownames(tmp) <- sub("^g[0-9]-", "", rownames(tmp))
    expect_equal(tmp,
                 do.call(rbind, se$readInfo)[c(1, 5, 3, 4, 2), ])
    tmp <- do.call(rbind, sere$QC)
    rownames(tmp) <- sub("^g[0-9]-", "", rownames(tmp))
    expect_equal(tmp,
                 do.call(rbind, se$QC)[c(1, 5, 3, 4, 2), ])
    expect_equal(colnames(sere), names(groups))
    expect_equal(rowRanges(se), rowRanges(sere))
    tmp <- do.call(c, unname(se$nucl))[c(1, 5, 3, 4, 2)]
    names(tmp) <- paste0("g", c(1, 1, 2, 3, 3), "-", names(tmp))
    expect_equal(tmp, do.call(c, unname(sere$nucl)))
    expect_equal(lengths(sere$nucl), c(g1 = 2, g2 = 1, g3 = 2))

    # ... identical results (with warning) if nonexistent reads are provided
    groups2 <- groups
    groups2$g1 <- c(groups2$g1, "missing1")
    groups2$g3 <- c(groups2$g3, "missing2")
    expect_warning({
        sere2 <- regroupReads(se, readGroups = groups2)
    }, "The following reads were not found")
    expect_identical(sere, sere2)

    # ... identical results also if non-read-level assays are present
    se2 <- flattenReadLevelAssay(se)
    sere2 <- regroupReads(se2, readGroups = groups)
    expect_identical(sere, sere2)

    # ... fail if there are no read-level assays
    se2 <- flattenReadLevelAssay(se, keepReads = FALSE)
    expect_error(regroupReads(se2, readGroups = groups),
                 "does not contain any read-level assays")

    # empty sample
    expect_warning({
        sere2 <- regroupReads(se, readGroups = c(groups, list(g4 = "missing")))
    }, "The following reads were not found")
    expect_identical(sere, sere2)

    # inconsistent modbase
    se2 <- se
    se2$modbase[2] <- "m"
    expect_error(regroupReads(se2, readGroups = groups),
                 "Some read groups correspond to reads with different modbases")

    # regroup by read annotation (variant label), across samples
    sere <- regroupReadsByColData(se, colNames = "variant_label",
                                  withinSample = FALSE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = do.call(rbind, se$readInfo)$variant_label)
    expect_equal(lapply(assay(sere, "mod_prob"), ncol),
                 list(`G-` = 3, GT = 2))
    expect_equal(colnames(as.matrix(assay(sere, "mod_prob"))),
                 paste0(rep(names(groups2), lengths(groups2)), "-",
                        unlist(groups2)))
    expect_equal(assayNames(sere), "mod_prob")
    expect_equal(unname(as.matrix(assay(sere, "mod_prob"))),
                 unname(as.matrix(assay(se, "mod_prob"))[, c(2, 4, 5, 1, 3)]))
    expect_equal(colnames(sere), names(groups2))
    expect_equal(rowRanges(se), rowRanges(sere))

    # ... within sample
    sere <- regroupReadsByColData(se, colNames = "variant_label",
                                  withinSample = TRUE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = paste0(rep(colnames(se), se$n_reads), "-",
                                do.call(rbind, se$readInfo)$variant_label))
    expect_equal(lapply(assay(sere, "mod_prob"), ncol),
                 list(`s1-G-` = 1, `s1-GT` = 2, `s2-G-` = 2))
    expect_equal(colnames(as.matrix(assay(sere, "mod_prob"))),
                 paste0(rep(names(groups2), lengths(groups2)), "-",
                        unlist(groups2)))
    expect_equal(assayNames(sere), "mod_prob")
    expect_equal(unname(as.matrix(assay(sere, "mod_prob"))),
                 unname(as.matrix(assay(se, "mod_prob"))[, c(2, 1, 3, 4, 5)]))
    expect_equal(colnames(sere), names(groups2))
    expect_equal(rowRanges(se), rowRanges(sere))

    # multiple annotation columns
    se2 <- se
    se2$readInfo <- lapply(se2$readInfo, function(ri) {
        ri$label2 <- ri$variant_label
        ri
    })
    # ... across samples
    sere <- regroupReadsByColData(se2, colNames = c("variant_label", "label2"),
                                  withinSample = FALSE)
    groups2 <- split(x = rownames(do.call(rbind, se2$readInfo)),
                     f = paste0(do.call(rbind, se2$readInfo)$variant_label, "-",
                                do.call(rbind, se2$readInfo)$label2))
    expect_equal(lapply(assay(sere, "mod_prob"), ncol),
                 list(`G--G-` = 3, `GT-GT` = 2))
    expect_equal(colnames(as.matrix(assay(sere, "mod_prob"))),
                 paste0(rep(names(groups2), lengths(groups2)), "-",
                        unlist(groups2)))
    expect_equal(assayNames(sere), "mod_prob")
    expect_equal(unname(as.matrix(assay(sere, "mod_prob"))),
                 unname(as.matrix(assay(se2, "mod_prob"))[, c(2, 4, 5, 1, 3)]))
    expect_equal(colnames(sere), names(groups2))
    expect_equal(rowRanges(se2), rowRanges(sere))

    # ... within sample
    sere <- regroupReadsByColData(se2, colNames = c("variant_label", "label2"),
                                  withinSample = TRUE)
    groups2 <- split(x = rownames(do.call(rbind, se$readInfo)),
                     f = paste0(rep(colnames(se), se2$n_reads), "-",
                                do.call(rbind, se2$readInfo)$variant_label, "-",
                                do.call(rbind, se2$readInfo)$label2))
    expect_equal(lapply(assay(sere, "mod_prob"), ncol),
                 list(`s1-G--G-` = 1, `s1-GT-GT` = 2, `s2-G--G-` = 2))
    expect_equal(colnames(as.matrix(assay(sere, "mod_prob"))),
                 paste0(rep(names(groups2), lengths(groups2)), "-",
                        unlist(groups2)))
    expect_equal(assayNames(sere), "mod_prob")
    expect_equal(unname(as.matrix(assay(sere, "mod_prob"))),
                 unname(as.matrix(assay(se2, "mod_prob"))[, c(2, 1, 3, 4, 5)]))
    expect_equal(colnames(sere), names(groups2))
    expect_equal(rowRanges(se2), rowRanges(sere))
})
