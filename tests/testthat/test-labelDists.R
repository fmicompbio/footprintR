suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
})

## -------------------------------------------------------------------------- ##
## Checks, labelDists
## -------------------------------------------------------------------------- ##
test_that("labelDists works", {
    # example data
    labs <- c("--AACACT-",
              "---ACCCT-", # 1vs2: 1/5
              "---ACAC--", # 1vs3: 0/4  2vs3: 1/4
              "-AAAC-TTT", # 1vs4: 2/6  2vs4: 2/5  3vs4: 2/4
              "---------", # 1vs5: 1    2vs5: 1    3vs5: 1    4vs5: 1
              "TTTTTTTTT") # 1vs6: 5/6  2vs6: 4/5  3vs6: 4/4  4vs6: 5/8  5vs6: 1
    nrand <- 100
    randlabs <- do.call(c, lapply(
        seq.int(nrand),
        \(i) paste(sample(c("A","C","G","T","-"), 50, TRUE), collapse = "")))
    randlabs <- c(randlabs, randlabs)

    # faulty arguments
    expect_error(labelDists(), "argument \"labels\" is missing")
    expect_error(labelDists(c("AA", "AAA")), "labels.2. .AAA. does not have 2 characters")

    # expected results
    resA <- labelDists(labs)
    expect_type(resA, "double")
    expect_identical(dim(resA), c(length(labs), length(labs)))
    expect_true(isSymmetric(resA))
    expect_identical(resA[upper.tri(resA)],
                     c(1/5,
                       0/4, 1/4,
                       2/6, 2/5, 2/4,
                       1,   1,   1,   1,
                       5/6, 4/5, 4/4, 5/8, 1))

    resB <- labelDists(labs, minOverlap = 4)
    expect_identical(resA, resB)

    resC <- labelDists(labs, minOverlap = 5)
    expect_type(resC, "double")
    expect_identical(dim(resC), c(length(labs), length(labs)))
    expect_true(isSymmetric(resC))
    expect_identical(resC[upper.tri(resC)],
                     c(1/5,
                       1,   1,
                       2/6, 2/5, 1,
                       1,   1,   1,   1,
                       5/6, 4/5, 1, 5/8, 1))

    resD <- labelDists(randlabs)
    expect_type(resD, "double")
    expect_identical(dim(resD), c(length(randlabs), length(randlabs)))
    expect_true(isSymmetric(resD))
    expect_identical(diag(resD), rep(0.0, length(randlabs)))
    expect_true(min(resD) >= 0 && max(resD) <= 1.0)
    expect_identical(resD[cbind(seq.int(nrand), seq.int(nrand) + nrand)],
                     rep(0.0, nrand))
})
