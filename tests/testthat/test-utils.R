library(testthat)

## -------------------------------------------------------------------------- ##
## Checks, .assertScalar
## -------------------------------------------------------------------------- ##
test_that(".assertScalar works", {
    expect_error(.assertScalar(1, type = TRUE))
    expect_error(.assertScalar(1, type = 1))
    expect_error(.assertScalar(1, type = c("numeric", "character")))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = 1))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = c(TRUE, FALSE)))

    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_true(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = 0))
    expect_true(.assertScalar(-1, type = "numeric", validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", validValues = c(-2, 0)))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_error(.assertScalar(NA, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 1)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 3:4))
    expect_true(.assertScalar(NULL, type = "numeric", allowNULL = TRUE))
    expect_error(.assertScalar(NULL, type = "numeric", allowNULL = FALSE))
    expect_error(.assertScalar(1, type = "character"))
    expect_error(.assertScalar("x", type = "numeric"))
    expect_error(.assertScalar(FALSE, type = "character"))
    expect_error(.assertScalar(c(1, 2), type = "numeric"))
    test <- "text"
    expect_error(.assertScalar(x = test, type = "numeric"),
                 "'test' must be of class 'numeric")
})

## -------------------------------------------------------------------------- ##
## Checks, .assertVector
## -------------------------------------------------------------------------- ##
test_that(".assertVector works", {
    expect_error(.assertVector(1, type = TRUE))
    expect_error(.assertVector(1, type = 1))
    expect_error(.assertVector(1, type = c("numeric", "character")))
    expect_error(.assertVector(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertVector(1, type = "numeric", allowNULL = 1))
    expect_error(.assertVector(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertVector(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertVector(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    expect_error(.assertVector(1, type = "numeric", len = TRUE))
    expect_error(.assertVector(1, type = "numeric", len = "rng"))
    expect_error(.assertVector(1, type = "numeric", len = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngLen = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngLen = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1:3))

    expect_true(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 1.5)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_error(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 2)))
    expect_error(.assertVector(c(1, 2), type = "numeric", len = 1))
    expect_true(.assertVector(c(1, 2), type = "numeric", len = 2))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngLen = c(3, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(2, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(1, 2)))
    expect_error(.assertVector(c("a", "b"), type = "character", validValues = c("A", "B")))
    expect_true(.assertVector(LETTERS[1:2], type = "character", validValues = LETTERS))
    test <- "text"
    expect_error(.assertVector(x = test, type = "numeric"),
                 "'test' must be of class 'numeric")
})

## -------------------------------------------------------------------------- ##
## Checks, .assertPackagesAvailable
## -------------------------------------------------------------------------- ##
test_that(".assertPackagesAvailable works", {
    testfunc <- function(...) .assertPackagesAvailable(...)
    expect_error(testfunc(1L))
    expect_error(testfunc("test", "error"))
    expect_error(testfunc("test", c(TRUE, FALSE)))

    expect_true(testfunc("base"))
    expect_true(testfunc("githubuser/base"))
    expect_true(testfunc(c("base", "methods")))
    expect_error(testfunc(c("error", "error2")), "BiocManager")
    expect_error(testfunc("error1", suggestInstallation = FALSE), "installed.\n$")
    rm(testfunc)
})

## -------------------------------------------------------------------------- ##
## Checks, .interpolateColumns
## -------------------------------------------------------------------------- ##
test_that(".interpolateColumns works", {
    naa <- SparseArray::NaArray(dim = c(10, 3), type = "double")
    naa[cbind(c(1:3, 8:10, 1, 10), rep(1:3, c(6, 0, 2)))] <-  rep(c(0,1,0,1), c(3, 3, 1, 1))
    colnames(naa) <- paste0("read_", seq.int(ncol(naa)))
    pos1 <- 1:10
    pos2 <- c(1:5, 11:15)

    # valid arguments
    expect_error(.interpolateColumns())
    expect_error(.interpolateColumns(assaydat = naa))
    expect_error(.interpolateColumns(assaydat = unname(naa), pos = pos2))
    expect_error(
        expect_warning(
            expect_warning(
                .interpolateColumns(assaydat = naa, pos = pos, maxgap = "error")
            )
        )
    )

    # expected results
    res1 <- .interpolateColumns(assaydat = naa, pos = pos1, maxgap = Inf)
    res2 <- .interpolateColumns(assaydat = naa, pos = pos2, maxgap = Inf)
    res3 <- .interpolateColumns(assaydat = naa, pos = pos2, maxgap = 10)
    res4 <- .interpolateColumns(assaydat = naa, pos = pos2, maxgap = 5)

    # ... structure
    expect_type(res1, "double")
    expect_type(res2, "double")
    expect_type(res3, "double")
    expect_type(res4, "double")

    expect_identical(dim(res1), c(length(pos1), ncol(naa)))
    expect_identical(dim(res1), c(length(pos2), ncol(naa)))
    expect_identical(dim(res1), c(length(pos2), ncol(naa)))
    expect_identical(dim(res1), c(length(pos2), ncol(naa)))

    expect_identical(colnames(res1), colnames(naa))
    expect_identical(colnames(res2), colnames(naa))
    expect_identical(colnames(res3), colnames(naa))
    expect_identical(colnames(res4), colnames(naa))

    expect_identical(attr(res1, "pos"), seq(min(pos1), max(pos1)))
    expect_identical(attr(res2, "pos"), seq(min(pos2), max(pos2)))
    expect_identical(attr(res3, "pos"), seq(min(pos2), max(pos2)))
    expect_identical(attr(res4, "pos"), seq(min(pos2), max(pos2)))

    # ... content of res1
    expect_identical(res1[attr(res1, "pos") %in% pos1, ][which(as.matrix(SparseArray::is_nonna(naa)), arr.ind = TRUE)],
                     SparseArray::nnavals(naa))
    expect_equal(res1[, "read_1"], c(0, 0, seq(0, 1, by = 0.2), 1, 1))
    expect_equal(res1[, "read_2"], rep(NA_real_, diff(range(pos1)) + 1))
    expect_equal(res1[, "read_3"], seq(0, 1, length.out = diff(range(pos1)) + 1))

    # ... content of res2
    expect_identical(res2[attr(res2, "pos") %in% pos2, ][which(as.matrix(SparseArray::is_nonna(naa)), arr.ind = TRUE)],
                     SparseArray::nnavals(naa))
    expect_equal(res2[, "read_1"], c(0, 0, seq(0, 1, by = 0.1), 1, 1))
    expect_equal(res2[, "read_2"], rep(NA_real_, diff(range(pos2)) + 1))
    expect_equal(res2[, "read_3"], seq(0, 1, length.out = diff(range(pos2)) + 1))

    # ... content of res3
    expect_identical(res3[attr(res3, "pos") %in% pos2, ][which(as.matrix(SparseArray::is_nonna(naa)), arr.ind = TRUE)],
                     SparseArray::nnavals(naa))
    expect_equal(res3[, "read_1"], c(0, 0, seq(0, 1, by = 0.1), 1, 1))
    expect_equal(res3[, "read_2"], rep(NA_real_, diff(range(pos2)) + 1))
    expect_equal(res3[, "read_3"], c(0, rep(NA_real_, length.out = diff(range(pos2)) - 1), 1))

    # ... content of res4
    expect_identical(res4[attr(res4, "pos") %in% pos2, ][which(as.matrix(SparseArray::is_nonna(naa)), arr.ind = TRUE)],
                     SparseArray::nnavals(naa))
    expect_equal(res4[, "read_1"], c(0, 0, 0, rep(NA_real_, diff(range(pos2)) - 5), 1, 1, 1))
    expect_equal(res4[, "read_2"], rep(NA_real_, diff(range(pos2)) + 1))
    expect_equal(res4[, "read_3"], c(0, rep(NA_real_, length.out = diff(range(pos2)) - 1), 1))
})
