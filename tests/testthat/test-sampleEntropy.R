suppressPackageStartupMessages({
    library(testthat)
    library(footprintR)
})

## -------------------------------------------------------------------------- ##
## Checks, sampleEntropy
## -------------------------------------------------------------------------- ##
test_that("sampleEntropy works", {
    # example data
    set.seed(1L)
    x <- rnorm(1000)
    y <- sin(seq(0, 2 * pi, length.out = 1000)) + rnorm(1000, sd = 0.1)
    z <- rep(1:100, each = 10) + rnorm(1000, sd = 0.05)
    allNA <- rep(NA, 100)

    # missing arguments
    expect_error(sampleEntropy())
    expect_error(sampleEntropy(data = x))
    expect_error(sampleEntropy(data = x, m = 2))
    expect_error(sampleEntropy(data = x, r = 0.2))

    # expected results
    #   calculated using pracma::sample_entropy, e.g.
    #   pracma::sample_entropy(ts = x, edim = 2, r = 0.5 * sd(x))
    expect_equal(sampleEntropy(x, 2, 0.5), 1.27978397948638)
    expect_equal(sampleEntropy(y, 3, 0.2), 0.496944555201614)
    expect_equal(sampleEntropy(z, 2, 0.2), 0.014575909171002, tolerance = 1e-2)
    expect_identical(sampleEntropy(allNA, 2, 0.2), 0)
})
