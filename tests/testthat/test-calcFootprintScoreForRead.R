test_that("calcFootprintScore works", {
    res <- calcFootprintScoreForRead(pos = c(10L, 13L, 21L, 22L, 27L),
                                     pmod = 1:5  / 5,
                                     wgt = c(0.1, .9, .1),
                                     minconf = 0.7,
                                     minweight = 0.05)
    expect_equal(res,
                 data.frame(
                     pos = 10:27,
                     pmod = c(-0.3, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                              NA, 0.3, NA, NA, NA, NA, 0.5),
                     score = c(NA, -0.03, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                               0.03, 0.27, 0.03, NA, NA, 0.05, NA)))
})
