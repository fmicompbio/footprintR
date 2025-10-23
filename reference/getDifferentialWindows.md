# Perform differential analysis on the rows of an assay

Given a `SummarizedExperiment` with at least one assay, use either the
likelihood ratio framework of `edgeR` or the linear model framework of
`limma` to fit a model for each row of the specified assay and extract
statistics corresponding to a specific contrast.

## Usage

``` r
getDifferentialWindows(
  se,
  assayName = "phasingScoreAbs",
  designMatrix,
  contrast,
  method = "limma",
  verbose = FALSE
)
```

## Arguments

- se:

  `SummarizedExperiment`, for example returned by `phasingScoreFourier`.

- assayName:

  Character scalar that gives the assay name in `se` containing the
  values to test.

- designMatrix:

  Numeric matrix providing the design matrix for the statistical
  modeling. Must have row names, corresponding to the sample (column)
  names of `se`.

- contrast:

  Numeric vector providing the contrast to test.

- method:

  Character scalar giving the method to use. Currently supported methods
  are "edgeR" and "limma".

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

The
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object constructed from the
[`topTags`](https://rdrr.io/pkg/edgeR/man/topTags.html) or
[`topTable`](https://rdrr.io/pkg/limma/man/toptable.html) output
obtained for the statistical analysis, with an additional summary column
named "dirNegLog10PValue" (the sign of the logFC multiplied with the
-log10(PValue)).

## Author

Charlotte Soneson

## Examples

``` r
library(GenomicRanges)
library(SummarizedExperiment)
modbamfiles <- system.file("extdata",
                           c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
                             "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
rng <- GRanges("chr1", IRanges(6940000, 6943600))
windowStep <- 180
windowSize <- 4 * windowStep
s <- seq(start(rng), end(rng) - windowSize + 1, by = windowStep)
windowgr <- GRanges(seqnames = seqnames(rng),
                    ranges = IRanges(start = s, width = windowSize))
se <- readModBam(bamfiles = modbamfiles, regions = rng, level = "summary",
                 modbase = "a", trim = TRUE,
                 BPPARAM = BiocParallel::SerialParam())
seFourier <- phasingScoreFourier(se = se, gr = windowgr, numCoef = 5)
seFourier$group <- c("group1", "group1", "group2", "group2")
mm <- model.matrix(~ group, data = colData(seFourier))
gr <- getDifferentialWindows(seFourier, assayName = "phasingScoreAbs",
                             designMatrix = mm,
                             contrast = c(0, 1),
                             method = "limma",
                             verbose = TRUE)
#> Warning: More than half of residual variances are exactly zero: eBayes unreliable
class(gr)
#> [1] "GRanges"
#> attr(,"package")
#> [1] "GenomicRanges"
head(gr)
#> GRanges object with 6 ranges and 7 metadata columns:
#>       seqnames          ranges strand |     logFC   AveExpr         t
#>          <Rle>       <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#>   [1]     chr1 6940000-6940719      * | -10.34650  16.93087  -3271.85
#>   [2]     chr1 6940180-6940899      * | -23.24387  15.51079  -7350.36
#>   [3]     chr1 6940360-6941079      * | -16.91449   8.87515  -5348.83
#>   [4]     chr1 6940540-6941259      * |  -9.32276   4.66138  -2948.12
#>   [5]     chr1 6940720-6941439      * |  -8.24789   4.12394  -2608.21
#>   [6]     chr1 6940900-6941619      * |  -4.59829   2.29915  -1454.11
#>            P.Value    adj.P.Val         B dirNegLog10PValue
#>          <numeric>    <numeric> <numeric>         <numeric>
#>   [1]  4.62534e-95  2.62103e-94   5352491          -94.3349
#>   [2] 5.17201e-107 8.79241e-106  27013842         -106.2863
#>   [3] 2.55469e-102 2.17149e-101  14304976         -101.5927
#>   [4]  1.59810e-93  6.79193e-93   4345680          -92.7964
#>   [5]  1.02911e-91  3.49897e-91   3401368          -90.9875
#>   [6]  4.36260e-83  1.23607e-82   1057203          -82.3603
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
