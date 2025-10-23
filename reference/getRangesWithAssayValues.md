# Create a `GRanges` from a `SummarizedExperiment`.

Create a `GRanges` object from `rowRanges` and `assay` of a
`SummarizedExperiment` object.

## Usage

``` r
getRangesWithAssayValues(se, assayName)
```

## Arguments

- se:

  `RangedSummarizedExperiment` object, for example returned by
  `quantifyWindowsInRegion`.

- assayName:

  Character scalar defining the assay to be extracted from `se` and set
  as `mcols(...)` of the return value. Thus, each column of the assay
  will become a column in `mcols(...)`. If missing, the first assay from
  `se` will be used.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object, with the selected assay in its metadata columns.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
modbamfiles <- system.file("extdata",
                           c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
                             "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
se <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                              region = "chr1:6940000-6955000", modbase = "a",
                              BPPARAM = BiocParallel::SerialParam())
gr <- getRangesWithAssayValues(se, assayName = "Nmod")
gr
#> GRanges object with 136 ranges and 4 metadata columns:
#>       seqnames          ranges strand |   Nmod.s1   Nmod.s2   Nmod.s3   Nmod.s4
#>          <Rle>       <IRanges>  <Rle> | <numeric> <numeric> <numeric> <numeric>
#>     1     chr1 6940000-6940023      * |         0         0         6         6
#>     2     chr1 6940012-6940035      * |         0         0         8         8
#>     3     chr1 6940024-6940047      * |         0         0         3         3
#>     4     chr1 6940036-6940059      * |         2         2         1         1
#>     5     chr1 6940048-6940071      * |         2         2         2         2
#>   ...      ...             ...    ... .       ...       ...       ...       ...
#>   132     chr1 6941572-6941595      * |         0         0         0         0
#>   133     chr1 6941584-6941607      * |         0         0         0         0
#>   134     chr1 6941596-6941619      * |         0         0         0         0
#>   135     chr1 6941608-6941631      * |         0         0         0         0
#>   136     chr1 6941620-6941643      * |         0         0         0         0
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
