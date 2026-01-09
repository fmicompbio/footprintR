# Filter reads

Filter reads

## Usage

``` r
filterReads(
  se,
  assayName = "mod_prob",
  readInfoCol = "readInfo",
  qcCol = "QC",
  minQscore = 0,
  maxEntropy = Inf,
  minSNR = -Inf,
  maxFracLowConf = 1,
  minReadLength = 0,
  minAlignedLength = 0,
  minAlignedFraction = 0,
  minCoveredFraction = 0,
  region = NULL,
  prune = TRUE,
  onlyStats = FALSE,
  removeAllNApos = TRUE
)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- assayName:

  A character scalar providing the name of a read-level assay in `se`.
  This assay will be used to extract read names, as well as to filter
  out any read that is not overlapping any of the positions in the
  object.

- readInfoCol:

  A character scalar providing the name of the column in `colData` that
  contains read info. Can be `NULL` if no such column exists.

- qcCol:

  A character scalar providing the name of the column in `colData` that
  contains quality metrics (calculated by `calcReadStats`). Can be
  `NULL` if no such column exists.

- minQscore:

  A numeric scalar representing the smallest acceptable read-level
  Qscore. Reads with Qscore below this value will be filtered out.

- maxEntropy:

  A numeric scalar representing the largest acceptable read-level
  entropy. Reads without modified-base calls or with entropy above this
  value will be filtered out.

- minSNR:

  A numeric scalar representing the smallest acceptable read Signal to
  Noise ratio (SNR). Reads with SNR below this value will be filtered
  out.

- maxFracLowConf:

  A numeric scalar representing the maximally acceptable fraction of
  low-confidence modified base calls in a read. Reads without
  modified-base calls or with a fraction of low confidence calls greater
  than this value will be filtered out.

- minReadLength:

  A numeric scalar representing the smallest acceptable read length.
  Reads that are shorter than this value will be filtered out.

- minAlignedLength:

  A numeric scalar representing the smallest acceptable aligned length.
  Reads with aligned length shorter than this value will be filtered
  out.

- minAlignedFraction:

  A numeric scalar representing the smallest acceptable aligned fraction
  of a read. Reads where the aligned fraction is smaller than this value
  will be filtered out.

- minCoveredFraction:

  A numeric scalar giving the minimal fraction of `region` that a read
  alignment needs to cover to be retained.

- region:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object with a single region to be used for the `minCoveredFraction`
  filter. Alternatively, the region can be specified as a character
  scalar (e.g. "chr1:1200-1300") that can be coerced into a `GRanges`
  object. If `NULL` (the default), all reads are retained.

- prune:

  A logical scalar. If `TRUE` (the default), samples for which the
  filtering retains none of the reads will be completely removed from
  the returned `SummarizedExperiment` (also from `colData` and from
  assays that do not store read-level data). If `FALSE`, such samples
  are retained (in the assays with read-level data as a zero-column
  `SparseMatrix`).

- onlyStats:

  A logical scalar. If `FALSE` (the default), the `SummarizedExperiment`
  object will be filtered according to the provided thresholds. If
  `TRUE`, the filter statistics are calculated and returned, but the
  object is not subset.

- removeAllNApos:

  A logical scalar. If `TRUE`, remove all positions for which all values
  in `assayName` are `NA` after the read filtering.

## Value

If `onlyStats` is `FALSE`, a filtered `SummarizedExperiment` object. The
metadata of this object contains a slot named `filteredOutReads`, which
tabulate all reads that are filtered out, together with the reason(s)
for exclusion. If `onlyStats` is `TRUE`, only this table is returned.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
library(SummarizedExperiment)
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "footprintR")
se <- readModBam(bamfile = modbamfile, regions = "chr1:6920000-6995000",
                 modbase = "a", verbose = TRUE,
                 BPPARAM = BiocParallel::SerialParam())
#> ℹ extracting base modifications from modBAM files
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [22ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ collapsed 29545 positions to 8439 unique ones
#> ✔ collapsed 29545 positions to 8439 unique ones [165ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
se <- addReadStats(se, name = "QC",
                   BPPARAM = BiocParallel::SerialParam())
#> Warning: Too few points to estimate noise floor (8); raw noise variances are used.

## Filter se
## ... by quality score and aligned length
sefilt <- filterReads(se, minQscore = 14, minAlignedLength = 10000)
## ... by covered fraction of the specified region
sefilt <- filterReads(se, minCoveredFraction = 1, region = "chr1:6930000-6935000")

## Only calculate filter stats
filtstats <- filterReads(se, minQscore = 14, minAlignedLength = 10000,
                         onlyStats = TRUE)
filtstats
#> $s1
#> <7 x 9 SparseMatrix> of type "logical" [nzcount=11 (17%)]:
#>                                          Qscore Entropy ...   AllNA
#> s1-fc4646ce-66f9-401f-b968-e9b0cda14d61   FALSE   FALSE   .   FALSE
#> s1-6cf74134-e550-4c02-bd2b-91385422ee25   FALSE   FALSE   .   FALSE
#> s1-5d45d8d2-d5f5-47ff-a9fa-f3fd6b7bd3c7    TRUE   FALSE   .   FALSE
#> s1-b6fea9db-c92d-4152-9d29-4d021bbc45e8    TRUE   FALSE   .   FALSE
#> s1-49c1e21e-8cb0-415a-aba9-92912219c4bb    TRUE   FALSE   .   FALSE
#> s1-b0b20f04-931f-4f60-b3e4-0ee1f5666a61   FALSE   FALSE   .   FALSE
#> s1-41ca0e97-11b3-454b-9741-bc373e29ef37    TRUE   FALSE   .   FALSE
#> 

## Visualize filter stats in UpSet plot, e.g. with ComplexUpset
if (require(ComplexUpset)) {
    ComplexUpset::upset(as.data.frame(filtstats$s1),
                        intersect = colnames(filtstats$s1))
}
#> Loading required package: ComplexUpset
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the ComplexUpset package.
#>   Please report the issue at
#>   <https://github.com/krassowski/complex-upset/issues>.

```
