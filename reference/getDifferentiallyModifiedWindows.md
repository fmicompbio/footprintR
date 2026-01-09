# Analyze counts for sequential windows in a single region

Given a `SummarizedExperiment` with modified and total base counts in
assays `"Nmod"` and `"Nvalid"`, perform a pairwise statistical test for
differential modification.

## Usage

``` r
getDifferentiallyModifiedWindows(
  se,
  assayNameMod = "Nmod",
  assayNameValid = "Nvalid",
  groupCol = "group",
  verbose = FALSE
)
```

## Arguments

- se:

  `SummarizedExperiment`, for example returned by
  `quantifyWindowsInRegion`. It is expected to at least contain assays
  for modified and total counts (given by `assayNameMod` and
  `assayNameValid`) and a `colData` column that defines the groups
  (given by `groupCol`).

- assayNameMod, assayNameValid:

  Character scalars that give the assay names in `se` containing the
  modified and total counts, respectively.

- groupCol:

  Character scalar giving the column in `colData(se)` that defines the
  groups of samples to be compared.

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

The
[`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object constructed from the
[`topTags`](https://rdrr.io/pkg/edgeR/man/topTags.html) output obtained
for the statistical analysis, with additional columns for the average
fraction modified counts in each group, and summary columns named
"dirNegLog10PValue" (the sign of the logFC multiplied with the
-log10(PValue)) and "DeltaFracMod" (the difference of the average
modification fraction in the two groups).

## Author

Panagiotis Papasaikas, Sebastien Smallwood, Charlotte Soneson, Michael
Stadler

## Examples

``` r
modbamfiles <- system.file("extdata",
                           c("6mA_1_10reads.bam", "6mA_1_10reads.bam",
                             "6mA_2_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
se <- quantifyWindowsInRegion(bamfiles = modbamfiles,
                              region = "chr1:6940000-6955000", modbase = "a",
                              BPPARAM = BiocParallel::SerialParam())
se$group <- c("group1", "group1", "group2", "group2")
gr <- getDifferentiallyModifiedWindows(se, groupCol = "group", verbose = TRUE)
#> ℹ calculating library sizes and normalization factors
#> ✔ calculating library sizes and normalization factors [24ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ creating design matrix
#> ✔ creating design matrix [15ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ testing for differential modifications (group2 - group1)
#> ✔ testing for differential modifications (group2 - group1) [124ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
class(gr)
#> [1] "GRanges"
#> attr(,"package")
#> [1] "GenomicRanges"
head(gr)
#> GRanges object with 6 ranges and 9 metadata columns:
#>       seqnames          ranges strand |     logFC    logCPM        LR
#>          <Rle>       <IRanges>  <Rle> | <numeric> <numeric> <numeric>
#>   [1]     chr1 6940000-6940023      * |  6.696961   12.9861 22.480673
#>   [2]     chr1 6940012-6940035      * |  7.018526   13.4242 29.825436
#>   [3]     chr1 6940024-6940047      * |  5.675384   13.7867 13.472295
#>   [4]     chr1 6940036-6940059      * |  0.585622   13.7867  0.166427
#>   [5]     chr1 6940048-6940071      * |  1.214138   13.6466  1.284355
#>   [6]     chr1 6940060-6940083      * | -1.165075   13.5237  1.256727
#>            PValue         FDR dirNegLog10PValue FracMod_group1 FracMod_group2
#>         <numeric>   <numeric>         <numeric>      <numeric>      <numeric>
#>   [1] 2.12269e-06 1.44343e-04          5.673114      0.0000000      0.5454545
#>   [2] 4.72749e-08 6.42938e-06          7.325370      0.0000000      0.5000000
#>   [3] 2.42112e-04 8.23181e-03          3.615984      0.0000000      0.1875000
#>   [4] 6.83307e-01 1.00000e+00          0.165384      0.0588235      0.0833333
#>   [5] 2.57091e-01 1.00000e+00          0.589913      0.0714286      0.1538462
#>   [6] 2.62272e-01 1.00000e+00         -0.581249      0.1666667      0.0769231
#>       DeltaFracMod
#>          <numeric>
#>   [1]    0.5454545
#>   [2]    0.5000000
#>   [3]    0.1875000
#>   [4]    0.0245098
#>   [5]    0.0824176
#>   [6]   -0.0897436
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
