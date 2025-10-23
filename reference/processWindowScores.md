# Process scores of sequential genomic windows.

Given scores or statistical estimates for sequential windows, process
them according to the `scoreAction` argument. Typical actions are to
smooth the scores and possibly identify regions of interest by fusing
consistent neighboring windows along the genome.

## Usage

``` r
processWindowScores(
  x,
  scoreCol = "dirNegLog10PValue",
  scoreAction = c("smoothFuse", "smooth", "pass", "select"),
  minperiod = 3,
  thresh = 3,
  maxGap = 50,
  verbose = FALSE
)
```

## Arguments

- x:

  `GRanges` object with window-based scores or estimates. Ranges
  correspond to windows, and columns in `mcols(x)` to scores.

- scoreCol:

  Character scalar giving the column name in `mcols(x)` to use for the
  analysis. For `scanForHighScoringRegions`, it can be a vector, in
  which case `processWindowScores` will be run for each of them and the
  results will be returned as a list.

- scoreAction:

  A character scalar indicating how to process the scores. Currently
  supported values are:

  pass

  :   : Do not modify window scores (pass-through). Note that if
      `mcols(x)` has multiple columns, they will all be passed through
      (i.e., `scoreCol` is ignored). To return only `scoreCol`, use
      `select`.

  select

  :   : Pass through only `scoreCol`.

  smooth

  :   : Smooth window scores (`minperiod` argument).

  smoothFuse (default)

  :   : Smooth window scores, identify high scoring elements (`thresh`
      argument) and fuse nearby high-scoring elements (`maxGap`
      argument).

- minperiod:

  Numeric scalar that defines the low-pass filter parameter used to
  smooth the scores for segmentation. `minperiod` gives the minimal
  period (in number of windows) for the critical frequency in the score
  signal that should be retained. The default value is suitable to
  retain signals that occur in three neighboring windows. The filtering
  is performed using
  [`filtfilt`](https://rdrr.io/pkg/signal/man/filtfilt.html) and thus
  requires the `signal` package to be installed.

- thresh:

  A numeric scalar giving the minimal absolute window score (after
  smoothing, see `minperiod` argument) defining a region of interest.
  Higher values make the region detection more stringent.

- maxGap:

  Numeric scalar giving the maximal gap between neighboring windows, in
  base pairs, from the end of the first to the start of the next, to be
  fused into a single region of interest.

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

A [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object. For `scoreAction` equal to `"pass"` or `"smooth"`, the returned
object has the same dimensions as the input object `x`. For
`scoreAction="smoothFuse"`, the returned object will only contain
thresholded, fused regions.

## Author

Sebastien Smallwood, Charlotte Soneson, Michael Stadler

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
gr <- getDifferentiallyModifiedWindows(se, groupCol = "group")
gr
#> GRanges object with 136 ranges and 9 metadata columns:
#>         seqnames          ranges strand |     logFC    logCPM          LR
#>            <Rle>       <IRanges>  <Rle> | <numeric> <numeric>   <numeric>
#>     [1]     chr1 6940000-6940023      * |  6.696961   12.9861   22.480673
#>     [2]     chr1 6940012-6940035      * |  7.018526   13.4242   29.825436
#>     [3]     chr1 6940024-6940047      * |  5.675384   13.7867   13.472295
#>     [4]     chr1 6940036-6940059      * |  0.585622   13.7867    0.166427
#>     [5]     chr1 6940048-6940071      * |  1.214138   13.6466    1.284355
#>     ...      ...             ...    ... .       ...       ...         ...
#>   [132]     chr1 6941572-6941595      * |   5.67334   11.9387 2.48245e-07
#>   [133]     chr1 6941584-6941607      * |   5.45565   11.8392 2.08965e-07
#>   [134]     chr1 6941596-6941619      * |   5.86244   12.0319 2.87441e-07
#>   [135]     chr1 6941608-6941631      * |   5.67334   11.9387 2.48245e-07
#>   [136]     chr1 6941620-6941643      * |   4.48815   11.4913 2.93492e-07
#>              PValue         FDR dirNegLog10PValue FracMod_group1 FracMod_group2
#>           <numeric>   <numeric>         <numeric>      <numeric>      <numeric>
#>     [1] 2.12269e-06 1.44343e-04          5.673114      0.0000000      0.5454545
#>     [2] 4.72749e-08 6.42938e-06          7.325370      0.0000000      0.5000000
#>     [3] 2.42112e-04 8.23181e-03          3.615984      0.0000000      0.1875000
#>     [4] 6.83307e-01 1.00000e+00          0.165384      0.0588235      0.0833333
#>     [5] 2.57091e-01 1.00000e+00          0.589913      0.0714286      0.1538462
#>     ...         ...         ...               ...            ...            ...
#>   [132]    0.999602           1       0.000172683              0            NaN
#>   [133]    0.999635           1       0.000158431              0            NaN
#>   [134]    0.999572           1       0.000185820              0            NaN
#>   [135]    0.999602           1       0.000172683              0            NaN
#>   [136]    0.999568           1       0.000187766              0            NaN
#>         DeltaFracMod
#>            <numeric>
#>     [1]    0.5454545
#>     [2]    0.5000000
#>     [3]    0.1875000
#>     [4]    0.0245098
#>     [5]    0.0824176
#>     ...          ...
#>   [132]          NaN
#>   [133]          NaN
#>   [134]          NaN
#>   [135]          NaN
#>   [136]          NaN
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

grFused <- processWindowScores(x = gr, scoreCol = "logFC", thresh = 5.0)
grFused
#> GRanges object with 5 ranges and 5 metadata columns:
#>       seqnames          ranges strand | logFCThresh numWindowsThresh direction
#>          <Rle>       <IRanges>  <Rle> |   <numeric>        <integer>  <factor>
#>   [1]     chr1 6940000-6940035      * |     6.88595                2  positive
#>   [2]     chr1 6940468-6940539      * |     6.24705                5  positive
#>   [3]     chr1 6940816-6940863      * |     5.87596                3  positive
#>   [4]     chr1 6940984-6941079      * |     5.76636                7  positive
#>   [5]     chr1 6941152-6941631      * |     5.52749               33  positive
#>           logFC numWindows
#>       <numeric>  <integer>
#>   [1]   6.88595          2
#>   [2]   6.24705          5
#>   [3]   5.87596          3
#>   [4]   5.76636          7
#>   [5]   5.38673         39
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
