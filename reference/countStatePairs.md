# Count pairs of modified bases by distance and modification state

For all pairs of bases with modification calls in a read, tabulate the
number of pairs at a given distance with a given modification state.

## Usage

``` r
countStatePairs(
  bamfile,
  regions = ".",
  modbase,
  threshUnmod = 0.5,
  threshMod = 0.5,
  windowSize = 200,
  minMapQ = 0,
  minAlignedLength = 0,
  BPPARAM = MulticoreParam(4L),
  verbose = FALSE
)
```

## Arguments

- bamfile:

  Character scalar with the path to a `modBAM` file, containing
  information about base modifications in `MM` and `ML` tags. The
  `bamfile` must have an index.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object specifying which genomic regions to extract the reads from.
  Alternatively, regions can be specified as a character vector (e.g.
  "chr1:1200-1300", "chr2:-6000", "chr1:10-", "chrM" or ".").

- modbase:

  Character scalar defining the modified base.

- threshUnmod, threshMod:

  Numeric scalars defining how to convert modification probabilities `p`
  to modification states. Bases with `p < threshUnmod` will be
  considered unmodified, and bases with `p >= threshMod` modified. All
  remaining bases will be ignored.

- windowSize:

  Numeric scalar giving the maximum window size covering pairs of
  modified bases to consider in pair-counting mode. A window size of 1
  corresponds to a single base, a size of 2 to directly adjacent bases,
  etc.

- minMapQ:

  Numeric scalar giving the minimal mapping quality to include
  alignments in pair-counting mode.

- minAlignedLength:

  Numeric scalar giving the minimal alignment length to include
  alignments in pair-counting mode.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object that controls the number of parallel CPU threads to use for
  decompressing bam records.

- verbose:

  A logical scalar. If `TRUE`, report on progress.

## Value

A `DataFrame` with `windowSize` rows and five columns.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "footprintR")
res <- countStatePairs(bamfile = modbamfile,
                       regions = "chr1",
                       modbase = "a", windowSize = 300,
                       BPPARAM = BiocParallel::SerialParam())
res
#> DataFrame with 300 rows and 5 columns
#>             S unmod_unmod unmod_mod mod_unmod   mod_mod
#>     <integer>   <numeric> <numeric> <numeric> <numeric>
#> 1           1       27084         0         0      2461
#> 2           2        8877       281       340       408
#> 3           3        9131       441       439       264
#> 4           4        8140       399       403       282
#> 5           5        8565       485       493       164
#> ...       ...         ...       ...       ...       ...
#> 296       296        7399       701       615        63
#> 297       297        7097       654       584        60
#> 298       298        7342       680       650        66
#> 299       299        7222       673       654        70
#> 300       300        7196       642       640        65
```
