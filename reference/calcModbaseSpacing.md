# Calculate distances between modified bases on individual reads

Calculate the frequencies of same-read modified base distances, for
example from read-level modification data to estimate nucleosome repeat
length. Distances are calculated separately for each sample (column in
`se`), but if needed they can be easily combined for estimating NRL on a
pool of samples by summing up observed counts (e.g. using
`Reduce("+", sampleDistsList)`. Distance calculations are implemented in
C++
([`calcAndCountDist`](https://fmicompbio.github.io/footprintR/reference/calcAndCountDist.md))
for efficiency.

## Usage

``` r
calcModbaseSpacing(
  se,
  assayName = "mod_prob",
  minModProb = 0.5,
  poolReads = TRUE,
  dmax = 1000L
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object with read-level footprinting data, for example returned by
  [`readModBam`](https://fmicompbio.github.io/footprintR/reference/readModBam.md).
  Rows should correspond to positions and columns to samples.

- assayName:

  A character scalar specifying the assay of `se` containing the
  read-level modification probabilities.

- minModProb:

  Numeric scalar giving the minimal modification probability for a
  modified base.

- poolReads:

  Logical scalar indicating if reads within a sample should be pooled.
  If `TRUE` (the default), distances from reads within a sample are
  combined and returned as a vector. If `FALSE`, distances obtained from
  each read in a sample are returned separately as columns in a matrix.

- dmax:

  Numeric scalar specifying the maximal distance between modified bases
  on the same read to count.

## Value

A named list of length `ncol(se)` (one element for each sample). If
`poolReads=TRUE`, the elements are `integer` vectors of length `dmax`,
with the value at position `d` giving the observed number of within-read
modified base pairs at distance `d`. If `poolReads=FALSE`, each list
element is a matrix with `dmax` rows and individual reads in columns,
with the value at row `d` and column `r` giving the observed number of
modified base pairs at distance `d` for read `r`.

## References

Phasograms were originally described in Valouev et al., Nature 2011
(doi:10.1038/nature10002). The implementation here differs in three ways
from the original algorithms:

1.  Instead of same strand alignment start positions, this function is
    adapted to single-molecule footprinting data and measures the
    distances between same-read modified base positions.

2.  It does not implement removing of positions that have been seen less
    than `n` times (referred to as a `n`-pile subset in the paper).

3.  It does allow to retain only alignments that fall into selected
    genomic intervals (`regions` argument).

## See also

[`estimateNRL`](https://fmicompbio.github.io/footprintR/reference/estimateNRL.md)
to estimate the nucleosome repeat length from a phasogram,
[`plotModbaseSpacing`](https://fmicompbio.github.io/footprintR/reference/plotModbaseSpacing.md)
to visualize an annotated phasogram,
[`calcAndCountDist`](https://fmicompbio.github.io/footprintR/reference/calcAndCountDist.md)
for low-level distance counting.

## Author

Michael Stadler

## Examples

``` r
modbamfiles <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
se <- readModBam(modbamfiles, "chr1:6940000-6955000", "a",
                 BPPARAM = BiocParallel::SerialParam())

# get distances
moddist <- calcModbaseSpacing(se)
str(moddist)
#> List of 1
#>  $ s1: Named num [1:1000] 175 110 103 70 83 99 111 93 77 68 ...
#>   ..- attr(*, "names")= chr [1:1000] "1" "2" "3" "4" ...
```
