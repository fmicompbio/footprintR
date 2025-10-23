# Write bam records from `infile` to `outfile` if they pass filter criteria.

For each bam file in `infiles`, parse alignments, calculate read
statistics and write the alignment to the corresponding output file from
`outfiles` if the read passes all criteria defined by the filtering
arguments. Filters are processed hierarchically: If a read does not pass
a given filter, the remaining filters will not be examined and the
processing continues with the next read. The filters are examined in
this order: `keepUnmapped`, `keepSecondary`, `keepSupplementary`,
`minReadLength`, `minAlignedLength`, `minAlignedFraction`, `minQscore`,
`minSNR`, `maxFracLowConf`, `maxEntropy`.

## Usage

``` r
filterReadsBam(
  infiles,
  outfiles,
  modbase,
  indexOutfiles = TRUE,
  overwriteOutfiles = FALSE,
  keepUnmapped = TRUE,
  keepSecondary = TRUE,
  keepSupplementary = TRUE,
  minReadLength = 0,
  minAlignedLength = 0,
  minAlignedFraction = 0,
  minQscore = 0,
  minSNR = -Inf,
  maxFracLowConf = 1,
  maxEntropy = Inf,
  noiseCoef = c(NA_real_, NA_real_),
  LowConf = 0.7,
  BPPARAM = MulticoreParam(4L, RNGseed = 42L),
  verbose = FALSE
)
```

## Arguments

- infiles:

  Character vector with name(s) of the input bam file(s).

- outfiles:

  Character vector with name(s) of the output bam file(s). Needs to have
  the same length as `infiles`.

- modbase:

  Character scalar defining the modified base to analyze (used by
  `maxEntropy` and `maxFracLowConf`).

- indexOutfiles:

  Logical scalar. If `TRUE` (the default) create a bam index file
  (`.bai` file) for each of the generated `outfiles`.

- overwriteOutfiles:

  Logical scalar. If `FALSE` (the default), existing `outfiles` will not
  be overwritten and the function will abort with an error message.

- keepUnmapped, keepSecondary, keepSupplementary:

  Logical scalars indicating whether to keep unmapped, secondary or
  supplementary alignments.

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

- minQscore:

  A numeric scalar representing the smallest acceptable read-level
  Qscore. Reads with Qscore below this value will be filtered out.

- minSNR:

  Numeric scalar. Minimum acceptable read Signal-to-Noise Ratio (SNR).
  Reads with SNR below this value are filtered out. Set to `-Inf` to
  disable SNR filtering. The SNR is computed per read as
  \\\log_2(\mathrm{SignalVar}/\mathrm{NoiseVar})\\ where
  \\\mathrm{NoiseVar} \approx 0.5\\\mathrm{Var}(\Delta x)\\ using
  adjacent methylation differences that can jump up to \\k=2\\ gaps, and
  \\\mathrm{SignalVar} = \max(\mathrm{Var}(x) - \mathrm{NoiseVar},
  \varepsilon)\\ with \\\varepsilon=10^{-3}\\.

- maxFracLowConf:

  A numeric scalar representing the maximally acceptable fraction of
  low-confidence modified base calls in a read. Reads without
  modified-base calls or with a fraction of low confidence calls greater
  than this value will be filtered out.

- maxEntropy:

  A numeric scalar representing the largest acceptable read-level
  entropy. Reads without modified-base calls or with entropy above this
  value will be filtered out. A value of `Inf` deactivates the entropy
  filter.

- noiseCoef:

  Numeric vector of length 2 giving the background noise model
  coefficients `(b0, b1)`. If provided, these are used to impose a floor
  on the estimated per-read noise variance (see also
  [`calcReadStats`](https://fmicompbio.github.io/footprintR/reference/calcReadStats.md)).

- LowConf:

  A numeric scalar with the minimum call confidence below which calls
  are considered "low confidence".

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object that controls the number of parallel CPU threads to use for
  some of the steps in `filterReadsBam()`. The default value is
  ([`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)`(4L, RNGseed = 42L)`).

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

`filterReadsBam` is called for its side effect of generating new bam
files containing the subset of bam records from input bam files that
pass all filtering criteria. In addition, it returns a `data.frame` with
one row per `infiles` giving the numbers of bam records that were read
in `total`, that were `retained` in the `outfiles` and that were
filtered-out by reason of exclusion.

## Author

Michael Stadler

## Examples

``` r
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
filtbamfiles <- tempfile(fileext = rep(".bam", length(modbamfiles)))
res <- filterReadsBam(infiles = modbamfiles, outfiles = filtbamfiles,
                      modbase = "a", indexOutfiles = FALSE, minReadLength = 6746,
                      minAlignedLength = 6896, minAlignedFraction = 0.56, minSNR=-0.768,
                      minQscore = 9.7, maxFracLowConf = 0.11, maxEntropy = 0.29,
                      BPPARAM = BiocParallel::SerialParam(), verbose = TRUE)
#> ℹ start filtering of /Users/runner/work/_temp/Library/footprintR/extdata/6mA_1_10reads.bam using 1 thread
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ merging 1 filtered chunks
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ done filtering: retained 6 of 10 records (60%)
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ start filtering of /Users/runner/work/_temp/Library/footprintR/extdata/6mA_2_10reads.bam using 1 thread
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ merging 1 filtered chunks
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
#> ℹ done filtering: retained 7 of 10 records (70%)
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [2ms]
res
#>   sample                                                                infile
#> 1     s1 /Users/runner/work/_temp/Library/footprintR/extdata/6mA_1_10reads.bam
#> 2     s2 /Users/runner/work/_temp/Library/footprintR/extdata/6mA_2_10reads.bam
#>                                                                             outfile
#> 1 /var/folders/q0/wmf37v850txck86cpnvwm_zw0000gn/T//Rtmp0GBwua/file36ae207e2fec.bam
#> 2 /var/folders/q0/wmf37v850txck86cpnvwm_zw0000gn/T//Rtmp0GBwua/file36ae3c7c3f9e.bam
#>   total retained filtered_unmapped filtered_secondary filtered_supplementary
#> 1    10        6                 0                  0                      0
#> 2    10        7                 0                  0                      0
#>   filtered_minReadLength filtered_minAlignedLength filtered_minAlignedFraction
#> 1                      0                         1                           1
#> 2                      1                         0                           0
#>   filtered_minQscore filtered_minSNR filtered_maxFracLowConf
#> 1                  0               1                       0
#> 2                  1               0                       1
#>   filtered_maxEntropy
#> 1                   1
#> 2                   0
unlink(filtbamfiles)
```
