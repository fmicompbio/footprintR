# Calculate or add summary statistics for read-level base modification data

`calcReadStats` calculates various per-read summary statistics on
modification probabilities or calls from a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with genomic positions in rows and samples in columns.
`addReadStats` adds them to the `colData` under `name`. See details for
more information on the statistics that are calculated.

## Usage

``` r
calcReadStats(
  se,
  assayName = "mod_prob",
  stats = NULL,
  regions = NULL,
  sequenceContext = NULL,
  minNobsPpos = 0,
  minNobsPread = 0,
  LowConf = 0.7,
  LagRange = c(12, 64),
  EntrControl = NULL,
  BPPARAM = MulticoreParam(4L, RNGseed = 42L),
  verbose = FALSE
)

addReadStats(se, ..., name = "QC")
```

## Arguments

- se:

  A
  [`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object with assay `assayName` typically returned by
  [`readModkitExtract`](https://fmicompbio.github.io/footprintR/reference/readModkitExtract.md)
  or
  [`readModBam`](https://fmicompbio.github.io/footprintR/reference/readModBam.md).

- assayName:

  A character scalar specifying the assay of `se` containing the
  read-level data to be summarized. Typically, this assay contains
  modification probabilities.

- stats:

  Character vector specifying which statistics to calculate. When set to
  `NULL` all available statistics except the sample entropy are
  calculated. See details for available read statistics.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object limiting the positions included in the calculations to the ones
  overlapping the corresponding genomic regions. Alternatively, regions
  can be specified as a character vector (e.g. "chr1:1200-1300") that
  can be coerced into a `GRanges` object.

- sequenceContext:

  A character vector with sequence context(s) to include in the
  calculations. Only positions that match one of the provided sequence
  contexts will be included. Sequence contexts can be provided using
  IUPAC redundancy codes. The sequence contexts of modified bases are
  obtained from `rowData(se)$sequenceContext` and thus requires that
  `se` contains the appropriate information, for example by setting the
  `sequenceContextWidth` and `sequenceReference` arguments of
  [`readModkitExtract`](https://fmicompbio.github.io/footprintR/reference/readModkitExtract.md)
  when it was generated, or by adding it using
  [`addSeqContext`](https://fmicompbio.github.io/footprintR/reference/addSeqContext.md).

- minNobsPpos:

  A numeric scalar value \>=1 indicating the minimum coverage on
  individual positions for them to be included in the calculations. In
  high coverage data this is an effective filter for removing spurious
  modbases, typically the result of erroneous basecalling.

- minNobsPread:

  A numeric scalar with the minimum number of observed modifiable bases
  per read for it to be included in the calculations. `NA` values are
  returned for the reads that do not pass this threshold.

- LowConf:

  A numeric scalar with the minimum call confidence below which calls
  are considered "low confidence".

- LagRange:

  A numeric vector of two values (minimum and maximum) defining the
  range of lags for the calculation of autocorrelation and partial
  autocorrelation (see details section).

- EntrControl:

  Optional named list with elements `m`, `r`, `maxStarts`, `nThreads` to
  control Sample Entropy calculation. Missing elements use defaults:
  m=2L, r=0.2, maxStarts=1000, nThreads=1L. See also
  [`sampleEntropy`](https://fmicompbio.github.io/footprintR/reference/sampleEntropy.md)

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object that controls the number of parallel CPU threads to use for
  calculating read statistics. The default value is
  ([`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)`(4L, RNGseed = 42L)`).

- verbose:

  If `TRUE`, report on progress.

- ...:

  For `addReadStats` only: Additional arguments passed on to
  `calcReadStats`.

- name:

  For `addReadStats` only: A character scalar specifying the name to be
  used to store the result in the
  [`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  of the output.

## Value

For `calcReadStats`, a `SimpleList` object with summary statistics for
the samples (columns) in `se`.

For `addReadStats`, a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with summary statistics added to the `name` column of
[`colData`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

## Details

`calcReadStats` calculates a collection of location/scatter statistics
and information theoretic/signal-processing metrics for the modification
probability, confidence or modification call value vectors across
individual reads (data in assay `assayName`). Only bases matching the
criteria given by`regions`, `sequenceContext`, `minNobsPpos` and
`minNobsPread` are included in the calculations. The values of these
filtering parameters are stored in the attribute of the output.

`stats` selects the summaries to be calculated. Currently available
values are:

- MeanModProb:

  : Mean modification probability across the read.

- FracMod:

  : Fraction of confidently called bases (either modified or
  unmodified), that are modified. By default, all bases are considered
  confidently called.

- MeanConf:

  : Mean call confidence across the read.

- MeanConfUnm:

  : Mean call confidence confined to unmodified bases (modifiable bases
  with modification probability \< 0.5).

- MeanConfMod:

  : Mean call confidence confined to modified bases (modifiable bases
  with modification probability \>= 0.5).

- FracLowConf:

  : Fraction of modifiable bases called with low confidence (call
  confidence \< `LowConf`).

- IQRModProb:

  : Interquartile range of modification probabilities across the read.

- sdModProb:

  : Standard deviation of modification probabilities across the read.

- SEntrModProb:

  : Sample entropy of the modification probability signal. Sample
  entropy is a metric assessing the complexity of one-dimensional
  physiological signals. The higher the sample entropy, the more
  irregular, unpredictable and therefore complex the signal. See
  [wikipedia:Sample_entropy](https://en.wikipedia.org/wiki/Sample_entropy)
  for more details.

- ACModProb:

  : Autocorrelation of the modification probability values for lags in
  the range `LagRange`. This range typically covers the signal of
  nucleosome periodicity.

- PACModProb:

  : Partial autocorrelation of the modification probability values for
  lags in the range `LagRange`. This range typically covers the signal
  of nucleosome periodicity.

- NoiseVar:

  : raw **Read Noise variance** estimated as \\0.5\\\mathrm{Var}(\Delta
  x)\\, where \\\Delta x\\ are lag-1 differences that may skip up to
  \\k\\ consecutive NAs. A floor is applied: \\\mathrm{noise} =
  \max(\mathrm{rawNoise},\\ b_0 + b_1 \bar{x})\\.

- SignalVar:

  : **Read Signal variance**
  \\\max(\mathrm{totalVar}-\mathrm{NoiseVar},\\\varepsilon)\\.

- SNR:

  : **Read Signal-to-Noise Ratio**
  \\\log_2(\mathrm{SignalVar}/\mathrm{NoiseVar})\\.

When SNR-related statistics are requested, a robust noise floor
`b0 + b1 * mean(x)` is fitted per sample (see *NoiseVar*). The fitted
coefficients are stored in the result metadata (see below). If
SNR-related statistics are not computed, `NA`/`NA` placeholders are
stored for uniformity.

## Author

Panagiotis Papapasaikas, Charlotte Soneson, Michael Stadler

## Examples

``` r
# load example data
library(SummarizedExperiment)
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "footprintR")
se <- readModBam(bamfile = modbamfile, regions = "chr1:6940000-6955000",
           modbase = "a", verbose = TRUE,
           BPPARAM = BiocParallel::SerialParam())
#> ℹ extracting base modifications from modBAM files
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [58ms]
#> 
#> ℹ collapsed 11300 positions to 4772 unique ones
#> ✔ collapsed 11300 positions to 4772 unique ones [315ms]
#> 

readStats <- calcReadStats(se, BPPARAM = BiocParallel::SerialParam())
#> Warning: Too few points to estimate noise floor (1); raw noise variances are used.
readStats$s1
#> DataFrame with 3 rows and 13 columns
#>                                         MeanModProb   FracMod  MeanConf
#>                                           <numeric> <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563   0.1652868 0.1498969  0.943914
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9   0.0835376 0.0646707  0.958320
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d   0.0556013 0.0347512  0.965854
#>                                         MeanConfUnm MeanConfMod FracLowConf
#>                                           <numeric>   <numeric>   <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563    0.957961    0.864255   0.0680724
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9    0.967633    0.823622   0.0470060
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d    0.971512    0.808703   0.0355852
#>                                         IQRModProb sdModProb
#>                                          <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563  0.1308594  0.312860
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9  0.0488281  0.215142
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d  0.0000000  0.164023
#>                                                                       ACModProb
#>                                                                          <list>
#> s1-233e48a7-f379-4dcf-9270-958231125563       0.1306906,0.0840176,0.0591518,...
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9       0.0619512,0.0479904,0.0392903,...
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d -0.00193260,-0.01562365,-0.00424019,...
#>                                                                      PACModProb
#>                                                                          <list>
#> s1-233e48a7-f379-4dcf-9270-958231125563    -0.0415954,-0.0216324,-0.0132802,...
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9 -0.00443190,-0.00209608,-0.00710440,...
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d -0.01637500, 0.00794205,-0.02053895,...
#>                                               SNR SignalVar  NoiseVar
#>                                         <numeric> <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563  0.698999 0.0605703 0.0373113
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9  0.154130 0.0243781 0.0219079
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d  0.307366 0.0148793 0.0120242

se_withReadStats <- addReadStats(se, name = "QC",
                                 BPPARAM = BiocParallel::SerialParam())
#> Warning: Too few points to estimate noise floor (1); raw noise variances are used.
se_withReadStats$QC$s1
#> DataFrame with 3 rows and 13 columns
#>                                         MeanModProb   FracMod  MeanConf
#>                                           <numeric> <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563   0.1652868 0.1498969  0.943914
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9   0.0835376 0.0646707  0.958320
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d   0.0556013 0.0347512  0.965854
#>                                         MeanConfUnm MeanConfMod FracLowConf
#>                                           <numeric>   <numeric>   <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563    0.957961    0.864255   0.0680724
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9    0.967633    0.823622   0.0470060
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d    0.971512    0.808703   0.0355852
#>                                         IQRModProb sdModProb
#>                                          <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563  0.1308594  0.312860
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9  0.0488281  0.215142
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d  0.0000000  0.164023
#>                                                                       ACModProb
#>                                                                          <list>
#> s1-233e48a7-f379-4dcf-9270-958231125563       0.1306906,0.0840176,0.0591518,...
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9       0.0619512,0.0479904,0.0392903,...
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d -0.00193260,-0.01562365,-0.00424019,...
#>                                                                      PACModProb
#>                                                                          <list>
#> s1-233e48a7-f379-4dcf-9270-958231125563    -0.0415954,-0.0216324,-0.0132802,...
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9 -0.00443190,-0.00209608,-0.00710440,...
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d -0.01637500, 0.00794205,-0.02053895,...
#>                                               SNR SignalVar  NoiseVar
#>                                         <numeric> <numeric> <numeric>
#> s1-233e48a7-f379-4dcf-9270-958231125563  0.698999 0.0605703 0.0373113
#> s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9  0.154130 0.0243781 0.0219079
#> s1-92e906ae-cddb-4347-a114-bf9137761a8d  0.307366 0.0148793 0.0120242
metadata(se_withReadStats$QC$s1)
#> list()
```
