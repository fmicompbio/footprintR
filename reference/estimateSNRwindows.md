# Quantify footprint coherence in windows by estimating the Signal to Noise Ratio (SNR)

Calculate a measure of footprint coherence and window "structureness" by
decomposing total aggregate signal variance into (i) `noise` and (ii)
`signal` variance components and returning
\\\log_2(\text{signal}/\text{noise})\\ as the SNR.

## Usage

``` r
estimateSNRwindows(
  se,
  gr,
  assayNameAgg = "FracMod",
  k = 2L,
  min_diffs = NULL,
  noise_pars = NULL
)
```

## Arguments

- se:

  A per-position `RangedSummarizedExperiment`, typically the output of
  `readModBam(..., level = "summary")`. Must contain the assay given by
  `assayNameAgg` and an assay `"Nvalid"` with coverage.

- gr:

  `GRanges` specifying the windows.

- assayNameAgg:

  Character scalar, assay to aggregate (default `"FracMod"`).

- k:

  Integer (default 2L). Two measurements are considered “adjacent” if
  they are at most `k` bases apart when estimating the noise.

- min_diffs:

  Integer or `NULL`. Minimal number of paired differences to compute a
  reliable noise estimate. If `NULL` the value is set to
  `max(16, floor(0.05 * n))`, where `n` is the number of valid positions
  in the window.

- noise_pars:

  Numeric vector as returned by the `coefficients` slot of the
  [`estimateNoiseParsWindows()`](https://fmicompbio.github.io/footprintR/reference/estimateNoiseParsWindows.md)
  result. If `NULL` the function looks for `metadata(se)$NoisePars`. If
  still `NULL` raw noise estimates are used.

## Value

A `SummarizedExperiment` with assays `"TotalVar"`, `"SignalVar"`,
`"NoiseVar"` and `"SNR"` (\\\log_2\\ signal/noise).

## Details

Noise variance estimation can be performed in two ways:

1.  **Parametric estimation**: If noise parameters (`noise_pars`) are
    supplied (typically estimated using
    [`estimateNoiseParsWindows()`](https://fmicompbio.github.io/footprintR/reference/estimateNoiseParsWindows.md),
    which fits a linear regression model to relate noise variance to
    average methylation level and coverage across random genomic
    windows), the noise variance returned is the expected noise given
    the window's average methylation and coverage.

2.  **Non-parametric (raw) estimation**: If noise parameters are not
    supplied, the function estimates noise directly from the data, based
    on short-range variation between neighboring positions. While
    unbiased for random signals, this method can *overestimate*
    noise—and thus underestimate the true SNR—if the biological signal
    is strongly correlated across reads (e.g., due to systematic
    accessibility patterns), because the structured variation is
    incorrectly interpreted as random noise.

## See also

[`estimateNoiseParsWindows`](https://fmicompbio.github.io/footprintR/reference/estimateNoiseParsWindows.md)
for estimation of the global noise model parameters.

## Author

Panagiotis Papasaikas

## Examples

``` r
library(GenomicRanges)
modbamfiles <- system.file("extdata",
                           c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
gr0 <- GRanges("chr1", IRanges(6920000, 6950000))
gr_tiles <- unlist(tile(gr0, width = 500))
gr_tiles
#> GRanges object with 61 ranges and 0 metadata columns:
#>        seqnames          ranges strand
#>           <Rle>       <IRanges>  <Rle>
#>    [1]     chr1 6920000-6920490      *
#>    [2]     chr1 6920491-6920982      *
#>    [3]     chr1 6920983-6921474      *
#>    [4]     chr1 6921475-6921966      *
#>    [5]     chr1 6921967-6922458      *
#>    ...      ...             ...    ...
#>   [57]     chr1 6947541-6948032      *
#>   [58]     chr1 6948033-6948524      *
#>   [59]     chr1 6948525-6949016      *
#>   [60]     chr1 6949017-6949508      *
#>   [61]     chr1 6949509-6950000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# estimate noise parameters (typically should be performed on a larger,
# independent sample of genomic windows)
resNoiseEst <- estimateNoiseParsWindows(modbamfiles, modbase = "a",
                                        windows = gr_tiles,
                                        return_data = FALSE,
                                        BPPARAM = BiocParallel::SerialParam())

se <- readModBam(bamfiles = modbamfiles,
                 regions = gr_tiles,
                 modbase = "a",
                 level = "summary",
                 trim = TRUE,
                 BPPARAM = BiocParallel::SerialParam()) |>
   filterPositions(filters = "coverage",
                   minCov=2, assayNameNA = NULL)
# parametric noise estimation
se_snr_scores1 <- estimateSNRwindows(se, gr = gr_tiles,
                                     noise_pars = resNoiseEst$coefficients)
# non-parametric noise estimation
se_snr_scores2 <- estimateSNRwindows(se, gr = gr_tiles)
```
