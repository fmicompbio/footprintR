# Summarize a read-level object to sample-level

This function will take a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with read-level footprinting data (e.g. returned by
[`readModkitExtract`](https://fmicompbio.github.io/footprintR/reference/readModkitExtract.md)
or
[`readModBam`](https://fmicompbio.github.io/footprintR/reference/readModBam.md))
and summarize reads in each sample, for instance to generate modified
and total counts at each position for each sample.

## Usage

``` r
flattenReadLevelAssay(
  se,
  assayName = "mod_prob",
  statistics = c("Nmod", "Nvalid", "FracMod"),
  modProbThreshold = 0.5,
  keepReads = TRUE,
  replaceExisting = TRUE,
  verbose = FALSE
)
```

## Arguments

- se:

  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object with read-level footprinting data. Rows should correspond to
  positions and columns to samples.

- assayName:

  A character scalar specifying the assay of `se` containing the
  read-level data to be summarized. Typically, this assay contains
  modification probabilities.

- statistics:

  Character vector specifying the type of statistics to be computed.
  Currently supported values are "Nmod" (number of values per row in the
  `assayName` assay that are greater than or equal to
  `modProbThreshold`), "Nvalid" (number of valid/non-NA values per row,
  typically the number of overlapping reads), "FracMod" (Nmod/Nvalid),
  "Pmod" (row-wise average values), "Mean" (equivalent to "Pmod"), "Sum"
  (row-wise sums of non-NA values), "AvgConf" (average confidence of
  (non-)modification probabilities, more precisely the row-wise averages
  of the largest of the observed values and 1 - the observed values).

- modProbThreshold:

  A numeric scalar, indicating the modification probability threshold to
  use to classify a base as 'modified' or 'unmodified'.

- keepReads:

  A logical scalar. If `TRUE` (the default), the read-level data from
  `assayName` will be retained in an assay of the same name.

- replaceExisting:

  A logical scalar. If `TRUE` (the default), any existing assays with
  the same name as the ones requested will be overwritten. Otherwise,
  existing assays will be retained and the corresponding summary
  statistic(s) will not be recalculated.

- verbose:

  If `TRUE`, report on progress.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with the same dimensions as `se` (positions in rows and samples
in columns), and added assays corresponding to the requested statistics.

## See also

[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for the returned object type,
[`readModkitExtract`](https://fmicompbio.github.io/footprintR/reference/readModkitExtract.md)
for the function used to read the input files.

## Author

Charlotte Soneson, Michael Stadler

## Examples

``` r
exfile <- system.file("extdata", "modkit_extract_rc_6mA_1.tsv.gz",
                      package = "footprintR")
se <- readModkitExtract(exfile, modbase = "a",
                        BPPARAM = BiocParallel::SerialParam())
se
#> class: RangedSummarizedExperiment 
#> dim: 8344 1 
#> metadata(3): modkit_threshold filter_threshold readLevelData
#> assays(1): mod_prob
#> rownames(8344): chr1:6925830:- chr1:6925834:- ... chr1:6941530:-
#>   chr1:6941531:-
#> rowData names(0):
#> colnames(1): s1
#> colData names(2): sample modbase

se_summary <- flattenReadLevelAssay(se)
se_summary
#> class: RangedSummarizedExperiment 
#> dim: 8344 1 
#> metadata(3): modkit_threshold filter_threshold readLevelData
#> assays(4): mod_prob Nmod Nvalid FracMod
#> rownames(8344): chr1:6925830:- chr1:6925834:- ... chr1:6941530:-
#>   chr1:6941531:-
#> rowData names(0):
#> colnames(1): s1
#> colData names(2): sample modbase
```
