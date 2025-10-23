# Filter positions

Filter positions based on any combination of sequence context, coverage,
repetition of the same position, and presence of non-NA values. Filters
are applied in the order they are specified to the `filters` argument.
Any filter type can be repeated an arbitrary number of times.

## Usage

``` r
filterPositions(
  se,
  filters = c("sequenceContext", "coverage", "all.na"),
  sequenceContext = NULL,
  assayNameCov = "Nvalid",
  minCov = 1,
  minNbrSamples = NULL,
  regions = NULL,
  seqinfo = NULL,
  assayNameAmbig = "Nvalid",
  assayNameNA = "mod_prob"
)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- filters:

  A character vector. All values must be one of `"sequenceContext"`,
  `"coverage"`, `"repeated.positions"`, `"regions"` and `"all.na"`.
  Filters are applied in the order specified by this vector.

- sequenceContext:

  A character vector with sequence contexts to retain. To apply this
  filter, the `"sequenceContext"` column must be present in
  `rowData(se)` (see `addSeqContext`).

- assayNameCov:

  A character scalar indicating the assay to use to define the coverage.
  If this is a read-level assay, coverage is first calculated using
  `flattenReadLevelAssay(..., statistics = "Nvalid")`.

- minCov:

  A numeric scalar indicating the lowest acceptable coverage in order to
  keep a position.

- minNbrSamples:

  A numeric scalar, or `NULL`. If `NULL` (default), the row sum of
  `assayNameCov` (i.e., the total coverage across all samples) is used
  for the coverage filtering. If not `NULL`, a position is required to
  have at least `minCov` coverage in at least `minNbrSamples` to be
  retained.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object specifying the genomic regions to restrict the output to.
  Alternatively, regions can be specified as a character vector (e.g.
  "chr1:1200-1300", "chr2:-6000" or "chrM") that will be coerced into a
  `GRanges` object. If the end coordinate is not provided (for example
  in "chr1:10-", which means "to the end of chr1"), it will be obtained
  from `seqinfo` or default to a large value if `seqinfo` is not
  provided.

- seqinfo:

  `NULL` or a
  [`Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html) object
  containing information about the set of genomic sequences
  (chromosomes). Alternatively, a named numeric vector with genomic
  sequence names and lengths.

- assayNameAmbig:

  A character scalar indicating the assay to use to decide which row to
  retain if multiple rows represent the same genomic position (on
  different strands). The row with the largest row sum in this assay is
  retained.

- assayNameNA:

  A character scalar indicating the assay to use as the basis for
  filtering out positions with NA values across all reads. This should
  be a read level assay.

## Value

A filtered `SummarizedExperiment`.

## Author

Charlotte Soneson

## Examples

``` r
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")

se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                 modbase = "a", verbose = FALSE,
                 BPPARAM = BiocParallel::SerialParam())
se <- flattenReadLevelAssay(se)
se <- addSeqContext(se, sequenceContextWidth = 3, sequenceReference = reffile)
sefilt <- filterPositions(se, c("sequenceContext", "coverage", "all.na"),
                          minCov = 5, sequenceContext = "TAG")
```
