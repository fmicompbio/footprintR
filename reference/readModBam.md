# Read base modifications from bam file(s)

Parse ML and MM tags (see
https://samtools.github.io/hts-specs/SAMtags.pdf, section 1.7) and
return a
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with information on modified bases. Implicitly called bases will
get a modification probability of zero. Secondary and supplementary
alignments are ignored.

## Usage

``` r
readModBam(
  bamfiles,
  regions = NULL,
  modbase,
  level = ifelse(nAlnsToSample == 0 & is.null(variantPositions), "quickread", "read"),
  sampleAnnot = NULL,
  nAlnsToSample = 0,
  seqnamesToSampleFrom = "chr19",
  seqinfo = NULL,
  sequenceContextWidth = 0,
  sequenceReference = NULL,
  variantPositions = NULL,
  modProbThreshold = 0.5,
  trim = FALSE,
  BPPARAM = MulticoreParam(4L, RNGseed = 42L),
  verbose = FALSE
)
```

## Arguments

- bamfiles:

  Character vector with one or several paths of `modBAM` files,
  containing information about base modifications in `MM` and `ML` tags.
  If `bamfiles` is a named vector, the names are used as sample names
  and prefixes for read names. Otherwise, the prefixes will be `s1`,
  ..., `sN`, where `N` is the length of `bamfiles`. All `bamfiles` must
  have an index.

- regions:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object specifying which genomic regions to extract the reads from.
  Alternatively, regions can be specified as a character vector (e.g.
  "chr1:1200-1300", "chr2:-6000" or "chrM") that will be coerced into a
  `GRanges` object. If the end coordinate is not provided (for example
  in "chr1:10-", which means "to the end of chr1"), it will be obtained
  from `seqinfo` or default to a large value if `seqinfo` is not
  provided. Note that unless `trim=TRUE`, the reads are not trimmed to
  the boundaries of the specified ranges. As a result, returned
  positions will typically extend out of the specified regions. If
  `nAlnsToSample` is set to a non-zero value, `regions` is ignored.

- modbase:

  Character vector defining the modified base for each sample. If
  `modbase` is a named vector, the names should correspond to the names
  of `bamfiles`. Otherwise, it will be assumed that the elements are in
  the same order as the files in `bamfiles`. If `modbase` has length 1,
  the same modified base will be used for all samples.

- level:

  Character scalar specifying the level of returned modification data.
  Supported values are:

  "read"

  :   : Extracts modification probabilities for individual reads into an
      assay called `"mod_prob"`, with each column (sample) consisting of
      a position-by-read
      [`NaMatrix`](https://rdrr.io/pkg/SparseArray/man/NaArray-class.html).
      This is the default if `nAlnsToSample` is non-zero or
      `variantPositions` is not `NULL`.

  "quickread"

  :   : Like "read", but runs faster, does not support read sampling,
      and does not return all annotations currently provided by "read".
      This is the default if `nAlnsToSample` is zero and
      `variantPositions` is `NULL`.

  "summary"

  :   : Counts the total and modified bases for each position and strand
      and returns them in assays named `"Nvalid"` and `"Nmod"`,
      respectively, as well as the `Nmod/Nvalid` ratio in the assay
      `"FracMod"`.

- sampleAnnot:

  A `data.frame` (or `NULL`) providing annotations for the samples. It
  must contain at least one column, named `"sample"`, which must contain
  all the values of `names(bamfiles)`. The provided annotations will be
  propagated to the returned `SummarizedExperiment` object.

- nAlnsToSample:

  A numeric scalar. If non-zero, `regions` is ignored and approximately
  `nAlnsToSample` randomly selected alignments on `seqnamesToSampleFrom`
  are read from each of the `bamfiles`. In order to make the results
  reproducible, make sure to set the `RNGseed` argument in the provided
  `BPPARAM` object (see below). Please note that secondary and
  supplementary alignments in `bamfiles` contribute to the total number
  of alignments but will not be sampled, thus the number of returned
  alignments may be lower than `nAlnsToSample`.

- seqnamesToSampleFrom:

  A character vector with one or several sequence names (chromosomes)
  from which to sample alignments from (only used if `nAlnsToSample` is
  greater than zero).

- seqinfo:

  `NULL` or a
  [`Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html) object
  containing information about the set of genomic sequences
  (chromosomes). Alternatively, a named numeric vector with genomic
  sequence names and lengths. Useful to set the sorting order of
  sequence names.

- sequenceContextWidth, sequenceReference:

  Define the sequence context to be extracted around modified bases. By
  default ( `sequenceContextWidth = 0`), no sequence context will be
  extracted, otherwise it will be returned in
  `rowData(x)$sequenceContext`. See
  [`addSeqContext`](https://fmicompbio.github.io/footprintR/reference/addSeqContext.md)
  for details.

- variantPositions:

  An optional `GPos` object with seqnames and coordinates of single
  nucleotide variant positions, to be used to construct read labels for
  allele-specific analysis. Ignored if `NULL` or `nAlnsToSample > 0`
  (sampling-mode).

- modProbThreshold:

  A numeric scalar, indicating the modification probability threshold to
  use to classify a base as 'modified' or 'unmodified'. Only used if
  `level` is `"summary"`.

- trim:

  A logical scalar. If `TRUE`, the returned `SummarizedExperiment`
  object will only contain the positions overlapping the specified
  `regions`. If `FALSE` (default), the object will be extended to all
  positions covered by the reads overlapping `regions`. In both cases,
  only reads overlapping the specified `regions` are included.

- BPPARAM:

  A
  [`BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  object that controls the number of parallel CPU threads to use for
  some of the steps in `readModBam()`. The default value is
  ([`MulticoreParam`](https://rdrr.io/pkg/BiocParallel/man/MulticoreParam-class.html)`(4L, RNGseed = 42L)`).
  If randomly sampling reads (`nAlnsToSample > 0`), make sure to set the
  `RNGseed` argument when constructing the `BPPARAM` object for
  reproducible results (see also
  [`vignette("Random_Numbers", package = "BiocParallel")`](https://bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Random_Numbers.html)).

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with genomic positions in rows and samples in columns. The assays
depend on the value of the `level` argument (see above).

## See also

https://samtools.github.io/hts-specs/SAMtags.pdf describing the SAM ML
and MM tags for base modifications.

## Author

Michael Stadler, Charlotte Soneson

## Examples

``` r
modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
                          package = "footprintR")
readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
           modbase = "a", verbose = TRUE,
           BPPARAM = BiocParallel::SerialParam())
#> ℹ extracting base modifications from modBAM files
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [51ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ collapsed 11300 positions to 4772 unique ones
#> ✔ collapsed 11300 positions to 4772 unique ones [388ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> class: RangedSummarizedExperiment 
#> dim: 4772 1 
#> metadata(3): readLevelData variantPositions filteredOutReads
#> assays(1): mod_prob
#> rownames(4772): chr1:6925830:- chr1:6925834:- ... chr1:6941622:-
#>   chr1:6941631:-
#> rowData names(0):
#> colnames(1): s1
#> colData names(4): sample modbase n_reads readInfo
```
