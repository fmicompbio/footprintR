# Extract data for one or more anchor regions

Extract assay values for each read (or sample) and each position in a
set of anchor regions of the same width. The anchor regions are defined
by their genomic midpoint coordinate and the width.

## Usage

``` r
getAnchorRegions(
  se,
  assayName = "mod_prob",
  regionMidpoints,
  regionWidth,
  anchorName = "anchor",
  prune = TRUE,
  ignore.strand = FALSE,
  reverseMinusStrandRegions = FALSE,
  verbose = FALSE
)
```

## Arguments

- se:

  A `SummarizedExperiment` object.

- assayName:

  Character vector, the assay(s) from which to extract values.

- regionMidpoints:

  Either a `GPos` object or a character vector that can be coerced into
  a `GPos` object, representing the midpoints of the desired anchor
  regions.

- regionWidth:

  Integer scalar, the desired width of the anchor regions. Must be an
  odd value.

- anchorName:

  Character scalar that gives the "sequence name" of the aligned anchor
  regions and will be used in generating the `rowRanges` of the return
  value.

- prune:

  Logical scalar. If `TRUE` (the default), samples for which there are
  no reads overlapping any of the anchor regions in any of the
  read-level assays in `assayName` will be completely removed from the
  returned `SummarizedExperiment` (also from `colData`). If `FALSE`,
  such samples are retained (in assays with read-level data as a
  zero-column `NAMatrix`, in other assays as a dense matrix with a
  single column of NA values).

- ignore.strand:

  Logical scalar, whether to ignore the strand information when matching
  anchor regions with observations. Will be passed on to
  `GenomicRanges::match()`. If `TRUE`, a pruning step will be applied to
  ensure that each genomic position is represented by at most one row in
  the object. If multiple rows are found corresponding to the same
  position (on different strands), the one with the highest number of
  supporting reads (determined by the `Nvalid`, `Nmod`, or `mod_prob`
  assay, in this preference order, depending on what is present in `se`)
  will be retained.

- reverseMinusStrandRegions:

  Logical scalar. If `TRUE`, data extracted from regions on the negative
  strand will be reversed before they are concatenated with the rest of
  the data, so that negative relative positions correspond to regions
  upstream of the region midpoint.

- verbose:

  Logical scalar. If `TRUE`, report on progress.

## Value

A `SummarizedExperiment` with rows representing relative positions
within an anchor region (the midpoint of the region corresponds to a
relative position of 0) and columns representing samples. Each column of
the assay is an `NaArray` (if `assayName` is a read-level assay) or a
dense matrix (otherwise), with columns representing read-anchor region
(or sample-anchor region) combinations. The `region` column of the
`colData` records which anchor region a given column corresponds to. In
addition, all assays will be designated as 'read-level' assays (each
column represents a sample, which is in turn represented by multiple
columns in the actual data, corresponding to the different regions).
This allows downstream analysis similar to that for read-level data,
including flattening and plotting.

## Author

Charlotte Soneson

## Examples

``` r
library(SummarizedExperiment)
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam", "6mA_2_10reads.bam"),
                           package = "footprintR")
se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6920000-6940000",
                 modbase = "a", verbose = FALSE,
                 BPPARAM = BiocParallel::SerialParam())
se <- flattenReadLevelAssay(se)
ar <- getAnchorRegions(se, assayName = c("mod_prob", "FracMod", "Nvalid"),
                       regionMidpoints = c("chr1:6929389:-", "chr1:6935630:-"),
                       regionWidth = 9)

## Modification probabilities
assay(ar)
#> DataFrame with 9 rows and 2 columns
#>                                        s1
#>                                <NaMatrix>
#> 1                            NA:NA:NA:...
#> 2                              0:NA:0:...
#> 3                              0:NA:0:...
#> 4 0.955078125:0.068359375:0.654296875:...
#> 5 0.455078125:0.111328125:0.994140625:...
#> 6                            NA:NA:NA:...
#> 7                            NA:NA:NA:...
#> 8                    NA:0:0.068359375:...
#> 9                            NA:NA:NA:...
#>                                        s2
#>                                <NaMatrix>
#> 1                            NA:NA:NA:...
#> 2          0.990234375:0.103515625:NA:...
#> 3 0.994140625:0.064453125:0.939453125:...
#> 4          0.990234375:0.056640625:NA:...
#> 5                   0.970703125:NA:NA:...
#> 6                            NA:NA:NA:...
#> 7                   0.955078125:NA:NA:...
#> 8          NA:0.166015625:0.423828125:...
#> 9                            NA:NA:NA:...

## Region assignment
colData(ar)
#> DataFrame with 2 rows and 4 columns
#>         sample
#>    <character>
#> s1          s1
#> s2          s2
#>                                                                                                                              region_mod_prob
#>                                                                                                                                       <List>
#> s1 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6929385-6929393..:chr1:6929385-6929393:-
#> s2 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6929385-6929393..:chr1:6929385-6929393:-
#>                                                                                 region_FracMod
#>                                                                                         <List>
#> s1 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6935626-6935634..:chr1:6935626-6935634:-
#> s2 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6935626-6935634..:chr1:6935626-6935634:-
#>                                                                                  region_Nvalid
#>                                                                                         <List>
#> s1 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6935626-6935634..:chr1:6935626-6935634:-
#> s2 chr1:6929385-6929393..:chr1:6929385-6929393:-,chr1:6935626-6935634..:chr1:6935626-6935634:-
```
