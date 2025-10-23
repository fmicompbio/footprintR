# Regroup reads

Regroup reads into new "samples", defined by `readGroups`.

## Usage

``` r
regroupReads(se, readGroups)

regroupReadsByColData(se, colNames, withinSample = FALSE)
```

## Arguments

- se:

  A
  [`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  object, e.g. obtained by
  [`readModBam`](https://fmicompbio.github.io/footprintR/reference/readModBam.md).
  Rows represent positions, and columns samples. The object should have
  at least one read-level assay.

- readGroups:

  A named list, where each element corresponds to a new group of reads
  (a "sample").

- colNames:

  A character vector corresponding to the names of columns in
  `colData(se)$readInfo`, the combination of which represent the desired
  grouping of the reads.

- withinSample:

  A logical scalar, indicating whether the regrouping should be done
  within each current sample (column of `se`) or not. If `FALSE`
  (default), reads are pooled across samples before being regrouped.

## Value

A
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object with columns corresponding to the names of `readGroups`. Only the
read-level assays will be returned (as summary-level assays are likely
no longer consistent with the regrouped read-level assays).

## Author

Charlotte Soneson

## Examples

``` r
library(SummarizedExperiment)
library(GenomicRanges)
modbamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                        "6mA_2_10reads.bam"),
                           package = "footprintR")
se <- readModBam(bamfiles = modbamfiles, regions = "chr1:6940000-6955000",
                 modbase = "a", verbose = TRUE,
                 variantPositions = GPos(seqnames = "chr1",
                                         pos = c(6940000, 6940500)),
                 BPPARAM = BiocParallel::SerialParam())
#> ℹ extracting base modifications from modBAM files
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ opening input file /Users/runner/work/_temp/Library/footprintR/extdata/6mA_1_10reads.bam using 1 thread
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ reading alignments overlapping 1 region
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ removed 1724 unaligned (e.g. soft-masked) of 25090 called bases
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ read 3 alignments
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ opening input file /Users/runner/work/_temp/Library/footprintR/extdata/6mA_2_10reads.bam using 1 thread
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ reading alignments overlapping 1 region
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ removed 320 unaligned (e.g. soft-masked) of 10174 called bases
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ read 2 alignments
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ finding unique genomic positions...
#> ✔ finding unique genomic positions... [32ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
#> ℹ collapsed 17739 positions to 7967 unique ones
#> ✔ collapsed 17739 positions to 7967 unique ones [256ms]
#> 
#> ⠙ 0.000 Mio. genomic positions processed (0.001 Mio./s) [1ms]
# number of reads per sample
lapply(assay(se, "mod_prob"), ncol)
#> $s1
#> [1] 3
#> 
#> $s2
#> [1] 2
#> 
# define groups
groups <- list(g1 = c("s1-233e48a7-f379-4dcf-9270-958231125563",
                      "s2-d03efe3b-a45b-430b-9cb6-7e5882e4faf8"),
               g2 = "s1-92e906ae-cddb-4347-a114-bf9137761a8d",
               g3 = c("s2-034b625e-6230-4f8d-a713-3a32cd96c298",
                      "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9"))
# regroup reads
sere <- regroupReads(se, readGroups = groups)
lapply(assay(sere, "mod_prob"), ncol)
#> $g1
#> [1] 2
#> 
#> $g2
#> [1] 1
#> 
#> $g3
#> [1] 2
#> 

# regroup by read annotation (variant label)
# ... across samples
sere <- regroupReadsByColData(se, colNames = "variant_label",
                              withinSample = FALSE)
sere
#> class: RangedSummarizedExperiment 
#> dim: 7967 2 
#> metadata(3): readLevelData variantPositions filteredOutReads
#> assays(1): mod_prob
#> rownames(7967): chr1:6925830:- chr1:6925834:- ... chr1:6941622:-
#>   chr1:6941631:-
#> rowData names(0):
#> colnames(2): G- GT
#> colData names(4): sample modbase n_reads readInfo
lapply(assay(sere, "mod_prob"), ncol)
#> $`G-`
#> [1] 3
#> 
#> $GT
#> [1] 2
#> 
# ... within sample
sere <- regroupReadsByColData(se, colNames = "variant_label",
                              withinSample = TRUE)
sere
#> class: RangedSummarizedExperiment 
#> dim: 7967 3 
#> metadata(3): readLevelData variantPositions filteredOutReads
#> assays(1): mod_prob
#> rownames(7967): chr1:6925830:- chr1:6925834:- ... chr1:6941622:-
#>   chr1:6941631:-
#> rowData names(0):
#> colnames(3): s1-G- s1-GT s2-G-
#> colData names(4): sample modbase n_reads readInfo
lapply(assay(sere, "mod_prob"), ncol)
#> $`s1-G-`
#> [1] 1
#> 
#> $`s1-GT`
#> [1] 2
#> 
#> $`s2-G-`
#> [1] 2
#> 
```
