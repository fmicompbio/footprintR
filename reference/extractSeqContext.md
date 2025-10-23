# Extract the sequence context around positions of interest

This function will extract a sequence context of `sequenceContextWidth`
bases around the center of the regions defined in `x` from
`sequenceReference`.

## Usage

``` r
extractSeqContext(x, sequenceContextWidth, sequenceReference)
```

## Arguments

- x:

  A
  [`GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  object defining the regions of interest. Alternatively, a
  [`RangedSummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/RangedSummarizedExperiment-class.html)
  object from which regions can be extracted using `rowRanges`. The
  extracted sequences will correspond to the regions defined as
  `resize(x, width = sequenceContextWidth, fix = "center"`.

- sequenceContextWidth:

  A numeric scalar giving the width of the sequence context to be
  extracted from the reference (`sequenceReference` argument). This must
  be an odd number so that the sequence can be centered on the modified
  base. If `sequenceContextWidth = 0` (the default), no sequence context
  will be extracted.

- sequenceReference:

  A
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object, a
  [`BSgenome`](https://rdrr.io/pkg/BSgenome/man/BSgenome-class.html)
  object, or a character scalar giving the path to a fasta formatted
  file with reference sequences. If a
  [`BSgenome`](https://rdrr.io/pkg/BSgenome/man/BSgenome-class.html)
  object is provided, it will be internally converted to a
  [`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object, since the latter allows for faster sequence retrieval.

## Value

A
[`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
object of the same length as `x` with extracted sequence context. All
elements are guaranteed to have identical length (if the sequence
context extends to before the start or beyond the end of a reference
sequence, it will be padded with 'N' bases.

## See also

[`resize`](https://rdrr.io/pkg/IRanges/man/intra-range-methods.html),
[`DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)

## Author

Michael Stadler

## Examples

``` r
# file with sequence in fasta format of length 6957060
reffile <- system.file("extdata", "reference.fa.gz", package = "footprintR")

# define some regions at the end of the reference sequence
regions <- GenomicRanges::GRanges(
    "chr1", IRanges::IRanges(start = 6957060 - c(4, 2, 0),
    width = 1, names = c("a", "b", "c")))

# extract sequence context (note the padding with N's)
extractSeqContext(regions, 7, reffile)
#> DNAStringSet object of length 3:
#>     width seq                                               names               
#> [1]     7 AAAGGGG                                           a
#> [2]     7 AGGGGAN                                           b
#> [3]     7 GGGANNN                                           c
```
