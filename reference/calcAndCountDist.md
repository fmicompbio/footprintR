# Count frequency of differences between values in integer vectors.

Given two ascendingly sorted integer vectors `query` and `reference`,
calculate and count the differences between their elements that are
greater than zero and less than `maxd`. The number of observed distances
`d` are reported in `cnt[d]`, and `maxd` corresponds to the
`length(cnt)`. The function is called by
[`calcModbaseSpacing`](https://fmicompbio.github.io/footprintR/reference/calcModbaseSpacing.md),
which provides a higher level, more convenient interface.

## Usage

``` r
calcAndCountDist(query, reference, cnt)
```

## Arguments

- query:

  first `integer` vector.

- reference:

  second `integer` vector. Distances are calculated from each element in
  `query` to each greater element in `reference`.

- cnt:

  `NumericVector` to store the result in. The length of `cnt` defines
  the maximal distance that will be included in the analysis, and new
  counts will be added to the values of `cnt` (repeated calls to
  `calcAndCountDist` will increment existing counts).

## Value

`numeric` vector `cnt`, where `cnt[d]` correspond to the number of
observed distances `d`.

## Author

Michael Stadler

## Examples

``` r
cnt <- c(0, 0, 0)
calcAndCountDist(c(1, 4, 8), c(2, 3, 4, 5, 8, 9, 10), cnt)
#> [1] 3 2 1
cnt
#> [1] 3 2 1
```
