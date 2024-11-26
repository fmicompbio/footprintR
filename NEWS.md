# footprintR 0.2.1

* Refactor `calcReadStats()` (non-user facing changes)
* Add genomic region track type to `plotRegion()` and expose more arguments for the individual track types
* Export helper functions to create individual `plotRegion()` tracks
* Modify `plotRegion()` to display exactly the provided input region (with minimal padding, and not restricting to the subregion with observed data)
* Add support for relative coordinate system in `plotRegion()`
* Use 0.5 instead of 1.0 as the label distance between non-overlapping reads
* Modify `getAnchorRegions()` to return a `RangedSummarizedExperiment` object, allow the user to set the anchor name

# footprintR 0.2.0

* Harmonize function and argument names (using camelCase)
* Use `BiocParallel::BPPARAM` in functions that support parallel execution
* Rename `addReadsSummary()` to `flattenReadLevelAssay()`

# footprintR 0.1.4

* Refactor `plotRegion()` to use a single `tracks` argument instead of `tracks.summary` and `tracks.reads`.
* Allow `filterReads()` to return the filter statistics table without subsetting the `SummarizedExperiment` object.
* Add `readModBam(..., variantPositions)` argument that adds labels to reads, enabling allele specific analyses.

# footprintR 0.1.3

* Add option to decompress bam records in `readModBam()` in parallel on multiple threads.

# footprintR 0.1.2

* Add a `NEWS.md` file to track changes to the package.
* Use `SparseArray::NaArray` assays to represent read-level modification probabilities.
