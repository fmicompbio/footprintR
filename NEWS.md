# footprintR 0.2.0

* Harmonization of function and argument names (using camelCase)
* Use `BiocParallel::BPPARAM` in functions that support parallel execution
* Rename `addReadsSummary()` to `flattenReadLevelAssay()`

# footprintR 0.1.4

* Refactored `plotRegion()` to use a single `tracks` argument instead of `tracks.summary` and `tracks.reads`.
* Allow `filterReads()` to return the filter statistics table without subsetting the `SummarizedExperiment` object.
* Add `readModBam(..., variantPositions)` argument that adds labels to reads, enabling allele specific analyses.

# footprintR 0.1.3

* Added option to decompress bam records in `readModBam()` in parallel on multiple threads.

# footprintR 0.1.2

* Added a `NEWS.md` file to track changes to the package.
* Use `SparseArray::NaArray` assays to represent read-level modification probabilities.
