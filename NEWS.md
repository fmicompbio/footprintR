# footprintR 0.3.3

* Add argument `arglistFootprints` replacing and extending `footprintColors` in `plotRegion()`
* Various smaller improvements and bug fixes 

# footprintR 0.3.2

* Add removeAllNApos argument to filterReads and subsetReads
* Add BigWig plot panel
* Add regroupReads function to regroup reads by an annotation
* Rewrite and expand genome scanning framework

# footprintR 0.3.1

* Explicitly ignore supplementary and secondary alignments in all reading modes in `readModBam`
* Various smaller improvements and bug fixes

# footprintR 0.3.0

* Enable use of more flexible region specification in `readModBam()`
* Add additional supported `statistics` in `flattenReadLevelAssay()` (`Sum` and `Mean`)
* Various smaller improvements and bug fixes

# footprintR 0.2.5

* Expand readModBam to allow direct extraction of summary values (Nmod, Nvalid), and to use the pileup approach for quick extraction of read-level values
* Add `filterReadsBam()` to filter-out reads from input bam files and write reads that pass all filters to output bam files
* Add functions to statistically identify differentially modified regions (`scanForHighScoringRegions()` and helpers)

# footprintR 0.2.4

* Speed up `plotRegion()` for data sets with many reads by replacing `do.call()` with a faster alternative
* Add option to adjust facet heights to the number of reads in the facet
* Add option to squish the reads (put more than one read in a row, if possible)
* Change the `orderReads` argument for read-level plots from a logical to a character indicating the type of ordering

# footprintR 0.2.3

* Add `addFootprints()` and helper functions `calcFootprintScores()` and
  `segmentFootprintScores()`
* Add `sampleAnnot` argument to reading functions to enable inclusion of additional sample annotations
* Add `trim` argument to `readModBam`, to restrict the returned object to the specified regions
* Enable plotting of read-specific (footprint) regions

# footprintR 0.2.2

* Rename the facetBySample argument to the read-level plots to facetBy, allowing facetting by an arbitrary sample column
* Add groupBy and colorBy argument to the summary plots, for specification of arbitrary sample columns to group and color values by
* Open up for using other smoothing methods, specified via the smoothMethod argument to the summary plots. Add 'rollingMean' as a first alternative to the smoothing splines
* Add progress bars and expand messages to track progress when reading modbam files

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
