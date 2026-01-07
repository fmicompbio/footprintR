# Package index

## Filtering BAM files

- [`filterReadsBam()`](https://fmicompbio.github.io/footprintR/reference/filterReadsBam.md)
  :

  Write bam records from `infile` to `outfile` if they pass filter
  criteria.

## Reading data into R

- [`readBedMethyl()`](https://fmicompbio.github.io/footprintR/reference/readBedMethyl.md)
  :

  Read collapsed single-molecule footprinting data from a `bedMethyl`
  file

- [`readModBam()`](https://fmicompbio.github.io/footprintR/reference/readModBam.md)
  : Read base modifications from bam file(s)

- [`readModkitExtract()`](https://fmicompbio.github.io/footprintR/reference/readModkitExtract.md)
  : Read modkit extract file(s)

- [`countStatePairs()`](https://fmicompbio.github.io/footprintR/reference/countStatePairs.md)
  : Count pairs of modified bases by distance and modification state

## Calculating on data in R

- [`calcReadStats()`](https://fmicompbio.github.io/footprintR/reference/calcReadStats.md)
  [`addReadStats()`](https://fmicompbio.github.io/footprintR/reference/calcReadStats.md)
  : Calculate or add summary statistics for read-level base modification
  data
- [`sampleEntropy()`](https://fmicompbio.github.io/footprintR/reference/sampleEntropy.md)
  : Sample Entropy of Time series signal
- [`addSeqContext()`](https://fmicompbio.github.io/footprintR/reference/addSeqContext.md)
  : Add sequence context around positions of interest to a
  SummarizedExperiment
- [`flattenReadLevelAssay()`](https://fmicompbio.github.io/footprintR/reference/flattenReadLevelAssay.md)
  : Summarize a read-level object to sample-level
- [`filterPositions()`](https://fmicompbio.github.io/footprintR/reference/filterPositions.md)
  : Filter positions
- [`filterReads()`](https://fmicompbio.github.io/footprintR/reference/filterReads.md)
  : Filter reads
- [`subsetReads()`](https://fmicompbio.github.io/footprintR/reference/subsetReads.md)
  : Subset the reads from read-level assays
- [`regroupReads()`](https://fmicompbio.github.io/footprintR/reference/regroupReads.md)
  [`regroupReadsByColData()`](https://fmicompbio.github.io/footprintR/reference/regroupReads.md)
  : Regroup reads

## Footprint identification

- [`addFootprints()`](https://fmicompbio.github.io/footprintR/reference/footprintScoring.md)
  [`calcFootprintScores()`](https://fmicompbio.github.io/footprintR/reference/footprintScoring.md)
  [`segmentFootprintScores()`](https://fmicompbio.github.io/footprintR/reference/footprintScoring.md)
  : Calculate and segment scores for footprints and add them to a
  SummarizedExeriment

## Score (differential) modifications genome-wide

- [`quantifyWindowsInRegion()`](https://fmicompbio.github.io/footprintR/reference/quantifyWindowsInRegion.md)
  : Generate counts for sequential windows in a single region

- [`sumNmodNvalid()`](https://fmicompbio.github.io/footprintR/reference/sumNmodNvalid.md)
  :

  Quantify windows by summing across positions in each of `Nmod` and
  `Nvalid`

- [`strandDiffFracMod()`](https://fmicompbio.github.io/footprintR/reference/strandDiffFracMod.md)
  : Quantify windows by calculating the difference between modification
  fractions for the positive and negative strand

- [`phasingScoreFourier()`](https://fmicompbio.github.io/footprintR/reference/phasingScoreFourier.md)
  : Quantify windows by calculating a Fourier-transform-based phasing
  score

- [`estimateNRLwindows()`](https://fmicompbio.github.io/footprintR/reference/estimateNRLwindows.md)
  : Quantify windows by estimating nucleosome repeat lengths (NRLs)

- [`estimateNoiseParsWindows()`](https://fmicompbio.github.io/footprintR/reference/estimateNoiseParsWindows.md)
  : Estimate a background noise model

- [`estimateSNRwindows()`](https://fmicompbio.github.io/footprintR/reference/estimateSNRwindows.md)
  : Quantify footprint coherence in windows by estimating the Signal to
  Noise Ratio (SNR)

- [`getDifferentiallyModifiedWindows()`](https://fmicompbio.github.io/footprintR/reference/getDifferentiallyModifiedWindows.md)
  : Analyze counts for sequential windows in a single region

- [`getDifferentialWindows()`](https://fmicompbio.github.io/footprintR/reference/getDifferentialWindows.md)
  : Perform differential analysis on the rows of an assay

- [`getRangesWithAssayValues()`](https://fmicompbio.github.io/footprintR/reference/getRangesWithAssayValues.md)
  :

  Create a `GRanges` from a `SummarizedExperiment`.

- [`processWindowScores()`](https://fmicompbio.github.io/footprintR/reference/processWindowScores.md)
  : Process scores of sequential genomic windows.

- [`scanForHighScoringRegions()`](https://fmicompbio.github.io/footprintR/reference/scanForHighScoringRegions.md)
  : Scan one or more chromosomes for high-scoring regions

## Visualization

- [`plotReadStats()`](https://fmicompbio.github.io/footprintR/reference/plotReadStats.md)
  : Plot distribution of QC statistics
- [`plotNoisePars()`](https://fmicompbio.github.io/footprintR/reference/plotNoisePars.md)
  : Plot background noise model fit diagnostics.
- [`plotRegion()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  [`plotBigWig()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  [`plotReadsLollipop()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  [`plotReadsHeatmap()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  [`plotSummaryPointSmooth()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  [`plotGenomicRegions()`](https://fmicompbio.github.io/footprintR/reference/plotRegion.md)
  : Plot single-molecule footprinting data for a single genomic region
- [`plotValsBySeqContext()`](https://fmicompbio.github.io/footprintR/reference/plotValsBySeqContext.md)
  : Plot assay values stratified by sequence context

## Helper functions

- [`extractSeqContext()`](https://fmicompbio.github.io/footprintR/reference/extractSeqContext.md)
  : Extract the sequence context around positions of interest
- [`getAnchorRegions()`](https://fmicompbio.github.io/footprintR/reference/getAnchorRegions.md)
  : Extract data for one or more anchor regions

## Nucleosome repeat length (NRL) analysis

Functions for calculating and visualizing distance distributions between
modified bases and estimate NRL.

- [`calcModbaseSpacing()`](https://fmicompbio.github.io/footprintR/reference/calcModbaseSpacing.md)
  : Calculate distances between modified bases on individual reads
- [`estimateNRL()`](https://fmicompbio.github.io/footprintR/reference/estimateNRL.md)
  : Estimate the nucleosome repeat length (NRL) from modified-base
  distances
- [`plotModbaseSpacing()`](https://fmicompbio.github.io/footprintR/reference/plotModbaseSpacing.md)
  : Plot annotated distances between modified bases
- [`calcAndCountDist()`](https://fmicompbio.github.io/footprintR/reference/calcAndCountDist.md)
  : Count frequency of differences between values in integer vectors.
