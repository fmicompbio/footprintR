#include <Rcpp.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <cli/progress.h>
#include "utils.h"
#include "sampleEntropy.h"

#define NMODS 5


//' Write records from \code{infile} to \code{outfile} if they pass filter criteria.
//'
//' Workhorse function for filterReadsBam. Parses records from a single
//' \code{infile}, calculate read statistics and writes the record to
//' a single \code{outfile} if the record passes all criteria defined by
//' the filtering arguments.
//' Filters are processed hierarchically: If a record does not pass a given
//' filter, the remaining filters will not be examined and the processing
//' continues with the next record.
//' The filter order is: keepUnmapped, keepSecondary, keepSupplementary,
//' minReadLength, minAlignedLength, minAlignedFraction, minQscore, maxEntropy,
//' maxFracLowConf.
//'
//' @param infile Character scalar with name of the input bam file.
//' @param outfile Character scalar with name of the output bam file.
//' @param modbase Character scalar defining the modified base to analyze
//'     (used by \code{maxEntropy} and \code{maxFracLowConf}).
//' @param region Character scalar specifying the region for which
//'     to extract overlapping reads, for example in the form
//'     \code{"chr:start-end"} (genomic interval), \code{"chr"} (all records
//'     on the given reference sequence) or \code{"."} (all records in the file).
//' @param includeBamHeader Logical scalar. If \code{TRUE} (the default), the
//'     bam header from \code{infile} will be read and written to \code{outfile}.
//'     If \code{FALSE}, no header will be written to \code{outfile}.
//' @param keepUnmapped,keepSecondary,keepSupplementary Logical scalars
//'     indicating whether to keep unmapped, secondary or supplementary
//'     alignments.
//' @param minReadLength A numeric scalar representing the smallest acceptable
//'     read length. Reads that are shorter than this value will be filtered
//'     out.
//' @param minAlignedLength A numeric scalar representing the smallest acceptable
//'     aligned length. Reads with aligned length shorter than this value will
//'     be filtered out.
//' @param minAlignedFraction A numeric scalar representing the smallest
//'     acceptable aligned fraction of a read. Reads where the aligned fraction
//'     is smaller than this value will be filtered out.
//' @param minQscore A numeric scalar representing the smallest acceptable
//'     read-level Qscore. Reads with Qscore below this value will be filtered
//'     out.
//' @param maxFracLowConf A numeric scalar representing the maximally acceptable
//'     fraction of low-confidence modified base calls in a read. Reads with
//'     a fraction of low confidence calls greater than this value will be
//'     filtered out.
//' @param maxEntropy A numeric scalar representing the largest acceptable
//'     read-level entropy. Reads with entropy above this value will be filtered
//'     out. A negative value deactivates the entropy filter.
//' @param LowConf A numeric scalar with the minimum call confidence below which
//'     calls are considered "low confidence".
//' @param nThreads Numeric scalar defining the number of threads to
//'     use for (de-)compressing bam records.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @return A named \code{numeric} vector with the numbers of filtered out
//'     records per reason for exclusion.
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector filter_modbam_cpp(std::string infile,
                                      std::string outfile,
                                      char modbase,
                                      std::string region = ".",
                                      bool includeBamHeader = true,
                                      bool keepUnmapped = true,
                                      bool keepSecondary = true,
                                      bool keepSupplementary = true,
                                      int minReadLength = 0,
                                      int minAlignedLength = 0,
                                      double minAlignedFraction = 0,
                                      double minQscore = 0.0,
                                      double maxFracLowConf = 1.0,
                                      double maxEntropy = -1.0,
                                      double LowConf = 0.7,
                                      int nThreads = 2,
                                      bool verbose = false) {
    // turn htslib logging off -> handle via Rcpp::warning or Rcpp::stop
    hts_set_log_level(HTS_LOG_OFF);

    // variable declarations
    const char *infile_c = infile.c_str(), *outfile_c = outfile.c_str();
    int c = 0, this_read_len = 0;
    char unmodbase = '0';
    bool had_error = false;
    int buffer_len = 2000;
    char buffer[2000];
    unsigned int alncnt = 0, outcnt = 0;
    char *qseq = NULL;
    int qseq_len = 0;
    Rcpp::NumericVector mod_probs = Rcpp::NumericVector(0);
    double fracLowConf = 0.0;

    // ... htslib
    bam1_t *bamdata = NULL;
    htsThreadPool tpool = {NULL, 0};
    hts_base_mod_state *ms = NULL;
    samFile *inbamfile = NULL, *outbamfile = NULL;
    sam_hdr_t *inbamhdr = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;

    // ... return values
    unsigned int nUnmapped = 0, nSecondary = 0, nSupplementary = 0,
        nMinReadLength = 0, nMinAlignedLength = 0, nMinAlignedFraction = 0,
        nMinQscore = 0, nMaxFracLowConf = 0, nMaxEntropy = 0;

    // ... cli progress bar
    Rcpp::RObject bar;

    // initialize
    if (!(bamdata = bam_init1())) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize bamdata\n");
        goto end; // # nocov end
    }
    if (!(ms = hts_base_mod_state_alloc())) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to allocate state memory\n");
        goto end; // # nocov end
    }

    // open input file
    if (!(inbamfile = sam_open(infile_c, "r"))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Could not open %s\n", infile_c);
        goto end; // # nocov end
    }

    // load index file
    if (!(idx = sam_index_load(inbamfile, infile_c))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len,
                 "Failed to load the index for %s\n", infile_c);
        goto end; // # nocov end
    }

    // open output file
    if (!(outbamfile = sam_open(outfile_c, "wb"))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Could not open %s\n", outfile_c);
        goto end; // # nocov end
    }

    // create a pool of nThreads threads...
    if (!(tpool.pool = hts_tpool_init(nThreads))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize the thread pool using {%d} threads\n", nThreads);
        goto end; // # nocov end
    }
    // ... and use it for both inbamfile and outbamfile
    if (hts_set_opt(inbamfile, HTS_OPT_THREAD_POOL, &tpool) < 0 ||
        hts_set_opt(outbamfile, HTS_OPT_THREAD_POOL, &tpool) < 0) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to set thread options\n");
        goto end; // # nocov end
    }

    // read and write header
    if (!(inbamhdr = sam_hdr_read(inbamfile))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to read header from %s\n", infile_c);
        goto end; // # nocov end
    }
    if (includeBamHeader) {
        if (sam_hdr_write(outbamfile, inbamhdr) == -1) {
            had_error = true; // # nocov start
            snprintf(buffer, buffer_len, "Failed to write header to %s\n", outfile_c);
            goto end; // # nocov end
        }
    }

    // get expected unmodified base corresponding to `modbase`
    unmodbase = get_unmodified_base(modbase);

    // create iterator
    if (!(iter = sam_itr_querys(idx, inbamhdr, region.c_str()))) {
        had_error = true;
        snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
        goto end;
    }

    // iterate over records
    if (verbose) {
        bar = cli_progress_bar(
            NA_REAL,
            Rcpp::List::create(
                Rcpp::_["clear"] = false,
                Rcpp::_["show_after"] = 0.25,
                Rcpp::_["format"] = "{cli::pb_spin} {sprintf(\"%.1f\", cli::pb_current / 1e3)} thousand records processed ({sprintf(\"%.1f /s\", cli::pb_rate_raw)}) [{cli::pb_elapsed}]"));
    }

    while ((c = sam_itr_next(inbamfile, iter, bamdata)) >= 0) {
        alncnt++;

        // extract read information
        this_read_len = bamdata->core.l_qseq;
        if (maxEntropy >= 0 || maxFracLowConf < 1.0) {
            // extract *forward* read sequence (qseq)
            //     (populates qseq and qseq_len)
            if (extract_forward_qseq(bamdata, qseq, qseq_len) != 0) {
                had_error = true; // # nocov start
                snprintf(buffer, buffer_len,
                         "Failed to extract forward read sequence (read %s)\n",
                         bam_get_qname(bamdata));
                goto end; // # nocov end
            }

            // extract modification probabilities for read
            mod_probs = Rcpp::NumericVector(0);
            if (extract_mod_probs(bamdata, modbase, unmodbase, &mod_probs, qseq,
                                  ms, buffer, buffer_len) < 0) {
                had_error = true; // # nocov start
                goto end; // # nocov end
            }
        }

        // calculate filter statistics
        // ... keepUnmapped
        if (!keepUnmapped && (bamdata->core.flag & BAM_FUNMAP)) {
            nUnmapped++;
            continue;
        }

        // ... keepSecondary
        if (!keepSecondary && (bamdata->core.flag & BAM_FSECONDARY)) {
            nSecondary++;
            continue;
        }

        // ... keepSupplementary
        if (!keepSupplementary && (bamdata->core.flag & BAM_FSUPPLEMENTARY)) {
            nSupplementary++;
            continue;
        }

        // ... minReadLength
        if (minReadLength > 0 && this_read_len < minReadLength) {
            nMinReadLength++;
            continue;
        }

        // ... minAlignedLength
        if (minAlignedLength > 0 &&
            calculate_aligned_bases(bamdata) < minAlignedLength) {
            nMinAlignedLength++;
            continue;
        }

        // ... minAlignedFraction
        if (minAlignedFraction > 0.0 &&
            (((double)calculate_aligned_bases(bamdata) / (double)this_read_len) < minAlignedFraction)) {
            nMinAlignedFraction++;
            continue;
        }

        // ... minQscore
        if (minQscore > 0.0 && extract_qscore(bamdata) < minQscore) {
            nMinQscore++;
            continue;
        }

        // ... maxFracLowConf
        if (maxFracLowConf < 1.0) {
            fracLowConf = Rcpp::sum(
                Rcpp::abs(0.5 - mod_probs) < (LowConf - 0.5)) / (double)mod_probs.size();
            if (fracLowConf > maxFracLowConf) {
                nMaxFracLowConf++;
                continue;
            }
        }

        // ... maxEntropy
        if (maxEntropy >= 0) {
            // calculate sample entropy
            if (sampleEntropy(mod_probs, 2, 0.2) > maxEntropy) {
                nMaxEntropy++;
                continue;
            }
        }

        // still here -> write to output
        outcnt++;
        if (sam_write1(outbamfile, inbamhdr, bamdata) < 0) {
            had_error = true; // # nocov start
            snprintf(buffer, buffer_len, "Failed to write output data to {.file %s}\n", outfile_c);
            goto end; // # nocov end
        }

        // report on progress
        if (verbose && CLI_SHOULD_TICK) {
            cli_progress_set(bar, (double)alncnt);
        }
        if (alncnt % 100 == 0) { // # nocov start
            R_CheckUserInterrupt();
        } // # nocov end
    }
    if (-1 == c) {
        // reached EOF
        had_error = false;
    } else {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Error in reading data\n");
        goto end; // # nocov end
    }

end:
    //clean up
    if (qseq) {
        free((void*) qseq);
        qseq = NULL;
    }
    if (inbamhdr) {
        sam_hdr_destroy(inbamhdr);
    }
    if (inbamfile) {
        sam_close(inbamfile);
    }
    if (outbamfile) {
        sam_close(outbamfile);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (ms) {
        hts_base_mod_state_free(ms);
    }
    if (iter) {
        sam_itr_destroy(iter);
    }
    if (idx) {
        hts_idx_destroy(idx);
    }
    if (tpool.pool) {
        hts_tpool_destroy(tpool.pool);
    }

    if (had_error) {
        // we encountered an error (message in `buffer`) --> stop
        // # nocov start
        Rcpp::stop(buffer);
        // # nocov end

    } else {
        // create return value
        Rcpp::NumericVector res = Rcpp::NumericVector::create(
            Rcpp::_["total"] = alncnt,
            Rcpp::_["retained"] = outcnt,
            Rcpp::_["filtered_unmapped"] = nUnmapped,
            Rcpp::_["filtered_secondary"] = nSecondary,
            Rcpp::_["filtered_supplementary"] = nSupplementary,
            Rcpp::_["filtered_minReadLength"] = nMinReadLength,
            Rcpp::_["filtered_minAlignedLength"] = nMinAlignedLength,
            Rcpp::_["filtered_minAlignedFraction"] = nMinAlignedFraction,
            Rcpp::_["filtered_minQscore"] = nMinQscore,
            Rcpp::_["filtered_maxFracLowConf"] = nMaxFracLowConf,
            Rcpp::_["filtered_maxEntropy"] = nMaxEntropy);

        return res;
    }
}

