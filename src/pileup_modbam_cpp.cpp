/*
 The pileup_modbam_cpp function is in part based on pileup_mod.c (distributed
 with htslib) and subject to the following copyright and permission notice:

    pileup_mod.c --  showcases the htslib api usage

    Copyright (C) 2023 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE

*/

#include <string>
#include <vector>
#include <unistd.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <Rcpp.h>
#include <cli/progress.h>
#include "utils.h"

typedef struct plpconf {
    char *inname;
    samFile *infile;
    sam_hdr_t *in_samhdr;
    hts_idx_t *idx;
    hts_itr_t *iter;
} plpconf;

//' Constructor for pileup data in bam_pileup_cd*
//'
//' @param data void* (client data)
//' @param b bam1_t* (bam being loaded)
//' @param cd bam_pileup_cd* (client data)
//'
//' @return An integer scalar (zero on success, non-zero on failure)
//'
//' @noRd
//' @keywords internal
int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    //plpconf *conf= (plpconf*)data; can use this to access anything required from the data in pileup init

    //when using cd, initialize and use as it will be reused after destructor
    cd->p = hts_base_mod_state_alloc();
    if (!cd->p) {
        // # nocov start
        Rcpp::stop("Failed to allocate base modification state\n");
        return 1;
        // # nocov end
    }

    //parse the bam data and gather modification data from MM tags
    return (-1 == bam_parse_basemod(b, (hts_base_mod_state*)cd->p)) ? 1 : 0;
}

//' Destructor for pileup data in bam_pileup_cd*
//'
//' @param data void* (client data)
//' @param b bam1_t* (bam being loaded)
//' @param cd bam_pileup_cd* (client data)
//'
//' @return An integer scalar (zero)
//'
//' @noRd
//' @keywords internal
int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    if (cd->p) {
        hts_base_mod_state_free((hts_base_mod_state *)cd->p);
        cd->p = NULL;
    }
    return 0;
}

//' Read alignment data for pileup operation
//'
//' @param data void* (client callback data holding alignment file handle)
//' @param b bam1_t* (aligned read)
//'
//' @return same as sam_read1
//'
//' @noRd
//' @keywords internal
int readdata(void *data, bam1_t *b)
{
    plpconf *conf = (plpconf*)data;
    if (!conf || !conf->infile) {
        // # nocov start
        return -2;  //cant read data
        // # nocov end
    }

    //read alignment and send
    // return sam_read1(conf->infile, conf->infile->bam_header, b);
    return sam_itr_next(conf->infile, conf->iter, b);
}

//' Read and pile-up base modifications from a bam file.
//'
//' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
//' section 1.7) and return a list of information on modified bases.
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}.
//'     The strings are interpreted by htslib, which understands:
//'     \describe{
//'         \item{"REF" or "REF:"}{: All reads with RNAME REF}
//'         \item{"REF:START"}{: Reads with RNAME REF overlapping START to end of REF}
//'         \item{"REF:-END"}{: Reads with RNAME REF overlapping start of REF to END}
//'         \item{"REF:START-END"}{: Reads with RNAME REF overlapping START to END}
//'         \item{"."}{: All reads from the start of the file}
//'         \item{"*"}{: Unmapped reads at the end of the file (RNAME '*' in SAM)}
//'     }
//' @param modbase Character scalar defining the modified base to extract.
//'     Only modifications corresponding to \code{modbase} and with the
//'     corresponding expected base in the read sequence will be extracted.
//' @param mod_prob_thresh Double scalar defining the minimal mod_prob
//'     of a base to be considered modified.
//' @param n_threads Integer scalar defining the number of threads to
//'     use for decompressing a sam record. Especially using in sampling mode
//'     (\code{n_alns_to_sample > 0}), where more time is spend reading and
//'     decompressing bam records than processing them.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @return A named list with elements \code{"chrom"} (chromosome name),
//'     \code{"ref_position"} (1-based coordinate on \code{"chrom"}),
//'     \code{"ref_mod_strand"} (the strand relative to the reference on which
//'     the modification was identified),  \code{"Nmod"} (number of modified
//'     bases) and \code{"Nvalid"} (number of total bases).
//'
//' @examples
//' modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
//' res <- pileup_modbam_cpp(modbamfile, "chr1", "a", 0.7, 1, TRUE)
//' str(res)
//'
//' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
//'     SAM ML and MM tags for base modifications.
//'     Helpful examples are available in
//'      https://github.com/samtools/htslib/blob/develop/samples/pileup_mod.c
//'
//' @author Michael Stadler
//'
//' @importFrom cli cli_progress_step cli_progress_done
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List pileup_modbam_cpp(std::string inname_str,
                             std::vector<std::string> regions,
                             char modbase,
                             double mod_prob_thresh = 0.5,
                             int n_threads = 2,
                             bool verbose = false) {
    // variable declarations
    bam1_t *bamdata = NULL;
    plpconf conf = {0};
    conf.inname = (char*)inname_str.c_str();
    bam_plp_t plpiter = NULL;
    int tid = -1, depth = -1, j = 0, modlen = 0;
    #define NMODS 5
    hts_base_mod mods[NMODS] = {{0}}; //ACGTN
    int refpos = -1;
    const bam_pileup1_t *plp = NULL;
    kstring_t insdata = KS_INITIALIZE;
    bool had_error = false;
    int buffer_len = 2000;
    char buffer[2000];
    char readbase = '0';
    uint64_t refposcount = 0;
    unsigned int regcnt = 0;
    char **regions_c = NULL;

    std::vector<std::string> chrom;
    std::vector<int> ref_position;
    std::vector<char> ref_mod_strand;
    std::vector<int> Nmod;
    std::vector<int> Nvalid;
    int curr_Nmod[2] = {0, 0};   // for +/- strand modification counts
    int curr_Nvalid[2] = {0, 0};
    int curr_strand = 0;

    // ... cli progress bar
    Rcpp::RObject bar;

    // get expected unmodified base corresponding to modbase
    char unmodbase = get_unmodified_base(modbase);
    char unmodbase_complement = complement(unmodbase);

    // initialize
    if (!(bamdata = bam_init1())) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize bamdata\n");
        goto end; // # nocov end
    }

    // open input files
    if (!(conf.infile = sam_open(conf.inname, "r"))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Could not open %s\n", conf.inname);
        goto end; // # nocov end
    }

    // load index file
    if (!(conf.idx = sam_index_load(conf.infile, conf.inname))) {
        // # nocov start
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Failed to load the index for %s\n", conf.inname);
        goto end;
        // # nocov end
    }

    // read header
    if (!(conf.in_samhdr = sam_hdr_read(conf.infile))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to read header from file!\n");
        goto end; // # nocov end
    }

    // convert regions to C arrays
    regcnt = (unsigned int) regions.size();
    regions_c = (char**) calloc(regcnt, sizeof(char*));
    for (int i = 0; i < (int) regcnt; i++) {
        regions_c[i] = (char*) regions[i].c_str();
    }

    // create multi-region iterator
    if (!(conf.iter = sam_itr_regarray(conf.idx, conf.in_samhdr, regions_c, regcnt))) {
        // # nocov start
        had_error = true;
        snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
        goto end;
        // # nocov end
    }

    // initialize pileup iterator
    if (!(plpiter = bam_plp_init(readdata, &conf))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize pileup data\n");
        goto end; // # nocov end
    }

    // set constructor and destructor callbacks
    bam_plp_constructor(plpiter, plpconstructor);
    bam_plp_destructor(plpiter, plpdestructor);

    if (verbose) {
        bar = cli_progress_bar(
            NA_REAL,
            Rcpp::List::create(
                Rcpp::_["clear"] = false,
                Rcpp::_["show_after"] = 0.25,
                Rcpp::_["format"] = "{cli::pb_spin} {sprintf(\"%.3f\", cli::pb_current / 1e6)} Mio. genomic positions processed ({sprintf(\"%.3f Mio./s\", cli::pb_rate_raw / 1e6)}) [{cli::pb_elapsed}]"));
    }

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &depth))) {
        memset(&mods, 0, sizeof(mods));
        curr_Nmod[0] = 0;
        curr_Nmod[1] = 0;
        curr_Nvalid[0] = 0;
        curr_Nvalid[1] = 0;

        // iterate over reads overlapping refpos
        for (j = 0; j < depth; ++j) {

            if (plp[j].is_del || plp[j].is_refskip ||
                (plp[j].b->core.flag & BAM_FSECONDARY) ||
                (plp[j].b->core.flag & BAM_FSUPPLEMENTARY)) {
                continue;
            }

            // check if read j has expected base
            readbase = toupper(seq_nt16_str[bam_seqi(bam_get_seq(plp[j].b),
                                                     plp[j].qpos)]);
            if (readbase != (bam_is_rev(plp[j].b) ? unmodbase_complement : unmodbase)) {
                continue;
            }

            // retrieve base modification
            if ((modlen = bam_mods_at_qpos(plp[j].b, plp[j].qpos,
                                           (hts_base_mod_state*)plp[j].cd.p,
                                           mods, NMODS)) == -1) {
                // # nocov start
                had_error = true;
                snprintf(buffer, buffer_len,
                         "Failed to get modifications from %s on %s (qpos=%d, refpos=%d)\n",
                         bam_get_qname(plp[j].b), sam_hdr_tid2name(conf.in_samhdr, tid),
                         plp[j].qpos, refpos);
                goto end;
                // # nocov end
            }

            // increment curr_Nmod[curr_strand] if base is modified and has the expected base
            if (modlen > 0) {
                curr_strand = bam_is_rev(plp[j].b) == mods[0].strand ? 0 : 1;
                curr_Nvalid[curr_strand]++;
                if ((((double) mods[0].qual + 0.5) / 256.0) >= mod_prob_thresh) {
                    curr_Nmod[curr_strand]++;
                }
            } else {
                curr_strand = bam_is_rev(plp[j].b) ? 1 : 0;
                curr_Nvalid[curr_strand]++;
            }

        }

        // add counters for refpos to return value vectors
        // ... plus strand
        if (curr_Nvalid[0] > 0) {
            chrom.push_back(sam_hdr_tid2name(conf.in_samhdr, tid));
            ref_position.push_back(refpos + 1);
            ref_mod_strand.push_back('+');
            Nmod.push_back(curr_Nmod[0]);
            Nvalid.push_back(curr_Nvalid[0]);

        }

        // ... minus strand
        if (curr_Nvalid[1] > 0) {
            chrom.push_back(sam_hdr_tid2name(conf.in_samhdr, tid));
            ref_position.push_back(refpos + 1);
            ref_mod_strand.push_back('-');
            Nmod.push_back(curr_Nmod[1]);
            Nvalid.push_back(curr_Nvalid[1]);
        }

        refposcount++;
        if (verbose && CLI_SHOULD_TICK) {
            // # nocov start
            cli_progress_set(bar, (double)refposcount);
            // # nocov end
        }
        if (refposcount % 1000000 == 0) { // # nocov start
            R_CheckUserInterrupt();
        } // # nocov end
    }

end:
    //clean up
    if (conf.in_samhdr) {
        sam_hdr_destroy(conf.in_samhdr);
    }
    if (conf.infile) {
        sam_close(conf.infile);
    }
    if (conf.iter) {
        sam_itr_destroy(conf.iter);
    }
    if (conf.idx) {
        hts_idx_destroy(conf.idx);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (plpiter) {
        bam_plp_destroy(plpiter);
    }
    ks_free(&insdata);

    if (had_error) {
        // we encountered an error (message in `buffer`) --> stop
        // # nocov start
        Rcpp::stop(buffer);
        // # nocov end

    } else {
        // create return list
        Rcpp::List res = Rcpp::List::create(
            Rcpp::_["chrom"] = chrom,
            Rcpp::_["ref_position"] = ref_position,
            Rcpp::_["ref_mod_strand"] = ref_mod_strand,
            Rcpp::_["Nmod"] = Nmod,
            Rcpp::_["Nvalid"] = Nvalid);

        return res;
    }
}
