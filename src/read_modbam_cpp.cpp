#include <htslib/sam.h>
#include <string>
#include <vector>
#include <set>
#include <cstdlib>
#include <cstdio>
#include <Rcpp.h>


// Helper functions
// -----------------------------------------------------------------------------
// convert 0-based read position to 0-based reference sequence position
// (a position of -1 means unaligned)
std::vector<int> read_to_reference_pos(const bam1_t *aln,
                                       const std::vector<int> &read_pos) {
    // variables
    size_t read_pos_index = 0; // index to elements of read_pos
    const uint32_t *cigar = bam_get_cigar(aln);  // cigar array
    int ref_pos = aln->core.pos;  // reference position (0-based)
    int read_seq_index = 0;  // current index in the read sequence

    // return value: 0-based reference positions, initialized to -1
    std::vector<int> ref_positions(read_pos.size(), -1);

    // iterate over the CIGAR operations i
    for (unsigned int i = 0; i < aln->core.n_cigar && read_pos_index < read_pos.size(); i++) {
        int op = bam_cigar_op(cigar[i]);  // operation type
        int op_len = bam_cigar_oplen(cigar[i]);  // operation length

        switch (op) {
        case BAM_CMATCH:  // match or mismatch (M)
        case BAM_CEQUAL:  // match (=)
        case BAM_CDIFF:   // mismatch (X)
            while (read_pos_index < read_pos.size() &&
                   read_seq_index + op_len > read_pos[read_pos_index]) {
                ref_positions[read_pos_index] = ref_pos + (read_pos[read_pos_index] - read_seq_index);
                read_pos_index++;
            }
            ref_pos += op_len;
            read_seq_index += op_len;
            break;

        case BAM_CINS:  // insertion (I)
            if (read_seq_index + op_len > read_pos[read_pos_index]) {
                // the current read position is within an insertion -->
                //     no corresponding reference position
                while (read_pos_index < read_pos.size() &&
                       read_seq_index + op_len > read_pos[read_pos_index]) {
                    ref_positions[read_pos_index] = -1;
                    read_pos_index++;
                }
            }
            read_seq_index += op_len;
            break;

        case BAM_CDEL:       // deletion (D)
        case BAM_CREF_SKIP:  // reference skip (N)
            ref_pos += op_len;
            break;

        case BAM_CSOFT_CLIP:  // soft clipping (S)
            if (read_seq_index + op_len > read_pos[read_pos_index]) {
                // the current read position is within a soft-clipped region -->
                //     no corresponding reference position
                while (read_pos_index < read_pos.size() &&
                       read_seq_index + op_len > read_pos[read_pos_index]) {
                    ref_positions[read_pos_index] = -1;
                    read_pos_index++;
                }
            }
            read_seq_index += op_len;
            break;

        case BAM_CHARD_CLIP:  // hard clipping (H)
        case BAM_CPAD:        // padding (P)
            // these do not consume any positions in the read or reference
            break;

        default: // # nocov start
            Rcpp::warning("Unknown CIGAR operation: %d", op);
        return ref_positions; // # nocov end
        }
    }

    return ref_positions;
}

// create the complement of a base
// [[Rcpp::export]]
char complement(char n) {
    switch(n) {
    case 'A':
    case 'a':
        return 'T';
    case 'T':
    case 't':
        return 'A';
    case 'G':
    case 'g':
        return 'C';
    case 'C':
    case 'c':
        return 'G';
    case 'N':
    case 'n':
    default:
        return 'N';
    }
}

// get unmodified base corresponding to a modified base `b`
// [[Rcpp::export]]
char get_unmodified_base(char b) {
    switch (b) {
    case 'm':
    case 'h':
    case 'f':
    case 'c':
    case 'C':
        return 'C';
    case 'g':
    case 'e':
    case 'b':
    case 'T':
        return 'T';
    case 'U':
        return 'U';
    case 'a':
    case 'A':
        return 'A';
    case 'o':
    case 'G':
        return 'G';
    case 'n':
    case 'N':
    default:
        return 'N';
    }
}

// calculate aligned bases (sum of 'M', '=', or 'X' operation lengths)
int calculate_aligned_bases(bam1_t *bamdata) {
    uint32_t *cigar = bam_get_cigar(bamdata);
    int aligned_bases = 0;

    // loop over CIGAR operations
    for (int i = 0; i < bamdata->core.n_cigar; i++) {
        uint32_t op = bam_cigar_op(cigar[i]);

        // only count 'M', '=', or 'X' operations
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            aligned_bases += bam_cigar_oplen(cigar[i]);
        }
    }

    return aligned_bases;
}

// process a single bam record:
// - increase alignment counter (passed by reference)
// - extract information from record (qscore, modification information, etc.)
// - add information to vectors (passed by reference) for later returning to R
int process_bam_record(bam1_t *bamdata,        // bam record
                       unsigned int &alncnt,   // alignment counter
                       char *&qseq,            // buffer for forward read sequence
                       int &qseq_len,          // allocated length of qseq
                       hts_base_mod_state *ms, // modification state struct
                       bool &had_error,        // error flag
                       char *buffer,           // buffer for message
                       int &buffer_len,        // allocated length of message buffer
                       char modbase,           // modified base to analyze
                       sam_hdr_t *in_samhdr,   // sam file header
                       int &n_unaligned,       // number of unaligned modified bases
                       int &n_total,           // total number of modified bases
                       // vectors for return values (per modification)
                       std::vector<std::string> &read_id,
                       std::vector<char> &call_code,
                       std::vector<char> &canonical_base,
                       std::vector<char> &ref_mod_strand,
                       std::vector<std::string> &chrom,
                       std::vector<int> &aligned_read_position,
                       std::vector<int> &forward_read_position,
                       std::vector<int> &ref_position,
                       std::vector<double> &mod_prob,
                       // vectors for return values (per alignment)
                       std::vector<std::string> &df_read_id,
                       std::vector<double> &df_qscore,
                       std::vector<int> &df_read_length,
                       std::vector<int> &df_aligned_length) {
    // allocate variable only used inside process_bam_record()
    uint8_t *data = NULL, *qs_data = NULL, *qual = NULL;
    unsigned int sum_qual = 0;
    double qs_value = -1;
    int i = 0, j = 0, strand = 0, impl = 0, pos = 0, r = 0;
    hts_base_mod mod[5] = {{0}};  //for ATCGN
    char canonical = '0', unmodbase = '0';
    int this_read_len = bamdata->core.l_qseq;
    size_t size_before_this_read = read_id.size();

    // get expected unmodified base corresponding to `modbase`
    unmodbase = get_unmodified_base(modbase);

    // process alignment
    alncnt++;

    // ... extract *forward* read sequence to char*
    data = bam_get_seq(bamdata);
    if (qseq_len < this_read_len) {
        if (qseq) // # nocov start
            free((void*) qseq); // # nocov end
        qseq = (char*) calloc(this_read_len + 1, sizeof(char));
        qseq_len = this_read_len;
    }
    if (bam_is_rev(bamdata)) {
        for (j = 0; j < this_read_len; j++) {
            qseq[this_read_len - 1 - j] = complement(seq_nt16_str[bam_seqi(data, j)]);
        }
    } else {
        for (j = 0; j < this_read_len; j++) {
            qseq[j] = seq_nt16_str[bam_seqi(data, j)];
        }
    }

    // ... parse base modifications
    if (bam_parse_basemod(bamdata, ms)) {
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Failed to parse the base mods (read %s)\n",
                 bam_get_qname(bamdata));
        return -1;
    }

    // ... process read if modifications of the right type are present
    //     bam_mods_query_type:
    //     - returns 0 on success, -1 if not found
    //     - also fills out `canonical`, `strand` and `impl`
    //       (`impl` is a boolean for whether unlisted positions should be
    //        implicitly assumed to be unmodified, or require an explicit
    //        score and should be considered as unknown)
    if (bam_mods_query_type(ms, modbase, &strand, &impl, &canonical) == 0) {
        // ... loop over sequence positions i
        for (i = 0; i < this_read_len; i++) {
            // i is the position in the aligned read (possibly reverse-complemented)
            // pos is the position in the original read (qseq)
            if (bam_is_rev(bamdata)) {
                pos = this_read_len - 1 - i;
            } else{
                pos = i;
            }

            // r: number of found modifications (>=1, 0 or -1 if failed)
            r = bam_mods_at_next_pos(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]));
            if (r <= -1) {
                had_error = true; // # nocov start
                snprintf(buffer, buffer_len,
                         "Failed to get modifications (read %s)\n",
                         bam_get_qname(bamdata));
                return -2; // # nocov end

            } else if (r > (int)(sizeof(mod) / sizeof(mod[0]))) {
                had_error = true;
                snprintf(buffer, buffer_len,
                         "More modifications than footprintR:::read_modbam_cpp can handle (read %s)\n",
                         bam_get_qname(bamdata));
                return -3;

            } else if (!r && impl) {
                // implied base without modification at position i
                if (qseq[pos] == unmodbase) {
                    // base of the right type -> add to results
                    read_id.push_back(bam_get_qname(bamdata));
                    aligned_read_position.push_back(i);
                    forward_read_position.push_back(pos);
                    chrom.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                    call_code.push_back('-');
                    canonical_base.push_back(canonical);
                    ref_mod_strand.push_back(bam_is_rev(bamdata) ? '-' : '+');
                    mod_prob.push_back(-1.0); // special value of -1.0 indicates inferred unmodified base
                }
            }
            // modifications
            for (j = 0; j < r; j++) {
                if (mod[j].modified_base == modbase) {
                    // found modified base of the right type -> add to results
                    read_id.push_back(bam_get_qname(bamdata));
                    aligned_read_position.push_back(i);
                    forward_read_position.push_back(pos);
                    chrom.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                    call_code.push_back((char) mod[j].modified_base);
                    canonical_base.push_back((char) mod[j].canonical_base);
                    ref_mod_strand.push_back(bam_is_rev(bamdata) == mod[j].strand ? '+' : '-');
                    // `qual` of N corresponds to call probability
                    //     in [N/256, (N+1)/256] -> store midpoint
                    mod_prob.push_back(((double) mod[j].qual + 0.5) / 256.0);
                }
            }
        }
    }

    // ... convert 0-based read positions to 0-based reference coordinates
    //     (a coordinate of -1 means unaligned, e.g. soft-masked)
    std::vector<int> aligned_read_position_converted =
        read_to_reference_pos(bamdata, aligned_read_position);
    ref_position.reserve(ref_position.size() +
        aligned_read_position_converted.size());
    ref_position.insert(ref_position.end(),
                        aligned_read_position_converted.begin(),
                        aligned_read_position_converted.end());
    aligned_read_position.clear();

    // ... remove unaligned (e.g. soft-masked) read-bases
    //     (iterate backwards to avoid messing up indices
    //      when removing elements)
    n_total += ref_position.size();
    for (size_t e = ref_position.size(); e-- > 0;) {
        if (ref_position[e] == -1) {
            n_unaligned++;
            read_id.erase(read_id.begin() + e);
            chrom.erase(chrom.begin() + e);
            forward_read_position.erase(forward_read_position.begin() + e);
            ref_position.erase(ref_position.begin() + e);
            call_code.erase(call_code.begin() + e);
            canonical_base.erase(canonical_base.begin() + e);
            ref_mod_strand.erase(ref_mod_strand.begin() + e);
            mod_prob.erase(mod_prob.begin() + e);
        }
    }

    // ... extract read-level information if the read had modified bases
    if (size_before_this_read < read_id.size()) {
        // ... extract qscore
        qs_data = bam_aux_get(bamdata, "qs");
        if (qs_data != NULL) {
            qs_value = bam_aux2f(qs_data);
        } else {
            // qs tag is missing
            //   --> calculate mean of base QUAL values
            qual = bam_get_qual(bamdata);
            sum_qual = 0;
            for (j = 0; j < this_read_len; j++) {
                sum_qual += qual[j];
            }
            qs_value = ((double) sum_qual) / this_read_len;
        }

        // ... add to read-level results
        df_read_id.push_back(bam_get_qname(bamdata));
        df_qscore.push_back(qs_value);
        df_read_length.push_back(this_read_len);
        df_aligned_length.push_back(calculate_aligned_bases(bamdata));
    }

    return 0;
}


//' Read base modifications from a bam file.
//'
//' Parse ML and MM tags (see https://samtools.github.io/hts-specs/SAMtags.pdf,
//' section 1.7) and return a list of information on modified bases.
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}
//' @param modbase Character scalar defining the modified base to extract.
//' @param n_alns_to_sample Integer defining the number of alignments
//'     to randomly sample.
//' @param tnames_for_sampling String vector with target names (chromosomes)
//'     from which to sample \code{n_alns_to_sample} alignments. Ignored if
//'     \code{n_alns_to_sample = 0}.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @return A named list with elements \code{"read_id"},
//'     \code{"forward_read_position"}, \code{"ref_position"},
//'     \code{"chrom"}, \code{"ref_mod_strand"}, \code{"call_code"},
//'     \code{"canonical_base"}, \code{"mod_prob"} and \code{"read_df"}.
//'     The meaning of these elements is described in https://nanoporetech.github.io/modkit/intro_extract.html,
//'     apart from \code{"mod_prob"}, which is equal to \code{call_prob} for
//'     modified bases and equal to \code{1 - call_prob} for unmodified bases
//'     (\code{call_code == "-"}), and \code{"read_df"}, which is a
//'      \code{data.frame} with one row per read and columns \code{"read_id"}
//'     (the read identifier), \code{"qscore"} (the read quality score recorded
//'     in the \code{qs} tag of each bam record), \code{"read_length"} (the
//'     total read length), and \code{"aligned_length"} (the number of
//'     aligned bases).
//'
//' @examples
//' modbamfile <- system.file("extdata", "6mA_1_10reads.bam", package = "footprintR")
//' res <- read_modbam_cpp(modbamfile, "chr1:6940000-6955000", "a", 0, "", TRUE)
//' str(res)
//'
//' @seealso https://samtools.github.io/hts-specs/SAMtags.pdf describing the
//'     SAM ML and MM tags for base modifications.
//'     Helpful examples are available in
//'      https://github.com/samtools/htslib/blob/develop/samples/modstate.c
//'     Documentation of the htslib C API is available in
//'      https://github.com/samtools/htslib/blob/develop/htslib/sam.h
//'     Description of returned values
//'      https://nanoporetech.github.io/modkit/intro_extract.html
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List read_modbam_cpp(std::string inname_str,
                           std::vector<std::string> regions,
                           char modbase,
                           int n_alns_to_sample,
                           std::vector<std::string> tnames_for_sampling,
                           bool verbose = false) {
    // turn htslib logging off -> handle via Rcpp::warning or Rcpp::stop
    hts_set_log_level(HTS_LOG_OFF);

    // variable declarations
    // ... general variables
    int c = 0, i = 0, success = 0;
    int n_unaligned = 0, n_total = 0;
    bool had_error = false;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    hts_base_mod_state *ms = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    unsigned int regcnt = 0, alncnt = 0;
    char **regions_c = NULL, *qseq = NULL;
    int qseq_len = 0;
    int buffer_len = 2000;
    char buffer[2000];

    // ... return values (one per modification)
    std::vector<std::string> read_id;
    std::vector<int> aligned_length;
    std::vector<char> call_code;
    std::vector<char> canonical_base;
    std::vector<char> ref_mod_strand;
    std::vector<std::string> chrom;
    std::vector<int> aligned_read_position;
    std::vector<int> forward_read_position;
    std::vector<int> ref_position;
    std::vector<double> mod_prob;

    // ... return values (one per aligned read)
    std::vector<std::string> df_read_id;
    std::vector<double> df_qscore;
    std::vector<int> df_read_length;
    std::vector<int> df_aligned_length;

    const char* inname = inname_str.c_str();

    // initialize bam data storage
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
    if (verbose) {
        snprintf(buffer, buffer_len, "    opening input file %s", inname);
        Rcpp::message(Rcpp::wrap(buffer));
    }
    if (!(infile = sam_open(inname, "r"))) {
        had_error = true;
        snprintf(buffer, buffer_len, "Could not open input file %s\n", inname);
        goto end;
    }

    // load index file
    if (!(idx = sam_index_load(infile, inname))) {
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Failed to load the index for %s\n", inname);
        goto end;
    }

    // read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len,
                 "Failed to read header from file %s\n", inname);
        goto end; // # nocov end
    }

    if (n_alns_to_sample > 0) {
        // random-sampling-based alignment reading
        // ---------------------------------------------------------------------

        // check if tnames_for_sampling exist and count alignments
        uint64_t mapped = 0, unmapped = 0, total_for_sampling = 0;
        std::set<std::string> tnames_for_sampling_set(tnames_for_sampling.begin(), tnames_for_sampling.end());
        std::set<std::string> tnames_existing;
        double rand_val = 0.0;
        regcnt = 0;
        regions_c = (char**) calloc((unsigned int) tnames_for_sampling.size(),
                                    sizeof(char*));
        for (i = 0; i < in_samhdr->n_targets; i++) {
            tnames_existing.insert(in_samhdr->target_name[i]);

            // for each target i that is in tnames_for_sampling_set,
            // get the number of mapped and unmapped records
            // and add it to regions_c
            if (tnames_for_sampling_set.find(in_samhdr->target_name[i]) !=
                  tnames_for_sampling_set.end() &&
                hts_idx_get_stat(idx, i, &mapped, &unmapped) == 0) {
                total_for_sampling += mapped;
                regions_c[regcnt] = in_samhdr->target_name[i];
                regcnt++;
            }
        }
        for (i = 0; i < (int)tnames_for_sampling.size(); i++) {
            if (tnames_existing.find(tnames_for_sampling[i]) == tnames_existing.end()) {
                Rcpp::warning("Ignoring unknown target name: %s",
                              tnames_for_sampling[i].c_str());
            }
        }

        // check if we have enough alignments to sample from
        if (total_for_sampling < (uint64_t)n_alns_to_sample) {
            had_error = true;
            snprintf(buffer, buffer_len,
                     "Cannot sample %d alignments from a total of %" PRIu64 "\n",
                     n_alns_to_sample, total_for_sampling);
            goto end;
        }
        double keep_aln_fraction = (double) n_alns_to_sample / total_for_sampling;
        if (verbose) {
            snprintf(buffer, buffer_len, "    sampling alignments with probability %g", keep_aln_fraction);
            Rcpp::message(Rcpp::wrap(buffer));
        }

        // create multi-region iterator
        if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
            had_error = true; // # nocov start
            snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
            goto end; // # nocov end
        }

        // iterate over regions
        if (verbose) {
            snprintf(buffer, buffer_len, "    reading alignments overlapping any of %u targets", regcnt);
            Rcpp::message(Rcpp::wrap(buffer));
        }
        // read overlapping alignments using iterator
        while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
            rand_val = R::runif(0, 1);
            if (!(bamdata->core.flag & BAM_FUNMAP) &&
                rand_val < keep_aln_fraction) {
                success = process_bam_record(bamdata,          // bam record
                                             alncnt,           // alignment counter
                                             qseq,             // buffer for forward read sequence
                                             qseq_len,         // allocated length of qseq
                                             ms,               // modification state struct
                                             had_error,        // error flag
                                             buffer,           // buffer for message
                                             buffer_len,       // allocated length of message buffer
                                             modbase,          // modified base to analyze
                                             in_samhdr,        // sam file header
                                             n_unaligned,      // number of unaligned modified bases
                                             n_total,          // total number of modified bases
                                             // vectors for return values (per modification)
                                             read_id,
                                             call_code,
                                             canonical_base,
                                             ref_mod_strand,
                                             chrom,
                                             aligned_read_position,
                                             forward_read_position,
                                             ref_position,
                                             mod_prob,
                                             // vectors for return values (per alignment)
                                             df_read_id,
                                             df_qscore,
                                             df_read_length,
                                             df_aligned_length);
                if (success != 0) { // # nocov start
                    goto end;
                } // # nocov end
            }
        }

    } else {
        // region-based alignment reading
        // ---------------------------------------------------------------------
        // convert regions to C arrays
        regcnt = (unsigned int) regions.size();
        regions_c = (char**) calloc(regcnt, sizeof(char*));
        for (i = 0; i < (int) regcnt; i++) {
            regions_c[i] = (char*) regions[i].c_str();
        }

        // create multi-region iterator
        if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
            had_error = true;
            snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
            goto end;
        }

        // iterate over regions
        if (verbose) {
            snprintf(buffer, buffer_len, "    reading alignments overlapping any of %u regions", regcnt);
            Rcpp::message(Rcpp::wrap(buffer));
        }
        // read overlapping alignments using iterator
        while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
            success = process_bam_record(bamdata,          // bam record
                                         alncnt,           // alignment counter
                                         qseq,             // buffer for forward read sequence
                                         qseq_len,         // allocated length of qseq
                                         ms,               // modification state struct
                                         had_error,        // error flag
                                         buffer,           // buffer for message
                                         buffer_len,       // allocated length of message buffer
                                         modbase,          // modified base to analyze
                                         in_samhdr,        // sam file header
                                         n_unaligned,      // number of unaligned modified bases
                                         n_total,          // total number of modified bases
                                         // vectors for return values (per modification)
                                         read_id,
                                         call_code,
                                         canonical_base,
                                         ref_mod_strand,
                                         chrom,
                                         aligned_read_position,
                                         forward_read_position,
                                         ref_position,
                                         mod_prob,
                                         // vectors for return values (per alignment)
                                         df_read_id,
                                         df_qscore,
                                         df_read_length,
                                         df_aligned_length);
            if (success != 0) {
                goto end;
            }
        }
    }

    if (c != -1) {
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Error while reading from %s - aborting\n", inname);
        goto end;
    }

    if (verbose) {
        snprintf(buffer, buffer_len,
                 "    removed %d unaligned (e.g. soft-masked) of %d called bases\n    read %u alignments",
                 n_unaligned, n_total, alncnt);
        Rcpp::message(Rcpp::wrap(buffer));
    }

    end:
        //cleanup
        if (qseq) { // # nocov start
            free((void*) qseq);
            qseq = NULL;
        } // # nocov end
        if (regions_c) {
            free((void*) regions_c);
            regions_c = NULL;
        }
        if (in_samhdr) {
            sam_hdr_destroy(in_samhdr);
        }
        if (infile) {
            sam_close(infile);
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

        if (had_error) {
            // we encountered an error (message in `buffer`) --> stop
            Rcpp::stop(buffer);

        } else {
            // create data.frame for read-level data
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::_["read_id"] = df_read_id,
                Rcpp::_["qscore"] = df_qscore,
                Rcpp::_["read_length"] = df_read_length,
                Rcpp::_["aligned_length"] = df_aligned_length
            );

            // create return list
            Rcpp::List res = Rcpp::List::create(
                Rcpp::_["read_id"] = read_id,
                Rcpp::_["forward_read_position"] = forward_read_position,
                Rcpp::_["ref_position"] = ref_position,
                Rcpp::_["chrom"] = chrom,
                Rcpp::_["ref_mod_strand"] = ref_mod_strand,
                Rcpp::_["call_code"] = call_code,
                Rcpp::_["canonical_base"] = canonical_base,
                Rcpp::_["mod_prob"] = mod_prob,
                Rcpp::_["read_df"] = df);

            return res;
        }
}
