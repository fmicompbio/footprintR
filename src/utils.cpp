#include <htslib/sam.h>
#include <string>
#include <vector>
#include <Rcpp.h>

//' Get unmodified base corresponding to a modified base
//'
//' @param b Modified base as a char
//'
//' @return The upper-case unmodified base corresponding to \code{b} as a
//'     \code{char}.
//' @noRd
//' @keywords internal
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

//' Create the complement of a base
//'
//' @param n single base as a char
//'
//' @return char (complement of \code{n})
//'
//' @noRd
//' @keywords internal
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

//' Calculate aligned bases (sum of 'M', '=', or 'X' operation lengths)
//'
//' @param bamdata A \code{bam1_t*} with the alignment.
//'
//' @return An \code{int} giving the number of aligned bases.
//'
//' @noRd
//' @keywords internal
int calculate_aligned_bases(bam1_t *bamdata) {
    uint32_t *cigar = bam_get_cigar(bamdata);
    int aligned_bases = 0;

    // loop over CIGAR operations
    for (uint32_t i = 0; i < bamdata->core.n_cigar; i++) {
        uint32_t op = bam_cigar_op(cigar[i]);

        // only count 'M', '=', or 'X' operations
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            aligned_bases += bam_cigar_oplen(cigar[i]);
        }
    }

    return aligned_bases;
}

//' Extract quality score (qscore)
//'
//' @param bamdata A \code{bam1_t*} with the alignment.
//'
//' @return A \code{double} corresponding to the value extracted from the "qs"
//'     tag, or in case that is missing, calculated as the mean of base quality
//'     values.
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
double extract_qscore(bam1_t *bamdata) {
    uint8_t *qual = NULL, *qs_data = bam_aux_get(bamdata, "qs");
    double qs_value = 0.0, sum_qual = 0.0;
    if (qs_data != NULL) {
        qs_value = bam_aux2f(qs_data);
    } else {
        // qs tag is missing --> calculate mean of base QUAL values
        qual = bam_get_qual(bamdata);
        sum_qual = 0;
        for (int j = 0; j < bamdata->core.l_qseq; j++) {
            sum_qual += qual[j];
        }
        qs_value = ((double) sum_qual) / bamdata->core.l_qseq;
    }
    return qs_value;
}

//' Get the forward read sequence form an alignment
//'
//' Extract the read sequence from a bam1_t corresponding to the plus-strand
//' of the read (thus reverse-complementing the read for an minus-strand
//' alignment) and write it to the char* array at qseq, allocating memory of
//' sufficient length if needed. The allocated space (without terminating null
//' character) is stored in qseq_len.
//'
//' @param bamdata A \code{bam1_t*} with the alignment.
//' @param qseq A \code{char**} (pointer to a character array) to which the
//'     extracted sequence will be written.
//' @param qseq_len A \code{int*} (pointer to int) in which the number of
//'     allocated characters at \code{qseq} are stored (escluding the
//'     terminating null character).
//'
//' @returns 0 if sucessful, -1 if memory allocation failed
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
int extract_forward_qseq(bam1_t *bamdata, // alignment
                         char *&qseq,     // buffer for forward read sequence
                         int &qseq_len) { // allocated length of qseq
    uint8_t *data = bam_get_seq(bamdata);
    int this_read_len = bamdata->core.l_qseq, j = 0;

    if (qseq_len < this_read_len) {
        if (qseq) // # nocov start
            free((void*) qseq); // # nocov end
        qseq = (char*) calloc(this_read_len + 1, sizeof(char));
        if (qseq == NULL) // # nocov start
            return -1; // # nocov end
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
    return 0;
}

//' Extract vector with modification probabilities from alignment
//'
//' Use htslib functions to parse the modification probabilities for
//' `modbase`.
//'
//' @param bamdata A \code{bam1_t*} with the alignment.
//' @param modbase A \code{char} with the modified base code for which to
//'     extract modification probabilities.
//' @param unmodbase A \code{char} with the unmodified base corresponding to
//'     \code{modbase}.
//' @param mod_probs A \code{Rcpp::NumericVector*} to which the extracted
//'     modification probabilities will be appended at the end.
//' @param qseq A \code{char*} pointing to the forward read sequence.
//' @param ms A \code{hts_base_mod_state*} (modification state struct) expected
//'     to be pre-initialized.
//' @param buffer A \code{char*} pointing to a pre-allocated character array
//'     to which an error message is written in case of a failure.
//' @param buffer_len An \code{int} giving the pre-allocated size of the array
//'     at \code{buffer} (excluding the terminating null).
//'
//' @returns An \code{int}, if greater or equal to zero giving the number of
//'     extracted probabilities, or less than zero if something failed. In
//'     that case, the error message is giving in \code{buffer}.
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
int extract_mod_probs(bam1_t *bamdata,
                      char modbase,
                      char unmodbase,
                      Rcpp::NumericVector *mod_probs,
                      char* qseq,
                      hts_base_mod_state *ms,
                      char* buffer,
                      int buffer_len) {
    // declare variables
    int i = 0, j = 0, strand = 0, impl = 0, pos = 0, r = 0;
    int this_read_len = bamdata->core.l_qseq, n_probs = 0;
    hts_base_mod mod[5] = {{0}};  //for ATCGN
    char canonical = '0';

    // parse base modifications
    if (bam_parse_basemod(bamdata, ms)) { // # nocov start
        snprintf(buffer, buffer_len, "Failed to parse the base mods (read %s)\n",
                 bam_get_qname(bamdata));
        return -1; // # nocov end
    }

    // process read if modifications of the right type are present
    // bam_mods_query_type:
    // - returns 0 on success, -1 if not found
    // - also fills out `canonical`, `strand` and `impl`
    //   (`impl` is a boolean for whether unlisted positions should be
    //    implicitly assumed to be unmodified, or require an explicit
    //    score and should be considered as unknown)
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
            if (r <= -1) { // # nocov start
                snprintf(buffer, buffer_len, "Failed to get modifications (read %s)\n",
                         bam_get_qname(bamdata));
                return -2; // # nocov end
            } else if (r > (int)(sizeof(mod) / sizeof(mod[0]))) { // # nocov start
                snprintf(buffer, buffer_len,
                         "More modifications than footprintR:::extract_mod_probs can handle (read %s)\n",
                         bam_get_qname(bamdata));
                return -3; // # nocov end
            } else if (!r && impl) {
                // implied base without modification at position i
                if (qseq[pos] == unmodbase) {
                    // base of the right type -> add to results
                    mod_probs->push_back(0.0);
                    n_probs++;
                }
            }
            // modifications
            for (j = 0; j < r; j++) {
                if (mod[j].modified_base == modbase) {
                    // found modified base of the right type -> add to results
                    // `qual` of N corresponds to call probability
                    //     in [N/256, (N+1)/256] -> store midpoint
                    mod_probs->push_back(((double) mod[j].qual + 0.5) / 256.0);
                    n_probs++;
                }
            }
        }
    }

    return n_probs;
}
