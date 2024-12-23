#include <htslib/sam.h>
#include <string>
#include <vector>

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
