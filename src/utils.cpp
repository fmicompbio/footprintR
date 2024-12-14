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

// calculate aligned bases (sum of 'M', '=', or 'X' operation lengths)
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
