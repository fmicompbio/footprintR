#include <htslib/sam.h>
#include <Rcpp.h>

// for description of the arguments see function definitins in utils.cpp

char get_unmodified_base(char);
char complement(char);
int calculate_aligned_bases(bam1_t*);
double extract_qscore(bam1_t*);
int extract_forward_qseq(bam1_t*, char*&, int&);
Rcpp::NumericVector extract_mod_probs(bam1_t*, char, char, char*, hts_base_mod_state*);
