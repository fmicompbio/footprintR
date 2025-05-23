#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <string>
#include <vector>
#include <Rcpp.h>

// for description of the arguments see function definitions in utils.cpp

std::string concatenate_files(std::vector<std::string>, const std::string);
std::string concatenate_hts_files(std::vector<std::string>, const std::string, int);
Rcpp::CharacterVector getChromosomeNamesFromBam(const std::string);
char get_unmodified_base(char);
char complement(char);
int calculate_aligned_bases(bam1_t*);
double extract_qscore(bam1_t*);
int extract_forward_qseq(bam1_t*, char*&, int&);
int extract_mod_probs(bam1_t*, char, char, Rcpp::NumericVector*, char*, hts_base_mod_state*, char*, int);
