#include <htslib/sam.h>

// for description of the arguments see function definitins in utils.cpp

char get_unmodified_base(char);
char complement(char);
int calculate_aligned_bases(bam1_t*);
double extract_qscore(bam1_t*);