#include <htslib/sam.h>

char get_unmodified_base(char);
char complement(char);
int calculate_aligned_bases(bam1_t*);
double extract_qscore(bam1_t*);