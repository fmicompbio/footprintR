#include <string>
#include <Rcpp.h>
#include <htslib/sam.h>

//' Create an index for a given bam file
//'
//' The bam file is expected to be already sorted by coordinate and
//' the index file name will be automatically determined by appending
//' \code{.bai} to the bam file name.
//'
//' @param infile A \code{std::string} with the path and name to the input
//'     bam file to be indexed.
//'
//' @returns A \code{std::string} with the name of the created index file.
//'
//'
//' @author Michael Stadler
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
std::string index_bam_cpp(std::string infile) {
    // turn htslib logging off -> handle via Rcpp::warning or Rcpp::stop
    hts_set_log_level(HTS_LOG_OFF);

    if (sam_index_build(infile.c_str(), 0) < 0) {
        Rcpp::stop("Failed creating index for bam file");
    }

    return (infile + ".bai");
}
