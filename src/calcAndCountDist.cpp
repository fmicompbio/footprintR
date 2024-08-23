#include <Rcpp.h>

//' @title Count frequency of differences between values in integer vectors.
//'
//' @description Given two ascendingly sorted integer vectors \code{query} and
//'     \code{reference}, calculate and count the differences between their
//'     elements that are greater than zero and less than \code{maxd}. The
//'     number of observed distances \code{d} are reported in \code{cnt[d]},
//'     and \code{maxd} corresponds to the \code{length(cnt)}. The function is
//'     called by \code{\link{calcModbaseSpacing}}, which provides a higher level,
//'     more convenient interface.
//'
//' @author Michael Stadler
//'
//' @param query first \code{integer} vector.
//' @param reference second \code{integer} vector. Distances are calculated from
//'     each element in \code{query} to each greater element in
//'     \code{reference}.
//' @param cnt \code{NumericVector} to store the result in. The length of
//'     \code{cnt} defines the maximal distance that will be included in the
//'     analysis, and new counts will be added to the values of \code{cnt}
//'     (repeated calls to \code{calcAndCountDist} will increment existing
//'     counts).
//'
//' @return \code{numeric} vector \code{cnt}, where \code{cnt[d]} correspond to
//'   the number of observed distances \code{d}.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector calcAndCountDist(std::vector<int> query,
                                     std::vector<int> reference,
                                     Rcpp::NumericVector cnt) {
    // declare variables
    unsigned int qi, ri, currquery, d, riold = 0, maxd = cnt.size();

    // check if input is sorted
    if (!std::is_sorted(query.begin(), query.end()) ||
        !std::is_sorted(reference.begin(), reference.end())) {
        Rcpp::stop("'query' and 'reference' must be sorted ascendingly");
    }

    // for each element of query
    for (qi=0; qi < query.size(); qi++) {
        currquery = query[qi];
        // fast forward ri until currquery < reference[ri]
        while (riold < reference.size() && currquery >= reference[riold]) {
            riold++;
        }
        for (ri = riold; ri < reference.size(); ri++) {
            d = reference[ri] - currquery;
            if (d <= maxd) {
                cnt[d-1]++; // correct index for zero-based C++ vectors
            } else {
                break;
            }
        }
    }
    return cnt;
}
