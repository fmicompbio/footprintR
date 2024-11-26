#include <Rcpp.h>
#include <string>
#include <algorithm>


//' @title Calculate pairwise distances between read labels
//'
//' @description
//' \code{labelDists} returns all pairwise distances among a set of strings
//'    (read labels) of identical length, consisting of A, C, G, T and -
//'    letters. The distance is in [0, 1] and corresponds to the fraction of
//'    differences in the overlap range, defined as the range excluding the
//'    maximal number of leading and trailing - letters in any of the two
//'    compared labels.
//'
//' @param labels  Character vector of equally sized strings.
//' @param minOverlap An integer scalar giving the minimal number of overlapping
//'     letters. If a pair of labels overlap by less than this number of letters,
//'     their distance is set to the maximal distance (1.0).
//'
//' @return A numeric, symmetric matrix with pairwise distances between the
//'     elements of \code{labels}.
//'
//' @examples
//' labelDists(c("--AACACT-", "---ACCCT-", "---ACAC--"))
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix labelDists(std::vector<std::string> labels,
                               int minOverlap = 2) {
    int n = (int)labels.size(), len = (int)labels[0].size();
    int i = 0, j = 0, k = 0, kfrom = 0, kto = 0;
    double ndiff = 0.0;
    Rcpp::NumericMatrix dist(n, n);

    // count number of leading and trailing dashes in each label[i]
     std::vector<int> nDashLeading(n), nDashTrailing(n);
     for (i = 0; i < n; i ++) {
         // check length
         if (labels[i].size() != len) {
             Rcpp::stop("labels[%u] (%s) does not have %u characters",
                        i + 1, labels[i].c_str(), len);
         }

         // count leading dashes
         for (j = 0; j < len; j++) {
            if (labels[i][j] == '-') {
                nDashLeading[i]++;
            } else {
                break;
            }
         }

         // count trailing dashes
         for (j = len - 1; j >= 0; j--) {
             if (labels[i][j] == '-') {
                 nDashTrailing[i]++;
             } else {
                 break;
             }
         }
     }

     // calculate the distances between labels[i] and labels[j]
     for (i = 0; i < n - 1; i ++) {
         for (j = i + 1; j < n; j ++) {
             kfrom = std::max(nDashLeading[i], nDashLeading[j]);
             kto = len - std::max(nDashTrailing[i], nDashTrailing[j]);
             if (kto - kfrom < minOverlap) {
                 dist(i, j) = dist(j, i) = 0.5;
             } else {
                 ndiff = 0.0;
                 for (k = kfrom; k < kto; k++) {
                     if (labels[i][k] != labels[j][k]) {
                         ndiff++;
                     }
                 }
                 dist(i, j) = dist(j, i) = ndiff / ((double) kto - kfrom);
             }
         }
     }

     return dist;
 }

