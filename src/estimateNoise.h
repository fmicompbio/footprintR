#include <Rcpp.h>

Rcpp::NumericVector estimateNoise(const Rcpp::NumericVector& probs,
                                  const Rcpp::IntegerVector& read_pos,
                                  int k = 2,
                                  int min_diffs = -1);