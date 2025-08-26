#include <Rcpp.h>

double compute_snr_na_gap_aware(const Rcpp::NumericVector& probs,
                                const Rcpp::IntegerVector& read_pos,
                                int k = 2,
                                int min_diffs = -1,
                                double eps = 1e-3,
                                double b0 = R_NaReal,     // optional background noise coefficients
                                double b1 = R_NaReal);    // optional background noise coefficients
