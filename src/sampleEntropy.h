#include <Rcpp.h>

double sampleEntropy(Rcpp::NumericVector data,
                     unsigned int m,
                     double r,
                     int nThreads = 1);
