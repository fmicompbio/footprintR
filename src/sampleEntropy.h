#include <Rcpp.h>

double sampleEntropy(Rcpp::NumericVector data,
                     unsigned int m,
                     double r,
                     int maxStarts = 1000,
                     int nThreads = 1);
