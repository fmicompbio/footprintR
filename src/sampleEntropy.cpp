#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


//' @title Sample Entropy of Time series signal
//'
//' @description
//' \code{sampleEntropy} returns the sample entropy of a time-series signal
//' see also https://en.wikipedia.org/wiki/Sample_entropy
//'
//' @details
//' This function calculates the Sample Entropy of a time-series vector
//' given as argument. Sample Entropy is  used to assess the complexity of physiological
//' time-series signals.
//' This C++  implementation is a modified version of:
//' https://gist.github.com/schochastics/e3684645763e93cbc2ed7d1b70ee5fe6
//'
//' @param data  Numeric vector
//' @param m  Integer, the embedding dimension, as for chaotic time series; a preferred value is 2.
//' @param r  Scaling parameter for the filtering factor. The filtering factor is r x standard deviation of the signal
//'
//' @return The Sample Entropy value of the time-series signal
//'
//' @examples
//' ts <- runif(100,0,1)
//' sampleEntropy(ts, m=2L, r=0.2)
//'
//' @seealso [wikipedia:Sample_entropy](https://en.wikipedia.org/wiki/Sample_entropy)
//' [Multiscale entropy of biological signals](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.71.021906)
//'
//'@export
// [[Rcpp::export]]
 double sampleEntropy(NumericVector data, int m, double r) {
     int N = data.size();

     double ssum = std::accumulate(std::begin(data), std::end(data), 0.0);
     double mm =  ssum / N;
     double accum = 0.0;
     std::for_each (std::begin(data), std::end(data), [&](const double d) {
         accum += (d - mm) * (d - mm);
     });
     double sd = sqrt(accum / (N-1));

     int Cm = 0, Cm1 = 0;
     double err = 0.0;

     err = sd * r;

     for (unsigned int i = 0; i < N - (m + 1) + 1; i++) {
         for (unsigned int j = i + 1; j < N - (m + 1) + 1; j++) {
             bool eq = true;
             //m - length series
             for (unsigned int k = 0; k < m; k++) {
                 if (std::abs(data[i+k] - data[j+k]) > err) {
                     eq = false;
                     break;
                 }
             }
             if (eq) Cm++;

             //m+1 - length series
             int k = m;
             if (eq && std::abs(data[i+k] - data[j+k]) <= err)
                 Cm1++;
         }
     }

     if (Cm > 0 && Cm1 > 0)
         return std::log((double)Cm / (double)Cm1);
     else
         return 0.0;
 }

