#include <Rcpp.h>
#include <cmath>
#include <numeric>
#ifdef _OPENMP
  #include <omp.h>
#endif
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
//' @param nThreads Integer giving the number of parallel OpenMP threads to use for calculation.
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
//' @noRd
//' @keywords internal
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
double sampleEntropy(NumericVector data,
                     unsigned int m,
                     double r,
                     int nThreads = 1) {
    const unsigned int N = data.size();
    if (N <= m + 1) return 0.0;

    // Compute mean and stddev in one pass (serial, negligible cost compared to O(N^2))
    double mean = 0.0, M2 = 0.0;
    for (unsigned int i = 0; i < N; i++) {
        double delta = data[i] - mean;
        mean += delta / (i + 1);
        M2   += delta * (data[i] - mean);
    }
    double sd = std::sqrt(M2 / (N - 1));
    double err = sd * r;

    // Global counters
    unsigned long long Cm  = 0;
    unsigned long long Cm1 = 0;

    // Parallelize the outer loop
#pragma omp parallel for num_threads(nThreads) reduction(+:Cm,Cm1) schedule(dynamic)
    for (unsigned int i = 0; i <= N - (m + 1); i++) {
        for (unsigned int j = i + 1; j <= N - (m + 1); j++) {
            const double* xi = &data[i];
            const double* xj = &data[j];
            bool eq = true;

            // Compare subsequences of length m
            for (unsigned int k = 0; k < m; k++) {
                if (std::fabs(xi[k] - xj[k]) > err) {
                    eq = false;
                    break;
                }
            }

            if (eq) {
                Cm++;
                if (std::fabs(xi[m] - xj[m]) <= err)
                    Cm1++;
            }
        }
    }

    if (Cm > 0 && Cm1 > 0)
        return std::log(static_cast<double>(Cm) / Cm1);
    else
        return 0.0;
}
