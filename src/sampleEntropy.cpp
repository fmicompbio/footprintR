#include <Rcpp.h>
#include <cmath>
#include <numeric>
#include <algorithm>
#ifdef _OPENMP
  #include <omp.h>
// [[Rcpp::plugins(openmp)]]
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
//' @param m Integer. The pattern (embedding) length: Larger m captures
//'     finer structure but sharply reduces the number of matches,
//'     so it requires longer data and increases variance. Common
//'     choices are m = 2 or 3 for physiological time series.
//' @param r  Scaling parameter for the filtering factor. The filtering factor
//'     is r x standard deviation of the signal
//' @param maxStarts Integer specifying the maximum number of signal start
//'     positions to evaluate (default \code{1000}). If the time series has more
//'     possible starts, an evenly spaced subset of size \code{maxStarts} is used.
//'     Use \code{-1} to include all possible starts.
//' @param nThreads Integer giving the number of parallel OpenMP threads to use for calculation.
//'
//' @return The Sample Entropy value of the time-series signal
//'
//' @examples
//' ts <- runif(100, 0, 1)
//' sampleEntropy(ts, m = 2L, r = 0.2)
//'
//' @seealso [wikipedia:Sample_entropy](https://en.wikipedia.org/wiki/Sample_entropy)
//' [Multiscale entropy of biological signals](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.71.021906)
//'
//' @export
// [[Rcpp::export]]
double sampleEntropy(NumericVector data,
                     unsigned int m,
                     double r,
                     int maxStarts = 1000,
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

    // Number of possible Starting positions.
    // Each subsequence is of length m, so the last possible start is at index N - m.
    const unsigned int S = N - m;

    // Determine how many starting positions to use.
    // If maxStarts is Inf or <=0, use all starts.
    unsigned int maxI;
    if (!std::isfinite(maxStarts) || maxStarts <= 0) {
        maxI = S;  // use all starts
    } else {
        maxI = std::min(S, (unsigned int) maxStarts);
    }

    // Build list of all start indices ([0..S-1] or evenly sampled)
    std::vector<unsigned int> starts(maxI);

    // If maxStarts capped pick evenly spaced start positions:
    if (maxI < S) {
        for (size_t t = 0; t < (size_t) maxI; t++) {
            size_t idx = (t * (size_t) S) / (size_t) maxI;
            starts[t] = (unsigned int) idx;
        }
    } else{
        std::iota(starts.begin(), starts.end(), 0);
    }

    // Build an array of (data[index], index) pairs (vals),
    // where value = first value of the subsequence.
    std::vector<std::pair<double, unsigned int>> vals;
    vals.reserve(S);
    for (unsigned int i = 0; i < S; i++) {
        vals.emplace_back(data[i], i);
    }
    // Sort this array (by value) so later we can quickly find all subsequences
    // whose first element lies within +/- err of another subsequences' first element.
    std::sort(vals.begin(), vals.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });
    // Keep in a separe array (sorted_values) just the  sorted *values*
    std::vector<double> sorted_values(S);
    for (unsigned int p = 0; p < S; p++) {
        sorted_values[p] = vals[p].first;
    }

    // Global counters
    unsigned long long Cm  = 0;
    unsigned long long Cm1 = 0;

    // Parallelize the outer loop
#pragma omp parallel for num_threads(nThreads) reduction(+:Cm,Cm1) schedule(dynamic)
    for (unsigned int t = 0; t < maxI; t++) {
        const unsigned int i = starts[t];   // pick a sampled start index
        const double* xi = &data[i];

        // Find all candidates j, i.e. those with values between data[i]-err and
        // data[i]+err in the sorted_values
        const double lo_val = xi[0] - err;
        const double hi_val = xi[0] + err;

        // Find the index (lo) in sorted_values whose start value is >= (data[i] - err).
        auto lo_it = std::lower_bound(sorted_values.begin(), sorted_values.end(), lo_val);
        const unsigned int lo = (unsigned int) (lo_it - sorted_values.begin());
        // Find the index (hi) in sorted_values whose start value is <= (data[i] + err).
        auto hi_it = std::upper_bound(sorted_values.begin(), sorted_values.end(), hi_val);
        const unsigned int hi = (unsigned int) (hi_it - sorted_values.begin());

        // Scan ONLY the j indices with values in range that are also > i
        for (unsigned int p = lo; p < hi; p++) {
            const unsigned int j = vals[p].second;
            if (j <= i || j >= S) continue; // need j > i

            const double* xj = &data[j];
            bool eq = true;

            // Compare subsequences of length m
            // (skipping k=0 which is already guaranteed to be <= err)
            for (unsigned int k = 1; k < m; k++) {
                if (std::fabs(xi[k] - xj[k]) > err) {
                    eq = false;
                    break;
                }
            }

            if (eq) {
                Cm++;
                if (std::fabs(xi[m] - xj[m]) <= err) {
                    Cm1++;
                }
            }
        }
    }

    if (Cm > 0 && Cm1 > 0) {
        return std::log((double) Cm / Cm1);
    } else {
        return 0.0;
    }
}
