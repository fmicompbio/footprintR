#include <algorithm>  // std::max
#include <cmath>      // std::floor, std::log2
#include <limits>     // std::numeric_limits
#include <Rcpp.h>     // Rcpp::NumericVector / IntegerVector

double compute_snr_na_gap_aware(const Rcpp::NumericVector& probs,
                                const Rcpp::IntegerVector& read_pos,
                                int k,
                                int min_diffs,
                                double eps) {
    const int n = probs.size();
    if (n < 2 || read_pos.size() != n) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (min_diffs < 0) {
        min_diffs = std::max(16, (int)std::floor(0.05 * (double)n));
    }
    double sumx = 0.0, sumx2 = 0.0;
    for (int i = 0; i < n; ++i) { sumx += probs[i]; sumx2 += probs[i]*probs[i]; }
    const double meanx  = sumx / (double)n;
    const double totalV = (sumx2 - (double)n * meanx * meanx) / (double)(n - 1);
    int nd = 0; double sumd = 0.0, sumd2 = 0.0;
    for (int i = 1; i < n; ++i) {
        const int gap = (read_pos[i] - read_pos[i - 1]) - 1;
        if (gap <= k) {
            const double d = probs[i] - probs[i - 1];
            sumd  += d; sumd2 += d*d; nd++;
        }
    }
    if (nd < min_diffs || nd < 2) return std::numeric_limits<double>::quiet_NaN();
    const double meand   = sumd / (double)nd;
    const double diffVar = (sumd2 - (double)nd * meand * meand) / (double)(nd - 1);
    double noiseV = 0.5 * diffVar;
    if (noiseV < eps) noiseV = eps;
    double signalV = totalV - noiseV;
    if (signalV < eps) signalV = eps;
    return std::log2(signalV / noiseV);
}
