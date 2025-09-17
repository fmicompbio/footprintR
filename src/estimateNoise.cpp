#include <algorithm>
#include <limits>
#include <Rcpp.h>

//' @title Estimation of noise variance for Time Series.
//' @description Given a Time Series signal, potentially with missing/unobserved values,
//' this function estimates the level of noise under the assumption of
//' local continuity/smoothness
//'
//' @details
//' **Noise variance** ≈ \code{0.5 * Var(Δx)} where Δx are lag-1 differences that may skip
//'  up to *k* missing values. This follows from error propagation and the assumption
//'  of low varying x in adjacent measurements.
//'
//' @param probs Numeric vector of (observed) Time Series measurements.
//' @param read_pos Integer vector of measurement positions in the Time Series (same length as probs).
//' @param k Integer, maximum gap size tolerated when computing Δx.
//' @param min_diffs Integer, minimum number of lag-1 differences required to produce an estimate (default -1 = auto).
//'
//' @return A vector with the following items:
//' \enumerate{
//'   \item **Mean signal**
//'   \item **Total variance**
//'   \item **Noise variance estimate**
//'   \item **Number of adjacent positions used for the estimate**
//' }
//'
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector estimateNoise(const Rcpp::NumericVector& probs,
                                  const Rcpp::IntegerVector& read_pos,
                                  int k,
                                  int min_diffs) {
    const int n = probs.size();
    Rcpp::NumericVector out(4, NA_REAL); // meanx, totalV, noiseV_raw, ndiffs
    out.names() = Rcpp::CharacterVector::create("mean", "total", "noise_raw", "ndiffs");

    if (n < 2 || read_pos.size() != n) {
        return out;
    }
    if (min_diffs < 0) {
        min_diffs = std::max(16, (int)std::floor(0.05 * (double)n));
    }
    if (k < 0) {
        k = 0;
    }

    // mean & total variance
    double sumx = 0.0, sumx2 = 0.0;
    for (int i = 0; i < n; i++) {
        sumx += probs[i];
        sumx2 += probs[i] * probs[i];
    }
    const double meanx = sumx / (double)n;
    const double totalV = (sumx2 - (double)n * meanx * meanx) / (double)(n - 1);
    out[0] = meanx;
    out[1] = totalV;

    // 0.5 * Var(diff) over gaps <= k
    int nd = 0;
    double sumd = 0.0, sumd2 = 0.0;
    for (int i = 1; i < n; i++) {
        const int gap = (read_pos[i] - read_pos[i - 1]) - 1;
        if (gap <= k) {
            const double d = probs[i] - probs[i - 1];
            sumd += d;
            sumd2 += d * d;
            nd++;
        }
    }

    if (nd < min_diffs || nd < 2) {
        return out;
    }
    const double meand = sumd / (double)nd;
    const double diffVar = (sumd2 - (double)nd * meand * meand) / (double)(nd - 1);

    out[2] = 0.5 * diffVar; // noise_raw
    out[3] = (double)nd;
    return out;
}
