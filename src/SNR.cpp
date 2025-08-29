#include <algorithm> // std::max
#include <cmath>     // std::floor, std::log2
#include <limits>    // std::numeric_limits
#include <Rcpp.h>    // Rcpp::NumericVector / IntegerVector
using namespace Rcpp;

//' @title Internal: compute per-read SNR (NA-gap aware)
//'
//' @description
//' \code{compute_snr_na_gap_aware} estimates signal-to-noise ratio (SNR),
//' signal variance, and noise variance from a vector of modification
//' probabilities, allowing for missing values.
//'
//' @details
//' Implements a gap-aware estimate of SNR:
//' \enumerate{
//'   \item **Total variance** = \code{var(x, na.rm = TRUE)}
//'   \item **Noise variance** ≈ \code{0.5 * Var(Δx)}, where Δx are lag-1
//'         differences that may skip up to \code{k} missing values
//'   \item A noise floor is imposed: \code{noise = pmax(noise, b0 + b1 * mean(x))},
//'         with \code{b0}, \code{b1} from a robust linear fit
//'   \item **Signal variance** = \code{pmax(total - noise, eps)}
//'   \item **SNR** = \code{log2(signal / noise)}
//' }
//'
//'
//' @param probs Numeric vector of modification probabilities per read position.
//' @param read_pos Integer vector of read positions (same length as \code{probs}).
//' @param k Integer, maximum gap size tolerated when computing Δx.
//' @param min_diffs Integer, minimum number of differences required (default -1 = no requirement).
//' @param eps Numeric, small positive floor for signal variance.
//' @param b0 Numeric intercept for noise floor.
//' @param b1 Numeric slope for noise floor.
//'
//' @return A numeric scalar giving the estimated SNR (log2 scale).
//' Returns \code{NaN} if input length < 2 or lengths do not match.
//'
//' @examples
//' ## returns NaN because input length < 2
//' compute_snr_na_gap_aware(0.5, 1L, k=2L, min_diffs=-1L, eps=1e-3, b0=0, b1=0)
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
double compute_snr_na_gap_aware(const Rcpp::NumericVector& probs,
                                const Rcpp::IntegerVector& read_pos,
                                int k,
                                int min_diffs,
                                double eps,
                                double b0,
                                double b1) {
    const int n = probs.size();
    if (n < 2 || read_pos.size() != n) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (min_diffs < 0) {
        min_diffs = std::max(16, (int)std::floor(0.05 * (double)n));
    }
    if (k < 0) k = 0;
    
    // total variance
    double sumx = 0.0, sumx2 = 0.0;
    for (int i = 0; i < n; ++i) { sumx += probs[i]; sumx2 += probs[i]*probs[i]; }
    const double meanx  = sumx / (double)n;
    const double totalV = (sumx2 - (double)n * meanx * meanx) / (double)(n - 1);
    
    // 0.5 * Var(diff) over gaps <= k
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
    
    // optional robust noise floor: noise := max(noise, b0 + b1 * meanx)
    if (std::isfinite(b0) && std::isfinite(b1)) {
        const double fitted = b0 + b1 * meanx;
        if (noiseV < fitted) noiseV = fitted;
    }
    
    if (noiseV < eps) noiseV = eps;
    double signalV = totalV - noiseV;
    if (signalV < eps) signalV = eps;
    
    return std::log2(signalV / noiseV);
}
