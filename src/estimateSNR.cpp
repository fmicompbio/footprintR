#include <Rcpp.h>
#include <algorithm> // std::max
#include <cmath>     // std::log2
#include <limits>
#include <string>
using namespace Rcpp;

static inline std::string to_lower(const std::string& s) {
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(), ::tolower);
    return r;
}


 //' @title Estimate a Time Series Signal to Noise Ratio (SNR) given
 //' total and noise variance components and (optionally) a background Noise model.
 //'
 //' @description
 //' Given Time series variance components (total variance, raw noise variance),
 //' compute the final noise, signal, and SNR. For noise 
 //' estimation use one of:
 //' \enumerate{
 //'   \item **A Raw noise variance** e.g the result of \code{estimateNoise}
 //'   \item **A Background noise model** based on linear predictors
 //'   \item **A Flooring rule**: final noise = max(raw, background)
 //' }
 //' The Background noise model is a *general* linear noise baseline:
 //' \code{baseline = sum_i beta[i] * feat[i]}.
 //'
 //' @param totalVar Total time series variance. Typically calculated with \code{estimateNoise}
 //' @param noiseRaw noise variance estimate. Typically calculated with \code{estimateNoise}
 //' @param eps Numeric. Lower bound applied to both noise and signal.
 //' @param betas,features Coefficients (\eqn{betas}) and predictors (\eqn{features})
 //'   for the background noise model. Ignored when \code{noise_mode = "raw"}.
 //'   Both vectors must have the same length when used. The predictors form a
 //'   *design vector*: include a constant \code{1} if the model has an intercept.
 //'   For example, for \eqn{baseline = b0 + b1 * mean}, use
 //'   \code{betas = c(b0, b1)}, \code{features = c(1, mean)}.
 //' @param noise_mode Character scalar: one of \code{"raw"}, \code{"model"}, \code{"floor"}.
 //'        - \code{"raw"}   → use \code{noise_raw}
 //'        - \code{"model"} → use \code{baseline = sum(betas * features)}
 //'        - \code{"floor"} → use \code{max(noise_raw, baseline)}
 //'
 //' @return Named numeric vector with elements:
 //'         - snr      (log2 scale)
 //'         - signal
 //'         - noise
 //'         - baseline (sum(betas * features); NA if not applied)
 //'
 //' @noRd
 //' @keywords internal
 // [[Rcpp::export]]
 Rcpp::NumericVector estimateSNR(double totalVar,
                                 double noiseRaw,
                                 double eps,
                                 const Rcpp::NumericVector& betas,
                                 const Rcpp::NumericVector& features,
                                 const std::string& noise_mode = "floor") {
     NumericVector out(4, NA_REAL);
     out.names() = CharacterVector::create("snr", "signal", "noise", "baseline");
     
     const std::string mode = to_lower(noise_mode);
     
     if (!R_finite(totalVar)) return out;
     
     // For raw mode, both totalVar and noiseRaw must be finite
     if (mode == "raw" && !R_finite(noiseRaw)) {
         return out;
     }
     
     // Compute baseline if needed
     double baseline = std::numeric_limits<double>::quiet_NaN();
     const bool need_baseline = (mode == "model" || mode == "floor");
     
     if (need_baseline) {
         if (betas.size() != features.size() || betas.size() == 0) {
             stop("For noise_mode='%s', 'betas' and 'features' must be the same non-zero length.", mode.c_str());
         }
         double accum = 0.0;
         for (int i = 0; i < betas.size(); ++i) {
             const double b = betas[i];
             const double z = features[i];
             if (!R_finite(b) || !R_finite(z)) {
                 baseline = std::numeric_limits<double>::quiet_NaN();
                 break;
             }
             accum += b * z;
         }
         baseline = accum;
         if (!R_finite(baseline)) {
             // cannot proceed if baseline required but not finite
             return out;
         }
     }
     
     // Select noise according to mode
     double noiseV;
     if (mode == "raw") {
         noiseV = noiseRaw;
     } else if (mode == "model") {
         noiseV = baseline;
     } else if (mode == "floor") {
         noiseV = std::max(noiseRaw, baseline);
     } else {
         stop("Invalid noise_mode: '%s'. Use 'raw', 'model', or 'floor'.", noise_mode.c_str());
     }
     
     // Apply eps to noise and ensure finite
     if (!R_finite(noiseV)) return out;
     if (noiseV < eps) noiseV = eps;
     
     // Signal and SNR
     double signalV = totalVar - noiseV;
     if (!R_finite(signalV)) return out;
     if (signalV < eps) signalV = eps;
     
     const double snr = std::log2(signalV / noiseV);
     
     out[0] = snr;
     out[1] = signalV;
     out[2] = noiseV;
     out[3] = need_baseline ? baseline : NA_REAL;
     return out;
 }
